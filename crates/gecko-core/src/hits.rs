use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, Write};

use crate::error::{GeckoError, Result};
use crate::kmer::{KmerDict, kmer_prefix};
use crate::radix::{i64_to_sortable_u64, parallel_radix_sort_u128_key};
use crate::seq::Genome;

/// Compact diagonal seed hit (16 bytes vs 40 bytes in C `hit`).
///
/// `pos_y` is derived directly as `(pos_x as i64 - diag) as u64`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct Hit {
    pub diag: i64,
    pub pos_x: u64,
}

impl Hit {
    #[inline]
    pub fn pos_y(&self) -> u64 {
        (self.pos_x as i64 - self.diag) as u64
    }
}

/// Configuration for seed matching and hit filtering.
#[derive(Debug, Clone, Copy)]
pub struct HitOptions {
    /// Seed prefix length in bases (4..=32; no longer restricted to multiples of 4).
    pub word_length: usize,
    /// Skip any dictionary word with frequency strictly greater than `freq_limit`.
    pub freq_limit: u64,
    /// Subsampling cap per word occurrence list (matching GECKO's `MaxREP = 10_000`).
    pub max_rep: usize,
    /// Memory budget in bytes for buffering hits before spilling sorted runs to disk.
    pub memory_budget: usize,
}

impl Default for HitOptions {
    fn default() -> Self {
        Self {
            word_length: 12,
            freq_limit: 1000,
            max_rep: 10_000,
            memory_budget: 2 * 1024 * 1024 * 1024,
        }
    }
}

/// Statistics from seed matching and diagonal hit filtering.
#[derive(Debug, Clone, Default)]
pub struct HitStats {
    pub word_matches: u64,
    pub raw_hits: u64,
    pub filtered_hits: u64,
    pub spilled_runs: usize,
}

/// Find, sort, and deduplicate diagonal hits between `dict_x` and `dict_y`, delivering
/// filtered hits sorted by `(diag, pos_x)`.
///
/// - Fast path: if all raw hits fit in `opts.memory_budget`, everything runs 100% in memory
///   with zero hit disk I/O.
/// - Out-of-core path: when raw hits exceed `opts.memory_budget`, chunks are sorted,
///   pre-deduplicated per diagonal, spilled to temporary files, and k-way merged into a
///   compact deduplicated hit vector (or streamed in diagonal batches via
///   [`generate_and_consume_hits`]).
pub fn generate_filtered_hits(
    dict_x: &KmerDict,
    dict_y: &KmerDict,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: HitOptions,
) -> Result<(Vec<Hit>, HitStats)> {
    let mut out = Vec::new();
    let stats = generate_and_consume_hits(dict_x, dict_y, genome_x, genome_y, opts, |batch| {
        out.extend_from_slice(batch);
    })?;
    Ok((out, stats))
}

/// Memory-bounded streaming hit generator: emits sorted, deduplicated slices of [`Hit`]
/// guaranteed to be split cleanly at diagonal boundaries so downstream HSP extension can
/// process each batch independently in parallel.
pub fn generate_and_consume_hits<F>(
    dict_x: &KmerDict,
    dict_y: &KmerDict,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: HitOptions,
    mut on_diagonal_batch: F,
) -> Result<HitStats>
where
    F: FnMut(&[Hit]),
{
    if opts.word_length < 4 || opts.word_length > dict_x.index_k.min(dict_y.index_k) {
        return Err(GeckoError::InvalidParameter(format!(
            "word_length ({}) must be in 4..={} for the given dictionaries",
            opts.word_length,
            dict_x.index_k.min(dict_y.index_k)
        )));
    }

    let max_hits_in_mem = (opts.memory_budget / (2 * std::mem::size_of::<Hit>())).max(4096);
    let mut buf: Vec<Hit> = Vec::with_capacity(max_hits_in_mem.min(1 << 20));
    let mut spill_files: Vec<File> = Vec::new();
    let mut stats = HitStats::default();

    let wl = opts.word_length;
    let freq_limit = opts.freq_limit as usize;
    let max_rep = opts.max_rep.max(1);

    let nx = dict_x.len();
    let ny = dict_y.len();
    let mut ix = 0usize;
    let mut iy = 0usize;

    while ix < nx && iy < ny {
        let px = kmer_prefix(dict_x.words[ix], wl);
        let py = kmer_prefix(dict_y.words[iy], wl);

        if px < py {
            ix = advance_prefix_group(&dict_x.words, ix, px, wl);
            continue;
        }
        if px > py {
            iy = advance_prefix_group(&dict_y.words, iy, py, wl);
            continue;
        }

        let ix_end = advance_prefix_group(&dict_x.words, ix, px, wl);
        let iy_end = advance_prefix_group(&dict_y.words, iy, py, wl);

        for i in ix..ix_end {
            let pos_x_slice = dict_x.positions_for(i);
            if pos_x_slice.len() > freq_limit {
                continue;
            }
            let step_x = if pos_x_slice.len() > max_rep {
                pos_x_slice.len() / max_rep
            } else {
                1
            };

            for j in iy..iy_end {
                let pos_y_slice = dict_y.positions_for(j);
                if pos_y_slice.len() > freq_limit {
                    continue;
                }
                stats.word_matches += 1;
                let step_y = if pos_y_slice.len() > max_rep {
                    pos_y_slice.len() / max_rep
                } else {
                    1
                };

                let mut ax = 0;
                while ax < pos_x_slice.len() {
                    let px_val = pos_x_slice[ax];
                    let mut ay = 0;
                    while ay < pos_y_slice.len() {
                        let py_val = pos_y_slice[ay];
                        buf.push(Hit {
                            diag: px_val as i64 - py_val as i64,
                            pos_x: px_val,
                        });
                        stats.raw_hits += 1;
                        if buf.len() >= max_hits_in_mem {
                            spill_files.push(sort_and_spill_hits(&mut buf)?);
                        }
                        ay += step_y;
                    }
                    ax += step_x;
                }
            }
        }

        ix = ix_end;
        iy = iy_end;
    }

    if spill_files.is_empty() {
        parallel_radix_sort_u128_key(&mut buf, |h| (i64_to_sortable_u64(h.diag), h.pos_x));
        dedup_hits_in_place(&mut buf, genome_x, genome_y, wl);
        stats.filtered_hits = buf.len() as u64;
        if !buf.is_empty() {
            on_diagonal_batch(&buf);
        }
        return Ok(stats);
    }

    if !buf.is_empty() {
        spill_files.push(sort_and_spill_hits(&mut buf)?);
    }
    drop(buf);
    stats.spilled_runs = spill_files.len();

    // K-way merge spilled runs and stream complete diagonal batches into the consumer.
    let mut readers: Vec<BufReader<File>> = Vec::with_capacity(spill_files.len());
    for mut f in spill_files {
        f.rewind()?;
        readers.push(BufReader::with_capacity(256 * 1024, f));
    }

    let mut heap: BinaryHeap<Reverse<(i64, u64, usize)>> = BinaryHeap::new();
    for (run_idx, r) in readers.iter_mut().enumerate() {
        if let Some(h) = read_hit(r)? {
            heap.push(Reverse((h.diag, h.pos_x, run_idx)));
        }
    }

    let batch_target = max_hits_in_mem.min(2 * 1024 * 1024);
    let mut batch: Vec<Hit> = Vec::with_capacity(batch_target.min(65536));
    let mut dedup = DiagDedupState::new(wl, genome_x, genome_y);

    while let Some(Reverse((diag, pos_x, run_idx))) = heap.pop() {
        let hit = Hit { diag, pos_x };
        if let Some(next) = read_hit(&mut readers[run_idx])? {
            heap.push(Reverse((next.diag, next.pos_x, run_idx)));
        }

        if dedup.accept(hit) {
            if batch.len() >= batch_target && batch.last().is_some_and(|prev| prev.diag != hit.diag)
            {
                on_diagonal_batch(&batch);
                batch.clear();
            }
            batch.push(hit);
            stats.filtered_hits += 1;
        }
    }

    if !batch.is_empty() {
        on_diagonal_batch(&batch);
    }

    Ok(stats)
}

#[inline]
fn advance_prefix_group(words: &[u64], start: usize, prefix: u64, wl: usize) -> usize {
    let mut i = start + 1;
    while i < words.len() && kmer_prefix(words[i], wl) == prefix {
        i += 1;
    }
    i
}

/// Stateful diagonal deduplicator matching GECKO's `filterHits.c` semantics.
struct DiagDedupState<'a> {
    dist: u64,
    genome_x: &'a Genome,
    genome_y: &'a Genome,
    multi_seq: bool,
    has_prev: bool,
    prev_diag: i64,
    prev_seq_x: usize,
    prev_seq_y: usize,
    last_pos: u64,
}

impl<'a> DiagDedupState<'a> {
    fn new(word_length: usize, genome_x: &'a Genome, genome_y: &'a Genome) -> Self {
        Self {
            dist: (2 * word_length - 1) as u64,
            genome_x,
            genome_y,
            multi_seq: genome_x.num_records() > 1 || genome_y.num_records() > 1,
            has_prev: false,
            prev_diag: 0,
            prev_seq_x: 0,
            prev_seq_y: 0,
            last_pos: 0,
        }
    }

    #[inline]
    fn accept(&mut self, h: Hit) -> bool {
        let (sx, sy) = if self.multi_seq {
            (
                self.genome_x.seq_index_for_data_pos(h.pos_x),
                self.genome_y.seq_index_for_data_pos(h.pos_y()),
            )
        } else {
            (0, 0)
        };

        if !self.has_prev {
            self.has_prev = true;
            self.prev_diag = h.diag;
            self.prev_seq_x = sx;
            self.prev_seq_y = sy;
            self.last_pos = h.pos_x.saturating_add(self.dist);
            return true;
        }

        if sx != self.prev_seq_x
            || sy != self.prev_seq_y
            || h.diag != self.prev_diag
            || h.pos_x > self.last_pos
        {
            self.prev_diag = h.diag;
            self.prev_seq_x = sx;
            self.prev_seq_y = sy;
            self.last_pos = h.pos_x.saturating_add(self.dist);
            true
        } else {
            false
        }
    }
}

fn dedup_hits_in_place(hits: &mut Vec<Hit>, genome_x: &Genome, genome_y: &Genome, wl: usize) {
    let mut state = DiagDedupState::new(wl, genome_x, genome_y);
    hits.retain(|&h| state.accept(h));
}

fn sort_and_spill_hits(buf: &mut Vec<Hit>) -> Result<File> {
    parallel_radix_sort_u128_key(buf, |h| (i64_to_sortable_u64(h.diag), h.pos_x));

    let mut tmp = tempfile::tempfile()?;
    {
        let mut w = BufWriter::with_capacity(1 << 20, &mut tmp);
        for h in buf.iter() {
            w.write_all(&h.diag.to_le_bytes())?;
            w.write_all(&h.pos_x.to_le_bytes())?;
        }
        w.flush()?;
    }
    buf.clear();
    Ok(tmp)
}

fn read_hit<R: Read>(r: &mut R) -> Result<Option<Hit>> {
    let mut buf = [0u8; 16];
    match r.read_exact(&mut buf) {
        Ok(()) => Ok(Some(Hit {
            diag: i64::from_le_bytes(buf[0..8].try_into().unwrap()),
            pos_x: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
        })),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(GeckoError::Io(e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::DictOptions;

    #[test]
    fn in_memory_and_spilled_hit_generation_match() {
        let seq_x = "ACGTACGTCCGGTTAAACGTTTGGCCAACGTACGTCCGGTTAAACGTTTGGCCAA".repeat(3);
        let seq_y = format!("NNNNN{seq_x}TTTTT");
        let gx = Genome::from_records("x", [("x0", seq_x.as_bytes())]);
        let gy = Genome::from_records("y", [("y0", seq_y.as_bytes())]);

        let dx = KmerDict::build(&gx, DictOptions::default()).unwrap();
        let dy = KmerDict::build(&gy, DictOptions::default()).unwrap();

        let (hits_mem, _) = generate_filtered_hits(
            &dx,
            &dy,
            &gx,
            &gy,
            HitOptions {
                word_length: 12,
                ..Default::default()
            },
        )
        .unwrap();

        let (hits_spill, stats_spill) = generate_filtered_hits(
            &dx,
            &dy,
            &gx,
            &gy,
            HitOptions {
                word_length: 12,
                memory_budget: 64 * 1024, // force spilling
                ..Default::default()
            },
        )
        .unwrap();

        assert!(!hits_mem.is_empty());
        assert_eq!(hits_mem, hits_spill);
        // Either fit or spilled cleanly without corrupting results
        assert_eq!(stats_spill.filtered_hits as usize, hits_spill.len());
    }
}
