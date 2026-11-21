use rayon::prelude::*;

use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::seq::{Genome, SEQ_SEP};

/// Default GECKO scoring parameters (+4 match, -4 mismatch, -1 ambiguous 'N').
#[derive(Debug, Clone, Copy)]
pub struct Scoring {
    pub match_score: i64,
    pub mismatch_penalty: i64,
    pub n_penalty: i64,
}

impl Default for Scoring {
    fn default() -> Self {
        Self {
            match_score: 4,
            mismatch_penalty: 4,
            n_penalty: 1,
        }
    }
}

/// Parameters controlling bidirectional ungapped HSP extension (`FragHits`).
#[derive(Debug, Clone, Copy)]
pub struct ExtendOptions {
    /// Minimum fragment length (in bp if `fixed_length == true`, or percentage of min seq len).
    pub min_length: u64,
    /// Minimum similarity threshold (0..=100, strict `>` matching GECKO).
    pub min_similarity: f32,
    /// Seed word length used to produce the input hits.
    pub word_length: usize,
    /// If true, `min_length` is in bp; if false, percentage of shorter sequence length.
    pub fixed_length: bool,
    /// Strand being compared (`'f'` for forward, `'r'` for reverse-complemented Y).
    pub strand: char,
    /// Alignment scoring scheme.
    pub scoring: Scoring,
}

/// 1000 x 100 length-by-similarity coverage histogram (matching GECKO's `.MAT` file,
/// with similarity clamped to `0..=99` to fix the original C off-by-one overflow at 100%).
#[derive(Debug, Clone)]
pub struct CoverageMatrix {
    pub bins: Vec<[u64; 100]>,
}

impl Default for CoverageMatrix {
    fn default() -> Self {
        Self::new()
    }
}

impl CoverageMatrix {
    pub fn new() -> Self {
        Self {
            bins: vec![[0u64; 100]; 1000],
        }
    }

    #[inline]
    pub fn record(&mut self, length: u64, similarity: f32) {
        let row = (length as usize).min(999);
        let col = (similarity.max(0.0) as usize).min(99);
        self.bins[row][col] += 1;
    }

    pub fn merge_from(&mut self, other: &Self) {
        for (dst_row, src_row) in self.bins.iter_mut().zip(other.bins.iter()) {
            for (d, &s) in dst_row.iter_mut().zip(src_row.iter()) {
                *d += s;
            }
        }
    }
}

/// Per-batch extension output.
#[derive(Debug, Default)]
pub struct ExtendBatchOutput {
    pub fragments: Vec<Fragment>,
    pub hits_total: u64,
    pub hits_used: u64,
    pub coverage: CoverageMatrix,
}

/// Extend a diagonal-sorted slice of [`Hit`]s in parallel across diagonal groups.
///
/// Assumes `hits` is sorted by `(diag, pos_x)` and contains complete diagonals (never splits
/// a single diagonal across separate calls).
pub fn extend_hits_parallel(
    hits: &[Hit],
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
) -> ExtendBatchOutput {
    if hits.is_empty() {
        return ExtendBatchOutput::default();
    }

    // Partition `hits` into coarse chunks aligned to diagonal boundaries so `rayon` threads
    // process independent diagonals with zero synchronization.
    let target_chunks = (rayon::current_num_threads() * 4).max(1);
    let min_chunk = (hits.len() / target_chunks).max(256);

    let mut ranges: Vec<(usize, usize)> = Vec::with_capacity(target_chunks + 1);
    let mut start = 0usize;
    while start < hits.len() {
        let mut end = (start + min_chunk).min(hits.len());
        while end < hits.len() && hits[end].diag == hits[end - 1].diag {
            end += 1;
        }
        ranges.push((start, end));
        start = end;
    }

    let partial: Vec<ExtendBatchOutput> = ranges
        .into_par_iter()
        .map(|(s, e)| extend_diagonal_slice(&hits[s..e], genome_x, genome_y, opts))
        .collect();

    let mut merged = ExtendBatchOutput::default();
    let total_frags: usize = partial.iter().map(|p| p.fragments.len()).sum();
    merged.fragments.reserve(total_frags);

    for mut p in partial {
        merged.hits_total += p.hits_total;
        merged.hits_used += p.hits_used;
        merged.coverage.merge_from(&p.coverage);
        merged.fragments.append(&mut p.fragments);
    }

    merged
}

fn extend_diagonal_slice(
    hits: &[Hit],
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
) -> ExtendBatchOutput {
    let mut out = ExtendBatchOutput::default();
    if hits.is_empty() {
        return out;
    }

    let mut last_diag = hits[0].diag;
    let mut last_offset: i64 = hits[0].pos_x as i64 - 1;

    for &h in hits {
        out.hits_total += 1;
        if last_diag != h.diag {
            last_diag = h.diag;
            last_offset = h.pos_x as i64 - 1;
        }
        if last_offset > h.pos_x as i64 {
            // Covered by a previously extended fragment on the same diagonal.
            continue;
        }

        out.hits_used += 1;
        if let Some(frag) = extend_single_hit(h, genome_x, genome_y, opts, &mut out.coverage) {
            last_offset = (h.pos_x + frag.length) as i64;
            out.fragments.push(frag);
        }
    }

    out
}

/// Bidirectional ungapped extension from a single seed hit, faithful to GECKO's `FragFromHit`.
fn extend_single_hit(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    coverage: &mut CoverageMatrix,
) -> Option<Fragment> {
    let data_x = &genome_x.data;
    let data_y = &genome_y.data;
    let nx = data_x.len() as u64;
    let ny = data_y.len() as u64;

    let pos_x = h.pos_x;
    let pos_y = h.pos_y();
    let wl = opts.word_length as u64;

    if pos_x + wl > nx || pos_y + wl > ny {
        return None;
    }

    let seq_x = genome_x.seq_index_for_data_pos(pos_x);
    let seq_y = genome_y.seq_index_for_data_pos(pos_y);

    let min_length: u64 = if opts.fixed_length {
        opts.min_length
    } else {
        let lx = genome_x.records[seq_x].len;
        let ly = genome_y.records[seq_y].len;
        (lx.min(ly) as f64 * (opts.min_length as f64 / 100.0)) as u64
    };

    let point = opts.scoring.match_score;
    let mismatch = opts.scoring.mismatch_penalty;
    let n_pen = opts.scoring.n_penalty;

    // Forward extension from right edge of seed
    let ldiag = (nx - pos_x).min(ny - pos_y);
    let mut xfil = pos_x + wl;
    let mut ycol = pos_y + wl;
    let mut frag_len = wl;

    // Note: initialize inclusive right endpoint to the last base of the seed (`pos + wl - 1`).
    let mut xfil_max = pos_x + wl - 1;
    let mut ycol_max = pos_y + wl - 1;
    let mut n_identities = wl;
    let mut max_identities = wl;
    let mut fscore = point * wl as i64;
    let mut fscore_max = fscore;

    while frag_len < ldiag {
        let vx = data_x[xfil as usize].to_ascii_uppercase();
        let vy = data_y[ycol as usize].to_ascii_uppercase();
        if vx == SEQ_SEP || vy == SEQ_SEP {
            break;
        }

        if vx == b'N' || vy == b'N' {
            fscore -= n_pen;
        } else if vx == vy {
            fscore += point;
            n_identities += 1;
            if fscore > fscore_max {
                fscore_max = fscore;
                xfil_max = xfil;
                ycol_max = ycol;
                max_identities = n_identities;
            }
        } else {
            fscore -= mismatch;
        }

        xfil += 1;
        ycol += 1;
        frag_len += 1;
        if fscore < 0 {
            break;
        }
    }

    // Backward extension from left edge of seed
    let ldiag2 = pos_x.min(pos_y);
    let mut xfil_max2 = pos_x;
    let mut ycol_max2 = pos_y;
    fscore = fscore_max;
    n_identities = max_identities;

    if pos_x > 0 && pos_y > 0 {
        let mut xfil2 = pos_x - 1;
        let mut ycol2 = pos_y - 1;
        let mut back_len = 0u64;

        while back_len < ldiag2 {
            let vx = data_x[xfil2 as usize].to_ascii_uppercase();
            let vy = data_y[ycol2 as usize].to_ascii_uppercase();
            if vx == SEQ_SEP || vy == SEQ_SEP {
                break;
            }

            if vx == b'N' || vy == b'N' {
                fscore -= n_pen;
            } else if vx == vy {
                fscore += point;
                n_identities += 1;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max2 = xfil2;
                    ycol_max2 = ycol2;
                    max_identities = n_identities;
                }
            } else {
                fscore -= mismatch;
            }

            back_len += 1;
            if fscore < 0 || xfil2 == 0 || ycol2 == 0 {
                break;
            }
            xfil2 -= 1;
            ycol2 -= 1;
        }
    }

    // Convert data coordinates (with '*' separators) to GECKO concatenated nucleotide coordinates.
    let x_start = xfil_max2 - seq_x as u64;
    let y_start = ycol_max2 - seq_y as u64;
    let x_end = xfil_max - seq_x as u64;
    let y_end = ycol_max - seq_y as u64;
    let length = x_end - x_start + 1;
    let score = fscore_max.max(0) as u64;
    let max_possible = (length as f64) * (point as f64);
    let similarity = ((score as f64) * 100.0 / max_possible) as f32;

    let out_seq_y = if opts.strand == 'f' {
        seq_y as u64
    } else {
        (genome_y.num_records() - 1 - seq_y) as u64
    };

    coverage.record(length, similarity);

    if length > min_length && similarity > opts.min_similarity {
        Some(Fragment {
            diag: h.diag,
            x_start,
            y_start,
            x_end,
            y_end,
            length,
            ident: max_identities,
            score,
            similarity,
            seq_x: seq_x as u64,
            seq_y: out_seq_y,
            block: 0,
            strand: opts.strand,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extends_perfect_and_mismatched_hsp_without_oob() {
        let seq = "ACGTACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGTACGTACGT";
        let gx = Genome::from_records("x", [("x0", seq)]);
        let gy = Genome::from_records("y", [("y0", seq)]);
        let hit = Hit { diag: 0, pos_x: 10 };

        let out = extend_hits_parallel(
            &[hit],
            &gx,
            &gy,
            ExtendOptions {
                min_length: 20,
                min_similarity: 60.0,
                word_length: 12,
                fixed_length: true,
                strand: 'f',
                scoring: Scoring::default(),
            },
        );

        assert_eq!(out.fragments.len(), 1);
        let f = &out.fragments[0];
        assert_eq!(f.x_start, 0);
        assert_eq!(f.x_end, (seq.len() - 1) as u64);
        assert_eq!(f.length, seq.len() as u64);
        assert_eq!(f.ident, seq.len() as u64);
        assert!((f.similarity - 100.0).abs() < 1e-4);
        // Verified no OOB panic on 100% similarity bin!
        assert_eq!(out.coverage.bins[seq.len()][99], 1);
    }
}
