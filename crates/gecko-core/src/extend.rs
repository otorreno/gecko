use rayon::prelude::*;

use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::seq::{Genome, PackedSeq, SEQ_SEP};

/// Default GECKO scoring parameters (+4 match, -4 mismatch, -1 ambiguous 'N').
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

impl Scoring {
    /// Whether this scoring matches GECKO's default `+4/-4/-1` scheme accelerated by 2-bit/AVX2.
    #[inline(always)]
    pub fn is_default_gecko(&self) -> bool {
        self.match_score == 4 && self.mismatch_penalty == 4 && self.n_penalty == 1
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
    /// If true, the seed window `[pos_x .. pos_x + word_length)` may contain mismatches at
    /// don't-care positions of a spaced seed, so initial identity/score is computed from 2-bit words.
    pub spaced_seed: bool,
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

    #[inline(always)]
    pub fn record(&mut self, length: u64, similarity: f32) {
        let row = (length as usize).min(999);
        let col = (similarity.max(0.0) as usize).min(99);
        debug_assert_eq!(self.bins.len(), 1000);
        // Safety: `bins` always has length 1000 and `row <= 999`, `col <= 99`.
        unsafe {
            *self.bins.get_unchecked_mut(row).get_unchecked_mut(col) += 1;
        }
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

/// Outcome of a single seed hit extension, including the forward peak `xfil_max` needed
/// by the streaming diagonal cache and exact rejected-span skip.
#[derive(Debug, Clone)]
pub(crate) struct ExtendOutcome {
    pub fragment: Option<Fragment>,
    pub xfil_max: u64,
}

/// Reusable scratch buffers for [`exact_rejected_skip`].
#[derive(Debug, Default)]
pub(crate) struct SkipScratch {
    blk_min: Vec<i64>,
    blk_max: Vec<i64>,
    suf_min: Vec<i64>,
}

/// Runtime check for AVX2 support.
#[inline(always)]
pub(crate) fn has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Runtime check for BMI2 (`pext`) support.
#[inline(always)]
pub(crate) fn has_bmi2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("bmi2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Runtime check for hardware `popcnt` support.
#[inline(always)]
pub(crate) fn has_popcnt() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("popcnt")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
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

    let use_avx2 = has_avx2() && opts.scoring.is_default_gecko();
    let use_bmi2 = has_bmi2();

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
        .map(|(s, e)| {
            extend_diagonal_slice(&hits[s..e], genome_x, genome_y, opts, use_avx2, use_bmi2)
        })
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
    use_avx2: bool,
    use_bmi2: bool,
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
        let outcome = extend_single_hit_outcome(
            h,
            genome_x,
            genome_y,
            opts,
            use_avx2,
            use_bmi2,
            &mut out.coverage,
        );
        if let Some(frag) = outcome.fragment {
            last_offset = (h.pos_x + frag.length) as i64;
            out.fragments.push(frag);
        }
    }

    out
}

/// Forward pass state returned by scalar / 2-bit AVX2 forward extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ForwardResult {
    xfil_max: u64,
    ycol_max: u64,
    fscore_max: i64,
    max_identities: u64,
}

/// Backward pass state returned by scalar / 2-bit AVX2 backward extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BackwardResult {
    xfil_max2: u64,
    ycol_max2: u64,
    fscore_max: i64,
    max_identities: u64,
}

/// Compute shared `special_runs` cursor (`runs[..c].end <= pos`).
#[inline(always)]
fn sp_cursor_for(runs: &[(u64, u64)], pos: u64) -> usize {
    if runs.is_empty() {
        0
    } else {
        runs.partition_point(|&(_, e)| e <= pos)
    }
}

/// Finalize `(fwd, bwd)` extension results into an [`ExtendOutcome`], deferring
/// `seq_index_for_data_pos` lookups until a fragment is accepted when `opts.fixed_length` is set.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn finalize_extend_outcome(
    h: Hit,
    pos_x: u64,
    pos_y: u64,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    fwd: ForwardResult,
    bwd: BackwardResult,
    coverage: &mut CoverageMatrix,
) -> ExtendOutcome {
    // `(fwd.xfil_max - seq_x) - (bwd.xfil_max2 - seq_x) + 1 == fwd.xfil_max - bwd.xfil_max2 + 1`.
    let length = fwd.xfil_max - bwd.xfil_max2 + 1;
    let score = bwd.fscore_max.max(0) as u64;
    let max_possible = (length as f64) * (opts.scoring.match_score as f64);
    let similarity = ((score as f64) * 100.0 / max_possible) as f32;

    coverage.record(length, similarity);

    let (min_length, cached_seqs) = if opts.fixed_length {
        (opts.min_length, None)
    } else {
        let sx = genome_x.seq_index_for_data_pos(pos_x);
        let sy = genome_y.seq_index_for_data_pos(pos_y);
        let lx = genome_x.records[sx].len;
        let ly = genome_y.records[sy].len;
        let ml = (lx.min(ly) as f64 * (opts.min_length as f64 / 100.0)) as u64;
        (ml, Some((sx, sy)))
    };

    if length > min_length && similarity > opts.min_similarity {
        let (seq_x, seq_y) = cached_seqs.unwrap_or_else(|| {
            (
                genome_x.seq_index_for_data_pos(pos_x),
                genome_y.seq_index_for_data_pos(pos_y),
            )
        });
        let x_start = bwd.xfil_max2 - seq_x as u64;
        let y_start = bwd.ycol_max2 - seq_y as u64;
        let x_end = fwd.xfil_max - seq_x as u64;
        let y_end = fwd.ycol_max - seq_y as u64;
        let out_seq_y = if opts.strand == 'f' {
            seq_y as u64
        } else {
            (genome_y.num_records() - 1 - seq_y) as u64
        };
        ExtendOutcome {
            fragment: Some(Fragment {
                diag: h.diag,
                x_start,
                y_start,
                x_end,
                y_end,
                length,
                ident: bwd.max_identities,
                score,
                similarity,
                seq_x: seq_x as u64,
                seq_y: out_seq_y,
                block: 0,
                strand: opts.strand,
            }),
            xfil_max: fwd.xfil_max,
        }
    } else {
        ExtendOutcome {
            fragment: None,
            xfil_max: fwd.xfil_max,
        }
    }
}

#[inline(always)]
fn seed_init_stats(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    scoring: Scoring,
    spaced_seed: bool,
) -> (u64, i64) {
    if !spaced_seed {
        return (wl, wl as i64 * scoring.match_score);
    }
    let span_mask = if wl >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * wl)) - 1
    };
    let diff = (seq_x.load_32(pos_x as usize) ^ seq_y.load_32(pos_y as usize)) & span_mask;
    let mm = ((diff | (diff >> 1)) & 0x5555_5555_5555_5555).count_ones() as u64;
    let ident = wl - mm;
    let score = ident as i64 * scoring.match_score - mm as i64 * scoring.mismatch_penalty;
    (ident, score)
}

/// Bidirectional ungapped extension from a single seed hit, faithful to GECKO's `FragFromHit`.
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_single_hit_outcome(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    use_avx2: bool,
    use_bmi2: bool,
    coverage: &mut CoverageMatrix,
) -> ExtendOutcome {
    let seq_x_pack = &genome_x.seq;
    let seq_y_pack = &genome_y.seq;
    let nx = seq_x_pack.len as u64;
    let ny = seq_y_pack.len as u64;

    let pos_x = h.pos_x;
    let pos_y = h.pos_y();
    let wl = opts.word_length as u64;

    if pos_x + wl > nx || pos_y + wl > ny {
        return ExtendOutcome {
            fragment: None,
            xfil_max: pos_x + wl.saturating_sub(1),
        };
    }

    let (seed_ident, seed_score) = seed_init_stats(
        seq_x_pack,
        seq_y_pack,
        pos_x,
        pos_y,
        wl,
        opts.scoring,
        opts.spaced_seed,
    );

    let (fwd, bwd) = if use_avx2 && opts.scoring.is_default_gecko() {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            let cx = sp_cursor_for(&seq_x_pack.special_runs, pos_x);
            let cy = sp_cursor_for(&seq_y_pack.special_runs, pos_y);
            if use_bmi2 && has_popcnt() {
                let f = extend_forward_2bit_avx2_bmi2(
                    seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score, cx, cy,
                );
                let b =
                    extend_backward_2bit_avx2_bmi2(seq_x_pack, seq_y_pack, pos_x, pos_y, f, cx, cy);
                (f, b)
            } else if has_popcnt() {
                let f = extend_forward_2bit_avx2_swar(
                    seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score, cx, cy,
                );
                let b =
                    extend_backward_2bit_avx2_swar(seq_x_pack, seq_y_pack, pos_x, pos_y, f, cx, cy);
                (f, b)
            } else {
                let f = extend_forward_2bit_portable(
                    seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score,
                );
                let b = extend_backward_2bit_portable(seq_x_pack, seq_y_pack, pos_x, pos_y, f);
                (f, b)
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = use_bmi2;
            let f = extend_forward_2bit_portable(
                seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score,
            );
            let b = extend_backward_2bit_portable(seq_x_pack, seq_y_pack, pos_x, pos_y, f);
            (f, b)
        }
    } else if opts.scoring.is_default_gecko() {
        let f = extend_forward_2bit_portable(
            seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score,
        );
        let b = extend_backward_2bit_portable(seq_x_pack, seq_y_pack, pos_x, pos_y, f);
        (f, b)
    } else {
        let f = extend_forward_scalar(
            seq_x_pack,
            seq_y_pack,
            pos_x,
            pos_y,
            wl,
            seed_ident,
            seed_score,
            opts.scoring,
        );
        let b = extend_backward_scalar(seq_x_pack, seq_y_pack, pos_x, pos_y, f, opts.scoring);
        (f, b)
    };

    finalize_extend_outcome(
        h, pos_x, pos_y, genome_x, genome_y, opts, fwd, bwd, coverage,
    )
}

/// Inlineable AVX2 + BMI2 + POPCNT single-hit extension with caller-supplied monotonic `sp_x` cursor.
///
/// # Safety
/// Caller must run inside a `#[target_feature(enable = "avx2,bmi2,popcnt")]` context and guarantee
/// `sp_x[..cx_x].end <= h.pos_x`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,bmi2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn extend_single_hit_outcome_avx2_bmi2(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    cx_x: usize,
    coverage: &mut CoverageMatrix,
) -> ExtendOutcome {
    let seq_x_pack = &genome_x.seq;
    let seq_y_pack = &genome_y.seq;
    let nx = seq_x_pack.len as u64;
    let ny = seq_y_pack.len as u64;

    let pos_x = h.pos_x;
    let pos_y = h.pos_y();
    let wl = opts.word_length as u64;

    if pos_x + wl > nx || pos_y + wl > ny {
        return ExtendOutcome {
            fragment: None,
            xfil_max: pos_x + wl.saturating_sub(1),
        };
    }

    let (seed_ident, seed_score) = seed_init_stats(
        seq_x_pack,
        seq_y_pack,
        pos_x,
        pos_y,
        wl,
        opts.scoring,
        opts.spaced_seed,
    );

    let (fwd, bwd) = if opts.scoring.is_default_gecko() {
        let cy = sp_cursor_for(&seq_y_pack.special_runs, pos_y);
        let f = unsafe {
            extend_forward_2bit_avx2_bmi2(
                seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score, cx_x, cy,
            )
        };
        let b = unsafe {
            extend_backward_2bit_avx2_bmi2(seq_x_pack, seq_y_pack, pos_x, pos_y, f, cx_x, cy)
        };
        (f, b)
    } else {
        let f = extend_forward_scalar(
            seq_x_pack,
            seq_y_pack,
            pos_x,
            pos_y,
            wl,
            seed_ident,
            seed_score,
            opts.scoring,
        );
        let b = extend_backward_scalar(seq_x_pack, seq_y_pack, pos_x, pos_y, f, opts.scoring);
        (f, b)
    };

    finalize_extend_outcome(
        h, pos_x, pos_y, genome_x, genome_y, opts, fwd, bwd, coverage,
    )
}

/// Inlineable AVX2 + POPCNT (SWAR pack) single-hit extension with caller-supplied monotonic `sp_x` cursor.
///
/// # Safety
/// Caller must run inside a `#[target_feature(enable = "avx2,popcnt")]` context and guarantee
/// `sp_x[..cx_x].end <= h.pos_x`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn extend_single_hit_outcome_avx2_swar(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    cx_x: usize,
    coverage: &mut CoverageMatrix,
) -> ExtendOutcome {
    let seq_x_pack = &genome_x.seq;
    let seq_y_pack = &genome_y.seq;
    let nx = seq_x_pack.len as u64;
    let ny = seq_y_pack.len as u64;

    let pos_x = h.pos_x;
    let pos_y = h.pos_y();
    let wl = opts.word_length as u64;

    if pos_x + wl > nx || pos_y + wl > ny {
        return ExtendOutcome {
            fragment: None,
            xfil_max: pos_x + wl.saturating_sub(1),
        };
    }

    let (seed_ident, seed_score) = seed_init_stats(
        seq_x_pack,
        seq_y_pack,
        pos_x,
        pos_y,
        wl,
        opts.scoring,
        opts.spaced_seed,
    );

    let (fwd, bwd) = if opts.scoring.is_default_gecko() {
        let cy = sp_cursor_for(&seq_y_pack.special_runs, pos_y);
        let f = unsafe {
            extend_forward_2bit_avx2_swar(
                seq_x_pack, seq_y_pack, pos_x, pos_y, wl, seed_ident, seed_score, cx_x, cy,
            )
        };
        let b = unsafe {
            extend_backward_2bit_avx2_swar(seq_x_pack, seq_y_pack, pos_x, pos_y, f, cx_x, cy)
        };
        (f, b)
    } else {
        let f = extend_forward_scalar(
            seq_x_pack,
            seq_y_pack,
            pos_x,
            pos_y,
            wl,
            seed_ident,
            seed_score,
            opts.scoring,
        );
        let b = extend_backward_scalar(seq_x_pack, seq_y_pack, pos_x, pos_y, f, opts.scoring);
        (f, b)
    };

    finalize_extend_outcome(
        h, pos_x, pos_y, genome_x, genome_y, opts, fwd, bwd, coverage,
    )
}

// ---------------------------------------------------------------------------
// Portable scalar extension (custom scoring + reference oracle for fuzz tests)
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn extend_forward_scalar(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    seed_ident: u64,
    seed_score: i64,
    scoring: Scoring,
) -> ForwardResult {
    let nx = seq_x.len as u64;
    let ny = seq_y.len as u64;
    let ldiag = (nx - pos_x).min(ny - pos_y);
    let mut xfil = pos_x + wl;
    let mut ycol = pos_y + wl;
    let mut frag_len = wl;

    let mut xfil_max = pos_x + wl - 1;
    let mut ycol_max = pos_y + wl - 1;
    let mut n_identities = seed_ident;
    let mut max_identities = seed_ident;
    let mut fscore = seed_score;
    let mut fscore_max = fscore;

    let point = scoring.match_score;
    let mismatch = scoring.mismatch_penalty;
    let n_pen = scoring.n_penalty;

    while frag_len < ldiag {
        let vx = seq_x.byte_upper_at(xfil as usize);
        let vy = seq_y.byte_upper_at(ycol as usize);
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

    ForwardResult {
        xfil_max,
        ycol_max,
        fscore_max,
        max_identities,
    }
}

fn extend_backward_scalar(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    fwd: ForwardResult,
    scoring: Scoring,
) -> BackwardResult {
    let ldiag2 = pos_x.min(pos_y);
    let mut xfil_max2 = pos_x;
    let mut ycol_max2 = pos_y;
    let mut fscore = fwd.fscore_max;
    let mut fscore_max = fwd.fscore_max;
    let mut n_identities = fwd.max_identities;
    let mut max_identities = fwd.max_identities;

    let point = scoring.match_score;
    let mismatch = scoring.mismatch_penalty;
    let n_pen = scoring.n_penalty;

    if ldiag2 > 0 {
        let mut xfil2 = pos_x - 1;
        let mut ycol2 = pos_y - 1;
        let mut back_len = 0u64;

        while back_len < ldiag2 {
            let vx = seq_x.byte_upper_at(xfil2 as usize);
            let vy = seq_y.byte_upper_at(ycol2 as usize);
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

    BackwardResult {
        xfil_max2,
        ycol_max2,
        fscore_max,
        max_identities,
    }
}

// ---------------------------------------------------------------------------
// 2-bit helpers (special-run distance queries + even-bit compaction)
// ---------------------------------------------------------------------------

/// Distance forward from `pos` until the next special (`N`/`SEQ_SEP`/IUPAC) position.
/// Returns `0` if `pos` itself lies inside a special run.
#[inline(always)]
fn safe_forward_dist(runs: &[(u64, u64)], pos: u64, cursor: &mut usize) -> u64 {
    while *cursor < runs.len() && runs[*cursor].1 <= pos {
        *cursor += 1;
    }
    if *cursor < runs.len() {
        runs[*cursor].0.saturating_sub(pos)
    } else {
        u64::MAX
    }
}

/// Distance backward from `pos` (inclusive) until the nearest special position at or before `pos`.
/// Returns `0` if `pos` itself lies inside a special run.
#[inline(always)]
fn safe_backward_dist(runs: &[(u64, u64)], pos: u64, cursor: &mut usize) -> u64 {
    while *cursor > 0 && runs[*cursor - 1].0 > pos {
        *cursor -= 1;
    }
    if *cursor > 0 {
        let e = runs[*cursor - 1].1;
        if pos < e { 0 } else { pos - e + 1 }
    } else {
        pos + 1
    }
}

/// Portable SWAR compaction of even bits (`0, 2, 4, ..., 62`) of `x` into bits `0..31`.
#[inline(always)]
fn pack_even_bits_swar(mut x: u64) -> u32 {
    x = (x | (x >> 1)) & 0x3333_3333_3333_3333;
    x = (x | (x >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x >> 4)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x >> 8)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x >> 16)) & 0x0000_0000_FFFF_FFFF;
    x as u32
}

// ---------------------------------------------------------------------------
// Portable 2-bit extension (used on non-AVX2 CPUs for default +4/-4/-1 scoring)
// ---------------------------------------------------------------------------

fn extend_forward_2bit_portable(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    seed_ident: u64,
    seed_score: i64,
) -> ForwardResult {
    let nx = seq_x.len as u64;
    let ny = seq_y.len as u64;
    let mut rem = (nx - pos_x).min(ny - pos_y) - wl;
    let mut xfil = pos_x + wl;
    let mut ycol = pos_y + wl;

    let mut xfil_max = pos_x + wl - 1;
    let mut ycol_max = pos_y + wl - 1;
    let mut n_identities = seed_ident;
    let mut max_identities = seed_ident;
    let mut fscore = seed_score;
    let mut fscore_max = fscore;

    let sp_x = &seq_x.special_runs;
    let sp_y = &seq_y.special_runs;
    let mut cx = sp_x.partition_point(|&(_, e)| e <= xfil);
    let mut cy = sp_y.partition_point(|&(_, e)| e <= ycol);

    while rem > 0 {
        let mut safe = rem
            .min(safe_forward_dist(sp_x, xfil, &mut cx))
            .min(safe_forward_dist(sp_y, ycol, &mut cy));

        if safe == 0 {
            let vx = seq_x.byte_upper_at(xfil as usize);
            let vy = seq_y.byte_upper_at(ycol as usize);
            if vx == SEQ_SEP || vy == SEQ_SEP {
                break;
            }
            if vx == b'N' || vy == b'N' {
                fscore -= 1;
            } else if vx == vy {
                fscore += 4;
                n_identities += 1;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max = xfil;
                    ycol_max = ycol;
                    max_identities = n_identities;
                }
            } else {
                fscore -= 4;
            }
            xfil += 1;
            ycol += 1;
            rem -= 1;
            if fscore < 0 {
                break;
            }
            continue;
        }

        while safe >= 32 {
            let diff = seq_x.load_32(xfil as usize) ^ seq_y.load_32(ycol as usize);
            if diff == 0 {
                fscore += 128;
                n_identities += 32;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max = xfil + 31;
                    ycol_max = ycol + 31;
                    max_identities = n_identities;
                }
                xfil += 32;
                ycol += 32;
                rem -= 32;
                safe -= 32;
                continue;
            }

            let mm_even = (diff | (diff >> 1)) & 0x5555_5555_5555_5555;
            let matches = 32 - mm_even.count_ones() as i64;
            if fscore >= 128 && fscore + 128 <= fscore_max {
                fscore += (2 * matches - 32) * 4;
                n_identities += matches as u64;
                xfil += 32;
                ycol += 32;
                rem -= 32;
                safe -= 32;
                continue;
            }

            // Step this 32-base chunk using `mm_even` directly (zero memory reads).
            let mut bits = mm_even;
            for _ in 0..32 {
                if (bits & 1) == 0 {
                    fscore += 4;
                    n_identities += 1;
                    if fscore > fscore_max {
                        fscore_max = fscore;
                        xfil_max = xfil;
                        ycol_max = ycol;
                        max_identities = n_identities;
                    }
                } else {
                    fscore -= 4;
                }
                xfil += 1;
                ycol += 1;
                rem -= 1;
                safe -= 1;
                bits >>= 2;
                if fscore < 0 {
                    return ForwardResult {
                        xfil_max,
                        ycol_max,
                        fscore_max,
                        max_identities,
                    };
                }
            }
        }

        // Clean tail (< 32 bases before special run or end of diagonal).
        while safe > 0 {
            if seq_x.base_2bit(xfil as usize) == seq_y.base_2bit(ycol as usize) {
                fscore += 4;
                n_identities += 1;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max = xfil;
                    ycol_max = ycol;
                    max_identities = n_identities;
                }
            } else {
                fscore -= 4;
            }
            xfil += 1;
            ycol += 1;
            rem -= 1;
            safe -= 1;
            if fscore < 0 {
                return ForwardResult {
                    xfil_max,
                    ycol_max,
                    fscore_max,
                    max_identities,
                };
            }
        }
    }

    ForwardResult {
        xfil_max,
        ycol_max,
        fscore_max,
        max_identities,
    }
}

fn extend_backward_2bit_portable(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    fwd: ForwardResult,
) -> BackwardResult {
    let mut rem = pos_x.min(pos_y);
    let mut xfil_max2 = pos_x;
    let mut ycol_max2 = pos_y;
    let mut fscore = fwd.fscore_max;
    let mut fscore_max = fwd.fscore_max;
    let mut n_identities = fwd.max_identities;
    let mut max_identities = fwd.max_identities;

    if rem == 0 {
        return BackwardResult {
            xfil_max2,
            ycol_max2,
            fscore_max,
            max_identities,
        };
    }

    let mut xfil2 = pos_x - 1;
    let mut ycol2 = pos_y - 1;

    let sp_x = &seq_x.special_runs;
    let sp_y = &seq_y.special_runs;
    let mut cx = sp_x.partition_point(|&(s, _)| s <= xfil2);
    let mut cy = sp_y.partition_point(|&(s, _)| s <= ycol2);

    while rem > 0 {
        let mut safe = rem
            .min(safe_backward_dist(sp_x, xfil2, &mut cx))
            .min(safe_backward_dist(sp_y, ycol2, &mut cy));

        if safe == 0 {
            let vx = seq_x.byte_upper_at(xfil2 as usize);
            let vy = seq_y.byte_upper_at(ycol2 as usize);
            if vx == SEQ_SEP || vy == SEQ_SEP {
                break;
            }
            if vx == b'N' || vy == b'N' {
                fscore -= 1;
            } else if vx == vy {
                fscore += 4;
                n_identities += 1;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max2 = xfil2;
                    ycol_max2 = ycol2;
                    max_identities = n_identities;
                }
            } else {
                fscore -= 4;
            }
            rem -= 1;
            if fscore < 0 || rem == 0 {
                break;
            }
            xfil2 -= 1;
            ycol2 -= 1;
            continue;
        }

        while safe >= 32 {
            let diff = seq_x.load_32_rev(xfil2 as usize) ^ seq_y.load_32_rev(ycol2 as usize);
            if diff == 0 {
                fscore += 128;
                n_identities += 32;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max2 = xfil2 - 31;
                    ycol_max2 = ycol2 - 31;
                    max_identities = n_identities;
                }
                rem -= 32;
                safe -= 32;
                if rem == 0 {
                    return BackwardResult {
                        xfil_max2,
                        ycol_max2,
                        fscore_max,
                        max_identities,
                    };
                }
                xfil2 -= 32;
                ycol2 -= 32;
                continue;
            }

            let mm_even = (diff | (diff >> 1)) & 0x5555_5555_5555_5555;
            let matches = 32 - mm_even.count_ones() as i64;
            if fscore >= 128 && fscore + 128 <= fscore_max {
                fscore += (2 * matches - 32) * 4;
                n_identities += matches as u64;
                rem -= 32;
                safe -= 32;
                if rem == 0 {
                    return BackwardResult {
                        xfil_max2,
                        ycol_max2,
                        fscore_max,
                        max_identities,
                    };
                }
                xfil2 -= 32;
                ycol2 -= 32;
                continue;
            }

            let mut bits = mm_even;
            for _ in 0..32 {
                if (bits & 1) == 0 {
                    fscore += 4;
                    n_identities += 1;
                    if fscore > fscore_max {
                        fscore_max = fscore;
                        xfil_max2 = xfil2;
                        ycol_max2 = ycol2;
                        max_identities = n_identities;
                    }
                } else {
                    fscore -= 4;
                }
                rem -= 1;
                safe -= 1;
                bits >>= 2;
                if fscore < 0 || rem == 0 {
                    return BackwardResult {
                        xfil_max2,
                        ycol_max2,
                        fscore_max,
                        max_identities,
                    };
                }
                xfil2 -= 1;
                ycol2 -= 1;
            }
        }

        while safe > 0 {
            if seq_x.base_2bit(xfil2 as usize) == seq_y.base_2bit(ycol2 as usize) {
                fscore += 4;
                n_identities += 1;
                if fscore > fscore_max {
                    fscore_max = fscore;
                    xfil_max2 = xfil2;
                    ycol_max2 = ycol2;
                    max_identities = n_identities;
                }
            } else {
                fscore -= 4;
            }
            rem -= 1;
            safe -= 1;
            if fscore < 0 || rem == 0 {
                return BackwardResult {
                    xfil_max2,
                    ycol_max2,
                    fscore_max,
                    max_identities,
                };
            }
            xfil2 -= 1;
            ycol2 -= 1;
        }
    }

    BackwardResult {
        xfil_max2,
        ycol_max2,
        fscore_max,
        max_identities,
    }
}

// ---------------------------------------------------------------------------
// 2-bit + AVX2 + POPCNT (+ optional BMI2) bit-exact bidirectional extension
// ---------------------------------------------------------------------------

/// Compute 32-byte inclusive `i8` prefix sum, horizontal min, and horizontal max for `step`
/// whose lanes are `+1` (match) or `-1` (mismatch).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn prefix_sum_minmax_32(
    step: core::arch::x86_64::__m256i,
) -> (core::arch::x86_64::__m256i, i64, i64) {
    use core::arch::x86_64::*;

    let mut p = step;
    p = _mm256_add_epi8(p, _mm256_slli_si256::<1>(p));
    p = _mm256_add_epi8(p, _mm256_slli_si256::<2>(p));
    p = _mm256_add_epi8(p, _mm256_slli_si256::<4>(p));
    p = _mm256_add_epi8(p, _mm256_slli_si256::<8>(p));

    // Propagate lane 15 (sum of low 16 bytes) across lanes 16..31.
    let lo128 = _mm256_castsi256_si128(p);
    let bcast_lo15 = _mm_shuffle_epi8(lo128, _mm_set1_epi8(15));
    let carry = _mm256_inserti128_si256::<1>(_mm256_setzero_si256(), bcast_lo15);
    p = _mm256_add_epi8(p, carry);

    // Horizontal min and max across all 32 signed bytes.
    let p_lo = _mm256_castsi256_si128(p);
    let p_hi = _mm256_extracti128_si256::<1>(p);
    let mut vmin = _mm_min_epi8(p_lo, p_hi);
    let mut vmax = _mm_max_epi8(p_lo, p_hi);

    vmin = _mm_min_epi8(vmin, _mm_srli_si128::<8>(vmin));
    vmax = _mm_max_epi8(vmax, _mm_srli_si128::<8>(vmax));
    vmin = _mm_min_epi8(vmin, _mm_srli_si128::<4>(vmin));
    vmax = _mm_max_epi8(vmax, _mm_srli_si128::<4>(vmax));
    vmin = _mm_min_epi8(vmin, _mm_srli_si128::<2>(vmin));
    vmax = _mm_max_epi8(vmax, _mm_srli_si128::<2>(vmax));
    vmin = _mm_min_epi8(vmin, _mm_srli_si128::<1>(vmin));
    vmax = _mm_max_epi8(vmax, _mm_srli_si128::<1>(vmax));

    let chunk_min = _mm_cvtsi128_si32(vmin) as i8 as i64;
    let chunk_max = _mm_cvtsi128_si32(vmax) as i8 as i64;
    (p, chunk_min, chunk_max)
}

/// Expand a 32-bit `mm_mask` (bit `k` = 1 if mismatch, 0 if match) into a 32-lane `i8` AVX2
/// step vector (`+1` for match, `-1` for mismatch) using a single 1-cycle `vpor`:
/// `0xFF | 0x01 = 0xFF (-1)` for mismatch, `0x00 | 0x01 = 0x01 (+1)` for match.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn mm_mask_to_step(
    mm_mask: u32,
    byte_shuf: core::arch::x86_64::__m256i,
    bit_mask: core::arch::x86_64::__m256i,
    pos_one: core::arch::x86_64::__m256i,
) -> core::arch::x86_64::__m256i {
    use core::arch::x86_64::*;

    let bcast = _mm256_set1_epi32(mm_mask as i32);
    let bytes = _mm256_shuffle_epi8(bcast, byte_shuf);
    let mm = _mm256_cmpeq_epi8(_mm256_and_si256(bytes, bit_mask), bit_mask);
    _mm256_or_si256(mm, pos_one)
}

/// Emit both `avx2,bmi2,popcnt` (inline `_pext_u64` + hardware `popcnt`) and `avx2,popcnt`
/// (SWAR pack + hardware `popcnt`) bit-exact 2-bit extension kernels with zero `vzeroupper`
/// boundaries inside the hot loop.
#[cfg(target_arch = "x86_64")]
macro_rules! impl_extend_2bit_avx2 {
    ($fwd_fn:ident, $bwd_fn:ident, $features:literal, $pack_expr:expr) => {
        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = $features)]
        #[inline]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $fwd_fn(
            seq_x: &PackedSeq,
            seq_y: &PackedSeq,
            pos_x: u64,
            pos_y: u64,
            wl: u64,
            seed_ident: u64,
            seed_score: i64,
            init_cx: usize,
            init_cy: usize,
        ) -> ForwardResult {
            use core::arch::x86_64::*;

            let nx = seq_x.len as u64;
            let ny = seq_y.len as u64;
            let mut rem = (nx - pos_x).min(ny - pos_y) - wl;
            let mut xfil = pos_x + wl;
            let mut ycol = pos_y + wl;

            let mut xfil_max = pos_x + wl - 1;
            let mut ycol_max = pos_y + wl - 1;
            let mut n_identities = seed_ident;
            let mut max_identities = seed_ident;
            let mut fscore = seed_score;
            let mut fscore_max = fscore;

            let sp_x = &seq_x.special_runs;
            let sp_y = &seq_y.special_runs;
            let mut cx = init_cx;
            let mut cy = init_cy;

            let byte_shuf = _mm256_setr_epi8(
                0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
                3, 3, 3, 3,
            );
            let bit_mask = _mm256_set1_epi64x(0x8040_2010_0804_0201u64 as i64);
            let pos_one = _mm256_set1_epi8(1);

            while rem > 0 {
                let mut safe = rem
                    .min(safe_forward_dist(sp_x, xfil, &mut cx))
                    .min(safe_forward_dist(sp_y, ycol, &mut cy));

                if safe == 0 {
                    let vx = seq_x.byte_upper_at(xfil as usize);
                    let vy = seq_y.byte_upper_at(ycol as usize);
                    if vx == SEQ_SEP || vy == SEQ_SEP {
                        break;
                    }
                    if vx == b'N' || vy == b'N' {
                        fscore -= 1;
                    } else if vx == vy {
                        fscore += 4;
                        n_identities += 1;
                        if fscore > fscore_max {
                            fscore_max = fscore;
                            xfil_max = xfil;
                            ycol_max = ycol;
                            max_identities = n_identities;
                        }
                    } else {
                        fscore -= 4;
                    }
                    xfil += 1;
                    ycol += 1;
                    rem -= 1;
                    if fscore < 0 {
                        break;
                    }
                    continue;
                }

                while safe >= 32 {
                    let diff = seq_x.load_32(xfil as usize) ^ seq_y.load_32(ycol as usize);

                    // Fast-path 1: all 32 bases match (1 scalar instruction).
                    if diff == 0 {
                        fscore += 128;
                        n_identities += 32;
                        if fscore > fscore_max {
                            fscore_max = fscore;
                            xfil_max = xfil + 31;
                            ycol_max = ycol + 31;
                            max_identities = n_identities;
                        }
                        xfil += 32;
                        ycol += 32;
                        rem -= 32;
                        safe -= 32;
                        continue;
                    }

                    let mm_even = (diff | (diff >> 1)) & 0x5555_5555_5555_5555;
                    let matches = 32 - mm_even.count_ones() as i64;

                    // Fast-path 2: long decay tail (hardware popcnt).
                    if fscore >= 128 && fscore + 128 <= fscore_max {
                        fscore += (2 * matches - 32) * 4;
                        n_identities += matches as u64;
                        xfil += 32;
                        ycol += 32;
                        rem -= 32;
                        safe -= 32;
                        continue;
                    }

                    // General path: bit-exact 32-byte i8 AVX2 prefix sum from mismatch bitmask.
                    let mm_mask: u32 = ($pack_expr)(mm_even);
                    let step = unsafe { mm_mask_to_step(mm_mask, byte_shuf, bit_mask, pos_one) };
                    let (p, chunk_min, chunk_max) = unsafe { prefix_sum_minmax_32(step) };

                    if fscore + 4 * chunk_min < 0 {
                        if fscore + 4 * chunk_max > fscore_max {
                            let mut bits = mm_mask;
                            for _ in 0..32 {
                                if (bits & 1) == 0 {
                                    fscore += 4;
                                    n_identities += 1;
                                    if fscore > fscore_max {
                                        fscore_max = fscore;
                                        xfil_max = xfil;
                                        ycol_max = ycol;
                                        max_identities = n_identities;
                                    }
                                } else {
                                    fscore -= 4;
                                    if fscore < 0 {
                                        break;
                                    }
                                }
                                xfil += 1;
                                ycol += 1;
                                bits >>= 1;
                            }
                        }
                        return ForwardResult {
                            xfil_max,
                            ycol_max,
                            fscore_max,
                            max_identities,
                        };
                    }

                    if fscore + 4 * chunk_max > fscore_max {
                        let vmax_bcast = _mm256_set1_epi8(chunk_max as i8);
                        let max_mask =
                            _mm256_movemask_epi8(_mm256_cmpeq_epi8(p, vmax_bcast)) as u32;
                        let best_idx = max_mask.trailing_zeros() as usize;
                        let ident_before = ((best_idx as i64 + 1 + chunk_max) >> 1) as u64;
                        fscore_max = fscore + 4 * chunk_max;
                        xfil_max = xfil + best_idx as u64;
                        ycol_max = ycol + best_idx as u64;
                        max_identities = n_identities + ident_before;
                    }

                    fscore += (2 * matches - 32) * 4;
                    n_identities += matches as u64;
                    xfil += 32;
                    ycol += 32;
                    rem -= 32;
                    safe -= 32;
                }

                // Clean tail (< 32 bases before next special run or diagonal end).
                while safe > 0 {
                    if seq_x.base_2bit(xfil as usize) == seq_y.base_2bit(ycol as usize) {
                        fscore += 4;
                        n_identities += 1;
                        if fscore > fscore_max {
                            fscore_max = fscore;
                            xfil_max = xfil;
                            ycol_max = ycol;
                            max_identities = n_identities;
                        }
                    } else {
                        fscore -= 4;
                        if fscore < 0 {
                            return ForwardResult {
                                xfil_max,
                                ycol_max,
                                fscore_max,
                                max_identities,
                            };
                        }
                    }
                    xfil += 1;
                    ycol += 1;
                    rem -= 1;
                    safe -= 1;
                }
            }

            ForwardResult {
                xfil_max,
                ycol_max,
                fscore_max,
                max_identities,
            }
        }

        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = $features)]
        #[inline]
        unsafe fn $bwd_fn(
            seq_x: &PackedSeq,
            seq_y: &PackedSeq,
            pos_x: u64,
            pos_y: u64,
            fwd: ForwardResult,
            init_cx: usize,
            init_cy: usize,
        ) -> BackwardResult {
            use core::arch::x86_64::*;

            let mut rem = pos_x.min(pos_y);
            let mut xfil_max2 = pos_x;
            let mut ycol_max2 = pos_y;
            let mut fscore = fwd.fscore_max;
            let mut fscore_max = fwd.fscore_max;
            let mut n_identities = fwd.max_identities;
            let mut max_identities = fwd.max_identities;

            if rem == 0 {
                return BackwardResult {
                    xfil_max2,
                    ycol_max2,
                    fscore_max,
                    max_identities,
                };
            }

            let mut xfil2 = pos_x - 1;
            let mut ycol2 = pos_y - 1;

            let sp_x = &seq_x.special_runs;
            let sp_y = &seq_y.special_runs;
            // Derive exact backward `partition_point(|&(s, _)| s <= pos - 1)` from shared `init_c`
            // (`partition_point(|&(_, e)| e <= pos)`) in O(1) with at most a +1 adjustment if
            // `pos` itself lies inside `runs[init_c]`.
            let mut cx = init_cx + usize::from(init_cx < sp_x.len() && sp_x[init_cx].0 <= xfil2);
            let mut cy = init_cy + usize::from(init_cy < sp_y.len() && sp_y[init_cy].0 <= ycol2);

            let byte_shuf = _mm256_setr_epi8(
                0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
                3, 3, 3, 3,
            );
            let bit_mask = _mm256_set1_epi64x(0x8040_2010_0804_0201u64 as i64);
            let pos_one = _mm256_set1_epi8(1);

            while rem > 0 {
                let mut safe = rem
                    .min(safe_backward_dist(sp_x, xfil2, &mut cx))
                    .min(safe_backward_dist(sp_y, ycol2, &mut cy));

                if safe == 0 {
                    let vx = seq_x.byte_upper_at(xfil2 as usize);
                    let vy = seq_y.byte_upper_at(ycol2 as usize);
                    if vx == SEQ_SEP || vy == SEQ_SEP {
                        break;
                    }
                    if vx == b'N' || vy == b'N' {
                        fscore -= 1;
                    } else if vx == vy {
                        fscore += 4;
                        n_identities += 1;
                        if fscore > fscore_max {
                            fscore_max = fscore;
                            xfil_max2 = xfil2;
                            ycol_max2 = ycol2;
                            max_identities = n_identities;
                        }
                    } else {
                        fscore -= 4;
                    }
                    rem -= 1;
                    if fscore < 0 || rem == 0 {
                        break;
                    }
                    xfil2 -= 1;
                    ycol2 -= 1;
                    continue;
                }

                while safe >= 32 {
                    // Load forward-order 32-base window `[xfil2 - 31 ..= xfil2]` and XOR directly;
                    // fast paths 1 and 2 are order-independent, avoiding all 64-bit bit reversals!
                    let diff =
                        seq_x.load_32((xfil2 - 31) as usize) ^ seq_y.load_32((ycol2 - 31) as usize);

                    // Fast-path 1: all 32 bases match.
                    if diff == 0 {
                        fscore += 128;
                        n_identities += 32;
                        if fscore > fscore_max {
                            fscore_max = fscore;
                            xfil_max2 = xfil2 - 31;
                            ycol_max2 = ycol2 - 31;
                            max_identities = n_identities;
                        }
                        rem -= 32;
                        safe -= 32;
                        if rem == 0 {
                            return BackwardResult {
                                xfil_max2,
                                ycol_max2,
                                fscore_max,
                                max_identities,
                            };
                        }
                        xfil2 -= 32;
                        ycol2 -= 32;
                        continue;
                    }

                    let mm_even = (diff | (diff >> 1)) & 0x5555_5555_5555_5555;
                    let matches = 32 - mm_even.count_ones() as i64;

                    // Fast-path 2: long decay tail (order-independent hardware popcnt).
                    if fscore >= 128 && fscore + 128 <= fscore_max {
                        fscore += (2 * matches - 32) * 4;
                        n_identities += matches as u64;
                        rem -= 32;
                        safe -= 32;
                        if rem == 0 {
                            return BackwardResult {
                                xfil_max2,
                                ycol_max2,
                                fscore_max,
                                max_identities,
                            };
                        }
                        xfil2 -= 32;
                        ycol2 -= 32;
                        continue;
                    }

                    // General path: pack 32-bit mismatch mask and reverse bits once so bit `j`
                    // corresponds to `xfil2 - j`.
                    let mm_mask_fwd: u32 = ($pack_expr)(mm_even);
                    let mm_mask = mm_mask_fwd.reverse_bits();
                    let step = unsafe { mm_mask_to_step(mm_mask, byte_shuf, bit_mask, pos_one) };
                    let (p, chunk_min, chunk_max) = unsafe { prefix_sum_minmax_32(step) };

                    if fscore + 4 * chunk_min < 0 {
                        if fscore + 4 * chunk_max > fscore_max {
                            let mut bits = mm_mask;
                            for _ in 0..32 {
                                if (bits & 1) == 0 {
                                    fscore += 4;
                                    n_identities += 1;
                                    if fscore > fscore_max {
                                        fscore_max = fscore;
                                        xfil_max2 = xfil2;
                                        ycol_max2 = ycol2;
                                        max_identities = n_identities;
                                    }
                                } else {
                                    fscore -= 4;
                                }
                                rem -= 1;
                                bits >>= 1;
                                if fscore < 0 || rem == 0 {
                                    break;
                                }
                                xfil2 -= 1;
                                ycol2 -= 1;
                            }
                        }
                        return BackwardResult {
                            xfil_max2,
                            ycol_max2,
                            fscore_max,
                            max_identities,
                        };
                    }

                    if fscore + 4 * chunk_max > fscore_max {
                        let vmax_bcast = _mm256_set1_epi8(chunk_max as i8);
                        let max_mask =
                            _mm256_movemask_epi8(_mm256_cmpeq_epi8(p, vmax_bcast)) as u32;
                        let best_idx = max_mask.trailing_zeros() as usize;
                        let ident_before = ((best_idx as i64 + 1 + chunk_max) >> 1) as u64;
                        fscore_max = fscore + 4 * chunk_max;
                        xfil_max2 = xfil2 - best_idx as u64;
                        ycol_max2 = ycol2 - best_idx as u64;
                        max_identities = n_identities + ident_before;
                    }

                    fscore += (2 * matches - 32) * 4;
                    n_identities += matches as u64;
                    rem -= 32;
                    safe -= 32;
                    if rem == 0 {
                        return BackwardResult {
                            xfil_max2,
                            ycol_max2,
                            fscore_max,
                            max_identities,
                        };
                    }
                    xfil2 -= 32;
                    ycol2 -= 32;
                }

                // Clean tail (< 32 bases).
                while safe > 0 {
                    if seq_x.base_2bit(xfil2 as usize) == seq_y.base_2bit(ycol2 as usize) {
                        fscore += 4;
                        n_identities += 1;
                        if fscore > fscore_max {
                            fscore_max = fscore;
                            xfil_max2 = xfil2;
                            ycol_max2 = ycol2;
                            max_identities = n_identities;
                        }
                    } else {
                        fscore -= 4;
                    }
                    rem -= 1;
                    safe -= 1;
                    if fscore < 0 || rem == 0 {
                        return BackwardResult {
                            xfil_max2,
                            ycol_max2,
                            fscore_max,
                            max_identities,
                        };
                    }
                    xfil2 -= 1;
                    ycol2 -= 1;
                }
            }

            BackwardResult {
                xfil_max2,
                ycol_max2,
                fscore_max,
                max_identities,
            }
        }
    };
}

#[cfg(target_arch = "x86_64")]
impl_extend_2bit_avx2!(
    extend_forward_2bit_avx2_bmi2,
    extend_backward_2bit_avx2_bmi2,
    "avx2,bmi2,popcnt",
    |x: u64| core::arch::x86_64::_pext_u64(x, 0x5555_5555_5555_5555) as u32
);

#[cfg(target_arch = "x86_64")]
impl_extend_2bit_avx2!(
    extend_forward_2bit_avx2_swar,
    extend_backward_2bit_avx2_swar,
    "avx2,popcnt",
    pack_even_bits_swar
);

#[cfg(all(test, target_arch = "x86_64"))]
unsafe fn extend_forward_2bit_avx2(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    use_bmi2: bool,
) -> ForwardResult {
    let cx = sp_cursor_for(&seq_x.special_runs, pos_x);
    let cy = sp_cursor_for(&seq_y.special_runs, pos_y);
    if use_bmi2 {
        unsafe {
            extend_forward_2bit_avx2_bmi2(seq_x, seq_y, pos_x, pos_y, wl, wl, 4 * wl as i64, cx, cy)
        }
    } else {
        unsafe {
            extend_forward_2bit_avx2_swar(seq_x, seq_y, pos_x, pos_y, wl, wl, 4 * wl as i64, cx, cy)
        }
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
unsafe fn extend_backward_2bit_avx2(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    fwd: ForwardResult,
    use_bmi2: bool,
) -> BackwardResult {
    let cx = sp_cursor_for(&seq_x.special_runs, pos_x);
    let cy = sp_cursor_for(&seq_y.special_runs, pos_y);
    if use_bmi2 {
        unsafe { extend_backward_2bit_avx2_bmi2(seq_x, seq_y, pos_x, pos_y, fwd, cx, cy) }
    } else {
        unsafe { extend_backward_2bit_avx2_swar(seq_x, seq_y, pos_x, pos_y, fwd, cx, cy) }
    }
}

/// Compute the exact safe skip boundary along `[pos_x ..= xfil_max]` after a long rejected
/// extension (`xfil_max - pos_x > 256`).
///
/// Partitions `[pos_x ..= xfil_max]` into blocks of size `wl`, computes each block's min/max
/// cumulative score, takes the backward suffix minimum `suf_min`, and advances `safe_b` while
/// `running_max(0..=b) <= suf_min[b + 2]` (for any exact `wl`-mer starting in block `b`,
/// `S(t) > S(p - 1)` holds automatically through block `b + 1`, so its extension is
/// mathematically guaranteed to reproduce the exact same rejected HSP).
///
/// Returns `max(pos_x + safe_b * wl - 1, pos_x + wl / 2)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn exact_rejected_skip(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    xfil_max: u64,
    wl: usize,
    scoring: Scoring,
    scratch: &mut SkipScratch,
) -> u64 {
    let min_skip = pos_x + (wl / 2) as u64;
    if xfil_max <= pos_x {
        return min_skip;
    }

    let span_len = (xfil_max - pos_x + 1) as usize;
    let n_blocks = span_len.div_ceil(wl);
    if n_blocks < 3 {
        return min_skip;
    }

    scratch.blk_min.clear();
    scratch.blk_max.clear();
    scratch.blk_min.reserve(n_blocks);
    scratch.blk_max.reserve(n_blocks);

    let point = scoring.match_score;
    let mismatch = scoring.mismatch_penalty;
    let n_pen = scoring.n_penalty;

    let base_x = pos_x as usize;
    let base_y = pos_y as usize;
    let mut s = 0i64;
    let mut off = 0usize;

    let sp_x = &seq_x.special_runs;
    let sp_y = &seq_y.special_runs;
    let mut cx = sp_x.partition_point(|&(_, e)| e <= pos_x);
    let mut cy = sp_y.partition_point(|&(_, e)| e <= pos_y);

    for _ in 0..n_blocks {
        let end = (off + wl).min(span_len);
        let mut b_max = s;
        let mut b_min = i64::MAX;
        while off < end {
            let px = pos_x + off as u64;
            let py = pos_y + off as u64;
            let safe = safe_forward_dist(sp_x, px, &mut cx)
                .min(safe_forward_dist(sp_y, py, &mut cy))
                .min((end - off) as u64) as usize;
            if safe > 0 {
                for _ in 0..safe {
                    if seq_x.base_2bit(base_x + off) == seq_y.base_2bit(base_y + off) {
                        s += point;
                    } else {
                        s -= mismatch;
                    }
                    if s < b_min {
                        b_min = s;
                    }
                    if s > b_max {
                        b_max = s;
                    }
                    off += 1;
                }
            } else {
                let vx = seq_x.byte_upper_at(base_x + off);
                let vy = seq_y.byte_upper_at(base_y + off);
                if vx == SEQ_SEP || vy == SEQ_SEP {
                    break;
                }
                if vx == b'N' || vy == b'N' {
                    s -= n_pen;
                } else if vx == vy {
                    s += point;
                } else {
                    s -= mismatch;
                }
                if s < b_min {
                    b_min = s;
                }
                if s > b_max {
                    b_max = s;
                }
                off += 1;
            }
        }
        scratch.blk_min.push(b_min);
        scratch.blk_max.push(b_max);
    }

    scratch.suf_min.resize(n_blocks, 0);
    scratch.suf_min[n_blocks - 1] = scratch.blk_min[n_blocks - 1];
    for b in (0..n_blocks - 1).rev() {
        scratch.suf_min[b] = scratch.blk_min[b].min(scratch.suf_min[b + 1]);
    }

    let mut run_max = i64::MIN;
    let mut safe_b = 0usize;
    for b in 0..n_blocks - 2 {
        run_max = run_max.max(scratch.blk_max[b]);
        if run_max <= scratch.suf_min[b + 2] {
            safe_b = b + 1;
        } else {
            break;
        }
    }

    let skip_end = if safe_b > 0 {
        (pos_x + (safe_b * wl) as u64).saturating_sub(1)
    } else {
        0
    };
    skip_end.max(min_skip)
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
                spaced_seed: false,
            },
        );

        assert_eq!(out.fragments.len(), 1);
        let f = &out.fragments[0];
        assert_eq!(f.x_start, 0);
        assert_eq!(f.x_end, (seq.len() - 1) as u64);
        assert_eq!(f.length, seq.len() as u64);
        assert_eq!(f.ident, seq.len() as u64);
        assert!((f.similarity - 100.0).abs() < 1e-4);
        assert_eq!(out.coverage.bins[seq.len()][99], 1);
    }

    #[test]
    fn two_bit_portable_and_avx2_extension_are_bit_exact_with_scalar() {
        let mut rng = 0xDEAD_BEEF_CAFE_BABEu64;
        let mut next_u64 = || {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            rng
        };

        let upper = *b"ACGT";
        let lower = *b"acgt";
        let scoring = Scoring::default();
        let avx2 = has_avx2();
        let bmi2 = has_bmi2();

        for trial in 0..1000 {
            let len = 40 + (next_u64() as usize % 600);
            let mismatch_pct = [0, 5, 15, 30, 48, 60, 80][trial % 7];
            let mut sx = Vec::with_capacity(len);
            let mut sy = Vec::with_capacity(len);

            for _ in 0..len {
                let r = next_u64();
                let b = if (r & 7) == 0 {
                    lower[((r >> 8) & 3) as usize]
                } else {
                    upper[((r >> 8) & 3) as usize]
                };
                sx.push(b);

                let r2 = next_u64();
                let by = if (r2 % 100) < 3 {
                    if (r2 & 1) == 0 { b'N' } else { b'n' }
                } else if (r2 % 1000) < 3 {
                    b'R' // test rare IUPAC ambiguity code
                } else if (r2 % 100) < mismatch_pct {
                    let alt = upper[(((r >> 8) + 1 + (r2 % 3)) & 3) as usize];
                    if (r2 & 4) == 0 {
                        alt.to_ascii_lowercase()
                    } else {
                        alt
                    }
                } else if (r2 & 8) == 0 {
                    b.to_ascii_lowercase()
                } else {
                    b
                };
                sy.push(by);
            }

            // Inject occasional sequence separators to test SEQ_SEP termination in both directions.
            if trial % 5 == 0 && len > 80 {
                let sep_pos = 10 + (next_u64() as usize % (len - 20));
                if (trial & 1) == 0 {
                    sx[sep_pos] = b'*';
                } else {
                    sy[sep_pos] = b'*';
                }
            }

            let wl = [4usize, 8, 12, 16, 24, 32][trial % 6];
            if len <= wl + 4 {
                continue;
            }
            let pos_x = (next_u64() as usize) % (len - wl);
            // Deliberately test mismatched mod-32 alignments between X and Y!
            let shift_y = (next_u64() as usize) % 32;
            let pos_y = if pos_x + shift_y + wl < len {
                pos_x + shift_y
            } else {
                pos_x
            };
            for k in 0..wl {
                let b = upper[(pos_x + k) & 3];
                sx[pos_x + k] = b;
                sy[pos_y + k] = b;
            }

            let gx = Genome::from_records("x", [("x0", sx.as_slice())]);
            let gy = Genome::from_records("y", [("y0", sy.as_slice())]);
            let px = pos_x as u64;
            let py = pos_y as u64;
            let wl_u64 = wl as u64;

            let f_scalar = extend_forward_scalar(
                &gx.seq,
                &gy.seq,
                px,
                py,
                wl_u64,
                wl_u64,
                scoring.match_score * wl_u64 as i64,
                scoring,
            );
            let f_port = extend_forward_2bit_portable(
                &gx.seq,
                &gy.seq,
                px,
                py,
                wl_u64,
                wl_u64,
                4 * wl_u64 as i64,
            );
            assert_eq!(f_port, f_scalar, "portable forward mismatch trial={trial}");

            let b_scalar = extend_backward_scalar(&gx.seq, &gy.seq, px, py, f_scalar, scoring);
            let b_port = extend_backward_2bit_portable(&gx.seq, &gy.seq, px, py, f_port);
            assert_eq!(b_port, b_scalar, "portable backward mismatch trial={trial}");

            #[cfg(target_arch = "x86_64")]
            if avx2 {
                let f_avx2 =
                    unsafe { extend_forward_2bit_avx2(&gx.seq, &gy.seq, px, py, wl_u64, bmi2) };
                assert_eq!(f_avx2, f_scalar, "avx2 forward mismatch trial={trial}");
                let b_avx2 =
                    unsafe { extend_backward_2bit_avx2(&gx.seq, &gy.seq, px, py, f_avx2, bmi2) };
                assert_eq!(b_avx2, b_scalar, "avx2 backward mismatch trial={trial}");
                // Also verify SWAR fallback path directly!
                let f_swar =
                    unsafe { extend_forward_2bit_avx2(&gx.seq, &gy.seq, px, py, wl_u64, false) };
                assert_eq!(f_swar, f_scalar, "avx2+swar forward mismatch trial={trial}");
            }
            let _ = (avx2, bmi2);
        }
    }

    #[test]
    fn exact_rejected_skip_skips_monotonic_span_and_stops_before_dip() {
        let wl = 16usize;
        let g_good = Genome::from_records("g", [("g0", vec![b'A'; 400])]);
        let mut scratch = SkipScratch::default();
        let skip = exact_rejected_skip(
            &g_good.seq,
            &g_good.seq,
            0,
            0,
            399,
            wl,
            Scoring::default(),
            &mut scratch,
        );
        assert!(
            skip >= 350,
            "expected large skip on monotonic span, got {skip}"
        );

        let gx = Genome::from_records("x", [("x0", vec![b'A'; 428])]);
        let mut sy = vec![b'A'; 428];
        sy[64..128].fill(b'C');
        let gy = Genome::from_records("y", [("y0", sy)]);
        let skip2 = exact_rejected_skip(
            &gx.seq,
            &gy.seq,
            0,
            0,
            427,
            wl,
            Scoring::default(),
            &mut scratch,
        );
        assert!(
            skip2 < 128,
            "skip must not step over a valley boundary where a new HSP begins (got {skip2})"
        );
        assert!(skip2 >= (wl / 2) as u64);
    }
}
