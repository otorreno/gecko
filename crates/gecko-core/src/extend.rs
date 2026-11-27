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

/// Exact integer reciprocal for `floor((u_max * 100) / length)` on `length in 2..=288`
/// (`t[l] = ceil(2^32 / l)`). Zero-division bit-exact match for [`CoverageMatrix::record`]
/// on short default-scoring (`+4/-4/-1`) extensions up to 128-bp bidirectional fast-reject.
const RECIP_100: [u32; 289] = {
    let mut t = [0u32; 289];
    let mut l = 2usize;
    while l <= 288 {
        t[l] = (1u64 << 32).div_ceil(l as u64) as u32;
        l += 1;
    }
    t
};

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

    /// Zero-division bit-exact coverage recording for short default-scoring (`+4/-4/-1`)
    /// extensions with `2 <= length <= 288` and `0 <= u_max <= length`.
    #[inline(always)]
    pub(crate) fn record_unit_exact(&mut self, length: usize, u_max: u32) {
        debug_assert!((2..=288).contains(&length));
        debug_assert_eq!(self.bins.len(), 1000);
        let recip = unsafe { *RECIP_100.get_unchecked(length) } as u64;
        let col = ((((u_max as u64) * 100) * recip) >> 32).min(99) as usize;
        // Safety: `bins` has length 1000, `length <= 288 < 1000`, `col <= 99`.
        unsafe {
            *self.bins.get_unchecked_mut(length).get_unchecked_mut(col) += 1;
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
    pub exceeded_min_len: bool,
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
/// `seq_index_for_data_pos` lookups until a fragment is accepted when `opts.fixed_length` is set
/// and deferring FP64 similarity division behind `length > min_length` when `coverage` is `None`.
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
    coverage: Option<&mut CoverageMatrix>,
) -> ExtendOutcome {
    // `(fwd.xfil_max - seq_x) - (bwd.xfil_max2 - seq_x) + 1 == fwd.xfil_max - bwd.xfil_max2 + 1`.
    let length = fwd.xfil_max - bwd.xfil_max2 + 1;
    let score = bwd.fscore_max.max(0) as u64;

    let computed_sim = if let Some(cov) = coverage {
        let max_possible = (length as f64) * (opts.scoring.match_score as f64);
        let sim = ((score as f64) * 100.0 / max_possible) as f32;
        cov.record(length, sim);
        Some(sim)
    } else {
        None
    };

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

    let exceeded_min_len = length > min_length;
    if exceeded_min_len {
        let similarity = computed_sim.unwrap_or_else(|| {
            let max_possible = (length as f64) * (opts.scoring.match_score as f64);
            ((score as f64) * 100.0 / max_possible) as f32
        });
        if similarity > opts.min_similarity {
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
            return ExtendOutcome {
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
                exceeded_min_len: true,
            };
        }
    }

    ExtendOutcome {
        fragment: None,
        xfil_max: fwd.xfil_max,
        exceeded_min_len,
    }
}

/// Compute initial seed `(identities, score)` from a precomputed 2-bit XOR `seed_diff`
/// (`seq_x.load_32(pos_x) ^ seq_y.load_32(pos_y)`), avoiding redundant loads on spaced seeds.
#[inline(always)]
pub(crate) fn seed_init_stats_from_diff(
    seed_diff: u64,
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
    let diff = seed_diff & span_mask;
    let mm = ((diff | (diff >> 1)) & 0x5555_5555_5555_5555).count_ones() as u64;
    let ident = wl - mm;
    let score = ident as i64 * scoring.match_score - mm as i64 * scoring.mismatch_penalty;
    (ident, score)
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
    let seed_diff = seq_x.load_32(pos_x as usize) ^ seq_y.load_32(pos_y as usize);
    seed_init_stats_from_diff(seed_diff, wl, scoring, true)
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
            exceeded_min_len: false,
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

    let (fwd, bwd) = if seed_score >= 0 && use_avx2 && opts.scoring.is_default_gecko() {
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
    } else if seed_score >= 0 && opts.scoring.is_default_gecko() {
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
        h,
        pos_x,
        pos_y,
        genome_x,
        genome_y,
        opts,
        fwd,
        bwd,
        Some(coverage),
    )
}

/// Portable single-hit extension with caller-supplied `(seed_ident, seed_score)` and optional
/// coverage matrix.
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_single_hit_with_seed_portable(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    seed_ident: u64,
    seed_score: i64,
    coverage: Option<&mut CoverageMatrix>,
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
            exceeded_min_len: false,
        };
    }

    let (fwd, bwd) = if seed_score >= 0 && opts.scoring.is_default_gecko() {
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

/// Inlineable AVX2 + BMI2 + POPCNT single-hit extension with caller-supplied `(seed_ident,
/// seed_score)`, monotonic `sp_x` cursor, and optional coverage matrix.
///
/// # Safety
/// Caller must run inside a `#[target_feature(enable = "avx2,bmi2,popcnt")]` context and guarantee
/// `sp_x[..cx_x].end <= h.pos_x`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,bmi2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn extend_single_hit_with_seed_avx2_bmi2(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    seed_ident: u64,
    seed_score: i64,
    cx_x: usize,
    coverage: Option<&mut CoverageMatrix>,
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
            exceeded_min_len: false,
        };
    }

    let (fwd, bwd) = if seed_score >= 0 && opts.scoring.is_default_gecko() {
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

/// Inlineable AVX2 + POPCNT (SWAR pack) single-hit extension with caller-supplied `(seed_ident,
/// seed_score)`, monotonic `sp_x` cursor, and optional coverage matrix.
///
/// # Safety
/// Caller must run inside a `#[target_feature(enable = "avx2,popcnt")]` context and guarantee
/// `sp_x[..cx_x].end <= h.pos_x`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn extend_single_hit_with_seed_avx2_swar(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    seed_ident: u64,
    seed_score: i64,
    cx_x: usize,
    coverage: Option<&mut CoverageMatrix>,
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
            exceeded_min_len: false,
        };
    }

    let (fwd, bwd) = if seed_score >= 0 && opts.scoring.is_default_gecko() {
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

/// If `pos` lies inside an `'N'`-run `[s, e)` of `seq`, return the remaining forward length
/// `e - pos` (`>= 1`); otherwise return `0`.
///
/// Precondition: `cx` was just updated by `safe_forward_dist(&seq.special_runs, pos, &mut cx)`.
#[inline(always)]
fn n_run_forward_dist(seq: &PackedSeq, pos: u64, cx: usize) -> u64 {
    let sp = &seq.special_runs;
    if cx >= sp.len() || sp[cx].0 > pos {
        return 0;
    }
    if seq.sep_pos.is_empty() && seq.iupac.is_empty() {
        return sp[cx].1 - pos;
    }
    let nr = &seq.n_runs;
    let idx = nr.partition_point(|&(s, _)| s <= pos);
    if idx > 0 && pos < nr[idx - 1].1 {
        nr[idx - 1].1 - pos
    } else {
        0
    }
}

/// If `pos` lies inside an `'N'`-run `[s, e)` of `seq`, return the remaining backward length
/// `pos - s + 1` (`>= 1`); otherwise return `0`.
///
/// Precondition: `cx` was just updated by `safe_backward_dist(&seq.special_runs, pos, &mut cx)`.
#[inline(always)]
fn n_run_backward_dist(seq: &PackedSeq, pos: u64, cx: usize) -> u64 {
    let sp = &seq.special_runs;
    if cx == 0 || pos >= sp[cx - 1].1 {
        return 0;
    }
    if seq.sep_pos.is_empty() && seq.iupac.is_empty() {
        return pos - sp[cx - 1].0 + 1;
    }
    let nr = &seq.n_runs;
    let idx = nr.partition_point(|&(s, _)| s <= pos);
    if idx > 0 && pos < nr[idx - 1].1 {
        pos - nr[idx - 1].0 + 1
    } else {
        0
    }
}

/// Forward distance from `pos` to the next `SEQ_SEP` (`0` if `pos` is `SEQ_SEP`, `u64::MAX` if none).
#[inline(always)]
fn sep_forward_dist(sep_pos: &[u64], pos: u64) -> u64 {
    if sep_pos.is_empty() {
        return u64::MAX;
    }
    let idx = sep_pos.partition_point(|&p| p < pos);
    if idx < sep_pos.len() {
        sep_pos[idx] - pos
    } else {
        u64::MAX
    }
}

/// Backward distance from `pos` (inclusive) to the nearest `SEQ_SEP` at or before `pos`
/// (`0` if `pos` is `SEQ_SEP`, `u64::MAX` if none).
#[inline(always)]
fn sep_backward_dist(sep_pos: &[u64], pos: u64) -> u64 {
    if sep_pos.is_empty() {
        return u64::MAX;
    }
    let idx = sep_pos.partition_point(|&p| p <= pos);
    if idx > 0 {
        pos - sep_pos[idx - 1]
    } else {
        u64::MAX
    }
}

/// Portable SWAR compaction of even bits (`0, 2, 4, ..., 62`) of `x` into bits `0..31`.
#[inline(always)]
pub(crate) fn pack_even_bits_swar(mut x: u64) -> u32 {
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
            let n_span =
                n_run_forward_dist(seq_x, xfil, cx).max(n_run_forward_dist(seq_y, ycol, cy));
            if n_span > 0 {
                let step = n_span
                    .min(rem)
                    .min(sep_forward_dist(&seq_x.sep_pos, xfil))
                    .min(sep_forward_dist(&seq_y.sep_pos, ycol));
                if step > 0 {
                    if fscore < 0 || step > fscore as u64 {
                        break;
                    }
                    fscore -= step as i64;
                    xfil += step;
                    ycol += step;
                    rem -= step;
                    continue;
                }
            }
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
            let n_span =
                n_run_backward_dist(seq_x, xfil2, cx).max(n_run_backward_dist(seq_y, ycol2, cy));
            if n_span > 0 {
                let step = n_span
                    .min(rem)
                    .min(sep_backward_dist(&seq_x.sep_pos, xfil2))
                    .min(sep_backward_dist(&seq_y.sep_pos, ycol2));
                if step > 0 {
                    if fscore < 0 || step > fscore as u64 {
                        break;
                    }
                    fscore -= step as i64;
                    rem -= step;
                    if rem == 0 {
                        break;
                    }
                    xfil2 -= step;
                    ycol2 -= step;
                    continue;
                }
            }
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
                    let n_span = n_run_forward_dist(seq_x, xfil, cx)
                        .max(n_run_forward_dist(seq_y, ycol, cy));
                    if n_span > 0 {
                        let step = n_span
                            .min(rem)
                            .min(sep_forward_dist(&seq_x.sep_pos, xfil))
                            .min(sep_forward_dist(&seq_y.sep_pos, ycol));
                        if step > 0 {
                            if fscore < 0 || step > fscore as u64 {
                                break;
                            }
                            fscore -= step as i64;
                            xfil += step;
                            ycol += step;
                            rem -= step;
                            continue;
                        }
                    }
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
                    let n_span = n_run_backward_dist(seq_x, xfil2, cx)
                        .max(n_run_backward_dist(seq_y, ycol2, cy));
                    if n_span > 0 {
                        let step = n_span
                            .min(rem)
                            .min(sep_backward_dist(&seq_x.sep_pos, xfil2))
                            .min(sep_backward_dist(&seq_y.sep_pos, ycol2));
                        if step > 0 {
                            if fscore < 0 || step > fscore as u64 {
                                break;
                            }
                            fscore -= step as i64;
                            rem -= step;
                            if rem == 0 {
                                break;
                            }
                            xfil2 -= step;
                            ycol2 -= step;
                            continue;
                        }
                    }
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

/// Compute the exact, 100% lossless skip boundary along `[pos_x ..= xfil_max]` after a rejected
/// extension (contiguous seeds only; returns `pos_x` immediately for spaced seeds or when no
/// full leading block is provably valley-free).
#[allow(clippy::too_many_arguments)]
pub(crate) fn exact_rejected_skip_checked(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    xfil_max: u64,
    wl: usize,
    scoring: Scoring,
    spaced_seed: bool,
    scratch: &mut SkipScratch,
) -> u64 {
    if xfil_max <= pos_x || spaced_seed || wl == 0 || scoring.match_score <= 0 {
        return pos_x;
    }

    // Block-2 lookahead (`suf_min[b + 2]`) holds iff a contiguous `wl`-mer (`+wl * match_score`)
    // cannot drop below `S(p - 1)` within block `b + 1` and cannot zero-tie before block `b + 2`,
    // which requires `wl >= 2` and `max(mismatch_penalty, n_penalty) <= match_score`.
    // Otherwise fall back to unconditionally sound block-1 lookahead (`suf_min[b + 1]`).
    let lookahead = if wl >= 2
        && scoring.mismatch_penalty <= scoring.match_score
        && scoring.n_penalty <= scoring.match_score
    {
        2usize
    } else {
        1usize
    };

    let span_len = (xfil_max - pos_x + 1) as usize;
    let n_blocks = span_len.div_ceil(wl);
    if n_blocks <= lookahead {
        return pos_x;
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
                    return pos_x;
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
    for b in 0..n_blocks - lookahead {
        run_max = run_max.max(scratch.blk_max[b]);
        if run_max <= scratch.suf_min[b + lookahead] {
            safe_b = b + 1;
        } else {
            break;
        }
    }

    if safe_b > 0 {
        (pos_x + (safe_b * wl) as u64).saturating_sub(1)
    } else {
        pos_x
    }
}

/// Contiguous-seed wrapper around [`exact_rejected_skip_checked`] (preserves call compatibility
/// with `minimizer.rs` while enforcing lossless skipping and scoring guards).
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
    exact_rejected_skip_checked(
        seq_x, seq_y, pos_x, pos_y, xfil_max, wl, scoring, false, scratch,
    )
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
        assert_eq!(
            skip2, 0,
            "safe_b == 0 must return pos_x with no +wl/2 clamp (got {skip2})"
        );

        // Spaced seed must always return `pos_x` immediately.
        let skip_spaced = exact_rejected_skip_checked(
            &g_good.seq,
            &g_good.seq,
            0,
            0,
            399,
            wl,
            Scoring::default(),
            true,
            &mut scratch,
        );
        assert_eq!(skip_spaced, 0);

        // High mismatch penalty (`mismatch_penalty > match_score`) uses lookahead-1 safely:
        // on a monotonic all-match span it still skips cleanly (up to n_blocks - 1 blocks),
        // and when block 1 dips below block 0's peak it returns `pos_x`.
        let harsh = Scoring {
            match_score: 4,
            mismatch_penalty: 8,
            n_penalty: 1,
        };
        let skip_harsh_good = exact_rejected_skip_checked(
            &g_good.seq,
            &g_good.seq,
            0,
            0,
            399,
            wl,
            harsh,
            false,
            &mut scratch,
        );
        assert!(skip_harsh_good >= 350);

        let mut sy_dip = vec![b'A'; 400];
        sy_dip[wl..wl + 8].fill(b'C'); // dips block 1 below block 0's peak under harsh -8
        let gy_dip = Genome::from_records("ydip", [("y0", sy_dip)]);
        let skip_harsh_dip = exact_rejected_skip_checked(
            &g_good.seq,
            &gy_dip.seq,
            0,
            0,
            399,
            wl,
            harsh,
            false,
            &mut scratch,
        );
        assert_eq!(skip_harsh_dip, 0);
    }

    #[test]
    fn exact_rejected_skip_recovers_valley_trimmed_hsp() {
        use crate::pipeline::{CompareConfig, SeedingMode, compare_genomes};

        let bases = *b"ACGT";
        let alt_base = |b: u8| match b {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        let mut rng = 0x1234_5678_9ABC_DEF0u64;
        let mut pseudo_dna = |n: usize| -> Vec<u8> {
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                rng = rng
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                out.push(bases[((rng >> 32) & 3) as usize]);
            }
            out
        };

        // Helper to build the two valley-trimming scenarios:
        // - Case 1 (74-bp MEM): Region B = 74 matches, Region C = 34 mismatches
        // - Case 2 (44-bp MEM): Region B = 44 matches, Region C = 33 mismatches
        for &(mem_len, valley_len, expected_x_end, expected_len, expected_score) in &[
            (74usize, 34usize, 1223u64, 1124u64, 3984u64),
            (44usize, 33usize, 1193u64, 1094u64, 3864u64),
        ] {
            let total = 100 + 1050 + mem_len + valley_len + 720 + 1200;
            let sx = pseudo_dna(total);
            let mut sy = sx.clone();

            // Prefix [0..100): all mismatches
            for i in 0..100 {
                sy[i] = alt_base(sx[i]);
            }
            // Region A [100..1150): exact 32-mer at [100..132), then 1 mismatch every 16 bp
            for i in (132..1150).step_by(16) {
                sy[i] = alt_base(sx[i]);
            }
            // Region B [1150..1150+mem_len): exact matches (MEM)
            let reg_c = 1150 + mem_len;
            // Region C [reg_c..reg_c+valley_len): consecutive mismatches
            for i in reg_c..(reg_c + valley_len) {
                sy[i] = alt_base(sx[i]);
            }
            // Region D [reg_d..reg_d+720): 60 x (7 match + 5 mismatch) -> low-sim right tail
            let reg_d = reg_c + valley_len;
            for blk in 0..60 {
                let start = reg_d + blk * 12;
                for m in 7..12 {
                    sy[start + m] = alt_base(sx[start + m]);
                }
            }
            // Region E [reg_d+720..total): all mismatches
            for i in (reg_d + 720)..total {
                sy[i] = alt_base(sx[i]);
            }

            let gx = Genome::from_records("x", [("x0", sx)]);
            let gy = Genome::from_records("y", [("y0", sy)]);
            let cfg = CompareConfig {
                word_length: 32,
                min_length: 1000,
                min_similarity: 75.0,
                freq_limit: 20,
                seeding_mode: SeedingMode::Exact,
                compare_reverse: false,
                soft_mask: false,
                ..CompareConfig::default()
            };
            let res = compare_genomes(&gx, &gy, &cfg, None).unwrap();
            assert_eq!(
                res.fragments.len(),
                1,
                "expected 1 valley-trimmed HSP for {mem_len}-bp MEM, got {:?}",
                res.fragments
            );
            let f = &res.fragments[0];
            assert_eq!(f.x_start, 100);
            assert_eq!(f.x_end, expected_x_end);
            assert_eq!(f.length, expected_len);
            assert_eq!(f.score, expected_score);
        }
    }

    #[test]
    fn bulk_n_run_decay_matches_scalar_all_corner_cases() {
        let scoring = Scoring::default();
        let avx2 = has_avx2();
        let bmi2 = has_bmi2();

        let check = |label: &str,
                     gx: &Genome,
                     gy: &Genome,
                     px: u64,
                     py: u64,
                     wl: u64,
                     seed_ident: u64,
                     seed_score: i64| {
            let f_scalar = extend_forward_scalar(
                &gx.seq, &gy.seq, px, py, wl, seed_ident, seed_score, scoring,
            );
            let f_port =
                extend_forward_2bit_portable(&gx.seq, &gy.seq, px, py, wl, seed_ident, seed_score);
            assert_eq!(f_port, f_scalar, "portable forward mismatch on {label}");

            let b_scalar = extend_backward_scalar(&gx.seq, &gy.seq, px, py, f_scalar, scoring);
            let b_port = extend_backward_2bit_portable(&gx.seq, &gy.seq, px, py, f_port);
            assert_eq!(b_port, b_scalar, "portable backward mismatch on {label}");

            #[cfg(target_arch = "x86_64")]
            if avx2 && has_popcnt() {
                let cx = sp_cursor_for(&gx.seq.special_runs, px);
                let cy = sp_cursor_for(&gy.seq.special_runs, py);
                if bmi2 {
                    let f_bmi2 = unsafe {
                        extend_forward_2bit_avx2_bmi2(
                            &gx.seq, &gy.seq, px, py, wl, seed_ident, seed_score, cx, cy,
                        )
                    };
                    assert_eq!(f_bmi2, f_scalar, "avx2+bmi2 forward mismatch on {label}");
                    let b_bmi2 = unsafe {
                        extend_backward_2bit_avx2_bmi2(&gx.seq, &gy.seq, px, py, f_bmi2, cx, cy)
                    };
                    assert_eq!(b_bmi2, b_scalar, "avx2+bmi2 backward mismatch on {label}");
                }
                let f_swar = unsafe {
                    extend_forward_2bit_avx2_swar(
                        &gx.seq, &gy.seq, px, py, wl, seed_ident, seed_score, cx, cy,
                    )
                };
                assert_eq!(f_swar, f_scalar, "avx2+swar forward mismatch on {label}");
                let b_swar = unsafe {
                    extend_backward_2bit_avx2_swar(&gx.seq, &gy.seq, px, py, f_swar, cx, cy)
                };
                assert_eq!(b_swar, b_scalar, "avx2+swar backward mismatch on {label}");
            }
            let _ = (avx2, bmi2);
        };

        let wl = 16u64;
        let seed_score = 4 * wl as i64; // 64

        // 1. Long N runs (50, 64 exact survival, 65 death, 200, 1000) on X only, Y only, and both.
        for &n_len in &[1usize, 50, 64, 65, 200, 1000] {
            for mode in 0..3 {
                let prefix = 120usize;
                let suffix = 128usize;
                let total = prefix + n_len + wl as usize + n_len + suffix;
                let mut sx = vec![b'A'; total];
                let mut sy = vec![b'A'; total];
                let pos = prefix + n_len;
                // Backward N-run immediately before seed, forward N-run immediately after seed.
                if mode == 0 || mode == 2 {
                    sx[prefix..pos].fill(b'N');
                    sx[pos + wl as usize..pos + wl as usize + n_len].fill(b'N');
                }
                if mode == 1 || mode == 2 {
                    sy[prefix..pos].fill(b'N');
                    sy[pos + wl as usize..pos + wl as usize + n_len].fill(b'N');
                }
                // Offset overlap for mode 2
                if mode == 2 && n_len >= 10 {
                    sy[pos + wl as usize + n_len..pos + wl as usize + n_len + 10].fill(b'N');
                }
                let gx = Genome::from_records("x", [("x0", sx)]);
                let gy = Genome::from_records("y", [("y0", sy)]);
                check(
                    &format!("n_run_len={n_len}_mode={mode}"),
                    &gx,
                    &gy,
                    pos as u64,
                    pos as u64,
                    wl,
                    wl,
                    seed_score,
                );
            }
        }

        // 2. N-run vs SEQ_SEP at offset 0, 17, and end (both X-N/Y-SEP and Y-N/X-SEP).
        for &sep_off in &[0usize, 17, 40] {
            for swap in [false, true] {
                let mut sx = vec![b'A'; 300];
                let mut sy = vec![b'A'; 300];
                let pos = 100usize;
                // Backward N run [60..100), forward N run [116..156)
                if !swap {
                    sx[60..100].fill(b'N');
                    sx[116..156].fill(b'N');
                    sy[116 + sep_off] = b'*';
                    sy[99 - sep_off.min(39)] = b'*';
                } else {
                    sy[60..100].fill(b'N');
                    sy[116..156].fill(b'N');
                    sx[116 + sep_off] = b'*';
                    sx[99 - sep_off.min(39)] = b'*';
                }
                let gx = Genome::from_records("x", [("x0", sx)]);
                let gy = Genome::from_records("y", [("y0", sy)]);
                check(
                    &format!("n_vs_sep_off={sep_off}_swap={swap}"),
                    &gx,
                    &gy,
                    pos as u64,
                    pos as u64,
                    wl,
                    wl,
                    seed_score,
                );
            }
        }

        // 3. N-run aligned with IUPAC 'R'/'Y' on the other sequence.
        {
            let mut sx = vec![b'A'; 260];
            let mut sy = vec![b'A'; 260];
            sx[60..100].fill(b'N');
            sx[116..156].fill(b'N');
            sy[75] = b'R';
            sy[120] = b'Y';
            sy[130] = b'R';
            let gx = Genome::from_records("x", [("x0", sx)]);
            let gy = Genome::from_records("y", [("y0", sy)]);
            check("n_vs_iupac", &gx, &gy, 100, 100, wl, wl, seed_score);
        }

        // 4. Backward N-run reaching index 0 (`rem == step`) -> verifies zero u64 underflow.
        for &(px, py) in &[(20u64, 20u64), (20u64, 50u64), (50u64, 20u64)] {
            let mut sx = vec![b'A'; 200];
            let mut sy = vec![b'A'; 200];
            sx[0..px as usize].fill(b'N');
            sy[0..py as usize].fill(b'N');
            let gx = Genome::from_records("x", [("x0", sx)]);
            let gy = Genome::from_records("y", [("y0", sy)]);
            check(
                &format!("bwd_rem_eq_step_px={px}_py={py}"),
                &gx,
                &gy,
                px,
                py,
                wl,
                wl,
                seed_score,
            );
        }

        // 5. Reverse-complemented Genome with N-runs + SEQ_SEP + IUPAC.
        {
            let mut sx = vec![b'G'; 300];
            let mut sy = vec![b'C'; 300];
            sx[80..110].fill(b'N');
            sx[150..200].fill(b'N');
            sy[90..130].fill(b'N');
            sy[220] = b'*';
            sy[140] = b'R';
            let gx = Genome::from_records("x", [("x0", sx)]);
            let gy_rc = Genome::from_records("y", [("y0", sy)]).reverse_complement();
            // Place a matching seed at pos 120 in rc coordinates.
            let mut rx = vec![b'G'; 300];
            for k in 0..wl as usize {
                rx[120 + k] = gy_rc.seq.byte_upper_at(120 + k);
            }
            rx[80..110].fill(b'N');
            rx[150..190].fill(b'N');
            let gx2 = Genome::from_records("x2", [("x0", rx)]);
            check(
                "rev_comp_n_sep_iupac",
                &gx2,
                &gy_rc,
                120,
                120,
                wl,
                wl,
                seed_score,
            );
            let _ = gx;
        }

        // 6. Negative initial seed_score = -4 immediately hitting an N-run.
        {
            let mut sx = vec![b'A'; 128];
            let mut sy = vec![b'A'; 128];
            sx[40 + wl as usize..40 + wl as usize + 30].fill(b'N');
            sy[20..40].fill(b'N');
            let gx = Genome::from_records("x", [("x0", sx)]);
            let gy = Genome::from_records("y", [("y0", sy)]);
            check("neg_seed_score_n_run", &gx, &gy, 40, 40, wl, 7, -4);
        }
    }

    #[test]
    fn record_unit_exact_matches_fp64_record_exhaustively() {
        let mut cov_exact = CoverageMatrix::new();
        let mut cov_fp = CoverageMatrix::new();
        for length in 2usize..=288 {
            for u_max in 0u32..=(length as u32) {
                cov_exact.record_unit_exact(length, u_max);
                let sim = (((4 * u_max) as f64) * 100.0 / ((length as f64) * 4.0)) as f32;
                cov_fp.record(length as u64, sim);
            }
        }
        assert_eq!(cov_exact.bins, cov_fp.bins);
    }

    #[test]
    fn negative_seed_score_spaced_seed_matches_scalar_oracle() {
        use crate::pipeline::{CompareConfig, compare_genomes};
        use crate::spaced::SpacedSeed;

        // Custom palindromic spaced seed `10000111100001` (weight=6, span=14):
        // when all 8 don't-care positions mismatch, `seed_score = 6*4 - 8*4 = -8 < 0`.
        // Previously, 2-bit/AVX2 Fast-Path 1 (`diff == 0`) omitted `if fscore < 0 { break; }`,
        // causing default `(4, 4, 1)` to extend past base 0 and diverge from 2x-scaled `(8, 8, 2)`.
        let seed = SpacedSeed::parse("10000111100001").unwrap();
        let wl = seed.span as u64;
        let alt = |b: u8| match b {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };

        let alphabet = *b"ACGT";
        let mut rng = 0x5EED_CAFE_0001_0002u64;
        let sx: Vec<u8> = (0..320)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();
        let mut sy: Vec<u8> = sx.iter().map(|&b| alt(b)).collect();

        // Plant a hit at pos=100 matching ONLY the 6 care positions of `10000111100001`
        // (`seed_score = -8`), followed immediately by a 46-bp exact match `[114..160)`
        // (and 40-bp exact match before pos=100 in a second test case to exercise backward!).
        for (i, ch) in "10000111100001".bytes().enumerate() {
            if ch == b'1' {
                sy[100 + i] = sx[100 + i];
            }
        }
        sy[114..160].copy_from_slice(&sx[114..160]);
        let gx = Genome::from_records("x", [("x0", sx.as_slice())]);
        let gy = Genome::from_records("y", [("y0", sy.as_slice())]);

        // 1. Direct entry-point parity: portable, AVX2+BMI2, AVX2+SWAR, and extend_single_hit_outcome
        //    vs exact scalar oracle across `seed_score ∈ {-8, -4}` (forward & backward clean matches).
        for &(seed_ident, seed_score, label) in &[
            (6u64, -8i64, "neg8_clean_fwd"),
            (7u64, -4i64, "neg4_clean_fwd"),
        ] {
            let opts = ExtendOptions {
                min_length: 40,
                min_similarity: 50.0,
                word_length: seed.span,
                fixed_length: true,
                strand: 'f',
                scoring: Scoring::default(),
                spaced_seed: true,
            };
            let hit = Hit {
                diag: 0,
                pos_x: 100,
            };
            let f_oracle = extend_forward_scalar(
                &gx.seq,
                &gy.seq,
                100,
                100,
                wl,
                seed_ident,
                seed_score,
                opts.scoring,
            );
            let b_oracle =
                extend_backward_scalar(&gx.seq, &gy.seq, 100, 100, f_oracle, opts.scoring);
            let mut cov_oracle = CoverageMatrix::new();
            let out_oracle = finalize_extend_outcome(
                hit,
                100,
                100,
                &gx,
                &gy,
                opts,
                f_oracle,
                b_oracle,
                Some(&mut cov_oracle),
            );

            let mut cov_port = CoverageMatrix::new();
            let out_port = extend_single_hit_with_seed_portable(
                hit,
                &gx,
                &gy,
                opts,
                seed_ident,
                seed_score,
                Some(&mut cov_port),
            );
            assert_eq!(
                out_port.fragment, out_oracle.fragment,
                "portable frag {label}"
            );
            assert_eq!(
                out_port.xfil_max, out_oracle.xfil_max,
                "portable xfil_max {label}"
            );
            assert_eq!(cov_port.bins, cov_oracle.bins, "portable cov {label}");

            if seed_score == -8 {
                // With seed_score = -8 and base 0 matching, fscore = -4 < 0 -> MUST stop at 100 + wl.
                assert_eq!(out_oracle.xfil_max, 100 + wl);
                assert!(out_oracle.fragment.is_none());

                // Also verify extend_single_hit_outcome (computes seed_score = -8 from seq_x/seq_y).
                let mut cov_entry = CoverageMatrix::new();
                let out_entry = extend_single_hit_outcome(
                    hit,
                    &gx,
                    &gy,
                    opts,
                    has_avx2(),
                    has_bmi2(),
                    &mut cov_entry,
                );
                assert_eq!(out_entry.fragment, out_oracle.fragment);
                assert_eq!(out_entry.xfil_max, out_oracle.xfil_max);
                assert_eq!(cov_entry.bins, cov_oracle.bins);
            }

            #[cfg(target_arch = "x86_64")]
            if has_avx2() && has_popcnt() {
                let cx = sp_cursor_for(&gx.seq.special_runs, 100);
                if has_bmi2() {
                    let mut cov_bmi2 = CoverageMatrix::new();
                    let out_bmi2 = unsafe {
                        extend_single_hit_with_seed_avx2_bmi2(
                            hit,
                            &gx,
                            &gy,
                            opts,
                            seed_ident,
                            seed_score,
                            cx,
                            Some(&mut cov_bmi2),
                        )
                    };
                    assert_eq!(out_bmi2.fragment, out_oracle.fragment, "bmi2 frag {label}");
                    assert_eq!(
                        out_bmi2.xfil_max, out_oracle.xfil_max,
                        "bmi2 xfil_max {label}"
                    );
                    assert_eq!(cov_bmi2.bins, cov_oracle.bins, "bmi2 cov {label}");
                }
                let mut cov_swar = CoverageMatrix::new();
                let out_swar = unsafe {
                    extend_single_hit_with_seed_avx2_swar(
                        hit,
                        &gx,
                        &gy,
                        opts,
                        seed_ident,
                        seed_score,
                        cx,
                        Some(&mut cov_swar),
                    )
                };
                assert_eq!(out_swar.fragment, out_oracle.fragment, "swar frag {label}");
                assert_eq!(
                    out_swar.xfil_max, out_oracle.xfil_max,
                    "swar xfil_max {label}"
                );
                assert_eq!(cov_swar.bins, cov_oracle.bins, "swar cov {label}");
            }
        }

        // 2. Full compare_genomes differential parity: default (4, 4, 1) vs 2x-scaled (8, 8, 2)
        //    across two isolated contigs:
        //    - contig 0: 40-bp exact prefix `[60..100)` + `seed_score=-8` at 100 + `[114..160)`
        //      (exercises backward extension starting from `fwd.fscore_max = -4 < 0`).
        //    - contig 1: exact minimal repro — `seed_score=-8` at 100 + 55-bp HSP `[105..160)`
        //      with mismatches `{109..112}` (must emit `[105..=159]` len=55, NOT `[100..=159]` len=60).
        let mut sy0 = sy.clone();
        sy0[60..100].copy_from_slice(&sx[60..100]);

        let mut sy1: Vec<u8> = sx.iter().map(|&b| alt(b)).collect();
        sy1[105..160].copy_from_slice(&sx[105..160]);
        sy1[100] = sx[100];
        for pos in 109..113 {
            sy1[pos] = alt(sx[pos]);
        }

        let gx_full =
            Genome::from_records("x_full", [("x0", sx.as_slice()), ("x1", sx.as_slice())]);
        let gy_full =
            Genome::from_records("y_full", [("y0", sy0.as_slice()), ("y1", sy1.as_slice())]);

        let cfg_1x = CompareConfig {
            word_length: seed.span,
            spaced_seed: Some(seed.clone()),
            min_length: 50,
            min_similarity: 50.0,
            freq_limit: 50,
            compare_reverse: true,
            scoring: Scoring {
                match_score: 4,
                mismatch_penalty: 4,
                n_penalty: 1,
            },
            collect_coverage: true,
            ..Default::default()
        };
        let cfg_2x = CompareConfig {
            scoring: Scoring {
                match_score: 8,
                mismatch_penalty: 8,
                n_penalty: 2,
            },
            ..cfg_1x.clone()
        };

        let r1 = compare_genomes(&gx_full, &gy_full, &cfg_1x, None).unwrap();
        let r2 = compare_genomes(&gx_full, &gy_full, &cfg_2x, None).unwrap();

        let c1_off = sx.len() as u64;
        assert!(
            r1.fragments.iter().any(|f| {
                f.seq_x == 1
                    && f.seq_y == 1
                    && f.x_start == c1_off + 105
                    && f.x_end == c1_off + 159
                    && f.length == 55
                    && f.ident == 51
                    && f.score == 188
            }),
            "expected exact 55-bp HSP [105..=159] on contig 1, got {:?}",
            r1.fragments
        );
        assert!(
            !r1.fragments.iter().any(|f| {
                f.seq_x == 1 && f.seq_y == 1 && f.x_start == c1_off + 100 && f.length == 60
            }),
            "must not emit buggy 60-bp HSP [100..=159] from negative seed_score, got {:?}",
            r1.fragments
        );
        assert_eq!(r1.fragments.len(), r2.fragments.len());
        for (f1, f2) in r1.fragments.iter().zip(r2.fragments.iter()) {
            assert_eq!(f1.diag, f2.diag);
            assert_eq!(f1.x_start, f2.x_start);
            assert_eq!(f1.y_start, f2.y_start);
            assert_eq!(f1.x_end, f2.x_end);
            assert_eq!(f1.y_end, f2.y_end);
            assert_eq!(f1.length, f2.length);
            assert_eq!(f1.ident, f2.ident);
            assert_eq!(f1.strand, f2.strand);
            assert_eq!(f1.score * 2, f2.score);
        }
        assert_eq!(
            r1.coverage.bins, r2.coverage.bins,
            "CoverageMatrix (.frags.MAT) must match 2x-scaled scalar oracle bit-for-bit"
        );
    }
}
