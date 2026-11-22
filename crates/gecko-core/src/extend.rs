use rayon::prelude::*;

use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::seq::{Genome, SEQ_SEP};

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
    /// Whether this scoring matches GECKO's default `+4/-4/-1` scheme accelerated by AVX2.
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
        .map(|(s, e)| extend_diagonal_slice(&hits[s..e], genome_x, genome_y, opts, use_avx2))
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
        let outcome =
            extend_single_hit_outcome(h, genome_x, genome_y, opts, use_avx2, &mut out.coverage);
        if let Some(frag) = outcome.fragment {
            last_offset = (h.pos_x + frag.length) as i64;
            out.fragments.push(frag);
        }
    }

    out
}

/// Forward pass state returned by scalar/AVX2 forward extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ForwardResult {
    xfil_max: u64,
    ycol_max: u64,
    fscore_max: i64,
    max_identities: u64,
}

/// Backward pass state returned by scalar/AVX2 backward extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BackwardResult {
    xfil_max2: u64,
    ycol_max2: u64,
    fscore_max: i64,
    max_identities: u64,
}

/// Bidirectional ungapped extension from a single seed hit, faithful to GECKO's `FragFromHit`.
#[inline]
pub(crate) fn extend_single_hit_outcome(
    h: Hit,
    genome_x: &Genome,
    genome_y: &Genome,
    opts: ExtendOptions,
    use_avx2: bool,
    coverage: &mut CoverageMatrix,
) -> ExtendOutcome {
    let data_x = &genome_x.data;
    let data_y = &genome_y.data;
    let nx = data_x.len() as u64;
    let ny = data_y.len() as u64;

    let pos_x = h.pos_x;
    let pos_y = h.pos_y();
    let wl = opts.word_length as u64;

    if pos_x + wl > nx || pos_y + wl > ny {
        return ExtendOutcome {
            fragment: None,
            xfil_max: pos_x + wl.saturating_sub(1),
        };
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

    let (fwd, bwd) = if use_avx2 && opts.scoring.is_default_gecko() {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            let f = extend_forward_avx2(data_x, data_y, pos_x, pos_y, wl);
            let b = extend_backward_avx2(data_x, data_y, pos_x, pos_y, f);
            (f, b)
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let f = extend_forward_scalar(data_x, data_y, pos_x, pos_y, wl, opts.scoring);
            let b = extend_backward_scalar(data_x, data_y, pos_x, pos_y, f, opts.scoring);
            (f, b)
        }
    } else {
        let f = extend_forward_scalar(data_x, data_y, pos_x, pos_y, wl, opts.scoring);
        let b = extend_backward_scalar(data_x, data_y, pos_x, pos_y, f, opts.scoring);
        (f, b)
    };

    // Convert data coordinates (with `SEQ_SEP` separators) to GECKO concatenated nucleotide coordinates.
    let x_start = bwd.xfil_max2 - seq_x as u64;
    let y_start = bwd.ycol_max2 - seq_y as u64;
    let x_end = fwd.xfil_max - seq_x as u64;
    let y_end = fwd.ycol_max - seq_y as u64;
    let length = x_end - x_start + 1;
    let score = bwd.fscore_max.max(0) as u64;
    let max_possible = (length as f64) * (opts.scoring.match_score as f64);
    let similarity = ((score as f64) * 100.0 / max_possible) as f32;

    let out_seq_y = if opts.strand == 'f' {
        seq_y as u64
    } else {
        (genome_y.num_records() - 1 - seq_y) as u64
    };

    coverage.record(length, similarity);

    if length > min_length && similarity > opts.min_similarity {
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
fn ascii_upper_fast(b: u8) -> u8 {
    b & 0xDF
}

fn extend_forward_scalar(
    data_x: &[u8],
    data_y: &[u8],
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    scoring: Scoring,
) -> ForwardResult {
    let nx = data_x.len() as u64;
    let ny = data_y.len() as u64;
    let ldiag = (nx - pos_x).min(ny - pos_y);
    let mut xfil = pos_x + wl;
    let mut ycol = pos_y + wl;
    let mut frag_len = wl;

    let mut xfil_max = pos_x + wl - 1;
    let mut ycol_max = pos_y + wl - 1;
    let mut n_identities = wl;
    let mut max_identities = wl;
    let mut fscore = scoring.match_score * wl as i64;
    let mut fscore_max = fscore;

    let point = scoring.match_score;
    let mismatch = scoring.mismatch_penalty;
    let n_pen = scoring.n_penalty;

    while frag_len < ldiag {
        let vx = ascii_upper_fast(data_x[xfil as usize]);
        let vy = ascii_upper_fast(data_y[ycol as usize]);
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
    data_x: &[u8],
    data_y: &[u8],
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
            let vx = ascii_upper_fast(data_x[xfil2 as usize]);
            let vy = ascii_upper_fast(data_y[ycol2 as usize]);
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

/// Bit-exact AVX2 forward ungapped extension for default `+4/-4/-1(N)` GECKO scoring.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn extend_forward_avx2(
    data_x: &[u8],
    data_y: &[u8],
    pos_x: u64,
    pos_y: u64,
    wl: u64,
) -> ForwardResult {
    use core::arch::x86_64::*;

    let nx = data_x.len() as u64;
    let ny = data_y.len() as u64;
    let mut rem = (nx - pos_x).min(ny - pos_y) - wl;
    let mut xfil = pos_x + wl;
    let mut ycol = pos_y + wl;

    let mut xfil_max = pos_x + wl - 1;
    let mut ycol_max = pos_y + wl - 1;
    let mut n_identities = wl;
    let mut max_identities = wl;
    let mut fscore = 4 * wl as i64;
    let mut fscore_max = fscore;

    let mask_df = _mm256_set1_epi8(0xDFu8 as i8);
    let zero = _mm256_setzero_si256();
    let vn = _mm256_set1_epi8(b'N' as i8);
    let pos_one = _mm256_set1_epi8(1);
    let neg_one = _mm256_set1_epi8(-1);

    let px = data_x.as_ptr();
    let py = data_y.as_ptr();

    while rem >= 32 {
        let vx = unsafe { _mm256_loadu_si256(px.add(xfil as usize) as *const __m256i) };
        let vy = unsafe { _mm256_loadu_si256(py.add(ycol as usize) as *const __m256i) };
        let ux = _mm256_and_si256(vx, mask_df);
        let uy = _mm256_and_si256(vy, mask_df);

        let sp_x = _mm256_or_si256(_mm256_cmpeq_epi8(ux, zero), _mm256_cmpeq_epi8(ux, vn));
        let sp_y = _mm256_or_si256(_mm256_cmpeq_epi8(uy, zero), _mm256_cmpeq_epi8(uy, vn));
        let sp = _mm256_movemask_epi8(_mm256_or_si256(sp_x, sp_y));

        if sp != 0 {
            // Handle this 32-byte chunk scalar (breaks immediately on SEQ_SEP or fscore < 0).
            for _ in 0..32 {
                let cx = ascii_upper_fast(unsafe { *px.add(xfil as usize) });
                let cy = ascii_upper_fast(unsafe { *py.add(ycol as usize) });
                if cx == SEQ_SEP || cy == SEQ_SEP {
                    return ForwardResult {
                        xfil_max,
                        ycol_max,
                        fscore_max,
                        max_identities,
                    };
                }
                if cx == b'N' || cy == b'N' {
                    fscore -= 1;
                } else if cx == cy {
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
                    return ForwardResult {
                        xfil_max,
                        ycol_max,
                        fscore_max,
                        max_identities,
                    };
                }
            }
            continue;
        }

        let eq = _mm256_cmpeq_epi8(ux, uy);
        let eq_mask = _mm256_movemask_epi8(eq) as u32;

        // Fast-path 1: all 32 match.
        if eq_mask == u32::MAX {
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
            continue;
        }

        // Fast-path 2: long decay tail where neither fscore < 0 nor fscore > fscore_max can occur.
        if fscore >= 128 && fscore + 128 <= fscore_max {
            let matches = eq_mask.count_ones() as i64;
            fscore += (2 * matches - 32) * 4;
            n_identities += matches as u64;
            xfil += 32;
            ycol += 32;
            rem -= 32;
            continue;
        }

        // General path: bit-exact 32-byte i8 SIMD prefix sum.
        let step = _mm256_blendv_epi8(neg_one, pos_one, eq);
        let (p, chunk_min, chunk_max) = unsafe { prefix_sum_minmax_32(step) };

        if fscore + 4 * chunk_min < 0 {
            // Terminating chunk: finish exact boundary in scalar tail below.
            break;
        }

        if fscore + 4 * chunk_max > fscore_max {
            let vmax_bcast = _mm256_set1_epi8(chunk_max as i8);
            let max_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(p, vmax_bcast)) as u32;
            let best_idx = max_mask.trailing_zeros() as usize;
            let ident_before = (eq_mask & (u32::MAX >> (31 - best_idx))).count_ones() as u64;
            fscore_max = fscore + 4 * chunk_max;
            xfil_max = xfil + best_idx as u64;
            ycol_max = ycol + best_idx as u64;
            max_identities = n_identities + ident_before;
        }

        let matches = eq_mask.count_ones() as i64;
        fscore += (2 * matches - 32) * 4;
        n_identities += matches as u64;
        xfil += 32;
        ycol += 32;
        rem -= 32;
    }

    // Scalar tail for remaining < 32 bases or the single terminating chunk.
    while rem > 0 {
        let cx = ascii_upper_fast(unsafe { *px.add(xfil as usize) });
        let cy = ascii_upper_fast(unsafe { *py.add(ycol as usize) });
        if cx == SEQ_SEP || cy == SEQ_SEP {
            break;
        }
        if cx == b'N' || cy == b'N' {
            fscore -= 1;
        } else if cx == cy {
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
    }

    ForwardResult {
        xfil_max,
        ycol_max,
        fscore_max,
        max_identities,
    }
}

/// Bit-exact AVX2 backward ungapped extension for default `+4/-4/-1(N)` GECKO scoring.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn extend_backward_avx2(
    data_x: &[u8],
    data_y: &[u8],
    pos_x: u64,
    pos_y: u64,
    fwd: ForwardResult,
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

    let mask_df = _mm256_set1_epi8(0xDFu8 as i8);
    let zero = _mm256_setzero_si256();
    let vn = _mm256_set1_epi8(b'N' as i8);
    let pos_one = _mm256_set1_epi8(1);
    let neg_one = _mm256_set1_epi8(-1);
    let rev_shuf = _mm256_setr_epi8(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6,
        5, 4, 3, 2, 1, 0,
    );

    let px = data_x.as_ptr();
    let py = data_y.as_ptr();

    while rem >= 32 {
        // Load [xfil2 - 31 ..= xfil2] and reverse so lane k corresponds to xfil2 - k.
        let vx_raw = unsafe { _mm256_loadu_si256(px.add((xfil2 - 31) as usize) as *const __m256i) };
        let vy_raw = unsafe { _mm256_loadu_si256(py.add((ycol2 - 31) as usize) as *const __m256i) };
        let vx_s = _mm256_shuffle_epi8(vx_raw, rev_shuf);
        let vy_s = _mm256_shuffle_epi8(vy_raw, rev_shuf);
        let vx = _mm256_permute2x128_si256::<0x01>(vx_s, vx_s);
        let vy = _mm256_permute2x128_si256::<0x01>(vy_s, vy_s);

        let ux = _mm256_and_si256(vx, mask_df);
        let uy = _mm256_and_si256(vy, mask_df);

        let sp_x = _mm256_or_si256(_mm256_cmpeq_epi8(ux, zero), _mm256_cmpeq_epi8(ux, vn));
        let sp_y = _mm256_or_si256(_mm256_cmpeq_epi8(uy, zero), _mm256_cmpeq_epi8(uy, vn));
        let sp = _mm256_movemask_epi8(_mm256_or_si256(sp_x, sp_y));

        if sp != 0 {
            for _ in 0..32 {
                let cx = ascii_upper_fast(unsafe { *px.add(xfil2 as usize) });
                let cy = ascii_upper_fast(unsafe { *py.add(ycol2 as usize) });
                if cx == SEQ_SEP || cy == SEQ_SEP {
                    return BackwardResult {
                        xfil_max2,
                        ycol_max2,
                        fscore_max,
                        max_identities,
                    };
                }
                if cx == b'N' || cy == b'N' {
                    fscore -= 1;
                } else if cx == cy {
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
            continue;
        }

        let eq = _mm256_cmpeq_epi8(ux, uy);
        let eq_mask = _mm256_movemask_epi8(eq) as u32;

        // Fast-path 1: all 32 match.
        if eq_mask == u32::MAX {
            fscore += 128;
            n_identities += 32;
            if fscore > fscore_max {
                fscore_max = fscore;
                xfil_max2 = xfil2 - 31;
                ycol_max2 = ycol2 - 31;
                max_identities = n_identities;
            }
            rem -= 32;
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

        // Fast-path 2: long decay tail.
        if fscore >= 128 && fscore + 128 <= fscore_max {
            let matches = eq_mask.count_ones() as i64;
            fscore += (2 * matches - 32) * 4;
            n_identities += matches as u64;
            rem -= 32;
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

        // General path: bit-exact 32-byte i8 SIMD prefix sum.
        let step = _mm256_blendv_epi8(neg_one, pos_one, eq);
        let (p, chunk_min, chunk_max) = unsafe { prefix_sum_minmax_32(step) };

        if fscore + 4 * chunk_min < 0 {
            break;
        }

        if fscore + 4 * chunk_max > fscore_max {
            let vmax_bcast = _mm256_set1_epi8(chunk_max as i8);
            let max_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(p, vmax_bcast)) as u32;
            let best_idx = max_mask.trailing_zeros() as usize;
            let ident_before = (eq_mask & (u32::MAX >> (31 - best_idx))).count_ones() as u64;
            fscore_max = fscore + 4 * chunk_max;
            xfil_max2 = xfil2 - best_idx as u64;
            ycol_max2 = ycol2 - best_idx as u64;
            max_identities = n_identities + ident_before;
        }

        let matches = eq_mask.count_ones() as i64;
        fscore += (2 * matches - 32) * 4;
        n_identities += matches as u64;
        rem -= 32;
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

    // Scalar tail for remaining < 32 bases or terminating chunk.
    while rem > 0 {
        let cx = ascii_upper_fast(unsafe { *px.add(xfil2 as usize) });
        let cy = ascii_upper_fast(unsafe { *py.add(ycol2 as usize) });
        if cx == SEQ_SEP || cy == SEQ_SEP {
            break;
        }
        if cx == b'N' || cy == b'N' {
            fscore -= 1;
        } else if cx == cy {
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
    }

    BackwardResult {
        xfil_max2,
        ycol_max2,
        fscore_max,
        max_identities,
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
    data_x: &[u8],
    data_y: &[u8],
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

    for _ in 0..n_blocks {
        let end = (off + wl).min(span_len);
        let mut b_max = s;
        let mut b_min = i64::MAX;
        while off < end {
            let vx = ascii_upper_fast(data_x[base_x + off]);
            let vy = ascii_upper_fast(data_y[base_y + off]);
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

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_extension_is_bit_exact_with_scalar() {
        if !has_avx2() {
            return;
        }

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

        for trial in 0..300 {
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
                    sx[sep_pos] = SEQ_SEP;
                } else {
                    sy[sep_pos] = SEQ_SEP;
                }
            }

            let wl = [4usize, 8, 12, 16, 24, 32][trial % 6];
            if len <= wl + 4 {
                continue;
            }
            let pos = (next_u64() as usize) % (len - wl);
            // Make sure the seed itself is an exact match (no SEQ_SEP/N).
            for k in 0..wl {
                let b = upper[(pos + k) & 3];
                sx[pos + k] = b;
                sy[pos + k] = b;
            }

            let pos_u64 = pos as u64;
            let wl_u64 = wl as u64;

            let f_scalar = extend_forward_scalar(&sx, &sy, pos_u64, pos_u64, wl_u64, scoring);
            let f_avx2 = unsafe { extend_forward_avx2(&sx, &sy, pos_u64, pos_u64, wl_u64) };
            assert_eq!(
                f_avx2, f_scalar,
                "forward mismatch at trial={trial} len={len} wl={wl} pos={pos}"
            );

            let b_scalar = extend_backward_scalar(&sx, &sy, pos_u64, pos_u64, f_scalar, scoring);
            let b_avx2 = unsafe { extend_backward_avx2(&sx, &sy, pos_u64, pos_u64, f_avx2) };
            assert_eq!(
                b_avx2, b_scalar,
                "backward mismatch at trial={trial} len={len} wl={wl} pos={pos}"
            );
        }
    }

    #[test]
    fn exact_rejected_skip_skips_monotonic_span_and_stops_before_dip() {
        let wl = 16usize;
        // 1) Monotonically increasing 400 bp match span: safe_b should advance nearly to the end.
        let good = vec![b'A'; 400];
        let mut scratch = SkipScratch::default();
        let skip = exact_rejected_skip(
            &good,
            &good,
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

        // 2) Prefix match (64 bp), deep mismatch valley (64 bp), then long match (300 bp):
        // The deep valley causes suf_min to drop below early blk_max, so skip must stop safely!
        let sx = vec![b'A'; 428];
        let mut sy = vec![b'A'; 428];
        sy[64..128].fill(b'C');
        let skip2 = exact_rejected_skip(&sx, &sy, 0, 0, 427, wl, Scoring::default(), &mut scratch);
        assert!(
            skip2 < 128,
            "skip must not step over a valley boundary where a new HSP begins (got {skip2})"
        );
        // Verify minimum skip bound when safe_b == 0.
        assert!(skip2 >= (wl / 2) as u64);
    }
}
