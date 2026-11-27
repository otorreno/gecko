//! Shift-invariant canonical minimizer super-seeding engine (`F9`, Wave 3 overhaul) with
//! 100% lossless exact-$k$ verification, high-frequency exact-$k$ overflow routing (`W3-F6`),
//! branchless AVX2/portable sliding-window extraction (`W3-F7`), flat 2-level split-tag radix
//! indexing + $O(w)$ rolling HF window checks + balanced multi-slab (`W3-F8`), and scalar
//! short-death fast-reject + `DiagCache32` + `ActiveDiagRing` + lossless valley-trimming (`W3-F9`).
//!
//! For contiguous seeds $k \in [18, 32]$, every exact $k$-mer window of length $k$ contains
//! $w = k - k' + 1$ sub-$k'$-mers (e.g. $k' = 21, w = 12$ for $k = 32$). Indexing canonical
//! minimizer positions of the reference (leftmost + rightmost ties for strict strand symmetry)
//! and probing the leftmost canonical minimizer per query window shrinks index size, build work,
//! and query probes by ~3x–6.5x while guaranteeing 100% bit-exact `.frags` parity with stride-1
//! `RefIndex` — including finite per-strand `freq_limit`, palindromes, soft-masking, multi-FASTA,
//! `--forward-only`, and hash-partitioned multi-slab memory budgets.
//!
//! Crucial hashing detail: the minimizer ordering score `minimizer_order_hash(canon)` is skewed
//! toward `0` (minimum of `w` hashes). Therefore, all partition/bucket/Bloom/tag indexing uses an
//! independent bijective secondary hash `idx_hash = index_hash(canon) = fmix64(canon ^ IDX_HASH_SALT)`
//! so all 256 partitions, CSR buckets, and L3 Blocked Bloom blocks remain 100% uniformly loaded.

use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use rayon::prelude::*;

use crate::error::{GeckoError, Result};
use crate::extend::{
    CoverageMatrix, ExtendOptions, SkipScratch, exact_rejected_skip,
    extend_single_hit_with_seed_portable, has_avx2, has_bmi2, has_popcnt,
};
#[cfg(target_arch = "x86_64")]
use crate::extend::{extend_single_hit_with_seed_avx2_bmi2, extend_single_hit_with_seed_avx2_swar};
use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::pipeline::{CompareConfig, CompareResult, CompareStats, profile_enabled};
use crate::ref_index::{
    BlockedBloom, BloomBlock, PosStrandVec, choose_bloom_params, count_valid_kmers, fmix64,
    insert_in_bloom_block, prefetch_t1, prefetch_t2, probe_bucket_binary_u32,
    probe_sorted_tags_scalar,
};
use crate::seq::{Genome, PackedSeq};
#[cfg(target_arch = "x86_64")]
use crate::stream::pack_mm32_bmi2;
use crate::stream::{
    ActiveDiagRing, COMPACT_DIAG_MAX_LEN, DiagPair32, DiagPair64, DiagPairOps, FastRejectOutcome,
    StrandDiagCache, dedup_fragments, pack_mm32_swar, try_fast_short_reject_auto,
};

/// Minimum contiguous seed length eligible for canonical minimizer super-seeding.
pub const MINIMIZER_MIN_K: usize = 18;

const NUM_PARTITIONS: usize = 256;
const MIN_QUERY_CHUNK_BYTES: usize = 1 << 20; // ~1 Mb dynamic Rayon chunks
const MIN_BATCH_SIZE: usize = 64;

/// High-frequency minimizer bucket routing cap (`W3-F6`).
///
/// Any canonical minimizer occurring more than `MINIMIZER_OVERFLOW_CAP` times in the reference
/// is replaced by a single sentinel in the normal minimizer CSR, and its surviving (non-`freq_limit`
/// dropped) exact `k`-mers are stored in a tiny exact-`k` `OverflowKmerCsr` with its own cache-line
/// Blocked Bloom filter. This eliminates quadratic repeat cross-products on unmasked DNA while
/// keeping exact-`k` overflow memory bounded to a tiny fraction of the reference.
const MINIMIZER_OVERFLOW_CAP: usize = 4;

/// Reserved `pos_strand` sentinel indicating that this minimizer was routed to `OverflowKmerCsr`.
/// Valid `(pos << 1) | strand` entries satisfy `pos <= seq.len - 1 < (1 << 31) - 2` when
/// `use_u32_ps` is true, so `u32::MAX` / `u64::MAX` can never collide with a real coordinate.
const OVERFLOW_SENTINEL_U32: u32 = u32::MAX;
const OVERFLOW_SENTINEL_U64: u64 = u64::MAX;
const OVERFLOW_SENTINEL_POS: u64 = u64::MAX;

/// Independent salt mixed into `canon` before `fmix64` to decorrelate partition/bucket/Bloom
/// hashing from the minimum-biased minimizer selection score `minimizer_order_hash(canon)`.
const IDX_HASH_SALT: u64 = 0x9E37_79B9_7F4A_7C15;

/// Odd golden-ratio 32-bit salt so poly-A (`canon == 0`) does NOT map to score `0`.
const ORDER_SALT: u32 = 0x9E37_79B9;
/// MurmurHash3 32-bit finalizer odd multiplier (strong bit avalanche in 1x 32x32->64 mul).
const ORDER_MUL: u64 = 0x85EB_CA6B;

// Block-SIMD sliding-window constants (`W3-F7`).
const POS_BITS: u32 = 9;
const POS_MASK: u32 = (1u32 << POS_BITS) - 1; // 0x1FF (511)
#[cfg(target_arch = "x86_64")]
const SCORE_MASK: u32 = !POS_MASK; // 0xFFFF_FE00 (top 23 bits)
const KP_CANON_MASK: u64 = (1u64 << 42) - 1; // 2 * max_kp (42 bits)
const BLOCK_WINS: usize = 256; // exact 4 x u64 words in emit_bm[0..4]
const MAX_W: usize = 16; // strictly > max window (12)
const BUF_LEN: usize = 288; // >= BLOCK_WINS + MAX_W + 8 (8-lane aligned pad)

// Flat 2-level radix Pass-1 staging constants (`W3-F8a`).
const STAGE_CAP: usize = 64;
const STAGE_MID_BYTES: usize = STAGE_CAP * 3 + 1;

// Peak memory estimation constants (`W3-F8c`).
const MIN_PEAK_BYTES_PER_ENTRY_U32: usize = 18;
const MIN_PEAK_BYTES_PER_ENTRY_U64: usize = 22;
const MIN_FIXED_OVERHEAD: usize = 64 * 1024;

/// Return `(kp, window)` for contiguous `word_length` in `18..=32`, or `None` if ineligible.
///
/// Guarantees `kp + window - 1 == word_length`, `kp in 14..=21`, and `window in 5..=12`.
#[inline]
pub fn minimizer_params(word_length: usize) -> Option<(usize, usize)> {
    if !(MINIMIZER_MIN_K..=32).contains(&word_length) {
        return None;
    }
    let kp = word_length.div_ceil(2).saturating_add(5).min(21);
    let window = word_length - kp + 1;
    debug_assert!((1..=MAX_W).contains(&window));
    Some((kp, window))
}

/// Fast 1-multiply 23-bit canonical `kp`-mer ordering score (`0 .. (1 << 23)`).
///
/// - Folds all `2 * kp <= 42` bits of LSB-first `canon` into 32 bits via odd shift `>> 21`
///   (avoids period-2/4/6/8 cancellation on DNA repeats).
/// - XORs `ORDER_SALT` before multiplying so poly-A (`canon == 0`) hashes to `~0.577`, NOT `0`.
/// - Single 32x32->64 multiply (`imul` scalar / `vpmuludq` 4-wide AVX2) + high/low 32-bit XOR fold.
#[inline(always)]
fn minimizer_order_hash(canon: u64) -> u32 {
    let folded = ((canon as u32) ^ ((canon >> 21) as u32) ^ ORDER_SALT) as u64;
    let prod = folded.wrapping_mul(ORDER_MUL);
    ((prod ^ (prod >> 32)) as u32) >> 9
}

/// Uniform bijective index hash decorrelated from the minimum-biased minimizer score.
#[inline(always)]
fn index_hash(canon: u64) -> u64 {
    fmix64(canon ^ IDX_HASH_SALT)
}

/// Bijective bit-permutation of `idx_hash` placing `[bucket_bits (high) | tag32 (mid) | rest (low)]`
/// so a single `sort_unstable` groups identical `canon` runs contiguously AND sorts buckets + 32-bit
/// tags simultaneously with zero information loss.
#[inline(always)]
fn idx_hash_to_sort_key(idx_hash: u64, hash_shift: u32) -> u64 {
    let mid_bits = hash_shift - 32;
    let top = idx_hash & (u64::MAX << hash_shift);
    let tag = idx_hash & 0xFFFF_FFFF;
    let mid = (idx_hash >> 32) & ((1u64 << mid_bits) - 1);
    top | (tag << mid_bits) | mid
}

/// Inverse of [`idx_hash_to_sort_key`] recovering the original `idx_hash`.
#[inline(always)]
fn sort_key_to_idx_hash(sort_key: u64, hash_shift: u32) -> u64 {
    let mid_bits = hash_shift - 32;
    let top = sort_key & (u64::MAX << hash_shift);
    let tag = (sort_key >> mid_bits) & 0xFFFF_FFFF;
    let mid = sort_key & ((1u64 << mid_bits) - 1);
    top | (mid << 32) | tag
}

// ---------------------------------------------------------------------------
// Branchless SIMD / portable sliding-window canonical minimizer extraction (W3-F7)
// ---------------------------------------------------------------------------

struct MinScratch {
    canon_kind: [u64; BUF_LEN], // bits 0..42 = LSB-first canon, bits 62..63 = kind (0/1/2)
    packed_l: [u32; BUF_LEN],   // (score23 << 9) | local_idx
    packed_r: [u32; BUF_LEN],   // packed_l ^ POS_MASK (used iff NEED_RIGHT)
    s2: [u32; BUF_LEN],
    s4: [u32; BUF_LEN],
    out_l: [u32; BUF_LEN + 1], // out_l[0] = u32::MAX sentinel; outputs at out_l[1..1+n_wins]
    out_r: [u32; BUF_LEN],     // used iff NEED_RIGHT
}

impl MinScratch {
    #[inline]
    fn new() -> Self {
        Self {
            canon_kind: [0u64; BUF_LEN],
            packed_l: [0u32; BUF_LEN],
            packed_r: [0u32; BUF_LEN],
            s2: [0u32; BUF_LEN],
            s4: [0u32; BUF_LEN],
            out_l: [0u32; BUF_LEN + 1],
            out_r: [0u32; BUF_LEN],
        }
    }
}

/// Branchless binary-step sliding-window minimum over `src[0..n_wins + w - 1]` (portable).
#[inline]
fn sliding_min_u32_portable(
    src: &[u32; BUF_LEN],
    s2: &mut [u32; BUF_LEN],
    s4: &mut [u32; BUF_LEN],
    out: &mut [u32],
    n_wins: usize,
    w: usize,
) {
    let n_kps = n_wins + w - 1;
    if w == 1 {
        out[..n_wins].copy_from_slice(&src[..n_wins]);
        return;
    }
    if w == 2 {
        for i in 0..n_wins {
            unsafe {
                *out.get_unchecked_mut(i) = (*src.get_unchecked(i)).min(*src.get_unchecked(i + 1));
            }
        }
        return;
    }
    // Pass 1: s2[i] = min(src[i..i+2])
    let len_s2 = n_kps - 1;
    for i in 0..len_s2 {
        unsafe {
            *s2.get_unchecked_mut(i) = (*src.get_unchecked(i)).min(*src.get_unchecked(i + 1));
        }
    }
    if w <= 4 {
        let d = w - 2;
        for i in 0..n_wins {
            unsafe {
                *out.get_unchecked_mut(i) = (*s2.get_unchecked(i)).min(*s2.get_unchecked(i + d));
            }
        }
        return;
    }
    // Pass 2: s4[i] = min(s2[i], s2[i+2]) = min(src[i..i+4])
    let len_s4 = n_kps - 3;
    for i in 0..len_s4 {
        unsafe {
            *s4.get_unchecked_mut(i) = (*s2.get_unchecked(i)).min(*s2.get_unchecked(i + 2));
        }
    }
    if w <= 8 {
        let d = w - 4;
        for i in 0..n_wins {
            unsafe {
                *out.get_unchecked_mut(i) = (*s4.get_unchecked(i)).min(*s4.get_unchecked(i + d));
            }
        }
    } else if w <= 12 {
        let d = w - 4;
        for i in 0..n_wins {
            unsafe {
                let ab = (*s4.get_unchecked(i)).min(*s4.get_unchecked(i + 4));
                *out.get_unchecked_mut(i) = ab.min(*s4.get_unchecked(i + d));
            }
        }
    } else {
        let len_s8 = n_kps - 7;
        for i in 0..len_s8 {
            unsafe {
                *s2.get_unchecked_mut(i) = (*s4.get_unchecked(i)).min(*s4.get_unchecked(i + 4));
            }
        }
        let d = w - 8;
        for i in 0..n_wins {
            unsafe {
                *out.get_unchecked_mut(i) = (*s2.get_unchecked(i)).min(*s2.get_unchecked(i + d));
            }
        }
    }
}

/// Portable branchless block minimizer extractor (identical bit-for-bit output to AVX2 kernel).
#[allow(clippy::too_many_arguments)]
fn extract_block_minimizers_portable<const NEED_RIGHT: bool, const FWD_ONLY: bool>(
    seq: &PackedSeq,
    base_pos: usize,
    start_idx: usize,
    n_wins: usize,
    kp: usize,
    w: usize,
    kp_mask: u64,
    sc: &mut MinScratch,
    emit_bm: &mut [u64; 5],
) {
    let n_kps = n_wins + w - 1;
    let rc_base = 2 * (32 - kp);

    let mut idx = start_idx;
    while idx < n_kps {
        let batch_n = (n_kps - idx).min(8);
        let word = seq.load_32(base_pos + idx);
        let rc = if FWD_ONLY {
            0
        } else {
            !crate::seq::reverse_32_2bit(word)
        };
        for j in 0..batch_n {
            let kf = (word >> (2 * j)) & kp_mask;
            let (canon, kind) = if FWD_ONLY {
                (kf, 0u8)
            } else {
                let kr = (rc >> (rc_base - 2 * j)) & kp_mask;
                let knd = u8::from(kr < kf) | (u8::from(kf == kr) << 1);
                (kf.min(kr), knd)
            };
            let s = minimizer_order_hash(canon);
            let t = idx + j;
            unsafe {
                *sc.canon_kind.get_unchecked_mut(t) = canon | ((kind as u64) << 62);
                *sc.packed_l.get_unchecked_mut(t) = (s << POS_BITS) | (t as u32);
            }
        }
        idx += batch_n;
    }

    sliding_min_u32_portable(
        &sc.packed_l,
        &mut sc.s2,
        &mut sc.s4,
        &mut sc.out_l[1..],
        n_wins,
        w,
    );
    sc.out_l[0] = u32::MAX;

    if NEED_RIGHT {
        for i in 0..n_kps {
            unsafe {
                *sc.packed_r.get_unchecked_mut(i) = *sc.packed_l.get_unchecked(i) ^ POS_MASK;
            }
        }
        sliding_min_u32_portable(
            &sc.packed_r,
            &mut sc.s2,
            &mut sc.s4,
            &mut sc.out_r,
            n_wins,
            w,
        );
    }

    for i in 0..n_wins {
        unsafe {
            let cur = *sc.out_l.get_unchecked(1 + i);
            let prev = *sc.out_l.get_unchecked(i);
            if cur != prev {
                let pl = (cur & POS_MASK) as usize;
                *emit_bm.get_unchecked_mut(pl >> 6) |= 1u64 << (pl & 63);
            }
            if NEED_RIGHT {
                let cur_r = *sc.out_r.get_unchecked(i);
                if (cur ^ POS_MASK) != cur_r {
                    let pr = ((cur_r ^ POS_MASK) & POS_MASK) as usize;
                    *emit_bm.get_unchecked_mut(pr >> 6) |= 1u64 << (pr & 63);
                }
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn sliding_min_u32_avx2(
    src: &[u32; BUF_LEN],
    s2: &mut [u32; BUF_LEN],
    s4: &mut [u32; BUF_LEN],
    out: &mut [u32],
    n_wins: usize,
    w: usize,
) {
    use core::arch::x86_64::*;

    let n_kps = n_wins + w - 1;
    let sp = src.as_ptr();
    let p2 = s2.as_mut_ptr();
    let p4 = s4.as_mut_ptr();
    let op = out.as_mut_ptr();

    if w == 1 {
        let mut i = 0;
        while i < n_wins {
            unsafe {
                _mm256_storeu_si256(
                    op.add(i) as *mut __m256i,
                    _mm256_loadu_si256(sp.add(i) as *const __m256i),
                );
            }
            i += 8;
        }
        return;
    }
    if w == 2 {
        let mut i = 0;
        while i < n_wins {
            unsafe {
                let a = _mm256_loadu_si256(sp.add(i) as *const __m256i);
                let b = _mm256_loadu_si256(sp.add(i + 1) as *const __m256i);
                _mm256_storeu_si256(op.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
            }
            i += 8;
        }
        return;
    }
    // Pass 1: s2[i] = min(src[i..i+2])
    let len_s2 = n_kps - 1;
    let mut i = 0;
    while i < len_s2 {
        unsafe {
            let a = _mm256_loadu_si256(sp.add(i) as *const __m256i);
            let b = _mm256_loadu_si256(sp.add(i + 1) as *const __m256i);
            _mm256_storeu_si256(p2.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
        }
        i += 8;
    }
    if w <= 4 {
        let d = w - 2;
        let mut i = 0;
        while i < n_wins {
            unsafe {
                let a = _mm256_loadu_si256(p2.add(i) as *const __m256i);
                let b = _mm256_loadu_si256(p2.add(i + d) as *const __m256i);
                _mm256_storeu_si256(op.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
            }
            i += 8;
        }
        return;
    }
    // Pass 2: s4[i] = min(s2[i], s2[i+2]) = min(src[i..i+4])
    let len_s4 = n_kps - 3;
    let mut i = 0;
    while i < len_s4 {
        unsafe {
            let a = _mm256_loadu_si256(p2.add(i) as *const __m256i);
            let b = _mm256_loadu_si256(p2.add(i + 2) as *const __m256i);
            _mm256_storeu_si256(p4.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
        }
        i += 8;
    }
    if w <= 8 {
        let d = w - 4;
        let mut i = 0;
        while i < n_wins {
            unsafe {
                let a = _mm256_loadu_si256(p4.add(i) as *const __m256i);
                let b = _mm256_loadu_si256(p4.add(i + d) as *const __m256i);
                _mm256_storeu_si256(op.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
            }
            i += 8;
        }
    } else if w <= 12 {
        // Fused pass for k=26..=32 (w=9..=12): min(min(s4[i], s4[i+4]), s4[i + w - 4])
        let d = w - 4;
        let mut i = 0;
        while i < n_wins {
            unsafe {
                let a = _mm256_loadu_si256(p4.add(i) as *const __m256i);
                let b = _mm256_loadu_si256(p4.add(i + 4) as *const __m256i);
                let c = _mm256_loadu_si256(p4.add(i + d) as *const __m256i);
                let ab = _mm256_min_epu32(a, b);
                _mm256_storeu_si256(op.add(i) as *mut __m256i, _mm256_min_epu32(ab, c));
            }
            i += 8;
        }
    } else {
        // Defensive w in 13..=16: reuse s2 as s8.
        let len_s8 = n_kps - 7;
        let mut i = 0;
        while i < len_s8 {
            unsafe {
                let a = _mm256_loadu_si256(p4.add(i) as *const __m256i);
                let b = _mm256_loadu_si256(p4.add(i + 4) as *const __m256i);
                _mm256_storeu_si256(p2.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
            }
            i += 8;
        }
        let d = w - 8;
        let mut i = 0;
        while i < n_wins {
            unsafe {
                let a = _mm256_loadu_si256(p2.add(i) as *const __m256i);
                let b = _mm256_loadu_si256(p2.add(i + d) as *const __m256i);
                _mm256_storeu_si256(op.add(i) as *mut __m256i, _mm256_min_epu32(a, b));
            }
            i += 8;
        }
    }
}

/// Pure leaf AVX2 block minimizer extractor (1 `vzeroupper` per 256 bases, 0 closure callbacks
/// inside AVX2 code): 8-wide `load_32` canonicalization + 4-wide `vpmuludq` 1-mul ordering hash
/// + 8-wide `vpminud` binary-step sliding-window min + `vmovmskps` change/tie detection.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn extract_block_minimizers_avx2<const NEED_RIGHT: bool, const FWD_ONLY: bool>(
    seq: &PackedSeq,
    base_pos: usize,
    start_idx: usize,
    n_wins: usize,
    kp: usize,
    w: usize,
    kp_mask: u64,
    sc: &mut MinScratch,
    emit_bm: &mut [u64; 5],
) {
    use core::arch::x86_64::*;

    let n_kps = n_wins + w - 1;
    let rc_base = 2 * (32 - kp);
    let can_pal = !FWD_ONLY && (kp & 1) == 0;

    // Stage 1: 8-wide canonical kp-mer extraction + 4-wide vpmuludq 1-mul ordering hash.
    let v_mask = _mm256_set1_epi64x(kp_mask as i64);
    let v_fsh0 = _mm256_setr_epi64x(0, 2, 4, 6);
    let v_fsh1 = _mm256_setr_epi64x(8, 10, 12, 14);
    let v_rsh0 = _mm256_setr_epi64x(
        rc_base as i64,
        (rc_base - 2) as i64,
        (rc_base - 4) as i64,
        (rc_base - 6) as i64,
    );
    let v_rsh1 = _mm256_setr_epi64x(
        (rc_base - 8) as i64,
        (rc_base - 10) as i64,
        (rc_base - 12) as i64,
        (rc_base - 14) as i64,
    );
    let v_bit62 = _mm256_set1_epi64x(1i64 << 62);
    let v_bit63 = _mm256_set1_epi64x(i64::MIN);
    let v_salt32 = _mm256_set1_epi64x(ORDER_SALT as i64);
    let v_mul32 = _mm256_set1_epi64x(ORDER_MUL as i64);
    let v_score_mask = _mm256_set1_epi32(SCORE_MASK as i32);
    let v_eight = _mm256_set1_epi32(8);

    let mut idx = start_idx;
    let mut v_idx = _mm256_add_epi32(
        _mm256_set1_epi32(idx as i32),
        _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7),
    );
    while idx + 8 <= n_kps {
        let word = seq.load_32(base_pos + idx);
        let vw = _mm256_set1_epi64x(word as i64);
        let kf0 = _mm256_and_si256(_mm256_srlv_epi64(vw, v_fsh0), v_mask);
        let kf1 = _mm256_and_si256(_mm256_srlv_epi64(vw, v_fsh1), v_mask);

        let (c0, c1, ck0, ck1) = if FWD_ONLY {
            (kf0, kf1, kf0, kf1)
        } else {
            let rc = !crate::seq::reverse_32_2bit(word);
            let vrc = _mm256_set1_epi64x(rc as i64);
            let kr0 = _mm256_and_si256(_mm256_srlv_epi64(vrc, v_rsh0), v_mask);
            let kr1 = _mm256_and_si256(_mm256_srlv_epi64(vrc, v_rsh1), v_mask);
            // Safe signed compare: kf, kr < 2^42 < 2^63.
            let gt0 = _mm256_cmpgt_epi64(kf0, kr0);
            let gt1 = _mm256_cmpgt_epi64(kf1, kr1);
            let c0 = _mm256_blendv_epi8(kf0, kr0, gt0);
            let c1 = _mm256_blendv_epi8(kf1, kr1, gt1);
            let mut ck0 = _mm256_or_si256(c0, _mm256_and_si256(gt0, v_bit62));
            let mut ck1 = _mm256_or_si256(c1, _mm256_and_si256(gt1, v_bit62));
            if can_pal {
                let eq0 = _mm256_cmpeq_epi64(kf0, kr0);
                let eq1 = _mm256_cmpeq_epi64(kf1, kr1);
                ck0 = _mm256_or_si256(ck0, _mm256_and_si256(eq0, v_bit63));
                ck1 = _mm256_or_si256(ck1, _mm256_and_si256(eq1, v_bit63));
            }
            (c0, c1, ck0, ck1)
        };
        unsafe {
            _mm256_storeu_si256(sc.canon_kind.as_mut_ptr().add(idx) as *mut __m256i, ck0);
            _mm256_storeu_si256(sc.canon_kind.as_mut_ptr().add(idx + 4) as *mut __m256i, ck1);
        }

        // 1-mul hash on c0 & c1 -> pack 8x u32 into packed_l[idx..idx+8].
        let f0 = _mm256_xor_si256(_mm256_xor_si256(c0, _mm256_srli_epi64(c0, 21)), v_salt32);
        let p0 = _mm256_mul_epu32(f0, v_mul32);
        let h0 = _mm256_xor_si256(p0, _mm256_srli_epi64(p0, 32));

        let f1 = _mm256_xor_si256(_mm256_xor_si256(c1, _mm256_srli_epi64(c1, 21)), v_salt32);
        let p1 = _mm256_mul_epu32(f1, v_mul32);
        let h1 = _mm256_xor_si256(p1, _mm256_srli_epi64(p1, 32));

        let shuf = _mm256_castps_si256(_mm256_shuffle_ps(
            _mm256_castsi256_ps(h0),
            _mm256_castsi256_ps(h1),
            0x88,
        ));
        let h8 = _mm256_permute4x64_epi64(shuf, 0xD8);
        let packed8 = _mm256_or_si256(_mm256_and_si256(h8, v_score_mask), v_idx);
        unsafe {
            _mm256_storeu_si256(sc.packed_l.as_mut_ptr().add(idx) as *mut __m256i, packed8);
        }

        v_idx = _mm256_add_epi32(v_idx, v_eight);
        idx += 8;
    }

    // Scalar tail (< 8 kp-mers) using identical load_32 + !reverse_32_2bit formula.
    if idx < n_kps {
        let word = seq.load_32(base_pos + idx);
        let rc = if FWD_ONLY {
            0
        } else {
            !crate::seq::reverse_32_2bit(word)
        };
        for j in 0..(n_kps - idx) {
            let kf = (word >> (2 * j)) & kp_mask;
            let (canon, kind) = if FWD_ONLY {
                (kf, 0u8)
            } else {
                let kr = (rc >> (rc_base - 2 * j)) & kp_mask;
                let knd = u8::from(kr < kf) | (u8::from(kf == kr) << 1);
                (kf.min(kr), knd)
            };
            let s = minimizer_order_hash(canon);
            let t = idx + j;
            unsafe {
                *sc.canon_kind.get_unchecked_mut(t) = canon | ((kind as u64) << 62);
                *sc.packed_l.get_unchecked_mut(t) = (s << POS_BITS) | (t as u32);
            }
        }
    }

    // Stage 2: 8-wide AVX2 vpminud sliding-window minimum (packed_l -> out_l[1..1+n_wins]).
    unsafe {
        sliding_min_u32_avx2(
            &sc.packed_l,
            &mut sc.s2,
            &mut sc.s4,
            &mut sc.out_l[1..],
            n_wins,
            w,
        );
    }
    sc.out_l[0] = u32::MAX;

    if NEED_RIGHT {
        let v_pmask = _mm256_set1_epi32(POS_MASK as i32);
        let mut i = 0usize;
        while i < n_kps {
            unsafe {
                let v = _mm256_loadu_si256(sc.packed_l.as_ptr().add(i) as *const __m256i);
                _mm256_storeu_si256(
                    sc.packed_r.as_mut_ptr().add(i) as *mut __m256i,
                    _mm256_xor_si256(v, v_pmask),
                );
            }
            i += 8;
        }
        unsafe {
            sliding_min_u32_avx2(
                &sc.packed_r,
                &mut sc.s2,
                &mut sc.s4,
                &mut sc.out_r,
                n_wins,
                w,
            );
        }
    }

    // Stage 3: 8-wide vmovmskps change detection (+ tie detection if NEED_RIGHT) -> emit_bm.
    let v_pmask = _mm256_set1_epi32(POS_MASK as i32);
    let mut i = 0usize;
    while i < n_wins {
        let valid_mask = if i + 8 <= n_wins {
            0xFFu32
        } else {
            (1u32 << (n_wins - i)) - 1
        };
        unsafe {
            let cur = _mm256_loadu_si256(sc.out_l.as_ptr().add(1 + i) as *const __m256i);
            let prev = _mm256_loadu_si256(sc.out_l.as_ptr().add(i) as *const __m256i);
            let eq = _mm256_cmpeq_epi32(cur, prev);
            let mut changed = (!(_mm256_movemask_ps(_mm256_castsi256_ps(eq)) as u32)) & valid_mask;
            while changed != 0 {
                let b = changed.trailing_zeros() as usize;
                changed &= changed - 1;
                let pl = (*sc.out_l.get_unchecked(1 + i + b) & POS_MASK) as usize;
                *emit_bm.get_unchecked_mut(pl >> 6) |= 1u64 << (pl & 63);
            }
            if NEED_RIGHT {
                let cur_r = _mm256_loadu_si256(sc.out_r.as_ptr().add(i) as *const __m256i);
                let eq_tie = _mm256_cmpeq_epi32(_mm256_xor_si256(cur, v_pmask), cur_r);
                let mut ties =
                    (!(_mm256_movemask_ps(_mm256_castsi256_ps(eq_tie)) as u32)) & valid_mask;
                while ties != 0 {
                    let b = ties.trailing_zeros() as usize;
                    ties &= ties - 1;
                    let pr = ((*sc.out_r.get_unchecked(i + b) ^ POS_MASK) & POS_MASK) as usize;
                    *emit_bm.get_unchecked_mut(pr >> 6) |= 1u64 << (pr & 63);
                }
            }
        }
        i += 8;
    }
}

/// Extract canonical minimizers from `[start .. end)` of `seq` using block-wide branchless
/// AVX2 (`vpminud` + `vpmuludq`) or portable binary-step sliding-window minimum.
///
/// - When `EMIT_TIES == true` and `!FWD_ONLY` (reference dual-strand indexing): emits both
///   `p_left` and `p_right` per valid `k`-mer window (strictly increasing `pos`, zero duplicates),
///   guaranteeing 100% strand-invariant recall and exact per-strand `freq_limit` even under ties.
/// - When `EMIT_TIES == false` or `FWD_ONLY == true`: emits the deduplicated leftmost argmin per
///   valid `k`-mer window in strictly increasing `pos` order.
///
/// Calls `f(pos, canon, score, kind, run_start, run_end)` where `canon` is LSB-first, `score` is
/// the 23-bit `minimizer_order_hash(canon)`, and `kind` is `0` (fwd), `1` (rev), or `2` (pal).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn for_each_minimizer<const EMIT_TIES: bool, const FWD_ONLY: bool, F>(
    seq: &PackedSeq,
    start: usize,
    end: usize,
    k: usize,
    kp: usize,
    w: usize,
    soft_mask: bool,
    mut f: F,
) where
    F: FnMut(u64, u64, u32, u8, usize, usize),
{
    debug_assert!((1..=MAX_W).contains(&w));
    debug_assert!(kp + w - 1 == k && k <= 32);
    let kp_mask = (1u64 << (2 * kp)) - 1;
    let overlap = w - 1;
    #[cfg(target_arch = "x86_64")]
    let use_avx2 = has_avx2();

    let mut sc = MinScratch::new();

    seq.for_each_valid_run(start, end, soft_mask, |rs, re| {
        if re - rs < k {
            return;
        }
        let total_wins = re - rs - k + 1;
        let mut emit_bm = [0u64; 5]; // [0..4] = finalized 0..256; [4] = tail 256..256+w-1
        let mut win_done = 0usize;
        let mut carried = 0usize;

        while win_done < total_wins {
            let n_wins = (total_wins - win_done).min(BLOCK_WINS);
            let base_pos = rs + win_done;

            #[cfg(target_arch = "x86_64")]
            if use_avx2 {
                unsafe {
                    if EMIT_TIES && !FWD_ONLY {
                        extract_block_minimizers_avx2::<true, FWD_ONLY>(
                            seq,
                            base_pos,
                            carried,
                            n_wins,
                            kp,
                            w,
                            kp_mask,
                            &mut sc,
                            &mut emit_bm,
                        );
                    } else {
                        extract_block_minimizers_avx2::<false, FWD_ONLY>(
                            seq,
                            base_pos,
                            carried,
                            n_wins,
                            kp,
                            w,
                            kp_mask,
                            &mut sc,
                            &mut emit_bm,
                        );
                    }
                }
            } else if EMIT_TIES && !FWD_ONLY {
                extract_block_minimizers_portable::<true, FWD_ONLY>(
                    seq,
                    base_pos,
                    carried,
                    n_wins,
                    kp,
                    w,
                    kp_mask,
                    &mut sc,
                    &mut emit_bm,
                );
            } else {
                extract_block_minimizers_portable::<false, FWD_ONLY>(
                    seq,
                    base_pos,
                    carried,
                    n_wins,
                    kp,
                    w,
                    kp_mask,
                    &mut sc,
                    &mut emit_bm,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            {
                if EMIT_TIES && !FWD_ONLY {
                    extract_block_minimizers_portable::<true, FWD_ONLY>(
                        seq,
                        base_pos,
                        carried,
                        n_wins,
                        kp,
                        w,
                        kp_mask,
                        &mut sc,
                        &mut emit_bm,
                    );
                } else {
                    extract_block_minimizers_portable::<false, FWD_ONLY>(
                        seq,
                        base_pos,
                        carried,
                        n_wins,
                        kp,
                        w,
                        kp_mask,
                        &mut sc,
                        &mut emit_bm,
                    );
                }
            }

            let is_last = win_done + n_wins == total_wins;
            let flush_words = if is_last {
                (n_wins + overlap).div_ceil(64)
            } else {
                4
            };

            for (wi, &word) in emit_bm.iter().enumerate().take(flush_words) {
                let mut m = word;
                while m != 0 {
                    let bit = m.trailing_zeros() as usize;
                    m &= m - 1;
                    let idx = (wi << 6) + bit;
                    let ck = unsafe { *sc.canon_kind.get_unchecked(idx) };
                    let canon = ck & KP_CANON_MASK;
                    let kind = (ck >> 62) as u8;
                    let score = unsafe { *sc.packed_l.get_unchecked(idx) } >> POS_BITS;
                    f((base_pos + idx) as u64, canon, score, kind, rs, re);
                }
            }

            if !is_last {
                emit_bm[0] = emit_bm[4];
                emit_bm[1] = 0;
                emit_bm[2] = 0;
                emit_bm[3] = 0;
                emit_bm[4] = 0;
                for j in 0..overlap {
                    unsafe {
                        *sc.canon_kind.get_unchecked_mut(j) =
                            *sc.canon_kind.get_unchecked(BLOCK_WINS + j);
                        *sc.packed_l.get_unchecked_mut(j) =
                            *sc.packed_l.get_unchecked(BLOCK_WINS + j) - BLOCK_WINS as u32;
                    }
                }
                carried = overlap;
            }
            win_done += n_wins;
        }
    });
}

/// Reverse the lowest `k` 2-bit bases (LSB-first) and complement each base (`^ 3`).
#[inline(always)]
fn reverse_2bit_comp(word: u64, k: usize) -> u64 {
    let x = word.swap_bytes();
    let x = ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((x & 0xF0F0_F0F0_F0F0_F0F0) >> 4);
    let x = ((x & 0x3333_3333_3333_3333) << 2) | ((x & 0xCCCC_CCCC_CCCC_CCCC) >> 2);
    let comp = !x;
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    (comp >> (2 * (32 - k))) & mask
}

/// Canonical `(key, strand)` for a contiguous `k`-mer at `pos` in LSB-first 2-bit space.
/// Palindromes (`fwd == rc`) map to strand `0` (matching `RefIndex` per-strand counting).
#[inline(always)]
fn load_kmer_canon_strand(seq: &PackedSeq, pos: usize, k: usize, forward_only: bool) -> (u64, u8) {
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let fwd = seq.load_32(pos) & mask;
    if forward_only {
        return (fwd, 0);
    }
    let rc = reverse_2bit_comp(fwd, k);
    if fwd <= rc { (fwd, 0) } else { (rc, 1) }
}

/// Compute the 23-bit minimizer ordering score `minimizer_order_hash(canon_kp)` for the LSB-first
/// `kp`-mer at `pos` in `seq`.
#[inline(always)]
fn load_kp_score(seq: &PackedSeq, pos: usize, kp: usize, kp_mask: u64, forward_only: bool) -> u32 {
    let kf = seq.load_32(pos) & kp_mask;
    let canon = if forward_only {
        kf
    } else {
        kf.min(reverse_2bit_comp(kf, kp))
    };
    minimizer_order_hash(canon)
}

/// Check whether the valid `k`-mer window `[r .. r + k)` has minimum canonical `kp`-mer score
/// strictly equal to `target_score` (using LSB-first `canon` + `minimizer_order_hash`).
#[cfg(test)]
#[inline(always)]
fn window_has_min_score(
    seq: &PackedSeq,
    r: usize,
    kp: usize,
    w: usize,
    kp_mask: u64,
    forward_only: bool,
    target_score: u32,
) -> bool {
    let word = seq.load_32(r);
    let rc = if forward_only {
        0
    } else {
        !crate::seq::reverse_32_2bit(word)
    };
    let rc_base = 2 * (32 - kp);
    let mut found = false;
    for j in 0..w {
        let kf = (word >> (2 * j)) & kp_mask;
        let canon = if forward_only {
            kf
        } else {
            let kr = (rc >> (rc_base - 2 * j)) & kp_mask;
            kf.min(kr)
        };
        let s = minimizer_order_hash(canon);
        if s < target_score {
            return false;
        }
        if s == target_score {
            found = true;
        }
    }
    found
}

/// $O(w)$ rolling interval computation (`W3-F8b`): given a valid `kp`-mer position `p_r` inside
/// maximal valid run `[rs, re)` (`re - rs >= k`) with score `min_score`, return the exact
/// inclusive range `(r_min, r_max)` of valid `k`-mer start positions `r in [r_lo..=r_hi]` whose
/// window `[r..r+w)` has minimum canonical `kp`-mer score equal to `min_score`.
///
/// Returns `(1, 0)` (empty range) if no valid `k`-mer window achieves `min_score`.
#[inline]
#[allow(clippy::too_many_arguments)]
fn minimizer_valid_kmer_range(
    seq: &PackedSeq,
    p_r: usize,
    rs: usize,
    re: usize,
    k: usize,
    kp: usize,
    w: usize,
    kp_mask: u64,
    forward_only: bool,
    min_score: u32,
) -> (usize, usize) {
    let r_lo = p_r.saturating_sub(w - 1).max(rs);
    let r_hi = p_r.min(re - k);
    if r_lo > r_hi {
        return (1, 0);
    }
    let rc_base = 2 * (32 - kp);
    let mut r_min = r_lo;
    let mut r_max = r_hi;

    // Right scan: q = p_r + 1 ..= r_hi + w - 1.
    // Rightmost base is at r_hi + w - 1 + kp - 1 = r_hi + k - 1 <= p_r + 31, so one `load_32(p_r)`
    // covers the entire right scan in LSB-first order.
    let right_end = r_hi + w - 1;
    if right_end > p_r {
        let right_steps = right_end - p_r;
        let w_right = seq.load_32(p_r);
        let rc_right = if forward_only {
            0
        } else {
            !crate::seq::reverse_32_2bit(w_right)
        };
        for j in 1..=right_steps {
            let kf = (w_right >> (2 * j)) & kp_mask;
            let canon = if forward_only {
                kf
            } else {
                let kr = (rc_right >> (rc_base - 2 * j)) & kp_mask;
                kf.min(kr)
            };
            if minimizer_order_hash(canon) < min_score {
                let q = p_r + j;
                if q < w {
                    return (1, 0);
                }
                r_max = r_max.min(q - w);
                break;
            }
        }
    }

    if r_min > r_max {
        return (1, 0);
    }

    // Left scan: q = p_r - 1 down to r_lo.
    // Rightmost base is at (p_r - 1) + kp - 1 <= r_lo + k - 2 <= r_lo + 30, so one `load_32(r_lo)`
    // covers the entire left scan in LSB-first order.
    if p_r > r_lo {
        let left_steps = p_r - r_lo;
        let w_left = seq.load_32(r_lo);
        let rc_left = if forward_only {
            0
        } else {
            !crate::seq::reverse_32_2bit(w_left)
        };
        for j in (0..left_steps).rev() {
            let kf = (w_left >> (2 * j)) & kp_mask;
            let canon = if forward_only {
                kf
            } else {
                let kr = (rc_left >> (rc_base - 2 * j)) & kp_mask;
                kf.min(kr)
            };
            if minimizer_order_hash(canon) < min_score {
                r_min = r_lo + j + 1;
                break;
            }
        }
    }

    (r_min, r_max)
}

/// Return the maximal `[left, right)` interval around valid position `pos` that does not overlap
/// any sorted non-overlapping invalid run in `runs`.
#[inline(always)]
fn valid_bounds_in_runs(runs: &[(u64, u64)], pos: usize, seq_len: usize) -> (usize, usize) {
    if runs.is_empty() {
        return (0, seq_len);
    }
    let pos64 = pos as u64;
    let idx = runs.partition_point(|&(_, e)| e <= pos64);
    let left = if idx > 0 { runs[idx - 1].1 as usize } else { 0 };
    let right = if idx < runs.len() {
        runs[idx].0 as usize
    } else {
        seq_len
    };
    (left, right)
}

/// Return the maximal valid seed run `[rs, re)` containing valid position `pos` in `seq`.
#[inline(always)]
fn valid_run_around(seq: &PackedSeq, pos: usize, soft_mask: bool) -> (usize, usize) {
    let (mut lo, mut hi) = valid_bounds_in_runs(&seq.special_runs, pos, seq.len);
    if soft_mask && !seq.soft_runs.is_empty() {
        let (slo, shi) = valid_bounds_in_runs(&seq.soft_runs, pos, seq.len);
        lo = lo.max(slo);
        hi = hi.min(shi);
    }
    (lo, hi)
}

#[inline]
fn choose_bucket_bits(total_entries: usize) -> u32 {
    let target = (total_entries / 4).max(NUM_PARTITIONS);
    target.next_power_of_two().trailing_zeros().clamp(8, 26)
}

#[inline]
fn make_chunk_ranges(seq_len: usize, k: usize, n_chunks: usize) -> Vec<(usize, usize)> {
    if seq_len < k {
        return Vec::new();
    }
    let max_starts = seq_len - k + 1;
    let chunk_starts = max_starts.div_ceil(n_chunks.max(1)).max(1);
    let mut ranges = Vec::with_capacity(n_chunks);
    let mut s = 0usize;
    while s < max_starts {
        let e_start = (s + chunk_starts).min(max_starts);
        let e = (e_start + k - 1).min(seq_len);
        ranges.push((s, e));
        s = e_start;
    }
    ranges
}

// ---------------------------------------------------------------------------
// High-Frequency Exact-k Overflow Hash-CSR (W3-F6)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct OverflowKmerCsr {
    hash_shift: u32,
    offsets: Vec<u32>,
    tags: Vec<u32>,
    pos_strand: PosStrandVec,
    bloom: Option<BlockedBloom>,
    kept_entries: usize,
}

impl OverflowKmerCsr {
    fn heap_bytes(&self) -> usize {
        self.offsets.capacity() * 4
            + self.tags.capacity() * 4
            + self.pos_strand.heap_bytes()
            + self.bloom.as_ref().map_or(0, BlockedBloom::heap_bytes)
    }

    fn build(mut items: Vec<(u64, u64)>, use_u32_ps: bool) -> Option<Self> {
        if items.is_empty() {
            return None;
        }
        let bucket_bits = (items.len() / 4)
            .max(64)
            .next_power_of_two()
            .trailing_zeros()
            .clamp(6, 24);
        let hash_shift = 64 - bucket_bits;
        let n_buckets = 1usize << bucket_bits;

        for e in &mut items {
            e.0 = idx_hash_to_sort_key(e.0, hash_shift);
        }
        items.par_sort_unstable();
        items.dedup(); // Deduplicate exact k-mers emitted via tied overflowed minimizers.

        let kept_entries = items.len();
        let mut offsets = vec![0u32; n_buckets + 1];
        let mut tags = Vec::with_capacity(kept_entries);
        let mut ps_u32 = if use_u32_ps {
            Vec::with_capacity(kept_entries)
        } else {
            Vec::new()
        };
        let mut ps_u64 = if use_u32_ps {
            Vec::new()
        } else {
            Vec::with_capacity(kept_entries)
        };

        // Optional flat cache-line Blocked Bloom (10 bits/entry, capped at 16 MiB).
        let (mut bloom_blocks, bloom_shift, bloom_k_h, bloom_sub_bits) = if kept_entries >= 1024 {
            let n_blocks = kept_entries
                .saturating_mul(10)
                .div_ceil(512)
                .next_power_of_two()
                .clamp(16, 1 << 18);
            let b_bits = n_blocks.trailing_zeros();
            let b_shift = 64 - b_bits;
            let bpe = (n_blocks * 512) as f64 / (kept_entries as f64);
            let kh = (bpe * std::f64::consts::LN_2).round().clamp(1.0, 7.0) as u8;
            (vec![BloomBlock([0u64; 8]); n_blocks], b_shift, kh, b_bits)
        } else {
            (Vec::new(), 64u32, 0u8, 0u32)
        };

        let mut prev_sk = u64::MAX;
        for &(sk, ps) in &items {
            let kh = sort_key_to_idx_hash(sk, hash_shift);
            let b = (kh >> hash_shift) as usize;
            let tag = kh as u32;
            offsets[b + 1] += 1;
            tags.push(tag);
            if use_u32_ps {
                ps_u32.push(ps as u32);
            } else {
                ps_u64.push(ps);
            }
            if !bloom_blocks.is_empty() && sk != prev_sk {
                let blk = (kh >> bloom_shift) as usize;
                insert_in_bloom_block(&mut bloom_blocks[blk], tag, bloom_k_h as usize);
                prev_sk = sk;
            }
        }
        for b in 0..n_buckets {
            offsets[b + 1] += offsets[b];
        }

        let bloom = if !bloom_blocks.is_empty() {
            Some(BlockedBloom {
                blocks: bloom_blocks,
                sub_bits: bloom_sub_bits,
                block_shift: bloom_shift,
                k_h: bloom_k_h,
            })
        } else {
            None
        };
        let pos_strand = if use_u32_ps {
            PosStrandVec::U32(ps_u32)
        } else {
            PosStrandVec::U64(ps_u64)
        };

        Some(Self {
            hash_shift,
            offsets,
            tags,
            pos_strand,
            bloom,
            kept_entries,
        })
    }
}

// ---------------------------------------------------------------------------
// Flat 2-Level Split-Tag Radix MinimizerIndex Builder (W3-F8a + W3-F6)
// ---------------------------------------------------------------------------

trait MinPackedPs: Copy + Default + Ord + Send + Sync + 'static {
    fn pack(pos: u64, strand_bit: u8) -> Self;
    fn pos(self) -> usize;
    fn sentinel() -> Self;
    fn into_ps_vec(v: Vec<Self>) -> PosStrandVec;
}

impl MinPackedPs for u32 {
    #[inline(always)]
    fn pack(pos: u64, strand_bit: u8) -> Self {
        ((pos as u32) << 1) | (strand_bit as u32)
    }

    #[inline(always)]
    fn pos(self) -> usize {
        (self >> 1) as usize
    }

    #[inline(always)]
    fn sentinel() -> Self {
        OVERFLOW_SENTINEL_U32
    }

    fn into_ps_vec(v: Vec<Self>) -> PosStrandVec {
        PosStrandVec::U32(v)
    }
}

impl MinPackedPs for u64 {
    #[inline(always)]
    fn pack(pos: u64, strand_bit: u8) -> Self {
        (pos << 1) | (strand_bit as u64)
    }

    #[inline(always)]
    fn pos(self) -> usize {
        (self >> 1) as usize
    }

    #[inline(always)]
    fn sentinel() -> Self {
        OVERFLOW_SENTINEL_U64
    }

    fn into_ps_vec(v: Vec<Self>) -> PosStrandVec {
        PosStrandVec::U64(v)
    }
}

#[repr(align(64))]
struct PaddedAtomicUsize(AtomicUsize);

#[derive(Clone, Copy)]
struct SendMutPtr<T>(*mut T);
unsafe impl<T> Send for SendMutPtr<T> {}
unsafe impl<T> Sync for SendMutPtr<T> {}

#[derive(Clone, Copy)]
struct MinOverflowEntry<P: MinPackedPs> {
    lp: u8,
    tag: u32,
    mid: [u8; 3],
    ps: P,
}

type MinOverflowAppend<P> = (usize, Vec<u32>, Vec<P>);

/// Insertion sort parallel `(keys, vals)` slices by lexicographic `(key, val)` so duplicate
/// `(h_rot, ps)` entries from overlapping Pass-1 chunks are guaranteed adjacent.
#[inline]
fn insertion_sort_pair_by_kv<P: MinPackedPs>(keys: &mut [u64], vals: &mut [P]) {
    let n = keys.len();
    for i in 1..n {
        let k = unsafe { *keys.get_unchecked(i) };
        let v = unsafe { *vals.get_unchecked(i) };
        let mut j = i;
        while j > 0 {
            let pk = unsafe { *keys.get_unchecked(j - 1) };
            let pv = unsafe { *vals.get_unchecked(j - 1) };
            if pk > k || (pk == k && pv > v) {
                unsafe {
                    *keys.get_unchecked_mut(j) = pk;
                    *vals.get_unchecked_mut(j) = pv;
                }
                j -= 1;
            } else {
                break;
            }
        }
        unsafe {
            *keys.get_unchecked_mut(j) = k;
            *vals.get_unchecked_mut(j) = v;
        }
    }
}

#[cold]
#[inline(never)]
fn sort_large_subbucket_by_kv<P: MinPackedPs>(
    keys: &mut [u64],
    vals: &mut [P],
    sort_buf: &mut Vec<(u64, P)>,
) {
    sort_buf.clear();
    sort_buf.reserve(keys.len());
    for (&k, &v) in keys.iter().zip(vals.iter()) {
        sort_buf.push((k, v));
    }
    sort_buf.sort_unstable();
    for (i, &(k, v)) in sort_buf.iter().enumerate() {
        unsafe {
            *keys.get_unchecked_mut(i) = k;
            *vals.get_unchecked_mut(i) = v;
        }
    }
}

#[inline(always)]
unsafe fn insert_min_bloom_const<const K_H: usize>(
    bloom_ptr: *mut BloomBlock,
    h_rot: u64,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
) {
    if K_H > 0 {
        let blk_local = ((h_rot as usize) >> bloom_sub_shift) & bloom_sub_mask;
        let tag = (h_rot >> 32) as u32;
        unsafe {
            insert_in_bloom_block(&mut *bloom_ptr.add(blk_local), tag, K_H);
        }
    }
}

struct Pass2Scratch<P: MinPackedPs> {
    sub_off: Vec<u32>,
    scratch_k: Vec<u64>,
    scratch_p: Vec<P>,
    sort_buf: Vec<(u64, P)>,
    cand_r: Vec<usize>,
    kmer_items: Vec<(u64, u8, usize)>,
    dropped_local: Vec<usize>,
}

impl<P: MinPackedPs> Pass2Scratch<P> {
    fn new(sub_buckets: usize) -> Self {
        Self {
            sub_off: vec![0u32; sub_buckets],
            scratch_k: Vec::new(),
            scratch_p: Vec::new(),
            sort_buf: Vec::new(),
            cand_r: Vec::new(),
            kmer_items: Vec::new(),
            dropped_local: Vec::new(),
        }
    }
}

struct PartCompactOut<P: MinPackedPs> {
    kept: usize,
    dropped_r: Vec<usize>,
    ov_kmers: Vec<(u64, u64)>,
    append: Option<MinOverflowAppend<P>>,
}

/// Shared Pass-2 sub-bucket sort, cross-chunk dedup, 3-way run handler (`W3-F6` + `W3-F8b`),
/// and compile-time-specialized Bloom insertion writing into `(tag_ptr, ps_ptr)`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn compact_min_buckets<P: MinPackedPs, const K_H: usize>(
    seq: &PackedSeq,
    word_length: usize,
    kp: usize,
    window: usize,
    kp_mask: u64,
    soft_mask: bool,
    forward_only: bool,
    freq_limit: u64,
    part_off: usize,
    sub_buckets: usize,
    bloom_ptr: *mut BloomBlock,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    poff: &mut [MaybeUninit<u32>],
    tag_ptr: *mut u32,
    ps_ptr: *mut P,
    sc: &mut Pass2Scratch<P>,
    dropped_r: &mut Vec<usize>,
    ov_kmers: &mut Vec<(u64, u64)>,
) -> usize {
    let mut write_idx = 0usize;
    let mut sub_start = 0usize;

    for (slot, &sub_end_u32) in poff.iter_mut().zip(sc.sub_off.iter()) {
        slot.write((part_off + write_idx) as u32);
        let sub_end = sub_end_u32 as usize;
        let cnt = sub_end - sub_start;
        if cnt == 0 {
            continue;
        }
        if cnt == 1 {
            // `freq_limit >= 1` and `MINIMIZER_OVERFLOW_CAP >= 1` always hold here.
            let h_rot = unsafe { *sc.scratch_k.get_unchecked(sub_start) };
            let ps = unsafe { *sc.scratch_p.get_unchecked(sub_start) };
            unsafe {
                *tag_ptr.add(write_idx) = (h_rot >> 32) as u32;
                *ps_ptr.add(write_idx) = ps;
                insert_min_bloom_const::<K_H>(bloom_ptr, h_rot, bloom_sub_shift, bloom_sub_mask);
            }
            write_idx += 1;
            sub_start = sub_end;
            continue;
        }

        if cnt <= 24 {
            insertion_sort_pair_by_kv(
                unsafe { sc.scratch_k.get_unchecked_mut(sub_start..sub_end) },
                unsafe { sc.scratch_p.get_unchecked_mut(sub_start..sub_end) },
            );
        } else {
            sort_large_subbucket_by_kv(
                unsafe { sc.scratch_k.get_unchecked_mut(sub_start..sub_end) },
                unsafe { sc.scratch_p.get_unchecked_mut(sub_start..sub_end) },
                &mut sc.sort_buf,
            );
        }

        let mut i = sub_start;
        while i < sub_end {
            let h_rot = unsafe { *sc.scratch_k.get_unchecked(i) };
            let mut j = i + 1;
            while j < sub_end && unsafe { *sc.scratch_k.get_unchecked(j) } == h_rot {
                j += 1;
            }

            // Compact duplicate `(h_rot, ps)` from overlapping Pass-1 chunks in-place.
            let mut b_dedup = i + 1;
            let mut prev_ps = unsafe { *sc.scratch_p.get_unchecked(i) };
            for idx in (i + 1)..j {
                let cur_ps = unsafe { *sc.scratch_p.get_unchecked(idx) };
                if cur_ps != prev_ps {
                    unsafe { *sc.scratch_p.get_unchecked_mut(b_dedup) = cur_ps };
                    b_dedup += 1;
                    prev_ps = cur_ps;
                }
            }
            let run_len = b_dedup - i;
            let tag = (h_rot >> 32) as u32;
            let before_write = write_idx;

            if (run_len as u64) <= freq_limit && run_len <= MINIMIZER_OVERFLOW_CAP {
                // Branch 1: Normal low-frequency minimizer bucket (common fast path).
                for idx in i..b_dedup {
                    unsafe {
                        *tag_ptr.add(write_idx) = tag;
                        *ps_ptr.add(write_idx) = *sc.scratch_p.get_unchecked(idx);
                    }
                    write_idx += 1;
                }
            } else if (run_len as u64) <= freq_limit {
                // Branch 2: `MINIMIZER_OVERFLOW_CAP < run_len <= freq_limit`.
                // No k-mer can exceed `freq_limit`; route surviving exact k-mers to OverflowKmerCsr.
                let p_first = unsafe { *sc.scratch_p.get_unchecked(i) }.pos();
                let min_score = load_kp_score(seq, p_first, kp, kp_mask, forward_only);

                sc.cand_r.clear();
                for idx in i..b_dedup {
                    let p_r = unsafe { *sc.scratch_p.get_unchecked(idx) }.pos();
                    let (rs, re) = valid_run_around(seq, p_r, soft_mask);
                    if re - rs < word_length {
                        continue;
                    }
                    let (r_min, r_max) = minimizer_valid_kmer_range(
                        seq,
                        p_r,
                        rs,
                        re,
                        word_length,
                        kp,
                        window,
                        kp_mask,
                        forward_only,
                        min_score,
                    );
                    if r_min <= r_max {
                        sc.cand_r.extend(r_min..=r_max);
                    }
                }
                sc.cand_r.sort_unstable();
                sc.cand_r.dedup();

                if !sc.cand_r.is_empty() {
                    unsafe {
                        *tag_ptr.add(write_idx) = tag;
                        *ps_ptr.add(write_idx) = P::sentinel();
                    }
                    write_idx += 1;
                    ov_kmers.reserve(sc.cand_r.len());
                    for &r in &sc.cand_r {
                        let (kc, ks) = load_kmer_canon_strand(seq, r, word_length, forward_only);
                        ov_kmers.push((fmix64(kc), ((r as u64) << 1) | (ks as u64)));
                    }
                }
            } else {
                // Branch 3: `run_len > freq_limit` — exact per-strand k-mer `freq_limit` check
                // using O(w) rolling interval computation + HF overflow routing if `> M_CAP`.
                let p_first = unsafe { *sc.scratch_p.get_unchecked(i) }.pos();
                let min_score = load_kp_score(seq, p_first, kp, kp_mask, forward_only);

                sc.cand_r.clear();
                for idx in i..b_dedup {
                    let p_r = unsafe { *sc.scratch_p.get_unchecked(idx) }.pos();
                    let (rs, re) = valid_run_around(seq, p_r, soft_mask);
                    if re - rs < word_length {
                        continue;
                    }
                    let (r_min, r_max) = minimizer_valid_kmer_range(
                        seq,
                        p_r,
                        rs,
                        re,
                        word_length,
                        kp,
                        window,
                        kp_mask,
                        forward_only,
                        min_score,
                    );
                    if r_min <= r_max {
                        sc.cand_r.extend(r_min..=r_max);
                    }
                }
                sc.cand_r.sort_unstable();
                sc.cand_r.dedup();

                sc.kmer_items.clear();
                sc.kmer_items.reserve(sc.cand_r.len());
                for &r in &sc.cand_r {
                    let (kc, ks) = load_kmer_canon_strand(seq, r, word_length, forward_only);
                    sc.kmer_items.push((kc, ks, r));
                }
                sc.kmer_items.sort_unstable();

                sc.dropped_local.clear();
                let route_overflow = run_len > MINIMIZER_OVERFLOW_CAP;
                let mut ov_emitted = false;
                let mut ka = 0usize;
                while ka < sc.kmer_items.len() {
                    let (kc, ks, _) = sc.kmer_items[ka];
                    let mut kb = ka + 1;
                    while kb < sc.kmer_items.len()
                        && sc.kmer_items[kb].0 == kc
                        && sc.kmer_items[kb].1 == ks
                    {
                        kb += 1;
                    }
                    if ((kb - ka) as u64) > freq_limit {
                        for &(_, _, r) in &sc.kmer_items[ka..kb] {
                            sc.dropped_local.push(r);
                            dropped_r.push(r);
                        }
                    } else if route_overflow {
                        let kh = fmix64(kc);
                        ov_kmers.reserve(kb - ka);
                        for &(_, _, r) in &sc.kmer_items[ka..kb] {
                            ov_kmers.push((kh, ((r as u64) << 1) | (ks as u64)));
                        }
                        ov_emitted = true;
                    }
                    ka = kb;
                }

                if route_overflow {
                    if ov_emitted {
                        unsafe {
                            *tag_ptr.add(write_idx) = tag;
                            *ps_ptr.add(write_idx) = P::sentinel();
                        }
                        write_idx += 1;
                    }
                } else if sc.dropped_local.is_empty() {
                    for idx in i..b_dedup {
                        unsafe {
                            *tag_ptr.add(write_idx) = tag;
                            *ps_ptr.add(write_idx) = *sc.scratch_p.get_unchecked(idx);
                        }
                        write_idx += 1;
                    }
                } else {
                    sc.dropped_local.sort_unstable();
                    for idx in i..b_dedup {
                        let ps = unsafe { *sc.scratch_p.get_unchecked(idx) };
                        let p_r = ps.pos();
                        let r_lo = p_r.saturating_sub(window - 1);
                        let s_idx = sc.cand_r.partition_point(|&r| r < r_lo);
                        let mut survives = false;
                        for &r in &sc.cand_r[s_idx..] {
                            if r > p_r {
                                break;
                            }
                            if sc.dropped_local.binary_search(&r).is_err() {
                                survives = true;
                                break;
                            }
                        }
                        if survives {
                            unsafe {
                                *tag_ptr.add(write_idx) = tag;
                                *ps_ptr.add(write_idx) = ps;
                            }
                            write_idx += 1;
                        }
                    }
                }
            }

            if write_idx > before_write {
                unsafe {
                    insert_min_bloom_const::<K_H>(
                        bloom_ptr,
                        h_rot,
                        bloom_sub_shift,
                        bloom_sub_mask,
                    );
                }
            }

            i = j;
        }

        sub_start = sub_end;
    }

    poff[sub_buckets].write((part_off + write_idx) as u32);
    write_idx
}

#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
unsafe fn process_min_partition_overflow<P: MinPackedPs, const K_H: usize>(
    seq: &PackedSeq,
    word_length: usize,
    kp: usize,
    window: usize,
    kp_mask: u64,
    soft_mask: bool,
    forward_only: bool,
    freq_limit: u64,
    lp: usize,
    part_off: usize,
    flat_len: usize,
    slot_cap: usize,
    pmid: &[u8],
    ptag: &mut [MaybeUninit<u32>],
    pp: &mut [MaybeUninit<P>],
    poff: &mut [MaybeUninit<u32>],
    ov: &[MinOverflowEntry<P>],
    bloom_ptr: *mut BloomBlock,
    sub_buckets: usize,
    sub_shift_hi: u32,
    sub_mask: usize,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    sc: &mut Pass2Scratch<P>,
) -> PartCompactOut<P> {
    let n_total = flat_len + ov.len();
    sc.sub_off.fill(0);
    let m_ptr = pmid.as_ptr();
    for i in 0..flat_len {
        let mid24 = u32::from_le(unsafe { (m_ptr.add(i * 3) as *const u32).read_unaligned() })
            & 0x00FF_FFFF;
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        unsafe { *sc.sub_off.get_unchecked_mut(sub) += 1 };
    }
    for e in ov {
        let mid24 = u32::from_le_bytes([e.mid[0], e.mid[1], e.mid[2], 0]);
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        unsafe { *sc.sub_off.get_unchecked_mut(sub) += 1 };
    }

    let mut acc = 0u32;
    for slot in &mut sc.sub_off {
        let c = *slot;
        *slot = acc;
        acc += c;
    }

    if sc.scratch_k.len() < n_total {
        sc.scratch_k.resize(n_total, 0);
        sc.scratch_p.resize(n_total, P::default());
    }
    for i in 0..flat_len {
        let mid24 = u32::from_le(unsafe { (m_ptr.add(i * 3) as *const u32).read_unaligned() })
            & 0x00FF_FFFF;
        let tag = unsafe { ptag.get_unchecked(i).assume_init() };
        let ps = unsafe { pp.get_unchecked(i).assume_init() };
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        let dst = unsafe {
            let s = sc.sub_off.get_unchecked_mut(sub);
            let d = *s as usize;
            *s = (d + 1) as u32;
            d
        };
        unsafe {
            *sc.scratch_k.get_unchecked_mut(dst) = ((tag as u64) << 32) | (mid24 as u64);
            *sc.scratch_p.get_unchecked_mut(dst) = ps;
        }
    }
    for e in ov {
        let mid24 = u32::from_le_bytes([e.mid[0], e.mid[1], e.mid[2], 0]);
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        let dst = unsafe {
            let s = sc.sub_off.get_unchecked_mut(sub);
            let d = *s as usize;
            *s = (d + 1) as u32;
            d
        };
        unsafe {
            *sc.scratch_k.get_unchecked_mut(dst) = ((e.tag as u64) << 32) | (mid24 as u64);
            *sc.scratch_p.get_unchecked_mut(dst) = e.ps;
        }
    }

    let mut tmp_tags = vec![0u32; n_total];
    let mut tmp_ps = vec![P::default(); n_total];
    let mut dropped_r = Vec::new();
    let mut ov_kmers = Vec::new();
    let write_idx = unsafe {
        compact_min_buckets::<P, K_H>(
            seq,
            word_length,
            kp,
            window,
            kp_mask,
            soft_mask,
            forward_only,
            freq_limit,
            0,
            sub_buckets,
            bloom_ptr,
            bloom_sub_shift,
            bloom_sub_mask,
            poff,
            tmp_tags.as_mut_ptr(),
            tmp_ps.as_mut_ptr(),
            sc,
            &mut dropped_r,
            &mut ov_kmers,
        )
    };

    if write_idx <= slot_cap {
        let t_ptr = ptag.as_mut_ptr() as *mut u32;
        let p_ptr = pp.as_mut_ptr() as *mut P;
        unsafe {
            std::ptr::copy_nonoverlapping(tmp_tags.as_ptr(), t_ptr, write_idx);
            std::ptr::copy_nonoverlapping(tmp_ps.as_ptr(), p_ptr, write_idx);
        }
        for slot in &mut ptag[write_idx..slot_cap] {
            slot.write(0);
        }
        for slot in &mut pp[write_idx..slot_cap] {
            slot.write(P::default());
        }
        let base_u32 = part_off as u32;
        for slot in poff.iter_mut() {
            let v = unsafe { slot.assume_init() };
            slot.write(v + base_u32);
        }
        PartCompactOut {
            kept: write_idx,
            dropped_r,
            ov_kmers,
            append: None,
        }
    } else {
        for slot in &mut ptag[..slot_cap] {
            slot.write(0);
        }
        for slot in &mut pp[..slot_cap] {
            slot.write(P::default());
        }
        tmp_tags.truncate(write_idx);
        tmp_ps.truncate(write_idx);
        PartCompactOut {
            kept: write_idx,
            dropped_r,
            ov_kmers,
            append: Some((lp, tmp_tags, tmp_ps)),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_min_pass2_const<P: MinPackedPs, const K_H: usize>(
    seq: &PackedSeq,
    word_length: usize,
    kp: usize,
    window: usize,
    kp_mask: u64,
    soft_mask: bool,
    forward_only: bool,
    freq_limit: u64,
    n_parts: usize,
    slot_cap: usize,
    sub_buckets: usize,
    sub_shift_hi: u32,
    sub_mask: usize,
    bloom_blocks_per_part: usize,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    flat_lens: &[usize],
    ov_by_part: &[Vec<MinOverflowEntry<P>>],
    coarse_mid: &[u8],
    tags: &mut [MaybeUninit<u32>],
    pos_strand: &mut [MaybeUninit<P>],
    offsets: &mut [MaybeUninit<u32>],
    bloom_blocks: &mut [BloomBlock],
) -> Vec<PartCompactOut<P>> {
    let mid_stride = slot_cap * 3;
    let off_stride = sub_buckets + 1;
    let bloom_base_ptr = SendMutPtr(bloom_blocks.as_mut_ptr());

    tags.par_chunks_mut(slot_cap)
        .zip(pos_strand.par_chunks_mut(slot_cap))
        .zip(offsets.par_chunks_mut(off_stride))
        .enumerate()
        .map_init(
            || Pass2Scratch::<P>::new(sub_buckets),
            |sc, (lp, ((ptag, pp), poff))| {
                let _ = n_parts;
                let bp = bloom_base_ptr;
                let part_off = lp * slot_cap;
                let flat_len = flat_lens[lp];
                let ov = &ov_by_part[lp];
                let bloom_ptr = if bloom_blocks_per_part > 0 {
                    unsafe { bp.0.add(lp * bloom_blocks_per_part) }
                } else {
                    std::ptr::null_mut()
                };

                if flat_len == 0 && ov.is_empty() {
                    for slot in poff.iter_mut() {
                        slot.write(part_off as u32);
                    }
                    for slot in &mut ptag[..slot_cap] {
                        slot.write(0);
                    }
                    for slot in &mut pp[..slot_cap] {
                        slot.write(P::default());
                    }
                    return PartCompactOut {
                        kept: 0,
                        dropped_r: Vec::new(),
                        ov_kmers: Vec::new(),
                        append: None,
                    };
                }

                let pmid = &coarse_mid[lp * mid_stride..lp * mid_stride + mid_stride + 1];

                if !ov.is_empty() {
                    return unsafe {
                        process_min_partition_overflow::<P, K_H>(
                            seq,
                            word_length,
                            kp,
                            window,
                            kp_mask,
                            soft_mask,
                            forward_only,
                            freq_limit,
                            lp,
                            part_off,
                            flat_len,
                            slot_cap,
                            pmid,
                            ptag,
                            pp,
                            poff,
                            ov,
                            bloom_ptr,
                            sub_buckets,
                            sub_shift_hi,
                            sub_mask,
                            bloom_sub_shift,
                            bloom_sub_mask,
                            sc,
                        )
                    };
                }

                // Normal in-place path: flat_len <= slot_cap, zero overflow.
                sc.sub_off.fill(0);
                let m_ptr = pmid.as_ptr();
                for i in 0..flat_len {
                    let mid24 =
                        u32::from_le(unsafe { (m_ptr.add(i * 3) as *const u32).read_unaligned() })
                            & 0x00FF_FFFF;
                    let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
                    unsafe { *sc.sub_off.get_unchecked_mut(sub) += 1 };
                }

                let mut acc = 0u32;
                for slot in &mut sc.sub_off {
                    let c = *slot;
                    *slot = acc;
                    acc += c;
                }

                if sc.scratch_k.len() < flat_len {
                    sc.scratch_k.resize(flat_len, 0);
                    sc.scratch_p.resize(flat_len, P::default());
                }
                for i in 0..flat_len {
                    let mid24 =
                        u32::from_le(unsafe { (m_ptr.add(i * 3) as *const u32).read_unaligned() })
                            & 0x00FF_FFFF;
                    let tag = unsafe { ptag.get_unchecked(i).assume_init() };
                    let ps = unsafe { pp.get_unchecked(i).assume_init() };
                    let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
                    let dst = unsafe {
                        let s = sc.sub_off.get_unchecked_mut(sub);
                        let d = *s as usize;
                        *s = (d + 1) as u32;
                        d
                    };
                    unsafe {
                        *sc.scratch_k.get_unchecked_mut(dst) =
                            ((tag as u64) << 32) | (mid24 as u64);
                        *sc.scratch_p.get_unchecked_mut(dst) = ps;
                    }
                }

                let mut dropped_r = Vec::new();
                let mut ov_kmers = Vec::new();
                let write_idx = unsafe {
                    compact_min_buckets::<P, K_H>(
                        seq,
                        word_length,
                        kp,
                        window,
                        kp_mask,
                        soft_mask,
                        forward_only,
                        freq_limit,
                        part_off,
                        sub_buckets,
                        bloom_ptr,
                        bloom_sub_shift,
                        bloom_sub_mask,
                        poff,
                        ptag.as_mut_ptr() as *mut u32,
                        pp.as_mut_ptr() as *mut P,
                        sc,
                        &mut dropped_r,
                        &mut ov_kmers,
                    )
                };

                for slot in &mut ptag[write_idx..slot_cap] {
                    slot.write(0);
                }
                for slot in &mut pp[write_idx..slot_cap] {
                    slot.write(P::default());
                }

                PartCompactOut {
                    kept: write_idx,
                    dropped_r,
                    ov_kmers,
                    append: None,
                }
            },
        )
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn build_minimizer_typed<P: MinPackedPs>(
    reference: &Genome,
    part_start: u32,
    part_end: u32,
    word_length: usize,
    kp: usize,
    window: usize,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    use_u32_ps: bool,
) -> Result<MinimizerIndex> {
    let seq = &reference.seq;
    let valid_kmers = count_valid_kmers(reference, word_length, soft_mask);
    if valid_kmers == 0 {
        return Ok(MinimizerIndex::empty(
            word_length,
            kp,
            window,
            freq_limit,
            soft_mask,
            forward_only,
            part_start,
            part_end,
            use_u32_ps,
        ));
    }

    let n_parts = (part_end - part_start) as usize;
    let part_filter = part_start > 0 || part_end < NUM_PARTITIONS as u32;
    let n_threads = rayon::current_num_threads().max(1);
    let n_chunks = if seq.len >= 16_384 {
        (seq.len / (1 << 20)).clamp(n_threads * 4, 64)
    } else {
        n_threads * 2
    };
    let chunks = make_chunk_ranges(seq.len, word_length, n_chunks);

    // Pre-size flat per-partition slots with ~16% slack + 6*sqrt(mean) so normal genomes hit
    // 0% overflow while adversarial homopolymer/tandem runs spill safely to thread-local Vecs.
    let est_total = (valid_kmers.saturating_mul(2)).div_ceil(window + 1).max(1);
    let est_slab = if part_filter {
        (est_total.saturating_mul(n_parts))
            .div_ceil(NUM_PARTITIONS)
            .max(1)
    } else {
        est_total
    };
    let mean = est_slab.div_ceil(n_parts.max(1));
    let sqrt_mean = (mean as f64).sqrt() as usize;
    let slack = (mean / 6).max(6usize.saturating_mul(sqrt_mean)).max(32);
    let slot_cap = mean.saturating_add(slack).next_multiple_of(64);
    let alloc_total = n_parts
        .checked_mul(slot_cap)
        .filter(|&n| n <= u32::MAX as usize)
        .ok_or_else(|| {
            GeckoError::InvalidParameter(format!(
                "minimizer slab allocation ({n_parts} * {slot_cap}) exceeds 32-bit CSR limit"
            ))
        })?;

    let mut tags: Vec<MaybeUninit<u32>> = Vec::with_capacity(alloc_total);
    unsafe { tags.set_len(alloc_total) };
    let mut pos_strand: Vec<MaybeUninit<P>> = Vec::with_capacity(alloc_total);
    unsafe { pos_strand.set_len(alloc_total) };
    let mut coarse_mid: Vec<u8> = vec![0u8; alloc_total * 3 + 1];
    let part_cursors: Vec<PaddedAtomicUsize> = (0..n_parts)
        .map(|_| PaddedAtomicUsize(AtomicUsize::new(0)))
        .collect();

    let tptr = SendMutPtr(tags.as_mut_ptr() as *mut u32);
    let mptr = SendMutPtr(coarse_mid.as_mut_ptr());
    let pptr = SendMutPtr(pos_strand.as_mut_ptr() as *mut P);

    #[cold]
    #[inline(never)]
    fn push_min_overflow<P: MinPackedPs>(
        ov: &mut Vec<MinOverflowEntry<P>>,
        lp: usize,
        tags_s: &[u32],
        mids_s: &[u8],
        ps_s: &[P],
    ) {
        let lp_u8 = lp as u8;
        for i in 0..tags_s.len() {
            ov.push(MinOverflowEntry {
                lp: lp_u8,
                tag: tags_s[i],
                mid: [mids_s[i * 3], mids_s[i * 3 + 1], mids_s[i * 3 + 2]],
                ps: ps_s[i],
            });
        }
    }

    // Pass 1: Parallel single-pass atomic split-tag scatter (zero serial transpose!).
    let chunk_ovs: Vec<Vec<MinOverflowEntry<P>>> = chunks
        .par_iter()
        .map_init(
            || {
                let buf_t = vec![[0u32; STAGE_CAP]; n_parts];
                let buf_m = vec![[0u8; STAGE_MID_BYTES]; n_parts];
                let buf_p = vec![[P::default(); STAGE_CAP]; n_parts];
                (buf_t, buf_m, buf_p)
            },
            |(buf_t, buf_m, buf_p), &(cs, ce)| {
                let (tp, mp, pp) = (tptr, mptr, pptr);
                let mut buf_len = [0u8; NUM_PARTITIONS];
                let mut local_ov: Vec<MinOverflowEntry<P>> = Vec::new();

                let flush_part =
                    |lp: usize,
                     len: usize,
                     bt: &[u32; STAGE_CAP],
                     bm: &[u8; STAGE_MID_BYTES],
                     bp: &[P; STAGE_CAP],
                     ov: &mut Vec<MinOverflowEntry<P>>| {
                        let local_dst = unsafe { part_cursors.get_unchecked(lp) }
                            .0
                            .fetch_add(len, Ordering::Relaxed);
                        let flat_avail = slot_cap.saturating_sub(local_dst);
                        let n_flat = len.min(flat_avail);
                        if n_flat > 0 {
                            let dst = lp * slot_cap + local_dst;
                            unsafe {
                                std::ptr::copy_nonoverlapping(bt.as_ptr(), tp.0.add(dst), n_flat);
                                std::ptr::copy_nonoverlapping(
                                    bm.as_ptr(),
                                    mp.0.add(dst * 3),
                                    n_flat * 3,
                                );
                                std::ptr::copy_nonoverlapping(bp.as_ptr(), pp.0.add(dst), n_flat);
                            }
                        }
                        if n_flat < len {
                            push_min_overflow(
                                ov,
                                lp,
                                &bt[n_flat..len],
                                &bm[n_flat * 3..len * 3],
                                &bp[n_flat..len],
                            );
                        }
                    };

                let mut emit =
                    |pos: u64, canon: u64, _score: u32, kind: u8, _rs: usize, _re: usize| {
                        let ih = index_hash(canon);
                        let p = (ih >> 56) as u32;
                        if part_filter && (p < part_start || p >= part_end) {
                            return;
                        }
                        let lp = (p - part_start) as usize;
                        let ps = P::pack(pos, u8::from(kind == 1));
                        unsafe {
                            let idx = *buf_len.get_unchecked(lp) as usize;
                            let bt = buf_t.get_unchecked_mut(lp);
                            let bm = buf_m.get_unchecked_mut(lp);
                            let bp = buf_p.get_unchecked_mut(lp);
                            *bt.get_unchecked_mut(idx) = ih as u32;
                            let mid = (ih >> 32) as u32;
                            (bm.as_mut_ptr().add(idx * 3) as *mut u32).write_unaligned(mid.to_le());
                            *bp.get_unchecked_mut(idx) = ps;
                            let next_len = idx + 1;
                            if next_len == STAGE_CAP {
                                flush_part(lp, STAGE_CAP, bt, bm, bp, &mut local_ov);
                                *buf_len.get_unchecked_mut(lp) = 0;
                            } else {
                                *buf_len.get_unchecked_mut(lp) = next_len as u8;
                            }
                        }
                    };

                if forward_only {
                    for_each_minimizer::<true, true, _>(
                        seq,
                        cs,
                        ce,
                        word_length,
                        kp,
                        window,
                        soft_mask,
                        &mut emit,
                    );
                } else {
                    for_each_minimizer::<true, false, _>(
                        seq,
                        cs,
                        ce,
                        word_length,
                        kp,
                        window,
                        soft_mask,
                        &mut emit,
                    );
                }

                for lp in 0..n_parts {
                    let rem = unsafe { *buf_len.get_unchecked(lp) } as usize;
                    if rem > 0 {
                        unsafe {
                            flush_part(
                                lp,
                                rem,
                                buf_t.get_unchecked(lp),
                                buf_m.get_unchecked(lp),
                                buf_p.get_unchecked(lp),
                                &mut local_ov,
                            );
                        }
                    }
                }
                local_ov
            },
        )
        .collect();

    let mut raw_entries = 0usize;
    let mut flat_lens = Vec::with_capacity(n_parts);
    for cur in &part_cursors {
        let total_p = cur.0.load(Ordering::Relaxed);
        raw_entries += total_p;
        flat_lens.push(total_p.min(slot_cap));
    }

    if raw_entries == 0 {
        return Ok(MinimizerIndex::empty(
            word_length,
            kp,
            window,
            freq_limit,
            soft_mask,
            forward_only,
            part_start,
            part_end,
            use_u32_ps,
        ));
    }

    let mut ov_by_part: Vec<Vec<MinOverflowEntry<P>>> = vec![Vec::new(); n_parts];
    for ov_vec in chunk_ovs {
        for e in ov_vec {
            ov_by_part[e.lp as usize].push(e);
        }
    }

    let bucket_bits = choose_bucket_bits(raw_entries);
    let sub_bits = bucket_bits - 8;
    let hash_shift = 64 - bucket_bits;
    let sub_buckets = 1usize << sub_bits;
    let sub_shift_hi = 24 - sub_bits;
    let sub_mask = sub_buckets - 1;

    let bloom_params = choose_bloom_params(raw_entries, n_parts, freq_limit);
    let bloom_sub_bits = bloom_params.map_or(0u32, |(sb, _)| sb);
    let bloom_k_h = bloom_params.map_or(0usize, |(_, kh)| kh as usize);
    let bloom_blocks_per_part = if bloom_params.is_some() {
        1usize << bloom_sub_bits
    } else {
        0
    };
    let bloom_sub_shift = 24 - bloom_sub_bits;
    let bloom_sub_mask = bloom_blocks_per_part.saturating_sub(1);

    let off_total = n_parts * (sub_buckets + 1);
    let mut offsets: Vec<MaybeUninit<u32>> = Vec::with_capacity(off_total);
    unsafe { offsets.set_len(off_total) };
    let mut bloom_blocks = if bloom_blocks_per_part > 0 {
        vec![BloomBlock([0u64; 8]); n_parts * bloom_blocks_per_part]
    } else {
        Vec::new()
    };

    let kp_mask = (1u64 << (2 * kp)) - 1;

    // Pass 2: Dispatch ONCE on compile-time `bloom_k_h` (0 indirect jumps!).
    macro_rules! dispatch_pass2 {
        ($kh:expr) => {
            run_min_pass2_const::<P, $kh>(
                seq,
                word_length,
                kp,
                window,
                kp_mask,
                soft_mask,
                forward_only,
                freq_limit,
                n_parts,
                slot_cap,
                sub_buckets,
                sub_shift_hi,
                sub_mask,
                bloom_blocks_per_part,
                bloom_sub_shift,
                bloom_sub_mask,
                &flat_lens,
                &ov_by_part,
                &coarse_mid,
                &mut tags,
                &mut pos_strand,
                &mut offsets,
                &mut bloom_blocks,
            )
        };
    }

    let part_outs = match bloom_k_h {
        0 => dispatch_pass2!(0),
        1 => dispatch_pass2!(1),
        2 => dispatch_pass2!(2),
        3 => dispatch_pass2!(3),
        4 => dispatch_pass2!(4),
        5 => dispatch_pass2!(5),
        6 => dispatch_pass2!(6),
        7 => dispatch_pass2!(7),
        _ => dispatch_pass2!(8),
    };

    // Release transient `coarse_mid` immediately (returns 3.45 B/entry via munmap on large slabs).
    drop(coarse_mid);
    drop(ov_by_part);

    // SAFETY: Every element of `tags`, `pos_strand`, and `offsets` has been initialized by Pass 2.
    let mut tags: Vec<u32> = unsafe {
        let ptr = tags.as_mut_ptr() as *mut u32;
        let len = tags.len();
        let cap = tags.capacity();
        std::mem::forget(tags);
        Vec::from_raw_parts(ptr, len, cap)
    };
    let mut pos_strand_raw: Vec<P> = unsafe {
        let ptr = pos_strand.as_mut_ptr() as *mut P;
        let len = pos_strand.len();
        let cap = pos_strand.capacity();
        std::mem::forget(pos_strand);
        Vec::from_raw_parts(ptr, len, cap)
    };
    let mut offsets: Vec<u32> = unsafe {
        let ptr = offsets.as_mut_ptr() as *mut u32;
        let len = offsets.len();
        let cap = offsets.capacity();
        std::mem::forget(offsets);
        Vec::from_raw_parts(ptr, len, cap)
    };

    let mut kept_entries = 0usize;
    let mut any_dropped = false;
    let mut total_ov_kmers = 0usize;
    for po in &part_outs {
        kept_entries += po.kept;
        any_dropped |= !po.dropped_r.is_empty();
        total_ov_kmers += po.ov_kmers.len();
    }

    let mut dropped_kmer_bm = if any_dropped {
        vec![0u64; seq.len.div_ceil(64)]
    } else {
        Vec::new()
    };
    let mut all_ov_kmers = Vec::with_capacity(total_ov_kmers);
    let off_stride = sub_buckets + 1;

    for po in part_outs {
        for r in po.dropped_r {
            dropped_kmer_bm[r >> 6] |= 1u64 << (r & 63);
        }
        all_ov_kmers.extend_from_slice(&po.ov_kmers);
        if let Some((lp, tmp_t, tmp_p)) = po.append {
            let base = tags.len();
            if base + tmp_t.len() > u32::MAX as usize {
                return Err(GeckoError::InvalidParameter(
                    "minimizer overflow append exceeds 32-bit CSR limit".into(),
                ));
            }
            let base_u32 = base as u32;
            tags.extend_from_slice(&tmp_t);
            pos_strand_raw.extend_from_slice(&tmp_p);
            for slot in &mut offsets[lp * off_stride..(lp + 1) * off_stride] {
                *slot += base_u32;
            }
        }
    }

    let overflow = OverflowKmerCsr::build(all_ov_kmers, use_u32_ps);

    let bloom = if let Some((sb, kh)) = bloom_params {
        Some(BlockedBloom {
            blocks: bloom_blocks,
            sub_bits: sb,
            block_shift: 56 - sb,
            k_h: kh,
        })
    } else {
        None
    };

    Ok(MinimizerIndex {
        word_length,
        kp,
        window,
        freq_limit,
        soft_mask,
        forward_only,
        part_start,
        part_end,
        sub_bits,
        hash_shift,
        offsets,
        tags,
        pos_strand: P::into_ps_vec(pos_strand_raw),
        bloom,
        overflow,
        dropped_kmer_bm,
        kept_entries,
        raw_entries,
    })
}

/// Compact canonical Hash-CSR + L3 Blocked Bloom index over shift-invariant reference minimizers,
/// with optional high-frequency exact-`k` `OverflowKmerCsr`.
#[derive(Debug, Clone)]
pub struct MinimizerIndex {
    pub word_length: usize,
    pub kp: usize,
    pub window: usize,
    pub freq_limit: u64,
    pub soft_mask: bool,
    pub forward_only: bool,
    pub part_start: u32,
    pub part_end: u32,
    pub sub_bits: u32,
    pub hash_shift: u32,
    offsets: Vec<u32>,
    tags: Vec<u32>,
    pos_strand: PosStrandVec,
    bloom: Option<BlockedBloom>,
    overflow: Option<OverflowKmerCsr>,
    /// 1-bit/bp bitmap over `reference`: bit `r == 1` iff the `k`-mer starting at forward
    /// reference coordinate `r` exceeded `freq_limit` and was dropped. Empty when no `k`-mer in
    /// the indexed slab(s) exceeded `freq_limit`.
    dropped_kmer_bm: Vec<u64>,
    pub kept_entries: usize,
    pub raw_entries: usize,
}

#[inline]
fn estimate_minimizer_entries(reference: &Genome, word_length: usize, soft_mask: bool) -> usize {
    let Some((_kp, w)) = minimizer_params(word_length) else {
        return 0;
    };
    let valid_kmers = count_valid_kmers(reference, word_length, soft_mask);
    if valid_kmers == 0 {
        0
    } else {
        (valid_kmers.saturating_mul(2)).div_ceil(w + 1).max(1)
    }
}

#[inline]
fn minimizer_bytes_per_entry(reference: &Genome) -> usize {
    if reference.seq.len < ((1usize << 31) - 1) {
        MIN_PEAK_BYTES_PER_ENTRY_U32
    } else {
        MIN_PEAK_BYTES_PER_ENTRY_U64
    }
}

impl MinimizerIndex {
    #[allow(clippy::too_many_arguments)]
    fn empty(
        word_length: usize,
        kp: usize,
        window: usize,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
        part_start: u32,
        part_end: u32,
        use_u32_ps: bool,
    ) -> Self {
        Self {
            word_length,
            kp,
            window,
            freq_limit,
            soft_mask,
            forward_only,
            part_start,
            part_end,
            sub_bits: 0,
            hash_shift: 56,
            offsets: vec![0, 0],
            tags: Vec::new(),
            pos_strand: if use_u32_ps {
                PosStrandVec::U32(Vec::new())
            } else {
                PosStrandVec::U64(Vec::new())
            },
            bloom: None,
            overflow: None,
            dropped_kmer_bm: Vec::new(),
            kept_entries: 0,
            raw_entries: 0,
        }
    }

    /// Approximate heap bytes held by this minimizer index (including overflow exact-`k` CSR).
    pub fn heap_bytes(&self) -> usize {
        self.offsets.capacity() * 4
            + self.tags.capacity() * 4
            + self.pos_strand.heap_bytes()
            + self.bloom.as_ref().map_or(0, BlockedBloom::heap_bytes)
            + self
                .overflow
                .as_ref()
                .map_or(0, OverflowKmerCsr::heap_bytes)
            + self.dropped_kmer_bm.capacity() * 8
    }

    /// Estimate whether a single-slab `MinimizerIndex` (plus `rev_ref`) fits inside `memory_budget`
    /// taking `soft_mask` and `'N'`-runs into account (`W3-F8c`).
    pub fn fits_in_budget_with_mask(
        reference: &Genome,
        word_length: usize,
        soft_mask: bool,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> bool {
        if minimizer_params(word_length).is_none() {
            return false;
        }
        let rev_bytes = if compare_reverse {
            reference.seq.heap_bytes()
        } else {
            0
        };
        let est_entries = estimate_minimizer_entries(reference, word_length, soft_mask);
        let peak_bytes = est_entries
            .saturating_mul(minimizer_bytes_per_entry(reference))
            .saturating_add(rev_bytes)
            .saturating_add(MIN_FIXED_OVERHEAD);
        peak_bytes <= memory_budget
    }

    /// Estimate whether a single-slab `MinimizerIndex` (plus `rev_ref`) fits inside `memory_budget`.
    /// Conservative unmasked wrapper over [`Self::fits_in_budget_with_mask`].
    pub fn fits_in_budget(
        reference: &Genome,
        word_length: usize,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> bool {
        Self::fits_in_budget_with_mask(
            reference,
            word_length,
            false,
            compare_reverse,
            memory_budget,
        )
    }

    /// Compute the maximum minimizer entries per slab for `memory_budget` (`W3-F8c`).
    pub fn max_entries_per_slab(
        reference: &Genome,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> usize {
        let rev_bytes = if compare_reverse {
            reference.seq.heap_bytes()
        } else {
            0
        };
        let avail = memory_budget
            .saturating_sub(rev_bytes)
            .saturating_sub(MIN_FIXED_OVERHEAD);
        (avail / minimizer_bytes_per_entry(reference)).max(512)
    }

    /// Build a full-reference (`0..256`) canonical minimizer index.
    pub fn build(
        reference: &Genome,
        word_length: usize,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
    ) -> Result<Self> {
        Self::build_partitions(
            reference,
            0,
            NUM_PARTITIONS as u32,
            word_length,
            freq_limit,
            soft_mask,
            forward_only,
        )
    }

    /// Build a canonical minimizer index covering coarse hash partitions `[part_start .. part_end)`.
    pub fn build_partitions(
        reference: &Genome,
        part_start: u32,
        part_end: u32,
        word_length: usize,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
    ) -> Result<Self> {
        let Some((kp, window)) = minimizer_params(word_length) else {
            return Err(GeckoError::InvalidParameter(format!(
                "minimizer seeding requires contiguous word_length in {MINIMIZER_MIN_K}..=32, got {word_length}"
            )));
        };
        if part_start >= part_end || part_end > NUM_PARTITIONS as u32 {
            return Err(GeckoError::InvalidParameter(format!(
                "invalid minimizer partition range [{part_start}..{part_end})"
            )));
        }
        let seq = &reference.seq;
        // Reserve `u32::MAX` strictly for `OVERFLOW_SENTINEL_U32`.
        let use_u32_ps = seq.len < ((1usize << 31) - 1);
        if seq.len < word_length || freq_limit == 0 {
            return Ok(Self::empty(
                word_length,
                kp,
                window,
                freq_limit,
                soft_mask,
                forward_only,
                part_start,
                part_end,
                use_u32_ps,
            ));
        }

        if use_u32_ps {
            build_minimizer_typed::<u32>(
                reference,
                part_start,
                part_end,
                word_length,
                kp,
                window,
                freq_limit,
                soft_mask,
                forward_only,
                true,
            )
        } else {
            build_minimizer_typed::<u64>(
                reference,
                part_start,
                part_end,
                word_length,
                kp,
                window,
                freq_limit,
                soft_mask,
                forward_only,
                false,
            )
        }
    }

    #[inline(always)]
    fn is_kmer_dropped(&self, fwd_r: usize) -> bool {
        if self.dropped_kmer_bm.is_empty() {
            return false;
        }
        ((unsafe { *self.dropped_kmer_bm.get_unchecked(fwd_r >> 6) }) >> (fwd_r & 63)) & 1 != 0
    }
}

/// Partition the 256 uniform coarse minimizer hash partitions into balanced equal-width slabs
/// without performing an extra reference scan (`W3-F8c`).
pub fn compute_minimizer_partition_slabs(
    reference: &Genome,
    word_length: usize,
    soft_mask: bool,
    max_entries_per_slab: usize,
    _forward_only: bool,
) -> Vec<(u32, u32)> {
    let est_entries = estimate_minimizer_entries(reference, word_length, soft_mask);
    if est_entries == 0 {
        return Vec::new();
    }
    let num_slabs = est_entries
        .div_ceil(max_entries_per_slab.max(1))
        .clamp(1, NUM_PARTITIONS);
    let mut slabs = Vec::with_capacity(num_slabs);
    for s in 0..num_slabs {
        let ps = (s * NUM_PARTITIONS / num_slabs) as u32;
        let pe = ((s + 1) * NUM_PARTITIONS / num_slabs) as u32;
        if ps < pe {
            slabs.push((ps, pe));
        }
    }
    slabs
}

// ---------------------------------------------------------------------------
// Minimizer streaming query engine (W3-F6 + W3-F9)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Default)]
struct MinQueryItem {
    pos_x: u64,
    run_start: usize,
    run_end: usize,
    score: u32,
    b_idx: u32,
    tag: u32,
    bloom_blk: u32,
    start: u32,
    end: u32,
    q_kind: u8,
}

#[derive(Clone, Copy)]
enum MinIsa {
    #[cfg(target_arch = "x86_64")]
    Avx2Bmi2Popcnt,
    #[cfg(target_arch = "x86_64")]
    Avx2Popcnt,
    Portable,
}

struct MinWorkerState<C: DiagPairOps> {
    fwd_cache: StrandDiagCache<C>,
    rev_cache: Option<StrandDiagCache<C>>,
    active_ring: ActiveDiagRing,
    skip_scratch: SkipScratch,
    last_chunk_end: usize,
    sp_x_cursor: usize,
    fragments: Vec<Fragment>,
    coverage: Option<CoverageMatrix>,
    word_matches: u64,
    raw_hits: u64,
    hits_used: u64,
    query_minimizers: u64,
    bloom_rejected: u64,
    exact_k_rejected: u64,
    fast_rejects: u64,
}

impl<C: DiagPairOps> MinWorkerState<C> {
    fn new(compare_reverse: bool, collect_coverage: bool) -> Self {
        Self {
            fwd_cache: StrandDiagCache::new(),
            rev_cache: if compare_reverse {
                Some(StrandDiagCache::new())
            } else {
                None
            },
            active_ring: ActiveDiagRing::new(),
            skip_scratch: SkipScratch::default(),
            last_chunk_end: 0,
            sp_x_cursor: 0,
            fragments: Vec::new(),
            coverage: if collect_coverage {
                Some(CoverageMatrix::new())
            } else {
                None
            },
            word_matches: 0,
            raw_hits: 0,
            hits_used: 0,
            query_minimizers: 0,
            bloom_rejected: 0,
            exact_k_rejected: 0,
            fast_rejects: 0,
        }
    }

    #[inline(always)]
    fn sp_x_cursor_for(&mut self, sp_x: &[(u64, u64)], pos_x: u64) -> usize {
        let mut c = self.sp_x_cursor;
        while c > 0 && unsafe { sp_x.get_unchecked(c - 1).1 } > pos_x {
            c -= 1;
        }
        while c < sp_x.len() && unsafe { sp_x.get_unchecked(c).1 } <= pos_x {
            c += 1;
        }
        self.sp_x_cursor = c;
        c
    }
}

#[inline(always)]
fn is_diag_covered<C: DiagPairOps>(
    state: &MinWorkerState<C>,
    is_rev: bool,
    diag: i64,
    pos_x: u64,
) -> bool {
    if state.active_ring.covers(diag, pos_x, is_rev as u8) {
        return true;
    }
    if !is_rev {
        state.fwd_cache.is_covered(diag, pos_x)
    } else {
        state.rev_cache.as_ref().unwrap().is_covered(diag, pos_x)
    }
}

#[inline(always)]
fn unpack_ps(ps: &PosStrandVec, idx: usize) -> (u64, u8) {
    match ps {
        PosStrandVec::U32(v) => {
            let raw = unsafe { *v.get_unchecked(idx) };
            ((raw >> 1) as u64, (raw & 1) as u8)
        }
        PosStrandVec::U64(v) => {
            let raw = unsafe { *v.get_unchecked(idx) };
            (raw >> 1, (raw & 1) as u8)
        }
    }
}

#[inline(always)]
fn unpack_ps_or_sentinel(ps: &PosStrandVec, idx: usize) -> (u64, u8) {
    match ps {
        PosStrandVec::U32(v) => {
            let raw = unsafe { *v.get_unchecked(idx) };
            if raw == OVERFLOW_SENTINEL_U32 {
                (OVERFLOW_SENTINEL_POS, 0)
            } else {
                ((raw >> 1) as u64, (raw & 1) as u8)
            }
        }
        PosStrandVec::U64(v) => {
            let raw = unsafe { *v.get_unchecked(idx) };
            if raw == OVERFLOW_SENTINEL_U64 {
                (OVERFLOW_SENTINEL_POS, 0)
            } else {
                (raw >> 1, (raw & 1) as u8)
            }
        }
    }
}

#[inline(always)]
fn prefetch_ps(ps: &PosStrandVec, idx: usize) {
    match ps {
        PosStrandVec::U32(v) => prefetch_t2(unsafe { v.as_ptr().add(idx) }),
        PosStrandVec::U64(v) => prefetch_t2(unsafe { v.as_ptr().add(idx) }),
    }
}

/// Const-specialized Stage-1 Blocked Bloom batch filter (0 indirect jumps!).
#[inline(always)]
unsafe fn filter_min_bloom_batch_const<const K_H: usize>(
    bloom: &BlockedBloom,
    batch: &[MinQueryItem],
    off_ptr: *const u32,
    survivors: &mut [u8; MIN_BATCH_SIZE],
) -> usize {
    let mut ns = 0usize;
    for (i, item) in batch.iter().enumerate() {
        if unsafe { bloom.contains_unchecked_const::<K_H>(item.bloom_blk as usize, item.tag) } {
            prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
            unsafe { *survivors.get_unchecked_mut(ns) = i as u8 };
            ns += 1;
        }
    }
    ns
}

#[cold]
#[inline(never)]
unsafe fn filter_min_bloom_batch_rare(
    bloom: &BlockedBloom,
    batch: &[MinQueryItem],
    off_ptr: *const u32,
    survivors: &mut [u8; MIN_BATCH_SIZE],
) -> usize {
    match bloom.k_h {
        3 => unsafe { filter_min_bloom_batch_const::<3>(bloom, batch, off_ptr, survivors) },
        2 => unsafe { filter_min_bloom_batch_const::<2>(bloom, batch, off_ptr, survivors) },
        1 => unsafe { filter_min_bloom_batch_const::<1>(bloom, batch, off_ptr, survivors) },
        7 => unsafe { filter_min_bloom_batch_const::<7>(bloom, batch, off_ptr, survivors) },
        _ => unsafe { filter_min_bloom_batch_const::<8>(bloom, batch, off_ptr, survivors) },
    }
}

#[inline(always)]
unsafe fn filter_min_bloom_batch(
    bloom: &BlockedBloom,
    batch: &[MinQueryItem],
    off_ptr: *const u32,
    survivors: &mut [u8; MIN_BATCH_SIZE],
) -> usize {
    let kh = bloom.k_h;
    if kh == 4 {
        unsafe { filter_min_bloom_batch_const::<4>(bloom, batch, off_ptr, survivors) }
    } else if kh == 6 {
        unsafe { filter_min_bloom_batch_const::<6>(bloom, batch, off_ptr, survivors) }
    } else if kh == 5 {
        unsafe { filter_min_bloom_batch_const::<5>(bloom, batch, off_ptr, survivors) }
    } else {
        unsafe { filter_min_bloom_batch_rare(bloom, batch, off_ptr, survivors) }
    }
}

#[cold]
#[inline(never)]
unsafe fn ov_bloom_contains_rare(bloom: &BlockedBloom, blk: usize, tag: u32) -> bool {
    match bloom.k_h {
        1 => unsafe { bloom.contains_unchecked_const::<1>(blk, tag) },
        2 => unsafe { bloom.contains_unchecked_const::<2>(blk, tag) },
        3 => unsafe { bloom.contains_unchecked_const::<3>(blk, tag) },
        _ => unsafe { bloom.contains_unchecked_const::<8>(blk, tag) },
    }
}

#[inline(always)]
unsafe fn ov_bloom_contains(bloom: &BlockedBloom, blk: usize, tag: u32) -> bool {
    let kh = bloom.k_h;
    if kh == 7 {
        unsafe { bloom.contains_unchecked_const::<7>(blk, tag) }
    } else if kh == 6 {
        unsafe { bloom.contains_unchecked_const::<6>(blk, tag) }
    } else if kh == 5 {
        unsafe { bloom.contains_unchecked_const::<5>(blk, tag) }
    } else if kh == 4 {
        unsafe { bloom.contains_unchecked_const::<4>(blk, tag) }
    } else {
        unsafe { ov_bloom_contains_rare(bloom, blk, tag) }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
#[inline]
unsafe fn try_fast_reject_min_bmi2<const COLLECT_COV: bool>(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    seed_score: i64,
    min_len: u64,
) -> Option<FastRejectOutcome> {
    unsafe {
        try_fast_short_reject_auto::<COLLECT_COV, _>(
            seq_x,
            seq_y,
            pos_x,
            pos_y,
            wl,
            seed_score,
            min_len,
            |mm_even| pack_mm32_bmi2(mm_even),
        )
    }
}

/// Cold AVX2/portable slow-path extension keeping `vzeroupper` completely off the scalar
/// fast-reject hot path, with lossless `exceeded_min_len`, `exact_rejected_skip_checked`,
/// `ActiveDiagRing`, and `update_noise(outcome.xfil_max)` (`W3-F9`).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn extend_verified_seed_slow<const COLLECT_COV: bool, C: DiagPairOps>(
    is_rev: bool,
    diag: i64,
    seed_x: u64,
    seed_y: u64,
    wl: u64,
    seed_score: i64,
    query: &Genome,
    target: &Genome,
    opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState<C>,
) -> u64 {
    let cx_x = state.sp_x_cursor_for(&query.seq.special_runs, seed_x);
    let cov_opt = if COLLECT_COV {
        state.coverage.as_mut()
    } else {
        None
    };
    let outcome = match isa {
        #[cfg(target_arch = "x86_64")]
        MinIsa::Avx2Bmi2Popcnt => unsafe {
            extend_single_hit_with_seed_avx2_bmi2(
                Hit {
                    diag,
                    pos_x: seed_x,
                },
                query,
                target,
                opts,
                wl,
                seed_score,
                cx_x,
                cov_opt,
            )
        },
        #[cfg(target_arch = "x86_64")]
        MinIsa::Avx2Popcnt => unsafe {
            extend_single_hit_with_seed_avx2_swar(
                Hit {
                    diag,
                    pos_x: seed_x,
                },
                query,
                target,
                opts,
                wl,
                seed_score,
                cx_x,
                cov_opt,
            )
        },
        MinIsa::Portable => extend_single_hit_with_seed_portable(
            Hit {
                diag,
                pos_x: seed_x,
            },
            query,
            target,
            opts,
            wl,
            seed_score,
            cov_opt,
        ),
    };

    if let Some(frag) = outcome.fragment {
        if outcome.xfil_max > seed_x {
            if !is_rev {
                state.fwd_cache.update_hsp(diag, outcome.xfil_max);
            } else {
                state
                    .rev_cache
                    .as_mut()
                    .unwrap()
                    .update_hsp(diag, outcome.xfil_max);
            }
            state
                .active_ring
                .record(diag, outcome.xfil_max, is_rev as u8);
        }
        state.fragments.push(frag);
        outcome.xfil_max
    } else if outcome.exceeded_min_len || outcome.xfil_max > seed_x + 256 {
        let skip_end = exact_rejected_skip(
            &query.seq,
            &target.seq,
            seed_x,
            seed_y,
            outcome.xfil_max,
            opts.word_length,
            opts.scoring,
            &mut state.skip_scratch,
        );
        // Critical guard: when `skip_end == seed_x` (valley `safe_b == 0` or `SEQ_SEP`),
        // never poison cache or active_ring!
        if skip_end > seed_x {
            if !is_rev {
                state.fwd_cache.update_hsp(diag, skip_end);
            } else {
                state.rev_cache.as_mut().unwrap().update_hsp(diag, skip_end);
            }
            state.active_ring.record(diag, skip_end, is_rev as u8);
        }
        skip_end
    } else {
        // Lossless short-reject noise bound (`outcome.xfil_max` never crosses `SEQ_SEP`).
        let skip = outcome.xfil_max;
        if skip > seed_x {
            if !is_rev {
                state.fwd_cache.update_noise(diag, skip);
            } else {
                state.rev_cache.as_mut().unwrap().update_noise(diag, skip);
            }
        }
        skip
    }
}

/// Scalar 64/128-bp short-death fast-reject followed by cold AVX2/portable extension.
/// Returns `advance_to`: inclusive query `x` coordinate covered by this seed's outcome.
#[inline]
#[allow(clippy::too_many_arguments)]
fn extend_verified_seed<const COLLECT_COV: bool, C: DiagPairOps>(
    is_rev: bool,
    diag: i64,
    seed_x: u64,
    seed_y: u64,
    query: &Genome,
    target: &Genome,
    opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState<C>,
) -> u64 {
    state.hits_used += 1;
    let k = opts.word_length;
    let wl = k as u64;
    let seed_score = k as i64 * opts.scoring.match_score;

    if opts.fixed_length && opts.scoring.is_default_gecko() {
        let fast_rej = match isa {
            #[cfg(target_arch = "x86_64")]
            MinIsa::Avx2Bmi2Popcnt => unsafe {
                try_fast_reject_min_bmi2::<COLLECT_COV>(
                    &query.seq,
                    &target.seq,
                    seed_x,
                    seed_y,
                    wl,
                    seed_score,
                    opts.min_length,
                )
            },
            _ => unsafe {
                try_fast_short_reject_auto::<COLLECT_COV, _>(
                    &query.seq,
                    &target.seq,
                    seed_x,
                    seed_y,
                    wl,
                    seed_score,
                    opts.min_length,
                    pack_mm32_swar,
                )
            },
        };
        if let Some(rej) = fast_rej {
            if COLLECT_COV && let Some(cov) = state.coverage.as_mut() {
                cov.record_unit_exact(rej.length, rej.u_max);
            }
            state.fast_rejects += 1;
            let skip = seed_x + wl - 1 + rej.fwd_best as u64;
            if !is_rev {
                state.fwd_cache.update_noise(diag, skip);
            } else {
                state.rev_cache.as_mut().unwrap().update_noise(diag, skip);
            }
            return skip;
        }
    }

    extend_verified_seed_slow::<COLLECT_COV, C>(
        is_rev, diag, seed_x, seed_y, wl, seed_score, query, target, opts, isa, state,
    )
}

/// Compare `query` against a pre-built [`MinimizerIndex`] with optional `.MAT` coverage matrix
/// collection, achieving 100% bit-exact `.frags` parity with stride-1 `RefIndex`.
pub fn compare_minimizer_streaming_with_coverage(
    query: &Genome,
    reference: &Genome,
    min_index: &MinimizerIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
    collect_coverage: bool,
) -> Result<CompareResult> {
    config.validate()?;
    let compare_reverse = config.compare_reverse && rev_ref.is_some();
    if min_index.forward_only && compare_reverse {
        return Err(GeckoError::InvalidParameter(
            "forward-only MinimizerIndex cannot be queried with compare_reverse=true".into(),
        ));
    }

    let k = min_index.word_length;
    if query.seq.len < k || min_index.kept_entries == 0 {
        return Ok(CompareResult {
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            stats: CompareStats::default(),
        });
    }

    let use_compact =
        query.seq.len <= COMPACT_DIAG_MAX_LEN && reference.seq.len <= COMPACT_DIAG_MAX_LEN;
    match (use_compact, collect_coverage) {
        (true, true) => compare_minimizer_streaming_typed::<true, DiagPair32>(
            query,
            reference,
            min_index,
            rev_ref,
            config,
            compare_reverse,
        ),
        (true, false) => compare_minimizer_streaming_typed::<false, DiagPair32>(
            query,
            reference,
            min_index,
            rev_ref,
            config,
            compare_reverse,
        ),
        (false, true) => compare_minimizer_streaming_typed::<true, DiagPair64>(
            query,
            reference,
            min_index,
            rev_ref,
            config,
            compare_reverse,
        ),
        (false, false) => compare_minimizer_streaming_typed::<false, DiagPair64>(
            query,
            reference,
            min_index,
            rev_ref,
            config,
            compare_reverse,
        ),
    }
}

fn compare_minimizer_streaming_typed<const COLLECT_COV: bool, C: DiagPairOps>(
    query: &Genome,
    reference: &Genome,
    min_index: &MinimizerIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
    compare_reverse: bool,
) -> Result<CompareResult> {
    let k = min_index.word_length;
    let t_query = Instant::now();
    let isa = {
        #[cfg(target_arch = "x86_64")]
        {
            if has_avx2() && has_popcnt() {
                if has_bmi2() {
                    MinIsa::Avx2Bmi2Popcnt
                } else {
                    MinIsa::Avx2Popcnt
                }
            } else {
                MinIsa::Portable
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (has_avx2(), has_bmi2(), has_popcnt());
            MinIsa::Portable
        }
    };

    let fwd_opts = ExtendOptions {
        min_length: config.min_length,
        min_similarity: config.min_similarity,
        word_length: k,
        fixed_length: config.fixed_length,
        strand: 'f',
        scoring: config.scoring,
        spaced_seed: false,
    };
    let rev_opts = ExtendOptions {
        strand: 'r',
        ..fwd_opts
    };

    let max_starts = query.seq.len - k + 1;
    let mut chunks = Vec::with_capacity(max_starts.div_ceil(MIN_QUERY_CHUNK_BYTES));
    let mut s = 0usize;
    while s < max_starts {
        let e = (s + MIN_QUERY_CHUNK_BYTES).min(max_starts);
        chunks.push((s, e));
        s = e;
    }

    let worker_results: Vec<MinWorkerState<C>> = chunks
        .into_par_iter()
        .fold(
            || MinWorkerState::<C>::new(compare_reverse, COLLECT_COV),
            |mut state, (chunk_start, chunk_end)| {
                if chunk_start < state.last_chunk_end {
                    state.fwd_cache.clear();
                    if let Some(rc) = state.rev_cache.as_mut() {
                        rc.clear();
                    }
                    state.active_ring.clear();
                    state.sp_x_cursor = if query.seq.special_runs.is_empty() {
                        0
                    } else {
                        query
                            .seq
                            .special_runs
                            .partition_point(|&(_, e)| e <= chunk_start as u64)
                    };
                } else if state.last_chunk_end == 0 && chunk_start > 0 {
                    state.sp_x_cursor = if query.seq.special_runs.is_empty() {
                        0
                    } else {
                        query
                            .seq
                            .special_runs
                            .partition_point(|&(_, e)| e <= chunk_start as u64)
                    };
                }
                scan_minimizer_chunk::<COLLECT_COV, C>(
                    query,
                    reference,
                    rev_ref,
                    min_index,
                    chunk_start,
                    chunk_end,
                    config.max_rep,
                    fwd_opts,
                    rev_opts,
                    isa,
                    &mut state,
                );
                state.last_chunk_end = chunk_end;
                state
            },
        )
        .collect();

    let query_secs = t_query.elapsed().as_secs_f64();
    let t_dedup = Instant::now();

    let mut fragments = Vec::new();
    let mut coverage = CoverageMatrix::new();
    let mut stats = CompareStats::default();
    let mut total_q_mins = 0u64;
    let mut total_bloom_rej = 0u64;
    let mut total_exact_rej = 0u64;
    let mut total_fast_rej = 0u64;

    for mut w in worker_results {
        stats.word_matches += w.word_matches;
        stats.raw_hits += w.raw_hits;
        stats.hits_used += w.hits_used;
        total_q_mins += w.query_minimizers;
        total_bloom_rej += w.bloom_rejected;
        total_exact_rej += w.exact_k_rejected;
        total_fast_rej += w.fast_rejects;
        if let Some(cov) = w.coverage {
            coverage.merge_from(&cov);
        }
        fragments.append(&mut w.fragments);
    }

    dedup_fragments(&mut fragments);
    for f in &fragments {
        if f.strand == 'f' {
            stats.forward_fragments += 1;
        } else {
            stats.reverse_fragments += 1;
        }
    }
    stats.filtered_hits = stats.hits_used;

    let dedup_secs = t_dedup.elapsed().as_secs_f64();
    if profile_enabled() {
        eprintln!(
            "[GECKO_PROFILE] minimizer_query k={} kp={} w={} parts=[{}..{}) ref_entries={} ov_entries={} q_mins={} bloom_reject={} exact_k_reject={} fast_reject={} hits_used={} frags={} query_ms={:.2} dedup_ms={:.2}",
            min_index.word_length,
            min_index.kp,
            min_index.window,
            min_index.part_start,
            min_index.part_end,
            min_index.kept_entries,
            min_index.overflow.as_ref().map_or(0, |o| o.kept_entries),
            total_q_mins,
            total_bloom_rej,
            total_exact_rej,
            total_fast_rej,
            stats.hits_used,
            fragments.len(),
            query_secs * 1e3,
            dedup_secs * 1e3,
        );
    }

    Ok(CompareResult {
        fragments,
        coverage,
        stats,
    })
}

#[allow(clippy::too_many_arguments)]
fn scan_minimizer_chunk<const COLLECT_COV: bool, C: DiagPairOps>(
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    min_index: &MinimizerIndex,
    chunk_start: usize,
    chunk_end: usize,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState<C>,
) {
    let k = min_index.word_length;
    let kp = min_index.kp;
    let w = min_index.window;
    let part_start = min_index.part_start;
    let part_end = min_index.part_end;
    let part_filter = part_start > 0 || part_end < NUM_PARTITIONS as u32;
    let sub_bits = min_index.sub_bits;
    let hash_shift = min_index.hash_shift;
    let part_bucket_base = part_start << sub_bits;
    let bloom_shift = min_index.bloom.as_ref().map_or(56u32, |b| b.block_shift);
    let bloom_base = min_index
        .bloom
        .as_ref()
        .map_or(0u32, |b| part_start << b.sub_bits);
    let scan_end = (chunk_end + k - 1).min(query.seq.len);

    let mut batch = [MinQueryItem::default(); MIN_BATCH_SIZE];
    let mut batch_len = 0usize;

    let mut handle = |pos_x: u64, canon: u64, score: u32, q_kind: u8, rs: usize, re: usize| {
        state.query_minimizers += 1;
        let ih = index_hash(canon);
        if part_filter {
            let p = (ih >> 56) as u32;
            if p < part_start || p >= part_end {
                return;
            }
        }
        let b_global = (ih >> hash_shift) as u32;
        let b_bucket = if part_filter {
            b_global - part_bucket_base
        } else {
            b_global
        };
        let b_idx = b_bucket + (b_bucket >> sub_bits);
        let bloom_blk = if part_filter {
            ((ih >> bloom_shift) as u32) - bloom_base
        } else {
            (ih >> bloom_shift) as u32
        };
        batch[batch_len] = MinQueryItem {
            pos_x,
            run_start: rs,
            run_end: re,
            score,
            b_idx,
            tag: ih as u32,
            bloom_blk,
            start: 0,
            end: 0,
            q_kind,
        };
        batch_len += 1;
        if batch_len == MIN_BATCH_SIZE {
            flush_minimizer_batch::<COLLECT_COV, C>(
                &mut batch[..MIN_BATCH_SIZE],
                query,
                reference,
                rev_ref,
                min_index,
                max_rep,
                fwd_opts,
                rev_opts,
                isa,
                state,
            );
            batch_len = 0;
        }
    };

    if min_index.forward_only {
        for_each_minimizer::<false, true, _>(
            &query.seq,
            chunk_start,
            scan_end,
            k,
            kp,
            w,
            min_index.soft_mask,
            &mut handle,
        );
    } else {
        for_each_minimizer::<false, false, _>(
            &query.seq,
            chunk_start,
            scan_end,
            k,
            kp,
            w,
            min_index.soft_mask,
            &mut handle,
        );
    }

    if batch_len > 0 {
        flush_minimizer_batch::<COLLECT_COV, C>(
            &mut batch[..batch_len],
            query,
            reference,
            rev_ref,
            min_index,
            max_rep,
            fwd_opts,
            rev_opts,
            isa,
            state,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn flush_minimizer_batch<const COLLECT_COV: bool, C: DiagPairOps>(
    batch: &mut [MinQueryItem],
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    min_index: &MinimizerIndex,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState<C>,
) {
    let off_ptr = min_index.offsets.as_ptr();
    let tags = &min_index.tags;
    let tags_ptr = tags.as_ptr();

    let mut survivors = [0u8; MIN_BATCH_SIZE];
    let n_surv: usize;

    if let Some(ref bloom) = min_index.bloom {
        let bloom_ptr = bloom.blocks.as_ptr();
        for item in batch.iter() {
            prefetch_t1(unsafe { bloom_ptr.add(item.bloom_blk as usize) });
        }
        let ns = unsafe { filter_min_bloom_batch(bloom, batch, off_ptr, &mut survivors) };
        state.bloom_rejected += (batch.len() - ns) as u64;
        if ns == 0 {
            return;
        }
        n_surv = ns;
    } else {
        for (i, item) in batch.iter().enumerate() {
            prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
            survivors[i] = i as u8;
        }
        n_surv = batch.len();
    }

    for &idx in &survivors[..n_surv] {
        let item = &mut batch[idx as usize];
        let b_idx = item.b_idx as usize;
        let s = unsafe { *off_ptr.add(b_idx) };
        let e = unsafe { *off_ptr.add(b_idx + 1) };
        item.start = s;
        item.end = e;
        if s < e {
            prefetch_t2(unsafe { tags_ptr.add(s as usize) });
        }
    }

    let mut hits = [0u8; MIN_BATCH_SIZE];
    let mut n_hits = 0usize;
    for &idx in &survivors[..n_surv] {
        let item = &mut batch[idx as usize];
        let s = item.start as usize;
        let e = item.end as usize;
        if s == e {
            continue;
        }
        let (lo, hi) = if e - s <= 16 {
            probe_sorted_tags_scalar(tags, s, e, item.tag)
        } else {
            probe_bucket_binary_u32(tags, s, e, item.tag)
        };
        item.start = lo as u32;
        item.end = hi as u32;
        if lo < hi {
            prefetch_ps(&min_index.pos_strand, lo);
            hits[n_hits] = idx;
            n_hits += 1;
        }
    }
    if n_hits == 0 {
        return;
    }

    let k = min_index.word_length;
    let kp = min_index.kp;
    let w = min_index.window;
    let kp_mask = (1u64 << (2 * kp)) - 1;
    let is_multi_slab = min_index.part_start > 0 || min_index.part_end < NUM_PARTITIONS as u32;
    let fwd_only = min_index.forward_only;

    for &hit_idx in &hits[..n_hits] {
        let item = batch[hit_idx as usize];
        let lo = item.start as usize;
        let hi = item.end as usize;
        let total = hi - lo;
        let step = if total > max_rep {
            (total / max_rep).max(1)
        } else {
            1
        };

        let mut hit_overflow = false;
        if step > 1 && min_index.overflow.is_some() {
            // Ensure `max_rep` subsampling can never jump over an overflow sentinel.
            for scan_idx in lo..hi {
                let (p_r, _) = unpack_ps_or_sentinel(&min_index.pos_strand, scan_idx);
                if p_r == OVERFLOW_SENTINEL_POS {
                    hit_overflow = true;
                    break;
                }
            }
        }

        let mut idx = lo;
        let mut matched_any = false;
        while idx < hi {
            let (p_r, r_strand) = unpack_ps_or_sentinel(&min_index.pos_strand, idx);
            idx += step;

            if p_r == OVERFLOW_SENTINEL_POS {
                hit_overflow = true;
                continue;
            }

            let (try_fwd, try_rev) = if fwd_only {
                (true, false)
            } else if item.q_kind == 2 {
                (true, rev_ref.is_some())
            } else if r_strand == item.q_kind {
                (true, false)
            } else {
                (false, rev_ref.is_some())
            };

            if try_fwd {
                matched_any = true;
                state.raw_hits += 1;
                verify_and_extend_minimizer_hit::<COLLECT_COV, C>(
                    false,
                    item.pos_x,
                    p_r,
                    item.score,
                    item.run_start,
                    item.run_end,
                    query,
                    reference,
                    reference,
                    min_index,
                    k,
                    kp,
                    w,
                    kp_mask,
                    is_multi_slab,
                    fwd_opts,
                    isa,
                    state,
                );
            }
            if try_rev && let Some(rr) = rev_ref {
                matched_any = true;
                state.raw_hits += 1;
                verify_and_extend_minimizer_hit::<COLLECT_COV, C>(
                    true,
                    item.pos_x,
                    p_r,
                    item.score,
                    item.run_start,
                    item.run_end,
                    query,
                    reference,
                    rr,
                    min_index,
                    k,
                    kp,
                    w,
                    kp_mask,
                    is_multi_slab,
                    rev_opts,
                    isa,
                    state,
                );
            }
        }

        if hit_overflow
            && let Some(ref ov) = min_index.overflow
            && probe_overflow_minimizer_window::<COLLECT_COV, C>(
                &item, ov, query, reference, rev_ref, k, w, fwd_only, max_rep, fwd_opts, rev_opts,
                isa, state,
            )
        {
            matched_any = true;
        }

        if matched_any {
            state.word_matches += 1;
        }
    }
}

/// Probe `OverflowKmerCsr` for all valid query exact `k`-mers `[q_lo..=q_hi]` containing the
/// overflowed minimizer at `item.pos_x` (`W3-F6`). Visits `sx` in strictly increasing order so
/// diagonal coverage and valley-trimming parity with stride-1 `Exact` are 100% preserved.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn probe_overflow_minimizer_window<const COLLECT_COV: bool, C: DiagPairOps>(
    item: &MinQueryItem,
    ov: &OverflowKmerCsr,
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    k: usize,
    w: usize,
    fwd_only: bool,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState<C>,
) -> bool {
    let q_mpos = item.pos_x as usize;
    let q_lo = q_mpos.saturating_sub(w - 1).max(item.run_start);
    let q_hi = q_mpos.min(item.run_end - k);
    let n_win = q_hi - q_lo + 1;
    debug_assert!(n_win <= MAX_W);
    let k_mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let ny = reference.seq.len as u64;
    let hash_shift = ov.hash_shift;
    let mut matched_any = false;

    // 3-stage prefetched stack batch over the <= 12 k-mers in `[q_lo..=q_hi]` so L3/DRAM misses
    // to `ov.bloom`, `ov.offsets`, and `ov.tags` overlap in parallel instead of stalling serially.
    let mut q_fwds = [0u64; MAX_W];
    let mut q_kinds = [0u8; MAX_W];
    let mut khs = [0u64; MAX_W];
    let mut surv = [0u8; MAX_W];
    let mut n_surv = 0usize;

    if let Some(ref bloom) = ov.bloom {
        let bloom_ptr = bloom.blocks.as_ptr();
        let b_shift = bloom.block_shift;
        for j in 0..n_win {
            let sx = q_lo + j;
            let q_fwd = query.seq.load_32(sx) & k_mask;
            let (q_canon, q_k_kind) = if fwd_only {
                (q_fwd, 0u8)
            } else {
                let q_rc = reverse_2bit_comp(q_fwd, k);
                if q_fwd < q_rc {
                    (q_fwd, 0u8)
                } else if q_rc < q_fwd {
                    (q_rc, 1u8)
                } else {
                    (q_fwd, 2u8)
                }
            };
            let kh = fmix64(q_canon);
            unsafe {
                *q_fwds.get_unchecked_mut(j) = q_fwd;
                *q_kinds.get_unchecked_mut(j) = q_k_kind;
                *khs.get_unchecked_mut(j) = kh;
                prefetch_t1(bloom_ptr.add((kh >> b_shift) as usize));
            }
        }
        let off_ptr = ov.offsets.as_ptr();
        for j in 0..n_win {
            let kh = unsafe { *khs.get_unchecked(j) };
            let blk = (kh >> b_shift) as usize;
            let tag = kh as u32;
            if unsafe { ov_bloom_contains(bloom, blk, tag) } {
                let b = (kh >> hash_shift) as usize;
                prefetch_t2(unsafe { off_ptr.add(b) });
                unsafe { *surv.get_unchecked_mut(n_surv) = j as u8 };
                n_surv += 1;
            }
        }
    } else {
        let off_ptr = ov.offsets.as_ptr();
        for j in 0..n_win {
            let sx = q_lo + j;
            let q_fwd = query.seq.load_32(sx) & k_mask;
            let (q_canon, q_k_kind) = if fwd_only {
                (q_fwd, 0u8)
            } else {
                let q_rc = reverse_2bit_comp(q_fwd, k);
                if q_fwd < q_rc {
                    (q_fwd, 0u8)
                } else if q_rc < q_fwd {
                    (q_rc, 1u8)
                } else {
                    (q_fwd, 2u8)
                }
            };
            let kh = fmix64(q_canon);
            let b = (kh >> hash_shift) as usize;
            unsafe {
                *q_fwds.get_unchecked_mut(j) = q_fwd;
                *q_kinds.get_unchecked_mut(j) = q_k_kind;
                *khs.get_unchecked_mut(j) = kh;
                prefetch_t2(off_ptr.add(b));
                *surv.get_unchecked_mut(j) = j as u8;
            }
        }
        n_surv = n_win;
    }

    if n_surv == 0 {
        return false;
    }

    let mut starts = [0u32; MAX_W];
    let mut ends = [0u32; MAX_W];
    let tags_ptr = ov.tags.as_ptr();
    for &j_u8 in &surv[..n_surv] {
        let j = j_u8 as usize;
        let b = (unsafe { *khs.get_unchecked(j) } >> hash_shift) as usize;
        let s = unsafe { *ov.offsets.get_unchecked(b) };
        let e = unsafe { *ov.offsets.get_unchecked(b + 1) };
        unsafe {
            *starts.get_unchecked_mut(j) = s;
            *ends.get_unchecked_mut(j) = e;
            if s < e {
                prefetch_t2(tags_ptr.add(s as usize));
            }
        }
    }

    for &j_u8 in &surv[..n_surv] {
        let j = j_u8 as usize;
        let s = unsafe { *starts.get_unchecked(j) } as usize;
        let e = unsafe { *ends.get_unchecked(j) } as usize;
        if s == e {
            continue;
        }
        let tag = unsafe { *khs.get_unchecked(j) } as u32;
        let (lo, hi) = if e - s <= 16 {
            probe_sorted_tags_scalar(&ov.tags, s, e, tag)
        } else {
            probe_bucket_binary_u32(&ov.tags, s, e, tag)
        };
        if lo >= hi {
            continue;
        }

        let sx = q_lo + j;
        let sx64 = sx as u64;
        let q_fwd = unsafe { *q_fwds.get_unchecked(j) };
        let q_k_kind = unsafe { *q_kinds.get_unchecked(j) };
        let total = hi - lo;
        let step = if total > max_rep {
            (total / max_rep).max(1)
        } else {
            1
        };

        let mut idx = lo;
        while idx < hi {
            let (fwd_r, r_strand) = unpack_ps(&ov.pos_strand, idx);
            idx += step;

            let (try_fwd, try_rev) = if fwd_only {
                (true, false)
            } else if q_k_kind == 2 {
                (true, rev_ref.is_some())
            } else if r_strand == q_k_kind {
                (true, false)
            } else {
                (false, rev_ref.is_some())
            };

            if try_fwd {
                let diag = sx as i64 - fwd_r as i64;
                if !is_diag_covered(state, false, diag, sx64) {
                    if ((q_fwd ^ reference.seq.load_32(fwd_r as usize)) & k_mask) == 0 {
                        matched_any = true;
                        state.raw_hits += 1;
                        let _ = extend_verified_seed::<COLLECT_COV, C>(
                            false, diag, sx64, fwd_r, query, reference, fwd_opts, isa, state,
                        );
                    } else {
                        state.exact_k_rejected += 1;
                    }
                }
            }
            if try_rev && let Some(rr) = rev_ref {
                let target_y = ny - fwd_r - k as u64;
                let diag = sx as i64 - target_y as i64;
                if !is_diag_covered(state, true, diag, sx64) {
                    if ((q_fwd ^ rr.seq.load_32(target_y as usize)) & k_mask) == 0 {
                        matched_any = true;
                        state.raw_hits += 1;
                        let _ = extend_verified_seed::<COLLECT_COV, C>(
                            true, diag, sx64, target_y, query, rr, rev_opts, isa, state,
                        );
                    } else {
                        state.exact_k_rejected += 1;
                    }
                }
            }
        }
    }

    matched_any
}

/// Verify that a minimizer match `(p_q, p_r)` belongs to a valid, non-freq-dropped exact `k`-mer
/// on diagonal `diag`, and extend uncovered surviving exact `k`-mers with lossless valley-trimming.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn verify_and_extend_minimizer_hit<const COLLECT_COV: bool, C: DiagPairOps>(
    is_rev: bool,
    p_q: u64,
    p_r: u64,
    min_score: u32,
    q_run_start: usize,
    q_run_end: usize,
    query: &Genome,
    reference: &Genome,
    target: &Genome,
    min_index: &MinimizerIndex,
    k: usize,
    kp: usize,
    w: usize,
    kp_mask: u64,
    is_multi_slab: bool,
    opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState<C>,
) {
    let ny = reference.seq.len as u64;
    let target_y = if !is_rev { p_r } else { ny - p_r - kp as u64 };
    let diag = p_q as i64 - target_y as i64;

    // Sound early skip: every containing `k`-mer has `seed_x <= p_q`, so if `p_q` is already
    // covered on `diag`, every containing `k`-mer is also covered.
    if is_diag_covered(state, is_rev, diag, p_q) {
        return;
    }

    let seq_x = &query.seq;
    let seq_y = &target.seq;
    let pq = p_q as usize;
    let ty = target_y as usize;

    // 1. Exact `kp`-mer 2-bit verification (filters 32-bit hash-tag collisions).
    if ((seq_x.load_32(pq) ^ seq_y.load_32(ty)) & kp_mask) != 0 {
        state.exact_k_rejected += 1;
        return;
    }

    // 2. Fast raw 2-bit expansion up to `w - 1` bp left/right BEFORE binary-searching reference
    //    `special_runs`/`soft_runs`. If raw 2-bit match length < `k`, reject immediately!
    let max_left_raw = (w - 1).min(pq - q_run_start).min(ty);
    let ext_left_raw = if max_left_raw == 0 {
        0
    } else {
        let diff = seq_x.load_32(pq - max_left_raw) ^ seq_y.load_32(ty - max_left_raw);
        let shifted = diff << (64 - 2 * max_left_raw);
        ((shifted.leading_zeros() >> 1) as usize).min(max_left_raw)
    };

    let max_right_raw = (w - 1)
        .min(q_run_end - (pq + kp))
        .min(seq_y.len - (ty + kp));
    let ext_right_raw = if max_right_raw == 0 {
        0
    } else {
        let diff = seq_x.load_32(pq + kp) ^ seq_y.load_32(ty + kp);
        ((diff.trailing_zeros() >> 1) as usize).min(max_right_raw)
    };

    if ext_left_raw + kp + ext_right_raw < k {
        state.exact_k_rejected += 1;
        return;
    }

    // 3. Clamp against valid seed-run bounds around `p_r` in forward `reference` (and `rev_ref`).
    let (ref_rs, ref_re) = valid_run_around(&reference.seq, p_r as usize, min_index.soft_mask);
    let (ty_rs, ty_re) = if !is_rev {
        (ref_rs, ref_re)
    } else {
        let ny_u = reference.seq.len;
        (ny_u - ref_re, ny_u - ref_rs)
    };

    let ext_left = ext_left_raw.min(ty - ty_rs);
    let ext_right = ext_right_raw.min(ty_re - (ty + kp));
    let match_len = ext_left + kp + ext_right;
    if match_len < k {
        state.exact_k_rejected += 1;
        return;
    }

    // 4. In multi-slab mode, compute the O(w) valid k-mer range around `p_r` ONCE so each
    //    candidate `fwd_r` is checked in O(1) instead of rescanning `w` kp-mers per offset.
    let (valid_r_min, valid_r_max) = if is_multi_slab {
        if ref_re - ref_rs < k {
            state.exact_k_rejected += 1;
            return;
        }
        minimizer_valid_kmer_range(
            &reference.seq,
            p_r as usize,
            ref_rs,
            ref_re,
            k,
            kp,
            w,
            kp_mask,
            min_index.forward_only,
            min_score,
        )
    } else {
        (0usize, usize::MAX)
    };

    // 5. Lossless valley-trimming candidate loop (`W3-F9`):
    //    - Common case (accepted HSP, fast-reject, or positive short reject): `advance_to >= sx + k - 1 > last_sx`
    //      exits immediately after 1 iteration (0 overhead).
    //    - Valley (`safe_b == 0`): `advance_to == sx` -> `off += 1`, trying every remaining
    //      uncovered seed inside `[first_sx..=last_sx]` identically to stride-1 `Exact`!
    let n_cands = match_len - k + 1;
    let first_sx = pq - ext_left;
    let first_sy = ty - ext_left;
    let check_dropped = !min_index.dropped_kmer_bm.is_empty();

    let mut off = 0usize;
    let mut matched_or_covered = false;
    while off < n_cands {
        let sx = (first_sx + off) as u64;
        if is_diag_covered(state, is_rev, diag, sx) {
            matched_or_covered = true;
            off += 1;
            continue;
        }
        let sy = first_sy + off;
        let fwd_r = if !is_rev {
            sy
        } else {
            reference.seq.len - sy - k
        };
        if check_dropped && min_index.is_kmer_dropped(fwd_r) {
            off += 1;
            continue;
        }
        if is_multi_slab && (fwd_r < valid_r_min || fwd_r > valid_r_max) {
            off += 1;
            continue;
        }

        matched_or_covered = true;
        let advance_to = extend_verified_seed::<COLLECT_COV, C>(
            is_rev, diag, sx, sy as u64, query, target, opts, isa, state,
        );
        let next_sx = advance_to.max(sx).saturating_add(1);
        off = (next_sx - first_sx as u64) as usize;
    }

    if !matched_or_covered {
        state.exact_k_rejected += 1;
    }
}

/// Full pairwise comparison using shift-invariant canonical minimizer super-seeding (single-slab
/// or memory-bounded hash-partitioned multi-slab).
pub fn compare_genomes_minimizer_with_coverage(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    collect_coverage: bool,
) -> Result<CompareResult> {
    config.validate()?;
    let forward_only = !config.compare_reverse
        && config.freq_limit.saturating_mul(2) <= config.max_rep.max(1) as u64;

    if MinimizerIndex::fits_in_budget_with_mask(
        reference,
        config.word_length,
        config.soft_mask,
        config.compare_reverse,
        config.memory_budget,
    ) {
        let t_idx = Instant::now();
        let min_index = MinimizerIndex::build(
            reference,
            config.word_length,
            config.freq_limit,
            config.soft_mask,
            forward_only,
        )?;
        let idx_secs = t_idx.elapsed().as_secs_f64();

        let t_rev = Instant::now();
        let rev_ref = if config.compare_reverse {
            Some(reference.reverse_complement())
        } else {
            None
        };
        let rev_secs = t_rev.elapsed().as_secs_f64();

        if profile_enabled() {
            eprintln!(
                "[GECKO_PROFILE] minimizer_build slabs=1 raw_entries={} kept_entries={} ov_entries={} heap_mib={:.2} index_ms={:.2} rev_ms={:.2}",
                min_index.raw_entries,
                min_index.kept_entries,
                min_index.overflow.as_ref().map_or(0, |o| o.kept_entries),
                min_index.heap_bytes() as f64 / (1024.0 * 1024.0),
                idx_secs * 1e3,
                rev_secs * 1e3,
            );
        }

        return compare_minimizer_streaming_with_coverage(
            query,
            reference,
            &min_index,
            rev_ref.as_ref(),
            config,
            collect_coverage,
        );
    }

    // Memory-bounded balanced hash-partitioned minimizer slab streaming.
    let max_entries = MinimizerIndex::max_entries_per_slab(
        reference,
        config.compare_reverse,
        config.memory_budget,
    );
    let slabs = compute_minimizer_partition_slabs(
        reference,
        config.word_length,
        config.soft_mask,
        max_entries,
        forward_only,
    );
    if slabs.is_empty() {
        return Ok(CompareResult {
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            stats: CompareStats::default(),
        });
    }

    let t_rev = Instant::now();
    let rev_ref = if config.compare_reverse {
        Some(reference.reverse_complement())
    } else {
        None
    };
    let rev_secs = t_rev.elapsed().as_secs_f64();
    let n_slabs = slabs.len();

    let mut merged = CompareResult {
        fragments: Vec::new(),
        coverage: CoverageMatrix::new(),
        stats: CompareStats::default(),
    };
    let mut total_idx_secs = 0.0f64;

    for (ps, pe) in slabs {
        let t_idx = Instant::now();
        let slab_index = MinimizerIndex::build_partitions(
            reference,
            ps,
            pe,
            config.word_length,
            config.freq_limit,
            config.soft_mask,
            forward_only,
        )?;
        total_idx_secs += t_idx.elapsed().as_secs_f64();

        let mut slab_res = compare_minimizer_streaming_with_coverage(
            query,
            reference,
            &slab_index,
            rev_ref.as_ref(),
            config,
            collect_coverage,
        )?;
        drop(slab_index);

        merged.stats.word_matches += slab_res.stats.word_matches;
        merged.stats.raw_hits += slab_res.stats.raw_hits;
        merged.stats.filtered_hits += slab_res.stats.filtered_hits;
        merged.stats.hits_used += slab_res.stats.hits_used;
        if collect_coverage {
            merged.coverage.merge_from(&slab_res.coverage);
        }
        merged.fragments.append(&mut slab_res.fragments);
    }

    dedup_fragments(&mut merged.fragments);
    for f in &merged.fragments {
        if f.strand == 'f' {
            merged.stats.forward_fragments += 1;
        } else {
            merged.stats.reverse_fragments += 1;
        }
    }

    if profile_enabled() {
        eprintln!(
            "[GECKO_PROFILE] minimizer_multislab slabs={} index_ms={:.2} rev_ms={:.2} frags={}",
            n_slabs,
            total_idx_secs * 1e3,
            rev_secs * 1e3,
            merged.fragments.len(),
        );
    }

    Ok(merged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::{CompareConfig, SeedingMode, compare_genomes};

    fn pseudo_dna(seed: u64, len: usize) -> String {
        const BASES: [u8; 4] = *b"ACGT";
        let mut x = seed;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            out.push(BASES[((x >> 32) & 3) as usize]);
        }
        String::from_utf8(out).unwrap()
    }

    fn revcomp_str(s: &str) -> String {
        s.bytes()
            .rev()
            .map(crate::seq::complement_base)
            .map(|b| b as char)
            .collect()
    }

    /// Naive O(w)-per-window reference oracle for verifying `extract_block_minimizers_{avx2,portable}`
    /// and `for_each_minimizer` across all edge cases.
    #[allow(clippy::too_many_arguments)]
    fn naive_minimizers(
        seq: &PackedSeq,
        start: usize,
        end: usize,
        k: usize,
        kp: usize,
        w: usize,
        soft_mask: bool,
        emit_ties: bool,
        fwd_only: bool,
    ) -> Vec<(u64, u64, u32, u8, usize, usize)> {
        let kp_mask = (1u64 << (2 * kp)) - 1;
        let mut out = Vec::new();
        seq.for_each_valid_run(start, end, soft_mask, |rs, re| {
            if re - rs < k {
                return;
            }
            let need_right = emit_ties && !fwd_only;
            let n_wins = re - rs - k + 1;
            let n_kps = re - rs - kp + 1;
            let mut kps = Vec::with_capacity(n_kps);
            for i in 0..n_kps {
                let p = rs + i;
                let kf = seq.load_32(p) & kp_mask;
                let (canon, kind) = if fwd_only {
                    (kf, 0u8)
                } else {
                    let kr = reverse_2bit_comp(kf, kp);
                    let knd = if kf < kr {
                        0u8
                    } else if kr < kf {
                        1u8
                    } else {
                        2u8
                    };
                    (kf.min(kr), knd)
                };
                let score = minimizer_order_hash(canon);
                kps.push((p as u64, canon, score, kind));
            }
            let mut marked = vec![false; n_kps];
            for win in 0..n_wins {
                let mut best_s = u32::MAX;
                let mut best_l = win;
                let mut best_r = win;
                for (j, item) in kps.iter().enumerate().skip(win).take(w) {
                    let s = item.2;
                    if s < best_s {
                        best_s = s;
                        best_l = j;
                        best_r = j;
                    } else if s == best_s {
                        best_r = j;
                    }
                }
                marked[best_l] = true;
                if need_right {
                    marked[best_r] = true;
                }
            }
            for (i, &m) in marked.iter().enumerate() {
                if m {
                    let (p, c, s, knd) = kps[i];
                    out.push((p, c, s, knd, rs, re));
                }
            }
        });
        out
    }

    fn run_portable_minimizers<const EMIT_TIES: bool, const FWD_ONLY: bool>(
        seq: &PackedSeq,
        start: usize,
        end: usize,
        k: usize,
        kp: usize,
        w: usize,
        soft_mask: bool,
    ) -> Vec<(u64, u64, u32, u8, usize, usize)> {
        let kp_mask = (1u64 << (2 * kp)) - 1;
        let overlap = w - 1;
        let mut sc = MinScratch::new();
        let mut out = Vec::new();
        seq.for_each_valid_run(start, end, soft_mask, |rs, re| {
            if re - rs < k {
                return;
            }
            let total_wins = re - rs - k + 1;
            let mut emit_bm = [0u64; 5];
            let mut win_done = 0usize;
            let mut carried = 0usize;
            while win_done < total_wins {
                let n_wins = (total_wins - win_done).min(BLOCK_WINS);
                let base_pos = rs + win_done;
                if EMIT_TIES && !FWD_ONLY {
                    extract_block_minimizers_portable::<true, FWD_ONLY>(
                        seq,
                        base_pos,
                        carried,
                        n_wins,
                        kp,
                        w,
                        kp_mask,
                        &mut sc,
                        &mut emit_bm,
                    );
                } else {
                    extract_block_minimizers_portable::<false, FWD_ONLY>(
                        seq,
                        base_pos,
                        carried,
                        n_wins,
                        kp,
                        w,
                        kp_mask,
                        &mut sc,
                        &mut emit_bm,
                    );
                }
                let is_last = win_done + n_wins == total_wins;
                let flush_words = if is_last {
                    (n_wins + overlap).div_ceil(64)
                } else {
                    4
                };
                for (wi, &word) in emit_bm.iter().enumerate().take(flush_words) {
                    let mut m = word;
                    while m != 0 {
                        let bit = m.trailing_zeros() as usize;
                        m &= m - 1;
                        let idx = (wi << 6) + bit;
                        let ck = sc.canon_kind[idx];
                        let canon = ck & KP_CANON_MASK;
                        let kind = (ck >> 62) as u8;
                        let score = sc.packed_l[idx] >> POS_BITS;
                        out.push(((base_pos + idx) as u64, canon, score, kind, rs, re));
                    }
                }
                if !is_last {
                    emit_bm[0] = emit_bm[4];
                    emit_bm[1] = 0;
                    emit_bm[2] = 0;
                    emit_bm[3] = 0;
                    emit_bm[4] = 0;
                    for j in 0..overlap {
                        sc.canon_kind[j] = sc.canon_kind[BLOCK_WINS + j];
                        sc.packed_l[j] = sc.packed_l[BLOCK_WINS + j] - BLOCK_WINS as u32;
                    }
                    carried = overlap;
                }
                win_done += n_wins;
            }
        });
        out
    }

    #[test]
    fn sort_key_bijection_roundtrips_all_bucket_bits() {
        for bucket_bits in 8u32..=26 {
            let hash_shift = 64 - bucket_bits;
            for i in 0..1024u64 {
                let ih = index_hash(i.wrapping_mul(0x517C_C1B7_2722_0A95));
                let sk = idx_hash_to_sort_key(ih, hash_shift);
                assert_eq!(sort_key_to_idx_hash(sk, hash_shift), ih);
            }
        }
    }

    #[test]
    fn reverse_2bit_comp_matches_string_revcomp() {
        let s = pseudo_dna(42, 64);
        let g = Genome::from_records("test", [("test", s.as_bytes())]);
        for k in [18usize, 20, 24, 31, 32] {
            let mask = if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            };
            let fwd = g.seq.load_32(7) & mask;
            let rc_fast = reverse_2bit_comp(fwd, k);
            let rc_s = revcomp_str(&s[7..7 + k]);
            let g_rc = Genome::from_records("rc", [("rc", rc_s.as_bytes())]);
            let rc_expected = g_rc.seq.load_32(0) & mask;
            assert_eq!(rc_fast, rc_expected, "mismatch at k={k}");
        }
    }

    #[test]
    fn avx2_matches_portable_and_naive_all_edge_cases() {
        let pal14 = format!("ACGTCCG{}", revcomp_str("ACGTCCG"));
        let pal16 = format!("ACGTCCGG{}", revcomp_str("ACGTCCGG"));
        let pal18 = format!("ACGTCCGGT{}", revcomp_str("ACGTCCGGT"));
        let pal20 = format!("ACGTCCGGTT{}", revcomp_str("ACGTCCGGTT"));

        let cases: Vec<(&str, Vec<(&str, String)>)> = vec![
            ("short_17", vec![("r", pseudo_dna(1, 17))]),
            ("exact_18", vec![("r", pseudo_dna(2, 18))]),
            ("exact_32", vec![("r", pseudo_dna(3, 32))]),
            ("exact_33", vec![("r", pseudo_dna(4, 33))]),
            ("boundary_255p32", vec![("r", pseudo_dna(5, 255 + 32))]),
            ("boundary_256p32", vec![("r", pseudo_dna(6, 256 + 32))]),
            ("boundary_257p32", vec![("r", pseudo_dna(7, 257 + 32))]),
            ("len_600", vec![("r", pseudo_dna(8, 600))]),
            ("homopoly_400", vec![("r", "A".repeat(400))]),
            ("dinuc_400", vec![("r", "AT".repeat(200))]),
            ("tetranuc_400", vec![("r", "ACGT".repeat(100))]),
            (
                "palindromes",
                vec![(
                    "r",
                    format!(
                        "{}{}{}{}{}{}{}{}{}",
                        pseudo_dna(9, 40),
                        pal14,
                        pseudo_dna(10, 40),
                        pal16,
                        pseudo_dna(11, 40),
                        pal18,
                        pseudo_dna(12, 40),
                        pal20,
                        pseudo_dna(13, 40),
                    ),
                )],
            ),
            (
                "gaps_and_multi_rec",
                vec![
                    (
                        "r1",
                        format!(
                            "{}NNN{}acgtacgtacgtacgtacgtacgtacgtacgtacgt{}",
                            pseudo_dna(14, 150),
                            pseudo_dna(15, 290),
                            pseudo_dna(16, 40),
                        ),
                    ),
                    ("r2", format!("{}N{}", "T".repeat(280), pseudo_dna(17, 300))),
                ],
            ),
        ];

        for (label, recs) in &cases {
            let rec_refs: Vec<(&str, &[u8])> =
                recs.iter().map(|(n, s)| (*n, s.as_bytes())).collect();
            let g = Genome::from_records(*label, rec_refs);
            for k in 18usize..=32 {
                let (kp, w) = minimizer_params(k).unwrap();
                for &soft_mask in &[true, false] {
                    for &(emit_ties, fwd_only) in
                        &[(false, false), (true, false), (false, true), (true, true)]
                    {
                        let oracle = naive_minimizers(
                            &g.seq, 0, g.seq.len, k, kp, w, soft_mask, emit_ties, fwd_only,
                        );
                        let mut actual = Vec::new();
                        let portable = match (emit_ties, fwd_only) {
                            (false, false) => {
                                for_each_minimizer::<false, false, _>(
                                    &g.seq,
                                    0,
                                    g.seq.len,
                                    k,
                                    kp,
                                    w,
                                    soft_mask,
                                    |p, c, s, knd, rs, re| actual.push((p, c, s, knd, rs, re)),
                                );
                                run_portable_minimizers::<false, false>(
                                    &g.seq, 0, g.seq.len, k, kp, w, soft_mask,
                                )
                            }
                            (true, false) => {
                                for_each_minimizer::<true, false, _>(
                                    &g.seq,
                                    0,
                                    g.seq.len,
                                    k,
                                    kp,
                                    w,
                                    soft_mask,
                                    |p, c, s, knd, rs, re| actual.push((p, c, s, knd, rs, re)),
                                );
                                run_portable_minimizers::<true, false>(
                                    &g.seq, 0, g.seq.len, k, kp, w, soft_mask,
                                )
                            }
                            (false, true) => {
                                for_each_minimizer::<false, true, _>(
                                    &g.seq,
                                    0,
                                    g.seq.len,
                                    k,
                                    kp,
                                    w,
                                    soft_mask,
                                    |p, c, s, knd, rs, re| actual.push((p, c, s, knd, rs, re)),
                                );
                                run_portable_minimizers::<false, true>(
                                    &g.seq, 0, g.seq.len, k, kp, w, soft_mask,
                                )
                            }
                            (true, true) => {
                                for_each_minimizer::<true, true, _>(
                                    &g.seq,
                                    0,
                                    g.seq.len,
                                    k,
                                    kp,
                                    w,
                                    soft_mask,
                                    |p, c, s, knd, rs, re| actual.push((p, c, s, knd, rs, re)),
                                );
                                run_portable_minimizers::<true, true>(
                                    &g.seq, 0, g.seq.len, k, kp, w, soft_mask,
                                )
                            }
                        };
                        assert_eq!(
                            actual, oracle,
                            "for_each_minimizer != naive on {label} k={k} soft={soft_mask} ties={emit_ties} fwd={fwd_only}"
                        );
                        assert_eq!(
                            portable, oracle,
                            "portable != naive on {label} k={k} soft={soft_mask} ties={emit_ties} fwd={fwd_only}"
                        );
                        // Verify strictly increasing pos (zero duplicates).
                        for pair in actual.windows(2) {
                            assert!(
                                pair[0].0 < pair[1].0,
                                "non-monotonic pos on {label} k={k}: {} >= {}",
                                pair[0].0,
                                pair[1].0
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn order_hash_consistency_with_load_kp_and_rolling_range() {
        assert_ne!(
            minimizer_order_hash(0),
            0,
            "poly-A (canon=0) must not map to global minimum 0"
        );
        let dna = format!(
            "{}{}{}{}",
            pseudo_dna(77, 180),
            "A".repeat(45),
            "AT".repeat(30),
            pseudo_dna(88, 180),
        );
        let g = Genome::from_records("t", [("t", dna.as_bytes())]);

        for k in 18usize..=32 {
            let (kp, w) = minimizer_params(k).unwrap();
            let kp_mask = (1u64 << (2 * kp)) - 1;
            for &fwd_only in &[false, true] {
                let mut check = |pos: u64,
                                 _canon: u64,
                                 score: u32,
                                 _kind: u8,
                                 rs: usize,
                                 re: usize| {
                    let p_r = pos as usize;
                    assert_eq!(
                        load_kp_score(&g.seq, p_r, kp, kp_mask, fwd_only),
                        score,
                        "load_kp_score mismatch at k={k} fwd={fwd_only} pos={p_r}"
                    );
                    let (r_min, r_max) = minimizer_valid_kmer_range(
                        &g.seq, p_r, rs, re, k, kp, w, kp_mask, fwd_only, score,
                    );
                    let r_lo = p_r.saturating_sub(w - 1).max(rs);
                    let r_hi = p_r.min(re - k);
                    for r in r_lo..=r_hi {
                        let expected =
                            window_has_min_score(&g.seq, r, kp, w, kp_mask, fwd_only, score);
                        let got = r_min <= r_max && r >= r_min && r <= r_max;
                        assert_eq!(
                            got, expected,
                            "rolling range mismatch at k={k} fwd={fwd_only} p_r={p_r} r={r} range=({r_min},{r_max})"
                        );
                    }
                };
                if fwd_only {
                    for_each_minimizer::<true, true, _>(
                        &g.seq, 0, g.seq.len, k, kp, w, true, &mut check,
                    );
                } else {
                    for_each_minimizer::<true, false, _>(
                        &g.seq, 0, g.seq.len, k, kp, w, true, &mut check,
                    );
                }
            }
        }
    }

    #[test]
    fn overflow_routing_palindromes_freq_limit_and_valley_trimming_parity() {
        // 1. High-frequency repeat copies (> MINIMIZER_OVERFLOW_CAP) sharing a 24-bp core
        //    with unique flanks, plus palindromic blocks and finite freq_limit = {0, 1, 20}.
        let core = pseudo_dna(901, 24);
        let kept_hsp = format!("{}{}{}", pseudo_dna(902, 80), core, pseudo_dna(903, 80));
        let rev_hsp = format!("{}{}{}", pseudo_dna(904, 80), core, pseudo_dna(905, 80));
        let half = "ACGTCCGGTTAAACGT";
        let pal32 = format!("{half}{}", revcomp_str(half));
        let pal_hsp = format!("{}{}{}", pseudo_dna(906, 60), pal32, pseudo_dna(907, 60));

        let mut r_seq = String::new();
        for i in 0..40u64 {
            r_seq.push_str(&pseudo_dna(1000 + i, 50));
            r_seq.push_str(&core);
            r_seq.push_str(&pseudo_dna(2000 + i, 50));
        }
        r_seq.push_str(&kept_hsp);
        r_seq.push_str(&pseudo_dna(3001, 200));
        r_seq.push_str(&revcomp_str(&rev_hsp));
        r_seq.push_str(&pseudo_dna(3002, 200));
        for _ in 0..10 {
            r_seq.push_str(&pal_hsp);
            r_seq.push_str(&pseudo_dna(3003, 40));
        }

        let q_seq = format!(
            "{}{}{}{}{}{}{}",
            pseudo_dna(4001, 200),
            kept_hsp,
            pseudo_dna(4002, 200),
            rev_hsp,
            pseudo_dna(4003, 200),
            pal_hsp,
            pseudo_dna(4004, 200),
        );

        let q = Genome::from_records("q", [("q", q_seq.as_bytes())]);
        let r = Genome::from_records("r", [("r", r_seq.as_bytes())]);

        // Verify overflow CSR actually populated when `core` repeats > 40 times.
        let idx = MinimizerIndex::build(&r, 32, 20, false, false).unwrap();
        assert!(
            idx.overflow.as_ref().is_some_and(|o| o.kept_entries > 0),
            "expected OverflowKmerCsr to populate for 40x repeated core"
        );

        for &fl in &[0u64, 1, 20] {
            for &fwd_only in &[false, true] {
                for &budget in &[64 * 1024usize, 64 * 1024 * 1024] {
                    let cfg_exact = CompareConfig {
                        word_length: 32,
                        min_length: 100,
                        min_similarity: 75.0,
                        freq_limit: fl,
                        soft_mask: false,
                        compare_reverse: !fwd_only,
                        seeding_mode: SeedingMode::Exact,
                        ..Default::default()
                    };
                    let cfg_min = CompareConfig {
                        seeding_mode: SeedingMode::Minimizer,
                        memory_budget: budget,
                        ..cfg_exact.clone()
                    };
                    let res_ex = compare_genomes(&q, &r, &cfg_exact, None).unwrap();
                    let res_mn = compare_genomes(&q, &r, &cfg_min, None).unwrap();
                    assert_eq!(
                        res_ex.fragments, res_mn.fragments,
                        "overflow parity mismatch at fl={fl} fwd_only={fwd_only} budget={budget}"
                    );
                    if fl == 0 {
                        assert!(res_mn.fragments.is_empty());
                    }
                }
            }
        }

        // 2. Valley-trimming 74-bp MEM HSP recovery test (`safe_b == 0`):
        //    500-bp high-identity prefix (400 matches / 100 mismatches), followed by a deep valley
        //    of 42 mismatches, followed by a 74-bp exact MEM + 1,000-bp high-identity HSP.
        let mut qv = pseudo_dna(5001, 200).into_bytes();
        let mut rv = pseudo_dna(5002, 200).into_bytes();
        let mem74 = pseudo_dna(5003, 74).into_bytes();
        let hsp_tail = pseudo_dna(5004, 1000).into_bytes();
        // Inject 74-bp MEM + 1000-bp HSP with 6.3% mismatches, preceded by a valley right before it.
        let mut rv_tail = hsp_tail.clone();
        for i in (15..1000).step_by(16) {
            rv_tail[i] = match rv_tail[i] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        qv.extend_from_slice(&mem74);
        qv.extend_from_slice(&hsp_tail);
        qv.extend_from_slice(pseudo_dna(5005, 200).as_bytes());
        rv.extend_from_slice(&mem74);
        rv.extend_from_slice(&rv_tail);
        rv.extend_from_slice(pseudo_dna(5006, 200).as_bytes());

        let qg = Genome::from_records("qv", [("qv", qv.as_slice())]);
        let rg = Genome::from_records("rv", [("rv", rv.as_slice())]);
        let cfg_ex = CompareConfig {
            word_length: 32,
            min_length: 1000,
            min_similarity: 75.0,
            freq_limit: 20,
            seeding_mode: SeedingMode::Exact,
            ..Default::default()
        };
        let cfg_mn = CompareConfig {
            seeding_mode: SeedingMode::Minimizer,
            ..cfg_ex.clone()
        };
        let r_ex = compare_genomes(&qg, &rg, &cfg_ex, None).unwrap();
        let r_mn = compare_genomes(&qg, &rg, &cfg_mn, None).unwrap();
        assert_eq!(r_ex.fragments, r_mn.fragments);
        assert_eq!(r_mn.fragments.len(), 1);
    }

    #[test]
    fn minimizer_matches_exact_across_k_palindromes_and_freq_limit() {
        let h1 = pseudo_dna(101, 200);
        let h2 = pseudo_dna(102, 180);
        let rep = pseudo_dna(103, 90);
        let half = "ACGTCCGGTTAAACGT";
        let pal32 = format!("{half}{}", revcomp_str(half));
        let pal_block = format!("{}{}{}", pseudo_dna(104, 40), pal32, pseudo_dna(105, 40));
        let tandem = "ACGTACGT".repeat(16);

        let q_seq = format!(
            "{}{}{}{}{}{}{}{}{}{}{}",
            pseudo_dna(201, 200),
            h1,
            pseudo_dna(202, 200),
            h2,
            pseudo_dna(203, 200),
            rep,
            pseudo_dna(204, 200),
            pal_block,
            pseudo_dna(205, 200),
            tandem,
            pseudo_dna(206, 200),
        );
        let r_seq = format!(
            "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}",
            pseudo_dna(301, 200),
            h1,
            pseudo_dna(302, 200),
            revcomp_str(&h2),
            pseudo_dna(303, 200),
            rep,
            pseudo_dna(304, 100),
            rep,
            pseudo_dna(305, 100),
            rep,
            pseudo_dna(306, 200),
            pal_block,
            pseudo_dna(307, 200),
            tandem,
            pseudo_dna(308, 200),
        );

        let q = Genome::from_records("q", [("q", q_seq.as_bytes())]);
        let r = Genome::from_records("r", [("r", r_seq.as_bytes())]);

        for &k in &[18usize, 20, 24, 32] {
            for &fl in &[2u64, 100] {
                for &fwd_only in &[false, true] {
                    let cfg_exact = CompareConfig {
                        word_length: k,
                        min_length: 60,
                        min_similarity: 60.0,
                        freq_limit: fl,
                        compare_reverse: !fwd_only,
                        seeding_mode: SeedingMode::Exact,
                        ..Default::default()
                    };
                    let cfg_min = CompareConfig {
                        seeding_mode: SeedingMode::Minimizer,
                        ..cfg_exact.clone()
                    };
                    let res_ex = compare_genomes(&q, &r, &cfg_exact, None).unwrap();
                    let res_mn = compare_genomes(&q, &r, &cfg_min, None).unwrap();
                    assert_eq!(
                        res_ex.fragments, res_mn.fragments,
                        "minimizer vs exact mismatch at k={k} fl={fl} fwd_only={fwd_only}"
                    );
                }
            }
        }
    }

    #[test]
    fn minimizer_multi_slab_matches_exact() {
        let h1 = pseudo_dna(501, 220);
        let h2 = pseudo_dna(502, 200);
        let rep = pseudo_dna(503, 80);
        let q_seq = format!(
            "{}{}{}{}{}{}{}",
            pseudo_dna(510, 4000),
            h1,
            pseudo_dna(511, 4000),
            h2,
            pseudo_dna(512, 4000),
            rep,
            pseudo_dna(513, 4000),
        );
        let r_seq = format!(
            "{}{}{}{}{}{}{}{}{}{}{}",
            pseudo_dna(520, 4000),
            h1,
            pseudo_dna(521, 4000),
            revcomp_str(&h2),
            pseudo_dna(522, 2000),
            rep,
            pseudo_dna(523, 2000),
            rep,
            pseudo_dna(524, 2000),
            rep,
            pseudo_dna(525, 4000),
        );
        let q = Genome::from_records("q", [("q", q_seq.as_bytes())]);
        let r = Genome::from_records("r", [("r", r_seq.as_bytes())]);

        let cfg_exact = CompareConfig {
            word_length: 32,
            min_length: 100,
            min_similarity: 70.0,
            freq_limit: 2,
            seeding_mode: SeedingMode::Exact,
            ..Default::default()
        };
        let cfg_min_tiny = CompareConfig {
            seeding_mode: SeedingMode::Minimizer,
            memory_budget: 64 * 1024,
            ..cfg_exact.clone()
        };
        let res_ex = compare_genomes(&q, &r, &cfg_exact, None).unwrap();
        let res_mn = compare_genomes(&q, &r, &cfg_min_tiny, None).unwrap();
        assert_eq!(res_ex.fragments, res_mn.fragments);
        assert_eq!(res_mn.fragments.len(), 2);
    }
}
