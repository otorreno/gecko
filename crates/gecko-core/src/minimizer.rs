//! Shift-invariant canonical minimizer super-seeding engine (`F9`) with 100% lossless exact-$k$
//! verification and exact per-strand `freq_limit` resolution.
//!
//! For contiguous seeds $k \in [18, 32]$, every exact $k$-mer window of length $k$ contains
//! $w = k - k' + 1$ sub-$k'$-mers (e.g. $k' = 21, w = 12$ for $k = 32$). Indexing all canonical
//! minimizer positions of the reference (including hash ties for strict strand symmetry) and
//! probing the leftmost canonical minimizer per query window shrinks index size, build work, and
//! query probes by ~3x–6.5x while guaranteeing 100% bit-exact `.frags` parity with stride-1
//! `RefIndex` — including finite per-strand `freq_limit`, palindromes, soft-masking, multi-FASTA,
//! `--forward-only`, and hash-partitioned multi-slab memory budgets.
//!
//! Crucial hashing detail: the minimizer ordering score `score = fmix64(canon)` is exponentially
//! skewed toward `0` (minimum of `w` uniform hashes). Therefore, all partition/bucket/Bloom/tag
//! indexing uses an independent bijective secondary hash `idx_hash = index_hash(canon)` so all 256
//! partitions, CSR buckets, and L3 Blocked Bloom blocks remain 100% uniformly loaded.

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
    BlockedBloom, BloomBlock, PosStrandVec, choose_bloom_params, fmix64, insert_in_bloom_block,
    prefetch_t1, prefetch_t2, probe_bucket_binary_u32, probe_sorted_tags_scalar,
};
use crate::seq::{Genome, PackedSeq};
use crate::stream::dedup_fragments;

/// Minimum contiguous seed length eligible for canonical minimizer super-seeding.
pub const MINIMIZER_MIN_K: usize = 18;

const NUM_PARTITIONS: usize = 256;
const MIN_QUERY_CHUNK_BYTES: usize = 1 << 20; // ~1 Mb dynamic Rayon chunks
const MIN_BATCH_SIZE: usize = 64;
const DEQUE_CAP: usize = 16; // power of 2 strictly > max window (12)
const DEQUE_MASK: usize = DEQUE_CAP - 1;

/// Independent salt mixed into `canon` before `fmix64` to decorrelate partition/bucket/Bloom
/// hashing from the minimum-biased minimizer selection score `fmix64(canon)`.
const IDX_HASH_SALT: u64 = 0x9E37_79B9_7F4A_7C15;

const HSP_CACHE_BITS: u32 = 17;
const HSP_CACHE_LEN: usize = 1 << HSP_CACHE_BITS;
const NOISE_CACHE_BITS: u32 = 16;
const NOISE_CACHE_LEN: usize = 1 << NOISE_CACHE_BITS;
const FIB_MULT: u64 = 0x9E37_79B9_7F4A_7C15;

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
    debug_assert!((2..DEQUE_CAP).contains(&window));
    Some((kp, window))
}

/// Uniform bijective index hash decorrelated from the minimum-biased minimizer score `fmix64(canon)`.
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

#[derive(Clone, Copy, Default)]
struct DequeEntry {
    score: u64,
    canon: u64,
    pos: u64,
    kind: u8,
}

/// Extract canonical minimizers from `[start .. end)` of `seq` using a zero-heap stack circular
/// monotone deque.
///
/// - When `EMIT_TIES == true` (reference indexing): emits every position achieving the window
///   minimum score in at least one valid `k`-mer window (strictly increasing `pos`, zero duplicates),
///   guaranteeing 100% strand-invariant recall even under hash ties and palindromic `k`-mers.
/// - When `EMIT_TIES == false` (query scanning): emits the leftmost argmin per valid `k`-mer
///   window (deduplicated across consecutive windows sharing the same minimizer position).
///
/// Calls `f(pos, canon, score, kind, run_start, run_end)` where `kind` is `0` (fwd), `1` (rev), or
/// `2` (palindrome).
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
    F: FnMut(u64, u64, u64, u8, usize, usize),
{
    let kp_mask = (1u64 << (2 * kp)) - 1;
    let kp_rc_shift = (2 * (kp - 1)) as u32;

    seq.for_each_valid_run(start, end, soft_mask, |rs, re| {
        if re - rs < k {
            return;
        }
        let mut deque = [DequeEntry::default(); DEQUE_CAP];
        let mut head = 0usize;
        let mut tail = 0usize;
        let mut last_emitted_pos = u64::MAX;

        let mut kf = 0u64;
        let mut kr = 0u64;
        let mut idx = rs;
        let prime_end = rs + kp - 1;
        while idx < prime_end {
            let w_idx = idx >> 5;
            let lane = idx & 31;
            let avail = (32 - lane).min(prime_end - idx);
            let mut word = unsafe { *seq.bits.get_unchecked(w_idx) } >> (lane << 1);
            for _ in 0..avail {
                let b = word & 3;
                word >>= 2;
                kf = ((kf << 2) | b) & kp_mask;
                if !FWD_ONLY {
                    kr = (kr >> 2) | ((b ^ 3) << kp_rc_shift);
                }
            }
            idx += avail;
        }

        let first_win_end_pos = (rs + w - 1) as u64;
        while idx < re {
            let w_idx = idx >> 5;
            let lane = idx & 31;
            let avail = (32 - lane).min(re - idx);
            let mut word = unsafe { *seq.bits.get_unchecked(w_idx) } >> (lane << 1);
            let start_pos = (idx + 1 - kp) as u64;
            for pos in start_pos..start_pos + avail as u64 {
                let b = word & 3;
                word >>= 2;
                kf = ((kf << 2) | b) & kp_mask;
                let (canon, kind) = if FWD_ONLY {
                    (kf, 0u8)
                } else {
                    kr = (kr >> 2) | ((b ^ 3) << kp_rc_shift);
                    if kf < kr {
                        (kf, 0u8)
                    } else if kr < kf {
                        (kr, 1u8)
                    } else {
                        (kf, 2u8)
                    }
                };
                let score = fmix64(canon);

                // Strictly `>` pop so score ties are preserved in the deque.
                while head != tail {
                    let back_idx = tail.wrapping_sub(1) & DEQUE_MASK;
                    if unsafe { deque.get_unchecked(back_idx).score } > score {
                        tail = tail.wrapping_sub(1);
                    } else {
                        break;
                    }
                }
                unsafe {
                    *deque.get_unchecked_mut(tail & DEQUE_MASK) = DequeEntry {
                        score,
                        canon,
                        pos,
                        kind,
                    };
                }
                tail = tail.wrapping_add(1);

                if pos >= first_win_end_pos {
                    let win_start = pos + 1 - w as u64;
                    while head != tail
                        && unsafe { deque.get_unchecked(head & DEQUE_MASK).pos } < win_start
                    {
                        head = head.wrapping_add(1);
                    }
                    let front = unsafe { *deque.get_unchecked(head & DEQUE_MASK) };
                    if EMIT_TIES {
                        let min_s = front.score;
                        let mut cur = head;
                        while cur != tail {
                            let e = unsafe { *deque.get_unchecked(cur & DEQUE_MASK) };
                            if e.score != min_s {
                                break;
                            }
                            if last_emitted_pos == u64::MAX || e.pos > last_emitted_pos {
                                f(e.pos, e.canon, e.score, e.kind, rs, re);
                                last_emitted_pos = e.pos;
                            }
                            cur = cur.wrapping_add(1);
                        }
                    } else if front.pos != last_emitted_pos {
                        f(front.pos, front.canon, front.score, front.kind, rs, re);
                        last_emitted_pos = front.pos;
                    }
                }
            }
            idx += avail;
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

/// Compute the minimizer ordering score `fmix64(canon_kp)` for the `kp`-mer at `pos` in `seq`.
#[inline(always)]
fn load_kp_score(
    seq: &PackedSeq,
    pos: usize,
    kp: usize,
    kp_mask: u64,
    kp_rc_shift: u32,
    forward_only: bool,
) -> u64 {
    let mut word = seq.load_32(pos);
    let mut kf = 0u64;
    let mut kr = 0u64;
    for _ in 0..kp {
        let b = word & 3;
        word >>= 2;
        kf = ((kf << 2) | b) & kp_mask;
        if !forward_only {
            kr = (kr >> 2) | ((b ^ 3) << kp_rc_shift);
        }
    }
    let canon = if forward_only { kf } else { kf.min(kr) };
    fmix64(canon)
}

/// Check whether the valid `k`-mer window `[r .. r + k)` has minimum canonical `kp`-mer score
/// strictly equal to `target_score`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn window_has_min_score(
    seq: &PackedSeq,
    r: usize,
    kp: usize,
    w: usize,
    kp_mask: u64,
    kp_rc_shift: u32,
    forward_only: bool,
    target_score: u64,
) -> bool {
    let mut word = seq.load_32(r);
    let mut kf = 0u64;
    let mut kr = 0u64;
    for _ in 0..kp - 1 {
        let b = word & 3;
        word >>= 2;
        kf = ((kf << 2) | b) & kp_mask;
        if !forward_only {
            kr = (kr >> 2) | ((b ^ 3) << kp_rc_shift);
        }
    }
    let mut found = false;
    for _ in 0..w {
        let b = word & 3;
        word >>= 2;
        kf = ((kf << 2) | b) & kp_mask;
        let canon = if forward_only {
            kf
        } else {
            kr = (kr >> 2) | ((b ^ 3) << kp_rc_shift);
            kf.min(kr)
        };
        let s = fmix64(canon);
        if s < target_score {
            return false;
        }
        if s == target_score {
            found = true;
        }
    }
    found
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

/// Compact canonical Hash-CSR + L3 Blocked Bloom index over shift-invariant reference minimizers.
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
    /// 1-bit/bp bitmap over `reference`: bit `r == 1` iff the `k`-mer starting at forward
    /// reference coordinate `r` exceeded `freq_limit` and was dropped. Empty when no `k`-mer in
    /// the indexed slab(s) exceeded `freq_limit`.
    dropped_kmer_bm: Vec<u64>,
    pub kept_entries: usize,
    pub raw_entries: usize,
}

impl MinimizerIndex {
    /// Approximate heap bytes held by this minimizer index.
    pub fn heap_bytes(&self) -> usize {
        self.offsets.capacity() * 4
            + self.tags.capacity() * 4
            + self.pos_strand.heap_bytes()
            + self.bloom.as_ref().map_or(0, BlockedBloom::heap_bytes)
            + self.dropped_kmer_bm.capacity() * 8
    }

    /// Estimate whether a single-slab `MinimizerIndex` (plus `rev_ref`) fits inside `memory_budget`.
    pub fn fits_in_budget(
        reference: &Genome,
        word_length: usize,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> bool {
        let Some((_kp, w)) = minimizer_params(word_length) else {
            return false;
        };
        let rev_bytes = if compare_reverse {
            reference.seq.len.div_ceil(4) + 128 * 1024
        } else {
            0
        };
        let est_entries = (reference.seq.len * 3) / (w + 1);
        let peak_bytes = est_entries.saturating_mul(28) + rev_bytes + 64 * 1024;
        peak_bytes <= memory_budget
    }

    /// Compute the maximum minimizer entries per slab for `memory_budget`.
    pub fn max_entries_per_slab(
        reference: &Genome,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> usize {
        let rev_bytes = if compare_reverse {
            reference.seq.len.div_ceil(4)
        } else {
            0
        };
        let avail = memory_budget.saturating_sub(rev_bytes);
        (avail / 28).max(512)
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
        let use_u32_ps = seq.len < (1usize << 31);
        if seq.len < word_length || freq_limit == 0 {
            return Ok(Self {
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
                dropped_kmer_bm: Vec::new(),
                kept_entries: 0,
                raw_entries: 0,
            });
        }

        let n_parts = (part_end - part_start) as usize;
        let part_filter = part_start > 0 || part_end < NUM_PARTITIONS as u32;
        let n_threads = rayon::current_num_threads().max(1);
        let chunks = make_chunk_ranges(seq.len, word_length, n_threads * 2);

        // Pass 1: Parallel chunk extraction of reference minimizers into uniform `index_hash(canon)`
        // coarse partitions (decorrelated from minimum-biased minimizer selection scores).
        let est_per_chunk_part =
            ((seq.len / (window + 1)) * 2 / (chunks.len().max(1) * NUM_PARTITIONS)).max(16);
        let chunk_vecs: Vec<Vec<Vec<(u64, u64)>>> = chunks
            .par_iter()
            .map(|&(cs, ce)| {
                let mut local: Vec<Vec<(u64, u64)>> = (0..n_parts)
                    .map(|_| Vec::with_capacity(est_per_chunk_part))
                    .collect();
                let mut emit =
                    |pos: u64, canon: u64, _score: u64, kind: u8, _rs: usize, _re: usize| {
                        let ih = index_hash(canon);
                        let p = (ih >> 56) as u32;
                        if part_filter && (p < part_start || p >= part_end) {
                            return;
                        }
                        let strand_bit = u64::from(kind == 1);
                        let ps = (pos << 1) | strand_bit;
                        unsafe {
                            local
                                .get_unchecked_mut((p - part_start) as usize)
                                .push((ih, ps));
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
                local
            })
            .collect();

        // Transpose into 1 contiguous Vec per partition.
        let mut part_entries: Vec<Vec<(u64, u64)>> = (0..n_parts)
            .map(|p_idx| {
                let cap: usize = chunk_vecs.iter().map(|cv| cv[p_idx].len()).sum();
                Vec::with_capacity(cap)
            })
            .collect();
        for mut cv in chunk_vecs {
            for (p_idx, mut v) in cv.drain(..).enumerate() {
                part_entries[p_idx].append(&mut v);
            }
        }

        let raw_entries: usize = part_entries.iter().map(Vec::len).sum();
        if raw_entries == 0 {
            return Ok(Self {
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
                dropped_kmer_bm: Vec::new(),
                kept_entries: 0,
                raw_entries: 0,
            });
        }

        let bucket_bits = choose_bucket_bits(raw_entries);
        let sub_bits = bucket_bits - 8;
        let hash_shift = 64 - bucket_bits;
        let buckets_per_part = 1usize << sub_bits;
        let bloom_params = choose_bloom_params(raw_entries, n_parts, freq_limit);
        let bloom_sub_bits = bloom_params.map_or(0u32, |(sb, _)| sb);
        let bloom_k_h = bloom_params.map_or(0usize, |(_, kh)| kh as usize);
        let bloom_blocks_per_part = if bloom_params.is_some() {
            1usize << bloom_sub_bits
        } else {
            0
        };
        let bloom_shift = 56 - bloom_sub_bits;

        let kp_mask = (1u64 << (2 * kp)) - 1;
        let kp_rc_shift = (2 * (kp - 1)) as u32;

        struct PartBuildOut {
            offsets: Vec<u32>,
            tags: Vec<u32>,
            ps: Vec<u64>,
            bloom_blocks: Vec<BloomBlock>,
            dropped_r: Vec<usize>,
        }

        // Pass 2: Parallel single-sort per partition using the bijective sort key
        // `[bucket_bits | tag32 | mid_bits]`, performing exact `freq_limit` resolution and
        // emitting sorted partition CSR + L3 Blocked Bloom in a single linear pass.
        let part_outs: Vec<PartBuildOut> = part_entries
            .into_par_iter()
            .enumerate()
            .map(|(p_idx, mut entries)| {
                for e in &mut entries {
                    e.0 = idx_hash_to_sort_key(e.0, hash_shift);
                }
                entries.sort_unstable();
                entries.dedup();

                let p_global = part_start + p_idx as u32;
                let part_bucket_base = (p_global << sub_bits) as usize;
                let part_bloom_base = (p_global << bloom_sub_bits) as usize;

                let mut bloom_blocks = if bloom_blocks_per_part > 0 {
                    vec![BloomBlock([0u64; 8]); bloom_blocks_per_part]
                } else {
                    Vec::new()
                };

                let mut part_offsets = vec![0u32; buckets_per_part + 1];
                let mut part_tags = Vec::with_capacity(entries.len());
                let mut part_ps = Vec::with_capacity(entries.len());
                let mut dropped_r: Vec<usize> = Vec::new();
                let mut cand_r: Vec<usize> = Vec::new();
                let mut kmer_items: Vec<(u64, u8, usize)> = Vec::new();
                let mut dropped_local: Vec<usize> = Vec::new();

                let mut i = 0usize;
                while i < entries.len() {
                    let sk = entries[i].0;
                    let mut j = i + 1;
                    while j < entries.len() && entries[j].0 == sk {
                        j += 1;
                    }
                    let run_len = j - i;
                    let ih = sort_key_to_idx_hash(sk, hash_shift);
                    let b_local = ((ih >> hash_shift) as usize) - part_bucket_base;
                    let tag = ih as u32;
                    let before_len = part_tags.len();

                    if (run_len as u64) <= freq_limit {
                        part_offsets[b_local + 1] += run_len as u32;
                        for &(_, ps) in &entries[i..j] {
                            part_tags.push(tag);
                            part_ps.push(ps);
                        }
                    } else {
                        // ~1-2% high-frequency minimizer bucket: exact `k`-mer `freq_limit` check.
                        let p_first = (entries[i].1 >> 1) as usize;
                        let min_score =
                            load_kp_score(seq, p_first, kp, kp_mask, kp_rc_shift, forward_only);

                        cand_r.clear();
                        for &(_, ps) in &entries[i..j] {
                            let p_r = (ps >> 1) as usize;
                            let (rs, re) = valid_run_around(seq, p_r, soft_mask);
                            if re - rs < word_length {
                                continue;
                            }
                            let r_lo = p_r.saturating_sub(window - 1).max(rs);
                            let r_hi = p_r.min(re - word_length);
                            for r in r_lo..=r_hi {
                                if window_has_min_score(
                                    seq,
                                    r,
                                    kp,
                                    window,
                                    kp_mask,
                                    kp_rc_shift,
                                    forward_only,
                                    min_score,
                                ) {
                                    cand_r.push(r);
                                }
                            }
                        }
                        cand_r.sort_unstable();
                        cand_r.dedup();

                        kmer_items.clear();
                        kmer_items.reserve(cand_r.len());
                        for &r in &cand_r {
                            let (kc, ks) =
                                load_kmer_canon_strand(seq, r, word_length, forward_only);
                            kmer_items.push((kc, ks, r));
                        }
                        kmer_items.sort_unstable();

                        dropped_local.clear();
                        let mut a = 0usize;
                        while a < kmer_items.len() {
                            let (kc, ks, _) = kmer_items[a];
                            let mut b = a + 1;
                            while b < kmer_items.len()
                                && kmer_items[b].0 == kc
                                && kmer_items[b].1 == ks
                            {
                                b += 1;
                            }
                            if ((b - a) as u64) > freq_limit {
                                for &(_, _, r) in &kmer_items[a..b] {
                                    dropped_local.push(r);
                                    dropped_r.push(r);
                                }
                            }
                            a = b;
                        }

                        if dropped_local.is_empty() {
                            part_offsets[b_local + 1] += run_len as u32;
                            for &(_, ps) in &entries[i..j] {
                                part_tags.push(tag);
                                part_ps.push(ps);
                            }
                        } else {
                            dropped_local.sort_unstable();
                            let mut kept_cnt = 0u32;
                            for &(_, ps) in &entries[i..j] {
                                let p_r = (ps >> 1) as usize;
                                let r_lo = p_r.saturating_sub(window - 1);
                                let start_idx = cand_r.partition_point(|&r| r < r_lo);
                                let mut survives = false;
                                for &r in &cand_r[start_idx..] {
                                    if r > p_r {
                                        break;
                                    }
                                    if dropped_local.binary_search(&r).is_err() {
                                        survives = true;
                                        break;
                                    }
                                }
                                if survives {
                                    part_tags.push(tag);
                                    part_ps.push(ps);
                                    kept_cnt += 1;
                                }
                            }
                            part_offsets[b_local + 1] += kept_cnt;
                        }
                    }

                    if bloom_blocks_per_part > 0 && part_tags.len() > before_len {
                        let blk_local = ((ih >> bloom_shift) as usize) - part_bloom_base;
                        insert_in_bloom_block(&mut bloom_blocks[blk_local], tag, bloom_k_h);
                    }
                    i = j;
                }

                for b in 0..buckets_per_part {
                    part_offsets[b + 1] += part_offsets[b];
                }

                PartBuildOut {
                    offsets: part_offsets,
                    tags: part_tags,
                    ps: part_ps,
                    bloom_blocks,
                    dropped_r,
                }
            })
            .collect();

        // Assemble flat partition-segmented CSR arrays matching `RefIndex`'s `b_idx = b + (b >> sub_bits)` layout.
        let kept_entries: usize = part_outs.iter().map(|p| p.tags.len()).sum();
        let any_dropped = part_outs.iter().any(|p| !p.dropped_r.is_empty());
        let mut dropped_kmer_bm = if any_dropped {
            vec![0u64; seq.len.div_ceil(64)]
        } else {
            Vec::new()
        };

        let mut offsets = Vec::with_capacity(n_parts * (buckets_per_part + 1));
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
        let mut bloom_blocks = if bloom_blocks_per_part > 0 {
            Vec::with_capacity(n_parts * bloom_blocks_per_part)
        } else {
            Vec::new()
        };

        let mut base = 0u32;
        for po in part_outs {
            for off in po.offsets {
                offsets.push(base + off);
            }
            base += po.tags.len() as u32;
            tags.extend_from_slice(&po.tags);
            if use_u32_ps {
                ps_u32.extend(po.ps.iter().map(|&v| v as u32));
            } else {
                ps_u64.extend_from_slice(&po.ps);
            }
            if bloom_blocks_per_part > 0 {
                bloom_blocks.extend_from_slice(&po.bloom_blocks);
            }
            for r in po.dropped_r {
                dropped_kmer_bm[r >> 6] |= 1u64 << (r & 63);
            }
        }

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

        let pos_strand = if use_u32_ps {
            PosStrandVec::U32(ps_u32)
        } else {
            PosStrandVec::U64(ps_u64)
        };

        Ok(Self {
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
            pos_strand,
            bloom,
            dropped_kmer_bm,
            kept_entries,
            raw_entries,
        })
    }

    #[inline(always)]
    fn is_kmer_dropped(&self, fwd_r: usize) -> bool {
        if self.dropped_kmer_bm.is_empty() {
            return false;
        }
        ((unsafe { *self.dropped_kmer_bm.get_unchecked(fwd_r >> 6) }) >> (fwd_r & 63)) & 1 != 0
    }
}

/// Partition the 256 coarse minimizer hash partitions into contiguous slabs each holding at most
/// `max_entries_per_slab` raw minimizer entries.
pub fn compute_minimizer_partition_slabs(
    reference: &Genome,
    word_length: usize,
    soft_mask: bool,
    max_entries_per_slab: usize,
    forward_only: bool,
) -> Vec<(u32, u32)> {
    let Some((kp, window)) = minimizer_params(word_length) else {
        return Vec::new();
    };
    let seq = &reference.seq;
    if seq.len < word_length {
        return Vec::new();
    }
    let n_threads = rayon::current_num_threads().max(1);
    let chunks = make_chunk_ranges(seq.len, word_length, n_threads * 2);
    let counts = chunks
        .par_iter()
        .map(|&(cs, ce)| {
            let mut c = [0usize; NUM_PARTITIONS];
            let mut emit =
                |_pos: u64, canon: u64, _score: u64, _kind: u8, _rs: usize, _re: usize| {
                    c[(index_hash(canon) >> 56) as usize] += 1;
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
            c
        })
        .reduce(
            || [0usize; NUM_PARTITIONS],
            |mut a, b| {
                for i in 0..NUM_PARTITIONS {
                    a[i] += b[i];
                }
                a
            },
        );

    let total: usize = counts.iter().sum();
    if total == 0 {
        return Vec::new();
    }
    let limit = max_entries_per_slab.max(1);
    let mut slabs = Vec::new();
    let mut s = 0usize;
    while s < NUM_PARTITIONS {
        let mut acc = counts[s];
        let mut e = s + 1;
        while e < NUM_PARTITIONS && acc + counts[e] <= limit {
            acc += counts[e];
            e += 1;
        }
        slabs.push((s as u32, e as u32));
        s = e;
    }
    slabs
}

// ---------------------------------------------------------------------------
// Minimizer streaming query engine (batched Bloom + CSR prefetch + exact-k verify)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct DiagSlot {
    diag: i64,
    max_x: u64,
}

impl DiagSlot {
    const EMPTY: Self = Self {
        diag: i64::MIN,
        max_x: 0,
    };
}

struct MinDiagCache {
    hsp: Vec<DiagSlot>,
    noise: Vec<DiagSlot>,
}

impl MinDiagCache {
    fn new() -> Self {
        Self {
            hsp: vec![DiagSlot::EMPTY; HSP_CACHE_LEN],
            noise: vec![DiagSlot::EMPTY; NOISE_CACHE_LEN],
        }
    }

    fn clear(&mut self) {
        self.hsp.fill(DiagSlot::EMPTY);
        self.noise.fill(DiagSlot::EMPTY);
    }

    #[inline(always)]
    fn is_covered(&self, diag: i64, pos_x: u64) -> bool {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let hs = unsafe {
            *self
                .hsp
                .get_unchecked((h >> (64 - HSP_CACHE_BITS)) as usize)
        };
        if hs.diag == diag && pos_x <= hs.max_x {
            return true;
        }
        let ns = unsafe {
            *self
                .noise
                .get_unchecked((h >> (64 - NOISE_CACHE_BITS)) as usize)
        };
        ns.diag == diag && pos_x <= ns.max_x
    }

    #[inline(always)]
    fn update_hsp(&mut self, diag: i64, max_x: u64) {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let slot = unsafe {
            self.hsp
                .get_unchecked_mut((h >> (64 - HSP_CACHE_BITS)) as usize)
        };
        if slot.diag == diag {
            if max_x > slot.max_x {
                slot.max_x = max_x;
            }
        } else {
            *slot = DiagSlot { diag, max_x };
        }
    }

    #[inline(always)]
    fn update_noise(&mut self, diag: i64, max_x: u64) {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let slot = unsafe {
            self.noise
                .get_unchecked_mut((h >> (64 - NOISE_CACHE_BITS)) as usize)
        };
        if slot.diag == diag {
            if max_x > slot.max_x {
                slot.max_x = max_x;
            }
        } else {
            *slot = DiagSlot { diag, max_x };
        }
    }
}

#[derive(Clone, Copy, Default)]
struct MinQueryItem {
    pos_x: u64,
    score: u64,
    run_start: usize,
    run_end: usize,
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

struct MinWorkerState {
    fwd_cache: MinDiagCache,
    rev_cache: Option<MinDiagCache>,
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
}

impl MinWorkerState {
    fn new(compare_reverse: bool, collect_coverage: bool) -> Self {
        Self {
            fwd_cache: MinDiagCache::new(),
            rev_cache: if compare_reverse {
                Some(MinDiagCache::new())
            } else {
                None
            },
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
        }
    }

    #[inline(always)]
    fn sp_x_cursor_for(&mut self, sp_x: &[(u64, u64)], pos_x: u64) -> usize {
        let mut c = self.sp_x_cursor;
        while c < sp_x.len() && unsafe { sp_x.get_unchecked(c).1 } <= pos_x {
            c += 1;
        }
        self.sp_x_cursor = c;
        c
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
fn prefetch_ps(ps: &PosStrandVec, idx: usize) {
    match ps {
        PosStrandVec::U32(v) => prefetch_t2(unsafe { v.as_ptr().add(idx) }),
        PosStrandVec::U64(v) => prefetch_t2(unsafe { v.as_ptr().add(idx) }),
    }
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

    let worker_results: Vec<MinWorkerState> = chunks
        .into_par_iter()
        .fold(
            || MinWorkerState::new(compare_reverse, collect_coverage),
            |mut state, (chunk_start, chunk_end)| {
                if chunk_start < state.last_chunk_end {
                    state.fwd_cache.clear();
                    if let Some(rc) = state.rev_cache.as_mut() {
                        rc.clear();
                    }
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
                scan_minimizer_chunk(
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

    for mut w in worker_results {
        stats.word_matches += w.word_matches;
        stats.raw_hits += w.raw_hits;
        stats.hits_used += w.hits_used;
        total_q_mins += w.query_minimizers;
        total_bloom_rej += w.bloom_rejected;
        total_exact_rej += w.exact_k_rejected;
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
            "[GECKO_PROFILE] minimizer_query k={} kp={} w={} parts=[{}..{}) ref_entries={} q_mins={} bloom_reject={} exact_k_reject={} hits_used={} frags={} query_ms={:.2} dedup_ms={:.2}",
            min_index.word_length,
            min_index.kp,
            min_index.window,
            min_index.part_start,
            min_index.part_end,
            min_index.kept_entries,
            total_q_mins,
            total_bloom_rej,
            total_exact_rej,
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
fn scan_minimizer_chunk(
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
    state: &mut MinWorkerState,
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

    let mut handle = |pos_x: u64, canon: u64, score: u64, q_kind: u8, rs: usize, re: usize| {
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
            score,
            run_start: rs,
            run_end: re,
            b_idx,
            tag: ih as u32,
            bloom_blk,
            start: 0,
            end: 0,
            q_kind,
        };
        batch_len += 1;
        if batch_len == MIN_BATCH_SIZE {
            flush_minimizer_batch(
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
        flush_minimizer_batch(
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
fn flush_minimizer_batch(
    batch: &mut [MinQueryItem],
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    min_index: &MinimizerIndex,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState,
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
        let mut ns = 0usize;
        for (i, item) in batch.iter().enumerate() {
            if unsafe { bloom.contains_unchecked(item.bloom_blk as usize, item.tag) } {
                prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
                survivors[ns] = i as u8;
                ns += 1;
            }
        }
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
    let kp_rc_shift = (2 * (kp - 1)) as u32;
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

        let mut idx = lo;
        let mut matched_any = false;
        while idx < hi {
            let (p_r, r_strand) = unpack_ps(&min_index.pos_strand, idx);
            idx += step;

            if fwd_only {
                matched_any = true;
                state.raw_hits += 1;
                verify_and_extend_minimizer_hit(
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
                    kp_rc_shift,
                    is_multi_slab,
                    fwd_opts,
                    isa,
                    state,
                );
            } else if item.q_kind == 2 {
                matched_any = true;
                state.raw_hits += 1;
                verify_and_extend_minimizer_hit(
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
                    kp_rc_shift,
                    is_multi_slab,
                    fwd_opts,
                    isa,
                    state,
                );
                if let Some(rr) = rev_ref {
                    state.raw_hits += 1;
                    verify_and_extend_minimizer_hit(
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
                        kp_rc_shift,
                        is_multi_slab,
                        rev_opts,
                        isa,
                        state,
                    );
                }
            } else if r_strand == item.q_kind {
                matched_any = true;
                state.raw_hits += 1;
                verify_and_extend_minimizer_hit(
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
                    kp_rc_shift,
                    is_multi_slab,
                    fwd_opts,
                    isa,
                    state,
                );
            } else if let Some(rr) = rev_ref {
                matched_any = true;
                state.raw_hits += 1;
                verify_and_extend_minimizer_hit(
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
                    kp_rc_shift,
                    is_multi_slab,
                    rev_opts,
                    isa,
                    state,
                );
            }
        }
        if matched_any {
            state.word_matches += 1;
        }
    }
}

/// Verify that a minimizer match `(p_q, p_r)` belongs to a valid, non-freq-dropped exact `k`-mer
/// on diagonal `diag`, and extend from the leftmost **uncovered** surviving exact `k`-mer.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn verify_and_extend_minimizer_hit(
    is_rev: bool,
    p_q: u64,
    p_r: u64,
    min_score: u64,
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
    kp_rc_shift: u32,
    is_multi_slab: bool,
    opts: ExtendOptions,
    isa: MinIsa,
    state: &mut MinWorkerState,
) {
    let ny = reference.seq.len as u64;
    let target_y = if !is_rev { p_r } else { ny - p_r - kp as u64 };
    let diag = p_q as i64 - target_y as i64;

    // Sound early skip: every containing `k`-mer has `seed_x <= p_q`, so if `p_q` is already
    // covered on `diag`, every containing `k`-mer is also covered.
    if !is_rev {
        if state.fwd_cache.is_covered(diag, p_q) {
            return;
        }
    } else if state.rev_cache.as_ref().unwrap().is_covered(diag, p_q) {
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

    // 2. Valid seed-run bounds around `p_r` in forward `reference` (and mapped to `rev_ref`).
    let (ref_rs, ref_re) = valid_run_around(&reference.seq, p_r as usize, min_index.soft_mask);
    let (ty_rs, ty_re) = if !is_rev {
        (ref_rs, ref_re)
    } else {
        let ny_u = reference.seq.len;
        (ny_u - ref_re, ny_u - ref_rs)
    };

    // 3. Exact 2-bit expansion up to `w - 1` bp left and right within valid unmasked bounds.
    let max_left = (w - 1).min(pq - q_run_start).min(ty - ty_rs);
    let ext_left = if max_left == 0 {
        0
    } else {
        let diff = seq_x.load_32(pq - max_left) ^ seq_y.load_32(ty - max_left);
        let shifted = diff << (64 - 2 * max_left);
        ((shifted.leading_zeros() >> 1) as usize).min(max_left)
    };

    let max_right = (w - 1).min(q_run_end - (pq + kp)).min(ty_re - (ty + kp));
    let ext_right = if max_right == 0 {
        0
    } else {
        let diff = seq_x.load_32(pq + kp) ^ seq_y.load_32(ty + kp);
        ((diff.trailing_zeros() >> 1) as usize).min(max_right)
    };

    let match_len = ext_left + kp + ext_right;
    if match_len < k {
        state.exact_k_rejected += 1;
        return;
    }

    // 4. Select the leftmost UNCOVERED valid surviving `k`-mer in
    //    `[pq - ext_left ..= pq + kp + ext_right - k]`. Checking `is_covered(diag, sx)` INSIDE
    //    this loop guarantees 100% parity with stride-1 GECKO when an earlier `update_noise`
    //    boundary lands partway through an exact match run.
    let n_cands = match_len - k + 1;
    let first_sx = pq - ext_left;
    let first_sy = ty - ext_left;
    let check_dropped = !min_index.dropped_kmer_bm.is_empty();

    let mut chosen_sx: Option<u64> = None;
    for off in 0..n_cands {
        let sx = (first_sx + off) as u64;
        if !is_rev {
            if state.fwd_cache.is_covered(diag, sx) {
                continue;
            }
        } else if state.rev_cache.as_ref().unwrap().is_covered(diag, sx) {
            continue;
        }
        let sy = first_sy + off;
        let fwd_r = if !is_rev {
            sy
        } else {
            reference.seq.len - sy - k
        };
        if check_dropped && min_index.is_kmer_dropped(fwd_r) {
            continue;
        }
        if is_multi_slab
            && !window_has_min_score(
                &reference.seq,
                fwd_r,
                kp,
                w,
                kp_mask,
                kp_rc_shift,
                min_index.forward_only,
                min_score,
            )
        {
            continue;
        }
        chosen_sx = Some(sx);
        break;
    }

    let Some(seed_x) = chosen_sx else {
        state.exact_k_rejected += 1;
        return;
    };
    let seed_y = (seed_x as i64 - diag) as u64;

    // 5. Launch exact GECKO bidirectional drop-to-zero HSP extension from `(seed_x, seed_y)`.
    state.hits_used += 1;
    let seed_ident = k as u64;
    let seed_score = k as i64 * opts.scoring.match_score;
    let cx_x = state.sp_x_cursor_for(&query.seq.special_runs, seed_x);
    let cov_opt = state.coverage.as_mut();
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
                seed_ident,
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
                seed_ident,
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
            seed_ident,
            seed_score,
            cov_opt,
        ),
    };

    if let Some(frag) = outcome.fragment {
        if !is_rev {
            state.fwd_cache.update_hsp(diag, outcome.xfil_max);
        } else {
            state
                .rev_cache
                .as_mut()
                .unwrap()
                .update_hsp(diag, outcome.xfil_max);
        }
        state.fragments.push(frag);
    } else if outcome.xfil_max > seed_x + 256 {
        let skip_end = exact_rejected_skip(
            &query.seq,
            &target.seq,
            seed_x,
            seed_y,
            outcome.xfil_max,
            k,
            opts.scoring,
            &mut state.skip_scratch,
        );
        if !is_rev {
            state.fwd_cache.update_hsp(diag, skip_end);
        } else {
            state.rev_cache.as_mut().unwrap().update_hsp(diag, skip_end);
        }
    } else {
        let skip = seed_x + 2 * k as u64 - 1;
        if !is_rev {
            state.fwd_cache.update_noise(diag, skip);
        } else {
            state.rev_cache.as_mut().unwrap().update_noise(diag, skip);
        }
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

    if MinimizerIndex::fits_in_budget(
        reference,
        config.word_length,
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
                "[GECKO_PROFILE] minimizer_build slabs=1 raw_entries={} kept_entries={} heap_mib={:.2} index_ms={:.2} rev_ms={:.2}",
                min_index.raw_entries,
                min_index.kept_entries,
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

    // Memory-bounded hash-partitioned minimizer slab streaming.
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
