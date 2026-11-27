use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::error::{GeckoError, Result};
use crate::seq::{Genome, PackedSeq, reverse_32_2bit};

const NUM_PARTITIONS: usize = 256;
const MAX_HASH_COARSE_PARTITIONS: usize = 1024;
const STAGE_CAP: usize = 64; // 192 KiB HashCsr<u32> / 128 KiB DirectCsr<u32> fits cleanly in 256 KiB L2; 4x fewer flushes/RFOs/dTLB misses than 16
const UNIQ_BUF_CAP: usize = 64;

const BLOOM_MIN_KMERS: usize = 4_096;
const BLOOM_MAX_BYTES: usize = 32 * 1024 * 1024; // 32 MiB full-index cap (128 KiB/part, fits in 45 MiB L3)
const BLOOM_TARGET_BITS_PER_KMER: usize = 8; // 16 MiB on 16M-kmer chr21 (k_h=6), 32 MiB on 50M-kmer chr16 (k_h=4)
const BLOOM_MAX_K_H: usize = 8; // Retained for BLOOM_SALTS, contains_unchecked, & unit tests
const BLOOM_CHOOSE_MAX_K_H: usize = 6; // Putze blocked-Bloom sweet spot (<=0.04% FPR delta, -14%..-25% CPU)

/// Odd 32-bit multiply-shift salts (Putze / Impala family) for 9-bit lane hashes
/// inside a 512-bit cache-line block. Eliminates double-hash even-stride periodicity.
const BLOOM_SALTS: [u32; BLOOM_MAX_K_H] = [
    0x47b6_137b,
    0x4497_4d91,
    0x8824_ad5b,
    0xa2b7_289d,
    0x7054_95c7,
    0x2df1_424b,
    0x9efc_4947,
    0x5c6b_fb31,
];

/// 2-bit DNA encoding lookup table with soft-masking (uppercase A/C/G/T/U only).
pub(crate) const BASE_LUT_SOFT_MASK: [u8; 256] = make_base_lut(true);

/// 2-bit DNA encoding lookup table without soft-masking (case-insensitive A/C/G/T/U).
pub(crate) const BASE_LUT_NO_MASK: [u8; 256] = make_base_lut(false);

const fn make_base_lut(soft_mask: bool) -> [u8; 256] {
    let mut lut = [0xFFu8; 256];
    lut[b'A' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b'U' as usize] = 3;
    if !soft_mask {
        lut[b'a' as usize] = 0;
        lut[b'c' as usize] = 1;
        lut[b'g' as usize] = 2;
        lut[b't' as usize] = 3;
        lut[b'u' as usize] = 3;
    }
    lut
}

/// Return the 2-bit base encoding lookup table for the given `soft_mask` setting.
#[inline(always)]
pub fn base_lut(soft_mask: bool) -> &'static [u8; 256] {
    if soft_mask {
        &BASE_LUT_SOFT_MASK
    } else {
        &BASE_LUT_NO_MASK
    }
}

/// Invertible MurmurHash3 64-bit finalizer (`fmix64`) used for uniform bucket assignment.
#[inline(always)]
pub fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^= k >> 33;
    k
}

/// Exact algebraic inverse of MurmurHash3 [`fmix64`].
#[cfg(test)]
#[inline(always)]
pub(crate) fn inv_fmix64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0x9cb4_b2f8_1293_37db);
    h ^= h >> 33;
    h = h.wrapping_mul(0x4f74_430c_22a5_4005);
    h ^= h >> 33;
    h
}

/// Returns `true` iff a canonical (possibly spaced-masked) `wl`-bp MSB-first key is a palindrome
/// (`sf == sr`). Safe for all `wl in 4..=32`.
#[inline(always)]
fn is_palindrome_key(canon: u64, wl: usize, seed_mask: u64) -> bool {
    let rc = (reverse_32_2bit(!canon) >> (64 - 2 * wl)) & seed_mask;
    rc == canon
}

/// Fast setup check: an odd-length contiguous k-mer (or odd spaced seed whose center position
/// `(wl - 1) / 2` is a care `'1'` base) can never equal its reverse complement.
#[inline(always)]
fn can_key_be_palindrome(wl: usize, spaced_mask: Option<u64>) -> bool {
    (wl & 1) == 0 || spaced_mask.is_some_and(|sm| ((sm >> (wl - 1)) & 3) == 0)
}

/// Software prefetch into L2 cache (`_MM_HINT_T1` on x86_64) for L3-resident Bloom blocks
/// without evicting L1d active-diagonal ring or batch stack buffers.
#[inline(always)]
pub(crate) fn prefetch_t1<T>(ptr: *const T) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        core::arch::x86_64::_mm_prefetch(ptr as *const i8, core::arch::x86_64::_MM_HINT_T1);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = ptr;
    }
}

/// Software prefetch into L3/L2 cache (`_MM_HINT_T2` on x86_64, matching C `__builtin_prefetch(..., 0, 1)`
/// so 64-wide random index prefetches do not evict hot diagonal-cache lines from L1d/L2).
#[inline(always)]
pub(crate) fn prefetch_t2<T>(ptr: *const T) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        core::arch::x86_64::_mm_prefetch(ptr as *const i8, core::arch::x86_64::_MM_HINT_T2);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = ptr;
    }
}

/// Single 64-byte cache-line-aligned 512-bit Blocked Bloom filter block.
#[derive(Debug, Clone, Copy)]
#[repr(C, align(64))]
pub(crate) struct BloomBlock(pub [u64; 8]);

/// Partition-aligned L3-sized 64-byte cache-line-blocked Bloom filter built inline in Pass 2
/// from the existing `fmix64(canon)` hash (zero extra hashing, zero atomics).
#[derive(Debug, Clone)]
pub(crate) struct BlockedBloom {
    pub blocks: Vec<BloomBlock>,
    pub sub_bits: u32,
    pub block_shift: u32, // 56 - sub_bits
    pub k_h: u8,          // 1..=8
}

impl BlockedBloom {
    #[inline]
    pub(crate) fn heap_bytes(&self) -> usize {
        self.blocks.capacity() * std::mem::size_of::<BloomBlock>()
    }

    /// Compile-time-specialized branchless single-cache-line probe (`K_H in 1..=8`).
    ///
    /// Eliminates the dynamic `match self.k_h` indirect jump table (`Bi`) and unrolls
    /// straight-line with immediate `BLOOM_SALTS[i]` constants.
    ///
    /// # Safety
    /// Caller must ensure `blk_idx < self.blocks.len()` and `1 <= K_H && K_H <= BLOOM_MAX_K_H`.
    #[inline(always)]
    pub(crate) unsafe fn contains_unchecked_const<const K_H: usize>(
        &self,
        blk_idx: usize,
        tag: u32,
    ) -> bool {
        debug_assert!((1..=BLOOM_MAX_K_H).contains(&K_H));
        let words = unsafe { &self.blocks.get_unchecked(blk_idx).0 };
        let mut hit = u64::MAX;
        let mut i = 0usize;
        while i < K_H {
            let h = tag.wrapping_mul(BLOOM_SALTS[i]);
            let word = (h >> 29) as usize;
            let bit = (h >> 23) & 63;
            hit &= unsafe { *words.get_unchecked(word) } >> bit;
            i += 1;
        }
        (hit & 1) != 0
    }

    /// Safe compile-time-specialized single-slab Bloom probe directly from a 64-bit hash `h`.
    /// For multi-slab or precomputed `(blk_idx, tag)` batch probes, use
    /// [`Self::contains_unchecked_const`].
    #[inline(always)]
    #[allow(dead_code)]
    pub(crate) fn contains_const<const K_H: usize>(&self, h: u64) -> bool {
        let blk_idx = (h >> self.block_shift) as usize;
        if blk_idx >= self.blocks.len() {
            return true;
        }
        unsafe { self.contains_unchecked_const::<K_H>(blk_idx, h as u32) }
    }

    /// Branchless single-cache-line probe. Zero inner-loop branches eliminates branch
    /// mispredictions on divergent missing k-mers.
    ///
    /// # Safety
    /// Caller must ensure `blk_idx < self.blocks.len()`.
    #[inline(always)]
    pub(crate) unsafe fn contains_unchecked(&self, blk_idx: usize, tag: u32) -> bool {
        match self.k_h {
            1 => unsafe { self.contains_unchecked_const::<1>(blk_idx, tag) },
            2 => unsafe { self.contains_unchecked_const::<2>(blk_idx, tag) },
            3 => unsafe { self.contains_unchecked_const::<3>(blk_idx, tag) },
            4 => unsafe { self.contains_unchecked_const::<4>(blk_idx, tag) },
            5 => unsafe { self.contains_unchecked_const::<5>(blk_idx, tag) },
            6 => unsafe { self.contains_unchecked_const::<6>(blk_idx, tag) },
            7 => unsafe { self.contains_unchecked_const::<7>(blk_idx, tag) },
            _ => unsafe { self.contains_unchecked_const::<8>(blk_idx, tag) },
        }
    }
}

#[inline(always)]
pub(crate) fn insert_in_bloom_block(block: &mut BloomBlock, tag: u32, k_h: usize) {
    macro_rules! set_bit {
        ($i:expr) => {{
            let bit = (tag.wrapping_mul(BLOOM_SALTS[$i]) >> 23) as usize;
            unsafe {
                *block.0.get_unchecked_mut(bit >> 6) |= 1u64 << (bit & 63);
            }
        }};
    }
    match k_h {
        1 => {
            set_bit!(0);
        }
        2 => {
            set_bit!(0);
            set_bit!(1);
        }
        3 => {
            set_bit!(0);
            set_bit!(1);
            set_bit!(2);
        }
        4 => {
            set_bit!(0);
            set_bit!(1);
            set_bit!(2);
            set_bit!(3);
        }
        5 => {
            set_bit!(0);
            set_bit!(1);
            set_bit!(2);
            set_bit!(3);
            set_bit!(4);
        }
        6 => {
            set_bit!(0);
            set_bit!(1);
            set_bit!(2);
            set_bit!(3);
            set_bit!(4);
            set_bit!(5);
        }
        7 => {
            set_bit!(0);
            set_bit!(1);
            set_bit!(2);
            set_bit!(3);
            set_bit!(4);
            set_bit!(5);
            set_bit!(6);
        }
        _ => {
            set_bit!(0);
            set_bit!(1);
            set_bit!(2);
            set_bit!(3);
            set_bit!(4);
            set_bit!(5);
            set_bit!(6);
            set_bit!(7);
        }
    }
}

/// Choose `(sub_bits, k_h)` for the partition-aligned Blocked Bloom filter, or `None` when
/// the reference is tiny, `freq_limit == 0`, or the filter would be over-saturated.
#[inline]
pub(crate) fn choose_bloom_params(
    total_kmers: usize,
    n_parts: usize,
    freq_limit: u64,
) -> Option<(u32, u8)> {
    if total_kmers < BLOOM_MIN_KMERS || freq_limit == 0 || n_parts == 0 {
        return None;
    }
    // Scale-invariant per-partition cap: at most `BLOOM_MAX_BYTES / 64 / NUM_PARTITIONS`
    // (= 2,048 blocks = 128 KiB per coarse partition, `max_sub_bits = 11`).
    // Multi-slab runs (`n_parts < 256`) get the exact same bits/k-mer and `k_h` as single-slab
    // (`10.75 MiB` at 3 slabs, `6.50 MiB` at 5 slabs, `k_h = 4`) instead of 2x-4x oversizing Bloom
    // RAM and inflating Pass-2 insert work (`k_h = 7`).
    let max_per_part = ((BLOOM_MAX_BYTES / 64) / NUM_PARTITIONS).max(1);
    let max_sub_bits = max_per_part.ilog2(); // 11

    let desired_total_blocks = total_kmers
        .saturating_mul(BLOOM_TARGET_BITS_PER_KMER)
        .div_ceil(512);
    let desired_per_part = desired_total_blocks.div_ceil(n_parts).max(1);
    let sub_bits = desired_per_part
        .next_power_of_two()
        .trailing_zeros()
        .min(max_sub_bits);

    let total_blocks = n_parts << sub_bits;
    let total_bits = total_blocks.saturating_mul(512);
    // Disable if over-saturated (< 0.7 bits/k-mer).
    if total_bits.saturating_mul(10) < total_kmers.saturating_mul(7) {
        return None;
    }
    let bits_per_kmer = (total_bits as f64) / (total_kmers as f64);
    let k_h = (bits_per_kmer * std::f64::consts::LN_2)
        .round()
        .clamp(1.0, BLOOM_CHOOSE_MAX_K_H as f64) as u8;
    Some((sub_bits, k_h))
}

/// Packed `(pos << 1) | strand_bit` storage (`u32` for references < 2^31 bp, `u64` fallback).
#[derive(Debug, Clone)]
pub enum PosStrandVec {
    U32(Vec<u32>),
    U64(Vec<u64>),
}

impl PosStrandVec {
    #[inline(always)]
    pub fn len(&self) -> usize {
        match self {
            Self::U32(v) => v.len(),
            Self::U64(v) => v.len(),
        }
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return `(pos_y, strand_bit)` at `index`.
    #[inline(always)]
    pub fn get_unpacked(&self, index: usize) -> (u64, u8) {
        match self {
            Self::U32(v) => {
                let raw = v[index];
                ((raw >> 1) as u64, (raw & 1) as u8)
            }
            Self::U64(v) => {
                let raw = v[index];
                (raw >> 1, (raw & 1) as u8)
            }
        }
    }

    pub fn heap_bytes(&self) -> usize {
        match self {
            Self::U32(v) => v.capacity() * std::mem::size_of::<u32>(),
            Self::U64(v) => v.capacity() * std::mem::size_of::<u64>(),
        }
    }
}

/// Storage backend for [`RefIndex`].
#[derive(Debug, Clone)]
pub(crate) enum RefIndexKind {
    /// Direct CSR indexed by canonical key directly (zero key storage).
    /// `offsets` has `n_parts * (sub_keys + 1)` entries so each coarse partition has its own
    /// sentinel slot (`b_idx = b_local + (b_local >> sub_bits)`).
    DirectCsr {
        sub_bits: u32,
        offsets: Vec<u32>,
        pos_strand: PosStrandVec,
    },
    /// 2-level parallel radix-partitioned Hash-CSR keyed by `fmix64(canonical_key)`.
    /// Stores 32-bit hash tags (`fmix64(canon) as u32`, padded with 8 trailing sentinel zeros for
    /// branchless 8-wide AVX2 probing), halving key memory vs `Vec<u64>` while preserving exact
    /// global `freq_limit` filtering and 100% bit-exact seed verification on uncovered hits.
    HashCsr {
        sub_bits: u32,
        hash_shift: u32,
        offsets: Vec<u32>,
        tags: Vec<u32>,
        pos_strand: PosStrandVec,
        bloom: Option<BlockedBloom>,
    },
}

/// Single-reference canonical (or forward-only) k-mer index supporting lock-free streaming queries.
///
/// Indexes only Reference Y once at `word_length` (`4..=32`) using canonical keys
/// `min(kf, kr)` and a 1-bit strand flag `(kr < kf)` (or direct forward keys `kf` when
/// `forward_only` is enabled), with per-strand `freq_limit` filtering.
/// Supports exact global-`freq_limit` multi-slab streaming by restricting to a contiguous range
/// of coarse hash/prefix partitions `[part_start .. part_end)` within `0..256`.
#[derive(Debug, Clone)]
pub struct RefIndex {
    pub word_length: usize,
    pub freq_limit: u64,
    pub soft_mask: bool,
    /// Optional 2-bit palindromic spaced-seed mask (`None` for contiguous seeds).
    pub spaced_mask: Option<u64>,
    pub ref_data_len: usize,
    pub total_kmers: usize,
    pub kept_entries: usize,
    pub part_start: u16,
    pub part_end: u16,
    pub(crate) kind: RefIndexKind,
    /// 1-bit/bp bitmap encoding whether a reference position's k-mer is valid, kept in this slab,
    /// non-palindromic (when canonical), and globally unique (`cnt == 1` across both strands).
    /// Read as `((unique_bm[pos >> 6] ^ unique_bm_xor) >> (pos & 63)) & 1`.
    pub(crate) unique_bm: Vec<u64>,
    /// Sense mask for `unique_bm`: `0` on `DirectCsr` (bit=1 means unique) and `u64::MAX` on
    /// `HashCsr` (bit=0 means unique, so Pass 2 only writes duplicates — 10x-145x fewer writes).
    pub(crate) unique_bm_xor: u64,
    /// Whether this index keys directly on forward k-mers `kf` (`--forward-only` specialization).
    pub(crate) forward_only: bool,
}

/// Whether Direct CSR mode should be used for `(word_length, ref_total_bases)`.
#[inline]
pub fn use_direct_csr(word_length: usize, ref_total_bases: u64) -> bool {
    word_length <= 10
        || (word_length <= 12
            && (1u64 << (2 * word_length)) <= 2u64.saturating_mul(ref_total_bases))
}

#[inline]
fn use_direct_csr_for(word_length: usize, ref_total_bases: u64, spaced_mask: Option<u64>) -> bool {
    spaced_mask.is_none() && use_direct_csr(word_length, ref_total_bases)
}

/// Choose Hash-CSR bucket bits `P` in `14..=26` targeting ~3–4 entries per bucket.
#[inline]
fn choose_bucket_bits(total_kmers: usize) -> u32 {
    let mut p = 14u32;
    while p < 26 && (total_kmers >> p) >= 4 {
        p += 1;
    }
    p
}

/// Choose Hash-CSR coarse partition bit width `p_bits in {8, 9, 10}` (`256 / 512 / 1024`
/// coarse partitions) so Pass 2's per-thread `sub_off` stays inside 256 KiB L2 cache.
#[inline]
fn choose_hash_part_bits(equiv_total: usize) -> u32 {
    if equiv_total > 67_108_864 {
        10 // 1024 coarse parts: bucket_bits=25 -> 128 KiB, bucket_bits=26 -> 256 KiB
    } else if equiv_total > 20_000_000 {
        9 // 512 coarse parts: bucket_bits=23 -> 64 KiB, bucket_bits=24 -> 128 KiB
    } else {
        8 // 256 coarse parts: bucket_bits<=23 -> <=128 KiB
    }
}

/// Compute cache-line-aligned per-partition flat slot capacity for Single-Pass Atomic Scatter.
#[inline]
fn choose_hash_slot_cap(expected_slab_kmers: usize, n_coarse: usize) -> usize {
    let mean = expected_slab_kmers.div_ceil(n_coarse.max(1));
    let slack = (mean / 64).max(6usize.saturating_mul(mean.isqrt())).max(16);
    mean.saturating_add(slack).next_multiple_of(64)
}

/// Count exact valid k-mers in `reference.seq` in O(#runs) without touching sequence bits.
pub fn count_valid_kmers(reference: &Genome, word_length: usize, soft_mask: bool) -> usize {
    count_valid_kmers_in_range(&reference.seq, 0, reference.seq.len, word_length, soft_mask)
}

/// Count exact valid k-mers whose start positions fall in `[slab_start .. slab_end)`.
pub fn count_valid_kmers_in_range(
    seq: &PackedSeq,
    slab_start: usize,
    slab_end: usize,
    word_length: usize,
    soft_mask: bool,
) -> usize {
    if seq.len < word_length || slab_start >= slab_end {
        return 0;
    }
    let max_starts = seq.len - word_length + 1;
    let end = slab_end.min(max_starts);
    if slab_start >= end {
        return 0;
    }
    let scan_end = (end + word_length - 1).min(seq.len);
    let mut count = 0usize;
    seq.for_each_valid_run(slab_start, scan_end, soft_mask, |rs, re| {
        if re - rs >= word_length {
            count += re - rs - word_length + 1;
        }
    });
    count
}

/// Partition valid k-mer start positions `[0 .. seq.len - wl + 1)` into disjoint reference slabs,
/// each containing at most `max_kmers_per_slab` valid k-mers, for zero-spill slab streaming.
pub fn compute_slab_ranges(
    reference: &Genome,
    word_length: usize,
    soft_mask: bool,
    max_kmers_per_slab: usize,
) -> Vec<(usize, usize)> {
    let seq = &reference.seq;
    if seq.len < word_length {
        return Vec::new();
    }
    let max_starts = seq.len - word_length + 1;
    let cap = max_kmers_per_slab.max(1);
    let mut slabs = Vec::new();
    let mut cur_slab_start = 0usize;
    let mut cur_kmers = 0usize;

    seq.for_each_valid_run(0, seq.len, soft_mask, |rs, re| {
        if re - rs < word_length {
            return;
        }
        let mut run_pos = rs;
        let run_last_start_excl = re - word_length + 1;
        while run_pos < run_last_start_excl {
            let room = cap - cur_kmers;
            let avail = run_last_start_excl - run_pos;
            if avail < room {
                cur_kmers += avail;
                run_pos = run_last_start_excl;
            } else {
                let split_end = run_pos + room;
                slabs.push((cur_slab_start, split_end));
                cur_slab_start = split_end;
                cur_kmers = 0;
                run_pos = split_end;
            }
        }
    });

    if cur_slab_start < max_starts && (cur_kmers > 0 || slabs.is_empty()) {
        slabs.push((cur_slab_start, max_starts));
    }
    slabs
}

/// Partition the 256 coarse k-mer partitions (`0..256`) into disjoint contiguous ranges
/// `[part_start .. part_end)` whose k-mer counts each stay within `max_kmers_per_slab`.
///
/// Because every occurrence of a canonical k-mer across the entire reference maps to the same
/// coarse partition, per-strand `freq_limit` filtering inside each slab is **100% globally exact**,
/// and query streaming skips `(S - 1) / S` of all k-mers in O(1) before touching CSR tables.
pub fn compute_partition_slabs(
    reference: &Genome,
    word_length: usize,
    soft_mask: bool,
    max_kmers_per_slab: usize,
) -> Vec<(u16, u16)> {
    compute_partition_slabs_with_mask(reference, word_length, None, soft_mask, max_kmers_per_slab)
}

/// Partition the 256 coarse k-mer partitions (`0..256`) with optional spaced-seed mask.
pub fn compute_partition_slabs_with_mask(
    reference: &Genome,
    word_length: usize,
    spaced_mask: Option<u64>,
    soft_mask: bool,
    max_kmers_per_slab: usize,
) -> Vec<(u16, u16)> {
    compute_partition_slabs_with_mask_and_strand(
        reference,
        word_length,
        spaced_mask,
        soft_mask,
        max_kmers_per_slab,
        false,
    )
}

/// Partition the 256 coarse k-mer partitions (`0..256`) with optional spaced-seed mask and
/// forward-only strand specialization.
pub fn compute_partition_slabs_with_mask_and_strand(
    reference: &Genome,
    word_length: usize,
    spaced_mask: Option<u64>,
    soft_mask: bool,
    max_kmers_per_slab: usize,
    forward_only: bool,
) -> Vec<(u16, u16)> {
    let total_kmers = count_valid_kmers(reference, word_length, soft_mask);
    if total_kmers == 0 {
        return Vec::new();
    }
    let cap = max_kmers_per_slab.max(1);
    if total_kmers <= cap {
        return vec![(0, NUM_PARTITIONS as u16)];
    }

    if !use_direct_csr_for(word_length, reference.total_bases, spaced_mask) {
        // Hash-CSR partitions by top 8 bits of `fmix64(key)`, which is uniform to <0.2%.
        let num_slabs = total_kmers.div_ceil(cap).clamp(1, NUM_PARTITIONS);
        let mut slabs = Vec::with_capacity(num_slabs);
        for s in 0..num_slabs {
            let ps = (s * NUM_PARTITIONS / num_slabs) as u16;
            let pe = ((s + 1) * NUM_PARTITIONS / num_slabs) as u16;
            if ps < pe {
                slabs.push((ps, pe));
            }
        }
        return slabs;
    }

    // Direct-CSR partitions by top 8 bits of `canon` (or `kf` in forward-only mode),
    // so compute the exact 256-bucket histogram and greedily pack contiguous ranges <= cap.
    let seq = &reference.seq;
    let sub_bits = (2 * word_length - 8) as u32;
    let max_starts = seq.len.saturating_sub(word_length - 1);
    let ranges = make_chunk_ranges(0, max_starts, seq.len, word_length);
    let part_totals = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            if forward_only {
                for_each_canonical_kmer_2bit::<true, _>(
                    seq,
                    start,
                    end,
                    word_length,
                    soft_mask,
                    |key, _, _| {
                        counts[(key >> sub_bits) as usize] += 1;
                    },
                );
            } else {
                for_each_canonical_kmer_2bit::<false, _>(
                    seq,
                    start,
                    end,
                    word_length,
                    soft_mask,
                    |key, _, _| {
                        counts[(key >> sub_bits) as usize] += 1;
                    },
                );
            }
            counts
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

    let mut slabs = Vec::new();
    let mut p_start = 0usize;
    let mut acc = 0usize;
    for (p, &cnt) in part_totals.iter().enumerate() {
        if cnt == 0 {
            continue;
        }
        if acc > 0 && acc + cnt > cap {
            slabs.push((p_start as u16, p as u16));
            p_start = p;
            acc = cnt;
        } else {
            acc += cnt;
        }
    }
    if p_start < NUM_PARTITIONS {
        slabs.push((p_start as u16, NUM_PARTITIONS as u16));
    }
    slabs
}

impl RefIndex {
    /// Estimate peak RAM in bytes for a known `total_kmers` count.
    fn estimate_bytes_for_kmers(
        reference: &Genome,
        total_kmers: usize,
        word_length: usize,
        spaced_mask: Option<u64>,
        compare_reverse: bool,
    ) -> usize {
        if total_kmers >= u32::MAX as usize {
            return usize::MAX;
        }
        let ps_bytes = if reference.seq.len < (1usize << 31) {
            std::mem::size_of::<u32>()
        } else {
            std::mem::size_of::<u64>()
        };
        let rev_bytes = if compare_reverse {
            reference.seq.heap_bytes()
        } else {
            0
        };
        let unique_bm_bytes = if total_kmers > 0 {
            reference.seq.len.div_ceil(64) * std::mem::size_of::<u64>()
        } else {
            0
        };
        let threads = rayon::current_num_threads().max(1);

        if use_direct_csr_for(word_length, reference.total_bases, spaced_mask) {
            let sub_keys = 1usize << (2 * word_length - 8);
            let offsets_bytes = NUM_PARTITIONS * (sub_keys + 1) * std::mem::size_of::<u32>();
            let ps_total = total_kmers.saturating_mul(ps_bytes);
            // Single-slab `rev_ref` is built sequentially AFTER `RefIndex::build` drops transient
            // `coarse_keys`, reusing freed heap space.
            let coarse_keys_bytes = total_kmers.saturating_mul(std::mem::size_of::<u32>());
            let scratch_per_part = (total_kmers / NUM_PARTITIONS).max(1) * ps_bytes;
            let scratch_bytes = threads.saturating_mul(scratch_per_part + sub_keys * 4);
            let transient_build = coarse_keys_bytes.saturating_add(scratch_bytes);
            offsets_bytes
                .saturating_add(ps_total)
                .saturating_add(unique_bm_bytes)
                .saturating_add(transient_build.max(rev_bytes))
        } else {
            let p = choose_bucket_bits(total_kmers);
            let p_bits = choose_hash_part_bits(total_kmers);
            let n_coarse = 1usize << p_bits;
            let sub_buckets = 1usize << (p - p_bits);
            let offsets_bytes = n_coarse * (sub_buckets + 1) * std::mem::size_of::<u32>();
            let bloom_bytes = choose_bloom_params(total_kmers, NUM_PARTITIONS, 1)
                .map_or(0, |(sb, _)| {
                    (NUM_PARTITIONS << sb) * std::mem::size_of::<BloomBlock>()
                });
            let slot_cap = choose_hash_slot_cap(total_kmers, n_coarse);
            let alloc_total = n_coarse.saturating_mul(slot_cap);
            // In-place split-tag build: `tags` (4 B) + `pos_strand` (4/8 B) are allocated once
            // and rewritten in-place; transient `coarse_mid` (3 B) + Pass-2 `scratch` are
            // `munmap`ped at the end of Pass 2 before `rev_ref` is built sequentially.
            let tags_bytes = alloc_total
                .saturating_add(8)
                .saturating_mul(std::mem::size_of::<u32>());
            let ps_total = alloc_total.saturating_mul(ps_bytes);
            let coarse_mid_bytes = alloc_total.saturating_mul(3).saturating_add(1);
            let scratch_per_part = slot_cap
                .saturating_add(slot_cap >> 3)
                .saturating_add(1024)
                .saturating_mul(std::mem::size_of::<u64>() + ps_bytes);
            let scratch_bytes = threads.saturating_mul(scratch_per_part + sub_buckets * 4);
            let transient_build = coarse_mid_bytes.saturating_add(scratch_bytes);
            offsets_bytes
                .saturating_add(tags_bytes)
                .saturating_add(ps_total)
                .saturating_add(bloom_bytes)
                .saturating_add(unique_bm_bytes)
                .saturating_add(transient_build.max(rev_bytes))
        }
    }

    /// Estimate exact peak RAM in bytes required to build and hold `RefIndex` (plus 2-bit `rev_ref`
    /// if `compare_reverse` is true). Returns `usize::MAX` if the reference has `>= 2^32` valid k-mers.
    pub fn estimate_bytes(
        reference: &Genome,
        word_length: usize,
        soft_mask: bool,
        compare_reverse: bool,
    ) -> usize {
        Self::estimate_bytes_with_mask(reference, word_length, None, soft_mask, compare_reverse)
    }

    pub fn estimate_bytes_with_mask(
        reference: &Genome,
        word_length: usize,
        spaced_mask: Option<u64>,
        soft_mask: bool,
        compare_reverse: bool,
    ) -> usize {
        let total_kmers = count_valid_kmers(reference, word_length, soft_mask);
        Self::estimate_bytes_for_kmers(
            reference,
            total_kmers,
            word_length,
            spaced_mask,
            compare_reverse,
        )
    }

    /// Check whether a single-slab `RefIndex` fits within `memory_budget`.
    pub fn fits_in_budget(
        reference: &Genome,
        word_length: usize,
        soft_mask: bool,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> bool {
        Self::fits_in_budget_with_mask(
            reference,
            word_length,
            None,
            soft_mask,
            compare_reverse,
            memory_budget,
        )
    }

    pub fn fits_in_budget_with_mask(
        reference: &Genome,
        word_length: usize,
        spaced_mask: Option<u64>,
        soft_mask: bool,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> bool {
        Self::estimate_bytes_with_mask(
            reference,
            word_length,
            spaced_mask,
            soft_mask,
            compare_reverse,
        ) <= memory_budget
    }

    /// Compute the maximum number of valid k-mers per reference slab that fits comfortably inside
    /// `memory_budget` (accounting for 2-bit `rev_ref`, Hash/Direct CSR arrays, and build scratch).
    /// Returns `None` if `memory_budget` is too small for efficient slab streaming (< 128 MiB net).
    pub fn max_kmers_per_slab_for_budget(
        reference: &Genome,
        word_length: usize,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> Option<usize> {
        Self::max_kmers_per_slab_for_budget_with_mask(
            reference,
            word_length,
            None,
            compare_reverse,
            memory_budget,
        )
    }

    pub fn max_kmers_per_slab_for_budget_with_mask(
        reference: &Genome,
        word_length: usize,
        spaced_mask: Option<u64>,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> Option<usize> {
        const MIN_SLAB_BUDGET: usize = 128 * 1024 * 1024;
        let rev_bytes = if compare_reverse {
            reference.seq.heap_bytes()
        } else {
            0
        };
        let slab_budget = memory_budget.saturating_sub(rev_bytes);
        if slab_budget < MIN_SLAB_BUDGET {
            return None;
        }
        let ps_bytes = if reference.seq.len < (1usize << 31) {
            4usize
        } else {
            8usize
        };
        // Peak build RAM per k-mer (multi-slab passes `build_unique_bm = false`, so `unique_bm`
        // is 0 B): Split-tag `tags` (4 B) + transient `coarse_mid` (3 B) + `pos_strand` (4/8 B) +
        // scale-invariant Bloom (~0.66-1.0 B/k-mer) + offsets/slack headroom (~2 B).
        let bytes_per_kmer = if use_direct_csr_for(word_length, reference.total_bases, spaced_mask)
        {
            ps_bytes + 5
        } else {
            10 + ps_bytes
        };
        let usable = slab_budget * 7 / 10;
        let max_kmers = (usable / bytes_per_kmer).min((u32::MAX as usize) - 1);
        if max_kmers < 1_000_000 {
            None
        } else {
            Some(max_kmers)
        }
    }

    /// Build a canonical reference index directly at `word_length` (`4..=32`) over all of `reference`.
    pub fn build(
        reference: &Genome,
        word_length: usize,
        freq_limit: u64,
        soft_mask: bool,
    ) -> Result<Self> {
        Self::build_with_mask(reference, word_length, None, freq_limit, soft_mask)
    }

    /// Build a canonical reference index with an optional 2-bit palindromic spaced-seed mask.
    pub fn build_with_mask(
        reference: &Genome,
        word_length: usize,
        spaced_mask: Option<u64>,
        freq_limit: u64,
        soft_mask: bool,
    ) -> Result<Self> {
        Self::build_with_mask_and_strand(
            reference,
            word_length,
            spaced_mask,
            freq_limit,
            soft_mask,
            false,
        )
    }

    /// Build a canonical (`forward_only = false`) or forward-only (`forward_only = true`)
    /// reference index with an optional 2-bit palindromic spaced-seed mask.
    pub fn build_with_mask_and_strand(
        reference: &Genome,
        word_length: usize,
        spaced_mask: Option<u64>,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
    ) -> Result<Self> {
        Self::build_partitions_with_mask_and_strand(
            reference,
            0,
            NUM_PARTITIONS as u16,
            word_length,
            spaced_mask,
            freq_limit,
            soft_mask,
            forward_only,
        )
    }

    /// Build a canonical reference index covering all valid k-mers of `reference` whose coarse
    /// partition falls in `[part_start .. part_end)` (within `0..256`), guaranteeing globally exact
    /// per-strand `freq_limit` filtering across multi-slab runs.
    pub fn build_partitions(
        reference: &Genome,
        part_start: u16,
        part_end: u16,
        word_length: usize,
        freq_limit: u64,
        soft_mask: bool,
    ) -> Result<Self> {
        Self::build_partitions_with_mask(
            reference,
            part_start,
            part_end,
            word_length,
            None,
            freq_limit,
            soft_mask,
        )
    }

    /// Build a partition-slab canonical reference index with an optional 2-bit spaced-seed mask.
    #[allow(clippy::too_many_arguments)]
    pub fn build_partitions_with_mask(
        reference: &Genome,
        part_start: u16,
        part_end: u16,
        word_length: usize,
        spaced_mask: Option<u64>,
        freq_limit: u64,
        soft_mask: bool,
    ) -> Result<Self> {
        Self::build_partitions_with_mask_and_strand(
            reference,
            part_start,
            part_end,
            word_length,
            spaced_mask,
            freq_limit,
            soft_mask,
            false,
        )
    }

    /// Build a partition-slab canonical or forward-only reference index with an optional 2-bit
    /// spaced-seed mask.
    #[allow(clippy::too_many_arguments)]
    pub fn build_partitions_with_mask_and_strand(
        reference: &Genome,
        part_start: u16,
        part_end: u16,
        word_length: usize,
        spaced_mask: Option<u64>,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
    ) -> Result<Self> {
        Self::build_partitions_with_options(
            reference,
            part_start,
            part_end,
            word_length,
            spaced_mask,
            freq_limit,
            soft_mask,
            forward_only,
            true,
        )
    }

    /// Internal partition-slab builder with configurable `unique_bm` construction
    /// (`build_unique_bm = false` avoids rebuilding a full-genome uniqueness bitmap on every
    /// slab in multi-slab streaming while remaining 100% bit-exact).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build_partitions_with_options(
        reference: &Genome,
        part_start: u16,
        part_end: u16,
        word_length: usize,
        spaced_mask: Option<u64>,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
        build_unique_bm: bool,
    ) -> Result<Self> {
        let max_starts = reference
            .seq
            .len
            .saturating_sub(word_length.saturating_sub(1));
        Self::build_range_partitions(
            reference,
            0,
            max_starts,
            part_start as usize,
            part_end as usize,
            word_length,
            spaced_mask,
            freq_limit,
            soft_mask,
            forward_only,
            build_unique_bm,
        )
    }

    /// Build a canonical reference index covering valid k-mers whose start positions lie in
    /// `[slab_start .. slab_end)`, preserving exact global reference coordinates.
    pub fn build_range(
        reference: &Genome,
        slab_start: usize,
        slab_end: usize,
        word_length: usize,
        freq_limit: u64,
        soft_mask: bool,
    ) -> Result<Self> {
        Self::build_range_partitions(
            reference,
            slab_start,
            slab_end,
            0,
            NUM_PARTITIONS,
            word_length,
            None,
            freq_limit,
            soft_mask,
            false,
            true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build_range_partitions(
        reference: &Genome,
        slab_start: usize,
        slab_end: usize,
        part_start: usize,
        part_end: usize,
        word_length: usize,
        spaced_mask: Option<u64>,
        freq_limit: u64,
        soft_mask: bool,
        forward_only: bool,
        build_unique_bm: bool,
    ) -> Result<Self> {
        if !(4..=32).contains(&word_length) {
            return Err(GeckoError::InvalidParameter(format!(
                "word_length must be between 4 and 32, got {word_length}"
            )));
        }
        if part_start >= part_end || part_end > NUM_PARTITIONS {
            return Err(GeckoError::InvalidParameter(format!(
                "invalid partition range [{part_start}..{part_end})"
            )));
        }

        let use_u32_ps = reference.seq.len < (1usize << 31);
        if use_direct_csr_for(word_length, reference.total_bases, spaced_mask) {
            if use_u32_ps {
                build_direct_csr_u32(
                    reference,
                    slab_start,
                    slab_end,
                    part_start,
                    part_end,
                    word_length,
                    freq_limit,
                    soft_mask,
                    forward_only,
                    build_unique_bm,
                )
            } else {
                build_direct_csr_u64(
                    reference,
                    slab_start,
                    slab_end,
                    part_start,
                    part_end,
                    word_length,
                    freq_limit,
                    soft_mask,
                    forward_only,
                    build_unique_bm,
                )
            }
        } else if use_u32_ps {
            build_hash_csr_u32(
                reference,
                slab_start,
                slab_end,
                part_start,
                part_end,
                word_length,
                spaced_mask,
                freq_limit,
                soft_mask,
                forward_only,
                build_unique_bm,
            )
        } else {
            build_hash_csr_u64(
                reference,
                slab_start,
                slab_end,
                part_start,
                part_end,
                word_length,
                spaced_mask,
                freq_limit,
                soft_mask,
                forward_only,
                build_unique_bm,
            )
        }
    }

    /// Approximate heap memory currently held by the built index.
    pub fn heap_bytes(&self) -> usize {
        let unique_bytes = self.unique_bm.capacity() * std::mem::size_of::<u64>();
        match &self.kind {
            RefIndexKind::DirectCsr {
                offsets,
                pos_strand,
                ..
            } => {
                offsets.capacity() * std::mem::size_of::<u32>()
                    + pos_strand.heap_bytes()
                    + unique_bytes
            }
            RefIndexKind::HashCsr {
                offsets,
                tags,
                pos_strand,
                bloom,
                ..
            } => {
                offsets.capacity() * std::mem::size_of::<u32>()
                    + tags.capacity() * std::mem::size_of::<u32>()
                    + pos_strand.heap_bytes()
                    + bloom.as_ref().map_or(0, BlockedBloom::heap_bytes)
                    + unique_bytes
            }
        }
    }

    /// Whether this index uses Direct CSR mode (zero key storage).
    #[inline(always)]
    pub fn is_direct(&self) -> bool {
        matches!(self.kind, RefIndexKind::DirectCsr { .. })
    }

    /// Whether this index was built with `--forward-only` specialization.
    #[inline(always)]
    pub fn is_forward_only(&self) -> bool {
        self.forward_only
    }

    #[inline(always)]
    pub(crate) fn pos_strand_vec(&self) -> &PosStrandVec {
        match &self.kind {
            RefIndexKind::DirectCsr { pos_strand, .. } => pos_strand,
            RefIndexKind::HashCsr { pos_strand, .. } => pos_strand,
        }
    }

    /// Compute the local `offsets` index `b_idx` for a canonical k-mer key, or `None` if the key
    /// belongs to a partition outside `[self.part_start .. self.part_end)`.
    #[inline(always)]
    pub(crate) fn bucket_index(&self, canon_key: u64) -> Option<usize> {
        let ps = self.part_start as usize;
        let pe = self.part_end as usize;
        match &self.kind {
            RefIndexKind::DirectCsr { sub_bits, .. } => {
                let b = canon_key as usize;
                let p = b >> *sub_bits;
                if p < ps || p >= pe {
                    return None;
                }
                let b_local = b - (ps << *sub_bits);
                Some(b_local + (b_local >> *sub_bits))
            }
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                ..
            } => {
                let h = fmix64(canon_key);
                let p = (h >> 56) as usize;
                if p < ps || p >= pe {
                    return None;
                }
                let b = (h >> *hash_shift) as usize;
                let b_local = b - (ps << (56 - *hash_shift));
                Some(b_local + (b_local >> *sub_bits))
            }
        }
    }

    /// Stage 2 load: return `[start, end)` bucket range for `b_idx`.
    #[inline(always)]
    pub(crate) fn bucket_range(&self, b_idx: usize) -> (usize, usize) {
        let offsets = match &self.kind {
            RefIndexKind::DirectCsr { offsets, .. } => offsets,
            RefIndexKind::HashCsr { offsets, .. } => offsets,
        };
        // Safety: `b_idx` is in `0 .. offsets.len() - 1` by construction of `bucket_index`.
        unsafe {
            (
                *offsets.get_unchecked(b_idx) as usize,
                *offsets.get_unchecked(b_idx + 1) as usize,
            )
        }
    }

    /// Stage 3 probe: given bucket range `[start, end)` and `canon_key`, return the sub-range
    /// `[lo, hi)` of matching entries in `pos_strand` (and reference to `pos_strand`).
    #[inline(always)]
    pub(crate) fn probe_bucket(
        &self,
        start: usize,
        end: usize,
        canon_key: u64,
        use_avx2: bool,
    ) -> (usize, usize, &PosStrandVec) {
        match &self.kind {
            RefIndexKind::DirectCsr { pos_strand, .. } => (start, end, pos_strand),
            RefIndexKind::HashCsr {
                tags, pos_strand, ..
            } => {
                let tag = fmix64(canon_key) as u32;
                let (lo, hi) = probe_sorted_tags(tags, start, end, tag, use_avx2);
                (lo, hi, pos_strand)
            }
        }
    }

    /// Convenience single-key lookup returning all matching `(pos_y, strand_bit)` entries.
    pub fn lookup(&self, canon_key: u64) -> Vec<(u64, u8)> {
        if let RefIndexKind::HashCsr {
            bloom: Some(bloom), ..
        } = &self.kind
        {
            let h = fmix64(canon_key);
            let p = (h >> 56) as usize;
            if p < self.part_start as usize || p >= self.part_end as usize {
                return Vec::new();
            }
            let blk = ((h >> bloom.block_shift) as usize)
                - ((self.part_start as usize) << bloom.sub_bits);
            if !unsafe { bloom.contains_unchecked(blk, h as u32) } {
                return Vec::new();
            }
        }
        let Some(b_idx) = self.bucket_index(canon_key) else {
            return Vec::new();
        };
        let (start, end) = self.bucket_range(b_idx);
        if start == end {
            return Vec::new();
        }
        #[cfg(target_arch = "x86_64")]
        let use_avx2 = std::is_x86_feature_detected!("avx2");
        #[cfg(not(target_arch = "x86_64"))]
        let use_avx2 = false;
        let (lo, hi, ps) = self.probe_bucket(start, end, canon_key, use_avx2);
        (lo..hi).map(|i| ps.get_unpacked(i)).collect()
    }
}

/// Probe a sorted `u32` tag bucket `tags[start..end]` for exact occurrences of `tag`.
///
/// Requires `tags` to have at least 8 trailing padding elements past `end` when `use_avx2` is true
/// (guaranteed by `build_hash_csr`, which allocates `total_kmers + 8` initialized elements).
#[inline(always)]
pub(crate) fn probe_sorted_tags(
    tags: &[u32],
    start: usize,
    end: usize,
    tag: u32,
    use_avx2: bool,
) -> (usize, usize) {
    let len = end - start;
    if len == 0 {
        return (start, start);
    }
    if len <= 32 {
        #[cfg(target_arch = "x86_64")]
        if use_avx2 && start + 8 <= tags.len() {
            // Safety: `use_avx2` verified by `is_x86_feature_detected!("avx2")` and `tags` padded.
            return unsafe { probe_bucket_u32_avx2(tags, start, end, tag) };
        }
        let _ = use_avx2;
        return probe_bucket_linear_u32(tags, start, end, tag);
    }

    probe_bucket_binary_u32(tags, start, end, tag)
}

/// Scalar sorted `u32` tag bucket probe (`<= 16` linear with early-exit, `> 16` binary search).
#[inline(always)]
pub(crate) fn probe_sorted_tags_scalar(
    tags: &[u32],
    start: usize,
    end: usize,
    tag: u32,
) -> (usize, usize) {
    let len = end - start;
    if len <= 16 {
        probe_bucket_linear_u32(tags, start, end, tag)
    } else {
        probe_bucket_binary_u32(tags, start, end, tag)
    }
}

#[inline(always)]
fn probe_bucket_linear_u32(tags: &[u32], start: usize, end: usize, tag: u32) -> (usize, usize) {
    let ptr = tags.as_ptr();
    let mut i = start;
    while i < end {
        let t = unsafe { *ptr.add(i) };
        if t >= tag {
            if t > tag {
                return (start, start);
            }
            let lo = i;
            i += 1;
            while i < end && unsafe { *ptr.add(i) } == tag {
                i += 1;
            }
            return (lo, i);
        }
        i += 1;
    }
    (start, start)
}

#[inline(always)]
pub(crate) fn probe_bucket_binary_u32(
    tags: &[u32],
    start: usize,
    end: usize,
    tag: u32,
) -> (usize, usize) {
    let ptr = tags.as_ptr();
    let mut lo = start;
    let mut hi = end;
    while lo < hi {
        let mid = lo + ((hi - lo) >> 1);
        if unsafe { *ptr.add(mid) } < tag {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo >= end || unsafe { *ptr.add(lo) } != tag {
        return (start, start);
    }
    let mut r_hi = lo + 1;
    while r_hi < end && unsafe { *ptr.add(r_hi) } == tag {
        r_hi += 1;
    }
    (lo, r_hi)
}

/// Inlineable AVX2 8-wide sorted `u32` tag probe.
///
/// For the common case (`len <= 8`, ~98% of all buckets), executes a single branchless 256-bit
/// load + `vpcmpeqd` + `vmovmskps` + length mask.
///
/// # Safety
/// Caller must run in a `#[target_feature(enable = "avx2")]` context and guarantee
/// `start <= end`, `end + 8 <= tags.len()` (or `start + 8 <= tags.len()` when `len <= 8`).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
pub(crate) unsafe fn probe_sorted_tags_avx2_inline(
    tags: &[u32],
    start: usize,
    end: usize,
    tag: u32,
) -> (usize, usize) {
    use core::arch::x86_64::*;

    let len = end - start;
    let ptr = tags.as_ptr();
    if len <= 8 {
        let vtag = _mm256_set1_epi32(tag as i32);
        let vcand = unsafe { _mm256_loadu_si256(ptr.add(start) as *const __m256i) };
        let cmp = _mm256_cmpeq_epi32(vcand, vtag);
        let valid_mask = (1u32 << len) - 1;
        let mask = (_mm256_movemask_ps(_mm256_castsi256_ps(cmp)) as u32) & valid_mask;
        if mask == 0 {
            return (start, start);
        }
        let lo = start + mask.trailing_zeros() as usize;
        let hi = start + (32 - mask.leading_zeros()) as usize;
        return (lo, hi);
    }
    if len <= 32 {
        let vtag = _mm256_set1_epi32(tag as i32);
        let mut j = start;
        while j + 8 <= end {
            let vcand = unsafe { _mm256_loadu_si256(ptr.add(j) as *const __m256i) };
            let cmp = _mm256_cmpeq_epi32(vcand, vtag);
            let mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp)) as u32;
            if mask != 0 {
                let lo = j + mask.trailing_zeros() as usize;
                if (mask & 0x80) == 0 {
                    let hi = j + (32 - mask.leading_zeros()) as usize;
                    return (lo, hi);
                }
                let mut hi = j + 8;
                while hi < end && unsafe { *ptr.add(hi) } == tag {
                    hi += 1;
                }
                return (lo, hi);
            }
            if unsafe { *ptr.add(j + 7) } > tag {
                return (start, start);
            }
            j += 8;
        }
        let rem = end - j;
        if rem > 0 {
            let vcand = unsafe { _mm256_loadu_si256(ptr.add(j) as *const __m256i) };
            let cmp = _mm256_cmpeq_epi32(vcand, vtag);
            let valid_mask = (1u32 << rem) - 1;
            let mask = (_mm256_movemask_ps(_mm256_castsi256_ps(cmp)) as u32) & valid_mask;
            if mask != 0 {
                let lo = j + mask.trailing_zeros() as usize;
                let hi = j + (32 - mask.leading_zeros()) as usize;
                return (lo, hi);
            }
        }
        return (start, start);
    }

    probe_bucket_binary_u32(tags, start, end, tag)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn probe_bucket_u32_avx2(
    tags: &[u32],
    start: usize,
    end: usize,
    tag: u32,
) -> (usize, usize) {
    unsafe { probe_sorted_tags_avx2_inline(tags, start, end, tag) }
}

// ---------------------------------------------------------------------------
// Parallel chunking & pointer helpers
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct SendMutPtr<T>(*mut T);

unsafe impl<T: Send> Send for SendMutPtr<T> {}
unsafe impl<T: Send> Sync for SendMutPtr<T> {}

/// Partition valid start positions `[slab_start .. slab_end)` into balanced chunks for Pass 1.
fn make_chunk_ranges(
    slab_start: usize,
    slab_end: usize,
    seq_len: usize,
    wl: usize,
) -> Vec<(usize, usize)> {
    let max_starts = seq_len.saturating_sub(wl - 1);
    let lo = slab_start.min(max_starts);
    let hi = slab_end.min(max_starts);
    if lo >= hi {
        return Vec::new();
    }
    let span = hi - lo;
    let threads = rayon::current_num_threads().max(1);
    let target_chunks = (threads * 4).clamp(1, 64);
    let chunk_len = span.div_ceil(target_chunks).max(65_536);

    let mut ranges = Vec::with_capacity(target_chunks + 1);
    let mut s = lo;
    while s < hi {
        let e = (s + chunk_len).min(hi);
        ranges.push((s, e));
        s = e;
    }
    ranges
}

/// Split `slice` into disjoint mutable sub-slices given `offsets` of length `N + 1`.
fn split_by_offsets_mut<'a, T>(slice: &'a mut [T], offsets: &[usize]) -> Vec<&'a mut [T]> {
    let n = offsets.len() - 1;
    let mut out = Vec::with_capacity(n);
    let mut rest = slice;
    for i in 0..n {
        let len = offsets[i + 1] - offsets[i];
        let (head, tail) = rest.split_at_mut(len);
        out.push(head);
        rest = tail;
    }
    out
}

/// Compute global partition offsets, per-partition totals, and per-chunk partition write offsets
/// for `n_parts` active partitions (`n_parts <= NUM_PARTITIONS`).
fn compute_partition_layout(
    chunk_counts: &[[usize; NUM_PARTITIONS]],
    n_parts: usize,
) -> (Vec<usize>, Vec<usize>, Vec<[usize; NUM_PARTITIONS]>, usize) {
    let mut part_totals = vec![0usize; n_parts];
    for counts in chunk_counts {
        for p in 0..n_parts {
            part_totals[p] += counts[p];
        }
    }

    let mut part_offsets = vec![0usize; n_parts + 1];
    for p in 0..n_parts {
        part_offsets[p + 1] = part_offsets[p] + part_totals[p];
    }
    let total = part_offsets[n_parts];

    let mut chunk_write = vec![[0usize; NUM_PARTITIONS]; chunk_counts.len()];
    let mut running = [0usize; NUM_PARTITIONS];
    running[..n_parts].copy_from_slice(&part_offsets[..n_parts]);
    for (c, counts) in chunk_counts.iter().enumerate() {
        chunk_write[c][..n_parts].copy_from_slice(&running[..n_parts]);
        for p in 0..n_parts {
            running[p] += counts[p];
        }
    }

    (part_offsets, part_totals, chunk_write, total)
}

#[inline(always)]
fn kmer_mask(wl: usize) -> u64 {
    if wl == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * wl)) - 1
    }
}

/// Iterate all valid canonical (or forward-only) k-mers whose start position lies in
/// `[start..end)` directly from 2-bit packed `u64` words.
#[inline(always)]
fn for_each_canonical_kmer_2bit<const FWD_ONLY: bool, F>(
    seq: &PackedSeq,
    start: usize,
    end: usize,
    wl: usize,
    soft_mask: bool,
    f: F,
) where
    F: FnMut(u64, u64, u64),
{
    if wl == 32 {
        for_each_canonical_kmer_2bit_impl::<true, FWD_ONLY, F>(seq, start, end, 32, soft_mask, f);
    } else {
        for_each_canonical_kmer_2bit_impl::<false, FWD_ONLY, F>(seq, start, end, wl, soft_mask, f);
    }
}

#[inline(always)]
fn for_each_canonical_kmer_2bit_opt<const FWD_ONLY: bool, F>(
    seq: &PackedSeq,
    start: usize,
    end: usize,
    wl: usize,
    spaced_mask: Option<u64>,
    soft_mask: bool,
    f: F,
) where
    F: FnMut(u64, u64, u64),
{
    if let Some(sm) = spaced_mask {
        for_each_canonical_spaced_kmer_2bit::<FWD_ONLY, F>(seq, start, end, wl, sm, soft_mask, f);
    } else {
        for_each_canonical_kmer_2bit::<FWD_ONLY, F>(seq, start, end, wl, soft_mask, f);
    }
}

/// Iterate canonical (or forward-only) spaced-seed k-mers using a palindromic 2-bit `spaced_mask`.
#[inline(always)]
fn for_each_canonical_spaced_kmer_2bit<const FWD_ONLY: bool, F>(
    seq: &PackedSeq,
    start: usize,
    end: usize,
    wl: usize,
    spaced_mask: u64,
    soft_mask: bool,
    mut f: F,
) where
    F: FnMut(u64, u64, u64),
{
    let span_mask = if wl == 32 { u64::MAX } else { kmer_mask(wl) };
    let rc_shift = (2 * (wl - 1)) as u32;
    let scan_end = (end + wl - 1).min(seq.len);

    seq.for_each_valid_run(start, scan_end, soft_mask, |rs, re| {
        if re - rs < wl {
            return;
        }
        let mut kf = 0u64;
        let mut kr = 0u64;
        let mut idx = rs;
        let prime_end = rs + wl - 1;
        while idx < prime_end {
            let w_idx = idx >> 5;
            let lane = idx & 31;
            let avail = (32 - lane).min(prime_end - idx);
            let mut word = unsafe { *seq.bits.get_unchecked(w_idx) } >> (lane << 1);
            for _ in 0..avail {
                let b = word & 3;
                word >>= 2;
                kf = ((kf << 2) | b) & span_mask;
                if !FWD_ONLY {
                    kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                }
            }
            idx += avail;
        }

        while idx < re {
            let w_idx = idx >> 5;
            let lane = idx & 31;
            let avail = (32 - lane).min(re - idx);
            let mut word = unsafe { *seq.bits.get_unchecked(w_idx) } >> (lane << 1);
            let start_pos = (idx + 1 - wl) as u64;
            for pos in start_pos..start_pos + avail as u64 {
                let b = word & 3;
                word >>= 2;
                kf = ((kf << 2) | b) & span_mask;
                if FWD_ONLY {
                    f(kf & spaced_mask, pos, 0);
                } else {
                    kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                    let sf = kf & spaced_mask;
                    let sr = kr & spaced_mask;
                    if sr < sf {
                        f(sr, pos, 1);
                    } else {
                        f(sf, pos, if sf == sr { 2 } else { 0 });
                    }
                }
            }
            idx += avail;
        }
    });
}

#[inline(always)]
fn for_each_canonical_kmer_2bit_impl<const WL32: bool, const FWD_ONLY: bool, F>(
    seq: &PackedSeq,
    start: usize,
    end: usize,
    wl: usize,
    soft_mask: bool,
    mut f: F,
) where
    F: FnMut(u64, u64, u64),
{
    let mask = if WL32 { u64::MAX } else { kmer_mask(wl) };
    let rc_shift = if WL32 { 62u32 } else { (2 * (wl - 1)) as u32 };
    let scan_end = (end + wl - 1).min(seq.len);

    seq.for_each_valid_run(start, scan_end, soft_mask, |rs, re| {
        if re - rs < wl {
            return;
        }
        // Prime first `wl - 1` bases of the clean run.
        let mut kf = 0u64;
        let mut kr = 0u64;
        let mut idx = rs;
        let prime_end = rs + wl - 1;
        while idx < prime_end {
            let w_idx = idx >> 5;
            let lane = idx & 31;
            let avail = (32 - lane).min(prime_end - idx);
            let mut word = unsafe { *seq.bits.get_unchecked(w_idx) } >> (lane << 1);
            for _ in 0..avail {
                let b = word & 3;
                word >>= 2;
                kf = if WL32 {
                    (kf << 2) | b
                } else {
                    ((kf << 2) | b) & mask
                };
                if !FWD_ONLY {
                    kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                }
            }
            idx += avail;
        }

        // Emit one canonical (or forward-only) k-mer per base through `re`.
        while idx < re {
            let w_idx = idx >> 5;
            let lane = idx & 31;
            let avail = (32 - lane).min(re - idx);
            let mut word = unsafe { *seq.bits.get_unchecked(w_idx) } >> (lane << 1);
            let start_pos = (idx + 1 - wl) as u64;
            for pos in start_pos..start_pos + avail as u64 {
                let b = word & 3;
                word >>= 2;
                kf = if WL32 {
                    (kf << 2) | b
                } else {
                    ((kf << 2) | b) & mask
                };
                if FWD_ONLY {
                    f(kf, pos, 0);
                } else {
                    kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                    if kr < kf {
                        f(kr, pos, 1);
                    } else {
                        f(kf, pos, if kf == kr { 2 } else { 0 });
                    }
                }
            }
            idx += avail;
        }
    });
}

// ---------------------------------------------------------------------------
// Lock-free prefetched `unique_bm` Pass 2 helpers & branchless sorting networks
// ---------------------------------------------------------------------------

/// Set all bits in half-open `[start, end)` inside `bm`.
fn set_bit_range(bm: &mut [u64], start: usize, end: usize) {
    if start >= end {
        return;
    }
    let w_s = start >> 6;
    let w_e = (end - 1) >> 6;
    let bit_s = start & 63;
    let bit_e = (end - 1) & 63;
    if w_s == w_e {
        bm[w_s] |= (u64::MAX << bit_s) & (u64::MAX >> (63 - bit_e));
    } else {
        bm[w_s] |= u64::MAX << bit_s;
        bm[w_s + 1..w_e].fill(u64::MAX);
        bm[w_e] |= u64::MAX >> (63 - bit_e);
    }
}

/// Initialize the inverted `non_unique_bm` (`1` = invalid / soft-masked / out-of-range) in
/// O(#runs) sequential time so Pass 2 only marks duplicate keys (`cnt >= 2`).
fn init_non_unique_bm(
    seq: &PackedSeq,
    slab_start: usize,
    slab_end: usize,
    wl: usize,
    soft_mask: bool,
) -> Vec<u64> {
    let mut bm = vec![0u64; seq.len.div_ceil(64)];
    let max_starts = seq.len.saturating_sub(wl - 1);
    let lo = slab_start.min(max_starts);
    let hi = slab_end.min(max_starts);
    set_bit_range(&mut bm, 0, lo);
    set_bit_range(&mut bm, hi, seq.len);
    for &(s, e) in &seq.special_runs {
        let rs = (s as usize).saturating_sub(wl - 1);
        let re = (e as usize).min(seq.len);
        set_bit_range(&mut bm, rs, re);
    }
    if soft_mask {
        for &(s, e) in &seq.soft_runs {
            let rs = (s as usize).saturating_sub(wl - 1);
            let re = (e as usize).min(seq.len);
            set_bit_range(&mut bm, rs, re);
        }
    }
    bm
}

#[cold]
#[inline(never)]
fn flush_unique_buf(buf: &[usize], unique_atomic: &[AtomicU64]) {
    let ptr = unique_atomic.as_ptr();
    for &pos_y in buf {
        prefetch_t2(unsafe { ptr.add(pos_y >> 6) });
    }
    for &pos_y in buf {
        unsafe { unique_atomic.get_unchecked(pos_y >> 6) }
            .fetch_or(1u64 << (pos_y & 63), Ordering::Relaxed);
    }
}

#[inline(always)]
fn push_unique_pos(
    buf: &mut [usize; UNIQ_BUF_CAP],
    len: &mut usize,
    pos_y: usize,
    unique_atomic: &[AtomicU64],
) {
    if unique_atomic.is_empty() {
        return;
    }
    unsafe { *buf.get_unchecked_mut(*len) = pos_y };
    *len += 1;
    if *len == UNIQ_BUF_CAP {
        flush_unique_buf(buf.as_slice(), unique_atomic);
        *len = 0;
    }
}

/// Compile-time-specialized Blocked Bloom insert (`K_H == 0` folds away completely;
/// `K_H in 1..=6` unrolls straight-line with immediate `BLOOM_SALTS[i]`, zero jump tables,
/// zero salt spills).
#[inline(always)]
unsafe fn insert_bloom_const<const K_H: usize>(
    bloom_ptr: *mut BloomBlock,
    h_rot: u64,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
) {
    if K_H == 0 {
        return;
    }
    let sub_blk = (((h_rot as u32) >> bloom_sub_shift) as usize) & bloom_sub_mask;
    let words = unsafe { &mut (*bloom_ptr.add(sub_blk)).0 };
    let tag = (h_rot >> 32) as u32;
    let mut i = 0usize;
    while i < K_H {
        let h = tag.wrapping_mul(BLOOM_SALTS[i]);
        let word = (h >> 29) as usize;
        let bit = (h >> 23) & 63;
        unsafe {
            *words.get_unchecked_mut(word) |= 1u64 << bit;
        }
        i += 1;
    }
}

// ---------------------------------------------------------------------------
// Direct CSR builder (generic over u32 / u64 packed pos_strand)
// ---------------------------------------------------------------------------

trait PackedPs: Copy + Default + Ord + std::fmt::Debug + Send + Sync + 'static {
    fn pack(pos: u64, strand: u64) -> Self;
    fn pos(self) -> usize;
    fn strand_bit(self) -> u64;
    fn to_u64(self) -> u64;
    fn into_vec(v: Vec<Self>) -> PosStrandVec;
}

impl PackedPs for u32 {
    #[inline(always)]
    fn pack(pos: u64, strand: u64) -> Self {
        ((pos << 1) | (strand & 1)) as u32
    }

    #[inline(always)]
    fn pos(self) -> usize {
        (self >> 1) as usize
    }

    #[inline(always)]
    fn strand_bit(self) -> u64 {
        (self & 1) as u64
    }

    #[inline(always)]
    fn to_u64(self) -> u64 {
        self as u64
    }

    fn into_vec(v: Vec<Self>) -> PosStrandVec {
        PosStrandVec::U32(v)
    }
}

impl PackedPs for u64 {
    #[inline(always)]
    fn pack(pos: u64, strand: u64) -> Self {
        (pos << 1) | (strand & 1)
    }

    #[inline(always)]
    fn pos(self) -> usize {
        (self >> 1) as usize
    }

    #[inline(always)]
    fn strand_bit(self) -> u64 {
        self & 1
    }

    #[inline(always)]
    fn to_u64(self) -> u64 {
        self
    }

    fn into_vec(v: Vec<Self>) -> PosStrandVec {
        PosStrandVec::U64(v)
    }
}

/// Branchless compare-and-swap of paired `(key, ps)` values using 128-bit comparison
/// (`cmp + sbb` + `cmov` on x86_64, zero conditional branches).
#[inline(always)]
fn cas_pair<P: PackedPs>(ka: &mut u64, pa: &mut P, kb: &mut u64, pb: &mut P) {
    let ak = *ka;
    let bk = *kb;
    let ap = *pa;
    let bp = *pb;
    let a128 = ((ak as u128) << 64) | (ap.to_u64() as u128);
    let b128 = ((bk as u128) << 64) | (bp.to_u64() as u128);
    let swap = a128 > b128;
    *ka = if swap { bk } else { ak };
    *kb = if swap { ak } else { bk };
    *pa = if swap { bp } else { ap };
    *pb = if swap { ap } else { bp };
}

/// Unrolled branchless bubble-insert of element at `POS` into already-sorted `[0..POS)`.
/// Uses only 4 live registers and zero inner-loop branches.
#[inline(always)]
unsafe fn bubble_insert_at<const POS: usize, P: PackedPs>(sk: *mut u64, sp: *mut P) {
    let mut k = unsafe { *sk.add(POS) };
    let mut p = unsafe { *sp.add(POS) };
    let mut j = POS;
    while j > 0 {
        let mut pk = unsafe { *sk.add(j - 1) };
        let mut pp = unsafe { *sp.add(j - 1) };
        cas_pair(&mut pk, &mut pp, &mut k, &mut p);
        unsafe {
            *sk.add(j) = k;
            *sp.add(j) = p;
        }
        k = pk;
        p = pp;
        j -= 1;
    }
    unsafe {
        *sk = k;
        *sp = p;
    }
}

/// Branchless cumulative compare-and-swap sort for `2 <= cnt <= 6` in-place on `(sk, sp)`.
/// Structurally impossible to lower to a jump table (`ind_jmp = 0`, `Bim = 0`), half the
/// static `cas_pair` count of independent networks (15 vs 30), and peak 8 live registers.
#[inline(always)]
unsafe fn sort_small_2_to_6<P: PackedPs>(sk: *mut u64, sp: *mut P, cnt: usize) {
    debug_assert!((2..=6).contains(&cnt));
    let mut k0 = unsafe { *sk };
    let mut p0 = unsafe { *sp };
    let mut k1 = unsafe { *sk.add(1) };
    let mut p1 = unsafe { *sp.add(1) };
    cas_pair(&mut k0, &mut p0, &mut k1, &mut p1);
    if cnt == 2 {
        unsafe {
            *sk = k0;
            *sp = p0;
            *sk.add(1) = k1;
            *sp.add(1) = p1;
        }
        return;
    }

    let mut k2 = unsafe { *sk.add(2) };
    let mut p2 = unsafe { *sp.add(2) };
    cas_pair(&mut k1, &mut p1, &mut k2, &mut p2);
    cas_pair(&mut k0, &mut p0, &mut k1, &mut p1);
    if cnt == 3 {
        unsafe {
            *sk = k0;
            *sp = p0;
            *sk.add(1) = k1;
            *sp.add(1) = p1;
            *sk.add(2) = k2;
            *sp.add(2) = p2;
        }
        return;
    }

    let mut k3 = unsafe { *sk.add(3) };
    let mut p3 = unsafe { *sp.add(3) };
    cas_pair(&mut k2, &mut p2, &mut k3, &mut p3);
    cas_pair(&mut k1, &mut p1, &mut k2, &mut p2);
    cas_pair(&mut k0, &mut p0, &mut k1, &mut p1);
    unsafe {
        *sk = k0;
        *sp = p0;
        *sk.add(1) = k1;
        *sp.add(1) = p1;
        *sk.add(2) = k2;
        *sp.add(2) = p2;
        *sk.add(3) = k3;
        *sp.add(3) = p3;
    }
    if cnt == 4 {
        return;
    }
    unsafe { bubble_insert_at::<4, P>(sk, sp) };
    if cnt == 6 {
        unsafe { bubble_insert_at::<5, P>(sk, sp) };
    }
}

#[cold]
#[inline(never)]
fn sort_large_bucket<P: PackedPs>(sk: &mut [u64], sp: &mut [P], sort_buf_k: &mut Vec<(u64, P)>) {
    let cnt = sk.len();
    sort_buf_k.clear();
    sort_buf_k.reserve(cnt);
    for idx in 0..cnt {
        unsafe {
            sort_buf_k.push((*sk.get_unchecked(idx), *sp.get_unchecked(idx)));
        }
    }
    sort_buf_k.sort_unstable();
    for (off, &(h_rot, ps)) in sort_buf_k.iter().enumerate() {
        unsafe {
            *sk.get_unchecked_mut(off) = h_rot;
            *sp.get_unchecked_mut(off) = ps;
        }
    }
}

#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
unsafe fn emit_filtered_bucket<P: PackedPs, const K_H: usize>(
    scratch_k: &[u64],
    scratch_p: &[P],
    sub_start: usize,
    sub_end: usize,
    freq_limit: u64,
    forward_only: bool,
    tag_ptr: *mut u32,
    ps_ptr: *mut P,
    write_idx: &mut usize,
    bloom_ptr: *mut BloomBlock,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    dup_buf: &mut [usize; UNIQ_BUF_CAP],
    dup_len: &mut usize,
    non_unique_atomic: &[AtomicU64],
) {
    if forward_only {
        // Forward-only exact 64-bit freq_limit filtering (all entries have strand 0).
        let mut r = sub_start;
        while r < sub_end {
            let h_rot = unsafe { *scratch_k.get_unchecked(r) };
            let mut r_end = r + 1;
            while r_end < sub_end && unsafe { *scratch_k.get_unchecked(r_end) } == h_rot {
                r_end += 1;
            }
            let run_len = r_end - r;
            if run_len > 1 {
                for &ps in unsafe { scratch_p.get_unchecked(r..r_end) } {
                    push_unique_pos(dup_buf, dup_len, ps.pos(), non_unique_atomic);
                }
            }
            if (run_len as u64) <= freq_limit {
                unsafe {
                    insert_bloom_const::<K_H>(bloom_ptr, h_rot, bloom_sub_shift, bloom_sub_mask);
                }
                let tag = (h_rot >> 32) as u32;
                for &ps in unsafe { scratch_p.get_unchecked(r..r_end) } {
                    unsafe {
                        *tag_ptr.add(*write_idx) = tag;
                        *ps_ptr.add(*write_idx) = ps;
                    }
                    *write_idx += 1;
                }
            }
            r = r_end;
        }
    } else {
        // Canonical exact 64-bit per-strand freq_limit filtering across runs of identical `h_rot`.
        let mut r = sub_start;
        while r < sub_end {
            let h_rot = unsafe { *scratch_k.get_unchecked(r) };
            let mut r_end = r + 1;
            let mut cnt0 = 0u64;
            let mut cnt1 = 0u64;
            if unsafe { scratch_p.get_unchecked(r).strand_bit() } == 0 {
                cnt0 += 1;
            } else {
                cnt1 += 1;
            }
            while r_end < sub_end && unsafe { *scratch_k.get_unchecked(r_end) } == h_rot {
                if unsafe { scratch_p.get_unchecked(r_end).strand_bit() } == 0 {
                    cnt0 += 1;
                } else {
                    cnt1 += 1;
                }
                r_end += 1;
            }
            if r_end - r > 1 {
                for &ps in unsafe { scratch_p.get_unchecked(r..r_end) } {
                    push_unique_pos(dup_buf, dup_len, ps.pos(), non_unique_atomic);
                }
            }
            let keep0 = cnt0 <= freq_limit;
            let keep1 = cnt1 <= freq_limit;
            if keep0 || keep1 {
                unsafe {
                    insert_bloom_const::<K_H>(bloom_ptr, h_rot, bloom_sub_shift, bloom_sub_mask);
                }
                let tag = (h_rot >> 32) as u32;
                if keep0 && keep1 {
                    for &ps in unsafe { scratch_p.get_unchecked(r..r_end) } {
                        unsafe {
                            *tag_ptr.add(*write_idx) = tag;
                            *ps_ptr.add(*write_idx) = ps;
                        }
                        *write_idx += 1;
                    }
                } else {
                    for &ps in unsafe { scratch_p.get_unchecked(r..r_end) } {
                        let is_keep = if ps.strand_bit() == 0 { keep0 } else { keep1 };
                        if is_keep {
                            unsafe {
                                *tag_ptr.add(*write_idx) = tag;
                                *ps_ptr.add(*write_idx) = ps;
                            }
                            *write_idx += 1;
                        }
                    }
                }
            }
            r = r_end;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn build_direct_csr_u32(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
) -> Result<RefIndex> {
    build_direct_csr::<u32>(
        reference,
        slab_start,
        slab_end,
        part_start,
        part_end,
        wl,
        freq_limit,
        soft_mask,
        forward_only,
        build_unique_bm,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_direct_csr_u64(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
) -> Result<RefIndex> {
    build_direct_csr::<u64>(
        reference,
        slab_start,
        slab_end,
        part_start,
        part_end,
        wl,
        freq_limit,
        soft_mask,
        forward_only,
        build_unique_bm,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_direct_csr<P: PackedPs>(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
) -> Result<RefIndex> {
    let seq = &reference.seq;
    let sub_bits = (2 * wl - 8) as u32;
    let sub_keys = 1usize << sub_bits;
    let sub_mask = sub_keys - 1;
    let n_parts = part_end - part_start;
    let full_parts = part_start == 0 && part_end == NUM_PARTITIONS;

    let ranges = make_chunk_ranges(slab_start, slab_end, seq.len, wl);
    if ranges.is_empty() {
        return Ok(RefIndex {
            word_length: wl,
            freq_limit,
            soft_mask,
            spaced_mask: None,
            ref_data_len: seq.len,
            total_kmers: 0,
            kept_entries: 0,
            part_start: part_start as u16,
            part_end: part_end as u16,
            kind: RefIndexKind::DirectCsr {
                sub_bits,
                offsets: vec![0u32; n_parts * (sub_keys + 1)],
                pos_strand: P::into_vec(Vec::new()),
            },
            unique_bm: Vec::new(),
            unique_bm_xor: 0,
            forward_only,
        });
    }

    // Pass 1a: parallel count into active coarse partitions `[part_start .. part_end)`.
    let chunk_counts: Vec<[usize; NUM_PARTITIONS]> = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            macro_rules! count_direct {
                ($fwd:expr) => {
                    if full_parts {
                        for_each_canonical_kmer_2bit::<$fwd, _>(
                            seq,
                            start,
                            end,
                            wl,
                            soft_mask,
                            |canon, _, _| {
                                counts[(canon >> sub_bits) as usize] += 1;
                            },
                        );
                    } else {
                        for_each_canonical_kmer_2bit::<$fwd, _>(
                            seq,
                            start,
                            end,
                            wl,
                            soft_mask,
                            |canon, _, _| {
                                let p = (canon >> sub_bits) as usize;
                                if p >= part_start && p < part_end {
                                    counts[p - part_start] += 1;
                                }
                            },
                        );
                    }
                };
            }
            if forward_only {
                count_direct!(true);
            } else {
                count_direct!(false);
            }
            counts
        })
        .collect();

    let (part_offsets, _part_totals, chunk_write, total_kmers) =
        compute_partition_layout(&chunk_counts, n_parts);
    if total_kmers >= u32::MAX as usize {
        return Err(GeckoError::InvalidParameter(
            "reference slab has >= 2^32 valid k-mers; use smaller --memory budget".into(),
        ));
    }

    // Pass 1b: parallel scatter with L2-resident per-thread staging buffers (`STAGE_CAP = 64`).
    let mut coarse_keys: Vec<u32> = Vec::with_capacity(total_kmers);
    let mut pos_strand: Vec<P> = Vec::with_capacity(total_kmers);
    let keys_ptr = SendMutPtr(coarse_keys.as_mut_ptr());
    let ps_ptr = SendMutPtr(pos_strand.as_mut_ptr());

    ranges.par_iter().zip(chunk_write.par_iter()).for_each_init(
        || {
            let buf_k: Box<[[u32; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[0u32; STAGE_CAP]; NUM_PARTITIONS]
                    .into_boxed_slice()
                    .try_into()
                    .unwrap();
            let buf_p: Box<[[P; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[P::default(); STAGE_CAP]; NUM_PARTITIONS]
                    .into_boxed_slice()
                    .try_into()
                    .unwrap();
            (buf_k, buf_p)
        },
        |(buf_k, buf_p), (&(start, end), &init_write)| {
            let kptr = keys_ptr;
            let pptr = ps_ptr;
            let mut write_pos = init_write;
            let mut buf_len = [0u8; NUM_PARTITIONS];

            macro_rules! scatter_direct {
                ($fwd:expr) => {
                    for_each_canonical_kmer_2bit::<$fwd, _>(
                        seq,
                        start,
                        end,
                        wl,
                        soft_mask,
                        |canon, pos, strand| {
                            let p_glob = (canon >> sub_bits) as usize;
                            if !full_parts && (p_glob < part_start || p_glob >= part_end) {
                                return;
                            }
                            let p = if full_parts {
                                p_glob
                            } else {
                                p_glob - part_start
                            };
                            unsafe {
                                let idx = *buf_len.get_unchecked(p) as usize;
                                *buf_k.get_unchecked_mut(p).get_unchecked_mut(idx) =
                                    (canon as usize & sub_mask) as u32;
                                *buf_p.get_unchecked_mut(p).get_unchecked_mut(idx) =
                                    P::pack(pos, strand);
                                let next_len = idx + 1;
                                if next_len == STAGE_CAP {
                                    let dst = *write_pos.get_unchecked(p);
                                    std::ptr::copy_nonoverlapping(
                                        buf_k.get_unchecked(p).as_ptr(),
                                        kptr.0.add(dst),
                                        STAGE_CAP,
                                    );
                                    std::ptr::copy_nonoverlapping(
                                        buf_p.get_unchecked(p).as_ptr(),
                                        pptr.0.add(dst),
                                        STAGE_CAP,
                                    );
                                    *write_pos.get_unchecked_mut(p) = dst + STAGE_CAP;
                                    *buf_len.get_unchecked_mut(p) = 0;
                                } else {
                                    *buf_len.get_unchecked_mut(p) = next_len as u8;
                                }
                            }
                        },
                    );
                };
            }
            if forward_only {
                scatter_direct!(true);
            } else {
                scatter_direct!(false);
            }

            for p in 0..n_parts {
                let rem = unsafe { *buf_len.get_unchecked(p) } as usize;
                if rem > 0 {
                    let dst = unsafe { *write_pos.get_unchecked(p) };
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            buf_k.get_unchecked(p).as_ptr(),
                            kptr.0.add(dst),
                            rem,
                        );
                        std::ptr::copy_nonoverlapping(
                            buf_p.get_unchecked(p).as_ptr(),
                            pptr.0.add(dst),
                            rem,
                        );
                    }
                }
            }
        },
    );

    // Safety: Pass 1b wrote every element in `0..total_kmers` across disjoint partition chunks.
    unsafe {
        coarse_keys.set_len(total_kmers);
        pos_strand.set_len(total_kmers);
    }

    // Pass 2: Rayon-parallel fine-bucketing + per-strand freq_limit filtering + in-place compaction
    // + lock-free prefetched `unique_bm` construction.
    let off_len = n_parts * (sub_keys + 1);
    let mut offsets: Vec<u32> = Vec::with_capacity(off_len);
    let mut unique_bm = if build_unique_bm && total_kmers > 0 && freq_limit >= 1 {
        vec![0u64; seq.len.div_ceil(64)]
    } else {
        Vec::new()
    };
    let unique_atomic: &[AtomicU64] = AtomicU64::from_mut_slice(&mut unique_bm);
    let need_pal_check = !forward_only && can_key_be_palindrome(wl, None);
    let mask = kmer_mask(wl);

    let part_k_slices = split_by_offsets_mut(&mut coarse_keys, &part_offsets);
    let part_p_slices = split_by_offsets_mut(&mut pos_strand, &part_offsets);
    let part_off_slices: Vec<&mut [MaybeUninit<u32>]> = offsets.spare_capacity_mut()[..off_len]
        .chunks_exact_mut(sub_keys + 1)
        .collect();

    let kept_entries: usize = (0..n_parts)
        .into_par_iter()
        .zip(part_k_slices.into_par_iter())
        .zip(part_p_slices.into_par_iter())
        .zip(part_off_slices.into_par_iter())
        .map_init(
            || (vec![0u32; sub_keys], Vec::<P>::new()),
            |(sub_off, scratch_ps), (((lp, pk), pp), poff)| {
                let part_off = part_offsets[lp];
                let n = pk.len();
                if n == 0 || freq_limit == 0 {
                    for slot in poff.iter_mut() {
                        slot.write(part_off as u32);
                    }
                    return 0usize;
                }

                sub_off.fill(0);
                for &sub in pk.iter() {
                    unsafe {
                        *sub_off.get_unchecked_mut(sub as usize) += 1;
                    }
                }

                // In-place exclusive prefix sum; after scatter below, `sub_off[s] == end(s)`.
                let mut acc = 0u32;
                for s in 0..sub_keys {
                    unsafe {
                        let c = *sub_off.get_unchecked(s);
                        *sub_off.get_unchecked_mut(s) = acc;
                        acc += c;
                    }
                }

                if scratch_ps.len() < n {
                    scratch_ps.resize(n, P::default());
                }
                for i in 0..n {
                    unsafe {
                        let sub = *pk.get_unchecked(i) as usize;
                        let dst = *sub_off.get_unchecked(sub) as usize;
                        *scratch_ps.get_unchecked_mut(dst) = *pp.get_unchecked(i);
                        *sub_off.get_unchecked_mut(sub) = (dst + 1) as u32;
                    }
                }

                let mut uniq_buf = [0usize; UNIQ_BUF_CAP];
                let mut uniq_len = 0usize;
                let mut write_idx = 0usize;
                let mut sub_start = 0usize;
                for (s, (slot, &sub_end_u32)) in poff.iter_mut().zip(sub_off.iter()).enumerate() {
                    slot.write((part_off + write_idx) as u32);
                    let sub_end = sub_end_u32 as usize;
                    let cnt = sub_end - sub_start;
                    if cnt == 0 {
                        continue;
                    }
                    if cnt == 1 {
                        let ps = unsafe { *scratch_ps.get_unchecked(sub_start) };
                        unsafe { *pp.get_unchecked_mut(write_idx) = ps };
                        write_idx += 1;
                        let canon = (((part_start + lp) as u64) << sub_bits) | (s as u64);
                        if !need_pal_check || !is_palindrome_key(canon, wl, mask) {
                            push_unique_pos(&mut uniq_buf, &mut uniq_len, ps.pos(), unique_atomic);
                        }
                    } else if (cnt as u64) <= freq_limit {
                        unsafe {
                            std::ptr::copy_nonoverlapping(
                                scratch_ps.as_ptr().add(sub_start),
                                pp.as_mut_ptr().add(write_idx),
                                cnt,
                            );
                        }
                        write_idx += cnt;
                    } else if !forward_only {
                        let mut cnt0 = 0u64;
                        let mut cnt1 = 0u64;
                        for &ps in &scratch_ps[sub_start..sub_end] {
                            if ps.strand_bit() == 0 {
                                cnt0 += 1;
                            } else {
                                cnt1 += 1;
                            }
                        }
                        let keep0 = cnt0 <= freq_limit;
                        let keep1 = cnt1 <= freq_limit;
                        if keep0 && keep1 {
                            unsafe {
                                std::ptr::copy_nonoverlapping(
                                    scratch_ps.as_ptr().add(sub_start),
                                    pp.as_mut_ptr().add(write_idx),
                                    cnt,
                                );
                            }
                            write_idx += cnt;
                        } else if keep0 || keep1 {
                            for &ps in &scratch_ps[sub_start..sub_end] {
                                let is_keep = if ps.strand_bit() == 0 { keep0 } else { keep1 };
                                if is_keep {
                                    unsafe { *pp.get_unchecked_mut(write_idx) = ps };
                                    write_idx += 1;
                                }
                            }
                        }
                    }
                    sub_start = sub_end;
                }
                if uniq_len > 0 {
                    flush_unique_buf(&uniq_buf[..uniq_len], unique_atomic);
                }
                poff[sub_keys].write((part_off + write_idx) as u32);
                write_idx
            },
        )
        .sum();

    // Safety: every active partition wrote all `sub_keys + 1` elements of its `poff` slice.
    unsafe { offsets.set_len(off_len) };

    drop(coarse_keys);

    Ok(RefIndex {
        word_length: wl,
        freq_limit,
        soft_mask,
        spaced_mask: None,
        ref_data_len: seq.len,
        total_kmers,
        kept_entries,
        part_start: part_start as u16,
        part_end: part_end as u16,
        kind: RefIndexKind::DirectCsr {
            sub_bits,
            offsets,
            pos_strand: P::into_vec(pos_strand),
        },
        unique_bm,
        unique_bm_xor: 0,
        forward_only,
    })
}

// ---------------------------------------------------------------------------
// 2-Level Parallel Radix-Partitioned Hash-CSR builder (single-pass split-tag)
// ---------------------------------------------------------------------------

#[repr(align(64))]
struct PaddedAtomicUsize(AtomicUsize);

#[derive(Clone, Copy)]
struct OverflowEntry<P: PackedPs> {
    lp: u16,
    tag: u32,
    mid: [u8; 3],
    ps: P,
}

#[derive(Clone, Copy)]
struct SendConstPtr<T>(*const T);
unsafe impl<T> Send for SendConstPtr<T> {}
unsafe impl<T> Sync for SendConstPtr<T> {}

type OverflowAppend<P> = (usize, Vec<u32>, Vec<P>);

/// Single-pass atomic split-tag scatter into final `tags`, final `pos_strand`, and contiguous
/// transient `coarse_mid` (`3 B/k-mer`), spilling rare partition overflows into lock-free
/// thread-local `Vec<OverflowEntry<P>>` buffers.
#[allow(clippy::too_many_arguments)]
fn scatter_hash_single_pass<P: PackedPs, const STAGE_CAP_H: usize, const STAGE_MID_BYTES: usize>(
    seq: &PackedSeq,
    ranges: &[(usize, usize)],
    wl: usize,
    spaced_mask: Option<u64>,
    soft_mask: bool,
    forward_only: bool,
    full_parts: bool,
    part_shift: u32,
    ip_start: usize,
    ip_end: usize,
    n_coarse: usize,
    slot_cap: usize,
    tptr: SendMutPtr<u32>,
    mptr: SendMutPtr<u8>,
    pptr: SendMutPtr<P>,
    part_cursors: &[PaddedAtomicUsize],
    non_unique_atomic: &[AtomicU64],
) -> Vec<Vec<OverflowEntry<P>>> {
    debug_assert_eq!(STAGE_MID_BYTES, STAGE_CAP_H * 3 + 1);
    let track_non_unique = !non_unique_atomic.is_empty();

    #[cold]
    #[inline(never)]
    fn push_overflow<P: PackedPs>(
        ov: &mut Vec<OverflowEntry<P>>,
        lp: usize,
        tags: &[u32],
        mids: &[u8],
        ps: &[P],
    ) {
        let lp_u16 = lp as u16;
        for i in 0..tags.len() {
            ov.push(OverflowEntry {
                lp: lp_u16,
                tag: tags[i],
                mid: [mids[i * 3], mids[i * 3 + 1], mids[i * 3 + 2]],
                ps: ps[i],
            });
        }
    }

    ranges
        .par_iter()
        .map_init(
            || {
                let buf_t = vec![[0u32; STAGE_CAP_H]; n_coarse];
                let buf_m = vec![[0u8; STAGE_MID_BYTES]; n_coarse];
                let buf_p = vec![[P::default(); STAGE_CAP_H]; n_coarse];
                (buf_t, buf_m, buf_p)
            },
            |(buf_t, buf_m, buf_p), &(start, end)| {
                let (tp, mp, pp) = (tptr, mptr, pptr);
                let mut buf_len = [0u8; MAX_HASH_COARSE_PARTITIONS];
                let mut local_ov: Vec<OverflowEntry<P>> = Vec::new();
                let mut out_word_idx = usize::MAX;
                let mut out_word_mask = 0u64;

                let flush_part = |lp: usize,
                                  len: usize,
                                  bt: &[u32; STAGE_CAP_H],
                                  bm: &[u8; STAGE_MID_BYTES],
                                  bp: &[P; STAGE_CAP_H],
                                  ov: &mut Vec<OverflowEntry<P>>| {
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
                        push_overflow(
                            ov,
                            lp,
                            &bt[n_flat..len],
                            &bm[n_flat * 3..len * 3],
                            &bp[n_flat..len],
                        );
                    }
                };

                macro_rules! scatter_hash {
                    ($fwd:expr) => {
                        for_each_canonical_kmer_2bit_opt::<$fwd, _>(
                            seq,
                            start,
                            end,
                            wl,
                            spaced_mask,
                            soft_mask,
                            |canon, pos, strand| {
                                let h = fmix64(canon);
                                let ip = (h >> part_shift) as usize;
                                if !full_parts && (ip < ip_start || ip >= ip_end) {
                                    if track_non_unique {
                                        let w = (pos as usize) >> 6;
                                        if w != out_word_idx {
                                            if out_word_idx != usize::MAX {
                                                unsafe {
                                                    non_unique_atomic.get_unchecked(out_word_idx)
                                                }
                                                .fetch_or(out_word_mask, Ordering::Relaxed);
                                            }
                                            out_word_idx = w;
                                            out_word_mask = 1u64 << (pos & 63);
                                        } else {
                                            out_word_mask |= 1u64 << (pos & 63);
                                        }
                                    }
                                    return;
                                }
                                if !$fwd && strand == 2 && track_non_unique {
                                    unsafe { non_unique_atomic.get_unchecked((pos as usize) >> 6) }
                                        .fetch_or(1u64 << (pos & 63), Ordering::Relaxed);
                                }
                                let lp = if full_parts { ip } else { ip - ip_start };
                                unsafe {
                                    let idx = *buf_len.get_unchecked(lp) as usize;
                                    let bt = buf_t.get_unchecked_mut(lp);
                                    let bm = buf_m.get_unchecked_mut(lp);
                                    let bp = buf_p.get_unchecked_mut(lp);
                                    *bt.get_unchecked_mut(idx) = h as u32;
                                    let mid = (h >> 32) as u32;
                                    (bm.as_mut_ptr().add(idx * 3) as *mut u32)
                                        .write_unaligned(mid.to_le());
                                    *bp.get_unchecked_mut(idx) = P::pack(pos, strand);
                                    let next_len = idx + 1;
                                    if next_len == STAGE_CAP_H {
                                        flush_part(lp, STAGE_CAP_H, bt, bm, bp, &mut local_ov);
                                        *buf_len.get_unchecked_mut(lp) = 0;
                                    } else {
                                        *buf_len.get_unchecked_mut(lp) = next_len as u8;
                                    }
                                }
                            },
                        );
                    };
                }
                if forward_only {
                    scatter_hash!(true);
                } else {
                    scatter_hash!(false);
                }

                if !full_parts && out_word_idx != usize::MAX {
                    unsafe { non_unique_atomic.get_unchecked(out_word_idx) }
                        .fetch_or(out_word_mask, Ordering::Relaxed);
                }

                for lp in 0..n_coarse {
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
        .collect()
}

/// Shared Pass-2 sub-bucket sort, `freq_limit` filter, Bloom insert, and duplicate-only
/// `non_unique_bm` compaction loop writing into `(tag_ptr, ps_ptr)`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn compact_hash_buckets<P: PackedPs, const K_H: usize>(
    part_off: usize,
    sub_buckets: usize,
    freq_limit: u64,
    forward_only: bool,
    bloom_ptr: *mut BloomBlock,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    non_unique_atomic: &[AtomicU64],
    sub_off: &[u32],
    poff: &mut [MaybeUninit<u32>],
    tag_ptr: *mut u32,
    ps_ptr: *mut P,
    scratch_k: &mut [u64],
    scratch_p: &mut [P],
    sort_buf_k: &mut Vec<(u64, P)>,
) -> usize {
    let sk_base = scratch_k.as_mut_ptr();
    let sp_base = scratch_p.as_mut_ptr();
    let mut dup_buf = [0usize; UNIQ_BUF_CAP];
    let mut dup_len = 0usize;
    let mut write_idx = 0usize;
    let mut sub_start = 0usize;

    for (slot, &sub_end_u32) in poff.iter_mut().zip(sub_off.iter()) {
        slot.write((part_off + write_idx) as u32);
        let sub_end = sub_end_u32 as usize;
        let cnt = sub_end - sub_start;
        if cnt == 0 {
            continue;
        }
        if cnt == 1 {
            let h_rot = unsafe { *sk_base.add(sub_start) };
            let ps = unsafe { *sp_base.add(sub_start) };
            unsafe {
                *tag_ptr.add(write_idx) = (h_rot >> 32) as u32;
                *ps_ptr.add(write_idx) = ps;
                insert_bloom_const::<K_H>(bloom_ptr, h_rot, bloom_sub_shift, bloom_sub_mask);
            }
            write_idx += 1;
            sub_start = sub_end;
            continue;
        }

        if cnt <= 6 {
            unsafe { sort_small_2_to_6(sk_base.add(sub_start), sp_base.add(sub_start), cnt) };
        } else if cnt <= 24 {
            insertion_sort_pair(
                unsafe { scratch_k.get_unchecked_mut(sub_start..sub_end) },
                unsafe { scratch_p.get_unchecked_mut(sub_start..sub_end) },
            );
        } else {
            sort_large_bucket(
                unsafe { scratch_k.get_unchecked_mut(sub_start..sub_end) },
                unsafe { scratch_p.get_unchecked_mut(sub_start..sub_end) },
                sort_buf_k,
            );
        }

        if (cnt as u64) <= freq_limit {
            let first_k = unsafe { *sk_base.add(sub_start) };
            let mut prev_k = !first_k; // Guaranteed != first_k for all u64!
            let mut prev_ps = P::default();
            let mut prev_marked = false;
            for idx in sub_start..sub_end {
                let h_rot = unsafe { *sk_base.add(idx) };
                let ps = unsafe { *sp_base.add(idx) };
                unsafe {
                    *tag_ptr.add(write_idx) = (h_rot >> 32) as u32;
                    *ps_ptr.add(write_idx) = ps;
                }
                write_idx += 1;
                if h_rot != prev_k {
                    unsafe {
                        insert_bloom_const::<K_H>(
                            bloom_ptr,
                            h_rot,
                            bloom_sub_shift,
                            bloom_sub_mask,
                        );
                    }
                    prev_k = h_rot;
                    prev_ps = ps;
                    prev_marked = false;
                } else {
                    if !prev_marked {
                        push_unique_pos(
                            &mut dup_buf,
                            &mut dup_len,
                            prev_ps.pos(),
                            non_unique_atomic,
                        );
                        prev_marked = true;
                    }
                    push_unique_pos(&mut dup_buf, &mut dup_len, ps.pos(), non_unique_atomic);
                }
            }
        } else {
            unsafe {
                emit_filtered_bucket::<P, K_H>(
                    scratch_k,
                    scratch_p,
                    sub_start,
                    sub_end,
                    freq_limit,
                    forward_only,
                    tag_ptr,
                    ps_ptr,
                    &mut write_idx,
                    bloom_ptr,
                    bloom_sub_shift,
                    bloom_sub_mask,
                    &mut dup_buf,
                    &mut dup_len,
                    non_unique_atomic,
                );
            }
        }

        sub_start = sub_end;
    }

    if dup_len > 0 {
        flush_unique_buf(&dup_buf[..dup_len], non_unique_atomic);
    }
    poff[sub_buckets].write((part_off + write_idx) as u32);
    write_idx
}

#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
unsafe fn process_hash_partition_overflow<P: PackedPs, const K_H: usize>(
    lp: usize,
    part_off: usize,
    flat_len: usize,
    slot_cap: usize,
    pmid: &[u8],
    pp: &mut [MaybeUninit<P>],
    ptag: &mut [MaybeUninit<u32>],
    poff: &mut [MaybeUninit<u32>],
    ov: &[OverflowEntry<P>],
    bloom_ptr: *mut BloomBlock,
    sub_buckets: usize,
    sub_shift_hi: u32,
    sub_mask: usize,
    freq_limit: u64,
    forward_only: bool,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    non_unique_atomic: &[AtomicU64],
    sub_off: &mut [u32],
    scratch_k: &mut Vec<u64>,
    scratch_p: &mut Vec<P>,
    sort_buf_k: &mut Vec<(u64, P)>,
) -> (usize, Option<OverflowAppend<P>>) {
    if freq_limit == 0 {
        for slot in poff.iter_mut() {
            slot.write(part_off as u32);
        }
        for slot in &mut ptag[..slot_cap] {
            slot.write(0);
        }
        for slot in &mut pp[..slot_cap] {
            slot.write(P::default());
        }
        return (0, None);
    }

    let n_total = flat_len + ov.len();
    sub_off.fill(0);
    let m_ptr = pmid.as_ptr();
    for i in 0..flat_len {
        let mid24 = u32::from_le(unsafe { (m_ptr.add(i * 3) as *const u32).read_unaligned() })
            & 0x00FF_FFFF;
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        unsafe { *sub_off.get_unchecked_mut(sub) += 1 };
    }
    for e in ov {
        let mid24 = (e.mid[0] as u32) | ((e.mid[1] as u32) << 8) | ((e.mid[2] as u32) << 16);
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        unsafe { *sub_off.get_unchecked_mut(sub) += 1 };
    }

    let mut acc = 0u32;
    for s in 0..sub_buckets {
        unsafe {
            let c = *sub_off.get_unchecked(s);
            *sub_off.get_unchecked_mut(s) = acc;
            acc += c;
        }
    }

    if scratch_k.len() < n_total {
        let new_len = n_total + (n_total >> 3) + 1024;
        scratch_k.resize(new_len, 0);
        scratch_p.resize(new_len, P::default());
    }

    let tag_ptr = ptag.as_mut_ptr() as *mut u32;
    let ps_ptr = pp.as_mut_ptr() as *mut P;
    for i in 0..flat_len {
        unsafe {
            let mid24 =
                u32::from_le((m_ptr.add(i * 3) as *const u32).read_unaligned()) & 0x00FF_FFFF;
            let tag = *tag_ptr.add(i);
            let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
            let dst = *sub_off.get_unchecked(sub) as usize;
            *scratch_k.get_unchecked_mut(dst) = ((tag as u64) << 32) | (mid24 as u64);
            *scratch_p.get_unchecked_mut(dst) = *ps_ptr.add(i);
            *sub_off.get_unchecked_mut(sub) = (dst + 1) as u32;
        }
    }
    for e in ov {
        let mid24 = (e.mid[0] as u32) | ((e.mid[1] as u32) << 8) | ((e.mid[2] as u32) << 16);
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        unsafe {
            let dst = *sub_off.get_unchecked(sub) as usize;
            *scratch_k.get_unchecked_mut(dst) = ((e.tag as u64) << 32) | (mid24 as u64);
            *scratch_p.get_unchecked_mut(dst) = e.ps;
            *sub_off.get_unchecked_mut(sub) = (dst + 1) as u32;
        }
    }

    // Compact into temporary Vecs with relative offsets (`part_off = 0`).
    let mut tmp_tags = vec![0u32; n_total];
    let mut tmp_ps = vec![P::default(); n_total];
    let write_idx = unsafe {
        compact_hash_buckets::<P, K_H>(
            0,
            sub_buckets,
            freq_limit,
            forward_only,
            bloom_ptr,
            bloom_sub_shift,
            bloom_sub_mask,
            non_unique_atomic,
            sub_off,
            poff,
            tmp_tags.as_mut_ptr(),
            tmp_ps.as_mut_ptr(),
            scratch_k,
            scratch_p,
            sort_buf_k,
        )
    };

    if write_idx <= slot_cap {
        // Fits back into the partition's in-place slot (common case: high-freq repeats dropped).
        unsafe {
            std::ptr::copy_nonoverlapping(tmp_tags.as_ptr(), tag_ptr, write_idx);
            std::ptr::copy_nonoverlapping(tmp_ps.as_ptr(), ps_ptr, write_idx);
        }
        for slot in &mut ptag[write_idx..slot_cap] {
            slot.write(0);
        }
        for slot in &mut pp[write_idx..slot_cap] {
            slot.write(P::default());
        }
        let base = part_off as u32;
        for slot in poff.iter_mut() {
            unsafe {
                let v = slot.assume_init();
                slot.write(v + base);
            }
        }
        (write_idx, None)
    } else {
        // Adversarial kept count > slot_cap: leave relative offsets in `poff` and append after Pass 2.
        for slot in &mut ptag[..slot_cap] {
            slot.write(0);
        }
        for slot in &mut pp[..slot_cap] {
            slot.write(P::default());
        }
        tmp_tags.truncate(write_idx);
        tmp_ps.truncate(write_idx);
        (write_idx, Some((lp, tmp_tags, tmp_ps)))
    }
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn process_hash_partition_impl<P: PackedPs, const K_H: usize>(
    lp: usize,
    part_off: usize,
    flat_len: usize,
    slot_cap: usize,
    pmid: &[u8],
    pp: &mut [MaybeUninit<P>],
    ptag: &mut [MaybeUninit<u32>],
    poff: &mut [MaybeUninit<u32>],
    ov: &[OverflowEntry<P>],
    bloom_ptr: *mut BloomBlock,
    sub_buckets: usize,
    sub_shift_hi: u32,
    sub_mask: usize,
    freq_limit: u64,
    forward_only: bool,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    non_unique_atomic: &[AtomicU64],
    sub_off: &mut [u32],
    scratch_k: &mut Vec<u64>,
    scratch_p: &mut Vec<P>,
    sort_buf_k: &mut Vec<(u64, P)>,
) -> (usize, Option<OverflowAppend<P>>) {
    if !ov.is_empty() {
        return unsafe {
            process_hash_partition_overflow::<P, K_H>(
                lp,
                part_off,
                flat_len,
                slot_cap,
                pmid,
                pp,
                ptag,
                poff,
                ov,
                bloom_ptr,
                sub_buckets,
                sub_shift_hi,
                sub_mask,
                freq_limit,
                forward_only,
                bloom_sub_shift,
                bloom_sub_mask,
                non_unique_atomic,
                sub_off,
                scratch_k,
                scratch_p,
                sort_buf_k,
            )
        };
    }

    if flat_len == 0 || freq_limit == 0 {
        for slot in poff.iter_mut() {
            slot.write(part_off as u32);
        }
        for slot in &mut ptag[..slot_cap] {
            slot.write(0);
        }
        for slot in &mut pp[..slot_cap] {
            slot.write(P::default());
        }
        return (0usize, None);
    }

    sub_off.fill(0);
    let m_ptr = pmid.as_ptr();
    for i in 0..flat_len {
        let mid24 = u32::from_le(unsafe { (m_ptr.add(i * 3) as *const u32).read_unaligned() })
            & 0x00FF_FFFF;
        let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
        unsafe {
            *sub_off.get_unchecked_mut(sub) += 1;
        }
    }

    // In-place exclusive prefix sum; after scatter below, `sub_off[s] == end(s)`.
    let mut acc = 0u32;
    for s in 0..sub_buckets {
        unsafe {
            let c = *sub_off.get_unchecked(s);
            *sub_off.get_unchecked_mut(s) = acc;
            acc += c;
        }
    }

    if scratch_k.len() < flat_len {
        let new_len = flat_len + (flat_len >> 3) + 1024;
        scratch_k.resize(new_len, 0);
        scratch_p.resize(new_len, P::default());
    }

    let tag_ptr = ptag.as_mut_ptr() as *mut u32;
    let ps_ptr = pp.as_mut_ptr() as *mut P;
    for i in 0..flat_len {
        unsafe {
            let mid24 =
                u32::from_le((m_ptr.add(i * 3) as *const u32).read_unaligned()) & 0x00FF_FFFF;
            let tag = *tag_ptr.add(i);
            let sub = ((mid24 >> sub_shift_hi) as usize) & sub_mask;
            let dst = *sub_off.get_unchecked(sub) as usize;
            *scratch_k.get_unchecked_mut(dst) = ((tag as u64) << 32) | (mid24 as u64);
            *scratch_p.get_unchecked_mut(dst) = *ps_ptr.add(i);
            *sub_off.get_unchecked_mut(sub) = (dst + 1) as u32;
        }
    }

    let write_idx = unsafe {
        compact_hash_buckets::<P, K_H>(
            part_off,
            sub_buckets,
            freq_limit,
            forward_only,
            bloom_ptr,
            bloom_sub_shift,
            bloom_sub_mask,
            non_unique_atomic,
            sub_off,
            poff,
            tag_ptr,
            ps_ptr,
            scratch_k,
            scratch_p,
            sort_buf_k,
        )
    };

    // Zero-initialize any slack entries at the tail of this partition so `tags` and `pos_strand`
    // are 100% initialized across `[0 .. alloc_total)`.
    for slot in &mut ptag[write_idx..slot_cap] {
        slot.write(0);
    }
    for slot in &mut pp[write_idx..slot_cap] {
        slot.write(P::default());
    }
    (write_idx, None)
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn process_hash_partition_scalar<P: PackedPs, const K_H: usize>(
    lp: usize,
    part_off: usize,
    flat_len: usize,
    slot_cap: usize,
    pmid: &[u8],
    pp: &mut [MaybeUninit<P>],
    ptag: &mut [MaybeUninit<u32>],
    poff: &mut [MaybeUninit<u32>],
    ov: &[OverflowEntry<P>],
    bloom_ptr: *mut BloomBlock,
    sub_buckets: usize,
    sub_shift_hi: u32,
    sub_mask: usize,
    freq_limit: u64,
    forward_only: bool,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    non_unique_atomic: &[AtomicU64],
    sub_off: &mut [u32],
    scratch_k: &mut Vec<u64>,
    scratch_p: &mut Vec<P>,
    sort_buf_k: &mut Vec<(u64, P)>,
) -> (usize, Option<OverflowAppend<P>>) {
    unsafe {
        process_hash_partition_impl::<P, K_H>(
            lp,
            part_off,
            flat_len,
            slot_cap,
            pmid,
            pp,
            ptag,
            poff,
            ov,
            bloom_ptr,
            sub_buckets,
            sub_shift_hi,
            sub_mask,
            freq_limit,
            forward_only,
            bloom_sub_shift,
            bloom_sub_mask,
            non_unique_atomic,
            sub_off,
            scratch_k,
            scratch_p,
            sort_buf_k,
        )
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn process_hash_partition_bmi2<P: PackedPs, const K_H: usize>(
    lp: usize,
    part_off: usize,
    flat_len: usize,
    slot_cap: usize,
    pmid: &[u8],
    pp: &mut [MaybeUninit<P>],
    ptag: &mut [MaybeUninit<u32>],
    poff: &mut [MaybeUninit<u32>],
    ov: &[OverflowEntry<P>],
    bloom_ptr: *mut BloomBlock,
    sub_buckets: usize,
    sub_shift_hi: u32,
    sub_mask: usize,
    freq_limit: u64,
    forward_only: bool,
    bloom_sub_shift: u32,
    bloom_sub_mask: usize,
    non_unique_atomic: &[AtomicU64],
    sub_off: &mut [u32],
    scratch_k: &mut Vec<u64>,
    scratch_p: &mut Vec<P>,
    sort_buf_k: &mut Vec<(u64, P)>,
) -> (usize, Option<OverflowAppend<P>>) {
    unsafe {
        process_hash_partition_impl::<P, K_H>(
            lp,
            part_off,
            flat_len,
            slot_cap,
            pmid,
            pp,
            ptag,
            poff,
            ov,
            bloom_ptr,
            sub_buckets,
            sub_shift_hi,
            sub_mask,
            freq_limit,
            forward_only,
            bloom_sub_shift,
            bloom_sub_mask,
            non_unique_atomic,
            sub_off,
            scratch_k,
            scratch_p,
            sort_buf_k,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn build_hash_csr_u32(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    spaced_mask: Option<u64>,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
) -> Result<RefIndex> {
    build_hash_csr::<u32>(
        reference,
        slab_start,
        slab_end,
        part_start,
        part_end,
        wl,
        spaced_mask,
        freq_limit,
        soft_mask,
        forward_only,
        build_unique_bm,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_hash_csr_u64(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    spaced_mask: Option<u64>,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
) -> Result<RefIndex> {
    build_hash_csr::<u64>(
        reference,
        slab_start,
        slab_end,
        part_start,
        part_end,
        wl,
        spaced_mask,
        freq_limit,
        soft_mask,
        forward_only,
        build_unique_bm,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_hash_csr<P: PackedPs>(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    spaced_mask: Option<u64>,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
) -> Result<RefIndex> {
    build_hash_csr_inner::<P>(
        reference,
        slab_start,
        slab_end,
        part_start,
        part_end,
        wl,
        spaced_mask,
        freq_limit,
        soft_mask,
        forward_only,
        build_unique_bm,
        None,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_hash_csr_inner<P: PackedPs>(
    reference: &Genome,
    slab_start: usize,
    slab_end: usize,
    part_start: usize,
    part_end: usize,
    wl: usize,
    spaced_mask: Option<u64>,
    freq_limit: u64,
    soft_mask: bool,
    forward_only: bool,
    build_unique_bm: bool,
    p_bits_override: Option<u32>,
    slot_cap_override: Option<usize>,
) -> Result<RefIndex> {
    let seq = &reference.seq;
    let n_parts = part_end - part_start;
    let full_parts = part_start == 0 && part_end == NUM_PARTITIONS;

    let ranges = make_chunk_ranges(slab_start, slab_end, seq.len, wl);
    let range_kmers = count_valid_kmers_in_range(seq, slab_start, slab_end, wl, soft_mask);
    if ranges.is_empty() || range_kmers == 0 {
        let bucket_bits = 14u32;
        let sub_bits = bucket_bits - 8;
        let sub_buckets = 1usize << sub_bits;
        return Ok(RefIndex {
            word_length: wl,
            freq_limit,
            soft_mask,
            spaced_mask,
            ref_data_len: seq.len,
            total_kmers: 0,
            kept_entries: 0,
            part_start: part_start as u16,
            part_end: part_end as u16,
            kind: RefIndexKind::HashCsr {
                sub_bits,
                hash_shift: 64 - bucket_bits,
                offsets: vec![0u32; n_parts * (sub_buckets + 1)],
                tags: vec![0u32; 8],
                pos_strand: P::into_vec(Vec::new()),
                bloom: None,
            },
            unique_bm: Vec::new(),
            unique_bm_xor: u64::MAX,
            forward_only,
        });
    }

    let bucket_bits = choose_bucket_bits(range_kmers);
    let p_bits = p_bits_override.unwrap_or_else(|| choose_hash_part_bits(range_kmers));
    debug_assert!((8..=10).contains(&p_bits));
    let fine_shift = p_bits - 8;
    let part_shift = 64 - p_bits;
    let ip_start = part_start << fine_shift;
    let ip_end = part_end << fine_shift;
    let n_coarse = ip_end - ip_start;
    debug_assert!(n_coarse <= MAX_HASH_COARSE_PARTITIONS);

    let sub_bits = bucket_bits - p_bits;
    let sub_buckets = 1usize << sub_bits;
    let sub_mask = sub_buckets - 1;
    let hash_shift = 64 - bucket_bits;
    let sub_shift_hi = 32 - bucket_bits;

    let expected_slab_kmers = range_kmers.saturating_mul(n_parts).div_ceil(NUM_PARTITIONS);
    let slot_cap = slot_cap_override
        .unwrap_or_else(|| choose_hash_slot_cap(expected_slab_kmers, n_coarse))
        .max(1);
    let alloc_total = n_coarse.checked_mul(slot_cap).ok_or_else(|| {
        GeckoError::InvalidParameter(
            "reference slab exceeds addressable Hash-CSR capacity; use smaller --memory budget"
                .into(),
        )
    })?;
    if alloc_total > (u32::MAX as usize) - 8 {
        return Err(GeckoError::InvalidParameter(
            "reference slab has >= 2^32 Hash-CSR slots; use smaller --memory budget".into(),
        ));
    }

    // Single-pass atomic split-tag buffers:
    // - `tags` (4 B/slot) & `pos_strand` (4/8 B/slot) are allocated once and compacted in-place.
    // - `coarse_mid` (3 B/slot + 1 sentinel byte) is a single contiguous `mmap`-backed buffer
    //   that is immediately `munmap`ped at the end of Pass 2 (zero glibc arena stranding!).
    let mut tags: Vec<u32> = Vec::with_capacity(alloc_total + 8);
    let mut coarse_mid: Vec<u8> = Vec::with_capacity(alloc_total * 3 + 1);
    let mut pos_strand: Vec<P> = Vec::with_capacity(alloc_total);
    let part_cursors: Vec<PaddedAtomicUsize> = (0..n_coarse)
        .map(|_| PaddedAtomicUsize(AtomicUsize::new(0)))
        .collect();

    // Allocate inverted `non_unique_bm` (`0` = unique valid in-slab non-palindrome, `1` = non-unique)
    // before Pass 1 so Pass 1 marks rare palindromes + out-of-slab positions sequentially, and
    // Pass 2 ONLY writes duplicate keys (`cnt >= 2`) with ZERO `inv_fmix64` calls!
    let mut non_unique_bm = if build_unique_bm && freq_limit >= 1 {
        init_non_unique_bm(seq, slab_start, slab_end, wl, soft_mask)
    } else {
        Vec::new()
    };
    let non_unique_atomic: &[AtomicU64] = AtomicU64::from_mut_slice(&mut non_unique_bm);

    let tptr = SendMutPtr(tags.as_mut_ptr());
    let mptr = SendMutPtr(coarse_mid.as_mut_ptr());
    let pptr = SendMutPtr(pos_strand.as_mut_ptr());

    // Pass 1: single-pass atomic split-tag scatter with adaptive L2-resident staging capacity
    // (`STAGE_CAP = 64 / 32 / 16` for `n_coarse <= 256 / 512 / 1024`).
    let thread_overflows: Vec<Vec<OverflowEntry<P>>> = if n_coarse <= 256 {
        scatter_hash_single_pass::<P, 64, 193>(
            seq,
            &ranges,
            wl,
            spaced_mask,
            soft_mask,
            forward_only,
            full_parts,
            part_shift,
            ip_start,
            ip_end,
            n_coarse,
            slot_cap,
            tptr,
            mptr,
            pptr,
            &part_cursors,
            non_unique_atomic,
        )
    } else if n_coarse <= 512 {
        scatter_hash_single_pass::<P, 32, 97>(
            seq,
            &ranges,
            wl,
            spaced_mask,
            soft_mask,
            forward_only,
            full_parts,
            part_shift,
            ip_start,
            ip_end,
            n_coarse,
            slot_cap,
            tptr,
            mptr,
            pptr,
            &part_cursors,
            non_unique_atomic,
        )
    } else {
        scatter_hash_single_pass::<P, 16, 49>(
            seq,
            &ranges,
            wl,
            spaced_mask,
            soft_mask,
            forward_only,
            full_parts,
            part_shift,
            ip_start,
            ip_end,
            n_coarse,
            slot_cap,
            tptr,
            mptr,
            pptr,
            &part_cursors,
            non_unique_atomic,
        )
    };

    let mut flat_lens = Vec::with_capacity(n_coarse);
    let mut total_kmers = 0usize;
    for cur in &part_cursors {
        let count = cur.0.load(Ordering::Relaxed);
        total_kmers = total_kmers.saturating_add(count);
        flat_lens.push(count.min(slot_cap));
    }
    if total_kmers >= u32::MAX as usize {
        return Err(GeckoError::InvalidParameter(
            "reference slab has >= 2^32 valid k-mers; use smaller --memory budget".into(),
        ));
    }

    // Group rare overflow entries by coarse partition `lp` (empty 99.9% of the time).
    let has_overflow = thread_overflows.iter().any(|v| !v.is_empty());
    let mut part_ov: Vec<Vec<OverflowEntry<P>>> = if has_overflow {
        let mut grouped: Vec<Vec<OverflowEntry<P>>> = (0..n_coarse).map(|_| Vec::new()).collect();
        for chunk_ov in thread_overflows {
            for e in chunk_ov {
                grouped[e.lp as usize].push(e);
            }
        }
        grouped
    } else {
        Vec::new()
    };

    // Write 1 sentinel zero byte past each partition's `flat_len * 3` slice (only when `fl < slot_cap`,
    // so we never clobber byte 0 of partition `lp + 1` when `fl == slot_cap`; if `fl == slot_cap`,
    // byte `(lp + 1) * slot_cap * 3` is already initialized by partition `lp + 1` or by the final
    // sentinel at `alloc_total * 3`). Guarantees every 4-byte unaligned load of `mid24` in Pass 2
    // is 100% in-bounds and initialized (Miri-clean).
    {
        let mid_spare = coarse_mid.spare_capacity_mut();
        for (lp, &fl) in flat_lens.iter().enumerate() {
            if fl < slot_cap {
                mid_spare[(lp * slot_cap + fl) * 3].write(0);
            }
        }
        mid_spare[alloc_total * 3].write(0);
    }

    // Pass 2: Rayon-parallel fine-bucketing, branchless cumulative sorting (`cnt = 2..=6`), exact
    // 64-bit per-strand freq_limit filtering, in-place 32-bit tag + pos_strand compaction,
    // const-specialized partition-owned Blocked Bloom construction (BMI2-dispatched on x86_64),
    // and duplicate-only `non_unique_bm` marking.
    let off_len = n_coarse * (sub_buckets + 1);
    let mut offsets: Vec<u32> = Vec::with_capacity(off_len);

    let bloom_params = choose_bloom_params(total_kmers, n_parts, freq_limit)
        .filter(|&(ext_sb, _)| ext_sb >= fine_shift);
    let (ext_bloom_sb, bloom_k_h_u8) = bloom_params.unwrap_or((0, 0));
    let bloom_k_h = bloom_k_h_u8 as usize;
    debug_assert!(bloom_k_h <= BLOOM_CHOOSE_MAX_K_H);
    let bloom_coarse_sub_bits = ext_bloom_sb.saturating_sub(fine_shift);
    let bloom_sub_shift = 24 - ext_bloom_sb;
    let bloom_sub_mask = (1usize << bloom_coarse_sub_bits) - 1;
    let mut bloom_blocks: Vec<BloomBlock> = if bloom_params.is_some() {
        vec![BloomBlock([0u64; 8]); n_parts << ext_bloom_sb]
    } else {
        Vec::new()
    };

    #[cfg(target_arch = "x86_64")]
    let use_bmi2 = std::arch::is_x86_feature_detected!("bmi2");

    let mid_base = SendConstPtr(coarse_mid.as_ptr());
    let part_p_slices: Vec<&mut [MaybeUninit<P>]> = pos_strand.spare_capacity_mut()[..alloc_total]
        .chunks_exact_mut(slot_cap)
        .collect();
    let part_tag_slices: Vec<&mut [MaybeUninit<u32>]> = tags.spare_capacity_mut()[..alloc_total]
        .chunks_exact_mut(slot_cap)
        .collect();
    let part_off_slices: Vec<&mut [MaybeUninit<u32>]> = offsets.spare_capacity_mut()[..off_len]
        .chunks_exact_mut(sub_buckets + 1)
        .collect();
    let part_bloom_slices: Vec<Option<&mut [BloomBlock]>> = if bloom_params.is_some() {
        bloom_blocks
            .chunks_exact_mut(1usize << bloom_coarse_sub_bits)
            .map(Some)
            .collect()
    } else {
        (0..n_coarse).map(|_| None).collect()
    };

    let part_results: Vec<(usize, Option<OverflowAppend<P>>)> = (0..n_coarse)
        .into_par_iter()
        .zip(part_p_slices.into_par_iter())
        .zip(part_tag_slices.into_par_iter())
        .zip(part_off_slices.into_par_iter())
        .zip(part_bloom_slices.into_par_iter())
        .map_init(
            || {
                (
                    vec![0u32; sub_buckets],
                    Vec::<u64>::new(),
                    Vec::<P>::new(),
                    Vec::<(u64, P)>::new(),
                )
            },
            |(sub_off, scratch_k, scratch_p, sort_buf_k), ((((lp, pp), ptag), poff), pbloom)| {
                let mb = mid_base;
                let part_off = lp * slot_cap;
                let flat_len = flat_lens[lp];
                // Safety: `[part_off * 3 .. (part_off + flat_len) * 3 + 1)` is 100% initialized
                // by Pass 1 + the per-partition sentinel zero bytes, and is read immutably here.
                let pmid =
                    unsafe { std::slice::from_raw_parts(mb.0.add(part_off * 3), flat_len * 3 + 1) };
                let ov: &[OverflowEntry<P>] = if part_ov.is_empty() {
                    &[]
                } else {
                    part_ov[lp].as_slice()
                };
                let bloom_ptr = pbloom.map_or(std::ptr::null_mut(), |s| s.as_mut_ptr());
                macro_rules! dispatch_kh {
                    ($fn_name:ident) => {
                        if bloom_k_h == 4 {
                            $fn_name::<P, 4>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        } else if bloom_k_h == 6 {
                            $fn_name::<P, 6>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        } else if bloom_k_h == 0 {
                            $fn_name::<P, 0>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        } else if bloom_k_h == 5 {
                            $fn_name::<P, 5>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        } else if bloom_k_h == 3 {
                            $fn_name::<P, 3>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        } else if bloom_k_h == 2 {
                            $fn_name::<P, 2>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        } else {
                            $fn_name::<P, 1>(
                                lp,
                                part_off,
                                flat_len,
                                slot_cap,
                                pmid,
                                pp,
                                ptag,
                                poff,
                                ov,
                                bloom_ptr,
                                sub_buckets,
                                sub_shift_hi,
                                sub_mask,
                                freq_limit,
                                forward_only,
                                bloom_sub_shift,
                                bloom_sub_mask,
                                non_unique_atomic,
                                sub_off,
                                scratch_k,
                                scratch_p,
                                sort_buf_k,
                            )
                        }
                    };
                }
                #[cfg(target_arch = "x86_64")]
                if use_bmi2 {
                    return unsafe { dispatch_kh!(process_hash_partition_bmi2) };
                }
                unsafe { dispatch_kh!(process_hash_partition_scalar) }
            },
        )
        .collect();

    // Immediately unmap transient `coarse_mid` (`3 B/k-mer`) and drop overflow staging before
    // any overflow append or subsequent `rev_ref` build!
    drop(coarse_mid);
    part_ov.clear();
    part_ov.shrink_to_fit();
    if total_kmers == 0 {
        non_unique_bm.clear();
    }

    // Safety: every partition wrote all `sub_buckets + 1` elements of its `poff` slice and all
    // `slot_cap` slots of both its `ptag` and `pp` slices.
    unsafe {
        offsets.set_len(off_len);
        tags.set_len(alloc_total);
        pos_strand.set_len(alloc_total);
    }

    let mut kept_entries = 0usize;
    for (kept, append_opt) in part_results {
        kept_entries += kept;
        if let Some((lp, tmp_tags, tmp_ps)) = append_opt {
            let new_base = tags.len();
            if new_base.saturating_add(tmp_tags.len()) > (u32::MAX as usize) - 8 {
                return Err(GeckoError::InvalidParameter(
                    "reference slab overflow exceeds 2^32 Hash-CSR entries; use smaller --memory budget"
                        .into(),
                ));
            }
            tags.extend_from_slice(&tmp_tags);
            pos_strand.extend_from_slice(&tmp_ps);
            let base_u32 = new_base as u32;
            let off_slice = &mut offsets[lp * (sub_buckets + 1)..(lp + 1) * (sub_buckets + 1)];
            for slot in off_slice {
                *slot += base_u32;
            }
        }
    }

    // Append the 8 trailing SIMD padding zeros past the final tag entry (zero reallocation on the
    // common non-append path since `tags` was allocated with capacity `alloc_total + 8`).
    tags.resize(tags.len() + 8, 0);

    let bloom = bloom_params.map(|(sb, kh)| BlockedBloom {
        blocks: bloom_blocks,
        sub_bits: sb,
        block_shift: 56 - sb,
        k_h: kh,
    });

    Ok(RefIndex {
        word_length: wl,
        freq_limit,
        soft_mask,
        spaced_mask,
        ref_data_len: seq.len,
        total_kmers,
        kept_entries,
        part_start: part_start as u16,
        part_end: part_end as u16,
        kind: RefIndexKind::HashCsr {
            sub_bits,
            hash_shift,
            offsets,
            tags,
            pos_strand: P::into_vec(pos_strand),
            bloom,
        },
        unique_bm: non_unique_bm,
        unique_bm_xor: u64::MAX,
        forward_only,
    })
}

/// Fast in-place insertion sort of paired `(keys, ps)` slices by `(key, ps)` for `n <= 24`.
#[inline(always)]
fn insertion_sort_pair<P: Copy + Ord>(keys: &mut [u64], ps: &mut [P]) {
    let n = keys.len();
    debug_assert_eq!(n, ps.len());
    let kptr = keys.as_mut_ptr();
    let pptr = ps.as_mut_ptr();
    for i in 1..n {
        unsafe {
            let k = *kptr.add(i);
            let p = *pptr.add(i);
            let mut j = i;
            while j > 0 {
                let prev_k = *kptr.add(j - 1);
                let prev_p = *pptr.add(j - 1);
                if prev_k > k || (prev_k == k && prev_p > p) {
                    *kptr.add(j) = prev_k;
                    *pptr.add(j) = prev_p;
                    j -= 1;
                } else {
                    break;
                }
            }
            *kptr.add(j) = k;
            *pptr.add(j) = p;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_kmer(s: &[u8]) -> (u64, u64) {
        let wl = s.len();
        let mask = kmer_mask(wl);
        let rc_shift = (2 * (wl - 1)) as u32;
        let mut kf = 0u64;
        let mut kr = 0u64;
        for &b in s {
            let v: u64 = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("invalid base"),
            };
            kf = ((kf << 2) | v) & mask;
            kr = (kr >> 2) | ((v ^ 3) << rc_shift);
        }
        (kf, kr)
    }

    #[test]
    fn direct_and_hash_csr_handle_strands_palindromes_and_freq_limit() {
        // Repeat a forward-only 16-mer 4 times, its revcomp 1 time, plus a palindrome 2 times.
        let fwd_k = "AAACCCGGGTTTACGT";
        let rev_k = "ACGTAAACCCGGGTTT"; // revcomp of fwd_k
        let pal_k = "ACGTACGTACGTACGT"; // palindrome

        let seq = format!(
            "{fwd_k}NNNN{fwd_k}NNNN{fwd_k}NNNN{fwd_k}NNNN{rev_k}NNNN{pal_k}NNNN{pal_k}NNNNacgtacgtacgtacgt"
        );
        let g = Genome::from_records("r", [("r0", seq.as_bytes())]);

        // Hash-CSR with freq_limit = 2:
        // fwd_k occurs 4 times (strand 0 -> dropped!), rev_k occurs 1 time (strand 1 -> kept!),
        // pal_k occurs 2 times (strand 0 -> kept!), lowercase tail skipped by soft_mask.
        let idx = RefIndex::build(&g, 16, 2, true).unwrap();
        assert!(!idx.is_direct());

        let (kf, kr) = encode_kmer(fwd_k.as_bytes());
        let canon = kf.min(kr);
        let hits = idx.lookup(canon);
        assert_eq!(hits.len(), 1, "only the single rev_k strand should remain");
        let expected_rev_pos = (4 * (16 + 4)) as u64;
        assert_eq!(hits[0].0, expected_rev_pos);
        assert_eq!(hits[0].1, if kr < kf { 0 } else { 1 });

        let (pkf, pkr) = encode_kmer(pal_k.as_bytes());
        assert_eq!(pkf, pkr);
        let pal_hits = idx.lookup(pkf);
        assert_eq!(pal_hits.len(), 2);
        assert_eq!(pal_hits[0].1, 0);
        assert_eq!(pal_hits[1].1, 0);

        // Direct CSR (wl = 8) with freq_limit = 2 on same sequence.
        let idx_dir = RefIndex::build(&g, 8, 2, true).unwrap();
        assert!(idx_dir.is_direct());
        let (dkf, dkr) = encode_kmer(b"ACGTACGT");
        assert_eq!(dkf, dkr);
        // "ACGTACGT" appears > 2 times inside pal_k (multiple offsets per copy), so filtered out!
        assert!(idx_dir.lookup(dkf).is_empty());

        // With higher freq_limit, it is retained.
        let idx_dir_hi = RefIndex::build(&g, 8, 100, true).unwrap();
        assert!(!idx_dir_hi.lookup(dkf).is_empty());
    }

    #[test]
    fn hash_csr_tag_probe_matches_linear_and_binary_search() {
        // Test bucket probing across sizes 1..=40 including branchless <=8 and looped <=32 AVX2 paths.
        let mut tags = Vec::new();
        for k in 0u32..14 {
            for _ in 0..3 {
                tags.push(k * 100 + 7);
            }
        }
        let n_valid = tags.len();
        tags.extend_from_slice(&[0u32; 8]); // SIMD padding

        #[cfg(target_arch = "x86_64")]
        let has_avx2 = std::is_x86_feature_detected!("avx2");
        #[cfg(not(target_arch = "x86_64"))]
        let has_avx2 = false;

        for end in 1..=n_valid {
            for target in [7u32, 107, 407, 907, 1307, 50, 9999] {
                let expected = probe_bucket_linear_u32(&tags, 0, end, target);
                let actual = probe_sorted_tags(&tags, 0, end, target, has_avx2);
                assert_eq!(actual, expected, "mismatch for end={end}, target={target}");
            }
        }
    }

    #[test]
    fn inv_fmix64_exact_roundtrip() {
        let edges = [
            0u64,
            1,
            2,
            3,
            u64::MAX,
            u64::MAX - 1,
            1u64 << 63,
            0x5555_5555_5555_5555,
            0xAAAA_AAAA_AAAA_AAAA,
            0x0123_4567_89AB_CDEF,
        ];
        for &x in &edges {
            assert_eq!(inv_fmix64(fmix64(x)), x);
            assert_eq!(fmix64(inv_fmix64(x)), x);
        }
        let mut rng = 0x1337_C0DE_DEAD_BEEFu64;
        for _ in 0..10_000 {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            assert_eq!(inv_fmix64(fmix64(rng)), rng);
        }
    }

    #[test]
    fn is_palindrome_key_contiguous_and_spaced() {
        use crate::spaced::SpacedSeed;

        let mut rng = 0x4242_4242_1234_5678u64;
        let alphabet = *b"ACGT";
        let mut rand_seq = |wl: usize| -> Vec<u8> {
            (0..wl)
                .map(|_| {
                    rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                    alphabet[((rng >> 32) & 3) as usize]
                })
                .collect()
        };

        // Contiguous wl = 4, 12, 15, 16, 32
        for wl in [4usize, 12, 15, 16, 32] {
            let mask = kmer_mask(wl);
            assert_eq!(can_key_be_palindrome(wl, None), (wl & 1) == 0);
            for _ in 0..200 {
                let s = rand_seq(wl);
                let (kf, kr) = encode_kmer(&s);
                let canon = kf.min(kr);
                assert_eq!(is_palindrome_key(canon, wl, mask), kf == kr);
            }
            if (wl & 1) == 0 {
                // Construct an explicit palindrome.
                let half = rand_seq(wl / 2);
                let mut pal = half.clone();
                for &b in half.iter().rev() {
                    pal.push(match b {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        _ => b'A',
                    });
                }
                let (kf, kr) = encode_kmer(&pal);
                assert_eq!(kf, kr);
                assert!(is_palindrome_key(kf, wl, mask));
            }
        }

        // Spaced seeds: even (`w14s20`), odd with center `'0'` (`11011`, `1100011`), odd with center `'1'` (`11111`)
        for (pat, expect_can_pal) in [
            ("w14s20", true),
            ("11011", true),
            ("1100011", true),
            ("1101011", false),
        ] {
            let sp = SpacedSeed::parse(pat).unwrap();
            let sm = sp.gapped_mask_2bit();
            assert_eq!(
                can_key_be_palindrome(sp.span, sm),
                expect_can_pal,
                "pat={pat}"
            );
            let eff_mask = sm.unwrap_or_else(|| kmer_mask(sp.span));
            for _ in 0..200 {
                let s = rand_seq(sp.span);
                let (kf, kr) = encode_kmer(&s);
                let sf = kf & eff_mask;
                let sr = kr & eff_mask;
                assert_eq!(
                    is_palindrome_key(sf.min(sr), sp.span, eff_mask),
                    sf == sr,
                    "pat={pat}"
                );
            }
        }
    }

    #[test]
    fn sorting_networks_01_principle_and_pass2_parity() {
        // Exhaustively verify all 2^n binary inputs (Knuth 0-1 principle) for n = 2..=6.
        for mask in 0u32..(1 << 5) {
            let mut k = [0u64; 5];
            let mut p = [0u32; 5];
            for i in 0..5 {
                k[i] = ((mask >> i) & 1) as u64;
                p[i] = (4 - i) as u32;
            }
            let mut expected: Vec<(u64, u32)> = k.iter().copied().zip(p.iter().copied()).collect();
            expected.sort();
            let (mut k0, mut k1, mut k2, mut k3, mut k4) = (k[0], k[1], k[2], k[3], k[4]);
            let (mut p0, mut p1, mut p2, mut p3, mut p4) = (p[0], p[1], p[2], p[3], p[4]);
            cas_pair(&mut k0, &mut p0, &mut k1, &mut p1);
            cas_pair(&mut k3, &mut p3, &mut k4, &mut p4);
            cas_pair(&mut k2, &mut p2, &mut k4, &mut p4);
            cas_pair(&mut k2, &mut p2, &mut k3, &mut p3);
            cas_pair(&mut k0, &mut p0, &mut k3, &mut p3);
            cas_pair(&mut k0, &mut p0, &mut k2, &mut p2);
            cas_pair(&mut k1, &mut p1, &mut k4, &mut p4);
            cas_pair(&mut k1, &mut p1, &mut k3, &mut p3);
            cas_pair(&mut k1, &mut p1, &mut k2, &mut p2);
            let actual = [(k0, p0), (k1, p1), (k2, p2), (k3, p3), (k4, p4)];
            assert_eq!(actual.as_slice(), expected.as_slice());
        }
        for mask in 0u32..(1 << 6) {
            let mut k = [0u64; 6];
            let mut p = [0u32; 6];
            for i in 0..6 {
                k[i] = ((mask >> i) & 1) as u64;
                p[i] = (5 - i) as u32;
            }
            let mut expected: Vec<(u64, u32)> = k.iter().copied().zip(p.iter().copied()).collect();
            expected.sort();
            let (mut k0, mut k1, mut k2, mut k3, mut k4, mut k5) =
                (k[0], k[1], k[2], k[3], k[4], k[5]);
            let (mut p0, mut p1, mut p2, mut p3, mut p4, mut p5) =
                (p[0], p[1], p[2], p[3], p[4], p[5]);
            cas_pair(&mut k1, &mut p1, &mut k2, &mut p2);
            cas_pair(&mut k4, &mut p4, &mut k5, &mut p5);
            cas_pair(&mut k0, &mut p0, &mut k2, &mut p2);
            cas_pair(&mut k3, &mut p3, &mut k5, &mut p5);
            cas_pair(&mut k0, &mut p0, &mut k1, &mut p1);
            cas_pair(&mut k3, &mut p3, &mut k4, &mut p4);
            cas_pair(&mut k2, &mut p2, &mut k5, &mut p5);
            cas_pair(&mut k0, &mut p0, &mut k3, &mut p3);
            cas_pair(&mut k1, &mut p1, &mut k4, &mut p4);
            cas_pair(&mut k2, &mut p2, &mut k4, &mut p4);
            cas_pair(&mut k1, &mut p1, &mut k3, &mut p3);
            cas_pair(&mut k2, &mut p2, &mut k3, &mut p3);
            let actual = [(k0, p0), (k1, p1), (k2, p2), (k3, p3), (k4, p4), (k5, p5)];
            assert_eq!(actual.as_slice(), expected.as_slice());
        }

        // Verify full Pass 2 HashCsr bucket ordering parity across buckets of size 1..=30.
        let alphabet = *b"ACGT";
        let mut rng = 0x9988_7766_5544_3322u64;
        let seq: Vec<u8> = (0..100_000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();
        let g = Genome::from_records("r", [("r0", seq.as_slice())]);
        let idx = RefIndex::build(&g, 16, 100, true).unwrap();
        let RefIndexKind::HashCsr {
            sub_bits,
            offsets,
            tags,
            pos_strand,
            ..
        } = &idx.kind
        else {
            panic!("expected HashCsr");
        };
        let sub_buckets = 1usize << *sub_bits;
        let mut max_cnt = 0usize;
        let mask16 = kmer_mask(16);
        let h_rot_at = |pos: usize| -> u64 {
            let mut kf = 0u64;
            let mut kr = 0u64;
            for k in 0..16usize {
                let b = g.seq.base_2bit(pos + k) as u64;
                kf = ((kf << 2) | b) & mask16;
                kr = (kr >> 2) | ((b ^ 3) << 30);
            }
            fmix64(kf.min(kr)).rotate_left(32)
        };
        for part_off in offsets.chunks_exact(sub_buckets + 1) {
            for w in part_off.windows(2) {
                let s = w[0] as usize;
                let e = w[1] as usize;
                if e > s {
                    max_cnt = max_cnt.max(e - s);
                    for i in s + 1..e {
                        assert!(
                            tags[i - 1] <= tags[i],
                            "tags not sorted in bucket [{s}..{e})"
                        );
                        let (p_prev, s_prev) = pos_strand.get_unpacked(i - 1);
                        let (p_cur, s_cur) = pos_strand.get_unpacked(i);
                        let a = (h_rot_at(p_prev as usize), (p_prev << 1) | s_prev as u64);
                        let b = (h_rot_at(p_cur as usize), (p_cur << 1) | s_cur as u64);
                        assert!(a < b, "Pass 2 pair ordering mismatch in bucket [{s}..{e})");
                    }
                }
            }
        }
        assert!(
            max_cnt >= 6,
            "expected multi-entry buckets >= 6, got {max_cnt}"
        );
    }

    #[test]
    fn blocked_bloom_zero_false_negatives_and_fpr_single_and_multislab() {
        // Test all k_h = 1..=8 directly on BlockedBloom.
        for k_h in 1u8..=8 {
            let mut blocks = vec![BloomBlock([0u64; 8]); 64];
            let mut inserted = Vec::new();
            for i in 0u32..500 {
                let h = fmix64(i as u64 + 1);
                let blk = (h >> 58) as usize;
                let tag = h as u32;
                insert_in_bloom_block(&mut blocks[blk], tag, k_h as usize);
                inserted.push((blk, tag));
            }
            let bloom = BlockedBloom {
                blocks,
                sub_bits: 0,
                block_shift: 58,
                k_h,
            };
            for &(blk, tag) in &inserted {
                assert!(
                    unsafe { bloom.contains_unchecked(blk, tag) },
                    "FN at k_h={k_h}"
                );
            }
            let mut fps = 0usize;
            for i in 10_000u32..15_000 {
                let h = fmix64(i as u64);
                if unsafe { bloom.contains_unchecked((h >> 58) as usize, h as u32) } {
                    fps += 1;
                }
            }
            assert!(fps < 1_500, "unexpectedly high FPR {fps}/5000 at k_h={k_h}");
        }

        // Build full HashCsr (> 4,096 kmers per 1/4 slab) with a dropped high-frequency repeat
        // and verify 0 FN on single-slab AND 4-slab multi-slab.
        let alphabet = *b"ACGT";
        let mut rng = 0xAA55_AA55_1122_3344u64;
        let mut seq: Vec<u8> = (0..24_000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();
        let drop_k = b"AAAACCCCGGGGTTTT";
        for k in 0..6 {
            seq[200 + k * 40..200 + k * 40 + 16].copy_from_slice(drop_k);
        }
        let g = Genome::from_records("b", [("b0", seq.as_slice())]);
        let (dkf, dkr) = encode_kmer(drop_k);
        let drop_canon = dkf.min(dkr);

        for (ps, pe) in [(0u16, 256u16), (0, 64), (64, 128), (128, 192), (192, 256)] {
            let idx = RefIndex::build_partitions_with_mask(&g, ps, pe, 16, None, 3, true).unwrap();
            let RefIndexKind::HashCsr {
                bloom: Some(ref bloom),
                ..
            } = idx.kind
            else {
                panic!("expected BlockedBloom to be enabled for [{ps}..{pe})");
            };
            // Every valid k-mer in this slab with count <= 3 MUST be present in Bloom & lookup.
            let mut kept_checked = 0usize;
            for_each_canonical_kmer_2bit::<false, _>(
                &g.seq,
                0,
                g.seq.len - 15,
                16,
                true,
                |canon, _, _| {
                    let h = fmix64(canon);
                    let p = (h >> 56) as u16;
                    if p < ps || p >= pe {
                        return;
                    }
                    if canon == drop_canon {
                        assert!(idx.lookup(canon).is_empty());
                        return;
                    }
                    let blk =
                        ((h >> bloom.block_shift) as usize) - ((ps as usize) << bloom.sub_bits);
                    let hits = idx.lookup(canon);
                    if !hits.is_empty() {
                        assert!(unsafe { bloom.contains_unchecked(blk, h as u32) });
                        kept_checked += 1;
                    }
                },
            );
            assert!(kept_checked > 2_000);
        }
    }

    #[test]
    fn unique_bm_exact_vs_naive_counts_across_all_modes() {
        use crate::spaced::SpacedSeed;
        use std::collections::HashMap;

        let alphabet = *b"ACGT";
        let mut rng = 0x55AA_1234_9876_FEDCu64;
        let mut seq: Vec<u8> = (0..10_000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();
        // Inject duplicates, reverse-complement duplicate, palindrome, N, and lowercase.
        seq[100..116].copy_from_slice(b"ACGTACGTACGTACGT"); // palindrome
        seq[300..316].copy_from_slice(b"AAACCCGGGTTTACGT");
        seq[500..516].copy_from_slice(b"ACGTAAACCCGGGTTT"); // revcomp of above
        seq[700..716].copy_from_slice(b"TTTTGGGGCCCCAAAA");
        seq[900..916].copy_from_slice(b"TTTTGGGGCCCCAAAA"); // fwd dup
        seq[1_100..1_110].fill(b'N');
        seq[1_300..1_320].copy_from_slice(b"acgtacgtacgtacgtacgt");

        let g = Genome::from_records("u", [("u0", seq.as_slice())]);
        let spaced = SpacedSeed::parse("w14s20").unwrap();

        let configs: Vec<(usize, Option<u64>, u16, u16, u64, bool)> = vec![
            (8, None, 0, 256, 2, false),   // DirectCsr canonical
            (8, None, 0, 256, 2, true),    // DirectCsr forward-only
            (16, None, 0, 256, 2, false),  // HashCsr canonical
            (16, None, 64, 192, 2, false), // HashCsr multi-slab canonical
            (16, None, 0, 256, 2, true),   // HashCsr forward-only
            (spaced.span, spaced.gapped_mask_2bit(), 0, 256, 2, false), // Spaced canonical
        ];

        for (wl, sm, ps, pe, fl, fwd_only) in configs {
            let idx = RefIndex::build_partitions_with_mask_and_strand(
                &g, ps, pe, wl, sm, fl, true, fwd_only,
            )
            .unwrap();
            let eff_mask = sm.unwrap_or_else(|| kmer_mask(wl));

            // Compute exact naive counts of every valid k-mer in this slab.
            let mut key_counts: HashMap<u64, usize> = HashMap::new();
            let mut pos_key: HashMap<usize, u64> = HashMap::new();
            let mut record = |key: u64, pos: u64, _strand: u64| {
                let p = if idx.is_direct() {
                    let sub_bits = (2 * wl - 8) as u32;
                    (key >> sub_bits) as u16
                } else {
                    (fmix64(key) >> 56) as u16
                };
                if p >= ps && p < pe {
                    *key_counts.entry(key).or_insert(0) += 1;
                    pos_key.insert(pos as usize, key);
                }
            };
            if fwd_only {
                for_each_canonical_kmer_2bit_opt::<true, _>(
                    &g.seq,
                    0,
                    g.seq.len - wl + 1,
                    wl,
                    sm,
                    true,
                    &mut record,
                );
            } else {
                for_each_canonical_kmer_2bit_opt::<false, _>(
                    &g.seq,
                    0,
                    g.seq.len - wl + 1,
                    wl,
                    sm,
                    true,
                    &mut record,
                );
            }

            for pos in 0..g.seq.len {
                let actual =
                    (((idx.unique_bm[pos >> 6] ^ idx.unique_bm_xor) >> (pos & 63)) & 1) != 0;
                let expected = match pos_key.get(&pos) {
                    Some(&k) => {
                        let cnt = key_counts[&k];
                        let pal = !fwd_only && is_palindrome_key(k, wl, eff_mask);
                        cnt == 1 && !pal
                    }
                    None => false,
                };
                assert_eq!(
                    actual,
                    expected,
                    "unique_bm mismatch at pos={pos}, wl={wl}, spaced={}, slab=[{ps}..{pe}), fwd_only={fwd_only}",
                    sm.is_some()
                );
            }
        }
    }

    #[test]
    fn freq_limit_zero_exact_direct_and_hash_canon_and_fwd() {
        let g = Genome::from_records(
            "z",
            [(
                "z0",
                b"ACGTACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGT".as_slice(),
            )],
        );
        for wl in [8usize, 16] {
            for fwd_only in [false, true] {
                let idx =
                    RefIndex::build_with_mask_and_strand(&g, wl, None, 0, true, fwd_only).unwrap();
                assert_eq!(idx.kept_entries, 0);
                assert!(idx.total_kmers > 0);
                assert!(idx.unique_bm.is_empty());
                if let RefIndexKind::HashCsr { ref bloom, .. } = idx.kind {
                    assert!(bloom.is_none());
                }
                let (kf, kr) = encode_kmer(&g.seq.to_ascii_vec()[0..wl]);
                let key = if fwd_only { kf } else { kf.min(kr) };
                assert!(idx.lookup(key).is_empty());
            }
        }
    }

    #[test]
    fn choose_bloom_params_scales_proportionally_for_multislab() {
        // 16.0 MiB (sb=10, k_h=6) on gorilla_chr21 (15.70M kmers, 256 parts)
        assert_eq!(choose_bloom_params(15_703_781, 256, 20), Some((10, 6)));
        // 32.0 MiB (sb=11, k_h=4) on mouse_chr16 single-slab (50.62M kmers, 256 parts)
        assert_eq!(choose_bloom_params(50_617_016, 256, 20), Some((11, 4)));
        // 10.75 MiB (sb=11, k_h=4) on mouse_chr16 3-slab (17.00M kmers, 86 parts)
        assert_eq!(choose_bloom_params(17_004_154, 86, 20), Some((11, 4)));
        // 6.50 MiB (sb=11, k_h=4) on mouse_chr16 5-slab (10.28M kmers, 52 parts)
        assert_eq!(choose_bloom_params(10_281_581, 52, 20), Some((11, 4)));
        // Disabled for tiny refs or freq_limit == 0
        assert_eq!(choose_bloom_params(1_000, 256, 20), None);
        assert_eq!(choose_bloom_params(50_617_016, 256, 0), None);
    }

    #[test]
    fn sort_small_2_to_6_exhaustive_permutations_and_ties() {
        fn next_perm(a: &mut [usize]) -> bool {
            let n = a.len();
            if n < 2 {
                return false;
            }
            let mut i = n - 2;
            while a[i] >= a[i + 1] {
                if i == 0 {
                    return false;
                }
                i -= 1;
            }
            let mut j = n - 1;
            while a[j] <= a[i] {
                j -= 1;
            }
            a.swap(i, j);
            a[i + 1..].reverse();
            true
        }

        // 1. Exhaustive n! permutations for n = 2..=6 (2+6+24+120+720 = 872 permutations)
        //    with both distinct keys and duplicate keys + reverse tie-break `ps`.
        let mut perm_count = 0usize;
        for n in 2..=6 {
            let mut perm: Vec<usize> = (0..n).collect();
            loop {
                perm_count += 1;
                // Distinct keys (including extreme u64 boundaries)
                let base_keys = [
                    0u64,
                    1,
                    0x7FFF_FFFF_FFFF_FFFF,
                    0x8000_0000_0000_0000,
                    u64::MAX - 1,
                    u64::MAX,
                ];
                let mut sk: Vec<u64> = perm.iter().map(|&idx| base_keys[idx]).collect();
                let mut sp: Vec<u32> = perm.iter().map(|&idx| ((n - idx) as u32) << 1).collect();
                let mut expected: Vec<(u64, u32)> =
                    sk.iter().copied().zip(sp.iter().copied()).collect();
                expected.sort();
                unsafe { sort_small_2_to_6(sk.as_mut_ptr(), sp.as_mut_ptr(), n) };
                let got: Vec<(u64, u32)> = sk.into_iter().zip(sp).collect();
                assert_eq!(got, expected, "distinct perm failed for n={n}");

                // Duplicate keys with distinct `ps` tie-breakers (tests 128-bit comparator)
                let dup_keys = [10u64, 10, 20, 20, 30, 30];
                let mut sk_d: Vec<u64> = perm.iter().map(|&idx| dup_keys[idx]).collect();
                let mut sp_d: Vec<u32> = perm.iter().map(|&idx| (idx as u32) * 3 + 1).collect();
                let mut expected_d: Vec<(u64, u32)> =
                    sk_d.iter().copied().zip(sp_d.iter().copied()).collect();
                expected_d.sort();
                unsafe { sort_small_2_to_6(sk_d.as_mut_ptr(), sp_d.as_mut_ptr(), n) };
                let got_d: Vec<(u64, u32)> = sk_d.into_iter().zip(sp_d).collect();
                assert_eq!(got_d, expected_d, "dup-tie perm failed for n={n}");

                if !next_perm(&mut perm) {
                    break;
                }
            }

            // 2. All 2^n Knuth 0-1 binary patterns with reverse tie-break `ps`.
            for mask in 0..(1usize << n) {
                let mut sk: Vec<u64> = (0..n).map(|i| ((mask >> i) & 1) as u64).collect();
                let mut sp: Vec<u32> = (0..n).map(|i| ((n - i) as u32) << 1).collect();
                let mut expected: Vec<(u64, u32)> =
                    sk.iter().copied().zip(sp.iter().copied()).collect();
                expected.sort();
                unsafe { sort_small_2_to_6(sk.as_mut_ptr(), sp.as_mut_ptr(), n) };
                let got: Vec<(u64, u32)> = sk.into_iter().zip(sp).collect();
                assert_eq!(got, expected, "0-1 pattern mask={mask:#b} failed for n={n}");
            }
        }
        assert_eq!(perm_count, 2 + 6 + 24 + 120 + 720);
    }

    #[test]
    fn blocked_bloom_contains_const_parity_all_kh() {
        for k_h in 1u8..=8 {
            let mut blocks = vec![BloomBlock([0u64; 8]); 64];
            let mut inserted = Vec::new();
            for i in 0u64..400 {
                let h = fmix64(i + 100);
                let blk = (h >> 58) as usize;
                let tag = h as u32;
                insert_in_bloom_block(&mut blocks[blk], tag, k_h as usize);
                inserted.push(h);
            }
            let bloom = BlockedBloom {
                blocks,
                sub_bits: 0,
                block_shift: 58,
                k_h,
            };
            let check_h = |h: u64| -> (bool, bool, bool) {
                let blk = (h >> 58) as usize;
                let tag = h as u32;
                let dyn_res = unsafe { bloom.contains_unchecked(blk, tag) };
                let (const_u, const_s) = match k_h {
                    1 => (
                        unsafe { bloom.contains_unchecked_const::<1>(blk, tag) },
                        bloom.contains_const::<1>(h),
                    ),
                    2 => (
                        unsafe { bloom.contains_unchecked_const::<2>(blk, tag) },
                        bloom.contains_const::<2>(h),
                    ),
                    3 => (
                        unsafe { bloom.contains_unchecked_const::<3>(blk, tag) },
                        bloom.contains_const::<3>(h),
                    ),
                    4 => (
                        unsafe { bloom.contains_unchecked_const::<4>(blk, tag) },
                        bloom.contains_const::<4>(h),
                    ),
                    5 => (
                        unsafe { bloom.contains_unchecked_const::<5>(blk, tag) },
                        bloom.contains_const::<5>(h),
                    ),
                    6 => (
                        unsafe { bloom.contains_unchecked_const::<6>(blk, tag) },
                        bloom.contains_const::<6>(h),
                    ),
                    7 => (
                        unsafe { bloom.contains_unchecked_const::<7>(blk, tag) },
                        bloom.contains_const::<7>(h),
                    ),
                    _ => (
                        unsafe { bloom.contains_unchecked_const::<8>(blk, tag) },
                        bloom.contains_const::<8>(h),
                    ),
                };
                (dyn_res, const_u, const_s)
            };
            for &h in &inserted {
                let (d, cu, cs) = check_h(h);
                assert!(d && cu && cs, "hit mismatch at k_h={k_h}");
            }
            for i in 10_000u64..14_000 {
                let h = fmix64(i);
                let (d, cu, cs) = check_h(h);
                assert_eq!(d, cu, "unchecked_const mismatch at k_h={k_h}");
                assert_eq!(d, cs, "contains_const mismatch at k_h={k_h}");
            }
        }
    }

    #[test]
    fn hash_csr_adaptive_pbits_8_9_10_single_and_multislab_parity() {
        use crate::pipeline::{CompareConfig, SeedingMode, compare_with_ref_index};
        use crate::spaced::SpacedSeed;

        let alphabet = *b"ACGT";
        let mut rng = 0x3141_5926_5358_9793u64;
        let mut seq_r: Vec<u8> = (0..32_000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();
        // Inject a repeated 16-mer (exceeding freq_limit=3) and a couple of matches for query.
        for k in 0..5 {
            seq_r[500 + k * 50..500 + k * 50 + 16].copy_from_slice(b"ACGTTGCAAACCTTGG");
        }
        let mut seq_q = seq_r.clone();
        // Mutate every 9th base outside two 120-bp conserved HSPs.
        for (i, b) in seq_q.iter_mut().enumerate() {
            if !(2_000..2_120).contains(&i) && !(18_000..18_120).contains(&i) && i % 9 == 0 {
                *b = match *b {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
        }
        let gr = Genome::from_records("r", [("r0", seq_r.as_slice())]);
        let gq = Genome::from_records("q", [("q0", seq_q.as_slice())]);
        let spaced = SpacedSeed::parse("w14s20").unwrap();

        let modes: Vec<(usize, Option<u64>, bool)> = vec![
            (16, None, false),
            (16, None, true),
            (spaced.span, spaced.gapped_mask_2bit(), false),
        ];

        for (wl, sm, fwd_only) in modes {
            let base = build_hash_csr_inner::<u32>(
                &gr,
                0,
                gr.seq.len,
                0,
                256,
                wl,
                sm,
                3,
                true,
                fwd_only,
                true,
                Some(8),
                None,
            )
            .unwrap();
            let cfg = CompareConfig {
                word_length: wl,
                min_length: 80,
                min_similarity: 70.0,
                freq_limit: 3,
                compare_reverse: !fwd_only,
                seeding_mode: SeedingMode::Exact,
                spaced_seed: if sm.is_some() {
                    Some(spaced.clone())
                } else {
                    None
                },
                ..CompareConfig::default()
            };
            let rev_r = if fwd_only {
                None
            } else {
                Some(gr.reverse_complement())
            };
            let base_cmp = compare_with_ref_index(&gq, &gr, &base, rev_r.as_ref(), &cfg).unwrap();

            for p_bits in [8u32, 9, 10] {
                let single = build_hash_csr_inner::<u32>(
                    &gr,
                    0,
                    gr.seq.len,
                    0,
                    256,
                    wl,
                    sm,
                    3,
                    true,
                    fwd_only,
                    true,
                    Some(p_bits),
                    None,
                )
                .unwrap();
                assert_eq!(single.total_kmers, base.total_kmers);
                assert_eq!(single.kept_entries, base.kept_entries);
                assert_eq!(single.unique_bm, base.unique_bm);

                // Verify lookup parity on 500 sampled k-mers.
                let mut checked = 0usize;
                for_each_canonical_kmer_2bit_opt::<false, _>(
                    &gr.seq,
                    0,
                    2_000,
                    wl,
                    sm,
                    true,
                    |canon, _, _| {
                        if !fwd_only && checked < 500 {
                            assert_eq!(single.lookup(canon), base.lookup(canon));
                            checked += 1;
                        }
                    },
                );

                // Verify streaming compare parity against p_bits=8 baseline.
                let cmp = compare_with_ref_index(&gq, &gr, &single, rev_r.as_ref(), &cfg).unwrap();
                assert_eq!(
                    cmp.fragments, base_cmp.fragments,
                    "single-slab fragment parity failed for p_bits={p_bits}, wl={wl}, fwd={fwd_only}"
                );

                // Verify 4-slab [0..64, 64..128, 128..192, 192..256] union parity.
                let mut slab_total = 0usize;
                let mut slab_kept = 0usize;
                let mut slabs = Vec::new();
                for (ps, pe) in [(0usize, 64usize), (64, 128), (128, 192), (192, 256)] {
                    let sidx = build_hash_csr_inner::<u32>(
                        &gr,
                        0,
                        gr.seq.len,
                        ps,
                        pe,
                        wl,
                        sm,
                        3,
                        true,
                        fwd_only,
                        true,
                        Some(p_bits),
                        None,
                    )
                    .unwrap();
                    slab_total += sidx.total_kmers;
                    slab_kept += sidx.kept_entries;
                    slabs.push(sidx);
                }
                assert_eq!(slab_total, base.total_kmers);
                assert_eq!(slab_kept, base.kept_entries);
                if !fwd_only {
                    let mut checked_s = 0usize;
                    for_each_canonical_kmer_2bit_opt::<false, _>(
                        &gr.seq,
                        0,
                        2_000,
                        wl,
                        sm,
                        true,
                        |canon, _, _| {
                            if checked_s < 400 {
                                let p = (fmix64(canon) >> 56) as usize;
                                let slab = &slabs[p / 64];
                                assert_eq!(slab.lookup(canon), base.lookup(canon));
                                checked_s += 1;
                            }
                        },
                    );
                }
            }
        }
    }

    #[test]
    fn hash_csr_single_pass_overflow_filtered_and_kept() {
        // Adversarial homopolymer + random tail: forces massive single-partition overflow into
        // both the filtered-in-place (`write_idx <= slot_cap`) and kept-append (`write_idx > slot_cap`)
        // branches of `process_hash_partition_overflow`.
        let mut seq = vec![b'A'; 6_000];
        let alphabet = *b"ACGT";
        let mut rng = 0xDEAD_1234_5678_BEEFu64;
        for _ in 0..4_000 {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            seq.push(alphabet[((rng >> 32) & 3) as usize]);
        }
        let g = Genome::from_records("ov", [("ov0", seq.as_slice())]);
        let (kf_a, kr_a) = encode_kmer(&[b'A'; 16]);
        let canon_a = kf_a.min(kr_a);

        // Case 1: freq_limit = 2 -> poly-A dropped during Pass 2 (`write_idx <= slot_cap`, in-place copy).
        let idx_filt_base = build_hash_csr_inner::<u32>(
            &g,
            0,
            g.seq.len,
            0,
            256,
            16,
            None,
            2,
            true,
            false,
            true,
            Some(8),
            Some(8_192),
        )
        .unwrap();
        let idx_filt_ov = build_hash_csr_inner::<u32>(
            &g,
            0,
            g.seq.len,
            0,
            256,
            16,
            None,
            2,
            true,
            false,
            true,
            Some(8),
            Some(16),
        )
        .unwrap();
        assert_eq!(idx_filt_ov.total_kmers, idx_filt_base.total_kmers);
        assert_eq!(idx_filt_ov.kept_entries, idx_filt_base.kept_entries);
        assert_eq!(idx_filt_ov.unique_bm, idx_filt_base.unique_bm);
        assert!(idx_filt_ov.lookup(canon_a).is_empty());

        // Case 2: freq_limit = 50_000 -> poly-A kept (`write_idx > slot_cap`, append + offset patch).
        let idx_kept_base = build_hash_csr_inner::<u32>(
            &g,
            0,
            g.seq.len,
            0,
            256,
            16,
            None,
            50_000,
            true,
            false,
            true,
            Some(8),
            Some(8_192),
        )
        .unwrap();
        let idx_kept_ov = build_hash_csr_inner::<u32>(
            &g,
            0,
            g.seq.len,
            0,
            256,
            16,
            None,
            50_000,
            true,
            false,
            true,
            Some(8),
            Some(16),
        )
        .unwrap();
        assert_eq!(idx_kept_ov.total_kmers, idx_kept_base.total_kmers);
        assert_eq!(idx_kept_ov.kept_entries, idx_kept_base.kept_entries);
        assert_eq!(idx_kept_ov.unique_bm, idx_kept_base.unique_bm);
        assert_eq!(idx_kept_ov.lookup(canon_a), idx_kept_base.lookup(canon_a));
        assert!(idx_kept_ov.lookup(canon_a).len() >= 5_980);

        // Spot-check random tail lookups across both overflow paths.
        for pos in (6_050..9_500).step_by(17) {
            let (kf, kr) = encode_kmer(&seq[pos..pos + 16]);
            let c = kf.min(kr);
            assert_eq!(idx_filt_ov.lookup(c), idx_filt_base.lookup(c));
            assert_eq!(idx_kept_ov.lookup(c), idx_kept_base.lookup(c));
        }
    }

    #[test]
    fn estimate_bytes_matches_split_tag_hash_csr_peak() {
        let dummy = Genome::from_records("e", [("e0", vec![b'A'; 1024])]);
        for &(kmers, min_mib, max_mib) in &[
            (15_700_000usize, 180usize, 240usize),
            (50_600_000usize, 600usize, 730usize),
            (90_100_000usize, 1_050usize, 1_240usize),
            (205_000_000usize, 2_350usize, 2_700usize),
        ] {
            let est = RefIndex::estimate_bytes_for_kmers(&dummy, kmers, 32, None, false);
            let mib = est / (1024 * 1024);
            assert!(
                (min_mib..=max_mib).contains(&mib),
                "estimate_bytes_for_kmers({kmers}) = {mib} MiB outside [{min_mib}..={max_mib}]"
            );
        }
    }
}
