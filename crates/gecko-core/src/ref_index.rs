use rayon::prelude::*;

use crate::error::{GeckoError, Result};
use crate::seq::Genome;

const NUM_PARTITIONS: usize = 256;
const STAGE_CAP: usize = 64;

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
pub(crate) fn base_lut(soft_mask: bool) -> &'static [u8; 256] {
    if soft_mask {
        &BASE_LUT_SOFT_MASK
    } else {
        &BASE_LUT_NO_MASK
    }
}

/// Invertible MurmurHash3 64-bit finalizer (`fmix64`) used for uniform bucket assignment.
#[inline(always)]
pub(crate) fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^= k >> 33;
    k
}

/// Software prefetch into L1 cache (`_MM_HINT_T0` on x86_64, no-op elsewhere).
#[inline(always)]
pub(crate) fn prefetch_t0<T>(ptr: *const T) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        core::arch::x86_64::_mm_prefetch(ptr as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = ptr;
    }
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

    #[inline(always)]
    pub(crate) fn prefetch_at(&self, index: usize) {
        match self {
            Self::U32(v) => {
                if index < v.len() {
                    prefetch_t0(unsafe { v.as_ptr().add(index) });
                }
            }
            Self::U64(v) => {
                if index < v.len() {
                    prefetch_t0(unsafe { v.as_ptr().add(index) });
                }
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
    /// `offsets` has `256 * (sub_keys + 1)` entries so each of the 256 partitions has its own
    /// sentinel slot (`b_idx = key + (key >> sub_bits)`).
    DirectCsr {
        sub_bits: u32,
        offsets: Vec<u32>,
        pos_strand: PosStrandVec,
    },
    /// 2-level parallel radix-partitioned Hash-CSR keyed by `fmix64(canonical_key)`.
    /// `offsets` has `256 * (2^(P - 8) + 1)` entries (`b_idx = b + (b >> (P - 8))`).
    HashCsr {
        sub_bits: u32,
        hash_shift: u32,
        offsets: Vec<u32>,
        keys: Vec<u64>,
        pos_strand: PosStrandVec,
    },
}

/// Single-reference canonical k-mer index supporting lock-free streaming queries.
///
/// Indexes only Reference Y once at `word_length` (`4..=32`) using canonical keys
/// `min(kf, kr)` and a 1-bit strand flag `(kr < kf)`, with per-strand `freq_limit` filtering.
#[derive(Debug, Clone)]
pub struct RefIndex {
    pub word_length: usize,
    pub freq_limit: u64,
    pub soft_mask: bool,
    pub ref_data_len: usize,
    pub total_kmers: usize,
    pub kept_entries: usize,
    pub(crate) kind: RefIndexKind,
}

/// Whether Direct CSR mode should be used for `(word_length, ref_total_bases)`.
#[inline]
pub fn use_direct_csr(word_length: usize, ref_total_bases: u64) -> bool {
    word_length <= 10
        || (word_length <= 12
            && (1u64 << (2 * word_length)) <= 2u64.saturating_mul(ref_total_bases))
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

/// Count exact valid k-mers in `reference.data` in parallel (used for RAM estimation).
pub fn count_valid_kmers(reference: &Genome, word_length: usize, soft_mask: bool) -> usize {
    let data = &reference.data;
    if data.len() < word_length {
        return 0;
    }
    let lut = base_lut(soft_mask);
    let ranges = make_chunk_ranges(data.len(), word_length);
    ranges
        .par_iter()
        .map(|&(start, end)| {
            let scan_end = (end + word_length - 1).min(data.len());
            let mut valid_run = 0usize;
            let mut count = 0usize;
            for &b in &data[start..scan_end] {
                if lut[b as usize] != 0xFF {
                    valid_run += 1;
                    if valid_run >= word_length {
                        count += 1;
                    }
                } else {
                    valid_run = 0;
                }
            }
            count
        })
        .sum()
}

impl RefIndex {
    /// Estimate peak RAM in bytes for a known `total_kmers` count.
    fn estimate_bytes_for_kmers(
        reference: &Genome,
        total_kmers: usize,
        word_length: usize,
        compare_reverse: bool,
    ) -> usize {
        if total_kmers >= u32::MAX as usize {
            return usize::MAX;
        }
        let ps_bytes = if reference.data.len() < (1usize << 31) {
            std::mem::size_of::<u32>()
        } else {
            std::mem::size_of::<u64>()
        };
        let rev_bytes = if compare_reverse {
            reference.data.len()
        } else {
            0
        };
        let threads = rayon::current_num_threads().max(1);

        if use_direct_csr(word_length, reference.total_bases) {
            let sub_keys = 1usize << (2 * word_length - 8);
            let offsets_bytes = NUM_PARTITIONS * (sub_keys + 1) * std::mem::size_of::<u32>();
            let ps_total = total_kmers.saturating_mul(ps_bytes);
            // Transient u32 coarse keys and per-thread scratch are freed before `rev_ref` is allocated.
            let coarse_keys_bytes = total_kmers.saturating_mul(std::mem::size_of::<u32>());
            let scratch_per_part = (total_kmers / NUM_PARTITIONS).max(1) * ps_bytes;
            let scratch_bytes = threads.saturating_mul(scratch_per_part + sub_keys * 8);
            let transient_build = coarse_keys_bytes.saturating_add(scratch_bytes);
            offsets_bytes
                .saturating_add(ps_total)
                .saturating_add(transient_build.max(rev_bytes))
        } else {
            let p = choose_bucket_bits(total_kmers);
            let sub_buckets = 1usize << (p - 8);
            let offsets_bytes = NUM_PARTITIONS * (sub_buckets + 1) * std::mem::size_of::<u32>();
            let keys_bytes = total_kmers.saturating_mul(std::mem::size_of::<u64>());
            let ps_total = total_kmers.saturating_mul(ps_bytes);
            let scratch_per_part =
                (total_kmers / NUM_PARTITIONS).max(1) * (std::mem::size_of::<u64>() + ps_bytes);
            let scratch_bytes = threads.saturating_mul(scratch_per_part + sub_buckets * 8);
            offsets_bytes
                .saturating_add(keys_bytes)
                .saturating_add(ps_total)
                .saturating_add(scratch_bytes.max(rev_bytes))
        }
    }

    /// Estimate exact peak RAM in bytes required to build and hold `RefIndex` (plus `rev_ref` if
    /// `compare_reverse` is true). Returns `usize::MAX` if the reference exceeds `u32` offset
    /// capacity (`>= u32::MAX` valid k-mers).
    pub fn estimate_bytes(
        reference: &Genome,
        word_length: usize,
        soft_mask: bool,
        compare_reverse: bool,
    ) -> usize {
        let total_kmers = count_valid_kmers(reference, word_length, soft_mask);
        Self::estimate_bytes_for_kmers(reference, total_kmers, word_length, compare_reverse)
    }

    /// Check whether `RefIndex` fits within `memory_budget`, short-circuiting in O(1) without
    /// scanning `reference.data` whenever the `total_bases` upper bound already fits.
    pub fn fits_in_budget(
        reference: &Genome,
        word_length: usize,
        soft_mask: bool,
        compare_reverse: bool,
        memory_budget: usize,
    ) -> bool {
        let upper_kmers = (reference.total_bases as usize).saturating_sub(word_length - 1);
        if Self::estimate_bytes_for_kmers(reference, upper_kmers, word_length, compare_reverse)
            <= memory_budget
        {
            return true;
        }
        Self::estimate_bytes(reference, word_length, soft_mask, compare_reverse) <= memory_budget
    }

    /// Build a canonical reference index directly at `word_length` (`4..=32`).
    pub fn build(
        reference: &Genome,
        word_length: usize,
        freq_limit: u64,
        soft_mask: bool,
    ) -> Result<Self> {
        if !(4..=32).contains(&word_length) {
            return Err(GeckoError::InvalidParameter(format!(
                "word_length must be between 4 and 32, got {word_length}"
            )));
        }

        let use_u32_ps = reference.data.len() < (1usize << 31);
        if use_direct_csr(word_length, reference.total_bases) {
            if use_u32_ps {
                build_direct_csr_u32(reference, word_length, freq_limit, soft_mask)
            } else {
                build_direct_csr_u64(reference, word_length, freq_limit, soft_mask)
            }
        } else if use_u32_ps {
            build_hash_csr_u32(reference, word_length, freq_limit, soft_mask)
        } else {
            build_hash_csr_u64(reference, word_length, freq_limit, soft_mask)
        }
    }

    /// Approximate heap memory currently held by the built index.
    pub fn heap_bytes(&self) -> usize {
        match &self.kind {
            RefIndexKind::DirectCsr {
                offsets,
                pos_strand,
                ..
            } => offsets.capacity() * std::mem::size_of::<u32>() + pos_strand.heap_bytes(),
            RefIndexKind::HashCsr {
                offsets,
                keys,
                pos_strand,
                ..
            } => {
                offsets.capacity() * std::mem::size_of::<u32>()
                    + keys.capacity() * std::mem::size_of::<u64>()
                    + pos_strand.heap_bytes()
            }
        }
    }

    /// Whether this index uses Direct CSR mode (zero key storage).
    #[inline(always)]
    pub fn is_direct(&self) -> bool {
        matches!(self.kind, RefIndexKind::DirectCsr { .. })
    }

    /// Compute the `offsets` index `b_idx` for a canonical k-mer key.
    #[inline(always)]
    pub(crate) fn bucket_index(&self, canon_key: u64) -> usize {
        match &self.kind {
            RefIndexKind::DirectCsr { sub_bits, .. } => {
                let b = canon_key as usize;
                b + (b >> *sub_bits)
            }
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                ..
            } => {
                let b = (fmix64(canon_key) >> *hash_shift) as usize;
                b + (b >> *sub_bits)
            }
        }
    }

    /// Stage 1 prefetch: prefetch `offsets[b_idx]`.
    #[inline(always)]
    pub(crate) fn prefetch_offset(&self, b_idx: usize) {
        let offsets = match &self.kind {
            RefIndexKind::DirectCsr { offsets, .. } => offsets,
            RefIndexKind::HashCsr { offsets, .. } => offsets,
        };
        prefetch_t0(unsafe { offsets.as_ptr().add(b_idx) });
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

    /// Stage 2 prefetch: prefetch `keys[start]` (Hash-CSR) or `pos_strand[start]` (Direct CSR).
    #[inline(always)]
    pub(crate) fn prefetch_bucket_data(&self, start: usize) {
        match &self.kind {
            RefIndexKind::DirectCsr { pos_strand, .. } => {
                pos_strand.prefetch_at(start);
            }
            RefIndexKind::HashCsr { keys, .. } => {
                if start < keys.len() {
                    prefetch_t0(unsafe { keys.as_ptr().add(start) });
                }
            }
        }
    }

    /// Stage 3 probe: given bucket range `[start, end)` and `canon_key`, return the sub-range
    /// `[lo, hi)` of exact matches in `pos_strand` (and reference to `pos_strand`).
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
                keys, pos_strand, ..
            } => {
                let (lo, hi) = probe_sorted_keys(keys, start, end, canon_key, use_avx2);
                (lo, hi, pos_strand)
            }
        }
    }

    /// Convenience single-key lookup returning all matching `(pos_y, strand_bit)` entries.
    pub fn lookup(&self, canon_key: u64) -> Vec<(u64, u8)> {
        let b_idx = self.bucket_index(canon_key);
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

/// Probe a sorted bucket `keys[start..end]` for exact occurrences of `key`.
#[inline(always)]
fn probe_sorted_keys(
    keys: &[u64],
    start: usize,
    end: usize,
    key: u64,
    use_avx2: bool,
) -> (usize, usize) {
    let len = end - start;
    if len == 0 {
        return (start, start);
    }
    if len <= 16 {
        #[cfg(target_arch = "x86_64")]
        if use_avx2 && len >= 4 {
            // Safety: `use_avx2` verified by `is_x86_feature_detected!("avx2")`.
            return unsafe { probe_bucket_avx2(keys, start, end, key) };
        }
        let _ = use_avx2;
        return probe_bucket_linear(keys, start, end, key);
    }

    // Binary search for larger buckets (> 16 entries).
    let slice = &keys[start..end];
    let rel_lo = slice.partition_point(|&k| k < key);
    if rel_lo >= slice.len() || slice[rel_lo] != key {
        return (start, start);
    }
    let mut rel_hi = rel_lo + 1;
    if rel_hi < slice.len() && slice[rel_hi] == key {
        rel_hi = rel_lo + slice[rel_lo..].partition_point(|&k| k == key);
    }
    (start + rel_lo, start + rel_hi)
}

#[inline(always)]
fn probe_bucket_linear(keys: &[u64], start: usize, end: usize, key: u64) -> (usize, usize) {
    let mut i = start;
    while i < end {
        let k = unsafe { *keys.get_unchecked(i) };
        if k >= key {
            if k > key {
                return (start, start);
            }
            let lo = i;
            i += 1;
            while i < end && unsafe { *keys.get_unchecked(i) } == key {
                i += 1;
            }
            return (lo, i);
        }
        i += 1;
    }
    (start, start)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn probe_bucket_avx2(keys: &[u64], start: usize, end: usize, key: u64) -> (usize, usize) {
    use core::arch::x86_64::*;

    let target = _mm256_set1_epi64x(key as i64);
    let mut i = start;
    let mut shift = 0u32;
    let mut mask = 0u32;

    while i + 4 <= end {
        let v = unsafe { _mm256_loadu_si256(keys.as_ptr().add(i) as *const __m256i) };
        let cmp = _mm256_cmpeq_epi64(v, target);
        let m = _mm256_movemask_pd(_mm256_castsi256_pd(cmp)) as u32;
        mask |= m << shift;
        i += 4;
        shift += 4;
    }

    while i < end {
        if unsafe { *keys.get_unchecked(i) } == key {
            mask |= 1u32 << shift;
        }
        i += 1;
        shift += 1;
    }

    if mask == 0 {
        (start, start)
    } else {
        let lo = start + mask.trailing_zeros() as usize;
        let hi = start + (32 - mask.leading_zeros()) as usize;
        (lo, hi)
    }
}

// ---------------------------------------------------------------------------
// Parallel chunking & pointer helpers
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct SendMutPtr<T>(*mut T);

unsafe impl<T: Send> Send for SendMutPtr<T> {}
unsafe impl<T: Send> Sync for SendMutPtr<T> {}

/// Partition `[0 .. data_len.saturating_sub(wl - 1))` into balanced chunks for Pass 1.
fn make_chunk_ranges(data_len: usize, wl: usize) -> Vec<(usize, usize)> {
    let max_starts = data_len.saturating_sub(wl - 1);
    if max_starts == 0 {
        return Vec::new();
    }
    let threads = rayon::current_num_threads().max(1);
    let target_chunks = (threads * 4).clamp(1, 64);
    let chunk_len = max_starts.div_ceil(target_chunks).max(65_536);

    let mut ranges = Vec::with_capacity(target_chunks + 1);
    let mut s = 0usize;
    while s < max_starts {
        let e = (s + chunk_len).min(max_starts);
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

/// Compute global partition offsets and per-chunk partition write offsets from chunk histograms.
fn compute_partition_layout(
    chunk_counts: &[[usize; NUM_PARTITIONS]],
) -> (Vec<usize>, Vec<[usize; NUM_PARTITIONS]>, usize) {
    let mut part_totals = [0usize; NUM_PARTITIONS];
    for counts in chunk_counts {
        for p in 0..NUM_PARTITIONS {
            part_totals[p] += counts[p];
        }
    }

    let mut part_offsets = vec![0usize; NUM_PARTITIONS + 1];
    for p in 0..NUM_PARTITIONS {
        part_offsets[p + 1] = part_offsets[p] + part_totals[p];
    }
    let total = part_offsets[NUM_PARTITIONS];

    let mut chunk_write = vec![[0usize; NUM_PARTITIONS]; chunk_counts.len()];
    let mut running = [0usize; NUM_PARTITIONS];
    running.copy_from_slice(&part_offsets[..NUM_PARTITIONS]);
    for (c, counts) in chunk_counts.iter().enumerate() {
        chunk_write[c].copy_from_slice(&running);
        for p in 0..NUM_PARTITIONS {
            running[p] += counts[p];
        }
    }

    (part_offsets, chunk_write, total)
}

#[inline(always)]
fn kmer_mask(wl: usize) -> u64 {
    if wl == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * wl)) - 1
    }
}

/// Inline helper: iterate valid canonical k-mers whose start position lies in `[start..end)`.
#[inline(always)]
fn for_each_canonical_kmer<F>(
    data: &[u8],
    start: usize,
    end: usize,
    wl: usize,
    lut: &[u8; 256],
    mut f: F,
) where
    F: FnMut(u64, u64, u64),
{
    let mask = kmer_mask(wl);
    let rc_shift = (2 * (wl - 1)) as u32;
    let scan_end = (end + wl - 1).min(data.len());
    let mut kf = 0u64;
    let mut kr = 0u64;
    let mut valid_run = 0usize;

    for (idx, &raw) in data[start..scan_end].iter().enumerate() {
        let b = lut[raw as usize];
        if b != 0xFF {
            kf = ((kf << 2) | (b as u64)) & mask;
            kr = (kr >> 2) | (((b ^ 3) as u64) << rc_shift);
            valid_run += 1;
            if valid_run >= wl {
                let pos = (start + idx + 1 - wl) as u64;
                if kr < kf {
                    f(kr, pos, 1);
                } else {
                    f(kf, pos, 0);
                }
            }
        } else {
            valid_run = 0;
            kf = 0;
            kr = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// Direct CSR builder (generic over u32 / u64 packed pos_strand)
// ---------------------------------------------------------------------------

trait PackedPs: Copy + Default + Ord + std::fmt::Debug + Send + Sync + 'static {
    fn pack(pos: u64, strand: u64) -> Self;
    fn strand_bit(self) -> u64;
    fn into_vec(v: Vec<Self>) -> PosStrandVec;
}

impl PackedPs for u32 {
    #[inline(always)]
    fn pack(pos: u64, strand: u64) -> Self {
        ((pos << 1) | strand) as u32
    }

    #[inline(always)]
    fn strand_bit(self) -> u64 {
        (self & 1) as u64
    }

    fn into_vec(v: Vec<Self>) -> PosStrandVec {
        PosStrandVec::U32(v)
    }
}

impl PackedPs for u64 {
    #[inline(always)]
    fn pack(pos: u64, strand: u64) -> Self {
        (pos << 1) | strand
    }

    #[inline(always)]
    fn strand_bit(self) -> u64 {
        self & 1
    }

    fn into_vec(v: Vec<Self>) -> PosStrandVec {
        PosStrandVec::U64(v)
    }
}

fn build_direct_csr_u32(
    reference: &Genome,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
) -> Result<RefIndex> {
    build_direct_csr::<u32>(reference, wl, freq_limit, soft_mask)
}

fn build_direct_csr_u64(
    reference: &Genome,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
) -> Result<RefIndex> {
    build_direct_csr::<u64>(reference, wl, freq_limit, soft_mask)
}

fn build_direct_csr<P: PackedPs>(
    reference: &Genome,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
) -> Result<RefIndex> {
    let data = &reference.data;
    let lut = base_lut(soft_mask);
    let sub_bits = (2 * wl - 8) as u32;
    let sub_keys = 1usize << sub_bits;
    let sub_mask = sub_keys - 1;

    let ranges = make_chunk_ranges(data.len(), wl);
    if ranges.is_empty() {
        return Ok(RefIndex {
            word_length: wl,
            freq_limit,
            soft_mask,
            ref_data_len: data.len(),
            total_kmers: 0,
            kept_entries: 0,
            kind: RefIndexKind::DirectCsr {
                sub_bits,
                offsets: vec![0u32; NUM_PARTITIONS * (sub_keys + 1)],
                pos_strand: P::into_vec(Vec::new()),
            },
        });
    }

    // Pass 1a: parallel count into 256 coarse partitions.
    let chunk_counts: Vec<[usize; NUM_PARTITIONS]> = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            for_each_canonical_kmer(data, start, end, wl, lut, |canon, _, _| {
                counts[(canon >> sub_bits) as usize] += 1;
            });
            counts
        })
        .collect();

    let (part_offsets, chunk_write, total_kmers) = compute_partition_layout(&chunk_counts);
    if total_kmers >= u32::MAX as usize {
        return Err(GeckoError::InvalidParameter(
            "reference has >= 2^32 valid k-mers; use out-of-core --memory fallback".into(),
        ));
    }

    // Pass 1b: parallel scatter with L1/L2-resident per-partition staging buffers.
    let mut coarse_keys = vec![0u32; total_kmers];
    let mut pos_strand = vec![P::default(); total_kmers];
    let keys_ptr = SendMutPtr(coarse_keys.as_mut_ptr());
    let ps_ptr = SendMutPtr(pos_strand.as_mut_ptr());

    ranges
        .par_iter()
        .zip(chunk_write.par_iter())
        .for_each(|(&(start, end), &init_write)| {
            let kptr = keys_ptr;
            let pptr = ps_ptr;
            let mut write_pos = init_write;
            let mut buf_k: Box<[[u32; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[0u32; STAGE_CAP]; NUM_PARTITIONS]
                    .into_boxed_slice()
                    .try_into()
                    .unwrap();
            let mut buf_p: Box<[[P; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[P::default(); STAGE_CAP]; NUM_PARTITIONS]
                    .into_boxed_slice()
                    .try_into()
                    .unwrap();
            let mut buf_len = [0u8; NUM_PARTITIONS];

            for_each_canonical_kmer(data, start, end, wl, lut, |canon, pos, strand| {
                let p = (canon >> sub_bits) as usize;
                let idx = buf_len[p] as usize;
                buf_k[p][idx] = (canon as usize & sub_mask) as u32;
                buf_p[p][idx] = P::pack(pos, strand);
                let next_len = idx + 1;
                if next_len == STAGE_CAP {
                    let dst = write_pos[p];
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            buf_k[p].as_ptr(),
                            kptr.0.add(dst),
                            STAGE_CAP,
                        );
                        std::ptr::copy_nonoverlapping(
                            buf_p[p].as_ptr(),
                            pptr.0.add(dst),
                            STAGE_CAP,
                        );
                    }
                    write_pos[p] = dst + STAGE_CAP;
                    buf_len[p] = 0;
                } else {
                    buf_len[p] = next_len as u8;
                }
            });

            for p in 0..NUM_PARTITIONS {
                let rem = buf_len[p] as usize;
                if rem > 0 {
                    let dst = write_pos[p];
                    unsafe {
                        std::ptr::copy_nonoverlapping(buf_k[p].as_ptr(), kptr.0.add(dst), rem);
                        std::ptr::copy_nonoverlapping(buf_p[p].as_ptr(), pptr.0.add(dst), rem);
                    }
                }
            }
        });

    // Pass 2: Rayon-parallel fine-bucketing + per-strand freq_limit filtering + in-place compaction.
    let mut offsets = vec![0u32; NUM_PARTITIONS * (sub_keys + 1)];
    let part_k_slices = split_by_offsets_mut(&mut coarse_keys, &part_offsets);
    let part_p_slices = split_by_offsets_mut(&mut pos_strand, &part_offsets);
    let part_off_slices: Vec<&mut [u32]> = offsets.chunks_exact_mut(sub_keys + 1).collect();

    let kept_entries: usize = (0..NUM_PARTITIONS)
        .into_par_iter()
        .zip(part_k_slices.into_par_iter())
        .zip(part_p_slices.into_par_iter())
        .zip(part_off_slices.into_par_iter())
        .map_init(
            || (vec![0u32; sub_keys], vec![0u32; sub_keys], Vec::<P>::new()),
            |(sub_counts, sub_write, scratch_ps), (((p, pk), pp), poff)| {
                let part_start = part_offsets[p];
                let n = pk.len();
                if n == 0 {
                    for slot in poff.iter_mut() {
                        *slot = part_start as u32;
                    }
                    return 0usize;
                }

                sub_counts.fill(0);
                for &sub in pk.iter() {
                    sub_counts[sub as usize] += 1;
                }

                let mut acc = 0u32;
                for s in 0..sub_keys {
                    sub_write[s] = acc;
                    acc += sub_counts[s];
                }

                scratch_ps.resize(n, P::default());
                for i in 0..n {
                    let sub = pk[i] as usize;
                    let dst = sub_write[sub] as usize;
                    scratch_ps[dst] = pp[i];
                    sub_write[sub] += 1;
                }

                // Entries within each sub-bucket are already in ascending `pos` order because
                // Pass 1 preserved chunk order and within-chunk scan order.
                let mut write_idx = 0usize;
                let mut sub_start = 0usize;
                for s in 0..sub_keys {
                    poff[s] = (part_start + write_idx) as u32;
                    let cnt = sub_counts[s] as usize;
                    if cnt == 0 {
                        continue;
                    }
                    let sub_end = sub_start + cnt;
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
                        pp[write_idx..write_idx + cnt]
                            .copy_from_slice(&scratch_ps[sub_start..sub_end]);
                        write_idx += cnt;
                    } else if keep0 || keep1 {
                        for &ps in &scratch_ps[sub_start..sub_end] {
                            let is_keep = if ps.strand_bit() == 0 { keep0 } else { keep1 };
                            if is_keep {
                                pp[write_idx] = ps;
                                write_idx += 1;
                            }
                        }
                    }
                    sub_start = sub_end;
                }
                poff[sub_keys] = (part_start + write_idx) as u32;
                write_idx
            },
        )
        .sum();

    drop(coarse_keys);

    Ok(RefIndex {
        word_length: wl,
        freq_limit,
        soft_mask,
        ref_data_len: data.len(),
        total_kmers,
        kept_entries,
        kind: RefIndexKind::DirectCsr {
            sub_bits,
            offsets,
            pos_strand: P::into_vec(pos_strand),
        },
    })
}

// ---------------------------------------------------------------------------
// 2-Level Parallel Radix-Partitioned Hash-CSR builder
// ---------------------------------------------------------------------------

fn build_hash_csr_u32(
    reference: &Genome,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
) -> Result<RefIndex> {
    build_hash_csr::<u32>(reference, wl, freq_limit, soft_mask)
}

fn build_hash_csr_u64(
    reference: &Genome,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
) -> Result<RefIndex> {
    build_hash_csr::<u64>(reference, wl, freq_limit, soft_mask)
}

fn build_hash_csr<P: PackedPs>(
    reference: &Genome,
    wl: usize,
    freq_limit: u64,
    soft_mask: bool,
) -> Result<RefIndex> {
    let data = &reference.data;
    let lut = base_lut(soft_mask);

    let ranges = make_chunk_ranges(data.len(), wl);
    if ranges.is_empty() {
        let bucket_bits = 14u32;
        let sub_bits = bucket_bits - 8;
        let sub_buckets = 1usize << sub_bits;
        return Ok(RefIndex {
            word_length: wl,
            freq_limit,
            soft_mask,
            ref_data_len: data.len(),
            total_kmers: 0,
            kept_entries: 0,
            kind: RefIndexKind::HashCsr {
                sub_bits,
                hash_shift: 64 - bucket_bits,
                offsets: vec![0u32; NUM_PARTITIONS * (sub_buckets + 1)],
                keys: Vec::new(),
                pos_strand: P::into_vec(Vec::new()),
            },
        });
    }

    // Pass 1a: parallel count into 256 coarse partitions (top 8 bits of fmix64).
    let chunk_counts: Vec<[usize; NUM_PARTITIONS]> = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            for_each_canonical_kmer(data, start, end, wl, lut, |canon, _, _| {
                let p = (fmix64(canon) >> 56) as usize;
                counts[p] += 1;
            });
            counts
        })
        .collect();

    let (part_offsets, chunk_write, total_kmers) = compute_partition_layout(&chunk_counts);
    if total_kmers >= u32::MAX as usize {
        return Err(GeckoError::InvalidParameter(
            "reference has >= 2^32 valid k-mers; use out-of-core --memory fallback".into(),
        ));
    }

    let bucket_bits = choose_bucket_bits(total_kmers);
    let sub_bits = bucket_bits - 8;
    let sub_buckets = 1usize << sub_bits;
    let sub_mask = sub_buckets - 1;
    let hash_shift = 64 - bucket_bits;

    // Pass 1b: parallel scatter into coarse partitions using 64-entry staging buffers.
    let mut keys = vec![0u64; total_kmers];
    let mut pos_strand = vec![P::default(); total_kmers];
    let keys_ptr = SendMutPtr(keys.as_mut_ptr());
    let ps_ptr = SendMutPtr(pos_strand.as_mut_ptr());

    ranges
        .par_iter()
        .zip(chunk_write.par_iter())
        .for_each(|(&(start, end), &init_write)| {
            let kptr = keys_ptr;
            let pptr = ps_ptr;
            let mut write_pos = init_write;
            let mut buf_k: Box<[[u64; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[0u64; STAGE_CAP]; NUM_PARTITIONS]
                    .into_boxed_slice()
                    .try_into()
                    .unwrap();
            let mut buf_p: Box<[[P; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[P::default(); STAGE_CAP]; NUM_PARTITIONS]
                    .into_boxed_slice()
                    .try_into()
                    .unwrap();
            let mut buf_len = [0u8; NUM_PARTITIONS];

            for_each_canonical_kmer(data, start, end, wl, lut, |canon, pos, strand| {
                let p = (fmix64(canon) >> 56) as usize;
                let idx = buf_len[p] as usize;
                buf_k[p][idx] = canon;
                buf_p[p][idx] = P::pack(pos, strand);
                let next_len = idx + 1;
                if next_len == STAGE_CAP {
                    let dst = write_pos[p];
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            buf_k[p].as_ptr(),
                            kptr.0.add(dst),
                            STAGE_CAP,
                        );
                        std::ptr::copy_nonoverlapping(
                            buf_p[p].as_ptr(),
                            pptr.0.add(dst),
                            STAGE_CAP,
                        );
                    }
                    write_pos[p] = dst + STAGE_CAP;
                    buf_len[p] = 0;
                } else {
                    buf_len[p] = next_len as u8;
                }
            });

            for p in 0..NUM_PARTITIONS {
                let rem = buf_len[p] as usize;
                if rem > 0 {
                    let dst = write_pos[p];
                    unsafe {
                        std::ptr::copy_nonoverlapping(buf_k[p].as_ptr(), kptr.0.add(dst), rem);
                        std::ptr::copy_nonoverlapping(buf_p[p].as_ptr(), pptr.0.add(dst), rem);
                    }
                }
            }
        });

    // Pass 2: Rayon-parallel fine-bucketing, tiny-bucket sort, per-strand freq_limit filter,
    // and in-place compaction per partition.
    let mut offsets = vec![0u32; NUM_PARTITIONS * (sub_buckets + 1)];
    let part_k_slices = split_by_offsets_mut(&mut keys, &part_offsets);
    let part_p_slices = split_by_offsets_mut(&mut pos_strand, &part_offsets);
    let part_off_slices: Vec<&mut [u32]> = offsets.chunks_exact_mut(sub_buckets + 1).collect();

    let kept_entries: usize = (0..NUM_PARTITIONS)
        .into_par_iter()
        .zip(part_k_slices.into_par_iter())
        .zip(part_p_slices.into_par_iter())
        .zip(part_off_slices.into_par_iter())
        .map_init(
            || {
                (
                    vec![0u32; sub_buckets],
                    vec![0u32; sub_buckets],
                    Vec::<u64>::new(),
                    Vec::<P>::new(),
                    Vec::<(u64, P)>::new(),
                )
            },
            |(sub_counts, sub_write, scratch_k, scratch_p, sort_buf), (((p, pk), pp), poff)| {
                let part_start = part_offsets[p];
                let n = pk.len();
                if n == 0 {
                    for slot in poff.iter_mut() {
                        *slot = part_start as u32;
                    }
                    return 0usize;
                }

                sub_counts.fill(0);
                for &k in pk.iter() {
                    let sub = ((fmix64(k) >> hash_shift) as usize) & sub_mask;
                    sub_counts[sub] += 1;
                }

                let mut acc = 0u32;
                for s in 0..sub_buckets {
                    sub_write[s] = acc;
                    acc += sub_counts[s];
                }

                scratch_k.resize(n, 0);
                scratch_p.resize(n, P::default());
                for i in 0..n {
                    let k = pk[i];
                    let sub = ((fmix64(k) >> hash_shift) as usize) & sub_mask;
                    let dst = sub_write[sub] as usize;
                    scratch_k[dst] = k;
                    scratch_p[dst] = pp[i];
                    sub_write[sub] += 1;
                }

                let mut write_idx = 0usize;
                let mut sub_start = 0usize;
                for s in 0..sub_buckets {
                    poff[s] = (part_start + write_idx) as u32;
                    let cnt = sub_counts[s] as usize;
                    if cnt == 0 {
                        continue;
                    }
                    let sub_end = sub_start + cnt;
                    if cnt == 1 {
                        if freq_limit >= 1 {
                            pk[write_idx] = scratch_k[sub_start];
                            pp[write_idx] = scratch_p[sub_start];
                            write_idx += 1;
                        }
                        sub_start = sub_end;
                        continue;
                    }

                    // Sort bucket by (key, pos_strand).
                    if cnt <= 24 {
                        insertion_sort_pair(
                            &mut scratch_k[sub_start..sub_end],
                            &mut scratch_p[sub_start..sub_end],
                        );
                    } else {
                        sort_buf.clear();
                        sort_buf.reserve(cnt);
                        for idx in sub_start..sub_end {
                            sort_buf.push((scratch_k[idx], scratch_p[idx]));
                        }
                        sort_buf.sort_unstable();
                        for (off, &(k, ps)) in sort_buf.iter().enumerate() {
                            scratch_k[sub_start + off] = k;
                            scratch_p[sub_start + off] = ps;
                        }
                    }

                    // Per-strand freq_limit filtering across runs of identical keys.
                    let mut r = sub_start;
                    while r < sub_end {
                        let key = scratch_k[r];
                        let mut r_end = r + 1;
                        let mut cnt0 = 0u64;
                        let mut cnt1 = 0u64;
                        if scratch_p[r].strand_bit() == 0 {
                            cnt0 += 1;
                        } else {
                            cnt1 += 1;
                        }
                        while r_end < sub_end && scratch_k[r_end] == key {
                            if scratch_p[r_end].strand_bit() == 0 {
                                cnt0 += 1;
                            } else {
                                cnt1 += 1;
                            }
                            r_end += 1;
                        }
                        let keep0 = cnt0 <= freq_limit;
                        let keep1 = cnt1 <= freq_limit;
                        if keep0 && keep1 {
                            let rlen = r_end - r;
                            pk[write_idx..write_idx + rlen].copy_from_slice(&scratch_k[r..r_end]);
                            pp[write_idx..write_idx + rlen].copy_from_slice(&scratch_p[r..r_end]);
                            write_idx += rlen;
                        } else if keep0 || keep1 {
                            for &ps in &scratch_p[r..r_end] {
                                let is_keep = if ps.strand_bit() == 0 { keep0 } else { keep1 };
                                if is_keep {
                                    pk[write_idx] = key;
                                    pp[write_idx] = ps;
                                    write_idx += 1;
                                }
                            }
                        }
                        r = r_end;
                    }

                    sub_start = sub_end;
                }
                poff[sub_buckets] = (part_start + write_idx) as u32;
                write_idx
            },
        )
        .sum();

    Ok(RefIndex {
        word_length: wl,
        freq_limit,
        soft_mask,
        ref_data_len: data.len(),
        total_kmers,
        kept_entries,
        kind: RefIndexKind::HashCsr {
            sub_bits,
            hash_shift,
            offsets,
            keys,
            pos_strand: P::into_vec(pos_strand),
        },
    })
}

/// Fast in-place insertion sort of paired `(keys, ps)` slices by `(key, ps)` for `n <= 24`.
#[inline(always)]
fn insertion_sort_pair<P: Copy + Ord>(keys: &mut [u64], ps: &mut [P]) {
    let n = keys.len();
    for i in 1..n {
        let k = keys[i];
        let p = ps[i];
        let mut j = i;
        while j > 0 && (keys[j - 1] > k || (keys[j - 1] == k && ps[j - 1] > p)) {
            keys[j] = keys[j - 1];
            ps[j] = ps[j - 1];
            j -= 1;
        }
        keys[j] = k;
        ps[j] = p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_kmer(s: &[u8]) -> (u64, u64) {
        let lut = base_lut(true);
        let wl = s.len();
        let mask = kmer_mask(wl);
        let rc_shift = (2 * (wl - 1)) as u32;
        let mut kf = 0u64;
        let mut kr = 0u64;
        for &b in s {
            let v = lut[b as usize];
            kf = ((kf << 2) | (v as u64)) & mask;
            kr = (kr >> 2) | (((v ^ 3) as u64) << rc_shift);
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
    fn hash_csr_probe_matches_linear_and_binary_search() {
        // Test bucket probing across sizes 1..=30 including AVX2 path.
        let mut keys = Vec::new();
        for k in 0u64..10 {
            for _ in 0..3 {
                keys.push(k * 100 + 7);
            }
        }
        #[cfg(target_arch = "x86_64")]
        let has_avx2 = std::is_x86_feature_detected!("avx2");
        #[cfg(not(target_arch = "x86_64"))]
        let has_avx2 = false;

        for end in 1..=keys.len() {
            for target in [7u64, 107, 407, 907, 50, 9999] {
                let expected = probe_bucket_linear(&keys, 0, end, target);
                let actual = probe_sorted_keys(&keys, 0, end, target, has_avx2);
                assert_eq!(actual, expected, "mismatch for end={end}, target={target}");
            }
        }
    }
}
