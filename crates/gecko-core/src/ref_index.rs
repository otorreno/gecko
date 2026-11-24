use std::mem::MaybeUninit;

use rayon::prelude::*;

use crate::error::{GeckoError, Result};
use crate::seq::{Genome, PackedSeq};

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
    },
}

/// Single-reference canonical k-mer index supporting lock-free streaming queries.
///
/// Indexes only Reference Y once at `word_length` (`4..=32`) using canonical keys
/// `min(kf, kr)` and a 1-bit strand flag `(kr < kf)`, with per-strand `freq_limit` filtering.
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
    let total_kmers = count_valid_kmers(reference, word_length, soft_mask);
    if total_kmers == 0 {
        return Vec::new();
    }
    let cap = max_kmers_per_slab.max(1);
    if total_kmers <= cap {
        return vec![(0, NUM_PARTITIONS as u16)];
    }

    if !use_direct_csr_for(word_length, reference.total_bases, spaced_mask) {
        // Hash-CSR partitions by top 8 bits of `fmix64(canon)`, which is uniform to <0.2%.
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

    // Direct-CSR partitions by top 8 bits of `min(kf, kr)` (non-uniform across 0..256),
    // so compute the exact 256-bucket histogram and greedily pack contiguous ranges <= cap.
    let seq = &reference.seq;
    let sub_bits = (2 * word_length - 8) as u32;
    let max_starts = seq.len.saturating_sub(word_length - 1);
    let ranges = make_chunk_ranges(0, max_starts, seq.len, word_length);
    let part_totals = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            for_each_canonical_kmer_2bit(seq, start, end, word_length, soft_mask, |canon, _, _| {
                counts[(canon >> sub_bits) as usize] += 1;
            });
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
        let threads = rayon::current_num_threads().max(1);

        if use_direct_csr_for(word_length, reference.total_bases, spaced_mask) {
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
            // Steady-state HashCsr stores 4-byte `tags: Vec<u32>`, while Pass 1b/2 holds transient
            // per-partition `u64` keys (freed partition-by-partition as `tags` pages are written).
            let peak_key_bytes = total_kmers.saturating_mul(std::mem::size_of::<u64>());
            let ps_total = total_kmers.saturating_mul(ps_bytes);
            let scratch_per_part =
                (total_kmers / NUM_PARTITIONS).max(1) * (std::mem::size_of::<u64>() + ps_bytes);
            let scratch_bytes = threads.saturating_mul(scratch_per_part + sub_buckets * 8);
            offsets_bytes
                .saturating_add(peak_key_bytes)
                .saturating_add(ps_total)
                .saturating_add(scratch_bytes.max(rev_bytes))
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
        // Peak build RAM per k-mer: Hash-CSR transient u64 keys (8 B, shrunk to 4 B tags) +
        // pos_strand + offsets/scratch headroom; Direct-CSR <= that.
        let bytes_per_kmer =
            if use_direct_csr_for(word_length, reference.total_bases, spaced_mask) {
                ps_bytes + 5
            } else {
                8 + ps_bytes + 3
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
        Self::build_partitions_with_mask(
            reference,
            0,
            NUM_PARTITIONS as u16,
            word_length,
            spaced_mask,
            freq_limit,
            soft_mask,
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
            )
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
                tags,
                pos_strand,
                ..
            } => {
                offsets.capacity() * std::mem::size_of::<u32>()
                    + tags.capacity() * std::mem::size_of::<u32>()
                    + pos_strand.heap_bytes()
            }
        }
    }

    /// Whether this index uses Direct CSR mode (zero key storage).
    #[inline(always)]
    pub fn is_direct(&self) -> bool {
        matches!(self.kind, RefIndexKind::DirectCsr { .. })
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
                let b_local = b - (ps << *sub_bits);
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
unsafe fn probe_bucket_u32_avx2(tags: &[u32], start: usize, end: usize, tag: u32) -> (usize, usize) {
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

/// Iterate all valid canonical k-mers whose start position lies in `[start..end)` directly from
/// 2-bit packed `u64` words (jumping over `N` and soft-masked gaps in O(1), zero LUT lookups).
#[inline(always)]
fn for_each_canonical_kmer_2bit<F>(
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
        for_each_canonical_kmer_2bit_impl::<true, F>(seq, start, end, 32, soft_mask, f);
    } else {
        for_each_canonical_kmer_2bit_impl::<false, F>(seq, start, end, wl, soft_mask, f);
    }
}

#[inline(always)]
fn for_each_canonical_kmer_2bit_opt<F>(
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
        for_each_canonical_spaced_kmer_2bit(seq, start, end, wl, sm, soft_mask, f);
    } else {
        for_each_canonical_kmer_2bit(seq, start, end, wl, soft_mask, f);
    }
}

/// Iterate canonical spaced-seed k-mers using a palindromic 2-bit `spaced_mask`.
#[inline(always)]
fn for_each_canonical_spaced_kmer_2bit<F>(
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
                kr = (kr >> 2) | ((b ^ 3) << rc_shift);
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
                kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                let sf = kf & spaced_mask;
                let sr = kr & spaced_mask;
                if sr < sf {
                    f(sr, pos, 1);
                } else {
                    f(sf, pos, 0);
                }
            }
            idx += avail;
        }
    });
}

#[inline(always)]
fn for_each_canonical_kmer_2bit_impl<const WL32: bool, F>(
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
                kr = (kr >> 2) | ((b ^ 3) << rc_shift);
            }
            idx += avail;
        }

        // Emit one canonical k-mer per base through `re`, loading 32 bases per u64 word.
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
                kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                if kr < kf {
                    f(kr, pos, 1);
                } else {
                    f(kf, pos, 0);
                }
            }
            idx += avail;
        }
    });
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
) -> Result<RefIndex> {
    build_direct_csr::<u32>(
        reference, slab_start, slab_end, part_start, part_end, wl, freq_limit, soft_mask,
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
) -> Result<RefIndex> {
    build_direct_csr::<u64>(
        reference, slab_start, slab_end, part_start, part_end, wl, freq_limit, soft_mask,
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
        });
    }

    // Pass 1a: parallel count into active coarse partitions `[part_start .. part_end)`.
    let chunk_counts: Vec<[usize; NUM_PARTITIONS]> = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            if full_parts {
                for_each_canonical_kmer_2bit(seq, start, end, wl, soft_mask, |canon, _, _| {
                    counts[(canon >> sub_bits) as usize] += 1;
                });
            } else {
                for_each_canonical_kmer_2bit(seq, start, end, wl, soft_mask, |canon, _, _| {
                    let p = (canon >> sub_bits) as usize;
                    if p >= part_start && p < part_end {
                        counts[p - part_start] += 1;
                    }
                });
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

    // Pass 1b: parallel scatter with L1/L2-resident per-thread staging buffers.
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

            for_each_canonical_kmer_2bit(seq, start, end, wl, soft_mask, |canon, pos, strand| {
                let p_glob = (canon >> sub_bits) as usize;
                if !full_parts && (p_glob < part_start || p_glob >= part_end) {
                    return;
                }
                let p = if full_parts { p_glob } else { p_glob - part_start };
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

            for p in 0..n_parts {
                let rem = buf_len[p] as usize;
                if rem > 0 {
                    let dst = write_pos[p];
                    unsafe {
                        std::ptr::copy_nonoverlapping(buf_k[p].as_ptr(), kptr.0.add(dst), rem);
                        std::ptr::copy_nonoverlapping(buf_p[p].as_ptr(), pptr.0.add(dst), rem);
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

    // Pass 2: Rayon-parallel fine-bucketing + per-strand freq_limit filtering + in-place compaction.
    let off_len = n_parts * (sub_keys + 1);
    let mut offsets: Vec<u32> = Vec::with_capacity(off_len);
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
            || (vec![0u32; sub_keys], vec![0u32; sub_keys], Vec::<P>::new()),
            |(sub_counts, sub_write, scratch_ps), (((lp, pk), pp), poff)| {
                let part_off = part_offsets[lp];
                let n = pk.len();
                if n == 0 {
                    for slot in poff.iter_mut() {
                        slot.write(part_off as u32);
                    }
                    return 0usize;
                }

                sub_counts.fill(0);
                for &sub in pk.iter() {
                    unsafe {
                        *sub_counts.get_unchecked_mut(sub as usize) += 1;
                    }
                }

                let mut acc = 0u32;
                for s in 0..sub_keys {
                    unsafe {
                        *sub_write.get_unchecked_mut(s) = acc;
                        acc += *sub_counts.get_unchecked(s);
                    }
                }

                if scratch_ps.len() < n {
                    scratch_ps.resize(n, P::default());
                }
                for i in 0..n {
                    unsafe {
                        let sub = *pk.get_unchecked(i) as usize;
                        let dst = *sub_write.get_unchecked(sub) as usize;
                        *scratch_ps.get_unchecked_mut(dst) = *pp.get_unchecked(i);
                        *sub_write.get_unchecked_mut(sub) = (dst + 1) as u32;
                    }
                }

                let mut write_idx = 0usize;
                let mut sub_start = 0usize;
                for (slot, &cnt_u32) in poff.iter_mut().zip(sub_counts.iter()) {
                    slot.write((part_off + write_idx) as u32);
                    let cnt = cnt_u32 as usize;
                    if cnt == 0 {
                        continue;
                    }
                    let sub_end = sub_start + cnt;
                    if (cnt as u64) <= freq_limit {
                        pp[write_idx..write_idx + cnt]
                            .copy_from_slice(&scratch_ps[sub_start..sub_end]);
                        write_idx += cnt;
                    } else {
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
                    }
                    sub_start = sub_end;
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
    })
}

// ---------------------------------------------------------------------------
// 2-Level Parallel Radix-Partitioned Hash-CSR builder (32-bit tag compressed)
// ---------------------------------------------------------------------------

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
) -> Result<RefIndex> {
    let seq = &reference.seq;
    let n_parts = part_end - part_start;
    let full_parts = part_start == 0 && part_end == NUM_PARTITIONS;

    let ranges = make_chunk_ranges(slab_start, slab_end, seq.len, wl);
    if ranges.is_empty() {
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
            },
        });
    }

    // Pass 1a: parallel count into active coarse partitions `[part_start .. part_end)` (top 8 bits of fmix64).
    let chunk_counts: Vec<[usize; NUM_PARTITIONS]> = ranges
        .par_iter()
        .map(|&(start, end)| {
            let mut counts = [0usize; NUM_PARTITIONS];
            if full_parts {
                for_each_canonical_kmer_2bit_opt(
                    seq,
                    start,
                    end,
                    wl,
                    spaced_mask,
                    soft_mask,
                    |canon, _, _| {
                        let p = (fmix64(canon) >> 56) as usize;
                        counts[p] += 1;
                    },
                );
            } else {
                for_each_canonical_kmer_2bit_opt(
                    seq,
                    start,
                    end,
                    wl,
                    spaced_mask,
                    soft_mask,
                    |canon, _, _| {
                        let p = (fmix64(canon) >> 56) as usize;
                        if p >= part_start && p < part_end {
                            counts[p - part_start] += 1;
                        }
                    },
                );
            }
            counts
        })
        .collect();

    let (part_offsets, part_totals, chunk_write, total_kmers) =
        compute_partition_layout(&chunk_counts, n_parts);
    if total_kmers >= u32::MAX as usize {
        return Err(GeckoError::InvalidParameter(
            "reference slab has >= 2^32 valid k-mers; use smaller --memory budget".into(),
        ));
    }

    // Scale `total_kmers` to its full-256-partition equivalent so per-sub-bucket load stays ~3-4
    // entries even when building a subset of partitions for multi-slab streaming.
    let equiv_total = total_kmers.saturating_mul(NUM_PARTITIONS) / n_parts.max(1);
    let bucket_bits = choose_bucket_bits(equiv_total);
    let sub_bits = bucket_bits - 8;
    let sub_buckets = 1usize << sub_bits;
    let sub_mask = sub_buckets - 1;
    let hash_shift = 64 - bucket_bits;

    let sub_shift_hi = 32 - bucket_bits;

    // Pass 1b: scatter rotated 64-bit hashes `h_rot = fmix64(canon).rotate_left(32)` (with the
    // 32-bit `tag` in the high 32 bits and bucket bits in the low 32 bits) into per-partition
    // `Vec<u64>`s. Because `fmix64` + rotate is a 64-bit bijection, Pass 2 needs ZERO `fmix64`
    // recomputations and only ONE sort per bucket to both enforce exact 64-bit `freq_limit` and
    // emit entries already sorted by ascending `tag`.
    let mut part_keys: Vec<Vec<u64>> = part_totals
        .iter()
        .map(|&n| Vec::with_capacity(n))
        .collect();
    let key_ptrs: Vec<SendMutPtr<u64>> = part_keys
        .iter_mut()
        .map(|v| SendMutPtr(v.as_mut_ptr()))
        .collect();
    let mut pos_strand: Vec<P> = Vec::with_capacity(total_kmers);
    let ps_ptr = SendMutPtr(pos_strand.as_mut_ptr());

    ranges.par_iter().zip(chunk_write.par_iter()).for_each_init(
        || {
            let buf_k: Box<[[u64; STAGE_CAP]; NUM_PARTITIONS]> =
                vec![[0u64; STAGE_CAP]; NUM_PARTITIONS]
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
            let pptr = ps_ptr;
            let mut write_pos = init_write;
            let mut buf_len = [0u8; NUM_PARTITIONS];

            for_each_canonical_kmer_2bit_opt(
                seq,
                start,
                end,
                wl,
                spaced_mask,
                soft_mask,
                |canon, pos, strand| {
                let h = fmix64(canon);
                let p_glob = (h >> 56) as usize;
                if !full_parts && (p_glob < part_start || p_glob >= part_end) {
                    return;
                }
                let lp = if full_parts {
                    p_glob
                } else {
                    p_glob - part_start
                };
                let idx = buf_len[lp] as usize;
                buf_k[lp][idx] = h.rotate_left(32);
                buf_p[lp][idx] = P::pack(pos, strand);
                let next_len = idx + 1;
                if next_len == STAGE_CAP {
                    let dst = write_pos[lp];
                    let local_dst = dst - part_offsets[lp];
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            buf_k[lp].as_ptr(),
                            key_ptrs[lp].0.add(local_dst),
                            STAGE_CAP,
                        );
                        std::ptr::copy_nonoverlapping(
                            buf_p[lp].as_ptr(),
                            pptr.0.add(dst),
                            STAGE_CAP,
                        );
                    }
                    write_pos[lp] = dst + STAGE_CAP;
                    buf_len[lp] = 0;
                } else {
                    buf_len[lp] = next_len as u8;
                }
            });

            for lp in 0..n_parts {
                let rem = buf_len[lp] as usize;
                if rem > 0 {
                    let dst = write_pos[lp];
                    let local_dst = dst - part_offsets[lp];
                    unsafe {
                        std::ptr::copy_nonoverlapping(
                            buf_k[lp].as_ptr(),
                            key_ptrs[lp].0.add(local_dst),
                            rem,
                        );
                        std::ptr::copy_nonoverlapping(buf_p[lp].as_ptr(), pptr.0.add(dst), rem);
                    }
                }
            }
        },
    );

    // Safety: Pass 1b wrote every element in every partition and in `0..total_kmers`.
    for (pk, &n) in part_keys.iter_mut().zip(part_totals.iter()) {
        unsafe { pk.set_len(n) };
    }
    unsafe { pos_strand.set_len(total_kmers) };

    // Pass 2: Rayon-parallel fine-bucketing, single sort by `(h_rot, ps)`, exact 64-bit
    // per-strand freq_limit filtering, direct 32-bit tag compaction, and immediate drop of `pk`.
    let off_len = n_parts * (sub_buckets + 1);
    let mut offsets: Vec<u32> = Vec::with_capacity(off_len);
    let mut tags: Vec<u32> = Vec::with_capacity(total_kmers + 8);
    let part_p_slices = split_by_offsets_mut(&mut pos_strand, &part_offsets);
    let part_tag_slices =
        split_by_offsets_mut(&mut tags.spare_capacity_mut()[..total_kmers], &part_offsets);
    let part_off_slices: Vec<&mut [MaybeUninit<u32>]> = offsets.spare_capacity_mut()[..off_len]
        .chunks_exact_mut(sub_buckets + 1)
        .collect();

    let kept_entries: usize = (0..n_parts)
        .into_par_iter()
        .zip(part_keys.into_par_iter())
        .zip(part_p_slices.into_par_iter())
        .zip(part_tag_slices.into_par_iter())
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
            |(sub_counts, sub_write, scratch_k, scratch_p, sort_buf_k),
             ((((lp, pk), pp), ptag), poff)| {
                let part_off = part_offsets[lp];
                let n = pk.len();
                if n == 0 {
                    for slot in poff.iter_mut() {
                        slot.write(part_off as u32);
                    }
                    drop(pk);
                    return 0usize;
                }

                sub_counts.fill(0);
                for &h_rot in pk.iter() {
                    let sub = (((h_rot as u32) >> sub_shift_hi) as usize) & sub_mask;
                    unsafe {
                        *sub_counts.get_unchecked_mut(sub) += 1;
                    }
                }

                let mut acc = 0u32;
                for s in 0..sub_buckets {
                    unsafe {
                        *sub_write.get_unchecked_mut(s) = acc;
                        acc += *sub_counts.get_unchecked(s);
                    }
                }

                if scratch_k.len() < n {
                    let new_len = n + (n >> 3) + 1024;
                    scratch_k.resize(new_len, 0);
                    scratch_p.resize(new_len, P::default());
                }
                for i in 0..n {
                    unsafe {
                        let h_rot = *pk.get_unchecked(i);
                        let sub = (((h_rot as u32) >> sub_shift_hi) as usize) & sub_mask;
                        let dst = *sub_write.get_unchecked(sub) as usize;
                        *scratch_k.get_unchecked_mut(dst) = h_rot;
                        *scratch_p.get_unchecked_mut(dst) = *pp.get_unchecked(i);
                        *sub_write.get_unchecked_mut(sub) = (dst + 1) as u32;
                    }
                }

                // Release this partition's 64-bit hashes immediately before writing 32-bit tags.
                drop(pk);

                let tag_ptr = ptag.as_mut_ptr() as *mut u32;
                let mut write_idx = 0usize;
                let mut sub_start = 0usize;
                for (slot, &cnt_u32) in poff.iter_mut().zip(sub_counts.iter()) {
                    slot.write((part_off + write_idx) as u32);
                    let cnt = cnt_u32 as usize;
                    if cnt == 0 {
                        continue;
                    }
                    let sub_end = sub_start + cnt;
                    if cnt == 1 {
                        if freq_limit >= 1 {
                            unsafe {
                                *tag_ptr.add(write_idx) =
                                    (*scratch_k.get_unchecked(sub_start) >> 32) as u32;
                                *pp.get_unchecked_mut(write_idx) =
                                    *scratch_p.get_unchecked(sub_start);
                            }
                            write_idx += 1;
                        }
                        sub_start = sub_end;
                        continue;
                    }

                    // Single sort by `(h_rot, pos_strand)` simultaneously groups identical 64-bit
                    // keys contiguously AND orders kept entries by ascending 32-bit `tag`.
                    if cnt <= 24 {
                        insertion_sort_pair(
                            &mut scratch_k[sub_start..sub_end],
                            &mut scratch_p[sub_start..sub_end],
                        );
                    } else {
                        sort_buf_k.clear();
                        sort_buf_k.reserve(cnt);
                        for idx in sub_start..sub_end {
                            unsafe {
                                sort_buf_k.push((
                                    *scratch_k.get_unchecked(idx),
                                    *scratch_p.get_unchecked(idx),
                                ));
                            }
                        }
                        sort_buf_k.sort_unstable();
                        for (off, &(h_rot, ps)) in sort_buf_k.iter().enumerate() {
                            unsafe {
                                *scratch_k.get_unchecked_mut(sub_start + off) = h_rot;
                                *scratch_p.get_unchecked_mut(sub_start + off) = ps;
                            }
                        }
                    }

                    if (cnt as u64) <= freq_limit {
                        for idx in sub_start..sub_end {
                            unsafe {
                                *tag_ptr.add(write_idx) =
                                    (*scratch_k.get_unchecked(idx) >> 32) as u32;
                                *pp.get_unchecked_mut(write_idx) = *scratch_p.get_unchecked(idx);
                            }
                            write_idx += 1;
                        }
                    } else {
                        // Exact 64-bit per-strand freq_limit filtering across runs of identical `h_rot`.
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
                            while r_end < sub_end
                                && unsafe { *scratch_k.get_unchecked(r_end) } == h_rot
                            {
                                if unsafe { scratch_p.get_unchecked(r_end).strand_bit() } == 0 {
                                    cnt0 += 1;
                                } else {
                                    cnt1 += 1;
                                }
                                r_end += 1;
                            }
                            let keep0 = cnt0 <= freq_limit;
                            let keep1 = cnt1 <= freq_limit;
                            if keep0 || keep1 {
                                let tag = (h_rot >> 32) as u32;
                                if keep0 && keep1 {
                                    for &ps in &scratch_p[r..r_end] {
                                        unsafe {
                                            *tag_ptr.add(write_idx) = tag;
                                            *pp.get_unchecked_mut(write_idx) = ps;
                                        }
                                        write_idx += 1;
                                    }
                                } else {
                                    for &ps in &scratch_p[r..r_end] {
                                        let is_keep =
                                            if ps.strand_bit() == 0 { keep0 } else { keep1 };
                                        if is_keep {
                                            unsafe {
                                                *tag_ptr.add(write_idx) = tag;
                                                *pp.get_unchecked_mut(write_idx) = ps;
                                            }
                                            write_idx += 1;
                                        }
                                    }
                                }
                            }
                            r = r_end;
                        }
                    }

                    sub_start = sub_end;
                }
                poff[sub_buckets].write((part_off + write_idx) as u32);
                // Zero-initialize any slack entries at the tail of this partition so `tags` is
                // 100% initialized across `[0 .. total_kmers + 8)` for wide AVX2 loads.
                for slot in &mut ptag[write_idx..n] {
                    slot.write(0);
                }
                write_idx
            },
        )
        .sum();

    // Initialize the 8 trailing SIMD padding slots past `total_kmers`.
    let tail_spare = &mut tags.spare_capacity_mut()[total_kmers..total_kmers + 8];
    for slot in tail_spare.iter_mut() {
        slot.write(0);
    }

    // Safety: every partition wrote all `sub_buckets + 1` elements of its `poff` slice and all
    // `n` slots of its `ptag` slice, plus the 8 trailing padding slots were initialized above.
    unsafe {
        offsets.set_len(off_len);
        tags.set_len(total_kmers + 8);
    }

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
        },
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
}
