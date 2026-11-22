use std::cmp::Reverse;

use rayon::prelude::*;

use crate::error::Result;
use crate::extend::{
    CoverageMatrix, ExtendOptions, SkipScratch, exact_rejected_skip, extend_single_hit_outcome,
    has_avx2, has_bmi2, has_popcnt,
};
#[cfg(target_arch = "x86_64")]
use crate::extend::{extend_single_hit_outcome_avx2_bmi2, extend_single_hit_outcome_avx2_swar};
use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::pipeline::{CompareConfig, CompareResult, CompareStats};
use crate::ref_index::{
    PosStrandVec, RefIndex, RefIndexKind, fmix64, prefetch_t2, probe_bucket_binary,
    probe_sorted_keys_scalar,
};
use crate::seq::Genome;

const QUERY_CHUNK_BYTES: usize = 1 << 20; // ~1 Mb dynamic Rayon chunks
const BATCH_SIZE: usize = 64;

const HSP_CACHE_BITS: u32 = 17;
const HSP_CACHE_LEN: usize = 1 << HSP_CACHE_BITS; // 131,072 slots = 2 MB
const NOISE_CACHE_BITS: u32 = 16;
const NOISE_CACHE_LEN: usize = 1 << NOISE_CACHE_BITS; // 65,536 slots = 1 MB
const FIB_MULT: u64 = 0x9E37_79B9_7F4A_7C15;

#[derive(Clone, Copy)]
#[repr(C)]
struct DiagSlot {
    diag: i64,
    max_x: u64,
}

const EMPTY_SLOT: DiagSlot = DiagSlot {
    diag: i64::MIN,
    max_x: 0,
};

/// Lock-free 2-tier per-thread, per-strand Fibonacci-hashed diagonal cache.
///
/// Separates accepted HSPs / exact reject-skip bounds (`hsp`, 2 MB) from short rejected hits
/// (`noise`, 1 MB) so random noise never evicts active HSP diagonals.
struct StrandDiagCache {
    hsp: Vec<DiagSlot>,
    noise: Vec<DiagSlot>,
}

impl StrandDiagCache {
    fn new() -> Self {
        Self {
            hsp: vec![EMPTY_SLOT; HSP_CACHE_LEN],
            noise: vec![EMPTY_SLOT; NOISE_CACHE_LEN],
        }
    }

    fn clear(&mut self) {
        self.hsp.fill(EMPTY_SLOT);
        self.noise.fill(EMPTY_SLOT);
    }

    #[inline(always)]
    fn is_covered(&self, diag: i64, pos_x: u64) -> bool {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let hi = (h >> (64 - HSP_CACHE_BITS)) as usize;
        // Safety: `hi < HSP_CACHE_LEN` by right-shift of 64-bit value.
        let hslot = unsafe { *self.hsp.get_unchecked(hi) };
        if hslot.diag == diag && pos_x <= hslot.max_x {
            return true;
        }
        let ni = (h >> (64 - NOISE_CACHE_BITS)) as usize;
        // Safety: `ni < NOISE_CACHE_LEN` by right-shift of 64-bit value.
        let nslot = unsafe { *self.noise.get_unchecked(ni) };
        nslot.diag == diag && pos_x <= nslot.max_x
    }

    #[inline(always)]
    fn update_hsp(&mut self, diag: i64, max_x: u64) {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let hi = (h >> (64 - HSP_CACHE_BITS)) as usize;
        // Safety: `hi < HSP_CACHE_LEN`.
        let slot = unsafe { self.hsp.get_unchecked_mut(hi) };
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
        let ni = (h >> (64 - NOISE_CACHE_BITS)) as usize;
        // Safety: `ni < NOISE_CACHE_LEN`.
        let slot = unsafe { self.noise.get_unchecked_mut(ni) };
        if slot.diag == diag {
            if max_x > slot.max_x {
                slot.max_x = max_x;
            }
        } else {
            *slot = DiagSlot { diag, max_x };
        }
    }
}

/// Typed `pos_strand` slice abstraction so `flush_hash` / `flush_direct` monomorphize once over
/// `u32` vs `u64` instead of matching on [`PosStrandVec`] per raw hit.
trait PsRead: Copy {
    unsafe fn unpack_at(self, idx: usize) -> (u64, u8);
    unsafe fn prefetch_at(self, idx: usize);
}

impl PsRead for *const u32 {
    #[inline(always)]
    unsafe fn unpack_at(self, idx: usize) -> (u64, u8) {
        let raw = unsafe { *self.add(idx) };
        ((raw >> 1) as u64, (raw & 1) as u8)
    }

    #[inline(always)]
    unsafe fn prefetch_at(self, idx: usize) {
        prefetch_t2(unsafe { self.add(idx) });
    }
}

impl PsRead for *const u64 {
    #[inline(always)]
    unsafe fn unpack_at(self, idx: usize) -> (u64, u8) {
        let raw = unsafe { *self.add(idx) };
        (raw >> 1, (raw & 1) as u8)
    }

    #[inline(always)]
    unsafe fn prefetch_at(self, idx: usize) {
        prefetch_t2(unsafe { self.add(idx) });
    }
}

/// 32-byte batch item storing both k-mer query metadata and Stage-2 bucket `[start, end)` range
/// in place (eliminating any separate stack range buffer or per-batch `memset`).
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct QueryItem {
    canon: u64,
    pos_x: u64,
    b_idx: u32,
    start: u32,
    end: u32,
    /// `0` = forward canonical (`kf < kr`), `1` = revcomp canonical (`kr < kf`), `2` = palindrome (`kf == kr`).
    q_kind: u8,
}

struct WorkerState {
    fwd_cache: StrandDiagCache,
    rev_cache: Option<StrandDiagCache>,
    skip_scratch: SkipScratch,
    last_chunk_end: usize,
    sp_x_cursor: usize,
    fragments: Vec<Fragment>,
    coverage: CoverageMatrix,
    word_matches: u64,
    raw_hits: u64,
    hits_used: u64,
}

impl WorkerState {
    fn new(compare_reverse: bool) -> Self {
        Self {
            fwd_cache: StrandDiagCache::new(),
            rev_cache: if compare_reverse {
                Some(StrandDiagCache::new())
            } else {
                None
            },
            skip_scratch: SkipScratch::default(),
            last_chunk_end: 0,
            sp_x_cursor: 0,
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            word_matches: 0,
            raw_hits: 0,
            hits_used: 0,
        }
    }

    /// Advance the monotonic `query.seq.special_runs` cursor so `sp_x[..cursor].end <= pos_x`.
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

/// Runtime ISA tier selected once per comparison (dispatched once per ~1 Mb chunk so bucket
/// probing, BMI2 `pext`, hardware `popcnt`, and 2-bit AVX2 extension inline with zero `vzeroupper`
/// boundaries on the hot path).
#[derive(Clone, Copy)]
enum IsaTier {
    #[cfg(target_arch = "x86_64")]
    Avx2Bmi2Popcnt,
    #[cfg(target_arch = "x86_64")]
    Avx2Popcnt,
    Portable,
}

/// Execute the streaming canonical-index comparison over `query` against a pre-built [`RefIndex`].
///
/// Zero hit materialization or global hit sorting is performed: query k-mers are streamed
/// directly from 2-bit packed `u64` words in ~1 Mb parallel chunks using 64-item software-prefetched
/// batches, filtered through lock-free per-thread 2-tier diagonal caches, and extended immediately
/// with bit-exact 2-bit + AVX2/BMI2 SIMD.
pub fn compare_streaming(
    query: &Genome,
    reference: &Genome,
    ref_index: &RefIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
) -> Result<CompareResult> {
    config.validate()?;

    let wl = config.word_length;
    if query.seq.len < wl || ref_index.kept_entries == 0 {
        return Ok(CompareResult {
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            stats: CompareStats::default(),
        });
    }

    let compare_reverse = config.compare_reverse && rev_ref.is_some();
    let isa = {
        #[cfg(target_arch = "x86_64")]
        {
            if has_avx2() && has_popcnt() {
                if has_bmi2() {
                    IsaTier::Avx2Bmi2Popcnt
                } else {
                    IsaTier::Avx2Popcnt
                }
            } else {
                IsaTier::Portable
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (has_avx2(), has_bmi2(), has_popcnt());
            IsaTier::Portable
        }
    };

    let mask = if wl == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * wl)) - 1
    };
    let rc_shift = (2 * (wl - 1)) as u32;
    let ny = reference.seq.len as u64;

    let fwd_opts = ExtendOptions {
        min_length: config.min_length,
        min_similarity: config.min_similarity,
        word_length: wl,
        fixed_length: config.fixed_length,
        strand: 'f',
        scoring: config.scoring,
    };
    let rev_opts = ExtendOptions {
        strand: 'r',
        ..fwd_opts
    };

    // Build ~1 Mb dynamic Rayon chunks partitioning all valid start positions `[0 .. max_starts)`.
    let max_starts = query.seq.len - wl + 1;
    let mut chunks = Vec::with_capacity(max_starts.div_ceil(QUERY_CHUNK_BYTES));
    let mut s = 0usize;
    while s < max_starts {
        let e = (s + QUERY_CHUNK_BYTES).min(max_starts);
        chunks.push((s, e));
        s = e;
    }

    // Process chunks in parallel, keeping one reusable `WorkerState` per Rayon thread.
    let worker_results: Vec<WorkerState> = chunks
        .into_par_iter()
        .fold(
            || WorkerState::new(compare_reverse),
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

                scan_query_chunk(
                    query,
                    reference,
                    rev_ref,
                    ref_index,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    config.soft_mask,
                    config.max_rep.max(1),
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

    // Merge thread-local outputs.
    let total_frags: usize = worker_results.iter().map(|w| w.fragments.len()).sum();
    let mut fragments = Vec::with_capacity(total_frags);
    let mut coverage = CoverageMatrix::new();
    let mut stats = CompareStats::default();

    for mut w in worker_results {
        stats.word_matches += w.word_matches;
        stats.raw_hits += w.raw_hits;
        stats.hits_used += w.hits_used;
        coverage.merge_from(&w.coverage);
        fragments.append(&mut w.fragments);
    }
    stats.filtered_hits = stats.hits_used;

    // Deduplicate chunk-boundary crossings and overlapping extensions on identical diagonals.
    dedup_fragments(&mut fragments);

    for f in &fragments {
        if f.strand == 'f' {
            stats.forward_fragments += 1;
        } else {
            stats.reverse_fragments += 1;
        }
    }

    Ok(CompareResult {
        fragments,
        coverage,
        stats,
    })
}

/// Inline AVX2 sorted-bucket probe expanded directly inside `flush_hash`'s Stage 3 loop so LLVM
/// never inserts a function call or `vzeroupper` boundary between bucket probing and diagonal-cache
/// filtering.
#[cfg(target_arch = "x86_64")]
macro_rules! probe_bucket_avx2_inline {
    ($keys:expr, $keys_ptr:expr, $start:expr, $end:expr, $key:expr) => {{
        use core::arch::x86_64::*;

        let start = $start;
        let end = $end;
        let key = $key;
        let len = end - start;
        'probe: {
            if len <= 16 {
                let ptr = $keys_ptr;
                let mut j = start;
                if len >= 4 {
                    let vkey = _mm256_set1_epi64x(key as i64);
                    while j + 4 <= end {
                        let vcand = unsafe { _mm256_loadu_si256(ptr.add(j) as *const __m256i) };
                        let cmp = _mm256_cmpeq_epi64(vcand, vkey);
                        let mask = _mm256_movemask_pd(_mm256_castsi256_pd(cmp)) as u32;
                        if mask != 0 {
                            let lo = j + mask.trailing_zeros() as usize;
                            if (mask & 8) == 0 {
                                let hi = j + (32 - mask.leading_zeros()) as usize;
                                break 'probe (lo, hi);
                            }
                            let mut hi = j + 4;
                            while hi < end && unsafe { *ptr.add(hi) } == key {
                                hi += 1;
                            }
                            break 'probe (lo, hi);
                        }
                        if unsafe { *ptr.add(j + 3) } > key {
                            break 'probe (start, start);
                        }
                        j += 4;
                    }
                }
                while j < end {
                    let k = unsafe { *ptr.add(j) };
                    if k >= key {
                        if k > key {
                            break 'probe (start, start);
                        }
                        let lo = j;
                        let mut hi = j + 1;
                        while hi < end && unsafe { *ptr.add(hi) } == key {
                            hi += 1;
                        }
                        break 'probe (lo, hi);
                    }
                    j += 1;
                }
                break 'probe (start, start);
            }
            probe_bucket_binary($keys, start, end, key)
        }
    }};
}

macro_rules! probe_bucket_scalar_inline {
    ($keys:expr, $_keys_ptr:expr, $start:expr, $end:expr, $key:expr) => {
        probe_sorted_keys_scalar($keys, $start, $end, $key)
    };
}

/// Generate a complete, single-ISA-domain chunk scanner + 3-stage prefetched batch flusher.
///
/// `flush_hash` / `flush_direct` are kept `#[inline(never)]` and monomorphized over `PsRead`
/// (`*const u32` vs `*const u64`) so:
/// 1. The k-mer generator loop stays compact and register-light.
/// 2. Bucket probing (`vpcmpeqq` with early-exit) and the 2-tier diagonal-cache check
///    (`is_covered`, which filters ~94% of raw hits) execute inline inside Stage 3 with ZERO
///    function calls and ZERO `vzeroupper` transitions on the hot path.
/// 3. Only the ~6% diagonal-cache misses call out to `miss_fwd` / `miss_rev`.
macro_rules! impl_stream_tier {
    (
        $scan_hash:ident,
        $scan_direct:ident,
        $flush_hash:ident,
        $flush_direct:ident,
        $miss_fwd:ident,
        $miss_rev:ident,
        [$($attr:meta),*],
        $probe_macro:ident,
        $extend_call:expr
    ) => {
        $(#[$attr])*
        #[allow(clippy::too_many_arguments)]
        unsafe fn $scan_hash<const WL32: bool, P: PsRead>(
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            keys: &[u64],
            ps_ptr: P,
            sub_bits: u32,
            hash_shift: u32,
            chunk_start: usize,
            chunk_end: usize,
            wl: usize,
            ny: u64,
            mask: u64,
            rc_shift: u32,
            soft_mask: bool,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState,
        ) {
            let seq_x = &query.seq;
            let scan_end = (chunk_end + wl - 1).min(seq_x.len);
            let mut batch = [QueryItem::default(); BATCH_SIZE];
            let mut batch_len = 0usize;

            seq_x.for_each_valid_run(chunk_start, scan_end, soft_mask, |rs, re| {
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
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    for _ in 0..avail {
                        let b = word & 3;
                        word >>= 2;
                        kf = if WL32 { (kf << 2) | b } else { ((kf << 2) | b) & mask };
                        kr = (kr >> 2) | ((b ^ 3) << if WL32 { 62 } else { rc_shift });
                    }
                    idx += avail;
                }

                while idx < re {
                    let w_idx = idx >> 5;
                    let lane = idx & 31;
                    let avail = (32 - lane).min(re - idx);
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    let start_pos = (idx + 1 - wl) as u64;
                    for pos_x in start_pos..start_pos + avail as u64 {
                        let b = word & 3;
                        word >>= 2;
                        kf = if WL32 { (kf << 2) | b } else { ((kf << 2) | b) & mask };
                        kr = (kr >> 2) | ((b ^ 3) << if WL32 { 62 } else { rc_shift });
                        let (canon, q_kind) = if kf < kr {
                            (kf, 0u8)
                        } else if kr < kf {
                            (kr, 1u8)
                        } else {
                            (kf, 2u8)
                        };
                        let b_bucket = (fmix64(canon) >> hash_shift) as u32;
                        let b_idx = b_bucket + (b_bucket >> sub_bits);
                        unsafe {
                            *batch.get_unchecked_mut(batch_len) = QueryItem {
                                canon,
                                pos_x,
                                b_idx,
                                start: 0,
                                end: 0,
                                q_kind,
                            };
                        }
                        batch_len += 1;
                        if batch_len == BATCH_SIZE {
                            unsafe {
                                $flush_hash(
                                    &mut batch[..BATCH_SIZE],
                                    query,
                                    reference,
                                    rev_ref,
                                    offsets,
                                    keys,
                                    ps_ptr,
                                    wl as u64,
                                    ny,
                                    max_rep,
                                    fwd_opts,
                                    rev_opts,
                                    state,
                                );
                            }
                            batch_len = 0;
                        }
                    }
                    idx += avail;
                }
            });

            if batch_len > 0 {
                unsafe {
                    $flush_hash(
                        &mut batch[..batch_len],
                        query,
                        reference,
                        rev_ref,
                        offsets,
                        keys,
                        ps_ptr,
                        wl as u64,
                        ny,
                        max_rep,
                        fwd_opts,
                        rev_opts,
                        state,
                    );
                }
            }
        }

        $(#[$attr])*
        #[allow(clippy::too_many_arguments)]
        unsafe fn $scan_direct<P: PsRead>(
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            ps_ptr: P,
            sub_bits: u32,
            chunk_start: usize,
            chunk_end: usize,
            wl: usize,
            ny: u64,
            mask: u64,
            rc_shift: u32,
            soft_mask: bool,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState,
        ) {
            let seq_x = &query.seq;
            let scan_end = (chunk_end + wl - 1).min(seq_x.len);
            let mut batch = [QueryItem::default(); BATCH_SIZE];
            let mut batch_len = 0usize;

            seq_x.for_each_valid_run(chunk_start, scan_end, soft_mask, |rs, re| {
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
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    for _ in 0..avail {
                        let b = word & 3;
                        word >>= 2;
                        kf = ((kf << 2) | b) & mask;
                        kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                    }
                    idx += avail;
                }

                while idx < re {
                    let w_idx = idx >> 5;
                    let lane = idx & 31;
                    let avail = (32 - lane).min(re - idx);
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    let start_pos = (idx + 1 - wl) as u64;
                    for pos_x in start_pos..start_pos + avail as u64 {
                        let b = word & 3;
                        word >>= 2;
                        kf = ((kf << 2) | b) & mask;
                        kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                        let (canon, q_kind) = if kf < kr {
                            (kf, 0u8)
                        } else if kr < kf {
                            (kr, 1u8)
                        } else {
                            (kf, 2u8)
                        };
                        let b_bucket = canon as u32;
                        let b_idx = b_bucket + (b_bucket >> sub_bits);
                        unsafe {
                            *batch.get_unchecked_mut(batch_len) = QueryItem {
                                canon,
                                pos_x,
                                b_idx,
                                start: 0,
                                end: 0,
                                q_kind,
                            };
                        }
                        batch_len += 1;
                        if batch_len == BATCH_SIZE {
                            unsafe {
                                $flush_direct(
                                    &mut batch[..BATCH_SIZE],
                                    query,
                                    reference,
                                    rev_ref,
                                    offsets,
                                    ps_ptr,
                                    wl as u64,
                                    ny,
                                    max_rep,
                                    fwd_opts,
                                    rev_opts,
                                    state,
                                );
                            }
                            batch_len = 0;
                        }
                    }
                    idx += avail;
                }
            });

            if batch_len > 0 {
                unsafe {
                    $flush_direct(
                        &mut batch[..batch_len],
                        query,
                        reference,
                        rev_ref,
                        offsets,
                        ps_ptr,
                        wl as u64,
                        ny,
                        max_rep,
                        fwd_opts,
                        rev_opts,
                        state,
                    );
                }
            }
        }

        $(#[$attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $flush_hash<P: PsRead>(
            batch: &mut [QueryItem],
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            keys: &[u64],
            ps_ptr: P,
            wl: u64,
            ny: u64,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState,
        ) {
            let off_ptr = offsets.as_ptr();
            let keys_ptr = keys.as_ptr();

            // Stage 1: L3/L2 prefetch `offsets[b_idx]` for all items in the batch.
            for item in batch.iter() {
                prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
            }

            // Stage 2: load `[start, end)` in place and prefetch `keys[start]`.
            for item in batch.iter_mut() {
                let b_idx = item.b_idx as usize;
                let s = unsafe { *off_ptr.add(b_idx) };
                let e = unsafe { *off_ptr.add(b_idx + 1) };
                item.start = s;
                item.end = e;
                if s < e {
                    prefetch_t2(unsafe { keys_ptr.add(s as usize) });
                }
            }

            // Stage 3: probe buckets inline, filter ~94% of raw hits via inline diagonal cache,
            // and dispatch only the ~6% uncovered hits to the out-of-line extender.
            for item in batch.iter() {
                let s = item.start as usize;
                let e = item.end as usize;
                if s == e {
                    continue;
                }
                let (lo, hi) = $probe_macro!(keys, keys_ptr, s, e, item.canon);
                if lo == hi {
                    continue;
                }

                let total = hi - lo;
                let step = if total > max_rep {
                    (total / max_rep).max(1)
                } else {
                    1
                };
                let pos_x = item.pos_x;
                let q_kind = item.q_kind;
                let mut matched_any = false;
                let mut idx = lo;

                while idx < hi {
                    let (pos_y, r_strand) = unsafe { ps_ptr.unpack_at(idx) };
                    idx += step;

                    if q_kind == 2 {
                        matched_any = true;
                        state.raw_hits += 1;
                        let diag_f = pos_x as i64 - pos_y as i64;
                        if !state.fwd_cache.is_covered(diag_f, pos_x) {
                            unsafe {
                                $miss_fwd(diag_f, pos_x, pos_y, query, reference, fwd_opts, state);
                            }
                        }
                        if let Some(rr) = rev_ref {
                            state.raw_hits += 1;
                            let rev_pos_y = ny - pos_y - wl;
                            let diag_r = pos_x as i64 - rev_pos_y as i64;
                            if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                                unsafe {
                                    $miss_rev(diag_r, pos_x, rev_pos_y, query, rr, rev_opts, state);
                                }
                            }
                        }
                    } else if r_strand == q_kind {
                        matched_any = true;
                        state.raw_hits += 1;
                        let diag_f = pos_x as i64 - pos_y as i64;
                        if !state.fwd_cache.is_covered(diag_f, pos_x) {
                            unsafe {
                                $miss_fwd(diag_f, pos_x, pos_y, query, reference, fwd_opts, state);
                            }
                        }
                    } else if let Some(rr) = rev_ref {
                        matched_any = true;
                        state.raw_hits += 1;
                        let rev_pos_y = ny - pos_y - wl;
                        let diag_r = pos_x as i64 - rev_pos_y as i64;
                        if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                            unsafe {
                                $miss_rev(diag_r, pos_x, rev_pos_y, query, rr, rev_opts, state);
                            }
                        }
                    }
                }

                if matched_any {
                    state.word_matches += 1;
                }
            }
        }

        $(#[$attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $flush_direct<P: PsRead>(
            batch: &mut [QueryItem],
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            ps_ptr: P,
            wl: u64,
            ny: u64,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState,
        ) {
            let off_ptr = offsets.as_ptr();

            for item in batch.iter() {
                prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
            }

            for item in batch.iter_mut() {
                let b_idx = item.b_idx as usize;
                let s = unsafe { *off_ptr.add(b_idx) };
                let e = unsafe { *off_ptr.add(b_idx + 1) };
                item.start = s;
                item.end = e;
                if s < e {
                    unsafe { ps_ptr.prefetch_at(s as usize) };
                }
            }

            for item in batch.iter() {
                let lo = item.start as usize;
                let hi = item.end as usize;
                if lo == hi {
                    continue;
                }
                let total = hi - lo;
                let step = if total > max_rep {
                    (total / max_rep).max(1)
                } else {
                    1
                };
                let pos_x = item.pos_x;
                let q_kind = item.q_kind;
                let mut matched_any = false;
                let mut idx = lo;

                while idx < hi {
                    let (pos_y, r_strand) = unsafe { ps_ptr.unpack_at(idx) };
                    idx += step;

                    if q_kind == 2 {
                        matched_any = true;
                        state.raw_hits += 1;
                        let diag_f = pos_x as i64 - pos_y as i64;
                        if !state.fwd_cache.is_covered(diag_f, pos_x) {
                            unsafe {
                                $miss_fwd(diag_f, pos_x, pos_y, query, reference, fwd_opts, state);
                            }
                        }
                        if let Some(rr) = rev_ref {
                            state.raw_hits += 1;
                            let rev_pos_y = ny - pos_y - wl;
                            let diag_r = pos_x as i64 - rev_pos_y as i64;
                            if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                                unsafe {
                                    $miss_rev(diag_r, pos_x, rev_pos_y, query, rr, rev_opts, state);
                                }
                            }
                        }
                    } else if r_strand == q_kind {
                        matched_any = true;
                        state.raw_hits += 1;
                        let diag_f = pos_x as i64 - pos_y as i64;
                        if !state.fwd_cache.is_covered(diag_f, pos_x) {
                            unsafe {
                                $miss_fwd(diag_f, pos_x, pos_y, query, reference, fwd_opts, state);
                            }
                        }
                    } else if let Some(rr) = rev_ref {
                        matched_any = true;
                        state.raw_hits += 1;
                        let rev_pos_y = ny - pos_y - wl;
                        let diag_r = pos_x as i64 - rev_pos_y as i64;
                        if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                            unsafe {
                                $miss_rev(diag_r, pos_x, rev_pos_y, query, rr, rev_opts, state);
                            }
                        }
                    }
                }

                if matched_any {
                    state.word_matches += 1;
                }
            }
        }

        $(#[$attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $miss_fwd(
            diag: i64,
            pos_x: u64,
            pos_y: u64,
            query: &Genome,
            reference: &Genome,
            opts: ExtendOptions,
            state: &mut WorkerState,
        ) {
            state.hits_used += 1;
            let cx_x = state.sp_x_cursor_for(&query.seq.special_runs, pos_x);
            let outcome = ($extend_call)(
                Hit { diag, pos_x },
                query,
                reference,
                opts,
                cx_x,
                &mut state.coverage,
            );
            if let Some(frag) = outcome.fragment {
                state.fwd_cache.update_hsp(diag, outcome.xfil_max);
                state.fragments.push(frag);
            } else if outcome.xfil_max > pos_x + 256 {
                let skip_end = exact_rejected_skip(
                    &query.seq,
                    &reference.seq,
                    pos_x,
                    pos_y,
                    outcome.xfil_max,
                    opts.word_length,
                    opts.scoring,
                    &mut state.skip_scratch,
                );
                state.fwd_cache.update_hsp(diag, skip_end);
            } else {
                state
                    .fwd_cache
                    .update_noise(diag, pos_x + 2 * opts.word_length as u64 - 1);
            }
        }

        $(#[$attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $miss_rev(
            diag: i64,
            pos_x: u64,
            pos_y: u64,
            query: &Genome,
            rev_ref: &Genome,
            opts: ExtendOptions,
            state: &mut WorkerState,
        ) {
            state.hits_used += 1;
            let cx_x = state.sp_x_cursor_for(&query.seq.special_runs, pos_x);
            let outcome = ($extend_call)(
                Hit { diag, pos_x },
                query,
                rev_ref,
                opts,
                cx_x,
                &mut state.coverage,
            );
            let rev_cache = state.rev_cache.as_mut().unwrap();
            if let Some(frag) = outcome.fragment {
                rev_cache.update_hsp(diag, outcome.xfil_max);
                state.fragments.push(frag);
            } else if outcome.xfil_max > pos_x + 256 {
                let skip_end = exact_rejected_skip(
                    &query.seq,
                    &rev_ref.seq,
                    pos_x,
                    pos_y,
                    outcome.xfil_max,
                    opts.word_length,
                    opts.scoring,
                    &mut state.skip_scratch,
                );
                rev_cache.update_hsp(diag, skip_end);
            } else {
                rev_cache.update_noise(diag, pos_x + 2 * opts.word_length as u64 - 1);
            }
        }
    };
}

#[cfg(target_arch = "x86_64")]
impl_stream_tier!(
    scan_hash_avx2_bmi2,
    scan_direct_avx2_bmi2,
    flush_hash_avx2_bmi2,
    flush_direct_avx2_bmi2,
    miss_fwd_avx2_bmi2,
    miss_rev_avx2_bmi2,
    [target_feature(enable = "avx2,bmi2,popcnt")],
    probe_bucket_avx2_inline,
    |h, gx, gy, opts, cx_x, cov| unsafe {
        extend_single_hit_outcome_avx2_bmi2(h, gx, gy, opts, cx_x, cov)
    }
);

#[cfg(target_arch = "x86_64")]
impl_stream_tier!(
    scan_hash_avx2_swar,
    scan_direct_avx2_swar,
    flush_hash_avx2_swar,
    flush_direct_avx2_swar,
    miss_fwd_avx2_swar,
    miss_rev_avx2_swar,
    [target_feature(enable = "avx2,popcnt")],
    probe_bucket_avx2_inline,
    |h, gx, gy, opts, cx_x, cov| unsafe {
        extend_single_hit_outcome_avx2_swar(h, gx, gy, opts, cx_x, cov)
    }
);

impl_stream_tier!(
    scan_hash_portable,
    scan_direct_portable,
    flush_hash_portable,
    flush_direct_portable,
    miss_fwd_portable,
    miss_rev_portable,
    [],
    probe_bucket_scalar_inline,
    |h, gx, gy, opts, _cx_x, cov| extend_single_hit_outcome(h, gx, gy, opts, false, false, cov)
);

#[allow(clippy::too_many_arguments)]
fn scan_query_chunk(
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    ref_index: &RefIndex,
    chunk_start: usize,
    chunk_end: usize,
    wl: usize,
    ny: u64,
    mask: u64,
    rc_shift: u32,
    soft_mask: bool,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    isa: IsaTier,
    state: &mut WorkerState,
) {
    match (&ref_index.kind, isa) {
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                offsets,
                keys,
                pos_strand,
            },
            IsaTier::Avx2Bmi2Popcnt,
        ) => match (wl == 32, pos_strand) {
            (true, PosStrandVec::U32(v)) => unsafe {
                scan_hash_avx2_bmi2::<true, *const u32>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    32,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (false, PosStrandVec::U32(v)) => unsafe {
                scan_hash_avx2_bmi2::<false, *const u32>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (true, PosStrandVec::U64(v)) => unsafe {
                scan_hash_avx2_bmi2::<true, *const u64>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    32,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (false, PosStrandVec::U64(v)) => unsafe {
                scan_hash_avx2_bmi2::<false, *const u64>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
        },
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                offsets,
                keys,
                pos_strand,
            },
            IsaTier::Avx2Popcnt,
        ) => match (wl == 32, pos_strand) {
            (true, PosStrandVec::U32(v)) => unsafe {
                scan_hash_avx2_swar::<true, *const u32>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    32,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (false, PosStrandVec::U32(v)) => unsafe {
                scan_hash_avx2_swar::<false, *const u32>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (true, PosStrandVec::U64(v)) => unsafe {
                scan_hash_avx2_swar::<true, *const u64>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    32,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (false, PosStrandVec::U64(v)) => unsafe {
                scan_hash_avx2_swar::<false, *const u64>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
        },
        (
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                offsets,
                keys,
                pos_strand,
            },
            IsaTier::Portable,
        ) => match (wl == 32, pos_strand) {
            (true, PosStrandVec::U32(v)) => unsafe {
                scan_hash_portable::<true, *const u32>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    32,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (false, PosStrandVec::U32(v)) => unsafe {
                scan_hash_portable::<false, *const u32>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (true, PosStrandVec::U64(v)) => unsafe {
                scan_hash_portable::<true, *const u64>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    32,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            (false, PosStrandVec::U64(v)) => unsafe {
                scan_hash_portable::<false, *const u64>(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    keys,
                    v.as_ptr(),
                    *sub_bits,
                    *hash_shift,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
        },
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::DirectCsr {
                sub_bits,
                offsets,
                pos_strand,
            },
            IsaTier::Avx2Bmi2Popcnt,
        ) => match pos_strand {
            PosStrandVec::U32(v) => unsafe {
                scan_direct_avx2_bmi2(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    v.as_ptr(),
                    *sub_bits,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            PosStrandVec::U64(v) => unsafe {
                scan_direct_avx2_bmi2(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    v.as_ptr(),
                    *sub_bits,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
        },
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::DirectCsr {
                sub_bits,
                offsets,
                pos_strand,
            },
            IsaTier::Avx2Popcnt,
        ) => match pos_strand {
            PosStrandVec::U32(v) => unsafe {
                scan_direct_avx2_swar(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    v.as_ptr(),
                    *sub_bits,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            PosStrandVec::U64(v) => unsafe {
                scan_direct_avx2_swar(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    v.as_ptr(),
                    *sub_bits,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
        },
        (
            RefIndexKind::DirectCsr {
                sub_bits,
                offsets,
                pos_strand,
            },
            IsaTier::Portable,
        ) => match pos_strand {
            PosStrandVec::U32(v) => unsafe {
                scan_direct_portable(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    v.as_ptr(),
                    *sub_bits,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
            PosStrandVec::U64(v) => unsafe {
                scan_direct_portable(
                    query,
                    reference,
                    rev_ref,
                    offsets,
                    v.as_ptr(),
                    *sub_bits,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    soft_mask,
                    max_rep,
                    fwd_opts,
                    rev_opts,
                    state,
                );
            },
        },
    }
}

/// Sort merged fragments and remove any fragment fully contained within a preceding fragment
/// on the same `(strand, diag, seq_x, seq_y)`.
pub fn dedup_fragments(fragments: &mut Vec<Fragment>) {
    if fragments.len() <= 1 {
        return;
    }
    fragments.sort_unstable_by_key(|f| {
        (
            f.strand,
            f.diag,
            f.seq_x,
            f.seq_y,
            f.x_start,
            Reverse(f.x_end),
            Reverse(f.score),
        )
    });

    let mut write = 1usize;
    for read in 1..fragments.len() {
        let prev = &fragments[write - 1];
        let cur = &fragments[read];
        let contained = cur.strand == prev.strand
            && cur.diag == prev.diag
            && cur.seq_x == prev.seq_x
            && cur.seq_y == prev.seq_y
            && cur.x_end <= prev.x_end;
        if !contained {
            if write != read {
                fragments[write] = fragments[read].clone();
            }
            write += 1;
        }
    }
    fragments.truncate(write);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dedup_removes_identical_and_contained_fragments_only() {
        let mk = |strand, diag, seq_x, seq_y, x_start, x_end, score| Fragment {
            diag,
            x_start,
            y_start: x_start,
            x_end,
            y_end: x_end,
            length: x_end - x_start + 1,
            ident: x_end - x_start + 1,
            score,
            similarity: 100.0,
            seq_x,
            seq_y,
            block: 0,
            strand,
        };

        let mut frags = vec![
            mk('f', 10, 0, 0, 100, 200, 400),
            mk('f', 10, 0, 0, 100, 200, 400), // exact duplicate
            mk('f', 10, 0, 0, 120, 180, 240), // strictly contained
            mk('f', 10, 0, 0, 150, 250, 400), // overlapping but extends further -> keep
            mk('f', 10, 1, 0, 120, 180, 240), // different seq_x -> keep
            mk('r', 10, 0, 0, 100, 200, 400), // different strand -> keep
        ];

        dedup_fragments(&mut frags);
        assert_eq!(frags.len(), 4);
    }
}
