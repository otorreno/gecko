use std::cmp::Reverse;

use rayon::prelude::*;

use crate::error::Result;
use crate::extend::{
    CoverageMatrix, ExtendOptions, SkipScratch, exact_rejected_skip, extend_single_hit_outcome,
    has_avx2, has_bmi2,
};
use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::pipeline::{CompareConfig, CompareResult, CompareStats};
use crate::ref_index::{PosStrandVec, RefIndex, RefIndexKind, fmix64};
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

#[derive(Clone, Copy, Default)]
struct QueryItem {
    canon: u64,
    pos_x: u64,
    b_idx: usize,
    /// `0` = forward canonical (`kf < kr`), `1` = revcomp canonical (`kr < kf`), `2` = palindrome (`kf == kr`).
    q_kind: u8,
}

struct WorkerState {
    fwd_cache: StrandDiagCache,
    rev_cache: Option<StrandDiagCache>,
    skip_scratch: SkipScratch,
    last_chunk_end: usize,
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
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            word_matches: 0,
            raw_hits: 0,
            hits_used: 0,
        }
    }
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
    let use_avx2 = has_avx2();
    let use_bmi2 = has_bmi2();
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
                    use_avx2,
                    use_bmi2,
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
    use_avx2: bool,
    use_bmi2: bool,
    state: &mut WorkerState,
) {
    if ref_index.is_direct() {
        scan_query_chunk_impl::<true>(
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
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            use_avx2,
            use_bmi2,
            state,
        );
    } else {
        scan_query_chunk_impl::<false>(
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
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            use_avx2,
            use_bmi2,
            state,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_query_chunk_impl<const IS_DIRECT: bool>(
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
    use_avx2: bool,
    use_bmi2: bool,
    state: &mut WorkerState,
) {
    let seq_x = &query.seq;
    let scan_end = (chunk_end + wl - 1).min(seq_x.len);

    let (sub_bits, hash_shift) = match &ref_index.kind {
        RefIndexKind::DirectCsr { sub_bits, .. } => (*sub_bits, 0u32),
        RefIndexKind::HashCsr {
            sub_bits,
            hash_shift,
            ..
        } => (*sub_bits, *hash_shift),
    };

    let mut batch = [QueryItem::default(); BATCH_SIZE];
    let mut batch_len = 0usize;

    seq_x.for_each_valid_run(chunk_start, scan_end, soft_mask, |rs, re| {
        if re - rs < wl {
            return;
        }
        // Prime first `wl - 1` bases of the clean run directly from 2-bit u64 words.
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

        // Stream 32 bases per u64 load, batching 64 k-mers for software-prefetched probing.
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
                let b_idx = if IS_DIRECT {
                    let b = canon as usize;
                    b + (b >> sub_bits)
                } else {
                    let b = (fmix64(canon) >> hash_shift) as usize;
                    b + (b >> sub_bits)
                };
                batch[batch_len] = QueryItem {
                    canon,
                    pos_x,
                    b_idx,
                    q_kind,
                };
                batch_len += 1;
                if batch_len == BATCH_SIZE {
                    flush_batch(
                        &batch[..BATCH_SIZE],
                        query,
                        reference,
                        rev_ref,
                        ref_index,
                        wl as u64,
                        ny,
                        max_rep,
                        fwd_opts,
                        rev_opts,
                        use_avx2,
                        use_bmi2,
                        state,
                    );
                    batch_len = 0;
                }
            }
            idx += avail;
        }
    });

    if batch_len > 0 {
        flush_batch(
            &batch[..batch_len],
            query,
            reference,
            rev_ref,
            ref_index,
            wl as u64,
            ny,
            max_rep,
            fwd_opts,
            rev_opts,
            use_avx2,
            use_bmi2,
            state,
        );
    }
}

/// 3-stage software-prefetched batch lookup and immediate hit extension.
#[allow(clippy::too_many_arguments)]
#[inline]
fn flush_batch(
    batch: &[QueryItem],
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    ref_index: &RefIndex,
    wl: u64,
    ny: u64,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    use_avx2: bool,
    use_bmi2: bool,
    state: &mut WorkerState,
) {
    let n = batch.len();

    // Stage 1: prefetch `offsets[b_idx]` for all items in the batch.
    for item in batch {
        ref_index.prefetch_offset(item.b_idx);
    }

    // Stage 2: load `[start, end)` ranges and prefetch `keys[start]` (or `pos_strand[start]`).
    let mut ranges = [(0usize, 0usize); BATCH_SIZE];
    for i in 0..n {
        let (s, e) = ref_index.bucket_range(batch[i].b_idx);
        ranges[i] = (s, e);
        if s < e {
            ref_index.prefetch_bucket_data(s);
        }
    }

    // Stage 3: probe buckets and immediately handle all matching hits.
    for i in 0..n {
        let (s, e) = ranges[i];
        if s == e {
            continue;
        }
        let item = batch[i];
        let (lo, hi, pos_strand) = ref_index.probe_bucket(s, e, item.canon, use_avx2);
        if lo == hi {
            continue;
        }

        handle_matched_range(
            item, lo, hi, pos_strand, query, reference, rev_ref, wl, ny, max_rep, fwd_opts,
            rev_opts, use_avx2, use_bmi2, state,
        );
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn handle_matched_range(
    item: QueryItem,
    lo: usize,
    hi: usize,
    pos_strand: &PosStrandVec,
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    wl: u64,
    ny: u64,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    use_avx2: bool,
    use_bmi2: bool,
    state: &mut WorkerState,
) {
    let total = hi - lo;
    let step = if total > max_rep {
        (total / max_rep).max(1)
    } else {
        1
    };

    let mut matched_any = false;
    let mut idx = lo;
    while idx < hi {
        let (pos_y, r_strand) = pos_strand.get_unpacked(idx);
        idx += step;

        if item.q_kind == 2 {
            // Palindrome (`kf == kr`): emit both forward and reverse hits.
            matched_any = true;
            process_hit_fwd(
                item.pos_x, pos_y, query, reference, fwd_opts, use_avx2, use_bmi2, state,
            );
            if let Some(rr) = rev_ref {
                let rev_pos_y = ny - pos_y - wl;
                process_hit_rev(
                    item.pos_x, rev_pos_y, query, rr, rev_opts, use_avx2, use_bmi2, state,
                );
            }
        } else if r_strand == item.q_kind {
            matched_any = true;
            process_hit_fwd(
                item.pos_x, pos_y, query, reference, fwd_opts, use_avx2, use_bmi2, state,
            );
        } else if let Some(rr) = rev_ref {
            matched_any = true;
            let rev_pos_y = ny - pos_y - wl;
            process_hit_rev(
                item.pos_x, rev_pos_y, query, rr, rev_opts, use_avx2, use_bmi2, state,
            );
        }
    }

    if matched_any {
        state.word_matches += 1;
    }
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn process_hit_fwd(
    pos_x: u64,
    pos_y: u64,
    query: &Genome,
    reference: &Genome,
    opts: ExtendOptions,
    use_avx2: bool,
    use_bmi2: bool,
    state: &mut WorkerState,
) {
    let diag = pos_x as i64 - pos_y as i64;
    state.raw_hits += 1;
    if state.fwd_cache.is_covered(diag, pos_x) {
        return;
    }
    state.hits_used += 1;
    let outcome = extend_single_hit_outcome(
        Hit { diag, pos_x },
        query,
        reference,
        opts,
        use_avx2,
        use_bmi2,
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

#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn process_hit_rev(
    pos_x: u64,
    pos_y: u64,
    query: &Genome,
    rev_ref: &Genome,
    opts: ExtendOptions,
    use_avx2: bool,
    use_bmi2: bool,
    state: &mut WorkerState,
) {
    let diag = pos_x as i64 - pos_y as i64;
    state.raw_hits += 1;
    let rev_cache = state.rev_cache.as_mut().unwrap();
    if rev_cache.is_covered(diag, pos_x) {
        return;
    }
    state.hits_used += 1;
    let outcome = extend_single_hit_outcome(
        Hit { diag, pos_x },
        query,
        rev_ref,
        opts,
        use_avx2,
        use_bmi2,
        &mut state.coverage,
    );
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
