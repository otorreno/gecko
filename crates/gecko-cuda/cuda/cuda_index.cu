/*
 * cuda_index.cu
 *
 * Stage 1: Word-level 2-bit GPU reverse-complement generation, warp-cooperative
 * GPU run-to-bitmap materialization (W3-F12b), CTA-cooperative coalesced 2-bit
 * canonical k-mer extraction (Direct-CSR uint32_t for wl <= 12 and Hash-CSR
 * uint64_t), bit-bounded CUB DeviceRadixSort::SortPairs, exact per-strand
 * freq_limit compaction + fused word-coalesced 1-bit dropped bitmap + inverted
 * 1-bit non-unique/palindromic ref k-mer bitmap (W3-F13a), fused
 * thrust::transform_iterator tag selection, and warp-aggregated pure-GPU CSR scan.
 */

#include <utility>
#include <cub/cub.cuh>
#include <thrust/iterator/transform_iterator.h>

#include "cuda_common.cuh"

struct ExtractTagOp {
    int tag_shift;
    __host__ __device__ __forceinline__ uint32_t operator()(uint64_t k) const {
        return (uint32_t)(k >> tag_shift);
    }
};

static __global__ void revcomp_2bit_kernel(
    const uint64_t * __restrict__ src_bits,
    uint64_t * __restrict__ dst_bits,
    uint64_t len,
    uint32_t n_words
) {
    uint32_t w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= n_words) return;

    uint64_t out_end = (uint64_t)(w + 1U) * 32ULL;
    if (out_end <= len) {
        uint32_t src_start = (uint32_t)(len - out_end);
        uint64_t fwd = load_32_2bit(src_bits, src_start);
        dst_bits[w] = ~reverse_32_2bit_dev(fwd);
    } else {
        uint32_t rem = (uint32_t)(len & 31ULL);
        uint64_t last = 0ULL;
        for (uint32_t k = 0; k < rem; k++) {
            uint32_t src_pos = rem - 1U - k;
            uint64_t comp = (uint64_t)(base_2bit_dev(src_bits, src_pos) ^ 3U);
            last |= comp << (2U * k);
        }
        dst_bits[w] = last;
    }
}

static __global__ void rev_bitmap_kernel(
    const uint32_t * __restrict__ src_bm,
    uint32_t * __restrict__ dst_bm,
    uint64_t len,
    uint32_t n_words
) {
    uint32_t w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= n_words) return;

    uint64_t out_end = (uint64_t)(w + 1U) * 32ULL;
    if (out_end <= len) {
        uint32_t src_start = (uint32_t)(len - out_end);
        uint32_t fwd = load_32_1bit(src_bm, src_start);
        dst_bm[w] = __brev(fwd);
    } else {
        uint32_t rem = (uint32_t)(len & 31ULL);
        uint32_t fwd = __ldg(&src_bm[0]);
        dst_bm[w] = __brev(fwd) >> (32U - rem);
    }
}

/* W3-F12b: Warp-cooperative GPU run-to-bitmap kernel.
 * Boundary partial words use atomicOr; strictly interior words belong exclusively
 * to this run (32-bp aligned splits + disjoint sorted runs) and use direct 0xFFFFFFFFU
 * stores. Long runs (> 16 interior words) are filled cooperatively across all 32 warp lanes. */
static __global__ void fill_runs_bitmap_kernel(
    const GeckoRun * __restrict__ runs,
    uint32_t n_runs,
    uint32_t len,
    uint32_t * __restrict__ bm
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t inner_lo = 0U, inner_hi = 0U;
    bool has_long = false;

    if (idx < n_runs) {
        GeckoRun r = runs[idx];
        uint32_t w_s, mask_s, w_e, mask_e;
        if (gecko_run_word_span(r.start, r.end, len,
                                w_s, mask_s, w_e, mask_e, inner_lo, inner_hi)) {
            atomicOr(&bm[w_s], mask_s);
            if (w_s != w_e) {
                atomicOr(&bm[w_e], mask_e);
                if (inner_hi > inner_lo) {
                    if (inner_hi - inner_lo <= 16U) {
                        for (uint32_t w = inner_lo; w < inner_hi; w++) {
                            bm[w] = 0xFFFFFFFFU;
                        }
                    } else {
                        has_long = true;
                    }
                }
            }
        }
    }

    /* sm_75+ ITS safety: ALL 32 lanes reach __ballot_sync/__shfl_sync — NO early return above! */
    uint32_t lane = threadIdx.x & 31U;
    uint32_t long_mask = __ballot_sync(0xFFFFFFFFU, has_long);
    while (long_mask != 0U) {
        int src = __ffs(long_mask) - 1;
        long_mask &= long_mask - 1U;
        uint32_t lo = __shfl_sync(0xFFFFFFFFU, inner_lo, src);
        uint32_t hi = __shfl_sync(0xFFFFFFFFU, inner_hi, src);
        for (uint32_t w = lo + lane; w < hi; w += 32U) {
            bm[w] = 0xFFFFFFFFU;
        }
    }
}

static __global__ void count_valid_kmers_kernel(
    GpuSeqView seq,
    int wl,
    uint32_t num_chunks,
    uint32_t * __restrict__ chunk_counts
) {
    uint32_t cid = blockIdx.x * blockDim.x + threadIdx.x;
    if (cid >= num_chunks) return;

    uint64_t c_start = (uint64_t)cid * GECKO_EXTRACT_CHUNK_BP;
    uint64_t c_end = c_start + GECKO_EXTRACT_CHUNK_BP;
    if (c_end > seq.len) c_end = seq.len;

    uint64_t scan_from = (c_start >= (uint64_t)(wl - 1)) ? (c_start - (uint64_t)(wl - 1)) : 0ULL;
    if (!seq.invalid_bm) {
        uint64_t first_emit = (c_start >= (uint64_t)(wl - 1)) ? c_start : (uint64_t)(wl - 1);
        chunk_counts[cid] = (c_end > first_emit) ? (uint32_t)(c_end - first_emit) : 0U;
        return;
    }

    int valid = 0;
    uint32_t cnt = 0;
    uint64_t i = scan_from;
    while (i < c_end) {
        uint32_t w_idx = (uint32_t)(i >> 5);
        uint32_t lane = (uint32_t)(i & 31ULL);
        uint32_t avail = (uint32_t)((32U - lane < c_end - i) ? (32U - lane) : (c_end - i));
        uint32_t inv_word = __ldg(&seq.invalid_bm[w_idx]) >> lane;
        uint32_t avail_mask = (avail == 32U) ? 0xFFFFFFFFU : ((1U << avail) - 1U);
        uint32_t inv_masked = inv_word & avail_mask;
        if (inv_masked == avail_mask) {
            valid = 0;
            i += avail;
            continue;
        }
        if (inv_masked == 0U) {
            int prev_valid = valid;
            valid += (int)avail;
            if (valid >= wl) {
                uint64_t reach = (prev_valid >= wl - 1) ? i : (i + (uint64_t)(wl - 1 - prev_valid));
                uint64_t emit_lo = (reach > c_start) ? reach : c_start;
                uint64_t slice_end = i + avail;
                if (slice_end > emit_lo) cnt += (uint32_t)(slice_end - emit_lo);
            }
            i += avail;
            continue;
        }
        for (uint32_t k = 0; k < avail; k++) {
            if (((inv_word >> k) & 1U) == 0U) {
                if (++valid >= wl && (i + k) >= c_start) cnt++;
            } else {
                valid = 0;
            }
        }
        i += avail;
    }
    chunk_counts[cid] = cnt;
}

static __global__ void count_valid_kmers_slab_kernel(
    GpuSeqView seq,
    int wl,
    uint64_t seed_mask,
    uint32_t part_start,
    uint32_t part_end,
    uint32_t num_chunks,
    uint32_t * __restrict__ chunk_counts
) {
    uint32_t cid = blockIdx.x * blockDim.x + threadIdx.x;
    if (cid >= num_chunks) return;

    uint64_t c_start = (uint64_t)cid * GECKO_EXTRACT_CHUNK_BP;
    uint64_t c_end = c_start + GECKO_EXTRACT_CHUNK_BP;
    if (c_end > seq.len) c_end = seq.len;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : mask;
    int rc_shift = 2 * (wl - 1);
    uint64_t scan_from = (c_start >= (uint64_t)(wl - 1)) ? (c_start - (uint64_t)(wl - 1)) : 0ULL;

    uint64_t kf = 0, kr = 0;
    int valid = 0;
    uint32_t cnt = 0;

    uint64_t i = scan_from;
    while (i < c_end) {
        uint32_t w_idx = (uint32_t)(i >> 5);
        uint32_t lane = (uint32_t)(i & 31ULL);
        uint32_t avail = (uint32_t)((32U - lane < c_end - i) ? (32U - lane) : (c_end - i));
        uint32_t inv_word = seq.invalid_bm ? (__ldg(&seq.invalid_bm[w_idx]) >> lane) : 0U;
        uint32_t avail_mask = (avail == 32U) ? 0xFFFFFFFFU : ((1U << avail) - 1U);
        if ((inv_word & avail_mask) == avail_mask) {
            valid = 0;
            kf = 0;
            kr = 0;
            i += avail;
            continue;
        }
        uint64_t bword = __ldg(&seq.bits[w_idx]) >> (lane << 1);
        for (uint32_t k = 0; k < avail; k++) {
            uint64_t b = bword & 3ULL;
            bword >>= 2;
            if (((inv_word >> k) & 1U) == 0U) {
                kf = ((kf << 2) | b) & mask;
                kr = (kr >> 2) | ((b ^ 3ULL) << rc_shift);
                if (++valid >= wl && (i + k) >= c_start) {
                    uint64_t sf = kf & eff_mask;
                    uint64_t sr = kr & eff_mask;
                    uint64_t key = (sf <= sr) ? sf : sr;
                    uint32_t p = (uint32_t)(fmix64(key) >> 56);
                    if (p >= part_start && p < part_end) cnt++;
                }
            } else {
                valid = 0;
                kf = 0;
                kr = 0;
            }
        }
        i += avail;
    }
    chunk_counts[cid] = cnt;
}

/* F15: CTA-cooperative coalesced k-mer extractor (1 CTA = 1 4,096-bp chunk,
 * 4 sequential 1,024-bp waves, compile-time r=0..3 register staging => 0 B stack,
 * 12.1 KiB smem for Hash-CSR / 8.1 KiB for Direct-CSR => 100% smem occupancy). */
template <typename KeyT, bool IS_DIRECT>
static __global__ void extract_kmers_coalesced_kernel(
    GpuSeqView seq,
    int wl,
    uint64_t eff_mask,
    uint32_t part_start,
    uint32_t part_end,
    const uint32_t * __restrict__ chunk_offsets,
    KeyT * __restrict__ out_keys,
    uint32_t * __restrict__ out_vals
) {
    uint32_t cid = blockIdx.x;
    uint32_t chunk_dst = chunk_offsets[cid];
    if (chunk_dst == chunk_offsets[cid + 1U]) return; /* uniform across CTA */

    __shared__ KeyT s_keys[1024];
    __shared__ uint32_t s_vals[1024];
    __shared__ uint32_t s_warp_totals[8];
    __shared__ uint32_t s_warp_prefix[8];
    __shared__ uint32_t s_wave_total;

    uint64_t c_start = (uint64_t)cid * GECKO_EXTRACT_CHUNK_BP;
    uint64_t c_end = c_start + GECKO_EXTRACT_CHUNK_BP;
    if (c_end > seq.len) c_end = seq.len;

    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint32_t valid_mask_wl = (wl == 32) ? 0xFFFFFFFFU : ((1U << wl) - 1U);
    int rev_rshift = 64 - 2 * wl;
    bool part_filter = !IS_DIRECT && (part_start != 0U || part_end != 256U);

    uint32_t tid = threadIdx.x;
    uint32_t warp_id = tid >> 5;
    uint32_t lane = tid & 31U;

    for (uint32_t wave = 0U; wave < 4U; wave++) {
        uint64_t w_start = c_start + (uint64_t)wave * 1024ULL;
        if (w_start >= c_end) break; /* 100% uniform across CTA */
        uint64_t w_end = w_start + 1024ULL;
        if (w_end > c_end) w_end = c_end;

        uint64_t t_start = w_start + (uint64_t)tid * 4ULL;
        KeyT lk[4];
        uint32_t lv[4];
        uint32_t ok_mask = 0U;

        #pragma unroll
        for (uint32_t r = 0U; r < 4U; r++) {
            uint64_t cur_i = t_start + r;
            if (cur_i < w_end && cur_i >= (uint64_t)(wl - 1)) {
                uint32_t pos = (uint32_t)(cur_i + 1ULL - (uint64_t)wl);
                if (!seq.invalid_bm || (load_32_1bit(seq.invalid_bm, pos) & valid_mask_wl) == 0U) {
                    uint64_t w_lsb = load_32_2bit(seq.bits, pos);
                    uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift;
                    uint64_t kr = (~w_lsb) & span_mask;
                    uint64_t sf = kf & eff_mask;
                    uint64_t sr = kr & eff_mask;
                    uint32_t strand = (sf <= sr) ? 0U : 1U;
                    uint64_t canon = strand ? sr : sf;
                    if constexpr (IS_DIRECT) {
                        lk[r] = (KeyT)canon;
                        lv[r] = (pos << 1) | strand;
                        ok_mask |= (1U << r);
                    } else {
                        uint64_t h = fmix64(canon);
                        if (!part_filter ||
                            (((uint32_t)(h >> 56) >= part_start) && ((uint32_t)(h >> 56) < part_end))) {
                            lk[r] = (KeyT)h;
                            lv[r] = (pos << 1) | strand;
                            ok_mask |= (1U << r);
                        }
                    }
                }
            }
        }

        uint32_t local_cnt = (uint32_t)__popc(ok_mask);

        uint32_t scan = local_cnt;
        #pragma unroll
        for (int off = 1; off < 32; off <<= 1) {
            uint32_t n = __shfl_up_sync(0xFFFFFFFFU, scan, off);
            if (lane >= (uint32_t)off) scan += n;
        }
        uint32_t warp_excl = scan - local_cnt;
        if (lane == 31U) s_warp_totals[warp_id] = scan;
        __syncthreads(); /* Sync #1: unconditional */

        if (warp_id == 0U) {
            uint32_t wt = (lane < 8U) ? s_warp_totals[lane] : 0U;
            uint32_t orig_wt = wt;
            #pragma unroll
            for (int off = 1; off < 8; off <<= 1) {
                uint32_t n = __shfl_up_sync(0xFFFFFFFFU, wt, off);
                if (lane >= (uint32_t)off) wt += n;
            }
            if (lane < 8U) s_warp_prefix[lane] = wt - orig_wt;
            if (lane == 7U) s_wave_total = wt;
        }
        __syncthreads(); /* Sync #2: unconditional */

        uint32_t my_off = s_warp_prefix[warp_id] + warp_excl;
        uint32_t wave_total = s_wave_total;

        uint32_t dst_s = my_off;
        #pragma unroll
        for (uint32_t r = 0U; r < 4U; r++) {
            if ((ok_mask >> r) & 1U) {
                s_keys[dst_s] = lk[r];
                s_vals[dst_s] = lv[r];
                dst_s++;
            }
        }
        __syncthreads(); /* Sync #3: unconditional */

        for (uint32_t idx = tid; idx < wave_total; idx += 256U) {
            out_keys[chunk_dst + idx] = s_keys[idx];
            out_vals[chunk_dst + idx] = s_vals[idx];
        }
        chunk_dst += wave_total;
        __syncthreads(); /* Sync #4: unconditional before next wave reuses smem */
    }
}

/* W3-F13a: Fused frequency filter + inverted ref_non_unique_bm + word-coalesced dropped_bitmap.
 * - Non-leaders of run_len >= 2 mark ref_non_unique_bm in O(1) parallel.
 * - ~90% unique non-palindromic k-mers (run_len == 1) do ZERO vals[i] loads and ZERO atomics.
 * - Accepts keep_flags == nullptr for gpu_build_global_dropped_bitmap (eliminates d_keep). */
template <typename KeyT, bool IS_DIRECT>
static __global__ void mark_freq_filter_kernel(
    const KeyT * __restrict__ keys,
    const uint32_t * __restrict__ vals,
    uint32_t n,
    uint32_t freq_limit,
    int wl,
    uint64_t eff_mask,
    uint8_t * __restrict__ keep_flags,         /* may be nullptr (global_dropped_bm mode) */
    uint32_t * __restrict__ ref_non_unique_bm, /* single-slab contiguous Hash-CSR only */
    uint32_t * __restrict__ dropped_bitmap,    /* non-null when building dropped bitmap */
    uint32_t * __restrict__ any_dropped        /* non-null only in gpu_build_global_dropped_bitmap */
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    KeyT k = keys[i];
    bool same_prev = (i > 0U && keys[i - 1U] == k);
    if (same_prev) {
        if (ref_non_unique_bm) {
            uint32_t pos = vals[i] >> 1;
            atomicOr(&ref_non_unique_bm[pos >> 5], 1U << (pos & 31U));
        }
        return;
    }

    bool same_next = (i + 1U < n && keys[i + 1U] == k);
    if (!same_next) {
        if (keep_flags) keep_flags[i] = 1;
        if (ref_non_unique_bm && ((wl & 1) == 0)) {
            bool is_pal = IS_DIRECT
                ? is_palindromic_canon((uint64_t)k, wl, eff_mask)
                : is_palindromic_canon(fmix64_inv((uint64_t)k), wl, eff_mask);
            if (is_pal) {
                uint32_t pos = vals[i] >> 1;
                atomicOr(&ref_non_unique_bm[pos >> 5], 1U << (pos & 31U));
            }
        }
        return;
    }

    /* Run leader of run_len >= 2: mark leader's own bit in ref_non_unique_bm */
    if (ref_non_unique_bm) {
        uint32_t pos = vals[i] >> 1;
        atomicOr(&ref_non_unique_bm[pos >> 5], 1U << (pos & 31U));
    }

    uint32_t step = 1U;
    uint32_t hi = i + 1U;
    while (hi < n && keys[hi] == k) {
        step <<= 1;
        uint32_t next = i + step;
        hi = (next < n && next > hi) ? next : n;
    }
    uint32_t lo = i + (step >> 1);
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1);
        if (keys[mid] == k) lo = mid + 1U;
        else hi = mid;
    }
    uint32_t run_end = lo;
    uint32_t run_len = run_end - i;

    if (run_len <= freq_limit) {
        if (keep_flags) {
            for (uint32_t m = i; m < run_end; m++) keep_flags[m] = 1;
        }
        return;
    }

    uint32_t c0 = 0U, c1 = 0U;
    for (uint32_t m = i; m < run_end; m++) {
        if (vals[m] & 1U) c1++;
        else c0++;
        /* break (NOT return) so dropped_bitmap + any_dropped are still marked for every element! */
        if (c0 > freq_limit && c1 > freq_limit) break;
    }
    uint8_t k0 = (c0 <= freq_limit) ? 1 : 0;
    uint8_t k1 = (c1 <= freq_limit) ? 1 : 0;

    if (dropped_bitmap && (!k0 || !k1)) {
        if (any_dropped && *any_dropped == 0U) {
            *any_dropped = 1U;
        }
        /* Single merged pass: set keep_flags (pre-zeroed) + word-coalesced dropped_bitmap REDG */
        uint32_t cur_w = 0xFFFFFFFFU;
        uint32_t cur_mask = 0U;
        for (uint32_t m = i; m < run_end; m++) {
            uint32_t v = vals[m];
            uint8_t keep = (v & 1U) ? k1 : k0;
            if (keep_flags && keep) keep_flags[m] = 1;
            if (!keep) {
                uint32_t pos = v >> 1;
                uint32_t w = pos >> 5;
                uint32_t bit = 1U << (pos & 31U);
                if (w == cur_w) {
                    cur_mask |= bit;
                } else {
                    if (cur_mask != 0U) atomicOr(&dropped_bitmap[cur_w], cur_mask);
                    cur_w = w;
                    cur_mask = bit;
                }
            }
        }
        if (cur_mask != 0U) atomicOr(&dropped_bitmap[cur_w], cur_mask);
    } else if (keep_flags && (k0 || k1)) {
        for (uint32_t m = i; m < run_end; m++) {
            uint8_t keep = (vals[m] & 1U) ? k1 : k0;
            if (keep) keep_flags[m] = 1;
        }
    }
}

/* F15 + W3-F13a: Warp-aggregated bucket count kernel via sm_70+ __match_any_sync.
 * Also writes 32-bit tags and inverted 1-bit ref_non_unique_bm inline on the zero-drop fast path. */
template <typename KeyT, bool IS_DIRECT>
static __global__ void count_kept_per_bucket_kernel(
    const KeyT * __restrict__ keys,
    const uint32_t * __restrict__ vals,       /* non-null ONLY when keep_flags == nullptr */
    const uint8_t * __restrict__ keep_flags,  /* non-null when freq_limit < n_valid */
    uint32_t n,
    int shift_bits,
    int tag_shift,
    uint32_t bucket_base,
    int wl,
    uint64_t eff_mask,
    uint32_t * __restrict__ bucket_counts,    /* length num_buckets + 1, pre-zeroed */
    uint32_t * __restrict__ out_tags,         /* non-null ONLY for Hash-CSR && no_freq_drops */
    uint32_t * __restrict__ ref_non_unique_bm /* non-null ONLY for single-slab contiguous Hash-CSR && no_freq_drops */
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (keep_flags && keep_flags[i] == 0) return;

    KeyT k = keys[i];
    uint32_t b = IS_DIRECT
                     ? (uint32_t)k
                     : ((uint32_t)((uint64_t)k >> shift_bits) - bucket_base);

    uint32_t active = __activemask();
    uint32_t same = __match_any_sync(active, b);
    int leader = __ffs(same) - 1;
    if ((int)(threadIdx.x & 31U) == leader) {
        atomicAdd(&bucket_counts[b], (uint32_t)__popc(same));
    }

    if constexpr (!IS_DIRECT) {
        if (out_tags) {
            out_tags[i] = (uint32_t)((uint64_t)k >> tag_shift);
        }
    }

    if (ref_non_unique_bm && vals) {
        bool non_unique = (i > 0U && keys[i - 1U] == k) || (i + 1U < n && keys[i + 1U] == k);
        if (!non_unique && ((wl & 1) == 0)) {
            non_unique = IS_DIRECT
                ? is_palindromic_canon((uint64_t)k, wl, eff_mask)
                : is_palindromic_canon(fmix64_inv((uint64_t)k), wl, eff_mask);
        }
        if (non_unique) {
            uint32_t pos = vals[i] >> 1;
            atomicOr(&ref_non_unique_bm[pos >> 5], 1U << (pos & 31U));
        }
    }
}

int gpu_fill_runs_bitmap(
    const GeckoRun *d_runs,
    uint32_t n_runs,
    uint64_t len,
    uint32_t *d_bm,
    cudaStream_t stream
) {
    if (!d_runs || !d_bm || n_runs == 0U || len == 0ULL) return 0;
    uint32_t grid = (n_runs + 255U) / 256U;
    fill_runs_bitmap_kernel<<<grid, 256, 0, stream>>>(d_runs, n_runs, (uint32_t)len, d_bm);
    CUDA_CHECK(cudaGetLastError());
    return 0;
}

void gpu_build_revcomp_2bit(
    const uint64_t *d_src_bits,
    uint64_t *d_dst_bits,
    const uint32_t *d_src_special_bm,
    uint32_t *d_dst_special_bm,
    uint64_t len,
    GeckoCudaProfile *prof
) {
    if (len == 0) return;
    CudaEventScope ev(prof, prof ? &prof->revcomp_ms : nullptr);
    uint32_t n_words = (uint32_t)((len + 31ULL) / 32ULL);
    uint32_t grid = (n_words + 255U) / 256U;
    revcomp_2bit_kernel<<<grid, 256>>>(d_src_bits, d_dst_bits, len, n_words);
    CUDA_CHECK_VOID(cudaGetLastError());
    if (d_src_special_bm && d_dst_special_bm) {
        rev_bitmap_kernel<<<grid, 256>>>(d_src_special_bm, d_dst_special_bm, len, n_words);
        CUDA_CHECK_VOID(cudaGetLastError());
    }
}

uint64_t gpu_count_valid_kmers(
    const GpuSeqView &seqY,
    int wl,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t *&out_d_chunk_offsets,
    uint32_t &out_num_chunks,
    GeckoCudaProfile *prof
) {
    out_d_chunk_offsets = nullptr;
    out_num_chunks = 0U;
    uint64_t lenY = seqY.len;
    if (lenY < (uint64_t)wl) return 0ULL;
    uint32_t num_chunks = (uint32_t)((lenY + GECKO_EXTRACT_CHUNK_BP - 1) / GECKO_EXTRACT_CHUNK_BP);
    if (num_chunks == 0U) return 0ULL;

    uint32_t *d_chunk_counts = nullptr;
    uint32_t *d_chunk_offsets = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_chunk_counts, (size_t)num_chunks * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_chunk_offsets, (size_t)(num_chunks + 1U) * sizeof(uint32_t)));

    size_t temp_bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, temp_bytes, d_chunk_counts, d_chunk_offsets, (int)num_chunks);
    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, temp_bytes);

    {
        CudaEventScope ev(prof, prof ? &prof->count_valid_ms : nullptr);
        count_valid_kmers_kernel<<<(num_chunks + 255U) / 256U, 256>>>(
            seqY, wl, num_chunks, d_chunk_counts);
        CUDA_CHECK_VOID(cudaGetLastError());
        cub::DeviceScan::ExclusiveSum(d_cub_temp, temp_bytes, d_chunk_counts, d_chunk_offsets, (int)num_chunks);
    }

    uint32_t last_off = 0U, last_cnt = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&last_off, d_chunk_offsets + (num_chunks - 1U), sizeof(uint32_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK_VOID(cudaMemcpy(&last_cnt, d_chunk_counts + (num_chunks - 1U), sizeof(uint32_t), cudaMemcpyDeviceToHost));
    uint32_t total = last_off + last_cnt;
    CUDA_CHECK_VOID(cudaMemcpy(d_chunk_offsets + num_chunks, &total, sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK_VOID(cudaFree(d_chunk_counts));

    out_d_chunk_offsets = d_chunk_offsets;
    out_num_chunks = num_chunks;
    return (uint64_t)total;
}

static inline int compute_slab_sort_end_bit(uint32_t part_start, uint32_t part_end) {
    if (part_start == 0U && part_end == 256U) return 64;
    if (part_end <= part_start) return 64;
    uint32_t diff_mask = part_start ^ (part_end - 1U);
    int shared_top_bits = (diff_mask == 0U) ? 8 : (__builtin_clz(diff_mask) - 24);
    if (shared_top_bits < 0) shared_top_bits = 0;
    if (shared_top_bits > 8) shared_top_bits = 8;
    return 64 - shared_top_bits;
}

uint32_t *gpu_build_global_dropped_bitmap(
    const GpuSeqView &seqY,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit,
    uint32_t num_slabs,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    GeckoCudaProfile *prof
) {
    uint64_t lenY = seqY.len;
    if (freq_limit == 0U || lenY < (uint64_t)wl || num_slabs <= 1U) return nullptr;

    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : span_mask;
    if (eff_mask != span_mask) return nullptr; /* spaced seeds always skip L/R */

    uint32_t num_chunks = (uint32_t)((lenY + GECKO_EXTRACT_CHUNK_BP - 1) / GECKO_EXTRACT_CHUNK_BP);
    if (num_chunks == 0U) return nullptr;

    size_t bm_words = (size_t)((lenY + 31ULL) / 32ULL) + 2;
    uint32_t *d_global_bm = nullptr;
    uint32_t *d_any_dropped = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_global_bm, bm_words * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMemset(d_global_bm, 0, bm_words * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_any_dropped, sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMemset(d_any_dropped, 0, sizeof(uint32_t)));

    uint32_t *d_chunk_counts = nullptr;
    uint32_t *d_chunk_offsets = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_chunk_counts, (size_t)num_chunks * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_chunk_offsets, (size_t)(num_chunks + 1U) * sizeof(uint32_t)));

    for (uint32_t slab = 0U; slab < num_slabs; slab++) {
        uint32_t part_start = (slab * 256U) / num_slabs;
        uint32_t part_end = ((slab + 1U) * 256U) / num_slabs;
        if (part_start >= part_end) continue;

        CudaEventScope ev_slab(prof, prof ? &prof->global_dropped_ms : nullptr);

        count_valid_kmers_slab_kernel<<<(num_chunks + 255U) / 256U, 256>>>(
            seqY, wl, eff_mask, part_start, part_end, num_chunks, d_chunk_counts);
        CUDA_CHECK_VOID(cudaGetLastError());

        size_t scan_bytes = 0;
        cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_chunk_counts, d_chunk_offsets, (int)num_chunks);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);
        cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_chunk_counts, d_chunk_offsets, (int)num_chunks);

        uint32_t last_off = 0U, last_cnt = 0U;
        CUDA_CHECK_VOID(cudaMemcpy(&last_off, d_chunk_offsets + (num_chunks - 1U), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        CUDA_CHECK_VOID(cudaMemcpy(&last_cnt, d_chunk_counts + (num_chunks - 1U), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        uint32_t n_valid = last_off + last_cnt;
        if (n_valid == 0U) continue;
        CUDA_CHECK_VOID(cudaMemcpy(d_chunk_offsets + num_chunks, &n_valid, sizeof(uint32_t), cudaMemcpyHostToDevice));

        uint64_t *d_keys_a = nullptr, *d_keys_b = nullptr;
        uint32_t *d_vals_a = nullptr, *d_vals_b = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_keys_a, (size_t)n_valid * sizeof(uint64_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_vals_a, (size_t)n_valid * sizeof(uint32_t)));

        extract_kmers_coalesced_kernel<uint64_t, false><<<num_chunks, 256>>>(
            seqY, wl, eff_mask, part_start, part_end, d_chunk_offsets, d_keys_a, d_vals_a);
        CUDA_CHECK_VOID(cudaGetLastError());

        CUDA_CHECK_VOID(cudaMalloc(&d_keys_b, (size_t)n_valid * sizeof(uint64_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_vals_b, (size_t)n_valid * sizeof(uint32_t)));
        cub::DoubleBuffer<uint64_t> db_keys(d_keys_a, d_keys_b);
        cub::DoubleBuffer<uint32_t> db_vals(d_vals_a, d_vals_b);

        if (n_valid > 1U) {
            int sort_end_bit = compute_slab_sort_end_bit(part_start, part_end);
            size_t sort_bytes = 0;
            cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes, db_keys, db_vals, (int)n_valid, 0, sort_end_bit);
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sort_bytes);
            cub::DeviceRadixSort::SortPairs(d_cub_temp, sort_bytes, db_keys, db_vals, (int)n_valid, 0, sort_end_bit);
        }

        uint64_t *d_sorted_keys = db_keys.Current();
        uint32_t *d_sorted_vals = db_vals.Current();
        CUDA_CHECK_VOID(cudaFree(db_keys.Alternate()));
        CUDA_CHECK_VOID(cudaFree(db_vals.Alternate()));

        /* W3-F13a: Fused single-pass global dropped bitmap marking (zero d_keep allocation!) */
        mark_freq_filter_kernel<uint64_t, false><<<(n_valid + 255U) / 256U, 256>>>(
            d_sorted_keys, d_sorted_vals, n_valid, freq_limit, wl, eff_mask,
            nullptr, nullptr, d_global_bm, d_any_dropped);
        CUDA_CHECK_VOID(cudaGetLastError());

        CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
        CUDA_CHECK_VOID(cudaFree(d_sorted_vals));
    }

    CUDA_CHECK_VOID(cudaFree(d_chunk_counts));
    CUDA_CHECK_VOID(cudaFree(d_chunk_offsets));

    uint32_t h_any_dropped = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&h_any_dropped, d_any_dropped, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK_VOID(cudaFree(d_any_dropped));

    if (h_any_dropped == 0U) {
        CUDA_CHECK_VOID(cudaFree(d_global_bm));
        return nullptr;
    }
    return d_global_bm;
}

GpuRefIndex gpu_build_ref_index(
    const GpuSeqView &seqY,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t part_start,
    uint32_t part_end,
    const uint32_t *precomp_chunk_offsets,
    uint32_t precomp_n_valid,
    bool build_dropped_bm,
    GeckoCudaProfile *prof
) {
    GpuRefIndex idx{};
    idx.part_start = part_start;
    idx.part_end = part_end;

    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : span_mask;
    bool is_spaced = (eff_mask != span_mask);
    bool part_filter = (part_start != 0U || part_end != 256U);
    uint32_t n_parts = (part_end > part_start) ? (part_end - part_start) : 1U;
    uint64_t lenY = seqY.len;
    idx.has_global_dropped_bm = (!part_filter || freq_limit == 0U || (uint64_t)freq_limit >= lenY);
    uint32_t num_chunks = (uint32_t)((lenY + GECKO_EXTRACT_CHUNK_BP - 1) / GECKO_EXTRACT_CHUNK_BP);

    const uint32_t *d_chunk_offsets_ro = nullptr;
    uint32_t *d_chunk_offsets_owned = nullptr;

    if (precomp_chunk_offsets != nullptr && !part_filter) {
        d_chunk_offsets_ro = precomp_chunk_offsets;
        idx.n_valid = precomp_n_valid;
    } else if (num_chunks > 0U) {
        uint32_t *d_chunk_counts = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_chunk_counts, (size_t)num_chunks * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_chunk_offsets_owned, (size_t)(num_chunks + 1U) * sizeof(uint32_t)));

        size_t scan_bytes = 0;
        cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_chunk_counts, d_chunk_offsets_owned, (int)num_chunks);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);

        {
            CudaEventScope ev(prof, prof ? &prof->count_valid_ms : nullptr);
            if (part_filter) {
                count_valid_kmers_slab_kernel<<<(num_chunks + 255U) / 256U, 256>>>(
                    seqY, wl, eff_mask, part_start, part_end, num_chunks, d_chunk_counts);
            } else {
                count_valid_kmers_kernel<<<(num_chunks + 255U) / 256U, 256>>>(
                    seqY, wl, num_chunks, d_chunk_counts);
            }
            CUDA_CHECK_VOID(cudaGetLastError());
            cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_chunk_counts, d_chunk_offsets_owned, (int)num_chunks);
        }

        uint32_t last_off = 0U, last_cnt = 0U;
        CUDA_CHECK_VOID(cudaMemcpy(&last_off, d_chunk_offsets_owned + (num_chunks - 1U), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        CUDA_CHECK_VOID(cudaMemcpy(&last_cnt, d_chunk_counts + (num_chunks - 1U), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        idx.n_valid = last_off + last_cnt;
        CUDA_CHECK_VOID(cudaMemcpy(d_chunk_offsets_owned + num_chunks, &idx.n_valid, sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK_VOID(cudaFree(d_chunk_counts));
        d_chunk_offsets_ro = d_chunk_offsets_owned;
    }

    /* Empty sequence / slab / freq_limit == 0 fast exit (count <= 0 => 0 k-mers kept) */
    if (freq_limit == 0U || idx.n_valid == 0U) {
        if (d_chunk_offsets_owned) CUDA_CHECK_VOID(cudaFree(d_chunk_offsets_owned));
        int bucket_bits = 14;
        uint32_t sub_bits = (uint32_t)(bucket_bits - 8);
        uint32_t bucket_base = part_start << sub_bits;
        uint32_t num_buckets = n_parts << sub_bits;
        CUDA_CHECK_VOID(cudaMalloc(&idx.d_bucket_off, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMemset(idx.d_bucket_off, 0, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        idx.d_tags = nullptr;
        idx.d_vals = nullptr;
        idx.d_ref_non_unique_bm = nullptr;
        idx.d_dropped_bitmap = nullptr;
        idx.owns_dropped_bitmap = false;
        idx.has_global_dropped_bm = true;
        idx.n_kept = 0U;
        idx.num_buckets = num_buckets;
        idx.bucket_bits = bucket_bits;
        idx.shift_bits = 64 - bucket_bits;
        idx.tag_shift = idx.shift_bits - 32;
        idx.bucket_base = bucket_base;
        idx.is_direct_csr = false;
        return idx;
    }

    bool no_freq_drops = (freq_limit >= idx.n_valid);

    /* F14: Direct-CSR for contiguous wl in 8..=12 (single-slab, non-empty).
     * Never check seed_mask == 0ULL (bridge normalizes seed_mask = span_mask)! */
    bool is_direct = !is_spaced && !part_filter && idx.n_valid > 0U
                     && wl >= 8 && wl <= 12
                     && (wl <= 10 || ((1ULL << (2 * wl)) <= 4ULL * (uint64_t)idx.n_valid));

    if (is_direct) {
        size_t n_alloc = (size_t)idx.n_valid;
        uint32_t *d_keys_a = nullptr, *d_keys_b = nullptr;
        uint32_t *d_vals_a = nullptr, *d_vals_b = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_keys_a, n_alloc * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_vals_a, n_alloc * sizeof(uint32_t)));

        {
            CudaEventScope ev(prof, prof ? &prof->extract_kmers_ms : nullptr);
            extract_kmers_coalesced_kernel<uint32_t, true><<<num_chunks, 256>>>(
                seqY, wl, eff_mask, 0U, 256U, d_chunk_offsets_ro, d_keys_a, d_vals_a);
            CUDA_CHECK_VOID(cudaGetLastError());
        }
        if (d_chunk_offsets_owned) {
            CUDA_CHECK_VOID(cudaFree(d_chunk_offsets_owned));
            d_chunk_offsets_owned = nullptr;
        }

        CUDA_CHECK_VOID(cudaMalloc(&d_keys_b, n_alloc * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_vals_b, n_alloc * sizeof(uint32_t)));
        cub::DoubleBuffer<uint32_t> db_keys(d_keys_a, d_keys_b);
        cub::DoubleBuffer<uint32_t> db_vals(d_vals_a, d_vals_b);

        if (idx.n_valid > 1U) {
            size_t sort_bytes = 0;
            cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes, db_keys, db_vals, (int)idx.n_valid, 0, 2 * wl);
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sort_bytes);
            CudaEventScope ev(prof, prof ? &prof->sort_pairs_ms : nullptr);
            cub::DeviceRadixSort::SortPairs(d_cub_temp, sort_bytes, db_keys, db_vals, (int)idx.n_valid, 0, 2 * wl);
        }

        uint32_t *d_sorted_keys = db_keys.Current();
        uint32_t *d_sorted_vals = db_vals.Current();
        uint32_t *d_alt_vals = db_vals.Alternate();
        CUDA_CHECK_VOID(cudaFree(db_keys.Alternate()));

        uint32_t num_buckets = 1U << (2 * wl);
        uint32_t *d_bucket_off = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_bucket_off, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMemset(d_bucket_off, 0, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));

        uint32_t *d_final_vals = nullptr;
        if (!no_freq_drops) {
            uint8_t *d_keep = nullptr;
            uint32_t *d_num_out = nullptr;
            CUDA_CHECK_VOID(cudaMalloc(&d_keep, n_alloc * sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMemset(d_keep, 0, n_alloc * sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMalloc(&d_num_out, sizeof(uint32_t)));

            if (build_dropped_bm && lenY > 0ULL) {
                size_t bm_words = (size_t)((lenY + 31ULL) / 32ULL) + 2;
                CUDA_CHECK_VOID(cudaMalloc(&idx.d_dropped_bitmap, bm_words * sizeof(uint32_t)));
                CUDA_CHECK_VOID(cudaMemset(idx.d_dropped_bitmap, 0, bm_words * sizeof(uint32_t)));
                idx.owns_dropped_bitmap = true;
            }

            {
                CudaEventScope ev(prof, prof ? &prof->freq_filter_ms : nullptr);
                mark_freq_filter_kernel<uint32_t, true><<<(idx.n_valid + 255U) / 256U, 256>>>(
                    d_sorted_keys, d_sorted_vals, idx.n_valid, freq_limit, wl, eff_mask,
                    d_keep, nullptr, idx.d_dropped_bitmap, nullptr);
                CUDA_CHECK_VOID(cudaGetLastError());
            }

            size_t sel_bytes = 0;
            cub::DeviceSelect::Flagged(nullptr, sel_bytes, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sel_bytes);
            {
                CudaEventScope ev(prof, prof ? &prof->select_vals_ms : nullptr);
                cub::DeviceSelect::Flagged(d_cub_temp, sel_bytes, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
            }
            CUDA_CHECK_VOID(cudaMemcpy(&idx.n_kept, d_num_out, sizeof(uint32_t), cudaMemcpyDeviceToHost));

            if (idx.d_dropped_bitmap && idx.n_kept == idx.n_valid) {
                CUDA_CHECK_VOID(cudaFree(idx.d_dropped_bitmap));
                idx.d_dropped_bitmap = nullptr;
                idx.owns_dropped_bitmap = false;
            }

            CUDA_CHECK_VOID(cudaFree(d_sorted_vals));
            CUDA_CHECK_VOID(cudaFree(d_num_out));
            d_final_vals = d_alt_vals;

            if (idx.n_kept > 0U) {
                size_t scan_bytes = 0;
                cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
                gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);
                {
                    CudaEventScope ev(prof, prof ? &prof->bucket_build_ms : nullptr);
                    count_kept_per_bucket_kernel<uint32_t, true><<<(idx.n_valid + 255U) / 256U, 256>>>(
                        d_sorted_keys, nullptr, d_keep, idx.n_valid,
                        0, 0, 0U, wl, eff_mask, d_bucket_off, nullptr, nullptr);
                    CUDA_CHECK_VOID(cudaGetLastError());
                    cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
                }
                CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
                CUDA_CHECK_VOID(cudaFree(d_keep));
            } else {
                CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
                CUDA_CHECK_VOID(cudaFree(d_keep));
            }
        } else {
            CUDA_CHECK_VOID(cudaFree(d_alt_vals));
            d_final_vals = d_sorted_vals;
            idx.n_kept = idx.n_valid;

            size_t scan_bytes = 0;
            cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);
            {
                CudaEventScope ev(prof, prof ? &prof->bucket_build_ms : nullptr);
                count_kept_per_bucket_kernel<uint32_t, true><<<(idx.n_valid + 255U) / 256U, 256>>>(
                    d_sorted_keys, nullptr, nullptr, idx.n_valid,
                    0, 0, 0U, wl, eff_mask, d_bucket_off, nullptr, nullptr);
                CUDA_CHECK_VOID(cudaGetLastError());
                cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
            }
            CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
        }

        idx.d_tags = nullptr;
        idx.d_vals = d_final_vals;
        idx.d_bucket_off = d_bucket_off;
        idx.d_ref_non_unique_bm = nullptr;
        idx.num_buckets = num_buckets;
        idx.bucket_bits = 2 * wl;
        idx.shift_bits = 0;
        idx.tag_shift = 0;
        idx.bucket_base = 0U;
        idx.is_direct_csr = true;
        return idx;
    }

    /* ---- Hash-CSR Path (uint64_t keys, 32-bit tags, W3-F13a inverted ref_non_unique_bm) ---- */
    size_t n_alloc = (size_t)idx.n_valid;
    uint64_t *d_keys_a = nullptr, *d_keys_b = nullptr;
    uint32_t *d_vals_a = nullptr, *d_vals_b = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_keys_a, n_alloc * sizeof(uint64_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_vals_a, n_alloc * sizeof(uint32_t)));

    {
        CudaEventScope ev(prof, prof ? &prof->extract_kmers_ms : nullptr);
        extract_kmers_coalesced_kernel<uint64_t, false><<<num_chunks, 256>>>(
            seqY, wl, eff_mask, part_start, part_end, d_chunk_offsets_ro, d_keys_a, d_vals_a);
        CUDA_CHECK_VOID(cudaGetLastError());
    }
    if (d_chunk_offsets_owned) {
        CUDA_CHECK_VOID(cudaFree(d_chunk_offsets_owned));
        d_chunk_offsets_owned = nullptr;
    }

    CUDA_CHECK_VOID(cudaMalloc(&d_keys_b, n_alloc * sizeof(uint64_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_vals_b, n_alloc * sizeof(uint32_t)));
    cub::DoubleBuffer<uint64_t> db_keys(d_keys_a, d_keys_b);
    cub::DoubleBuffer<uint32_t> db_vals(d_vals_a, d_vals_b);

    if (idx.n_valid > 1U) {
        int sort_end_bit = compute_slab_sort_end_bit(part_start, part_end);
        size_t sort_bytes = 0;
        cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes, db_keys, db_vals, (int)idx.n_valid, 0, sort_end_bit);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sort_bytes);
        CudaEventScope ev(prof, prof ? &prof->sort_pairs_ms : nullptr);
        cub::DeviceRadixSort::SortPairs(d_cub_temp, sort_bytes, db_keys, db_vals, (int)idx.n_valid, 0, sort_end_bit);
    }

    uint64_t *d_sorted_keys = db_keys.Current();
    uint32_t *d_sorted_vals = db_vals.Current();
    uint32_t *d_alt_vals = db_vals.Alternate();
    /* F15: Free alternate 64-bit key buffer immediately — we NEVER compact 64-bit keys! */
    CUDA_CHECK_VOID(cudaFree(db_keys.Alternate()));

    /* W3-F13a: Allocate inverted 1-bit non-unique/palindromic ref k-mer bitmap STRICTLY for
     * single-slab contiguous Hash-CSR (!part_filter guarantees every valid ref k-mer is present
     * in d_sorted_keys, so out-of-slab positions can never be falsely treated as unique). */
    if (!is_spaced && !part_filter && lenY > 0ULL) {
        size_t bm_words = (size_t)((lenY + 31ULL) / 32ULL) + 2;
        CUDA_CHECK_VOID(cudaMalloc(&idx.d_ref_non_unique_bm, bm_words * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMemset(idx.d_ref_non_unique_bm, 0, bm_words * sizeof(uint32_t)));
    }

    uint32_t *d_final_vals = nullptr;
    uint32_t *d_bucket_off = nullptr;
    uint32_t *d_tags = nullptr;
    int bucket_bits = 14;
    uint32_t num_buckets = 0U;
    int shift_bits = 0;
    int tag_shift = 0;
    uint32_t bucket_base = 0U;

    if (!no_freq_drops) {
        uint8_t *d_keep = nullptr;
        uint32_t *d_num_out = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_keep, n_alloc * sizeof(uint8_t)));
        CUDA_CHECK_VOID(cudaMemset(d_keep, 0, n_alloc * sizeof(uint8_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_num_out, sizeof(uint32_t)));

        if (build_dropped_bm && !part_filter && !is_spaced && lenY > 0ULL) {
            size_t bm_words = (size_t)((lenY + 31ULL) / 32ULL) + 2;
            CUDA_CHECK_VOID(cudaMalloc(&idx.d_dropped_bitmap, bm_words * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMemset(idx.d_dropped_bitmap, 0, bm_words * sizeof(uint32_t)));
            idx.owns_dropped_bitmap = true;
        }

        {
            CudaEventScope ev(prof, prof ? &prof->freq_filter_ms : nullptr);
            mark_freq_filter_kernel<uint64_t, false><<<(idx.n_valid + 255U) / 256U, 256>>>(
                d_sorted_keys, d_sorted_vals, idx.n_valid, freq_limit, wl, eff_mask,
                d_keep, idx.d_ref_non_unique_bm, idx.d_dropped_bitmap, nullptr);
            CUDA_CHECK_VOID(cudaGetLastError());
        }

        size_t sel_val_bytes = 0;
        cub::DeviceSelect::Flagged(nullptr, sel_val_bytes, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sel_val_bytes);
        {
            CudaEventScope ev(prof, prof ? &prof->select_vals_ms : nullptr);
            cub::DeviceSelect::Flagged(d_cub_temp, sel_val_bytes, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
        }
        CUDA_CHECK_VOID(cudaMemcpy(&idx.n_kept, d_num_out, sizeof(uint32_t), cudaMemcpyDeviceToHost));

        if (idx.d_dropped_bitmap && idx.n_kept == idx.n_valid) {
            CUDA_CHECK_VOID(cudaFree(idx.d_dropped_bitmap));
            idx.d_dropped_bitmap = nullptr;
            idx.owns_dropped_bitmap = false;
        }

        /* Release uncompacted vals BEFORE allocating d_bucket_off & d_tags */
        CUDA_CHECK_VOID(cudaFree(d_sorted_vals));
        d_final_vals = d_alt_vals;

        uint64_t equiv_kept = ((uint64_t)idx.n_kept * 256ULL) / (uint64_t)n_parts;
        while (bucket_bits < 26 && (1ULL << bucket_bits) * 4ULL < equiv_kept) bucket_bits++;
        uint32_t sub_bits = (uint32_t)(bucket_bits - 8);
        bucket_base = part_start << sub_bits;
        num_buckets = n_parts << sub_bits;
        shift_bits = 64 - bucket_bits;
        tag_shift = shift_bits - 32;

        CUDA_CHECK_VOID(cudaMalloc(&d_bucket_off, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMemset(d_bucket_off, 0, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        size_t tag_alloc = (idx.n_kept > 0U) ? (size_t)idx.n_kept : 1;
        CUDA_CHECK_VOID(cudaMalloc(&d_tags, tag_alloc * sizeof(uint32_t)));

        if (idx.n_kept > 0U) {
            auto tag_iter = thrust::make_transform_iterator(d_sorted_keys, ExtractTagOp{tag_shift});
            size_t sel_tag_bytes = 0;
            cub::DeviceSelect::Flagged(nullptr, sel_tag_bytes, tag_iter, d_keep, d_tags, d_num_out, (int)idx.n_valid);
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sel_tag_bytes);
            {
                CudaEventScope ev(prof, prof ? &prof->select_tags_ms : nullptr);
                cub::DeviceSelect::Flagged(d_cub_temp, sel_tag_bytes, tag_iter, d_keep, d_tags, d_num_out, (int)idx.n_valid);
            }

            size_t scan_bytes = 0;
            cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);
            {
                CudaEventScope ev(prof, prof ? &prof->bucket_build_ms : nullptr);
                count_kept_per_bucket_kernel<uint64_t, false><<<(idx.n_valid + 255U) / 256U, 256>>>(
                    d_sorted_keys, nullptr, d_keep, idx.n_valid,
                    shift_bits, tag_shift, bucket_base, wl, eff_mask,
                    d_bucket_off, nullptr, nullptr);
                CUDA_CHECK_VOID(cudaGetLastError());
                cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
            }

            CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
            CUDA_CHECK_VOID(cudaFree(d_keep));
            CUDA_CHECK_VOID(cudaFree(d_num_out));
        } else {
            CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
            CUDA_CHECK_VOID(cudaFree(d_keep));
            CUDA_CHECK_VOID(cudaFree(d_num_out));
            if (idx.d_ref_non_unique_bm) {
                CUDA_CHECK_VOID(cudaFree(idx.d_ref_non_unique_bm));
                idx.d_ref_non_unique_bm = nullptr;
            }
        }
    } else {
        CUDA_CHECK_VOID(cudaFree(d_alt_vals));
        d_final_vals = d_sorted_vals;
        idx.n_kept = idx.n_valid;

        uint64_t equiv_kept = ((uint64_t)idx.n_kept * 256ULL) / (uint64_t)n_parts;
        while (bucket_bits < 26 && (1ULL << bucket_bits) * 4ULL < equiv_kept) bucket_bits++;
        uint32_t sub_bits = (uint32_t)(bucket_bits - 8);
        bucket_base = part_start << sub_bits;
        num_buckets = n_parts << sub_bits;
        shift_bits = 64 - bucket_bits;
        tag_shift = shift_bits - 32;

        CUDA_CHECK_VOID(cudaMalloc(&d_bucket_off, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMemset(d_bucket_off, 0, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_tags, (size_t)idx.n_kept * sizeof(uint32_t)));

        size_t scan_bytes = 0;
        cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);
        {
            CudaEventScope ev(prof, prof ? &prof->bucket_build_ms : nullptr);
            count_kept_per_bucket_kernel<uint64_t, false><<<(idx.n_valid + 255U) / 256U, 256>>>(
                d_sorted_keys, d_final_vals, nullptr, idx.n_valid,
                shift_bits, tag_shift, bucket_base, wl, eff_mask,
                d_bucket_off, d_tags, idx.d_ref_non_unique_bm);
            CUDA_CHECK_VOID(cudaGetLastError());
            cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_bucket_off, d_bucket_off, (int)(num_buckets + 1U));
        }
        CUDA_CHECK_VOID(cudaFree(d_sorted_keys));
    }

    idx.d_tags = d_tags;
    idx.d_vals = d_final_vals;
    idx.d_bucket_off = d_bucket_off;
    idx.num_buckets = num_buckets;
    idx.bucket_bits = bucket_bits;
    idx.shift_bits = shift_bits;
    idx.tag_shift = tag_shift;
    idx.bucket_base = bucket_base;
    idx.is_direct_csr = false;
    return idx;
}

void gpu_free_ref_index(GpuRefIndex &idx) {
    if (idx.d_tags) (void)cudaFree(idx.d_tags);
    if (idx.d_vals) (void)cudaFree(idx.d_vals);
    if (idx.d_bucket_off) (void)cudaFree(idx.d_bucket_off);
    if (idx.owns_dropped_bitmap && idx.d_dropped_bitmap) (void)cudaFree(idx.d_dropped_bitmap);
    if (idx.d_ref_non_unique_bm) (void)cudaFree(idx.d_ref_non_unique_bm);
    idx = {};
}
