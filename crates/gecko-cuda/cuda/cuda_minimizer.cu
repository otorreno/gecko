/*
 * cuda_minimizer.cu
 *
 * W3-F14: 100% Lossless GPU Canonical Minimizer Super-Seeding with High-Freq
 * Exact-k Overflow Routing (M_CAP), 48-bit Bijective mix48 Hash-CSR (0.000%
 * kp-mer tag collisions, bit-31 overflow sentinel), 0-stack CTA shared-memory
 * sliding-window extraction, CTA-privatized L2-partitioned query probe, and
 * 4-instruction 2-bit exact-k expansion + L/R-maximal filter.
 */

#include <algorithm>
#include <cstdlib>
#include <cub/cub.cuh>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include "cuda_common.cuh"

/* 512 core positions/wave x 8 waves = 4096 bp/CTA (GECKO_EXTRACT_CHUNK_BP).
 * Keeps s_scores (544 u64 = 4.25 KiB) + stage/hist (<= 6.1 KiB) well under 16 KiB
 * and avoids local arrays so every kernel achieves 0 B stack & 100% T4 occupancy. */
static constexpr uint32_t MIN_WAVE_BP       = 512U;
static constexpr uint32_t MIN_WAVES_PER_CTA = 8U;
static constexpr uint32_t MIN_SCORES_LEN    = 544U; /* >= 512 + 2*11 (max halo for wl<=32), 32-aligned */
static constexpr uint32_t MIN_L2_MAX_PARTS  = 1024U;

struct CastU8ToU32Op {
    __host__ __device__ __forceinline__ uint32_t operator()(uint8_t x) const {
        return (uint32_t)x;
    }
};

struct IsNonZeroU8ToU32Op {
    __host__ __device__ __forceinline__ uint32_t operator()(uint8_t x) const {
        return (x != 0U) ? 1U : 0U;
    }
};

/* ---- Device helper: canonical kp-mer score / (h48, strand, is_pal) ---- */

__device__ __forceinline__ uint64_t kp_score_at(
    const GpuSeqView &seq,
    uint32_t pos,
    int kp,
    uint32_t valid_mask_kp,
    uint64_t kp_mask,
    int rev_rshift_kp
) {
    if ((uint64_t)pos + (uint64_t)kp > seq.len) return ~0ULL;
    if (seq.invalid_bm && (load_32_1bit(seq.invalid_bm, pos) & valid_mask_kp) != 0U) {
        return ~0ULL;
    }
    uint64_t w_lsb = load_32_2bit(seq.bits, pos);
    uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift_kp;
    uint64_t kr = (~w_lsb) & kp_mask;
    uint64_t canon = (kf <= kr) ? kf : kr;
    return fmix64(canon);
}

__device__ __forceinline__ void kp_extract_canon(
    const GpuSeqView &seq,
    uint32_t pos,
    uint64_t kp_mask,
    int rev_rshift_kp,
    uint64_t &out_h48,
    uint32_t &out_strand,
    uint32_t &out_is_pal
) {
    uint64_t w_lsb = load_32_2bit(seq.bits, pos);
    uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift_kp;
    uint64_t kr = (~w_lsb) & kp_mask;
    out_strand = (kf <= kr) ? 0U : 1U;
    out_is_pal = (kf == kr) ? 1U : 0U;
    uint64_t canon = out_strand ? kr : kf;
    out_h48 = gecko_mix48(canon);
}

__device__ __forceinline__ void kmer_extract_canon_hash(
    const GpuSeqView &seq,
    uint32_t pos,
    uint64_t span_mask_wl,
    int rev_rshift_wl,
    uint64_t &out_hash,
    uint32_t &out_strand
) {
    uint64_t w_lsb = load_32_2bit(seq.bits, pos);
    uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift_wl;
    uint64_t kr = (~w_lsb) & span_mask_wl;
    out_strand = (kf <= kr) ? 0U : 1U;
    uint64_t canon = out_strand ? kr : kf;
    out_hash = fmix64(canon);
}

/* Scan left/right around ref position p_r directly in seqY (0 smem, 0 stack). */
__device__ __forceinline__ bool scan_ref_window_direct(
    const GpuSeqView &seqY,
    uint32_t p_r,
    int kp,
    int w,
    uint32_t valid_mask_kp,
    uint64_t kp_mask,
    int rev_rshift_kp,
    bool emit_ties,
    uint32_t &out_L,
    uint32_t &out_R
) {
    uint64_t s = kp_score_at(seqY, p_r, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
    if (s == ~0ULL) {
        out_L = 0U;
        out_R = 0U;
        return false;
    }
    uint32_t max_d = (uint32_t)(w - 1);
    uint32_t L = 0U;
    for (uint32_t d = 1U; d <= max_d; d++) {
        if (d > p_r) break;
        uint64_t ls = kp_score_at(seqY, p_r - d, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
        if (ls == ~0ULL || (emit_ties ? (ls < s) : (ls <= s))) break;
        L = d;
    }
    uint32_t R = 0U;
    for (uint32_t d = 1U; d <= max_d; d++) {
        uint64_t rs = kp_score_at(seqY, p_r + d, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
        if (rs == ~0ULL || rs < s) break;
        R = d;
    }
    out_L = L;
    out_R = R;
    return (L + R >= max_d);
}

/* ============================================================================
 * PART 1: GPU Ref Minimizer Index Build Kernels
 * ============================================================================ */

/* Step 1a: Count canonical minimizer emissions per 4096-bp chunk (EMIT_TIES = true). */
static __global__ __launch_bounds__(256, 4) void count_ref_minimizers_kernel(
    GpuSeqView seqY,
    int kp,
    int w,
    uint32_t * __restrict__ chunk_counts
) {
    __shared__ uint64_t s_scores[MIN_SCORES_LEN];
    __shared__ uint32_t s_warp_sums[8];

    uint32_t tid = threadIdx.x;
    uint32_t chunk_id = blockIdx.x;
    uint64_t chunk_start = (uint64_t)chunk_id * (uint64_t)GECKO_EXTRACT_CHUNK_BP;
    uint64_t total_kp = (seqY.len >= (uint64_t)kp) ? (seqY.len - (uint64_t)kp + 1ULL) : 0ULL;
    if (chunk_start >= total_kp) {
        if (tid == 0U) chunk_counts[chunk_id] = 0U;
        return;
    }

    uint32_t H = (uint32_t)(w - 1);
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;
    uint32_t local_cnt = 0U;

    for (uint32_t wave = 0U; wave < MIN_WAVES_PER_CTA; wave++) {
        uint64_t w_start = chunk_start + (uint64_t)wave * (uint64_t)MIN_WAVE_BP;
        uint32_t wave_len = 0U;
        if (w_start < total_kp) {
            uint64_t rem = total_kp - w_start;
            wave_len = (rem < (uint64_t)MIN_WAVE_BP) ? (uint32_t)rem : MIN_WAVE_BP;
        }
        uint32_t fill_len = (wave_len > 0U) ? (wave_len + 2U * H) : 0U;

        for (uint32_t i = tid; i < fill_len; i += 256U) {
            int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
            if (gp < 0 || (uint64_t)gp >= total_kp) {
                s_scores[i] = ~0ULL;
            } else {
                s_scores[i] = kp_score_at(seqY, (uint32_t)gp, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
            }
        }
        __syncthreads();

        uint32_t p0 = tid;
        if (p0 < wave_len) {
            int L, R;
            if (gecko_min_window_eval(s_scores, (int)(p0 + H), w, true, L, R)) local_cnt++;
        }
        uint32_t p1 = tid + 256U;
        if (p1 < wave_len) {
            int L, R;
            if (gecko_min_window_eval(s_scores, (int)(p1 + H), w, true, L, R)) local_cnt++;
        }
        __syncthreads();
    }

    uint32_t lane = tid & 31U;
    uint32_t warp = tid >> 5;
    uint32_t wsum = local_cnt;
    #pragma unroll
    for (int off = 16; off > 0; off >>= 1) {
        wsum += __shfl_down_sync(0xFFFFFFFFU, wsum, off);
    }
    if (lane == 0U) s_warp_sums[warp] = wsum;
    __syncthreads();

    if (tid == 0U) {
        uint32_t tot = 0U;
        #pragma unroll
        for (int i = 0; i < 8; i++) tot += s_warp_sums[i];
        chunk_counts[chunk_id] = tot;
    }
}

/* Step 1b: Extract ref minimizers (h48, (p_r << 1) | r_strand) with smem union staging. */
static __global__ __launch_bounds__(256, 4) void extract_ref_minimizers_kernel(
    GpuSeqView seqY,
    int kp,
    int w,
    const uint32_t * __restrict__ chunk_offsets,
    uint64_t * __restrict__ out_keys,
    uint32_t * __restrict__ out_vals
) {
    __shared__ union {
        uint64_t scores[MIN_SCORES_LEN];
        struct {
            uint64_t keys[MIN_WAVE_BP];
            uint32_t vals[MIN_WAVE_BP];
        } stage;
    } smem;
    __shared__ uint32_t s_warp_totals[8];
    __shared__ uint32_t s_wave_total;

    uint32_t tid = threadIdx.x;
    uint32_t lane = tid & 31U;
    uint32_t warp = tid >> 5;
    uint32_t chunk_id = blockIdx.x;
    uint64_t chunk_start = (uint64_t)chunk_id * (uint64_t)GECKO_EXTRACT_CHUNK_BP;
    uint64_t total_kp = (seqY.len >= (uint64_t)kp) ? (seqY.len - (uint64_t)kp + 1ULL) : 0ULL;
    if (chunk_start >= total_kp) return;

    uint32_t chunk_dst = chunk_offsets[chunk_id];
    uint32_t H = (uint32_t)(w - 1);
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;

    for (uint32_t wave = 0U; wave < MIN_WAVES_PER_CTA; wave++) {
        uint64_t w_start = chunk_start + (uint64_t)wave * (uint64_t)MIN_WAVE_BP;
        uint32_t wave_len = 0U;
        if (w_start < total_kp) {
            uint64_t rem = total_kp - w_start;
            wave_len = (rem < (uint64_t)MIN_WAVE_BP) ? (uint32_t)rem : MIN_WAVE_BP;
        }
        uint32_t fill_len = (wave_len > 0U) ? (wave_len + 2U * H) : 0U;

        for (uint32_t i = tid; i < fill_len; i += 256U) {
            int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
            if (gp < 0 || (uint64_t)gp >= total_kp) {
                smem.scores[i] = ~0ULL;
            } else {
                smem.scores[i] = kp_score_at(seqY, (uint32_t)gp, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
            }
        }
        __syncthreads();

        uint64_t lk0 = 0ULL, lk1 = 0ULL;
        uint32_t lv0 = 0U, lv1 = 0U;
        uint32_t ok0 = 0U, ok1 = 0U;

        uint32_t p0 = tid;
        if (p0 < wave_len) {
            int L, R;
            if (gecko_min_window_eval(smem.scores, (int)(p0 + H), w, true, L, R)) {
                uint32_t gpos = (uint32_t)(w_start + (uint64_t)p0);
                uint32_t st, pal;
                kp_extract_canon(seqY, gpos, kp_mask, rev_rshift_kp, lk0, st, pal);
                lv0 = (gpos << 1) | st;
                ok0 = 1U;
            }
        }
        uint32_t p1 = tid + 256U;
        if (p1 < wave_len) {
            int L, R;
            if (gecko_min_window_eval(smem.scores, (int)(p1 + H), w, true, L, R)) {
                uint32_t gpos = (uint32_t)(w_start + (uint64_t)p1);
                uint32_t st, pal;
                kp_extract_canon(seqY, gpos, kp_mask, rev_rshift_kp, lk1, st, pal);
                lv1 = (gpos << 1) | st;
                ok1 = 1U;
            }
        }
        __syncthreads(); /* Done reading smem.scores before reusing as smem.stage */

        uint32_t my_cnt = ok0 + ok1;
        uint32_t incl = my_cnt;
        #pragma unroll
        for (int off = 1; off < 32; off <<= 1) {
            uint32_t prev = __shfl_up_sync(0xFFFFFFFFU, incl, off);
            if (lane >= (uint32_t)off) incl += prev;
        }
        if (lane == 31U) s_warp_totals[warp] = incl;
        __syncthreads();

        if (tid == 0U) {
            uint32_t acc = 0U;
            #pragma unroll
            for (int i = 0; i < 8; i++) {
                uint32_t c = s_warp_totals[i];
                s_warp_totals[i] = acc;
                acc += c;
            }
            s_wave_total = acc;
        }
        __syncthreads();

        uint32_t wave_total = s_wave_total;
        if (wave_total > 0U) {
            uint32_t dst_s = s_warp_totals[warp] + (incl - my_cnt);
            if (ok0) {
                smem.stage.keys[dst_s] = lk0;
                smem.stage.vals[dst_s] = lv0;
                dst_s++;
            }
            if (ok1) {
                smem.stage.keys[dst_s] = lk1;
                smem.stage.vals[dst_s] = lv1;
            }
        }
        __syncthreads();

        for (uint32_t idx = tid; idx < wave_total; idx += 256U) {
            out_keys[chunk_dst + idx] = smem.stage.keys[idx];
            out_vals[chunk_dst + idx] = smem.stage.vals[idx];
        }
        chunk_dst += wave_total;
        __syncthreads();
    }
}

/* Step 3a: O(log(freq_limit)) parallel HF run detector (zero serial homopolymer bottleneck). */
static __global__ void mark_hf_entries_kernel(
    const uint64_t * __restrict__ sorted_keys,
    uint32_t n,
    uint32_t freq_limit,
    uint8_t * __restrict__ is_hf
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    uint64_t k = sorted_keys[i];
    bool same_prev = (i > 0U && sorted_keys[i - 1U] == k);
    bool same_next = (i + 1U < n && sorted_keys[i + 1U] == k);
    if (!same_prev && !same_next) {
        is_hf[i] = 0U;
        return;
    }
    if ((i >= freq_limit && sorted_keys[i - freq_limit] == k) ||
        ((uint64_t)i + (uint64_t)freq_limit < (uint64_t)n && sorted_keys[i + freq_limit] == k)) {
        is_hf[i] = 1U;
        return;
    }
    uint32_t lo = (i > freq_limit) ? (i - freq_limit) : 0U;
    uint32_t hi = i;
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1);
        if (sorted_keys[mid] < k) lo = mid + 1U;
        else hi = mid;
    }
    uint32_t run_start = lo;
    is_hf[i] = ((uint64_t)run_start + (uint64_t)freq_limit < (uint64_t)n &&
                sorted_keys[run_start + freq_limit] == k) ? 1U : 0U;
}

/* Step 3b-1: Count exact k-mers owned by each HF minimizer entry (strict L_strict / R_ge partition). */
static __global__ void count_hf_cand_kmers_kernel(
    GpuSeqView seqY,
    const uint32_t * __restrict__ sorted_vals,
    const uint32_t * __restrict__ hf_idx,
    uint32_t n_hf,
    int kp,
    int w,
    uint32_t * __restrict__ out_r_lo,
    uint8_t  * __restrict__ out_r_cnt
) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_hf) return;

    uint32_t p_r = sorted_vals[hf_idx[t]] >> 1;
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;

    uint32_t L = 0U, R = 0U;
    if (scan_ref_window_direct(seqY, p_r, kp, w, valid_mask_kp, kp_mask, rev_rshift_kp, false, L, R)) {
        uint32_t r_lo = p_r - L;
        uint32_t cnt = L + R - (uint32_t)(w - 1) + 1U;
        out_r_lo[t] = r_lo;
        out_r_cnt[t] = (uint8_t)cnt;
    } else {
        out_r_lo[t] = 0U;
        out_r_cnt[t] = 0U;
    }
}

/* Step 3b-2: Extract candidate exact k-mers (fmix64(k_canon), (r << 1) | k_strand) for HF entries. */
static __global__ void extract_hf_cand_kmers_kernel(
    GpuSeqView seqY,
    const uint32_t * __restrict__ hf_r_lo,
    const uint8_t  * __restrict__ hf_r_cnt,
    const uint32_t * __restrict__ hf_off,
    uint32_t n_hf,
    int wl,
    uint64_t * __restrict__ out_keys,
    uint32_t * __restrict__ out_vals
) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_hf) return;

    uint32_t cnt = (uint32_t)hf_r_cnt[t];
    if (cnt == 0U) return;

    uint32_t r_lo = hf_r_lo[t];
    uint32_t dst = hf_off[t];
    uint64_t span_mask_wl = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rev_rshift_wl = 64 - 2 * wl;

    for (uint32_t d = 0U; d < cnt; d++) {
        uint32_t r = r_lo + d;
        uint64_t h;
        uint32_t st;
        kmer_extract_canon_hash(seqY, r, span_mask_wl, rev_rshift_wl, h, st);
        out_keys[dst + d] = h;
        out_vals[dst + d] = (r << 1) | st;
    }
}

/* Step 3c: Exact per-strand freq_limit marking into 1-bit d_dropped_bitmap (word-coalesced REDG). */
static __global__ void mark_hf_kmer_dropped_kernel(
    const uint64_t * __restrict__ keys,
    const uint32_t * __restrict__ vals,
    uint32_t n,
    uint32_t freq_limit,
    uint32_t * __restrict__ dropped_bitmap,
    uint32_t * __restrict__ any_dropped
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    uint64_t k = keys[i];
    if (i > 0U && keys[i - 1U] == k) return;
    if (i + 1U >= n || keys[i + 1U] != k) return;

    uint32_t step = 1U;
    uint32_t rem = n - i;
    while (step < rem && keys[i + step] == k) {
        step <<= 1;
    }
    uint32_t lo = i + (step >> 1);
    uint32_t hi = (step >= rem) ? n : (i + step);
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1);
        if (keys[mid] == k) lo = mid + 1U;
        else hi = mid;
    }
    uint32_t run_end = lo;
    uint32_t run_len = run_end - i;
    if (run_len <= freq_limit) return;

    uint32_t c0 = 0U, c1 = 0U;
    for (uint32_t m = i; m < run_end; m++) {
        if (vals[m] & 1U) c1++;
        else c0++;
        if (c0 > freq_limit && c1 > freq_limit) break;
    }
    uint8_t k0 = (c0 <= freq_limit) ? 1U : 0U;
    uint8_t k1 = (c1 <= freq_limit) ? 1U : 0U;
    if (k0 && k1) return;

    if (*any_dropped == 0U) *any_dropped = 1U;
    uint32_t cur_w = 0xFFFFFFFFU;
    uint32_t cur_mask = 0U;
    for (uint32_t m = i; m < run_end; m++) {
        uint32_t v = vals[m];
        uint8_t keep = (v & 1U) ? k1 : k0;
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
}

/* Step 3d-1: Evaluate HF minimizer survival using L_ge / R_ge (REQUIRED for reverse-strand tie symmetry!). */
static __global__ void eval_hf_min_survival_kernel(
    GpuSeqView seqY,
    const uint32_t * __restrict__ sorted_vals,
    const uint32_t * __restrict__ hf_idx,
    uint32_t n_hf,
    int kp,
    int w,
    const uint32_t * __restrict__ dropped_bitmap,
    uint8_t * __restrict__ min_keep
) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_hf) return;

    uint32_t idx = hf_idx[t];
    uint32_t p_r = sorted_vals[idx] >> 1;
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;

    uint32_t L = 0U, R = 0U;
    if (!scan_ref_window_direct(seqY, p_r, kp, w, valid_mask_kp, kp_mask, rev_rshift_kp, true, L, R)) {
        min_keep[idx] = 0U;
        return;
    }
    uint32_t r_lo = p_r - L;
    uint32_t cnt = L + R - (uint32_t)(w - 1) + 1U;
    uint32_t db = load_32_1bit(dropped_bitmap, r_lo);
    uint32_t mask = (1U << cnt) - 1U;
    min_keep[idx] = ((db & mask) != mask) ? 1U : 0U;
}

/* Step 3d-2: O(log(run_len)) parallel M_CAP overflow routing classification.
 * Routes any surviving minimizer bucket with kept_cnt > m_cap to a single bit-31 sentinel
 * entry (min_keep[run_start] = 2) and marks all HF entries in that run for exact-k overflow CSR. */
static __global__ void classify_overflow_entries_kernel(
    const uint64_t * __restrict__ sorted_keys,
    const uint32_t * __restrict__ keep_prefix,
    const uint32_t * __restrict__ hf_idx,
    uint32_t n_raw_min,
    uint32_t n_hf,
    uint32_t m_cap,
    uint8_t  * __restrict__ min_keep,
    uint8_t  * __restrict__ hf_is_ovf,
    uint32_t * __restrict__ ovf_key_cnt
) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_hf) return;

    uint32_t idx = hf_idx[t];
    uint64_t k = sorted_keys[idx];

    uint32_t run_start = idx;
    if (idx > 0U && sorted_keys[idx - 1U] == k) {
        uint32_t step = 1U;
        while (step < idx && sorted_keys[idx - step] == k) {
            step <<= 1;
        }
        uint32_t lo = (step >= idx) ? 0U : (idx - step);
        uint32_t hi = idx - (step >> 1);
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (sorted_keys[mid] < k) lo = mid + 1U;
            else hi = mid;
        }
        run_start = lo;
    }

    uint32_t run_end = idx + 1U;
    if (run_end < n_raw_min && sorted_keys[run_end] == k) {
        uint32_t step = 1U;
        uint32_t rem = n_raw_min - idx;
        while (step < rem && sorted_keys[idx + step] == k) {
            step <<= 1;
        }
        uint32_t lo = idx + (step >> 1);
        uint32_t hi = (step >= rem) ? n_raw_min : (idx + step);
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (sorted_keys[mid] == k) lo = mid + 1U;
            else hi = mid;
        }
        run_end = lo;
    }

    uint32_t kept_cnt = keep_prefix[run_end] - keep_prefix[run_start];
    if (kept_cnt > m_cap) {
        hf_is_ovf[t] = 1U;
        if (idx == run_start) {
            min_keep[idx] = 2U;
            atomicAdd(ovf_key_cnt, 1U);
        } else {
            min_keep[idx] = 0U;
        }
    }
}

/* Step 3e-1: Count non-dropped exact k-mers to route into the tiny exact-k overflow CSR. */
static __global__ void count_ovf_kmers_kernel(
    const uint8_t  * __restrict__ hf_is_ovf,
    const uint32_t * __restrict__ hf_r_lo,
    const uint8_t  * __restrict__ hf_r_cnt,
    uint32_t n_hf,
    const uint32_t * __restrict__ dropped_bitmap,
    uint8_t * __restrict__ out_ovf_cnt
) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_hf) return;

    if (hf_is_ovf[t] == 0U) {
        out_ovf_cnt[t] = 0U;
        return;
    }
    uint32_t cnt = (uint32_t)hf_r_cnt[t];
    if (cnt == 0U) {
        out_ovf_cnt[t] = 0U;
        return;
    }
    if (!dropped_bitmap) {
        out_ovf_cnt[t] = (uint8_t)cnt;
        return;
    }
    uint32_t r_lo = hf_r_lo[t];
    uint32_t db = load_32_1bit(dropped_bitmap, r_lo);
    uint32_t mask = (1U << cnt) - 1U;
    out_ovf_cnt[t] = (uint8_t)__popc((~db) & mask);
}

/* Step 3e-2: Extract non-dropped exact k-mers into overflow CSR sort buffers. */
static __global__ void extract_ovf_kmers_kernel(
    GpuSeqView seqY,
    const uint8_t  * __restrict__ hf_is_ovf,
    const uint32_t * __restrict__ hf_r_lo,
    const uint8_t  * __restrict__ hf_r_cnt,
    const uint32_t * __restrict__ ovf_off,
    uint32_t n_hf,
    int wl,
    const uint32_t * __restrict__ dropped_bitmap,
    uint64_t * __restrict__ out_keys,
    uint32_t * __restrict__ out_vals
) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_hf) return;
    if (hf_is_ovf[t] == 0U) return;

    uint32_t cnt = (uint32_t)hf_r_cnt[t];
    if (cnt == 0U) return;

    uint32_t r_lo = hf_r_lo[t];
    uint32_t dst = ovf_off[t];
    uint32_t db = dropped_bitmap ? load_32_1bit(dropped_bitmap, r_lo) : 0U;
    uint32_t mask = (1U << cnt) - 1U;
    uint32_t keep_mask = (~db) & mask;
    if (keep_mask == 0U) return;

    uint64_t span_mask_wl = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rev_rshift_wl = 64 - 2 * wl;

    while (keep_mask != 0U) {
        uint32_t d = (uint32_t)(__ffs((int)keep_mask) - 1);
        keep_mask &= keep_mask - 1U;
        uint32_t r = r_lo + d;
        uint64_t h;
        uint32_t st;
        kmer_extract_canon_hash(seqY, r, span_mask_wl, rev_rshift_wl, h, st);
        out_keys[dst] = h;
        out_vals[dst] = (r << 1) | st;
        dst++;
    }
}

/* Step 3e-3: Build overflow exact-k CSR tags & warp-aggregated bucket counts. */
static __global__ void scatter_ovf_csr_kernel(
    const uint64_t * __restrict__ sorted_ovf_keys,
    uint32_t n,
    int ovf_shift_bits,
    int ovf_tag_shift,
    uint32_t * __restrict__ out_ovf_tags,
    uint32_t * __restrict__ ovf_bucket_counts
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    uint64_t k = sorted_ovf_keys[i];
    out_ovf_tags[i] = (uint32_t)(k >> ovf_tag_shift);
    uint32_t b = (uint32_t)(k >> ovf_shift_bits);

    uint32_t active = __activemask();
    uint32_t same = __match_any_sync(active, b);
    int leader = __ffs((int)same) - 1;
    if ((int)(threadIdx.x & 31U) == leader) {
        atomicAdd(&ovf_bucket_counts[b], (uint32_t)__popc(same));
    }
}

/* Step 4: Fused single-pass main minimizer Hash-CSR scatter (tags + vals + warp-aggregated counts).
 * Sets bit 31 (GECKO_MIN_OVF_SENTINEL) on d_tags and 0xFFFFFFFFU on d_vals for overflow sentinels. */
static __global__ void scatter_min_csr_kernel(
    const uint64_t * __restrict__ sorted_keys,
    const uint32_t * __restrict__ sorted_vals,
    const uint8_t  * __restrict__ min_keep,
    const uint32_t * __restrict__ out_pos,
    uint32_t n_raw_min,
    int shift_bits,
    uint32_t tag_mask,
    uint32_t * __restrict__ out_tags,
    uint32_t * __restrict__ out_vals,
    uint32_t * __restrict__ bucket_counts
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_raw_min) return;

    uint8_t act = min_keep[i];
    if (act == 0U) return;

    uint32_t dst = out_pos[i];
    uint64_t h48 = sorted_keys[i];
    uint32_t bkt = (uint32_t)(h48 >> shift_bits);
    uint32_t tag = ((uint32_t)h48) & tag_mask;

    if (act == 2U) {
        out_tags[dst] = tag | GECKO_MIN_OVF_SENTINEL;
        out_vals[dst] = GECKO_MIN_OVF_VAL_SENTINEL;
    } else {
        out_tags[dst] = tag;
        out_vals[dst] = sorted_vals[i];
    }

    uint32_t active = __activemask();
    uint32_t same = __match_any_sync(active, bkt);
    int leader = __ffs((int)same) - 1;
    if ((int)(threadIdx.x & 31U) == leader) {
        atomicAdd(&bucket_counts[bkt], (uint32_t)__popc(same));
    }
}

/* ============================================================================
 * PART 2: GPU Minimizer Query Kernels (Monolithic + L2-Partitioned + Overflow)
 * ============================================================================ */

__device__ __forceinline__ void stage_or_emit_hit(
    uint64_t key,
    uint64_t *s_hits,
    uint32_t *s_n,
    uint64_t *out_hits,
    uint32_t *hit_count,
    uint32_t max_hits
) {
    uint32_t slot = atomicAdd(s_n, 1U);
    if (slot < GECKO_HIT_STAGE_CAP) {
        s_hits[slot] = key;
    } else {
        uint32_t gslot = atomicAdd(hit_count, 1U);
        if (gslot < max_hits) out_hits[gslot] = key;
    }
}

/* Expand a matched canonical kp-mer (p_q, p_r) using 2-bit word ops and emit exact L/R-maximal
 * hits within [r_lo .. r_hi] (respecting dropped_bitmap gaps on BOTH forward and reverse strands). */
template <bool IS_REV>
__device__ __forceinline__ void min_expand_and_emit(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const uint32_t * __restrict__ dropped_bitmap,
    int wl,
    int kp,
    int w,
    uint32_t p_q,
    uint32_t p_r,
    uint32_t r_lo,
    uint32_t r_hi,
    HitKeyLayout key_layout,
    uint64_t *s_hits,
    uint32_t *s_n,
    uint64_t *out_hits,
    uint32_t *hit_count,
    uint32_t max_hits
) {
    uint32_t need_left = p_q - r_lo;
    uint32_t need_right = (r_hi + (uint32_t)wl) - (p_q + (uint32_t)kp);

    uint32_t ext_left = 0U, ext_right = 0U;
    gecko_min_expand_extents<IS_REV>(
        seqX.bits, seqX.invalid_bm, seqX.len,
        seqY.bits, seqY.invalid_bm, seqY.len,
        kp, p_q, p_r, need_left, need_right, ext_left, ext_right);

    if (ext_left + (uint32_t)kp + ext_right < (uint32_t)wl) return;

    uint32_t m_lo = p_q - ext_left;
    uint32_t m_hi = p_q + (uint32_t)kp + ext_right - (uint32_t)wl;
    uint32_t e_lo = (m_lo > r_lo) ? m_lo : r_lo;
    uint32_t e_hi = (m_hi < r_hi) ? m_hi : r_hi;
    if (e_lo > e_hi) return;

    /* Map [e_lo .. e_hi] to original seqY start coordinates [y_min .. y_max].
     * For IS_REV == true: y_orig(posX) = p_r + (p_q - posX) - (w - 1), strictly decreasing
     * with posX, so y_min = y_orig(e_hi) and y_max = y_orig(e_lo). */
    uint32_t y_min = (!IS_REV)
        ? (p_r - (p_q - e_lo))
        : (p_r + (p_q - e_hi) - (uint32_t)(w - 1));
    uint32_t y_max = y_min + (e_hi - e_lo);
    uint32_t base_y = (y_min > 0U) ? (y_min - 1U) : 0U;
    uint32_t dbits = dropped_bitmap ? load_32_1bit(dropped_bitmap, base_y) : 0U;
    uint32_t dmask = (1U << (y_max + 2U - base_y)) - 1U;

    if ((dbits & dmask) == 0U) {
        if (e_lo == m_lo) {
            uint32_t ry = (!IS_REV)
                ? (p_r - (p_q - e_lo))
                : (p_r + (p_q - e_lo) - (uint32_t)(w - 1));
            uint32_t out_y = (!IS_REV) ? ry : (uint32_t)(seqY.len - (uint64_t)ry - (uint64_t)wl);
            stage_or_emit_hit(
                pack_hit_key_tile(key_layout, IS_REV ? 1U : 0U, e_lo, out_y),
                s_hits, s_n, out_hits, hit_count, max_hits);
        }
        if (e_hi == m_hi && (e_hi != e_lo || e_lo != m_lo)) {
            uint32_t ry = (!IS_REV)
                ? (p_r - (p_q - e_hi))
                : (p_r + (p_q - e_hi) - (uint32_t)(w - 1));
            uint32_t out_y = (!IS_REV) ? ry : (uint32_t)(seqY.len - (uint64_t)ry - (uint64_t)wl);
            stage_or_emit_hit(
                pack_hit_key_tile(key_layout, IS_REV ? 1U : 0U, e_hi, out_y),
                s_hits, s_n, out_hits, hit_count, max_hits);
        }
    } else {
        for (uint32_t posX = e_lo; posX <= e_hi; posX++) {
            uint32_t ry = (!IS_REV)
                ? (p_r - (p_q - posX))
                : (p_r + (p_q - posX) - (uint32_t)(w - 1));
            if (((dbits >> (ry - base_y)) & 1U) != 0U) continue;
            bool left_ext = (posX > m_lo);
            if (left_ext) {
                /* IS_REV: stepping posX left (posX - 1) steps y_orig RIGHT (ry + 1) */
                uint32_t ly = (!IS_REV) ? (ry - 1U) : (ry + 1U);
                left_ext = (((dbits >> (ly - base_y)) & 1U) == 0U);
            }
            bool right_ext = (posX < m_hi);
            if (right_ext) {
                /* IS_REV: stepping posX right (posX + 1) steps y_orig LEFT (ry - 1) */
                uint32_t rry = (!IS_REV) ? (ry + 1U) : (ry - 1U);
                right_ext = (((dbits >> (rry - base_y)) & 1U) == 0U);
            }
            if (!left_ext || !right_ext) {
                uint32_t out_y = (!IS_REV) ? ry : (uint32_t)(seqY.len - (uint64_t)ry - (uint64_t)wl);
                stage_or_emit_hit(
                    pack_hit_key_tile(key_layout, IS_REV ? 1U : 0U, posX, out_y),
                    s_hits, s_n, out_hits, hit_count, max_hits);
            }
        }
    }
}

/* Exact-k Overflow CSR probe helper for deferred overflow queue items (0-stack). */
static __global__ __launch_bounds__(256, 4) void min_probe_overflow_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ ovf_bucket_off,
    const uint32_t * __restrict__ ovf_tags,
    const uint32_t * __restrict__ ovf_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int ovf_shift_bits,
    int ovf_tag_shift,
    int wl,
    bool compare_reverse,
    HitKeyLayout key_layout,
    const uint32_t * __restrict__ ovf_q_r_lo,
    const uint8_t  * __restrict__ ovf_q_span,
    uint32_t n_ovf_q,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits
) {
    __shared__ uint64_t s_hits[GECKO_HIT_STAGE_CAP];
    __shared__ uint32_t s_n;
    __shared__ uint32_t s_base;

    if (threadIdx.x == 0U) s_n = 0U;
    __syncthreads();

    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n_ovf_q) {
        uint32_t r_lo = ovf_q_r_lo[idx];
        uint32_t span = (uint32_t)ovf_q_span[idx];
        uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
        int rev_rshift = 64 - 2 * wl;

        for (uint32_t d = 0U; d <= span; d++) {
            uint32_t posX = r_lo + d;
            uint64_t w_lsb = load_32_2bit(seqX.bits, posX);
            uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift;
            uint64_t kr = (~w_lsb) & span_mask;
            uint32_t x_strand = (kf <= kr) ? 0U : 1U;
            bool is_pal = (kf == kr);
            uint64_t canon = x_strand ? kr : kf;
            uint64_t h = fmix64(canon);

            uint32_t ob = (uint32_t)(h >> ovf_shift_bits);
            uint32_t otag = (uint32_t)(h >> ovf_tag_shift);
            uint32_t lo = __ldg(&ovf_bucket_off[ob]);
            uint32_t hi = __ldg(&ovf_bucket_off[ob + 1U]);
            if (lo >= hi) continue;

            if (hi - lo > 8U) {
                uint32_t l = lo, r = hi;
                while (l < r) {
                    uint32_t mid = l + ((r - l) >> 1);
                    if (__ldg(&ovf_tags[mid]) < otag) l = mid + 1U;
                    else r = mid;
                }
                lo = l;
            }

            uint64_t x_seed_lsb = w_lsb & span_mask;
            bool prev_x_valid = (posX > 0U && !is_bit_set(seqX.invalid_bm, posX - 1U));
            uint32_t prev_xb = prev_x_valid ? base_2bit_dev(seqX.bits, posX - 1U) : 0U;
            bool next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len &&
                                 !is_bit_set(seqX.invalid_bm, posX + (uint32_t)wl));
            uint32_t next_xb = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 0U;

            for (uint32_t j = lo; j < hi; j++) {
                uint32_t t = __ldg(&ovf_tags[j]);
                if (t > otag) break;
                if (t < otag) continue;

                uint32_t ps = __ldg(&ovf_vals[j]);
                uint32_t posY = ps >> 1;
                uint32_t y_strand = ps & 1U;
                uint64_t y_word = load_32_2bit(seqY.bits, posY);

                if ((is_pal || x_strand == y_strand) &&
                    (((y_word ^ x_seed_lsb) & span_mask) == 0ULL)) {
                    bool left_ext = prev_x_valid && posY > 0U &&
                                    !is_bit_set(seqY.invalid_bm, posY - 1U) &&
                                    base_2bit_dev(seqY.bits, posY - 1U) == prev_xb &&
                                    !is_bit_set(dropped_bitmap, posY - 1U);
                    bool emit = !left_ext;
                    if (!emit) {
                        bool right_ext = next_x_valid &&
                                         ((uint64_t)posY + (uint64_t)wl < seqY.len) &&
                                         !is_bit_set(seqY.invalid_bm, posY + (uint32_t)wl) &&
                                         base_2bit_dev(seqY.bits, posY + (uint32_t)wl) == next_xb &&
                                         !is_bit_set(dropped_bitmap, posY + 1U);
                        emit = !right_ext;
                    }
                    if (emit) {
                        stage_or_emit_hit(
                            pack_hit_key_tile(key_layout, 0U, posX, posY),
                            s_hits, &s_n, out_hits, hit_count, max_hits);
                    }
                }

                if (compare_reverse && (is_pal || x_strand != y_strand)) {
                    uint64_t y_rc_lsb = (~reverse_32_2bit_dev(y_word)) >> rev_rshift;
                    if (((y_rc_lsb ^ x_seed_lsb) & span_mask) == 0ULL) {
                        uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)posY - (uint64_t)wl);
                        bool left_ext = prev_x_valid && posY_rc > 0U &&
                                        !is_bit_set(seqY.invalid_bm, posY + (uint32_t)wl) &&
                                        (base_2bit_dev(seqY.bits, posY + (uint32_t)wl) ^ 3U) == prev_xb &&
                                        !is_bit_set(dropped_bitmap, posY + 1U);
                        bool emit = !left_ext;
                        if (!emit) {
                            bool right_ext = next_x_valid && posY > 0U &&
                                             !is_bit_set(seqY.invalid_bm, posY - 1U) &&
                                             (base_2bit_dev(seqY.bits, posY - 1U) ^ 3U) == next_xb &&
                                             !is_bit_set(dropped_bitmap, posY - 1U);
                            emit = !right_ext;
                        }
                        if (emit) {
                            stage_or_emit_hit(
                                pack_hit_key_tile(key_layout, 1U, posX, posY_rc),
                                s_hits, &s_n, out_hits, hit_count, max_hits);
                        }
                    }
                }
            }
        }
    }

    __syncthreads();
    uint32_t flush_cnt = s_n;
    if (flush_cnt > GECKO_HIT_STAGE_CAP) flush_cnt = GECKO_HIT_STAGE_CAP;
    if (flush_cnt > 0U) {
        if (threadIdx.x == 0U) {
            s_base = atomicAdd(hit_count, flush_cnt);
        }
        __syncthreads();
        uint32_t base = s_base;
        if (base < max_hits) {
            uint32_t rem = max_hits - base;
            uint32_t n_write = (flush_cnt < rem) ? flush_cnt : rem;
            for (uint32_t k = threadIdx.x; k < n_write; k += 256U) {
                out_hits[base + k] = s_hits[k];
            }
        }
    }
}

/* Probe one query minimizer against the main 48-bit CSR; if it hits a bit-31 overflow sentinel,
 * enqueue (r_lo, span) into the deferred overflow queue (keeps probe kernel <= 56 regs & 0 stack!). */
template <bool HAS_OVF>
__device__ __forceinline__ void min_probe_main_entry(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const uint32_t * __restrict__ idx_bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int shift_bits,
    uint32_t tag_mask,
    int wl,
    int kp,
    int w,
    bool compare_reverse,
    HitKeyLayout key_layout,
    uint64_t h48,
    uint32_t p_q,
    uint32_t r_lo,
    uint32_t r_hi,
    uint32_t q_strand,
    uint32_t q_is_pal,
    uint64_t *s_hits,
    uint32_t *s_n,
    uint64_t *out_hits,
    uint32_t *hit_count,
    uint32_t max_hits,
    uint32_t *ovf_q_r_lo,
    uint8_t  *ovf_q_span,
    uint32_t *ovf_q_count,
    uint32_t ovf_q_cap
) {
    uint32_t bkt = (uint32_t)(h48 >> shift_bits);
    uint32_t tag = ((uint32_t)h48) & tag_mask;
    uint32_t lo = __ldg(&idx_bucket_off[bkt]);
    uint32_t hi = __ldg(&idx_bucket_off[bkt + 1U]);
    if (lo >= hi) return;

    if (hi - lo > 8U) {
        uint32_t l = lo, r = hi;
        while (l < r) {
            uint32_t mid = l + ((r - l) >> 1);
            uint32_t mt = __ldg(&idx_tags[mid]);
            if constexpr (HAS_OVF) mt &= GECKO_MIN_TAG_MASK;
            if (mt < tag) l = mid + 1U;
            else r = mid;
        }
        lo = l;
    }

    for (uint32_t j = lo; j < hi; j++) {
        uint32_t raw_t = __ldg(&idx_tags[j]);
        uint32_t t = HAS_OVF ? (raw_t & GECKO_MIN_TAG_MASK) : raw_t;
        if (t > tag) break;
        if (t < tag) continue;

        if constexpr (HAS_OVF) {
            if ((raw_t & GECKO_MIN_OVF_SENTINEL) != 0U) {
                uint32_t qslot = atomicAdd(ovf_q_count, 1U);
                if (qslot < ovf_q_cap) {
                    ovf_q_r_lo[qslot] = r_lo;
                    ovf_q_span[qslot] = (uint8_t)(r_hi - r_lo);
                }
                break;
            }
        }

        uint32_t ps = __ldg(&idx_vals[j]);
        uint32_t p_r = ps >> 1;
        uint32_t r_strand = ps & 1U;

        if (q_is_pal || q_strand == r_strand) {
            min_expand_and_emit<false>(
                seqX, seqY, dropped_bitmap, wl, kp, w,
                p_q, p_r, r_lo, r_hi, key_layout,
                s_hits, s_n, out_hits, hit_count, max_hits);
        }
        if (compare_reverse && (q_is_pal || q_strand != r_strand)) {
            min_expand_and_emit<true>(
                seqX, seqY, dropped_bitmap, wl, kp, w,
                p_q, p_r, r_lo, r_hi, key_layout,
                s_hits, s_n, out_hits, hit_count, max_hits);
        }
    }
}

/* Query Mode 1: Monolithic minimizer query kernel (small/L2-resident index or fallback).
 * 8.27 KiB smem (s_scores[544] + s_hits[512]), 0 B stack, 100% T4 occupancy. */
template <bool HAS_OVF>
static __global__ __launch_bounds__(256, 4) void min_query_monolithic_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ idx_bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int shift_bits,
    uint32_t tag_mask,
    int wl,
    int kp,
    int w,
    bool compare_reverse,
    uint32_t tile_min_x,
    uint32_t tile_max_x,
    uint32_t p_begin,
    uint32_t p_end,
    HitKeyLayout key_layout,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits,
    uint32_t * __restrict__ ovf_q_r_lo,
    uint8_t  * __restrict__ ovf_q_span,
    uint32_t * __restrict__ ovf_q_count,
    uint32_t ovf_q_cap
) {
    __shared__ uint64_t s_scores[MIN_SCORES_LEN];
    __shared__ uint64_t s_hits[GECKO_HIT_STAGE_CAP];
    __shared__ uint32_t s_n;
    __shared__ uint32_t s_base;

    uint32_t tid = threadIdx.x;
    uint32_t chunk_start = p_begin + blockIdx.x * GECKO_EXTRACT_CHUNK_BP;
    if (chunk_start >= p_end) return;

    uint32_t H = (uint32_t)(w - 1);
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;
    uint64_t total_kp = (seqX.len >= (uint64_t)kp) ? (seqX.len - (uint64_t)kp + 1ULL) : 0ULL;

    for (uint32_t wave = 0U; wave < MIN_WAVES_PER_CTA; wave++) {
        uint32_t w_start = chunk_start + wave * MIN_WAVE_BP;
        uint32_t wave_len = 0U;
        if (w_start < p_end) {
            uint32_t rem = p_end - w_start;
            wave_len = (rem < MIN_WAVE_BP) ? rem : MIN_WAVE_BP;
        }
        uint32_t fill_len = (wave_len > 0U) ? (wave_len + 2U * H) : 0U;

        if (tid == 0U) s_n = 0U;
        for (uint32_t i = tid; i < fill_len; i += 256U) {
            int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
            if (gp < 0 || (uint64_t)gp >= total_kp) {
                s_scores[i] = ~0ULL;
            } else {
                s_scores[i] = kp_score_at(seqX, (uint32_t)gp, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
            }
        }
        __syncthreads();

        #pragma unroll
        for (uint32_t step = 0U; step < 2U; step++) {
            uint32_t local_i = tid + step * 256U;
            if (local_i < wave_len) {
                int L, R;
                if (gecko_min_window_eval(s_scores, (int)(local_i + H), w, false, L, R)) {
                    uint32_t p_q = w_start + local_i;
                    uint32_t r_lo = p_q - (uint32_t)L;
                    uint32_t r_hi = p_q + (uint32_t)R - H;
                    if (r_lo < tile_min_x) r_lo = tile_min_x;
                    if (r_hi > tile_max_x) r_hi = tile_max_x;
                    if (r_lo <= r_hi) {
                        uint64_t h48;
                        uint32_t q_strand, q_is_pal;
                        kp_extract_canon(seqX, p_q, kp_mask, rev_rshift_kp, h48, q_strand, q_is_pal);
                        min_probe_main_entry<HAS_OVF>(
                            seqX, seqY, idx_bucket_off, idx_tags, idx_vals, dropped_bitmap,
                            shift_bits, tag_mask, wl, kp, w, compare_reverse, key_layout,
                            h48, p_q, r_lo, r_hi, q_strand, q_is_pal,
                            s_hits, &s_n, out_hits, hit_count, max_hits,
                            ovf_q_r_lo, ovf_q_span, ovf_q_count, ovf_q_cap);
                    }
                }
            }
        }
        __syncthreads();

        uint32_t flush_cnt = s_n;
        if (flush_cnt > GECKO_HIT_STAGE_CAP) flush_cnt = GECKO_HIT_STAGE_CAP;
        if (flush_cnt > 0U) {
            if (tid == 0U) s_base = atomicAdd(hit_count, flush_cnt);
            __syncthreads();
            uint32_t base = s_base;
            if (base < max_hits) {
                uint32_t rem = max_hits - base;
                uint32_t n_write = (flush_cnt < rem) ? flush_cnt : rem;
                for (uint32_t k = tid; k < n_write; k += 256U) {
                    out_hits[base + k] = s_hits[k];
                }
            }
        }
        __syncthreads();
    }
}

/* Query Mode 2 (L2-Partitioned) — Stage 1: CTA-privatized partition histogram (ATOMS). */
static __global__ __launch_bounds__(256, 4) void min_l2_count_parts_kernel(
    GpuSeqView seqX,
    int shift_bits,
    uint32_t P,
    int part_shift,
    int kp,
    int w,
    uint32_t tile_min_x,
    uint32_t tile_max_x,
    uint32_t p_begin,
    uint32_t p_end,
    uint32_t * __restrict__ part_counts
) {
    __shared__ uint64_t s_scores[MIN_SCORES_LEN];
    __shared__ uint32_t s_hist[MIN_L2_MAX_PARTS];

    uint32_t tid = threadIdx.x;
    for (uint32_t p = tid; p < P; p += 256U) s_hist[p] = 0U;
    __syncthreads();

    uint32_t chunk_start = p_begin + blockIdx.x * GECKO_EXTRACT_CHUNK_BP;
    uint32_t H = (uint32_t)(w - 1);
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;
    uint64_t total_kp = (seqX.len >= (uint64_t)kp) ? (seqX.len - (uint64_t)kp + 1ULL) : 0ULL;

    if (chunk_start < p_end) {
        for (uint32_t wave = 0U; wave < MIN_WAVES_PER_CTA; wave++) {
            uint32_t w_start = chunk_start + wave * MIN_WAVE_BP;
            uint32_t wave_len = 0U;
            if (w_start < p_end) {
                uint32_t rem = p_end - w_start;
                wave_len = (rem < MIN_WAVE_BP) ? rem : MIN_WAVE_BP;
            }
            uint32_t fill_len = (wave_len > 0U) ? (wave_len + 2U * H) : 0U;

            for (uint32_t i = tid; i < fill_len; i += 256U) {
                int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
                if (gp < 0 || (uint64_t)gp >= total_kp) {
                    s_scores[i] = ~0ULL;
                } else {
                    s_scores[i] = kp_score_at(seqX, (uint32_t)gp, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
                }
            }
            __syncthreads();

            #pragma unroll
            for (uint32_t step = 0U; step < 2U; step++) {
                uint32_t local_i = tid + step * 256U;
                if (local_i < wave_len) {
                    int L, R;
                    if (gecko_min_window_eval(s_scores, (int)(local_i + H), w, false, L, R)) {
                        uint32_t p_q = w_start + local_i;
                        uint32_t r_lo = p_q - (uint32_t)L;
                        uint32_t r_hi = p_q + (uint32_t)R - H;
                        if (r_lo < tile_min_x) r_lo = tile_min_x;
                        if (r_hi > tile_max_x) r_hi = tile_max_x;
                        if (r_lo <= r_hi) {
                            uint64_t h48;
                            uint32_t st, pal;
                            kp_extract_canon(seqX, p_q, kp_mask, rev_rshift_kp, h48, st, pal);
                            uint32_t part = (uint32_t)(h48 >> (shift_bits + part_shift));
                            atomicAdd(&s_hist[part], 1U);
                        }
                    }
                }
            }
            __syncthreads();
        }
    }

    for (uint32_t p = tid; p < P; p += 256U) {
        uint32_t c = s_hist[p];
        if (c > 0U) atomicAdd(&part_counts[p], c);
    }
}

/* Query Mode 2 — Stage 2: Single-CTA prefix sum over P <= 1024 partitions. */
static __global__ void min_l2_prefix_sum_parts_kernel(
    const uint32_t * __restrict__ part_counts,
    uint32_t P,
    uint32_t * __restrict__ part_offsets,
    uint32_t * __restrict__ part_cursor,
    uint32_t * __restrict__ part_block_off,
    uint32_t * __restrict__ totals
) {
    using BlockScan = cub::BlockScan<uint32_t, 256>;
    __shared__ typename BlockScan::TempStorage temp_q;
    __shared__ typename BlockScan::TempStorage temp_b;

    uint32_t tid = threadIdx.x;
    uint32_t q_arr[4];
    uint32_t b_arr[4];
    #pragma unroll
    for (int r = 0; r < 4; r++) {
        uint32_t p = tid * 4U + (uint32_t)r;
        uint32_t c = (p < P) ? part_counts[p] : 0U;
        q_arr[r] = c;
        b_arr[r] = (c + 255U) >> 8;
    }

    uint32_t total_q = 0U;
    uint32_t total_b = 0U;
    BlockScan(temp_q).ExclusiveSum(q_arr, q_arr, total_q);
    __syncthreads();
    BlockScan(temp_b).ExclusiveSum(b_arr, b_arr, total_b);

    #pragma unroll
    for (int r = 0; r < 4; r++) {
        uint32_t p = tid * 4U + (uint32_t)r;
        if (p < P) {
            part_offsets[p] = q_arr[r];
            part_cursor[p] = q_arr[r];
            part_block_off[p] = b_arr[r];
        }
    }
    if (tid == 0U) {
        part_offsets[P] = total_q;
        part_block_off[P] = total_b;
        totals[0] = total_q;
        totals[1] = total_b;
    }
}

/* Query Mode 2 — Stage 3: 2-phase CTA-privatized L2 scatter (ATOMS > 0, 8.25 KiB smem).
 * Packs (h48 | dl<<48 | span<<52 | q_strand<<56 | q_is_pal<<57) into 8 B + p_q into 4 B = 12 B. */
static __global__ __launch_bounds__(256, 4) void min_l2_scatter_parts_kernel(
    GpuSeqView seqX,
    int shift_bits,
    uint32_t P,
    int part_shift,
    int kp,
    int w,
    uint32_t tile_min_x,
    uint32_t tile_max_x,
    uint32_t p_begin,
    uint32_t p_end,
    uint32_t * __restrict__ part_cursor,
    uint64_t * __restrict__ out_q_hash_meta,
    uint32_t * __restrict__ out_q_pos
) {
    __shared__ uint64_t s_scores[MIN_SCORES_LEN];
    __shared__ uint32_t s_cnt[MIN_L2_MAX_PARTS];

    uint32_t tid = threadIdx.x;
    for (uint32_t p = tid; p < P; p += 256U) s_cnt[p] = 0U;
    __syncthreads();

    uint32_t chunk_start = p_begin + blockIdx.x * GECKO_EXTRACT_CHUNK_BP;
    if (chunk_start >= p_end) return;

    uint32_t H = (uint32_t)(w - 1);
    uint32_t valid_mask_kp = (1U << kp) - 1U;
    uint64_t kp_mask = (1ULL << (2 * kp)) - 1ULL;
    int rev_rshift_kp = 64 - 2 * kp;
    uint64_t total_kp = (seqX.len >= (uint64_t)kp) ? (seqX.len - (uint64_t)kp + 1ULL) : 0ULL;

    /* Phase 1: CTA-local partition count (ATOMS) */
    for (uint32_t wave = 0U; wave < MIN_WAVES_PER_CTA; wave++) {
        uint32_t w_start = chunk_start + wave * MIN_WAVE_BP;
        uint32_t wave_len = 0U;
        if (w_start < p_end) {
            uint32_t rem = p_end - w_start;
            wave_len = (rem < MIN_WAVE_BP) ? rem : MIN_WAVE_BP;
        }
        uint32_t fill_len = (wave_len > 0U) ? (wave_len + 2U * H) : 0U;

        for (uint32_t i = tid; i < fill_len; i += 256U) {
            int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
            if (gp < 0 || (uint64_t)gp >= total_kp) {
                s_scores[i] = ~0ULL;
            } else {
                s_scores[i] = kp_score_at(seqX, (uint32_t)gp, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
            }
        }
        __syncthreads();

        #pragma unroll
        for (uint32_t step = 0U; step < 2U; step++) {
            uint32_t local_i = tid + step * 256U;
            if (local_i < wave_len) {
                int L, R;
                if (gecko_min_window_eval(s_scores, (int)(local_i + H), w, false, L, R)) {
                    uint32_t p_q = w_start + local_i;
                    uint32_t r_lo = p_q - (uint32_t)L;
                    uint32_t r_hi = p_q + (uint32_t)R - H;
                    if (r_lo < tile_min_x) r_lo = tile_min_x;
                    if (r_hi > tile_max_x) r_hi = tile_max_x;
                    if (r_lo <= r_hi) {
                        uint64_t h48;
                        uint32_t st, pal;
                        kp_extract_canon(seqX, p_q, kp_mask, rev_rshift_kp, h48, st, pal);
                        uint32_t part = (uint32_t)(h48 >> (shift_bits + part_shift));
                        atomicAdd(&s_cnt[part], 1U);
                    }
                }
            }
        }
        __syncthreads();
    }

    /* Phase 1.5: Reserve 1 contiguous global slice per non-empty partition per CTA */
    for (uint32_t p = tid; p < P; p += 256U) {
        uint32_t c = s_cnt[p];
        if (c > 0U) s_cnt[p] = atomicAdd(&part_cursor[p], c);
    }
    __syncthreads();

    /* Phase 2: CTA-local scatter using ATOMS */
    for (uint32_t wave = 0U; wave < MIN_WAVES_PER_CTA; wave++) {
        uint32_t w_start = chunk_start + wave * MIN_WAVE_BP;
        uint32_t wave_len = 0U;
        if (w_start < p_end) {
            uint32_t rem = p_end - w_start;
            wave_len = (rem < MIN_WAVE_BP) ? rem : MIN_WAVE_BP;
        }
        uint32_t fill_len = (wave_len > 0U) ? (wave_len + 2U * H) : 0U;

        for (uint32_t i = tid; i < fill_len; i += 256U) {
            int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
            if (gp < 0 || (uint64_t)gp >= total_kp) {
                s_scores[i] = ~0ULL;
            } else {
                s_scores[i] = kp_score_at(seqX, (uint32_t)gp, kp, valid_mask_kp, kp_mask, rev_rshift_kp);
            }
        }
        __syncthreads();

        #pragma unroll
        for (uint32_t step = 0U; step < 2U; step++) {
            uint32_t local_i = tid + step * 256U;
            if (local_i < wave_len) {
                int L, R;
                if (gecko_min_window_eval(s_scores, (int)(local_i + H), w, false, L, R)) {
                    uint32_t p_q = w_start + local_i;
                    uint32_t r_lo = p_q - (uint32_t)L;
                    uint32_t r_hi = p_q + (uint32_t)R - H;
                    if (r_lo < tile_min_x) r_lo = tile_min_x;
                    if (r_hi > tile_max_x) r_hi = tile_max_x;
                    if (r_lo <= r_hi) {
                        uint64_t h48;
                        uint32_t st, pal;
                        kp_extract_canon(seqX, p_q, kp_mask, rev_rshift_kp, h48, st, pal);
                        uint32_t part = (uint32_t)(h48 >> (shift_bits + part_shift));
                        uint32_t slot = atomicAdd(&s_cnt[part], 1U);
                        uint64_t dl = (uint64_t)(p_q - r_lo);
                        uint64_t span = (uint64_t)(r_hi - r_lo);
                        out_q_hash_meta[slot] = h48
                            | (dl << 48)
                            | (span << 52)
                            | ((uint64_t)st << 56)
                            | ((uint64_t)pal << 57);
                        out_q_pos[slot] = p_q;
                    }
                }
            }
        }
        __syncthreads();
    }
}

/* Query Mode 2 — Stage 4: Wave-scheduled L2-resident minimizer probe (4.02 KiB smem, 100% T4 occ). */
template <bool HAS_OVF>
static __global__ __launch_bounds__(256, 4) void min_l2_probe_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ idx_bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int shift_bits,
    uint32_t tag_mask,
    uint32_t P,
    int wl,
    int kp,
    int w,
    bool compare_reverse,
    HitKeyLayout key_layout,
    const uint32_t * __restrict__ part_offsets,
    const uint32_t * __restrict__ part_block_off,
    const uint64_t * __restrict__ q_hash_meta,
    const uint32_t * __restrict__ q_pos,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits,
    uint32_t * __restrict__ ovf_q_r_lo,
    uint8_t  * __restrict__ ovf_q_span,
    uint32_t * __restrict__ ovf_q_count,
    uint32_t ovf_q_cap
) {
    __shared__ uint64_t s_hits[GECKO_HIT_STAGE_CAP];
    __shared__ uint32_t s_n;
    __shared__ uint32_t s_base;
    __shared__ uint32_t s_part;

    if (threadIdx.x == 0U) {
        s_n = 0U;
        uint32_t bid = blockIdx.x;
        uint32_t lo = 0U, hi = P;
        while (lo + 1U < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (__ldg(&part_block_off[mid]) <= bid) lo = mid;
            else hi = mid;
        }
        s_part = lo;
    }
    __syncthreads();

    uint32_t p = s_part;
    uint32_t local_blk = blockIdx.x - __ldg(&part_block_off[p]);
    uint32_t q_idx = __ldg(&part_offsets[p]) + local_blk * 256U + threadIdx.x;
    uint32_t q_end = __ldg(&part_offsets[p + 1U]);

    if (q_idx < q_end) {
        uint64_t hm = q_hash_meta[q_idx];
        uint32_t p_q = q_pos[q_idx];
        uint64_t h48 = hm & GECKO_MASK48;
        uint32_t dl = (uint32_t)((hm >> 48) & 15ULL);
        uint32_t span = (uint32_t)((hm >> 52) & 15ULL);
        uint32_t q_strand = (uint32_t)((hm >> 56) & 1ULL);
        uint32_t q_is_pal = (uint32_t)((hm >> 57) & 1ULL);
        uint32_t r_lo = p_q - dl;
        uint32_t r_hi = r_lo + span;

        min_probe_main_entry<HAS_OVF>(
            seqX, seqY, idx_bucket_off, idx_tags, idx_vals, dropped_bitmap,
            shift_bits, tag_mask, wl, kp, w, compare_reverse, key_layout,
            h48, p_q, r_lo, r_hi, q_strand, q_is_pal,
            s_hits, &s_n, out_hits, hit_count, max_hits,
            ovf_q_r_lo, ovf_q_span, ovf_q_count, ovf_q_cap);
    }

    __syncthreads();
    uint32_t flush_cnt = s_n;
    if (flush_cnt > GECKO_HIT_STAGE_CAP) flush_cnt = GECKO_HIT_STAGE_CAP;
    if (flush_cnt > 0U) {
        if (threadIdx.x == 0U) s_base = atomicAdd(hit_count, flush_cnt);
        __syncthreads();
        uint32_t base = s_base;
        if (base < max_hits) {
            uint32_t rem = max_hits - base;
            uint32_t n_write = (flush_cnt < rem) ? flush_cnt : rem;
            for (uint32_t k = threadIdx.x; k < n_write; k += 256U) {
                out_hits[base + k] = s_hits[k];
            }
        }
    }
}

/* ============================================================================
 * PART 3: Host Orchestrators (gpu_build_min_index, gpu_free_min_index, gpu_query_tile_minimizer)
 * ============================================================================ */

void gpu_free_min_index(GpuMinIndex &idx) {
    if (idx.d_tags)           (void)cudaFree(idx.d_tags);
    if (idx.d_vals)           (void)cudaFree(idx.d_vals);
    if (idx.d_bucket_off)     (void)cudaFree(idx.d_bucket_off);
    if (idx.d_dropped_bitmap) (void)cudaFree(idx.d_dropped_bitmap);
    if (idx.d_ovf_tags)       (void)cudaFree(idx.d_ovf_tags);
    if (idx.d_ovf_vals)       (void)cudaFree(idx.d_ovf_vals);
    if (idx.d_ovf_bucket_off) (void)cudaFree(idx.d_ovf_bucket_off);
    idx = {};
}

bool gpu_build_min_index(
    const GpuSeqView &seqY,
    int wl,
    uint32_t freq_limit,
    uint64_t owned_scratch_bytes,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    GpuMinIndex &out_idx,
    GeckoCudaProfile *prof
) {
    out_idx = {};
    int kp = 0, w = 0;
    if (!gecko_minimizer_params(wl, kp, w) || freq_limit == 0U || seqY.len < (uint64_t)wl) {
        return false;
    }

    out_idx.kp = kp;
    out_idx.w = w;
    out_idx.m_cap = gecko_minimizer_m_cap(freq_limit);

    uint64_t total_kp = seqY.len - (uint64_t)kp + 1ULL;
    uint32_t num_chunks = (uint32_t)((total_kp + (uint64_t)GECKO_EXTRACT_CHUNK_BP - 1ULL) /
                                     (uint64_t)GECKO_EXTRACT_CHUNK_BP);
    if (num_chunks == 0U) return true;

    /* Step 1a: Count canonical ref minimizers per 4096-bp chunk */
    uint32_t *d_chunk_offsets = nullptr;
    if (cudaMalloc(&d_chunk_offsets, (size_t)(num_chunks + 1U) * sizeof(uint32_t)) != cudaSuccess) {
        (void)cudaGetLastError();
        return false;
    }
    if (cudaMemset(d_chunk_offsets + num_chunks, 0, sizeof(uint32_t)) != cudaSuccess) {
        (void)cudaFree(d_chunk_offsets);
        (void)cudaGetLastError();
        return false;
    }

    size_t scan_bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, scan_bytes, d_chunk_offsets, d_chunk_offsets, (int)(num_chunks + 1U));
    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, scan_bytes);
    {
        CudaEventScope ev(prof, prof ? &prof->count_valid_ms : nullptr);
        count_ref_minimizers_kernel<<<num_chunks, 256>>>(seqY, kp, w, d_chunk_offsets);
        CUDA_CHECK_VOID(cudaGetLastError());
        cub::DeviceScan::ExclusiveSum(d_cub_temp, scan_bytes, d_chunk_offsets, d_chunk_offsets, (int)(num_chunks + 1U));
    }

    uint32_t n_raw_min = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&n_raw_min, d_chunk_offsets + num_chunks, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    out_idx.n_raw_min = n_raw_min;
    if (n_raw_min == 0U) {
        (void)cudaFree(d_chunk_offsets);
        return true;
    }

    /* VRAM & GECKO_CUDA_MAX_SLAB_KMERS guard before allocating sort buffers */
    if (const char *env = std::getenv("GECKO_CUDA_MAX_SLAB_KMERS")) {
        unsigned long long parsed = std::strtoull(env, nullptr, 10);
        if (parsed >= 100000ULL && (uint64_t)n_raw_min > (uint64_t)parsed) {
            (void)cudaFree(d_chunk_offsets);
            return false;
        }
    }
    {
        size_t free_vram = 0, total_vram = 0;
        if (cudaMemGetInfo(&free_vram, &total_vram) == cudaSuccess) {
            size_t base_reserve = std::max<size_t>(768ULL << 20, free_vram / 5ULL);
            size_t reserve = (base_reserve > owned_scratch_bytes + (256ULL << 20))
                ? (base_reserve - (size_t)owned_scratch_bytes)
                : (256ULL << 20);
            if (free_vram <= reserve || (uint64_t)n_raw_min * 36ULL > (uint64_t)(free_vram - reserve)) {
                (void)cudaFree(d_chunk_offsets);
                return false;
            }
        }
    }

    /* Step 1b: Extract ref minimizers into d_keys_a / d_vals_a */
    uint64_t *d_keys_a = nullptr, *d_keys_b = nullptr;
    uint32_t *d_vals_a = nullptr, *d_vals_b = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_keys_a, (size_t)n_raw_min * sizeof(uint64_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_vals_a, (size_t)n_raw_min * sizeof(uint32_t)));
    {
        CudaEventScope ev(prof, prof ? &prof->min_extract_ms : nullptr);
        extract_ref_minimizers_kernel<<<num_chunks, 256>>>(
            seqY, kp, w, d_chunk_offsets, d_keys_a, d_vals_a);
        CUDA_CHECK_VOID(cudaGetLastError());
    }
    (void)cudaFree(d_chunk_offsets);
    d_chunk_offsets = nullptr;

    /* Step 2: 48-bit CUB::DeviceRadixSort::SortPairs (6-arg CCCL 3.0 DoubleBuffer) */
    CUDA_CHECK_VOID(cudaMalloc(&d_keys_b, (size_t)n_raw_min * sizeof(uint64_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_vals_b, (size_t)n_raw_min * sizeof(uint32_t)));
    cub::DoubleBuffer<uint64_t> db_keys(d_keys_a, d_keys_b);
    cub::DoubleBuffer<uint32_t> db_vals(d_vals_a, d_vals_b);

    if (n_raw_min > 1U) {
        size_t sort_bytes = 0;
        cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes, db_keys, db_vals, (int)n_raw_min, 0, 48);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sort_bytes);
        CudaEventScope ev(prof, prof ? &prof->sort_pairs_ms : nullptr);
        cub::DeviceRadixSort::SortPairs(d_cub_temp, sort_bytes, db_keys, db_vals, (int)n_raw_min, 0, 48);
    }

    uint64_t *d_sorted_keys = db_keys.Current();
    uint32_t *d_sorted_vals = db_vals.Current();
    (void)cudaFree(db_keys.Alternate());
    (void)cudaFree(db_vals.Alternate());

    /* Step 3: Exact per-strand freq_limit resolution + M_CAP exact-k overflow routing */
    uint8_t *d_min_keep = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_min_keep, (size_t)(n_raw_min + 1U) * sizeof(uint8_t)));
    CUDA_CHECK_VOID(cudaMemset(d_min_keep, 1, (size_t)n_raw_min * sizeof(uint8_t)));
    CUDA_CHECK_VOID(cudaMemset(d_min_keep + n_raw_min, 0, sizeof(uint8_t)));

    if (freq_limit < n_raw_min) {
        /* Step 3a: Parallel O(log F) HF entry detection */
        uint8_t *d_is_hf = nullptr;
        uint32_t *d_hf_idx = nullptr;
        uint32_t *d_num_hf = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_is_hf, (size_t)n_raw_min * sizeof(uint8_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_hf_idx, (size_t)n_raw_min * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_num_hf, sizeof(uint32_t)));

        thrust::counting_iterator<uint32_t> count_it(0U);
        size_t sel_bytes = 0;
        cub::DeviceSelect::Flagged(nullptr, sel_bytes, count_it, d_is_hf, d_hf_idx, d_num_hf, (int)n_raw_min);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sel_bytes);
        {
            CudaEventScope ev(prof, prof ? &prof->min_hf_ms : nullptr);
            mark_hf_entries_kernel<<<(n_raw_min + 255U) / 256U, 256>>>(
                d_sorted_keys, n_raw_min, freq_limit, d_is_hf);
            CUDA_CHECK_VOID(cudaGetLastError());
            cub::DeviceSelect::Flagged(d_cub_temp, sel_bytes, count_it, d_is_hf, d_hf_idx, d_num_hf, (int)n_raw_min);
        }
        uint32_t n_hf = 0U;
        CUDA_CHECK_VOID(cudaMemcpy(&n_hf, d_num_hf, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        (void)cudaFree(d_is_hf);
        (void)cudaFree(d_num_hf);

        if (n_hf == 0U) {
            (void)cudaFree(d_hf_idx);
        } else {
            /* Step 3b: Emit candidate exact k-mers owned by HF minimizer entries (exact-once) */
            uint32_t *d_hf_r_lo = nullptr;
            uint8_t  *d_hf_r_cnt = nullptr;
            uint32_t *d_hf_off = nullptr;
            CUDA_CHECK_VOID(cudaMalloc(&d_hf_r_lo, (size_t)n_hf * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMalloc(&d_hf_r_cnt, (size_t)(n_hf + 1U) * sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMemset(d_hf_r_cnt + n_hf, 0, sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMalloc(&d_hf_off, (size_t)(n_hf + 1U) * sizeof(uint32_t)));

            auto cnt_u32_it = thrust::make_transform_iterator(d_hf_r_cnt, CastU8ToU32Op{});
            size_t hf_scan_bytes = 0;
            cub::DeviceScan::ExclusiveSum(nullptr, hf_scan_bytes, cnt_u32_it, d_hf_off, (int)(n_hf + 1U));
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, hf_scan_bytes);
            {
                CudaEventScope ev(prof, prof ? &prof->min_hf_ms : nullptr);
                count_hf_cand_kmers_kernel<<<(n_hf + 255U) / 256U, 256>>>(
                    seqY, d_sorted_vals, d_hf_idx, n_hf, kp, w, d_hf_r_lo, d_hf_r_cnt);
                CUDA_CHECK_VOID(cudaGetLastError());
                cub::DeviceScan::ExclusiveSum(d_cub_temp, hf_scan_bytes, cnt_u32_it, d_hf_off, (int)(n_hf + 1U));
            }

            uint32_t n_cand = 0U;
            CUDA_CHECK_VOID(cudaMemcpy(&n_cand, d_hf_off + n_hf, sizeof(uint32_t), cudaMemcpyDeviceToHost));

            /* Pathological homopolymer VRAM guard before allocating n_cand sort buffers */
            if (n_cand > 0U) {
                size_t free_vram = 0, total_vram = 0;
                if (cudaMemGetInfo(&free_vram, &total_vram) == cudaSuccess) {
                    if (free_vram <= (512ULL << 20) ||
                        (uint64_t)n_cand * 30ULL > (uint64_t)(free_vram - (512ULL << 20))) {
                        (void)cudaFree(d_hf_off);
                        (void)cudaFree(d_hf_r_cnt);
                        (void)cudaFree(d_hf_r_lo);
                        (void)cudaFree(d_hf_idx);
                        (void)cudaFree(d_min_keep);
                        (void)cudaFree(d_sorted_keys);
                        (void)cudaFree(d_sorted_vals);
                        gpu_free_min_index(out_idx);
                        return false;
                    }
                }

                /* Step 3c: Exact per-strand freq_limit -> d_dropped_bitmap */
                uint64_t *d_cand_keys_a = nullptr, *d_cand_keys_b = nullptr;
                uint32_t *d_cand_vals_a = nullptr, *d_cand_vals_b = nullptr;
                CUDA_CHECK_VOID(cudaMalloc(&d_cand_keys_a, (size_t)n_cand * sizeof(uint64_t)));
                CUDA_CHECK_VOID(cudaMalloc(&d_cand_vals_a, (size_t)n_cand * sizeof(uint32_t)));
                CUDA_CHECK_VOID(cudaMalloc(&d_cand_keys_b, (size_t)n_cand * sizeof(uint64_t)));
                CUDA_CHECK_VOID(cudaMalloc(&d_cand_vals_b, (size_t)n_cand * sizeof(uint32_t)));

                size_t bm_words = (size_t)((seqY.len + 31ULL) / 32ULL) + 2;
                uint32_t *d_any_dropped = nullptr;
                CUDA_CHECK_VOID(cudaMalloc(&out_idx.d_dropped_bitmap, bm_words * sizeof(uint32_t)));
                CUDA_CHECK_VOID(cudaMemset(out_idx.d_dropped_bitmap, 0, bm_words * sizeof(uint32_t)));
                CUDA_CHECK_VOID(cudaMalloc(&d_any_dropped, sizeof(uint32_t)));
                CUDA_CHECK_VOID(cudaMemset(d_any_dropped, 0, sizeof(uint32_t)));

                cub::DoubleBuffer<uint64_t> db_ck(d_cand_keys_a, d_cand_keys_b);
                cub::DoubleBuffer<uint32_t> db_cv(d_cand_vals_a, d_cand_vals_b);
                size_t csort_bytes = 0;
                if (n_cand > 1U) {
                    cub::DeviceRadixSort::SortPairs(nullptr, csort_bytes, db_ck, db_cv, (int)n_cand, 0, 64);
                    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, csort_bytes);
                }

                {
                    CudaEventScope ev(prof, prof ? &prof->min_hf_ms : nullptr);
                    extract_hf_cand_kmers_kernel<<<(n_hf + 255U) / 256U, 256>>>(
                        seqY, d_hf_r_lo, d_hf_r_cnt, d_hf_off, n_hf, wl, d_cand_keys_a, d_cand_vals_a);
                    CUDA_CHECK_VOID(cudaGetLastError());
                    if (n_cand > 1U) {
                        cub::DeviceRadixSort::SortPairs(d_cub_temp, csort_bytes, db_ck, db_cv, (int)n_cand, 0, 64);
                    }
                    mark_hf_kmer_dropped_kernel<<<(n_cand + 255U) / 256U, 256>>>(
                        db_ck.Current(), db_cv.Current(), n_cand, freq_limit,
                        out_idx.d_dropped_bitmap, d_any_dropped);
                    CUDA_CHECK_VOID(cudaGetLastError());
                }

                (void)cudaFree(d_cand_keys_a);
                (void)cudaFree(d_cand_keys_b);
                (void)cudaFree(d_cand_vals_a);
                (void)cudaFree(d_cand_vals_b);

                uint32_t h_any_dropped = 0U;
                CUDA_CHECK_VOID(cudaMemcpy(&h_any_dropped, d_any_dropped, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                (void)cudaFree(d_any_dropped);
                if (h_any_dropped == 0U) {
                    (void)cudaFree(out_idx.d_dropped_bitmap);
                    out_idx.d_dropped_bitmap = nullptr;
                }
            }
            (void)cudaFree(d_hf_off);
            d_hf_off = nullptr;

            /* Step 3d: Survival check (L_ge/R_ge) + parallel M_CAP overflow classification */
            uint32_t *d_keep_prefix = nullptr;
            uint8_t  *d_hf_is_ovf = nullptr;
            uint32_t *d_ovf_key_cnt = nullptr;
            CUDA_CHECK_VOID(cudaMalloc(&d_keep_prefix, (size_t)(n_raw_min + 1U) * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMalloc(&d_hf_is_ovf, (size_t)n_hf * sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMemset(d_hf_is_ovf, 0, (size_t)n_hf * sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMalloc(&d_ovf_key_cnt, sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMemset(d_ovf_key_cnt, 0, sizeof(uint32_t)));

            auto keep_u32_it = thrust::make_transform_iterator(d_min_keep, CastU8ToU32Op{});
            size_t kp_scan_bytes = 0;
            cub::DeviceScan::ExclusiveSum(nullptr, kp_scan_bytes, keep_u32_it, d_keep_prefix, (int)(n_raw_min + 1U));
            gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, kp_scan_bytes);

            {
                CudaEventScope ev(prof, prof ? &prof->min_hf_ms : nullptr);
                if (out_idx.d_dropped_bitmap != nullptr) {
                    eval_hf_min_survival_kernel<<<(n_hf + 255U) / 256U, 256>>>(
                        seqY, d_sorted_vals, d_hf_idx, n_hf, kp, w, out_idx.d_dropped_bitmap, d_min_keep);
                    CUDA_CHECK_VOID(cudaGetLastError());
                }
                cub::DeviceScan::ExclusiveSum(d_cub_temp, kp_scan_bytes, keep_u32_it, d_keep_prefix, (int)(n_raw_min + 1U));
                classify_overflow_entries_kernel<<<(n_hf + 255U) / 256U, 256>>>(
                    d_sorted_keys, d_keep_prefix, d_hf_idx, n_raw_min, n_hf, out_idx.m_cap,
                    d_min_keep, d_hf_is_ovf, d_ovf_key_cnt);
                CUDA_CHECK_VOID(cudaGetLastError());
            }

            (void)cudaFree(d_keep_prefix);
            uint32_t h_ovf_routed_keys = 0U;
            CUDA_CHECK_VOID(cudaMemcpy(&h_ovf_routed_keys, d_ovf_key_cnt, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            (void)cudaFree(d_ovf_key_cnt);
            out_idx.ovf_routed_keys = h_ovf_routed_keys;

            /* Step 3e: Build tiny exact-k Overflow Hash-CSR iff h_ovf_routed_keys > 0 */
            if (h_ovf_routed_keys > 0U) {
                uint8_t  *d_ovf_cnt = nullptr;
                uint32_t *d_ovf_off = nullptr;
                CUDA_CHECK_VOID(cudaMalloc(&d_ovf_cnt, (size_t)(n_hf + 1U) * sizeof(uint8_t)));
                CUDA_CHECK_VOID(cudaMemset(d_ovf_cnt + n_hf, 0, sizeof(uint8_t)));
                CUDA_CHECK_VOID(cudaMalloc(&d_ovf_off, (size_t)(n_hf + 1U) * sizeof(uint32_t)));

                auto ovf_u32_it = thrust::make_transform_iterator(d_ovf_cnt, CastU8ToU32Op{});
                size_t oscan_bytes = 0;
                cub::DeviceScan::ExclusiveSum(nullptr, oscan_bytes, ovf_u32_it, d_ovf_off, (int)(n_hf + 1U));
                gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, oscan_bytes);
                {
                    CudaEventScope ev(prof, prof ? &prof->min_ovf_ms : nullptr);
                    count_ovf_kmers_kernel<<<(n_hf + 255U) / 256U, 256>>>(
                        d_hf_is_ovf, d_hf_r_lo, d_hf_r_cnt, n_hf, out_idx.d_dropped_bitmap, d_ovf_cnt);
                    CUDA_CHECK_VOID(cudaGetLastError());
                    cub::DeviceScan::ExclusiveSum(d_cub_temp, oscan_bytes, ovf_u32_it, d_ovf_off, (int)(n_hf + 1U));
                }
                (void)cudaFree(d_ovf_cnt);

                uint32_t n_ovf_kmers = 0U;
                CUDA_CHECK_VOID(cudaMemcpy(&n_ovf_kmers, d_ovf_off + n_hf, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                out_idx.ovf_n_kept = n_ovf_kmers;

                if (n_ovf_kmers > 0U) {
                    uint64_t *d_ok_a = nullptr, *d_ok_b = nullptr;
                    uint32_t *d_ov_a = nullptr, *d_ov_b = nullptr;
                    CUDA_CHECK_VOID(cudaMalloc(&d_ok_a, (size_t)n_ovf_kmers * sizeof(uint64_t)));
                    CUDA_CHECK_VOID(cudaMalloc(&d_ov_a, (size_t)n_ovf_kmers * sizeof(uint32_t)));
                    CUDA_CHECK_VOID(cudaMalloc(&d_ok_b, (size_t)n_ovf_kmers * sizeof(uint64_t)));
                    CUDA_CHECK_VOID(cudaMalloc(&d_ov_b, (size_t)n_ovf_kmers * sizeof(uint32_t)));

                    cub::DoubleBuffer<uint64_t> db_ok(d_ok_a, d_ok_b);
                    cub::DoubleBuffer<uint32_t> db_ov(d_ov_a, d_ov_b);
                    size_t osort_bytes = 0;
                    if (n_ovf_kmers > 1U) {
                        cub::DeviceRadixSort::SortPairs(nullptr, osort_bytes, db_ok, db_ov, (int)n_ovf_kmers, 0, 64);
                        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, osort_bytes);
                    }

                    int ob_bits = 10;
                    while (ob_bits < 24 && (1ULL << ob_bits) * 4ULL < (uint64_t)n_ovf_kmers) ob_bits++;
                    uint32_t ob_cnt = 1U << ob_bits;
                    int o_shift = 64 - ob_bits;
                    int o_tag_shift = o_shift - 32;

                    out_idx.ovf_bucket_bits = ob_bits;
                    out_idx.ovf_num_buckets = ob_cnt;
                    out_idx.ovf_shift_bits = o_shift;
                    out_idx.ovf_tag_shift = o_tag_shift;

                    CUDA_CHECK_VOID(cudaMalloc(&out_idx.d_ovf_tags, (size_t)n_ovf_kmers * sizeof(uint32_t)));
                    CUDA_CHECK_VOID(cudaMalloc(&out_idx.d_ovf_bucket_off, (size_t)(ob_cnt + 1U) * sizeof(uint32_t)));
                    CUDA_CHECK_VOID(cudaMemset(out_idx.d_ovf_bucket_off, 0, (size_t)(ob_cnt + 1U) * sizeof(uint32_t)));

                    size_t ob_scan_bytes = 0;
                    cub::DeviceScan::ExclusiveSum(
                        nullptr, ob_scan_bytes, out_idx.d_ovf_bucket_off, out_idx.d_ovf_bucket_off, (int)(ob_cnt + 1U));
                    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, ob_scan_bytes);

                    {
                        CudaEventScope ev(prof, prof ? &prof->min_ovf_ms : nullptr);
                        extract_ovf_kmers_kernel<<<(n_hf + 255U) / 256U, 256>>>(
                            seqY, d_hf_is_ovf, d_hf_r_lo, d_hf_r_cnt, d_ovf_off, n_hf, wl,
                            out_idx.d_dropped_bitmap, d_ok_a, d_ov_a);
                        CUDA_CHECK_VOID(cudaGetLastError());
                        if (n_ovf_kmers > 1U) {
                            cub::DeviceRadixSort::SortPairs(d_cub_temp, osort_bytes, db_ok, db_ov, (int)n_ovf_kmers, 0, 64);
                        }
                        scatter_ovf_csr_kernel<<<(n_ovf_kmers + 255U) / 256U, 256>>>(
                            db_ok.Current(), n_ovf_kmers, o_shift, o_tag_shift,
                            out_idx.d_ovf_tags, out_idx.d_ovf_bucket_off);
                        CUDA_CHECK_VOID(cudaGetLastError());
                        cub::DeviceScan::ExclusiveSum(
                            d_cub_temp, ob_scan_bytes, out_idx.d_ovf_bucket_off, out_idx.d_ovf_bucket_off, (int)(ob_cnt + 1U));
                    }

                    out_idx.d_ovf_vals = db_ov.Current();
                    (void)cudaFree(db_ov.Alternate());
                    (void)cudaFree(d_ok_a);
                    (void)cudaFree(d_ok_b);
                }
                (void)cudaFree(d_ovf_off);
            }

            (void)cudaFree(d_hf_is_ovf);
            (void)cudaFree(d_hf_r_cnt);
            (void)cudaFree(d_hf_r_lo);
            (void)cudaFree(d_hf_idx);
        }
    }

    /* Step 4: Fused single-pass main minimizer Hash-CSR build (bucket_bits in [17 .. 26]) */
    uint32_t *d_out_pos = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_out_pos, (size_t)(n_raw_min + 1U) * sizeof(uint32_t)));
    auto nz_it = thrust::make_transform_iterator(d_min_keep, IsNonZeroU8ToU32Op{});
    size_t out_scan_bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, out_scan_bytes, nz_it, d_out_pos, (int)(n_raw_min + 1U));
    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, out_scan_bytes);
    {
        CudaEventScope ev(prof, prof ? &prof->bucket_build_ms : nullptr);
        cub::DeviceScan::ExclusiveSum(d_cub_temp, out_scan_bytes, nz_it, d_out_pos, (int)(n_raw_min + 1U));
    }

    uint32_t n_kept_main = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&n_kept_main, d_out_pos + n_raw_min, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    out_idx.n_kept_main = n_kept_main;

    /* Clamp bucket_bits in [17 .. 26] so shift_bits = 48 - bucket_bits in [22 .. 31],
     * guaranteeing bit 31 of d_tags is strictly 0 for normal tags (0.000% kp-mer collisions!). */
    int bucket_bits = 17;
    while (bucket_bits < 26 && (1ULL << bucket_bits) * 4ULL < (uint64_t)n_kept_main) {
        bucket_bits++;
    }
    uint32_t num_buckets = 1U << bucket_bits;
    int shift_bits = 48 - bucket_bits;
    uint32_t tag_mask = (1U << shift_bits) - 1U;

    out_idx.bucket_bits = bucket_bits;
    out_idx.num_buckets = num_buckets;
    out_idx.shift_bits = shift_bits;
    out_idx.tag_mask = tag_mask;

    if (n_kept_main > 0U) {
        CUDA_CHECK_VOID(cudaMalloc(&out_idx.d_tags, (size_t)n_kept_main * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&out_idx.d_vals, (size_t)n_kept_main * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMalloc(&out_idx.d_bucket_off, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));
        CUDA_CHECK_VOID(cudaMemset(out_idx.d_bucket_off, 0, (size_t)(num_buckets + 1U) * sizeof(uint32_t)));

        size_t bscan_bytes = 0;
        cub::DeviceScan::ExclusiveSum(
            nullptr, bscan_bytes, out_idx.d_bucket_off, out_idx.d_bucket_off, (int)(num_buckets + 1U));
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, bscan_bytes);

        {
            CudaEventScope ev(prof, prof ? &prof->bucket_build_ms : nullptr);
            scatter_min_csr_kernel<<<(n_raw_min + 255U) / 256U, 256>>>(
                d_sorted_keys, d_sorted_vals, d_min_keep, d_out_pos, n_raw_min,
                shift_bits, tag_mask, out_idx.d_tags, out_idx.d_vals, out_idx.d_bucket_off);
            CUDA_CHECK_VOID(cudaGetLastError());
            cub::DeviceScan::ExclusiveSum(
                d_cub_temp, bscan_bytes, out_idx.d_bucket_off, out_idx.d_bucket_off, (int)(num_buckets + 1U));
        }
    }

    (void)cudaFree(d_out_pos);
    (void)cudaFree(d_min_keep);
    (void)cudaFree(d_sorted_keys);
    (void)cudaFree(d_sorted_vals);
    return true;
}

struct MinL2Plan {
    uint32_t P;
    int p_bits;
    int part_shift;
};

static MinL2Plan choose_min_l2_part_plan(const GpuMinIndex &idx, int l2_cache_bytes) {
    MinL2Plan plan{1U, 0, idx.bucket_bits};
    if (idx.n_kept_main == 0U || idx.num_buckets == 0U) return plan;

    size_t probe_bytes = ((size_t)idx.num_buckets + 1ULL) * sizeof(uint32_t) +
                         (size_t)idx.n_kept_main * sizeof(uint32_t);
    size_t l2 = (l2_cache_bytes > 0) ? (size_t)l2_cache_bytes : (4ULL << 20);
    size_t target = (l2 * 3ULL) / 4ULL;
    if (probe_bytes <= target) return plan;

    int p_bits = 1;
    while (p_bits < 10 && p_bits < idx.bucket_bits && (probe_bytes >> p_bits) > target) {
        p_bits++;
    }
    plan.P = 1U << p_bits;
    plan.p_bits = p_bits;
    plan.part_shift = idx.bucket_bits - p_bits;
    return plan;
}

uint32_t gpu_query_tile_minimizer(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuMinIndex &idx,
    int wl,
    bool compare_reverse,
    uint64_t tile_start,
    uint64_t tile_end,
    HitKeyLayout key_layout,
    int l2_cache_bytes,
    int query_mode,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint32_t *d_hit_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed,
    uint64_t *&out_sorted_hits,
    uint64_t *&out_alt_buf,
    GeckoCudaProfile *prof
) {
    out_overflowed = false;
    out_raw_leftmax = 0U;
    out_sorted_hits = d_hits;
    out_alt_buf = d_hits_sorted;
    if (idx.n_kept_main == 0U || tile_start >= tile_end || tile_end < (uint64_t)wl) return 0U;

    uint32_t tile_min_x = (tile_start >= (uint64_t)(wl - 1))
        ? (uint32_t)(tile_start - (uint64_t)(wl - 1))
        : 0U;
    uint32_t tile_max_x = (uint32_t)(tile_end - (uint64_t)wl);
    if (tile_min_x > tile_max_x) return 0U;

    uint32_t p_begin = tile_min_x;
    uint32_t p_end = tile_max_x + (uint32_t)idx.w;
    uint32_t num_p = p_end - p_begin;
    uint32_t grid_chunks = (num_p + GECKO_EXTRACT_CHUNK_BP - 1U) / GECKO_EXTRACT_CHUNK_BP;
    bool has_ovf = (idx.ovf_n_kept > 0U && idx.d_ovf_bucket_off != nullptr);

    CUDA_CHECK_VOID(cudaMemset(d_hit_count, 0, sizeof(uint32_t)));

    size_t scratch_bytes = (size_t)max_hits * sizeof(uint64_t);
    auto align256 = [](size_t n) -> size_t { return (n + 255ULL) & ~255ULL; };
    uint8_t *base_ptr = reinterpret_cast<uint8_t *>(d_hits_sorted);

    /* Optional deferred overflow queue carved at the tail of idle d_hits_sorted */
    uint32_t *d_ovf_q_r_lo = nullptr;
    uint8_t  *d_ovf_q_span = nullptr;
    uint32_t *d_ovf_q_count = nullptr;
    uint32_t ovf_q_cap = 0U;

    bool ran_l2_part = false;
    if (query_mode != GECKO_QMODE_MONOLITHIC) {
        MinL2Plan plan = choose_min_l2_part_plan(idx, l2_cache_bytes);
        if (plan.P >= 2U && (uint64_t)num_p >= (uint64_t)plan.P * 2048ULL) {
            size_t meta_bytes =
                align256((size_t)plan.P * sizeof(uint32_t)) +
                align256(((size_t)plan.P + 1ULL) * sizeof(uint32_t)) +
                align256((size_t)plan.P * sizeof(uint32_t)) +
                align256(((size_t)plan.P + 1ULL) * sizeof(uint32_t)) +
                align256(2ULL * sizeof(uint32_t)) +
                align256(sizeof(uint32_t));

            if (meta_bytes + 4096ULL < scratch_bytes) {
                /* Reserve top metadata block at start of scratch_bytes for count/prefix */
                size_t off_counts = 0ULL;
                size_t off_offs   = align256(off_counts + (size_t)plan.P * sizeof(uint32_t));
                size_t off_cursor = align256(off_offs   + ((size_t)plan.P + 1ULL) * sizeof(uint32_t));
                size_t off_blkoff = align256(off_cursor + (size_t)plan.P * sizeof(uint32_t));
                size_t off_totals = align256(off_blkoff + ((size_t)plan.P + 1ULL) * sizeof(uint32_t));
                size_t off_ovf_cnt = align256(off_totals + 2ULL * sizeof(uint32_t));
                size_t off_payload = align256(off_ovf_cnt + sizeof(uint32_t));

                uint32_t *d_part_counts    = reinterpret_cast<uint32_t *>(base_ptr + off_counts);
                uint32_t *d_part_offsets   = reinterpret_cast<uint32_t *>(base_ptr + off_offs);
                uint32_t *d_part_cursor    = reinterpret_cast<uint32_t *>(base_ptr + off_cursor);
                uint32_t *d_part_block_off = reinterpret_cast<uint32_t *>(base_ptr + off_blkoff);
                uint32_t *d_totals         = reinterpret_cast<uint32_t *>(base_ptr + off_totals);
                uint32_t *d_ovf_cnt_slot   = reinterpret_cast<uint32_t *>(base_ptr + off_ovf_cnt);

                CUDA_CHECK_VOID(cudaMemset(d_part_counts, 0, (size_t)plan.P * sizeof(uint32_t)));
                {
                    CudaEventScope ev(prof, prof ? &prof->l2_count_ms : nullptr);
                    min_l2_count_parts_kernel<<<grid_chunks, 256>>>(
                        seqX, idx.shift_bits, plan.P, plan.part_shift,
                        idx.kp, idx.w, tile_min_x, tile_max_x, p_begin, p_end, d_part_counts);
                    CUDA_CHECK_VOID(cudaGetLastError());
                    min_l2_prefix_sum_parts_kernel<<<1, 256>>>(
                        d_part_counts, plan.P, d_part_offsets, d_part_cursor, d_part_block_off, d_totals);
                    CUDA_CHECK_VOID(cudaGetLastError());
                }

                uint32_t h_totals[2] = {0U, 0U};
                CUDA_CHECK_VOID(cudaMemcpy(h_totals, d_totals, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                if (h_totals[0] == 0U) {
                    if (prof && prof->enabled) prof->l2_part_P = plan.P;
                    return 0U;
                }

                uint32_t total_q = h_totals[0];
                size_t off_q_hash = off_payload;
                size_t off_q_pos  = align256(off_q_hash + (size_t)total_q * sizeof(uint64_t));
                size_t off_ovf_lo = align256(off_q_pos  + (size_t)total_q * sizeof(uint32_t));
                size_t ovf_need_bytes = has_ovf
                    ? (align256((size_t)total_q * sizeof(uint32_t)) + align256((size_t)total_q))
                    : 0ULL;

                /* Ensure payload + optional full-capacity overflow queue fits inside d_hits_sorted.
                 * If not (e.g. pathological 64 Mb pure homopolymer), fall back to monolithic. */
                if (off_ovf_lo + ovf_need_bytes <= scratch_bytes) {
                    if (prof && prof->enabled) prof->l2_part_P = plan.P;

                    uint64_t *d_q_hash_meta = reinterpret_cast<uint64_t *>(base_ptr + off_q_hash);
                    uint32_t *d_q_pos       = reinterpret_cast<uint32_t *>(base_ptr + off_q_pos);

                    if (has_ovf) {
                        d_ovf_q_r_lo = reinterpret_cast<uint32_t *>(base_ptr + off_ovf_lo);
                        d_ovf_q_span = reinterpret_cast<uint8_t *>(
                            base_ptr + align256(off_ovf_lo + (size_t)total_q * sizeof(uint32_t)));
                        d_ovf_q_count = d_ovf_cnt_slot;
                        ovf_q_cap = total_q;
                        CUDA_CHECK_VOID(cudaMemset(d_ovf_q_count, 0, sizeof(uint32_t)));
                    }

                    {
                        CudaEventScope ev(prof, prof ? &prof->l2_scatter_ms : nullptr);
                        min_l2_scatter_parts_kernel<<<grid_chunks, 256>>>(
                            seqX, idx.shift_bits, plan.P, plan.part_shift,
                            idx.kp, idx.w, tile_min_x, tile_max_x, p_begin, p_end,
                            d_part_cursor, d_q_hash_meta, d_q_pos);
                        CUDA_CHECK_VOID(cudaGetLastError());
                    }

                    {
                        CudaEventScope ev(prof, prof ? &prof->l2_probe_ms : nullptr);
                        if (has_ovf) {
                            min_l2_probe_kernel<true><<<h_totals[1], 256>>>(
                                seqX, seqY, idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap,
                                idx.shift_bits, idx.tag_mask, plan.P, wl, idx.kp, idx.w, compare_reverse,
                                key_layout, d_part_offsets, d_part_block_off, d_q_hash_meta, d_q_pos,
                                d_hits, d_hit_count, max_hits,
                                d_ovf_q_r_lo, d_ovf_q_span, d_ovf_q_count, ovf_q_cap);
                        } else {
                            min_l2_probe_kernel<false><<<h_totals[1], 256>>>(
                                seqX, seqY, idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap,
                                idx.shift_bits, idx.tag_mask, plan.P, wl, idx.kp, idx.w, compare_reverse,
                                key_layout, d_part_offsets, d_part_block_off, d_q_hash_meta, d_q_pos,
                                d_hits, d_hit_count, max_hits,
                                nullptr, nullptr, nullptr, 0U);
                        }
                        CUDA_CHECK_VOID(cudaGetLastError());
                    }

                    if (has_ovf) {
                        uint32_t n_ovf_q = 0U;
                        CUDA_CHECK_VOID(cudaMemcpy(&n_ovf_q, d_ovf_q_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        if (n_ovf_q > ovf_q_cap) {
                            out_overflowed = true;
                            return 0U;
                        }
                        if (n_ovf_q > 0U) {
                            CudaEventScope ev(prof, prof ? &prof->l2_probe_ms : nullptr);
                            min_probe_overflow_kernel<<<(n_ovf_q + 255U) / 256U, 256>>>(
                                seqX, seqY, idx.d_ovf_bucket_off, idx.d_ovf_tags, idx.d_ovf_vals,
                                idx.d_dropped_bitmap, idx.ovf_shift_bits, idx.ovf_tag_shift,
                                wl, compare_reverse, key_layout,
                                d_ovf_q_r_lo, d_ovf_q_span, n_ovf_q,
                                d_hits, d_hit_count, max_hits);
                            CUDA_CHECK_VOID(cudaGetLastError());
                        }
                    }

                    ran_l2_part = true;
                }
            }
        }
    }

    if (!ran_l2_part) {
        if (has_ovf) {
            /* Carve deferred overflow queue directly inside idle d_hits_sorted (5 B/item) */
            size_t cap = (scratch_bytes - 512ULL) / 5ULL;
            if (cap > (size_t)num_p) cap = (size_t)num_p;
            if (cap == 0ULL) cap = 1ULL;
            ovf_q_cap = (uint32_t)cap;
            d_ovf_q_count = reinterpret_cast<uint32_t *>(base_ptr);
            d_ovf_q_r_lo  = reinterpret_cast<uint32_t *>(base_ptr + 256ULL);
            d_ovf_q_span  = reinterpret_cast<uint8_t *>(
                base_ptr + align256(256ULL + (size_t)ovf_q_cap * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMemset(d_ovf_q_count, 0, sizeof(uint32_t)));
        }

        {
            CudaEventScope ev(prof, prof ? &prof->mono_probe_ms : nullptr);
            if (has_ovf) {
                min_query_monolithic_kernel<true><<<grid_chunks, 256>>>(
                    seqX, seqY, idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap,
                    idx.shift_bits, idx.tag_mask, wl, idx.kp, idx.w, compare_reverse,
                    tile_min_x, tile_max_x, p_begin, p_end, key_layout,
                    d_hits, d_hit_count, max_hits,
                    d_ovf_q_r_lo, d_ovf_q_span, d_ovf_q_count, ovf_q_cap);
            } else {
                min_query_monolithic_kernel<false><<<grid_chunks, 256>>>(
                    seqX, seqY, idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap,
                    idx.shift_bits, idx.tag_mask, wl, idx.kp, idx.w, compare_reverse,
                    tile_min_x, tile_max_x, p_begin, p_end, key_layout,
                    d_hits, d_hit_count, max_hits,
                    nullptr, nullptr, nullptr, 0U);
            }
            CUDA_CHECK_VOID(cudaGetLastError());
        }

        if (has_ovf) {
            uint32_t n_ovf_q = 0U;
            CUDA_CHECK_VOID(cudaMemcpy(&n_ovf_q, d_ovf_q_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            if (n_ovf_q > ovf_q_cap) {
                out_overflowed = true;
                return 0U;
            }
            if (n_ovf_q > 0U) {
                CudaEventScope ev(prof, prof ? &prof->mono_probe_ms : nullptr);
                min_probe_overflow_kernel<<<(n_ovf_q + 255U) / 256U, 256>>>(
                    seqX, seqY, idx.d_ovf_bucket_off, idx.d_ovf_tags, idx.d_ovf_vals,
                    idx.d_dropped_bitmap, idx.ovf_shift_bits, idx.ovf_tag_shift,
                    wl, compare_reverse, key_layout,
                    d_ovf_q_r_lo, d_ovf_q_span, n_ovf_q,
                    d_hits, d_hit_count, max_hits);
                CUDA_CHECK_VOID(cudaGetLastError());
            }
        }
    }

    uint32_t h_nhits = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nhits, d_hit_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nhits > max_hits) {
        out_overflowed = true;
        return 0U;
    }
    out_raw_leftmax = h_nhits;
    if (h_nhits <= 1U) {
        out_sorted_hits = d_hits;
        out_alt_buf = d_hits_sorted;
        return h_nhits;
    }

    /* 6-arg CCCL 3.0 DoubleBuffer bit-bounded SortKeys (identical to exact path) */
    int end_bit = (int)key_layout.end_bit;
    cub::DoubleBuffer<uint64_t> d_keys_buf(d_hits, d_hits_sorted);
    size_t sort_bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, sort_bytes, d_keys_buf, (int)h_nhits, 0, end_bit);
    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sort_bytes);
    {
        CudaEventScope ev(prof, prof ? &prof->sort_hits_ms : nullptr);
        cub::DeviceRadixSort::SortKeys(d_cub_temp, sort_bytes, d_keys_buf, (int)h_nhits, 0, end_bit);
    }

    out_sorted_hits = d_keys_buf.Current();
    out_alt_buf = d_keys_buf.Alternate();
    return h_nhits;
}
