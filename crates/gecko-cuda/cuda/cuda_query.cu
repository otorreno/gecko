/*
 * cuda_query.cu
 *
 * Stage 2: Tiled 2-bit GPU Query + Lock-Free Left/Right-Maximal Seed Filtering:
 *  - F14: Direct-CSR query path (wl <= 12, zero fmix64/d_tags/2-bit verify),
 *    block-level shared-memory hit staging (STAGE_CAP=512, lockstep wave flush,
 *    lock-free overflow spill), lazy right_ext evaluation, and tile-relative
 *    bit-bounded CUB SortKeys via HitKeyLayout.
 *  - F15: Multi-slab skip_lr fix via global dropped_bitmap / freq_limit == 0.
 *  - F16: 1-bit ref_unique_bm active-diagonal bypass (STRIDE=64) AND auto-tuned
 *    L2-resident wave-scheduled hash-partitioned probe aliasing idle d_hits_sorted.
 */

#include <cub/cub.cuh>
#include "cuda_common.cuh"

static __device__ __forceinline__ void stage_hit_dev(
    uint64_t *s_hits,
    uint32_t *s_cnt,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits,
    uint64_t packed_key
) {
    uint32_t lslot = atomicAdd(s_cnt, 1U);
    if (lslot < GECKO_HIT_STAGE_CAP) {
        s_hits[lslot] = packed_key;
    } else {
        /* Lock-free global overflow spill: ZERO __syncthreads() -> 100% deadlock-free! */
        if (*(volatile const uint32_t *)hit_count <= max_hits) {
            uint32_t gslot = atomicAdd(hit_count, 1U);
            if (gslot < max_hits) {
                out_hits[gslot] = packed_key;
            }
        }
    }
}

/* ---- F14 + F16: Monolithic / Direct-CSR / Unique-Bypass Tiled Query Kernel ---- */
template <bool DIRECT_CSR, bool USE_UNIQUE_BM>
static __global__ void query_left_maximal_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    const uint32_t * __restrict__ ref_unique_bm,
    int shift_bits,
    int tag_shift,
    uint32_t part_start,
    uint32_t part_end,
    uint32_t bucket_base,
    int wl,
    uint64_t seed_mask,
    bool compare_reverse,
    bool has_global_dropped_bm,
    uint64_t tile_start,
    uint64_t tile_end,
    uint32_t num_waves,
    HitKeyLayout key_layout,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits,
    uint32_t * __restrict__ stats_buf
) {
    __shared__ uint64_t s_hits[GECKO_HIT_STAGE_CAP];
    __shared__ uint32_t s_cnt;
    __shared__ uint32_t s_base;
    __shared__ uint32_t s_valid_q;
    __shared__ uint32_t s_bypass_q;

    if (threadIdx.x == 0) {
        s_cnt = 0U;
    }
    __syncthreads();

    const uint64_t stride_bp = (uint64_t)num_waves * 16ULL;
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t t_start = tile_start + gtid * stride_bp;
    const bool active = (t_start < tile_end);
    const uint64_t t_end = active
        ? ((t_start + stride_bp < tile_end) ? (t_start + stride_bp) : tile_end)
        : t_start;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : mask;
    bool is_spaced = (eff_mask != mask);
    bool part_filter = !DIRECT_CSR && (part_start != 0U || part_end != 256U);
    bool skip_lr = is_spaced || (part_filter && !has_global_dropped_bm);
    int rc_shift = 2 * (wl - 1);
    int rev_rshift = 64 - 2 * wl;

    uint64_t i = active
        ? ((t_start >= (uint64_t)(wl - 1)) ? (t_start - (uint64_t)(wl - 1)) : 0ULL)
        : t_end;
    uint64_t kf = 0ULL, kr = 0ULL;
    int valid = 0;

    /* F16: Per-thread active-diagonal bypass state (zero smem, 3 registers) */
    bool act_valid = false;
    uint32_t act_strand = 0U;
    uint32_t act_posY_prev = 0U;
    uint32_t valid_q_local = 0U;
    uint32_t bypassed_local = 0U;

    for (uint32_t wave = 0U; wave < num_waves; ++wave) {
        if (active) {
            uint64_t sub_end = t_start + (uint64_t)(wave + 1U) * 16ULL;
            if (sub_end > t_end) sub_end = t_end;

            while (i < sub_end) {
                uint32_t w_idx = (uint32_t)(i >> 5);
                uint32_t lane = (uint32_t)(i & 31ULL);
                uint32_t avail = (uint32_t)((32U - lane < sub_end - i) ? (32U - lane) : (sub_end - i));
                uint32_t inv_word = seqX.invalid_bm ? (__ldg(&seqX.invalid_bm[w_idx]) >> lane) : 0U;
                uint32_t avail_mask = (avail == 32U) ? 0xFFFFFFFFU : ((1U << avail) - 1U);
                if ((inv_word & avail_mask) == avail_mask) {
                    valid = 0;
                    kf = 0ULL;
                    kr = 0ULL;
                    act_valid = false;
                    i += avail;
                    continue;
                }
                uint64_t bword = __ldg(&seqX.bits[w_idx]) >> (lane << 1);
                for (uint32_t k = 0U; k < avail; k++) {
                    uint64_t b = bword & 3ULL;
                    bword >>= 2;
                    if (((inv_word >> k) & 1U) == 0U) {
                        kf = ((kf << 2) | b) & mask;
                        kr = (kr >> 2) | ((b ^ 3ULL) << rc_shift);
                        uint64_t cur_i = i + k;
                        if (++valid >= wl && cur_i >= t_start) {
                            uint32_t posX = (uint32_t)(cur_i + 1ULL - (uint64_t)wl);
                            uint64_t sf = kf & eff_mask;
                            uint64_t sr = kr & eff_mask;
                            uint32_t x_strand = (sf <= sr) ? 0U : 1U;
                            uint64_t canon = x_strand ? sr : sf;
                            bool is_pal = (sf == sr);
                            if (stats_buf) valid_q_local++;

                            /* F16: 1-bit ref_unique_bm active-diagonal exact-probe bypass */
                            if constexpr (USE_UNIQUE_BM) {
                                if (act_valid) {
                                    bool ext_ok = false;
                                    uint32_t cand_y = 0U;
                                    if (act_strand == 0U) {
                                        cand_y = act_posY_prev + 1U;
                                        uint64_t new_y = (uint64_t)cand_y + (uint64_t)(wl - 1);
                                        if (new_y < seqY.len &&
                                            !is_bit_set(seqY.invalid_bm, (uint32_t)new_y) &&
                                            base_2bit_dev(seqY.bits, (uint32_t)new_y) == (uint32_t)b) {
                                            ext_ok = true;
                                        }
                                    } else {
                                        if (act_posY_prev > 0U) {
                                            cand_y = act_posY_prev - 1U;
                                            if (!is_bit_set(seqY.invalid_bm, cand_y) &&
                                                (base_2bit_dev(seqY.bits, cand_y) ^ 3U) == (uint32_t)b) {
                                                ext_ok = true;
                                            }
                                        }
                                    }

                                    if (ext_ok) {
                                        act_posY_prev = cand_y;
                                        if (!is_pal && is_bit_set(ref_unique_bm, cand_y)) {
                                            bool emit = skip_lr;
                                            if (!emit) {
                                                if (act_strand == 0U) {
                                                    bool left_ext = !is_bit_set(dropped_bitmap, cand_y - 1U);
                                                    if (!left_ext) {
                                                        emit = true;
                                                    } else {
                                                        bool next_x_valid;
                                                        uint32_t next_x_2bit;
                                                        if (k + 1U < avail) {
                                                            next_x_valid = (((inv_word >> (k + 1U)) & 1U) == 0U);
                                                            next_x_2bit = (uint32_t)(bword & 3ULL);
                                                        } else {
                                                            next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                                                           !is_bit_set(seqX.invalid_bm, posX + (uint32_t)wl);
                                                            next_x_2bit = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 4U;
                                                        }
                                                        bool right_ext = next_x_valid &&
                                                                         ((uint64_t)cand_y + (uint64_t)wl < seqY.len) &&
                                                                         !is_bit_set(seqY.invalid_bm, cand_y + (uint32_t)wl) &&
                                                                         (base_2bit_dev(seqY.bits, cand_y + (uint32_t)wl) == next_x_2bit) &&
                                                                         !is_bit_set(dropped_bitmap, cand_y + 1U);
                                                        emit = !right_ext;
                                                    }
                                                } else {
                                                    bool left_ext = !is_bit_set(dropped_bitmap, cand_y + 1U);
                                                    if (!left_ext) {
                                                        emit = true;
                                                    } else {
                                                        bool next_x_valid;
                                                        uint32_t next_x_2bit;
                                                        if (k + 1U < avail) {
                                                            next_x_valid = (((inv_word >> (k + 1U)) & 1U) == 0U);
                                                            next_x_2bit = (uint32_t)(bword & 3ULL);
                                                        } else {
                                                            next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                                                           !is_bit_set(seqX.invalid_bm, posX + (uint32_t)wl);
                                                            next_x_2bit = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 4U;
                                                        }
                                                        bool right_ext = next_x_valid && (cand_y > 0U) &&
                                                                         !is_bit_set(seqY.invalid_bm, cand_y - 1U) &&
                                                                         ((base_2bit_dev(seqY.bits, cand_y - 1U) ^ 3U) == next_x_2bit) &&
                                                                         !is_bit_set(dropped_bitmap, cand_y - 1U);
                                                        emit = !right_ext;
                                                    }
                                                }
                                            }
                                            if (emit) {
                                                if (act_strand == 0U) {
                                                    stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                                  pack_hit_key_tile(key_layout, 0U, posX, cand_y));
                                                } else {
                                                    uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)cand_y - (uint64_t)wl);
                                                    stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                                  pack_hit_key_tile(key_layout, 1U, posX, posY_rc));
                                                }
                                            }
                                            if (stats_buf) bypassed_local++;
                                            continue;
                                        }
                                    } else {
                                        act_valid = false;
                                    }
                                }
                            }

                            uint32_t st, en, j;
                            uint32_t tag = 0U;
                            if constexpr (DIRECT_CSR) {
                                uint32_t bkt = (uint32_t)canon;
                                st = __ldg(&bucket_off[bkt]);
                                en = __ldg(&bucket_off[bkt + 1U]);
                                if (st == en) continue;
                                j = st;
                            } else {
                                uint64_t h = fmix64(canon);
                                if (part_filter) {
                                    uint32_t p = (uint32_t)(h >> 56);
                                    if (p < part_start || p >= part_end) continue;
                                }
                                uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                                st = __ldg(&bucket_off[bkt]);
                                en = __ldg(&bucket_off[bkt + 1U]);
                                if (st == en) continue;

                                tag = (uint32_t)(h >> tag_shift);
                                j = st;
                                if (en - st <= 16U) {
                                    while (j < en && __ldg(&idx_tags[j]) < tag) j++;
                                } else {
                                    uint32_t lo = st, hi = en;
                                    while (lo < hi) {
                                        uint32_t mid = lo + ((hi - lo) >> 1);
                                        if (__ldg(&idx_tags[mid]) < tag) lo = mid + 1U;
                                        else hi = mid;
                                    }
                                    j = lo;
                                }
                                if (j >= en || __ldg(&idx_tags[j]) != tag) continue;
                            }

                            uint64_t x_seed_lsb = 0ULL;
                            if constexpr (!DIRECT_CSR) {
                                x_seed_lsb = load_32_2bit(seqX.bits, posX) & eff_mask;
                            }

                            bool prev_x_valid = !skip_lr && (posX > 0U) && !is_bit_set(seqX.invalid_bm, posX - 1U);
                            uint32_t prev_x_2bit = prev_x_valid ? base_2bit_dev(seqX.bits, posX - 1U) : 4U;
                            bool next_x_valid = false;
                            uint32_t next_x_2bit = 4U;
                            if (!skip_lr) {
                                if (k + 1U < avail) {
                                    next_x_valid = (((inv_word >> (k + 1U)) & 1U) == 0U);
                                    next_x_2bit = (uint32_t)(bword & 3ULL);
                                } else {
                                    next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                                   !is_bit_set(seqX.invalid_bm, posX + (uint32_t)wl);
                                    next_x_2bit = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 4U;
                                }
                            }

                            while (j < en && (DIRECT_CSR || __ldg(&idx_tags[j]) == tag)) {
                                uint32_t ps = __ldg(&idx_vals[j]);
                                uint32_t posY_orig = ps >> 1;
                                uint32_t y_strand = ps & 1U;
                                uint64_t y_word = 0ULL;
                                if constexpr (!DIRECT_CSR) {
                                    y_word = load_32_2bit(seqY.bits, posY_orig);
                                }

                                if (is_pal || (x_strand == y_strand)) {
                                    bool match_ok = DIRECT_CSR || (((y_word ^ x_seed_lsb) & eff_mask) == 0ULL);
                                    if (match_ok) {
                                        if constexpr (USE_UNIQUE_BM) {
                                            if (!act_valid) {
                                                act_valid = true;
                                                act_strand = 0U;
                                                act_posY_prev = posY_orig;
                                            }
                                        }
                                        bool emit = skip_lr;
                                        if (!emit) {
                                            bool left_ext = prev_x_valid && (posY_orig > 0U) &&
                                                            !is_bit_set(seqY.invalid_bm, posY_orig - 1U) &&
                                                            (base_2bit_dev(seqY.bits, posY_orig - 1U) == prev_x_2bit) &&
                                                            !is_bit_set(dropped_bitmap, posY_orig - 1U);
                                            if (!left_ext) {
                                                emit = true;
                                            } else {
                                                bool right_ext = next_x_valid && ((uint64_t)posY_orig + (uint64_t)wl < seqY.len) &&
                                                                 !is_bit_set(seqY.invalid_bm, posY_orig + (uint32_t)wl) &&
                                                                 (base_2bit_dev(seqY.bits, posY_orig + (uint32_t)wl) == next_x_2bit) &&
                                                                 !is_bit_set(dropped_bitmap, posY_orig + 1U);
                                                emit = !right_ext;
                                            }
                                        }
                                        if (emit) {
                                            stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                          pack_hit_key_tile(key_layout, 0U, posX, posY_orig));
                                        }
                                    }
                                }
                                if (compare_reverse && (is_pal || (x_strand != y_strand))) {
                                    bool match_ok = true;
                                    if constexpr (!DIRECT_CSR) {
                                        uint64_t y_rc_lsb = (~reverse_32_2bit_dev(y_word)) >> rev_rshift;
                                        match_ok = (((y_rc_lsb ^ x_seed_lsb) & eff_mask) == 0ULL);
                                    }
                                    if (match_ok) {
                                        if constexpr (USE_UNIQUE_BM) {
                                            if (!act_valid) {
                                                act_valid = true;
                                                act_strand = 1U;
                                                act_posY_prev = posY_orig;
                                            }
                                        }
                                        uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)posY_orig - (uint64_t)wl);
                                        bool emit = skip_lr;
                                        if (!emit) {
                                            bool left_ext = prev_x_valid && (posY_rc > 0U) &&
                                                            !is_bit_set(seqY.invalid_bm, posY_orig + (uint32_t)wl) &&
                                                            ((base_2bit_dev(seqY.bits, posY_orig + (uint32_t)wl) ^ 3U) == prev_x_2bit) &&
                                                            !is_bit_set(dropped_bitmap, posY_orig + 1U);
                                            if (!left_ext) {
                                                emit = true;
                                            } else {
                                                bool right_ext = next_x_valid && (posY_orig > 0U) &&
                                                                 !is_bit_set(seqY.invalid_bm, posY_orig - 1U) &&
                                                                 ((base_2bit_dev(seqY.bits, posY_orig - 1U) ^ 3U) == next_x_2bit) &&
                                                                 !is_bit_set(dropped_bitmap, posY_orig - 1U);
                                                emit = !right_ext;
                                            }
                                        }
                                        if (emit) {
                                            stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                          pack_hit_key_tile(key_layout, 1U, posX, posY_rc));
                                        }
                                    }
                                }
                                j++;
                            }
                        }
                    } else {
                        valid = 0;
                        kf = 0ULL;
                        kr = 0ULL;
                        act_valid = false;
                    }
                }
                i += avail;
            }
        }

        /* Uniform lockstep block flush (100% CTA-uniform, race-free barriers!) */
        __syncthreads();
        uint32_t raw_cnt = s_cnt;
        if (raw_cnt > 0U) {
            uint32_t flush_cnt = (raw_cnt < GECKO_HIT_STAGE_CAP) ? raw_cnt : GECKO_HIT_STAGE_CAP;
            if (threadIdx.x == 0) {
                if (*(volatile const uint32_t *)hit_count <= max_hits) {
                    s_base = atomicAdd(hit_count, flush_cnt);
                } else {
                    s_base = max_hits;
                }
            }
            __syncthreads();
            uint32_t base = s_base;
            if (base < max_hits) {
                uint32_t rem = max_hits - base;
                uint32_t n_write = (flush_cnt < rem) ? flush_cnt : rem;
                for (uint32_t k = threadIdx.x; k < n_write; k += blockDim.x) {
                    out_hits[base + k] = s_hits[k];
                }
            }
        }
        if (wave + 1U < num_waves) {
            if (raw_cnt > 0U && threadIdx.x == 0) {
                s_cnt = 0U;
            }
            __syncthreads();
        }
    }

    if (stats_buf) {
        if (threadIdx.x == 0) {
            s_valid_q = 0U;
            s_bypass_q = 0U;
        }
        __syncthreads();
        #pragma unroll
        for (int off = 16; off > 0; off >>= 1) {
            valid_q_local  += __shfl_xor_sync(0xFFFFFFFFU, valid_q_local,  off);
            bypassed_local += __shfl_xor_sync(0xFFFFFFFFU, bypassed_local, off);
        }
        if ((threadIdx.x & 31U) == 0U) {
            if (valid_q_local)  atomicAdd(&s_valid_q,  valid_q_local);
            if (bypassed_local) atomicAdd(&s_bypass_q, bypassed_local);
        }
        __syncthreads();
        if (threadIdx.x == 0) {
            if (s_valid_q)  atomicAdd(&stats_buf[0], s_valid_q);
            if (s_bypass_q) atomicAdd(&stats_buf[1], s_bypass_q);
        }
    }
}

/* ---- F16: L2-Resident Wave-Scheduled Hash-Partitioned Probe Kernels ---- */

struct L2PartPlan {
    uint32_t P;          /* 1 = disabled (monolithic), else power of 2 in [2..1024] */
    uint32_t p_bits;     /* log2(P), in [1..10] */
    uint32_t part_shift; /* bkt_bits - p_bits; part = bkt >> part_shift */
};

static inline L2PartPlan choose_l2_part_plan(const GpuRefIndex &idx, int l2_cache_bytes) {
    L2PartPlan plan{1U, 0U, 0U};
    if (idx.is_direct_csr || !idx.d_tags || idx.n_kept == 0U || idx.num_buckets == 0U || l2_cache_bytes <= 0) {
        return plan;
    }
    size_t probe_bytes = ((size_t)idx.num_buckets + 1ULL) * sizeof(uint32_t)
                       + (size_t)idx.n_kept * sizeof(uint32_t);
    size_t target_bytes = ((size_t)l2_cache_bytes * 3ULL) / 4ULL;
    if (probe_bytes <= target_bytes) return plan;

    uint32_t p_bits = 1U;
    while (p_bits < 10U && (probe_bytes >> p_bits) > target_bytes) {
        p_bits++;
    }
    uint32_t bkt_bits = (uint32_t)ceil_log2_u64((uint64_t)idx.num_buckets);
    if (p_bits > bkt_bits) p_bits = bkt_bits;
    if (p_bits == 0U) return plan;

    plan.P = 1U << p_bits;
    plan.p_bits = p_bits;
    plan.part_shift = bkt_bits - p_bits;
    return plan;
}

static __global__ void l2_count_query_parts_kernel(
    GpuSeqView seqX,
    int shift_bits,
    uint32_t part_start,
    uint32_t part_end,
    uint32_t bucket_base,
    uint32_t P,
    uint32_t part_shift,
    int wl,
    uint64_t seed_mask,
    uint64_t tile_start,
    uint64_t tile_end,
    uint32_t * __restrict__ d_part_counts
) {
    __shared__ uint32_t s_hist[1024];
    for (uint32_t p = threadIdx.x; p < P; p += blockDim.x) {
        s_hist[p] = 0U;
    }
    __syncthreads();

    const uint64_t STRIDE = 64ULL;
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t t_start = tile_start + gtid * STRIDE;
    if (t_start < tile_end) {
        uint64_t t_end = t_start + STRIDE;
        if (t_end > tile_end) t_end = tile_end;

        uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
        uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : mask;
        bool part_filter = (part_start != 0U || part_end != 256U);
        int rc_shift = 2 * (wl - 1);
        uint64_t i = (t_start >= (uint64_t)(wl - 1)) ? (t_start - (uint64_t)(wl - 1)) : 0ULL;

        uint64_t kf = 0ULL, kr = 0ULL;
        int valid = 0;

        while (i < t_end) {
            uint32_t w_idx = (uint32_t)(i >> 5);
            uint32_t lane = (uint32_t)(i & 31ULL);
            uint32_t avail = (uint32_t)((32U - lane < t_end - i) ? (32U - lane) : (t_end - i));
            uint32_t inv_word = seqX.invalid_bm ? (__ldg(&seqX.invalid_bm[w_idx]) >> lane) : 0U;
            uint32_t avail_mask = (avail == 32U) ? 0xFFFFFFFFU : ((1U << avail) - 1U);
            if ((inv_word & avail_mask) == avail_mask) {
                valid = 0;
                kf = 0ULL;
                kr = 0ULL;
                i += avail;
                continue;
            }
            uint64_t bword = __ldg(&seqX.bits[w_idx]) >> (lane << 1);
            for (uint32_t k = 0U; k < avail; k++) {
                uint64_t b = bword & 3ULL;
                bword >>= 2;
                if (((inv_word >> k) & 1U) == 0U) {
                    kf = ((kf << 2) | b) & mask;
                    kr = (kr >> 2) | ((b ^ 3ULL) << rc_shift);
                    uint64_t cur_i = i + k;
                    if (++valid >= wl && cur_i >= t_start) {
                        uint64_t sf = kf & eff_mask;
                        uint64_t sr = kr & eff_mask;
                        uint64_t canon = (sf <= sr) ? sf : sr;
                        uint64_t h = fmix64(canon);
                        if (part_filter) {
                            uint32_t cp = (uint32_t)(h >> 56);
                            if (cp < part_start || cp >= part_end) continue;
                        }
                        uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                        uint32_t part = bkt >> part_shift;
                        atomicAdd(&s_hist[part], 1U);
                    }
                } else {
                    valid = 0;
                    kf = 0ULL;
                    kr = 0ULL;
                }
            }
            i += avail;
        }
    }

    __syncthreads();
    for (uint32_t p = threadIdx.x; p < P; p += blockDim.x) {
        uint32_t c = s_hist[p];
        if (c > 0U) {
            atomicAdd(&d_part_counts[p], c);
        }
    }
}

static __global__ void l2_prefix_sum_parts_kernel(
    const uint32_t * __restrict__ d_part_counts,
    uint32_t P,
    uint32_t * __restrict__ d_part_offsets,
    uint32_t * __restrict__ d_part_cursor,
    uint32_t * __restrict__ d_part_block_off,
    uint32_t * __restrict__ d_totals
) {
    using BlockScan = cub::BlockScan<uint32_t, 256>;
    __shared__ typename BlockScan::TempStorage temp_storage;

    uint32_t cnt[4];
    uint32_t blk[4];
    #pragma unroll
    for (uint32_t k = 0U; k < 4U; k++) {
        uint32_t p = threadIdx.x * 4U + k;
        uint32_t c = (p < P) ? d_part_counts[p] : 0U;
        cnt[k] = c;
        blk[k] = (c + 255U) >> 8;
    }

    uint32_t off[4];
    uint32_t total_q = 0U;
    BlockScan(temp_storage).ExclusiveSum(cnt, off, total_q);
    #pragma unroll
    for (uint32_t k = 0U; k < 4U; k++) {
        uint32_t p = threadIdx.x * 4U + k;
        if (p < P) {
            d_part_offsets[p] = off[k];
            d_part_cursor[p] = off[k];
        }
    }
    if (threadIdx.x == 0) {
        d_part_offsets[P] = total_q;
        d_totals[0] = total_q;
    }

    __syncthreads();

    uint32_t boff[4];
    uint32_t total_blks = 0U;
    BlockScan(temp_storage).ExclusiveSum(blk, boff, total_blks);
    #pragma unroll
    for (uint32_t k = 0U; k < 4U; k++) {
        uint32_t p = threadIdx.x * 4U + k;
        if (p < P) {
            d_part_block_off[p] = boff[k];
        }
    }
    if (threadIdx.x == 0) {
        d_part_block_off[P] = total_blks;
        d_totals[1] = total_blks;
    }
}

static __global__ void l2_scatter_query_parts_kernel(
    GpuSeqView seqX,
    int shift_bits,
    uint32_t part_start,
    uint32_t part_end,
    uint32_t bucket_base,
    uint32_t part_shift,
    int wl,
    uint64_t seed_mask,
    uint64_t tile_start,
    uint64_t tile_end,
    HitKeyLayout key_layout,
    uint32_t * __restrict__ d_part_cursor,
    uint64_t * __restrict__ d_q_hash,
    uint32_t * __restrict__ d_q_meta
) {
    const uint64_t STRIDE = 64ULL;
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t t_start = tile_start + gtid * STRIDE;
    if (t_start >= tile_end) return;
    uint64_t t_end = t_start + STRIDE;
    if (t_end > tile_end) t_end = tile_end;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : mask;
    bool part_filter = (part_start != 0U || part_end != 256U);
    int rc_shift = 2 * (wl - 1);
    uint64_t i = (t_start >= (uint64_t)(wl - 1)) ? (t_start - (uint64_t)(wl - 1)) : 0ULL;

    uint64_t kf = 0ULL, kr = 0ULL;
    int valid = 0;

    while (i < t_end) {
        uint32_t w_idx = (uint32_t)(i >> 5);
        uint32_t lane = (uint32_t)(i & 31ULL);
        uint32_t avail = (uint32_t)((32U - lane < t_end - i) ? (32U - lane) : (t_end - i));
        uint32_t inv_word = seqX.invalid_bm ? (__ldg(&seqX.invalid_bm[w_idx]) >> lane) : 0U;
        uint32_t avail_mask = (avail == 32U) ? 0xFFFFFFFFU : ((1U << avail) - 1U);
        if ((inv_word & avail_mask) == avail_mask) {
            valid = 0;
            kf = 0ULL;
            kr = 0ULL;
            i += avail;
            continue;
        }
        uint64_t bword = __ldg(&seqX.bits[w_idx]) >> (lane << 1);
        for (uint32_t k = 0U; k < avail; k++) {
            uint64_t b = bword & 3ULL;
            bword >>= 2;
            if (((inv_word >> k) & 1U) == 0U) {
                kf = ((kf << 2) | b) & mask;
                kr = (kr >> 2) | ((b ^ 3ULL) << rc_shift);
                uint64_t cur_i = i + k;
                if (++valid >= wl && cur_i >= t_start) {
                    uint32_t posX = (uint32_t)(cur_i + 1ULL - (uint64_t)wl);
                    uint64_t sf = kf & eff_mask;
                    uint64_t sr = kr & eff_mask;
                    uint32_t x_strand = (sf <= sr) ? 0U : 1U;
                    uint64_t canon = x_strand ? sr : sf;
                    bool is_pal = (sf == sr);
                    uint64_t h = fmix64(canon);
                    if (part_filter) {
                        uint32_t cp = (uint32_t)(h >> 56);
                        if (cp < part_start || cp >= part_end) continue;
                    }
                    uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                    uint32_t part = bkt >> part_shift;
                    uint32_t slot = atomicAdd(&d_part_cursor[part], 1U);
                    uint32_t rel_x = posX - key_layout.min_x;
                    d_q_hash[slot] = h;
                    d_q_meta[slot] = (rel_x << 2) | ((uint32_t)is_pal << 1) | x_strand;
                }
            } else {
                valid = 0;
                kf = 0ULL;
                kr = 0ULL;
            }
        }
        i += avail;
    }
}

static __global__ void l2_probe_partitioned_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int shift_bits,
    int tag_shift,
    uint32_t part_start,
    uint32_t part_end,
    uint32_t bucket_base,
    uint32_t P,
    int wl,
    uint64_t seed_mask,
    bool compare_reverse,
    bool has_global_dropped_bm,
    HitKeyLayout key_layout,
    const uint32_t * __restrict__ d_part_offsets,
    const uint32_t * __restrict__ d_part_block_off,
    const uint64_t * __restrict__ d_q_hash,
    const uint32_t * __restrict__ d_q_meta,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits
) {
    __shared__ uint64_t s_hits[GECKO_HIT_STAGE_CAP];
    __shared__ uint32_t s_cnt;
    __shared__ uint32_t s_base;
    __shared__ uint32_t s_q_start;
    __shared__ uint32_t s_q_end;

    if (threadIdx.x == 0) {
        uint32_t bid = blockIdx.x;
        uint32_t lo = 0U, hi = P;
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (__ldg(&d_part_block_off[mid + 1U]) <= bid) lo = mid + 1U;
            else hi = mid;
        }
        uint32_t part = lo;
        uint32_t local_blk = bid - __ldg(&d_part_block_off[part]);
        uint32_t pq_start = __ldg(&d_part_offsets[part]);
        uint32_t pq_end   = __ldg(&d_part_offsets[part + 1U]);
        uint32_t qs = pq_start + (local_blk << 8);
        s_q_start = qs;
        s_q_end   = (qs + 256U < pq_end) ? (qs + 256U) : pq_end;
        s_cnt = 0U;
    }
    __syncthreads();

    uint32_t q_idx = s_q_start + threadIdx.x;
    if (q_idx < s_q_end) {
        uint64_t h = d_q_hash[q_idx];
        uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
        uint32_t st = __ldg(&bucket_off[bkt]);
        uint32_t en = __ldg(&bucket_off[bkt + 1U]);
        if (st != en) {
            uint32_t tag = (uint32_t)(h >> tag_shift);
            uint32_t j = st;
            if (en - st <= 16U) {
                while (j < en && __ldg(&idx_tags[j]) < tag) j++;
            } else {
                uint32_t lo = st, hi = en;
                while (lo < hi) {
                    uint32_t mid = lo + ((hi - lo) >> 1);
                    if (__ldg(&idx_tags[mid]) < tag) lo = mid + 1U;
                    else hi = mid;
                }
                j = lo;
            }
            if (j < en && __ldg(&idx_tags[j]) == tag) {
                uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
                uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : mask;
                bool is_spaced = (eff_mask != mask);
                bool part_filter = (part_start != 0U || part_end != 256U);
                bool skip_lr = is_spaced || (part_filter && !has_global_dropped_bm);
                int rev_rshift = 64 - 2 * wl;

                uint32_t meta = d_q_meta[q_idx];
                uint32_t posX = key_layout.min_x + (meta >> 2);
                bool is_pal = ((meta >> 1) & 1U) != 0U;
                uint32_t x_strand = meta & 1U;

                uint64_t x_seed_lsb = load_32_2bit(seqX.bits, posX) & eff_mask;
                bool prev_x_valid = !skip_lr && (posX > 0U) && !is_bit_set(seqX.invalid_bm, posX - 1U);
                uint32_t prev_x_2bit = prev_x_valid ? base_2bit_dev(seqX.bits, posX - 1U) : 4U;
                bool next_x_valid = !skip_lr && ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                    !is_bit_set(seqX.invalid_bm, posX + (uint32_t)wl);
                uint32_t next_x_2bit = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 4U;

                while (j < en && __ldg(&idx_tags[j]) == tag) {
                    uint32_t ps = __ldg(&idx_vals[j]);
                    uint32_t posY_orig = ps >> 1;
                    uint32_t y_strand = ps & 1U;
                    uint64_t y_word = load_32_2bit(seqY.bits, posY_orig);

                    if (is_pal || (x_strand == y_strand)) {
                        if (((y_word ^ x_seed_lsb) & eff_mask) == 0ULL) {
                            bool emit = skip_lr;
                            if (!emit) {
                                bool left_ext = prev_x_valid && (posY_orig > 0U) &&
                                                !is_bit_set(seqY.invalid_bm, posY_orig - 1U) &&
                                                (base_2bit_dev(seqY.bits, posY_orig - 1U) == prev_x_2bit) &&
                                                !is_bit_set(dropped_bitmap, posY_orig - 1U);
                                if (!left_ext) {
                                    emit = true;
                                } else {
                                    bool right_ext = next_x_valid && ((uint64_t)posY_orig + (uint64_t)wl < seqY.len) &&
                                                     !is_bit_set(seqY.invalid_bm, posY_orig + (uint32_t)wl) &&
                                                     (base_2bit_dev(seqY.bits, posY_orig + (uint32_t)wl) == next_x_2bit) &&
                                                     !is_bit_set(dropped_bitmap, posY_orig + 1U);
                                    emit = !right_ext;
                                }
                            }
                            if (emit) {
                                stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                              pack_hit_key_tile(key_layout, 0U, posX, posY_orig));
                            }
                        }
                    }
                    if (compare_reverse && (is_pal || (x_strand != y_strand))) {
                        uint64_t y_rc_lsb = (~reverse_32_2bit_dev(y_word)) >> rev_rshift;
                        if (((y_rc_lsb ^ x_seed_lsb) & eff_mask) == 0ULL) {
                            uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)posY_orig - (uint64_t)wl);
                            bool emit = skip_lr;
                            if (!emit) {
                                bool left_ext = prev_x_valid && (posY_rc > 0U) &&
                                                !is_bit_set(seqY.invalid_bm, posY_orig + (uint32_t)wl) &&
                                                ((base_2bit_dev(seqY.bits, posY_orig + (uint32_t)wl) ^ 3U) == prev_x_2bit) &&
                                                !is_bit_set(dropped_bitmap, posY_orig + 1U);
                                if (!left_ext) {
                                    emit = true;
                                } else {
                                    bool right_ext = next_x_valid && (posY_orig > 0U) &&
                                                     !is_bit_set(seqY.invalid_bm, posY_orig - 1U) &&
                                                     ((base_2bit_dev(seqY.bits, posY_orig - 1U) ^ 3U) == next_x_2bit) &&
                                                     !is_bit_set(dropped_bitmap, posY_orig - 1U);
                                    emit = !right_ext;
                                }
                            }
                            if (emit) {
                                stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                              pack_hit_key_tile(key_layout, 1U, posX, posY_rc));
                            }
                        }
                    }
                    j++;
                }
            }
        }
    }

    __syncthreads();
    uint32_t raw_cnt = s_cnt;
    if (raw_cnt > 0U) {
        uint32_t flush_cnt = (raw_cnt < GECKO_HIT_STAGE_CAP) ? raw_cnt : GECKO_HIT_STAGE_CAP;
        if (threadIdx.x == 0) {
            if (*(volatile const uint32_t *)hit_count <= max_hits) {
                s_base = atomicAdd(hit_count, flush_cnt);
            } else {
                s_base = max_hits;
            }
        }
        __syncthreads();
        uint32_t base = s_base;
        if (base < max_hits) {
            uint32_t rem = max_hits - base;
            uint32_t n_write = (flush_cnt < rem) ? flush_cnt : rem;
            for (uint32_t k = threadIdx.x; k < n_write; k += blockDim.x) {
                out_hits[base + k] = s_hits[k];
            }
        }
    }
}

static __global__ void filter_sorted_hits_kernel(
    const uint64_t * __restrict__ sorted_hits,
    uint32_t n,
    HitKeyLayout key_layout,
    uint32_t min_gap,
    uint8_t * __restrict__ keep_flags
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (i == 0U) {
        keep_flags[0] = 1;
        return;
    }
    uint64_t cur = sorted_hits[i];
    uint64_t prev = sorted_hits[i - 1U];
    if (hit_strand_diag_tile(key_layout, cur) != hit_strand_diag_tile(key_layout, prev)) {
        keep_flags[i] = 1;
    } else {
        uint32_t cur_rx  = (uint32_t)(cur  & (uint64_t)key_layout.x_mask);
        uint32_t prev_rx = (uint32_t)(prev & (uint64_t)key_layout.x_mask);
        keep_flags[i] = (cur_rx > prev_rx + min_gap) ? 1 : 0;
    }
}

uint32_t gpu_query_tile_left_maximal(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuRefIndex &idx,
    int wl,
    uint64_t seed_mask,
    bool compare_reverse,
    uint64_t tile_start,
    uint64_t tile_end,
    HitKeyLayout key_layout,
    int l2_cache_bytes,
    int query_mode,
    uint32_t *d_stats_buf,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint8_t *d_hit_flags,
    uint32_t *d_hit_count,
    uint32_t *d_filt_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed,
    uint64_t *out_valid_q,
    uint64_t *out_bypassed_q
) {
    out_overflowed = false;
    out_raw_leftmax = 0U;
    if (out_valid_q) *out_valid_q = 0ULL;
    if (out_bypassed_q) *out_bypassed_q = 0ULL;
    if (idx.n_kept == 0U || tile_start >= tile_end) return 0U;

    CUDA_CHECK_VOID(cudaMemset(d_hit_count, 0, sizeof(uint32_t)));

    uint64_t tile_len = tile_end - tile_start;
    bool ran_l2_part = false;

    /* F16: L2-Resident Wave-Scheduled Hash-Partitioned Probe (aliases idle d_hits_sorted) */
    if (query_mode == GECKO_QMODE_L2_PART) {
        L2PartPlan plan = choose_l2_part_plan(idx, l2_cache_bytes);
        if (plan.P >= 2U &&
            tile_len <= (1ULL << 24) &&
            tile_len >= (uint64_t)plan.P * 4096ULL) {
            auto align256 = [](size_t n) -> size_t { return (n + 255ULL) & ~255ULL; };
            size_t max_q = (size_t)tile_len;
            size_t off_hash   = 0ULL;
            size_t off_meta   = align256(off_hash   + max_q * sizeof(uint64_t));
            size_t off_counts = align256(off_meta   + max_q * sizeof(uint32_t));
            size_t off_offs   = align256(off_counts + (size_t)plan.P * sizeof(uint32_t));
            size_t off_cursor = align256(off_offs   + ((size_t)plan.P + 1ULL) * sizeof(uint32_t));
            size_t off_blkoff = align256(off_cursor + (size_t)plan.P * sizeof(uint32_t));
            size_t off_totals = align256(off_blkoff + ((size_t)plan.P + 1ULL) * sizeof(uint32_t));
            size_t total_need = align256(off_totals + 2ULL * sizeof(uint32_t));
            size_t scratch_bytes = (size_t)max_hits * sizeof(uint64_t);

            if (total_need <= scratch_bytes) {
                uint8_t *base_ptr = reinterpret_cast<uint8_t *>(d_hits_sorted);
                uint64_t *d_q_hash         = reinterpret_cast<uint64_t *>(base_ptr + off_hash);
                uint32_t *d_q_meta         = reinterpret_cast<uint32_t *>(base_ptr + off_meta);
                uint32_t *d_part_counts    = reinterpret_cast<uint32_t *>(base_ptr + off_counts);
                uint32_t *d_part_offsets   = reinterpret_cast<uint32_t *>(base_ptr + off_offs);
                uint32_t *d_part_cursor    = reinterpret_cast<uint32_t *>(base_ptr + off_cursor);
                uint32_t *d_part_block_off = reinterpret_cast<uint32_t *>(base_ptr + off_blkoff);
                uint32_t *d_totals         = reinterpret_cast<uint32_t *>(base_ptr + off_totals);

                CUDA_CHECK_VOID(cudaMemset(d_part_counts, 0, (size_t)plan.P * sizeof(uint32_t)));
                uint64_t n_threads_p1 = (tile_len + 63ULL) / 64ULL;
                uint32_t grid_p1 = (uint32_t)((n_threads_p1 + 255ULL) / 256ULL);

                l2_count_query_parts_kernel<<<grid_p1, 256>>>(
                    seqX, idx.shift_bits, idx.part_start, idx.part_end, idx.bucket_base,
                    plan.P, plan.part_shift, wl, seed_mask, tile_start, tile_end, d_part_counts);
                CUDA_CHECK_VOID(cudaGetLastError());

                l2_prefix_sum_parts_kernel<<<1, 256>>>(
                    d_part_counts, plan.P, d_part_offsets, d_part_cursor, d_part_block_off, d_totals);
                CUDA_CHECK_VOID(cudaGetLastError());

                uint32_t h_totals[2] = {0U, 0U};
                CUDA_CHECK_VOID(cudaMemcpy(h_totals, d_totals, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                if (h_totals[0] == 0U || h_totals[1] == 0U) {
                    return 0U;
                }

                l2_scatter_query_parts_kernel<<<grid_p1, 256>>>(
                    seqX, idx.shift_bits, idx.part_start, idx.part_end, idx.bucket_base,
                    plan.part_shift, wl, seed_mask, tile_start, tile_end, key_layout,
                    d_part_cursor, d_q_hash, d_q_meta);
                CUDA_CHECK_VOID(cudaGetLastError());

                l2_probe_partitioned_kernel<<<h_totals[1], 256>>>(
                    seqX, seqY,
                    idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap,
                    idx.shift_bits, idx.tag_shift, idx.part_start, idx.part_end, idx.bucket_base,
                    plan.P, wl, seed_mask, compare_reverse, idx.has_global_dropped_bm,
                    key_layout, d_part_offsets, d_part_block_off, d_q_hash, d_q_meta,
                    d_hits, d_hit_count, max_hits);
                CUDA_CHECK_VOID(cudaGetLastError());
                ran_l2_part = true;
            }
        }
    }

    if (!ran_l2_part) {
        uint32_t num_waves = (tile_len >= (512ULL << 10)) ? 4U : 1U;
        uint64_t stride_bp = (uint64_t)num_waves * 16ULL;
        uint64_t n_threads = (tile_len + stride_bp - 1ULL) / stride_bp;
        uint32_t grid = (uint32_t)((n_threads + 255ULL) / 256ULL);

        uint32_t *use_stats = (out_valid_q && out_bypassed_q && d_stats_buf) ? d_stats_buf : nullptr;
        if (use_stats) {
            CUDA_CHECK_VOID(cudaMemset(use_stats, 0, 2 * sizeof(uint32_t)));
        }

        if (idx.is_direct_csr) {
            query_left_maximal_kernel<true, false><<<grid, 256>>>(
                seqX, seqY,
                idx.d_bucket_off, nullptr, idx.d_vals, idx.d_dropped_bitmap, nullptr,
                0, 0, 0U, 256U, 0U,
                wl, seed_mask, compare_reverse, idx.has_global_dropped_bm,
                tile_start, tile_end, num_waves, key_layout,
                d_hits, d_hit_count, max_hits, use_stats);
        } else if (idx.d_ref_unique_bm != nullptr) {
            query_left_maximal_kernel<false, true><<<grid, 256>>>(
                seqX, seqY,
                idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap, idx.d_ref_unique_bm,
                idx.shift_bits, idx.tag_shift, idx.part_start, idx.part_end, idx.bucket_base,
                wl, seed_mask, compare_reverse, idx.has_global_dropped_bm,
                tile_start, tile_end, num_waves, key_layout,
                d_hits, d_hit_count, max_hits, use_stats);
        } else {
            query_left_maximal_kernel<false, false><<<grid, 256>>>(
                seqX, seqY,
                idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap, nullptr,
                idx.shift_bits, idx.tag_shift, idx.part_start, idx.part_end, idx.bucket_base,
                wl, seed_mask, compare_reverse, idx.has_global_dropped_bm,
                tile_start, tile_end, num_waves, key_layout,
                d_hits, d_hit_count, max_hits, use_stats);
        }
        CUDA_CHECK_VOID(cudaGetLastError());

        if (use_stats) {
            uint32_t h_stats[2] = {0U, 0U};
            CUDA_CHECK_VOID(cudaMemcpy(h_stats, use_stats, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
            *out_valid_q = (uint64_t)h_stats[0];
            *out_bypassed_q = (uint64_t)h_stats[1];
        }
    }

    uint32_t h_nhits = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nhits, d_hit_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nhits > max_hits) {
        out_overflowed = true;
        return 0U;
    }
    out_raw_leftmax = h_nhits;
    if (h_nhits <= 1U) return h_nhits;

    /* F14: Tile-relative bit-bounded CUB SortKeys (saves 1-4 full radix passes) */
    int end_bit = (int)key_layout.end_bit;
    size_t sort_bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, sort_bytes, d_hits, d_hits_sorted, (int)h_nhits, 0, end_bit);
    size_t sel_bytes = 0;
    cub::DeviceSelect::Flagged(nullptr, sel_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);
    size_t need_bytes = (sort_bytes > sel_bytes) ? sort_bytes : sel_bytes;
    gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, need_bytes);

    cub::DeviceRadixSort::SortKeys(d_cub_temp, sort_bytes, d_hits, d_hits_sorted, (int)h_nhits, 0, end_bit);

    filter_sorted_hits_kernel<<<(h_nhits + 255U) / 256U, 256>>>(
        d_hits_sorted, h_nhits, key_layout, (uint32_t)(wl / 2), d_hit_flags);
    CUDA_CHECK_VOID(cudaGetLastError());

    cub::DeviceSelect::Flagged(d_cub_temp, sel_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);

    uint32_t h_nfilt = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nfilt, d_filt_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    return h_nfilt;
}
