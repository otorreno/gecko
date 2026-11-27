/*
 * cuda_query.cu
 *
 * Stage 2: Tiled 2-bit GPU Query + Lock-Free Left/Right-Maximal Seed Filtering:
 *  - F14: Direct-CSR query path (wl <= 12, zero fmix64/d_tags/2-bit verify),
 *    block-level shared-memory hit staging (STAGE_CAP=512, lockstep wave flush,
 *    lock-free overflow spill), lazy right_ext evaluation, and tile-relative
 *    bit-bounded CUB SortKeys via HitKeyLayout.
 *  - F15: Multi-slab skip_lr fix via global dropped_bitmap / freq_limit == 0.
 *  - F16 + W3-F13a/b: Inverted 1-bit ref_non_unique_bm active-diagonal bypass,
 *    hybrid anchor-seeded pre-bypass (l2_pre_bypass_kernel + d_q_bypass_bm), and
 *    two-phase CTA-privatized shared-memory L2-partitioned scatter/probe.
 *  - W3-F11b: Lossless Oracle parity (removed filter_sorted_hits_kernel min_gap=wl/2
 *    heuristic + DeviceSelect::Flagged; uses CUB DoubleBuffer SortKeys pointer swap).
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
    const uint32_t * __restrict__ ref_non_unique_bm,
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
                                        /* W3-F13a: Inverted polarity — USE_UNIQUE_BM guarantees ref_non_unique_bm != nullptr */
                                        if (!is_pal && !is_bit_set(ref_non_unique_bm, cand_y)) {
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

/* ---- W3-F13b: Hybrid Anchor-Seeded Pre-Bypass Kernel before L2_PART ----
 * Each thread owns 4 consecutive 32-bp words (128 bp) of d_q_bypass_bm:
 *  - Zero atomicOr and zero cudaMemset on d_q_bypass_bm (1 coalesced u32 write/wave).
 *  - Probes 1 anchor every A=16 bp when !act_valid; once on a unique active diagonal,
 *    bypasses consecutive matches using only seqY + ref_non_unique_bm + dropped_bitmap
 *    (0 Hash-CSR probes, 0 rolling kf/kr).
 *  - Sets d_q_bypass_bm bit = 1 for: invalid seqX k-mers, anchor misses, and unique
 *    non-palindromic hits emitted directly here. Leaves bit = 0 for all remaining
 *    queries to be handled cleanly by L2_PART. */
static __global__ void l2_pre_bypass_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    const uint32_t * __restrict__ ref_non_unique_bm,
    int shift_bits,
    int tag_shift,
    uint32_t bucket_base,
    int wl,
    bool compare_reverse,
    uint32_t num_pos,
    uint32_t bp_words,
    HitKeyLayout key_layout,
    uint32_t * __restrict__ d_q_bypass_bm,
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

    uint32_t gtid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t w_base = gtid * 4U;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint32_t valid_mask_wl = (wl == 32) ? 0xFFFFFFFFU : ((1U << wl) - 1U);
    int rev_rshift = 64 - 2 * wl;

    bool act_valid = false;
    uint32_t act_strand = 0U;
    uint32_t act_posY_prev = 0U;
    uint32_t valid_q_local = 0U;
    uint32_t bypassed_local = 0U;

    for (uint32_t wave = 0U; wave < 4U; wave++) {
        uint32_t w_out = w_base + wave;
        if (w_out < bp_words) {
            uint32_t w_rel0 = w_out << 5;
            uint32_t w_cnt = (num_pos - w_rel0 < 32U) ? (num_pos - w_rel0) : 32U;
            uint32_t pos0 = key_layout.min_x + w_rel0;

            uint32_t inv_lo = seqX.invalid_bm ? load_32_1bit(seqX.invalid_bm, pos0) : 0U;
            uint32_t inv_hi = (seqX.invalid_bm && ((uint64_t)pos0 + 32ULL < seqX.len))
                                  ? load_32_1bit(seqX.invalid_bm, pos0 + 32U)
                                  : 0U;
            uint64_t inv64 = ((uint64_t)inv_hi << 32) | (uint64_t)inv_lo;

            uint32_t bp_word = (w_cnt < 32U) ? (~0U << w_cnt) : 0U;

            for (uint32_t s = 0U; s < w_cnt; s++) {
                uint32_t posX = pos0 + s;
                /* Case 1: Invalid k-mer in seqX — mark done so L2_PART skips completely */
                if (((uint32_t)(inv64 >> s) & valid_mask_wl) != 0U) {
                    act_valid = false;
                    bp_word |= (1U << s);
                    continue;
                }
                if (stats_buf) valid_q_local++;

                /* Case 2 & 3: Active-diagonal 1-bp extension */
                if (act_valid) {
                    uint32_t b = base_2bit_dev(seqX.bits, posX + (uint32_t)(wl - 1));
                    bool ext_ok = false;
                    uint32_t cand_y = 0U;
                    if (act_strand == 0U) {
                        cand_y = act_posY_prev + 1U;
                        uint64_t new_y = (uint64_t)cand_y + (uint64_t)(wl - 1);
                        if (new_y < seqY.len &&
                            !is_bit_set(seqY.invalid_bm, (uint32_t)new_y) &&
                            base_2bit_dev(seqY.bits, (uint32_t)new_y) == b) {
                            ext_ok = true;
                        }
                    } else {
                        if (act_posY_prev > 0U) {
                            cand_y = act_posY_prev - 1U;
                            if (!is_bit_set(seqY.invalid_bm, cand_y) &&
                                (base_2bit_dev(seqY.bits, cand_y) ^ 3U) == b) {
                                ext_ok = true;
                            }
                        }
                    }

                    if (ext_ok) {
                        act_posY_prev = cand_y;
                        /* ref_non_unique_bm already marks all run_len >= 2 AND all palindromes with 1! */
                        if (!is_bit_set(ref_non_unique_bm, cand_y)) {
                            /* Case 2: Unique non-palindromic kept match on active diagonal */
                            bool emit = false;
                            if (act_strand == 0U) {
                                bool left_ext = !is_bit_set(dropped_bitmap, cand_y - 1U);
                                if (!left_ext) {
                                    emit = true;
                                } else {
                                    bool next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                                        (((inv64 >> (s + (uint32_t)wl)) & 1ULL) == 0ULL);
                                    bool right_ext = next_x_valid &&
                                                     ((uint64_t)cand_y + (uint64_t)wl < seqY.len) &&
                                                     !is_bit_set(seqY.invalid_bm, cand_y + (uint32_t)wl) &&
                                                     (base_2bit_dev(seqY.bits, cand_y + (uint32_t)wl) ==
                                                      base_2bit_dev(seqX.bits, posX + (uint32_t)wl)) &&
                                                     !is_bit_set(dropped_bitmap, cand_y + 1U);
                                    emit = !right_ext;
                                }
                                if (emit) {
                                    stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                  pack_hit_key_tile(key_layout, 0U, posX, cand_y));
                                }
                            } else {
                                bool left_ext = !is_bit_set(dropped_bitmap, cand_y + 1U);
                                if (!left_ext) {
                                    emit = true;
                                } else {
                                    bool next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                                        (((inv64 >> (s + (uint32_t)wl)) & 1ULL) == 0ULL);
                                    bool right_ext = next_x_valid && (cand_y > 0U) &&
                                                     !is_bit_set(seqY.invalid_bm, cand_y - 1U) &&
                                                     ((base_2bit_dev(seqY.bits, cand_y - 1U) ^ 3U) ==
                                                      base_2bit_dev(seqX.bits, posX + (uint32_t)wl)) &&
                                                     !is_bit_set(dropped_bitmap, cand_y - 1U);
                                    emit = !right_ext;
                                }
                                if (emit) {
                                    uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)cand_y - (uint64_t)wl);
                                    stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                  pack_hit_key_tile(key_layout, 1U, posX, posY_rc));
                                }
                            }
                            bp_word |= (1U << s);
                            if (stats_buf) bypassed_local++;
                        }
                        /* Case 3 (non-unique or palindrome): keep act_valid=true, leave bp_word bit 0 */
                        continue;
                    }
                    act_valid = false;
                }

                /* Case 4: Strided anchor probe every A = 16 bp when !act_valid */
                if ((s & 15U) == 0U) {
                    uint64_t w_lsb = load_32_2bit(seqX.bits, posX);
                    uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift;
                    uint64_t kr = (~w_lsb) & mask;
                    uint32_t x_strand = (kf <= kr) ? 0U : 1U;
                    uint64_t canon = x_strand ? kr : kf;
                    bool is_pal = (kf == kr);
                    uint64_t h = fmix64(canon);

                    uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                    uint32_t st = __ldg(&bucket_off[bkt]);
                    uint32_t en = __ldg(&bucket_off[bkt + 1U]);

                    bool found = false;
                    uint32_t m_posY = 0U;
                    uint32_t m_strand = 0U;

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
                            uint64_t x_seed_lsb = w_lsb & mask;
                            while (j < en && __ldg(&idx_tags[j]) == tag) {
                                uint32_t ps = __ldg(&idx_vals[j]);
                                uint32_t posY_orig = ps >> 1;
                                uint32_t y_strand = ps & 1U;
                                uint64_t y_word = load_32_2bit(seqY.bits, posY_orig);
                                if (is_pal || (x_strand == y_strand)) {
                                    if (((y_word ^ x_seed_lsb) & mask) == 0ULL) {
                                        found = true;
                                        m_posY = posY_orig;
                                        m_strand = 0U;
                                        break;
                                    }
                                }
                                if (compare_reverse && (is_pal || (x_strand != y_strand))) {
                                    uint64_t y_rc_lsb = (~reverse_32_2bit_dev(y_word)) >> rev_rshift;
                                    if (((y_rc_lsb ^ x_seed_lsb) & mask) == 0ULL) {
                                        found = true;
                                        m_posY = posY_orig;
                                        m_strand = 1U;
                                        break;
                                    }
                                }
                                j++;
                            }
                        }
                    }

                    if (!found) {
                        /* Subcase 4a: Proven anchor miss — skip in L2_PART! */
                        bp_word |= (1U << s);
                    } else {
                        /* Seed active diagonal for posX + 1 */
                        act_valid = true;
                        act_strand = m_strand;
                        act_posY_prev = m_posY;

                        if (!is_pal && !is_bit_set(ref_non_unique_bm, m_posY)) {
                            /* Subcase 4b: Unique non-palindromic anchor match — FULL left_ext & right_ext! */
                            bool prev_x_valid = (posX > 0U) && !is_bit_set(seqX.invalid_bm, posX - 1U);
                            uint32_t prev_x_2bit = prev_x_valid ? base_2bit_dev(seqX.bits, posX - 1U) : 4U;
                            bool next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                                (((inv64 >> (s + (uint32_t)wl)) & 1ULL) == 0ULL);
                            uint32_t next_x_2bit = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 4U;

                            if (m_strand == 0U) {
                                bool left_ext = prev_x_valid && (m_posY > 0U) &&
                                                !is_bit_set(seqY.invalid_bm, m_posY - 1U) &&
                                                (base_2bit_dev(seqY.bits, m_posY - 1U) == prev_x_2bit) &&
                                                !is_bit_set(dropped_bitmap, m_posY - 1U);
                                bool emit = !left_ext;
                                if (!emit) {
                                    bool right_ext = next_x_valid && ((uint64_t)m_posY + (uint64_t)wl < seqY.len) &&
                                                     !is_bit_set(seqY.invalid_bm, m_posY + (uint32_t)wl) &&
                                                     (base_2bit_dev(seqY.bits, m_posY + (uint32_t)wl) == next_x_2bit) &&
                                                     !is_bit_set(dropped_bitmap, m_posY + 1U);
                                    emit = !right_ext;
                                }
                                if (emit) {
                                    stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                  pack_hit_key_tile(key_layout, 0U, posX, m_posY));
                                }
                            } else {
                                uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)m_posY - (uint64_t)wl);
                                bool left_ext = prev_x_valid && (posY_rc > 0U) &&
                                                !is_bit_set(seqY.invalid_bm, m_posY + (uint32_t)wl) &&
                                                ((base_2bit_dev(seqY.bits, m_posY + (uint32_t)wl) ^ 3U) == prev_x_2bit) &&
                                                !is_bit_set(dropped_bitmap, m_posY + 1U);
                                bool emit = !left_ext;
                                if (!emit) {
                                    bool right_ext = next_x_valid && (m_posY > 0U) &&
                                                     !is_bit_set(seqY.invalid_bm, m_posY - 1U) &&
                                                     ((base_2bit_dev(seqY.bits, m_posY - 1U) ^ 3U) == next_x_2bit) &&
                                                     !is_bit_set(dropped_bitmap, m_posY - 1U);
                                    emit = !right_ext;
                                }
                                if (emit) {
                                    stage_hit_dev(s_hits, &s_cnt, out_hits, hit_count, max_hits,
                                                  pack_hit_key_tile(key_layout, 1U, posX, posY_rc));
                                }
                            }
                            bp_word |= (1U << s);
                            if (stats_buf) bypassed_local++;
                        }
                        /* Subcase 4c (palindromic or non-unique): emit 0, leave bp_word bit 0 for L2_PART */
                    }
                }
                /* Case 5 (non-anchor with !act_valid): leave bp_word bit 0 for L2_PART */
            }
            d_q_bypass_bm[w_out] = bp_word;
        }

        /* Uniform lockstep block flush across all 4 waves */
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
        if (wave < 3U) {
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

template <bool USE_BYPASS_BM>
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
    HitKeyLayout key_layout,
    const uint32_t * __restrict__ d_q_bypass_bm,
    uint32_t bp_words,
    uint32_t * __restrict__ d_part_counts
) {
    __shared__ uint32_t s_hist[1024];
    for (uint32_t p = threadIdx.x; p < P; p += blockDim.x) {
        s_hist[p] = 0U;
    }
    __syncthreads();

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);

    if constexpr (USE_BYPASS_BM) {
        int rev_rshift = 64 - 2 * wl;
        uint32_t gtid = blockIdx.x * blockDim.x + threadIdx.x;
        uint32_t w0 = gtid * 2U;
        for (uint32_t u = 0U; u < 2U; u++) {
            uint32_t w = w0 + u;
            if (w >= bp_words) break;
            uint32_t todo = ~__ldg(&d_q_bypass_bm[w]);
            while (todo != 0U) {
                uint32_t s = (uint32_t)(__ffs((int)todo) - 1);
                todo &= (todo - 1U);
                uint32_t rel_x = (w << 5) + s;
                uint32_t posX = key_layout.min_x + rel_x;
                uint64_t w_lsb = load_32_2bit(seqX.bits, posX);
                uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift;
                uint64_t kr = (~w_lsb) & mask;
                uint64_t canon = (kf <= kr) ? kf : kr;
                uint64_t h = fmix64(canon);
                uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                uint32_t part = bkt >> part_shift;
                atomicAdd(&s_hist[part], 1U);
            }
        }
    } else {
        const uint64_t STRIDE = 64ULL;
        uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
        uint64_t t_start = tile_start + gtid * STRIDE;
        if (t_start < tile_end) {
            uint64_t t_end = t_start + STRIDE;
            if (t_end > tile_end) t_end = tile_end;

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

/* W3-F13b: Two-phase CTA-privatized shared-memory scatter (4,096 B smem -> 100% T4 occ).
 * Phase 1 histograms CTA queries via smem ATOMS; Phase 1.5 reserves 1 contiguous global slice
 * per non-empty partition per CTA (~32x fewer global ATOMG); Phase 2 scatters into contiguous
 * cache-line bursts via smem ATOMS. NO early return before __syncthreads()! */
template <bool USE_BYPASS_BM>
static __global__ void l2_scatter_query_parts_kernel(
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
    HitKeyLayout key_layout,
    const uint32_t * __restrict__ d_q_bypass_bm,
    uint32_t bp_words,
    uint32_t * __restrict__ d_part_cursor,
    uint64_t * __restrict__ d_q_hash,
    uint32_t * __restrict__ d_q_meta
) {
    __shared__ uint32_t s_cnt[1024];
    for (uint32_t p = threadIdx.x; p < P; p += blockDim.x) {
        s_cnt[p] = 0U;
    }
    __syncthreads();

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);

    if constexpr (USE_BYPASS_BM) {
        int rev_rshift = 64 - 2 * wl;
        uint32_t gtid = blockIdx.x * blockDim.x + threadIdx.x;
        uint32_t w0 = gtid * 2U;

        /* Phase 1: CTA-local histogram in smem (ATOMS) */
        for (uint32_t u = 0U; u < 2U; u++) {
            uint32_t w = w0 + u;
            if (w >= bp_words) break;
            uint32_t todo = ~__ldg(&d_q_bypass_bm[w]);
            while (todo != 0U) {
                uint32_t s = (uint32_t)(__ffs((int)todo) - 1);
                todo &= (todo - 1U);
                uint32_t rel_x = (w << 5) + s;
                uint32_t posX = key_layout.min_x + rel_x;
                uint64_t w_lsb = load_32_2bit(seqX.bits, posX);
                uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift;
                uint64_t kr = (~w_lsb) & mask;
                uint64_t canon = (kf <= kr) ? kf : kr;
                uint64_t h = fmix64(canon);
                uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                uint32_t part = bkt >> part_shift;
                atomicAdd(&s_cnt[part], 1U);
            }
        }
        __syncthreads();

        /* Phase 1.5: Reserve 1 contiguous global slice per non-empty partition per CTA */
        for (uint32_t p = threadIdx.x; p < P; p += blockDim.x) {
            uint32_t c = s_cnt[p];
            if (c > 0U) {
                s_cnt[p] = atomicAdd(&d_part_cursor[p], c);
            }
        }
        __syncthreads();

        /* Phase 2: CTA-local scatter using smem ATOMS into contiguous global slices */
        for (uint32_t u = 0U; u < 2U; u++) {
            uint32_t w = w0 + u;
            if (w >= bp_words) break;
            uint32_t todo = ~__ldg(&d_q_bypass_bm[w]);
            while (todo != 0U) {
                uint32_t s = (uint32_t)(__ffs((int)todo) - 1);
                todo &= (todo - 1U);
                uint32_t rel_x = (w << 5) + s;
                uint32_t posX = key_layout.min_x + rel_x;
                uint64_t w_lsb = load_32_2bit(seqX.bits, posX);
                uint64_t kf = reverse_32_2bit_dev(w_lsb) >> rev_rshift;
                uint64_t kr = (~w_lsb) & mask;
                uint32_t x_strand = (kf <= kr) ? 0U : 1U;
                uint64_t canon = x_strand ? kr : kf;
                bool is_pal = (kf == kr);
                uint64_t h = fmix64(canon);
                uint32_t bkt = (uint32_t)(h >> shift_bits) - bucket_base;
                uint32_t part = bkt >> part_shift;
                uint32_t slot = atomicAdd(&s_cnt[part], 1U);
                d_q_hash[slot] = h;
                d_q_meta[slot] = (rel_x << 2) | ((uint32_t)is_pal << 1) | x_strand;
            }
        }
    } else {
        const uint64_t STRIDE = 64ULL;
        uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
        uint64_t t_start = tile_start + gtid * STRIDE;
        bool active = (t_start < tile_end);
        uint64_t t_end = active
            ? ((t_start + STRIDE < tile_end) ? (t_start + STRIDE) : tile_end)
            : t_start;

        uint64_t eff_mask = (seed_mask != 0ULL) ? seed_mask : mask;
        bool part_filter = (part_start != 0U || part_end != 256U);
        int rc_shift = 2 * (wl - 1);
        uint64_t scan_from = (t_start >= (uint64_t)(wl - 1)) ? (t_start - (uint64_t)(wl - 1)) : 0ULL;

        /* Phase 1: CTA-local histogram in smem (ATOMS) */
        if (active) {
            uint64_t i = scan_from;
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
                            atomicAdd(&s_cnt[part], 1U);
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

        /* Phase 1.5: Reserve 1 contiguous global slice per non-empty partition per CTA */
        for (uint32_t p = threadIdx.x; p < P; p += blockDim.x) {
            uint32_t c = s_cnt[p];
            if (c > 0U) {
                s_cnt[p] = atomicAdd(&d_part_cursor[p], c);
            }
        }
        __syncthreads();

        /* Phase 2: CTA-local scatter using smem ATOMS into contiguous global slices */
        if (active) {
            uint64_t i = scan_from;
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
                            uint32_t slot = atomicAdd(&s_cnt[part], 1U);
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

/* W3-F15b: 0-stack BLAST/LASTZ sorted-diagonal two-hit filter kernel.
 * Checks backward on the same (strand, diag) for any prior non-overlapping hit within
 * [wl, wl + two_hit_window], plus a 4-instruction backward self-2hit check for contiguous
 * L/R-maximal MEMs (L >= 2*wl) whose interior seeds were suppressed. */
static __global__ void two_hit_filter_kernel(
    const uint64_t * __restrict__ sorted_hits,
    uint32_t n_hits,
    HitKeyLayout lay,
    int wl,
    uint32_t two_hit_window,
    bool check_self_prev,
    const uint64_t * __restrict__ seqX_bits,
    const uint32_t * __restrict__ seqX_invalid_bm,
    const uint64_t * __restrict__ seqY_bits,
    const uint32_t * __restrict__ seqY_invalid_bm,
    const uint64_t * __restrict__ seqY_rc_bits,
    const uint32_t * __restrict__ dropped_bitmap,
    uint64_t lenY,
    uint8_t * __restrict__ out_flags
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_hits) return;

    uint64_t wl_u64 = (uint64_t)wl;
    uint64_t max_dist = wl_u64 + (uint64_t)two_hit_window;
    uint8_t keep = gecko_two_hit_scan_backward(sorted_hits, i, lay, wl_u64, max_dist);
    if (!keep && check_self_prev) {
        keep = gecko_two_hit_self_prev(
            sorted_hits[i], lay, wl,
            seqX_bits, seqX_invalid_bm,
            seqY_bits, seqY_invalid_bm, seqY_rc_bits,
            dropped_bitmap, lenY);
    }
    out_flags[i] = keep;
}

uint32_t gpu_query_tile_left_maximal(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuSeqView &seqY_rc,
    const GpuRefIndex &idx,
    int wl,
    uint64_t seed_mask,
    bool compare_reverse,
    uint32_t two_hit_window,
    uint64_t tile_start,
    uint64_t tile_end,
    HitKeyLayout key_layout,
    int l2_cache_bytes,
    int query_mode,
    uint32_t *d_stats_buf,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    GpuFrag *d_out_frags,
    uint32_t *d_hit_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed,
    uint64_t *&out_sorted_hits,
    uint64_t *&out_alt_buf,
    uint64_t *out_valid_q,
    uint64_t *out_bypassed_q,
    uint32_t *out_two_hit_dropped,
    GeckoCudaProfile *prof
) {
    out_overflowed = false;
    out_raw_leftmax = 0U;
    out_sorted_hits = d_hits;
    out_alt_buf = d_hits_sorted;
    if (out_valid_q) *out_valid_q = 0ULL;
    if (out_bypassed_q) *out_bypassed_q = 0ULL;
    if (out_two_hit_dropped) *out_two_hit_dropped = 0U;
    if (idx.n_kept == 0U || tile_start >= tile_end) return 0U;

    uint64_t first_cur_i = (tile_start >= (uint64_t)(wl - 1)) ? tile_start : (uint64_t)(wl - 1);
    uint32_t num_pos = (tile_end > first_cur_i) ? (uint32_t)(tile_end - first_cur_i) : 0U;
    if (num_pos == 0U) return 0U;

    /* Zero d_hit_count ONCE per tile (both l2_pre_bypass and l2_probe append to d_hits) */
    CUDA_CHECK_VOID(cudaMemset(d_hit_count, 0, sizeof(uint32_t)));

    uint64_t tile_len = tile_end - tile_start;
    bool ran_l2_part = false;

    /* F16 + W3-F13b: L2-Resident Wave-Scheduled Hash-Partitioned Probe (optional Hybrid pre-bypass) */
    if (query_mode == GECKO_QMODE_L2_PART || query_mode == GECKO_QMODE_HYBRID) {
        L2PartPlan plan = choose_l2_part_plan(idx, l2_cache_bytes);
        if (plan.P >= 2U &&
            tile_len <= (1ULL << 24) &&
            tile_len >= (uint64_t)plan.P * 4096ULL) {
            bool use_pre_bypass = (query_mode == GECKO_QMODE_HYBRID) && (idx.d_ref_non_unique_bm != nullptr);
            uint32_t bp_words = (num_pos + 31U) >> 5;
            size_t bypass_bytes = use_pre_bypass ? ((size_t)bp_words * sizeof(uint32_t)) : 0ULL;

            auto align256 = [](size_t n) -> size_t { return (n + 255ULL) & ~255ULL; };
            size_t max_q = (size_t)num_pos;
            size_t off_hash   = 0ULL;
            size_t off_meta   = align256(off_hash   + max_q * sizeof(uint64_t));
            size_t off_bypass = align256(off_meta   + max_q * sizeof(uint32_t));
            size_t off_counts = align256(off_bypass + bypass_bytes);
            size_t off_offs   = align256(off_counts + (size_t)plan.P * sizeof(uint32_t));
            size_t off_cursor = align256(off_offs   + ((size_t)plan.P + 1ULL) * sizeof(uint32_t));
            size_t off_blkoff = align256(off_cursor + (size_t)plan.P * sizeof(uint32_t));
            size_t off_totals = align256(off_blkoff + ((size_t)plan.P + 1ULL) * sizeof(uint32_t));
            size_t total_need = align256(off_totals + 2ULL * sizeof(uint32_t));
            size_t scratch_bytes = (size_t)max_hits * sizeof(uint64_t);

            if (total_need <= scratch_bytes) {
                if (prof && prof->enabled) prof->l2_part_P = plan.P;

                uint8_t *base_ptr = reinterpret_cast<uint8_t *>(d_hits_sorted);
                uint64_t *d_q_hash         = reinterpret_cast<uint64_t *>(base_ptr + off_hash);
                uint32_t *d_q_meta         = reinterpret_cast<uint32_t *>(base_ptr + off_meta);
                uint32_t *d_q_bypass_bm    = use_pre_bypass
                    ? reinterpret_cast<uint32_t *>(base_ptr + off_bypass)
                    : nullptr;
                uint32_t *d_part_counts    = reinterpret_cast<uint32_t *>(base_ptr + off_counts);
                uint32_t *d_part_offsets   = reinterpret_cast<uint32_t *>(base_ptr + off_offs);
                uint32_t *d_part_cursor    = reinterpret_cast<uint32_t *>(base_ptr + off_cursor);
                uint32_t *d_part_block_off = reinterpret_cast<uint32_t *>(base_ptr + off_blkoff);
                uint32_t *d_totals         = reinterpret_cast<uint32_t *>(base_ptr + off_totals);

                if (use_pre_bypass) {
                    uint32_t *use_stats = (out_valid_q && out_bypassed_q && d_stats_buf) ? d_stats_buf : nullptr;
                    if (use_stats) {
                        CUDA_CHECK_VOID(cudaMemset(use_stats, 0, 2 * sizeof(uint32_t)));
                    }
                    uint32_t n_threads_bp = (bp_words + 3U) >> 2;
                    uint32_t grid_bp = (n_threads_bp + 255U) >> 8;
                    {
                        CudaEventScope ev(prof, prof ? &prof->l2_pre_bypass_ms : nullptr);
                        l2_pre_bypass_kernel<<<grid_bp, 256>>>(
                            seqX, seqY,
                            idx.d_bucket_off, idx.d_tags, idx.d_vals,
                            idx.d_dropped_bitmap, idx.d_ref_non_unique_bm,
                            idx.shift_bits, idx.tag_shift, idx.bucket_base,
                            wl, compare_reverse, num_pos, bp_words, key_layout,
                            d_q_bypass_bm, d_hits, d_hit_count, max_hits, use_stats);
                        CUDA_CHECK_VOID(cudaGetLastError());
                    }
                    if (use_stats) {
                        uint32_t h_stats[2] = {0U, 0U};
                        CUDA_CHECK_VOID(cudaMemcpy(h_stats, use_stats, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        *out_valid_q = (uint64_t)h_stats[0];
                        *out_bypassed_q = (uint64_t)h_stats[1];
                    }
                }

                CUDA_CHECK_VOID(cudaMemset(d_part_counts, 0, (size_t)plan.P * sizeof(uint32_t)));
                uint32_t grid_p1 = use_pre_bypass
                    ? (((bp_words + 1U) >> 1) + 255U) >> 8
                    : (uint32_t)((((tile_len + 63ULL) / 64ULL) + 255ULL) / 256ULL);

                {
                    CudaEventScope ev(prof, prof ? &prof->l2_count_ms : nullptr);
                    if (use_pre_bypass) {
                        l2_count_query_parts_kernel<true><<<grid_p1, 256>>>(
                            seqX, idx.shift_bits, idx.part_start, idx.part_end, idx.bucket_base,
                            plan.P, plan.part_shift, wl, seed_mask, tile_start, tile_end,
                            key_layout, d_q_bypass_bm, bp_words, d_part_counts);
                    } else {
                        l2_count_query_parts_kernel<false><<<grid_p1, 256>>>(
                            seqX, idx.shift_bits, idx.part_start, idx.part_end, idx.bucket_base,
                            plan.P, plan.part_shift, wl, seed_mask, tile_start, tile_end,
                            key_layout, nullptr, 0U, d_part_counts);
                    }
                    CUDA_CHECK_VOID(cudaGetLastError());

                    l2_prefix_sum_parts_kernel<<<1, 256>>>(
                        d_part_counts, plan.P, d_part_offsets, d_part_cursor, d_part_block_off, d_totals);
                    CUDA_CHECK_VOID(cudaGetLastError());
                }

                uint32_t h_totals[2] = {0U, 0U};
                CUDA_CHECK_VOID(cudaMemcpy(h_totals, d_totals, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));

                /* CRITICAL: When h_totals[0] == 0U (e.g. 100% bypassed by l2_pre_bypass_kernel),
                 * skip scatter & probe and FALL THROUGH to read d_hit_count + SortKeys! */
                if (h_totals[0] > 0U && h_totals[1] > 0U) {
                    {
                        CudaEventScope ev(prof, prof ? &prof->l2_scatter_ms : nullptr);
                        if (use_pre_bypass) {
                            l2_scatter_query_parts_kernel<true><<<grid_p1, 256>>>(
                                seqX, idx.shift_bits, idx.part_start, idx.part_end, idx.bucket_base,
                                plan.P, plan.part_shift, wl, seed_mask, tile_start, tile_end,
                                key_layout, d_q_bypass_bm, bp_words, d_part_cursor, d_q_hash, d_q_meta);
                        } else {
                            l2_scatter_query_parts_kernel<false><<<grid_p1, 256>>>(
                                seqX, idx.shift_bits, idx.part_start, idx.part_end, idx.bucket_base,
                                plan.P, plan.part_shift, wl, seed_mask, tile_start, tile_end,
                                key_layout, nullptr, 0U, d_part_cursor, d_q_hash, d_q_meta);
                        }
                        CUDA_CHECK_VOID(cudaGetLastError());
                    }

                    {
                        CudaEventScope ev(prof, prof ? &prof->l2_probe_ms : nullptr);
                        l2_probe_partitioned_kernel<<<h_totals[1], 256>>>(
                            seqX, seqY,
                            idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap,
                            idx.shift_bits, idx.tag_shift, idx.part_start, idx.part_end, idx.bucket_base,
                            plan.P, wl, seed_mask, compare_reverse, idx.has_global_dropped_bm,
                            key_layout, d_part_offsets, d_part_block_off, d_q_hash, d_q_meta,
                            d_hits, d_hit_count, max_hits);
                        CUDA_CHECK_VOID(cudaGetLastError());
                    }
                }
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

        {
            CudaEventScope ev(prof, prof ? &prof->mono_probe_ms : nullptr);
            if (idx.is_direct_csr) {
                query_left_maximal_kernel<true, false><<<grid, 256>>>(
                    seqX, seqY,
                    idx.d_bucket_off, nullptr, idx.d_vals, idx.d_dropped_bitmap, nullptr,
                    0, 0, 0U, 256U, 0U,
                    wl, seed_mask, compare_reverse, idx.has_global_dropped_bm,
                    tile_start, tile_end, num_waves, key_layout,
                    d_hits, d_hit_count, max_hits, use_stats);
            } else if (idx.d_ref_non_unique_bm != nullptr) {
                query_left_maximal_kernel<false, true><<<grid, 256>>>(
                    seqX, seqY,
                    idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap, idx.d_ref_non_unique_bm,
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
        }

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
    if (h_nhits == 0U) {
        out_sorted_hits = d_hits;
        out_alt_buf = d_hits_sorted;
        return 0U;
    }
    if (h_nhits == 1U && two_hit_window == 0U) {
        out_sorted_hits = d_hits;
        out_alt_buf = d_hits_sorted;
        return 1U;
    }

    uint64_t *sorted_buf = d_hits;
    uint64_t *alt_buf = d_hits_sorted;

    /* F14 + W3-F11b: Tile-relative bit-bounded CUB DoubleBuffer SortKeys (CCCL 3.0 safe 6-arg).
     * Hands Current() to extend and Alternate() as d_long_seeds with zero extra VRAM and zero
     * min_gap heuristic loss. */
    if (h_nhits > 1U) {
        int end_bit = (int)key_layout.end_bit;
        cub::DoubleBuffer<uint64_t> d_keys_buf(d_hits, d_hits_sorted);
        size_t sort_bytes = 0;
        cub::DeviceRadixSort::SortKeys(nullptr, sort_bytes, d_keys_buf, (int)h_nhits, 0, end_bit);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sort_bytes);
        {
            CudaEventScope ev(prof, prof ? &prof->sort_hits_ms : nullptr);
            cub::DeviceRadixSort::SortKeys(d_cub_temp, sort_bytes, d_keys_buf, (int)h_nhits, 0, end_bit);
        }
        sorted_buf = d_keys_buf.Current();
        alt_buf = d_keys_buf.Alternate();
    }

    /* W3-F15b: Optional sorted-diagonal --two-hit [WINDOW] filter using idle d_out_frags as
     * 1-byte/hit flag buffer (0 B extra permanent VRAM; frags_cap * 24 >= 3 * max_hits bytes). */
    if (two_hit_window > 0U) {
        uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
        bool is_spaced = (seed_mask != 0ULL && seed_mask != span_mask);
        bool part_filter = (idx.part_start != 0U || idx.part_end != 256U);
        bool check_self_prev = !is_spaced && (!part_filter || idx.has_global_dropped_bm);
        uint8_t *d_flags = reinterpret_cast<uint8_t *>(d_out_frags);

        size_t sel_bytes = 0;
        cub::DeviceSelect::Flagged(
            nullptr, sel_bytes, sorted_buf, d_flags, alt_buf, d_hit_count, (int)h_nhits);
        gpu_ensure_cub_temp(d_cub_temp, cub_temp_cap, sel_bytes);

        {
            CudaEventScope ev(prof, prof ? &prof->two_hit_ms : nullptr);
            two_hit_filter_kernel<<<(h_nhits + 255U) / 256U, 256>>>(
                sorted_buf, h_nhits, key_layout, wl, two_hit_window, check_self_prev,
                seqX.bits, seqX.invalid_bm,
                seqY.bits, seqY.invalid_bm, seqY_rc.bits,
                idx.d_dropped_bitmap, seqY.len, d_flags);
            CUDA_CHECK_VOID(cudaGetLastError());

            cub::DeviceSelect::Flagged(
                d_cub_temp, sel_bytes, sorted_buf, d_flags, alt_buf, d_hit_count, (int)h_nhits);
        }

        uint32_t h_kept = 0U;
        CUDA_CHECK_VOID(cudaMemcpy(&h_kept, d_hit_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        if (out_two_hit_dropped) {
            *out_two_hit_dropped = h_nhits - h_kept;
        }
        out_sorted_hits = alt_buf;
        out_alt_buf = sorted_buf;
        return h_kept;
    }

    out_sorted_hits = sorted_buf;
    out_alt_buf = alt_buf;
    return h_nhits;
}
