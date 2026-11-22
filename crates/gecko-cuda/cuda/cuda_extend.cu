/*
 * cuda_extend.cu
 *
 * Stage 3: Two-tier 2-bit GPU bidirectional ungapped extension matching GECKO's
 * exact scoring model (+match, -mismatch, -n_penalty, SEQ_SEP, rare IUPAC):
 *  - Tier 1 (Short, 1 thread/seed, 512 bp cap): 32-base word-loaded register stepping.
 *  - Tier 2 (Long, 1 full 32-thread warp/seed): single warp-broadcast 64-bit word load
 *    per 32 bp + 1-cycle 32-bp exact-match fast path + zero-divergence warp prefix sums.
 */

#include <climits>
#include "cuda_common.cuh"

static __device__ __forceinline__ int score_special_pair(
    const GpuSeqView &seqX,
    uint32_t px,
    bool sp_x,
    const GpuSeqView &dy,
    uint32_t py,
    bool sp_y,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    bool &is_sep,
    int &is_match
) {
    uint8_t cx = sp_x ? decode_special_byte(seqX, px) : base_2bit_to_ascii(base_2bit_dev(seqX.bits, px));
    uint8_t cy = sp_y ? decode_special_byte(dy, py) : base_2bit_to_ascii(base_2bit_dev(dy.bits, py));
    if (cx == 0U || cy == 0U) {
        is_sep = true;
        is_match = 0;
        return 0;
    }
    if (cx == (uint8_t)'N' || cy == (uint8_t)'N') {
        is_match = 0;
        return -n_penalty;
    }
    if (cx == cy) {
        is_match = 1;
        return match_score;
    }
    is_match = 0;
    return -mismatch_penalty;
}

static __global__ void extend_short_kernel(
    const uint64_t * __restrict__ hits,
    uint32_t n_hits,
    GpuSeqView seqX,
    GpuSeqView seqY,
    GpuSeqView seqY_rc,
    int wl,
    uint32_t min_len,
    float sim_th,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    uint64_t * __restrict__ long_seeds,
    uint32_t * __restrict__ long_count,
    GpuFrag * __restrict__ out_frags,
    uint32_t * __restrict__ frag_count,
    uint32_t max_frags
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_hits) return;

    uint64_t key = hits[idx];
    uint32_t is_rev, posX, posY;
    int64_t diag;
    unpack_hit_key(key, is_rev, posX, posY, diag);

    const GpuSeqView &dy = is_rev ? seqY_rc : seqY;

    int64_t xfil = (int64_t)posX + wl;
    int64_t ycol = (int64_t)posY + wl;
    int64_t rem_fwd = ((int64_t)seqX.len - xfil < (int64_t)dy.len - ycol)
                          ? ((int64_t)seqX.len - xfil)
                          : ((int64_t)dy.len - ycol);
    int64_t lim_fwd = (rem_fwd < GECKO_SHORT_CAP) ? rem_fwd : GECKO_SHORT_CAP;

    int64_t xfilmax = (int64_t)posX + wl - 1;
    int64_t ycolmax = (int64_t)posY + wl - 1;
    int fscore = match_score * wl;
    int fscoreMax = fscore;
    int nIdent = wl, maxIdent = wl;

    bool fwd_terminated = false;
    for (int64_t off = 0; off < lim_fwd && !fwd_terminated; off += 32) {
        int chunk = (int)((lim_fwd - off < 32) ? (lim_fwd - off) : 32);
        uint32_t px0 = (uint32_t)(xfil + off);
        uint32_t py0 = (uint32_t)(ycol + off);
        uint64_t wx = load_32_2bit(seqX.bits, px0);
        uint64_t wy = load_32_2bit(dy.bits, py0);
        uint32_t spx = load_32_1bit(seqX.special_bm, px0);
        uint32_t spy = load_32_1bit(dy.special_bm, py0);
        uint32_t sp_any = spx | spy;

        if (chunk == 32 && sp_any == 0U && wx == wy) {
            fscore += 32 * match_score;
            nIdent += 32;
            if (fscoreMax < fscore) {
                fscoreMax = fscore;
                xfilmax = xfil + off + 31;
                ycolmax = ycol + off + 31;
                maxIdent = nIdent;
            }
            continue;
        }

        if (sp_any == 0U) {
            for (int k = 0; k < chunk; k++) {
                if ((wx & 3ULL) == (wy & 3ULL)) {
                    fscore += match_score;
                    nIdent++;
                    if (fscoreMax < fscore) {
                        fscoreMax = fscore;
                        xfilmax = xfil + off + k;
                        ycolmax = ycol + off + k;
                        maxIdent = nIdent;
                    }
                } else {
                    fscore -= mismatch_penalty;
                    if (fscore < 0) {
                        fwd_terminated = true;
                        break;
                    }
                }
                wx >>= 2;
                wy >>= 2;
            }
        } else {
            for (int k = 0; k < chunk; k++) {
                if (((sp_any >> k) & 1U) == 0U) {
                    if ((wx & 3ULL) == (wy & 3ULL)) {
                        fscore += match_score;
                        nIdent++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax = xfil + off + k;
                            ycolmax = ycol + off + k;
                            maxIdent = nIdent;
                        }
                    } else {
                        fscore -= mismatch_penalty;
                    }
                } else {
                    bool is_sep = false;
                    int is_match = 0;
                    int delta = score_special_pair(
                        seqX, px0 + (uint32_t)k, ((spx >> k) & 1U) != 0U,
                        dy, py0 + (uint32_t)k, ((spy >> k) & 1U) != 0U,
                        match_score, mismatch_penalty, n_penalty, is_sep, is_match);
                    if (is_sep) {
                        fwd_terminated = true;
                        break;
                    }
                    fscore += delta;
                    if (is_match) {
                        nIdent++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax = xfil + off + k;
                            ycolmax = ycol + off + k;
                            maxIdent = nIdent;
                        }
                    }
                }
                if (fscore < 0) {
                    fwd_terminated = true;
                    break;
                }
                wx >>= 2;
                wy >>= 2;
            }
        }
    }
    if (!fwd_terminated && rem_fwd > GECKO_SHORT_CAP) {
        uint32_t slot = atomicAdd(long_count, 1U);
        long_seeds[slot] = key;
        return;
    }

    int64_t xfil2 = (int64_t)posX - 1;
    int64_t ycol2 = (int64_t)posY - 1;
    int64_t rem_bwd = (posX < posY) ? (int64_t)posX : (int64_t)posY;
    int64_t lim_bwd = (rem_bwd < GECKO_SHORT_CAP) ? rem_bwd : GECKO_SHORT_CAP;

    int64_t xfilmax2 = (int64_t)posX;
    int64_t ycolmax2 = (int64_t)posY;
    fscore = fscoreMax;
    nIdent = maxIdent;

    bool bwd_terminated = false;
    for (int64_t off = 0; off < lim_bwd && !bwd_terminated; off += 32) {
        int64_t rem_chunk = lim_bwd - off;
        if (rem_chunk >= 32) {
            uint32_t bx0 = (uint32_t)(xfil2 - off - 31);
            uint32_t by0 = (uint32_t)(ycol2 - off - 31);
            uint64_t wx = reverse_32_2bit_dev(load_32_2bit(seqX.bits, bx0));
            uint64_t wy = reverse_32_2bit_dev(load_32_2bit(dy.bits, by0));
            uint32_t spx = __brev(load_32_1bit(seqX.special_bm, bx0));
            uint32_t spy = __brev(load_32_1bit(dy.special_bm, by0));
            uint32_t sp_any = spx | spy;

            if (sp_any == 0U && wx == wy) {
                fscore += 32 * match_score;
                nIdent += 32;
                if (fscoreMax < fscore) {
                    fscoreMax = fscore;
                    xfilmax2 = xfil2 - off - 31;
                    ycolmax2 = ycol2 - off - 31;
                    maxIdent = nIdent;
                }
                continue;
            }

            if (sp_any == 0U) {
                for (int k = 0; k < 32; k++) {
                    if ((wx & 3ULL) == (wy & 3ULL)) {
                        fscore += match_score;
                        nIdent++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax2 = xfil2 - off - k;
                            ycolmax2 = ycol2 - off - k;
                            maxIdent = nIdent;
                        }
                    } else {
                        fscore -= mismatch_penalty;
                        if (fscore < 0) {
                            bwd_terminated = true;
                            break;
                        }
                    }
                    wx >>= 2;
                    wy >>= 2;
                }
            } else {
                for (int k = 0; k < 32; k++) {
                    if (((sp_any >> k) & 1U) == 0U) {
                        if ((wx & 3ULL) == (wy & 3ULL)) {
                            fscore += match_score;
                            nIdent++;
                            if (fscoreMax < fscore) {
                                fscoreMax = fscore;
                                xfilmax2 = xfil2 - off - k;
                                ycolmax2 = ycol2 - off - k;
                                maxIdent = nIdent;
                            }
                        } else {
                            fscore -= mismatch_penalty;
                        }
                    } else {
                        bool is_sep = false;
                        int is_match = 0;
                        int delta = score_special_pair(
                            seqX, (uint32_t)(xfil2 - off - k), ((spx >> k) & 1U) != 0U,
                            dy, (uint32_t)(ycol2 - off - k), ((spy >> k) & 1U) != 0U,
                            match_score, mismatch_penalty, n_penalty, is_sep, is_match);
                        if (is_sep) {
                            bwd_terminated = true;
                            break;
                        }
                        fscore += delta;
                        if (is_match) {
                            nIdent++;
                            if (fscoreMax < fscore) {
                                fscoreMax = fscore;
                                xfilmax2 = xfil2 - off - k;
                                ycolmax2 = ycol2 - off - k;
                                maxIdent = nIdent;
                            }
                        }
                    }
                    if (fscore < 0) {
                        bwd_terminated = true;
                        break;
                    }
                    wx >>= 2;
                    wy >>= 2;
                }
            }
        } else {
            int chunk = (int)rem_chunk;
            for (int k = 0; k < chunk; k++) {
                uint32_t px = (uint32_t)(xfil2 - off - k);
                uint32_t py = (uint32_t)(ycol2 - off - k);
                bool sp_x = is_bit_set(seqX.special_bm, px);
                bool sp_y = is_bit_set(dy.special_bm, py);
                if (!sp_x && !sp_y) {
                    if (base_2bit_dev(seqX.bits, px) == base_2bit_dev(dy.bits, py)) {
                        fscore += match_score;
                        nIdent++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax2 = px;
                            ycolmax2 = py;
                            maxIdent = nIdent;
                        }
                    } else {
                        fscore -= mismatch_penalty;
                    }
                } else {
                    bool is_sep = false;
                    int is_match = 0;
                    int delta = score_special_pair(
                        seqX, px, sp_x, dy, py, sp_y,
                        match_score, mismatch_penalty, n_penalty, is_sep, is_match);
                    if (is_sep) {
                        bwd_terminated = true;
                        break;
                    }
                    fscore += delta;
                    if (is_match) {
                        nIdent++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax2 = px;
                            ycolmax2 = py;
                            maxIdent = nIdent;
                        }
                    }
                }
                if (fscore < 0) {
                    bwd_terminated = true;
                    break;
                }
            }
        }
    }
    if (!bwd_terminated && rem_bwd > GECKO_SHORT_CAP) {
        uint32_t slot = atomicAdd(long_count, 1U);
        long_seeds[slot] = key;
        return;
    }

    uint32_t length = (uint32_t)(xfilmax - xfilmax2 + 1);
    float sim = (float)((double)fscoreMax * 100.0 / (double)((uint64_t)length * (uint64_t)match_score));
    if (length > min_len && sim > sim_th) {
        uint32_t slot = atomicAdd(frag_count, 1U);
        if (slot < max_frags) {
            out_frags[slot].diag = diag;
            out_frags[slot].xStart = (uint32_t)xfilmax2;
            out_frags[slot].yStart = (uint32_t)ycolmax2;
            out_frags[slot].xEnd = (uint32_t)xfilmax;
            out_frags[slot].yEnd = (uint32_t)ycolmax;
            out_frags[slot].ident = (uint32_t)maxIdent;
            out_frags[slot].score = (uint32_t)fscoreMax;
            out_frags[slot].strand = (uint8_t)is_rev;
        }
    }
}

static __device__ __forceinline__ int warp_inclusive_sum(int val, int lane) {
    #pragma unroll
    for (int off = 1; off < 32; off <<= 1) {
        int n = __shfl_up_sync(0xFFFFFFFFU, val, off);
        if (lane >= off) val += n;
    }
    return val;
}

static __device__ __forceinline__ int warp_max(int val) {
    #pragma unroll
    for (int off = 16; off > 0; off >>= 1) {
        int other = __shfl_xor_sync(0xFFFFFFFFU, val, off);
        if (other > val) val = other;
    }
    return val;
}

static __global__ void extend_long_warp_kernel(
    const uint64_t * __restrict__ long_seeds,
    uint32_t n_long,
    GpuSeqView seqX,
    GpuSeqView seqY,
    GpuSeqView seqY_rc,
    int wl,
    uint32_t min_len,
    float sim_th,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    GpuFrag * __restrict__ out_frags,
    uint32_t * __restrict__ frag_count,
    uint32_t max_frags
) {
    uint32_t warp_id = ((uint32_t)blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    int lane = (int)(threadIdx.x & 31U);
    if (warp_id >= n_long) return;

    uint64_t key = long_seeds[warp_id];
    uint32_t is_rev, posX, posY;
    int64_t diag;
    unpack_hit_key(key, is_rev, posX, posY, diag);

    const GpuSeqView &dy = is_rev ? seqY_rc : seqY;

    int64_t xfil = (int64_t)posX + wl;
    int64_t ycol = (int64_t)posY + wl;
    int64_t rem_fwd = ((int64_t)seqX.len - xfil < (int64_t)dy.len - ycol)
                          ? ((int64_t)seqX.len - xfil)
                          : ((int64_t)dy.len - ycol);

    int64_t xfilmax = (int64_t)posX + wl - 1;
    int64_t ycolmax = (int64_t)posY + wl - 1;
    int fscore = match_score * wl;
    int fscoreMax = fscore;
    int nIdent = wl, maxIdent = wl;

    for (int64_t off = 0; off < rem_fwd; off += 32) {
        int64_t rem_chunk = rem_fwd - off;
        uint64_t wx = 0ULL, wy = 0ULL;
        uint32_t spx = 0U, spy = 0U;
        if (lane == 0) {
            uint32_t px0 = (uint32_t)(xfil + off);
            uint32_t py0 = (uint32_t)(ycol + off);
            wx = load_32_2bit(seqX.bits, px0);
            wy = load_32_2bit(dy.bits, py0);
            spx = load_32_1bit(seqX.special_bm, px0);
            spy = load_32_1bit(dy.special_bm, py0);
        }
        wx = __shfl_sync(0xFFFFFFFFU, wx, 0);
        wy = __shfl_sync(0xFFFFFFFFU, wy, 0);
        spx = __shfl_sync(0xFFFFFFFFU, spx, 0);
        spy = __shfl_sync(0xFFFFFFFFU, spy, 0);
        uint32_t sp_any = spx | spy;

        /* 1-cycle fast path when all 32 bases are clean and identical */
        if (rem_chunk >= 32 && sp_any == 0U && wx == wy) {
            fscore += 32 * match_score;
            nIdent += 32;
            if (fscore > fscoreMax) {
                fscoreMax = fscore;
                xfilmax = xfil + off + 31;
                ycolmax = ycol + off + 31;
                maxIdent = nIdent;
            }
            continue;
        }

        bool in_bounds = (lane < rem_chunk);
        int delta = 0, is_match = 0;
        bool is_sep = !in_bounds;

        if (in_bounds) {
            if (((sp_any >> lane) & 1U) == 0U) {
                uint32_t bx = (uint32_t)((wx >> (lane << 1)) & 3ULL);
                uint32_t by = (uint32_t)((wy >> (lane << 1)) & 3ULL);
                if (bx == by) { delta = match_score; is_match = 1; }
                else { delta = -mismatch_penalty; }
            } else {
                delta = score_special_pair(
                    seqX, (uint32_t)(xfil + off + lane), ((spx >> lane) & 1U) != 0U,
                    dy, (uint32_t)(ycol + off + lane), ((spy >> lane) & 1U) != 0U,
                    match_score, mismatch_penalty, n_penalty, is_sep, is_match);
            }
        }

        uint32_t sep_mask = __ballot_sync(0xFFFFFFFFU, is_sep);
        int sep_lane = sep_mask ? (__ffs(sep_mask) - 1) : 32;
        if (lane >= sep_lane) { delta = 0; is_match = 0; }

        int psum = warp_inclusive_sum(delta, lane);
        int pmatch = warp_inclusive_sum(is_match, lane);
        int score_here = fscore + psum;

        uint32_t neg_mask = __ballot_sync(0xFFFFFFFFU, (lane < sep_lane) && (score_here < 0));
        int neg_lane = neg_mask ? (__ffs(neg_mask) - 1) : 32;
        int active_end = (sep_lane < neg_lane) ? sep_lane : neg_lane;

        int cand_score = (lane < active_end) ? score_here : INT_MIN;
        int chunk_max = warp_max(cand_score);
        if (chunk_max > fscoreMax) {
            uint32_t max_mask = __ballot_sync(0xFFFFFFFFU, (lane < active_end) && (score_here == chunk_max));
            int best_lane = __ffs(max_mask) - 1;
            int best_pmatch = __shfl_sync(0xFFFFFFFFU, pmatch, best_lane);
            fscoreMax = chunk_max;
            xfilmax = xfil + off + best_lane;
            ycolmax = ycol + off + best_lane;
            maxIdent = nIdent + best_pmatch;
        }

        if (sep_mask || neg_mask) break;

        fscore = __shfl_sync(0xFFFFFFFFU, score_here, 31);
        nIdent += __shfl_sync(0xFFFFFFFFU, pmatch, 31);
    }

    int64_t xfil2 = (int64_t)posX - 1;
    int64_t ycol2 = (int64_t)posY - 1;
    int64_t rem_bwd = (posX < posY) ? (int64_t)posX : (int64_t)posY;
    int64_t xfilmax2 = (int64_t)posX;
    int64_t ycolmax2 = (int64_t)posY;
    fscore = fscoreMax;
    nIdent = maxIdent;

    for (int64_t off = 0; off < rem_bwd; off += 32) {
        int64_t rem_chunk = rem_bwd - off;
        int delta = 0, is_match = 0;
        bool is_sep = false;

        if (rem_chunk >= 32) {
            uint64_t wx = 0ULL, wy = 0ULL;
            uint32_t spx = 0U, spy = 0U;
            if (lane == 0) {
                uint32_t bx0 = (uint32_t)(xfil2 - off - 31);
                uint32_t by0 = (uint32_t)(ycol2 - off - 31);
                wx = reverse_32_2bit_dev(load_32_2bit(seqX.bits, bx0));
                wy = reverse_32_2bit_dev(load_32_2bit(dy.bits, by0));
                spx = __brev(load_32_1bit(seqX.special_bm, bx0));
                spy = __brev(load_32_1bit(dy.special_bm, by0));
            }
            wx = __shfl_sync(0xFFFFFFFFU, wx, 0);
            wy = __shfl_sync(0xFFFFFFFFU, wy, 0);
            spx = __shfl_sync(0xFFFFFFFFU, spx, 0);
            spy = __shfl_sync(0xFFFFFFFFU, spy, 0);
            uint32_t sp_any = spx | spy;

            if (sp_any == 0U && wx == wy) {
                fscore += 32 * match_score;
                nIdent += 32;
                if (fscore > fscoreMax) {
                    fscoreMax = fscore;
                    xfilmax2 = xfil2 - off - 31;
                    ycolmax2 = ycol2 - off - 31;
                    maxIdent = nIdent;
                }
                continue;
            }

            if (((sp_any >> lane) & 1U) == 0U) {
                uint32_t bx = (uint32_t)((wx >> (lane << 1)) & 3ULL);
                uint32_t by = (uint32_t)((wy >> (lane << 1)) & 3ULL);
                if (bx == by) { delta = match_score; is_match = 1; }
                else { delta = -mismatch_penalty; }
            } else {
                delta = score_special_pair(
                    seqX, (uint32_t)(xfil2 - off - lane), ((spx >> lane) & 1U) != 0U,
                    dy, (uint32_t)(ycol2 - off - lane), ((spy >> lane) & 1U) != 0U,
                    match_score, mismatch_penalty, n_penalty, is_sep, is_match);
            }
        } else {
            bool in_bounds = (lane < rem_chunk);
            is_sep = !in_bounds;
            if (in_bounds) {
                uint32_t px = (uint32_t)(xfil2 - off - lane);
                uint32_t py = (uint32_t)(ycol2 - off - lane);
                bool sp_x = is_bit_set(seqX.special_bm, px);
                bool sp_y = is_bit_set(dy.special_bm, py);
                if (!sp_x && !sp_y) {
                    if (base_2bit_dev(seqX.bits, px) == base_2bit_dev(dy.bits, py)) {
                        delta = match_score;
                        is_match = 1;
                    } else {
                        delta = -mismatch_penalty;
                    }
                } else {
                    delta = score_special_pair(
                        seqX, px, sp_x, dy, py, sp_y,
                        match_score, mismatch_penalty, n_penalty, is_sep, is_match);
                }
            }
        }

        uint32_t sep_mask = __ballot_sync(0xFFFFFFFFU, is_sep);
        int sep_lane = sep_mask ? (__ffs(sep_mask) - 1) : 32;
        if (lane >= sep_lane) { delta = 0; is_match = 0; }

        int psum = warp_inclusive_sum(delta, lane);
        int pmatch = warp_inclusive_sum(is_match, lane);
        int score_here = fscore + psum;

        uint32_t neg_mask = __ballot_sync(0xFFFFFFFFU, (lane < sep_lane) && (score_here < 0));
        int neg_lane = neg_mask ? (__ffs(neg_mask) - 1) : 32;
        int active_end = (sep_lane < neg_lane) ? sep_lane : neg_lane;

        int cand_score = (lane < active_end) ? score_here : INT_MIN;
        int chunk_max = warp_max(cand_score);
        if (chunk_max > fscoreMax) {
            uint32_t max_mask = __ballot_sync(0xFFFFFFFFU, (lane < active_end) && (score_here == chunk_max));
            int best_lane = __ffs(max_mask) - 1;
            int best_pmatch = __shfl_sync(0xFFFFFFFFU, pmatch, best_lane);
            fscoreMax = chunk_max;
            xfilmax2 = xfil2 - off - best_lane;
            ycolmax2 = ycol2 - off - best_lane;
            maxIdent = nIdent + best_pmatch;
        }

        if (sep_mask || neg_mask) break;

        fscore = __shfl_sync(0xFFFFFFFFU, score_here, 31);
        nIdent += __shfl_sync(0xFFFFFFFFU, pmatch, 31);
    }

    if (lane == 0) {
        uint32_t length = (uint32_t)(xfilmax - xfilmax2 + 1);
        float sim = (float)((double)fscoreMax * 100.0 / (double)((uint64_t)length * (uint64_t)match_score));
        if (length > min_len && sim > sim_th) {
            uint32_t slot = atomicAdd(frag_count, 1U);
            if (slot < max_frags) {
                out_frags[slot].diag = diag;
                out_frags[slot].xStart = (uint32_t)xfilmax2;
                out_frags[slot].yStart = (uint32_t)ycolmax2;
                out_frags[slot].xEnd = (uint32_t)xfilmax;
                out_frags[slot].yEnd = (uint32_t)ycolmax;
                out_frags[slot].ident = (uint32_t)maxIdent;
                out_frags[slot].score = (uint32_t)fscoreMax;
                out_frags[slot].strand = (uint8_t)is_rev;
            }
        }
    }
}

uint32_t gpu_extend_tile_two_tier(
    const uint64_t *d_filtered_hits,
    uint32_t n_hits,
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuSeqView &seqY_rc,
    int wl,
    uint32_t min_len,
    float sim_th,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    uint64_t *d_long_seeds,
    uint32_t *d_long_count,
    GpuFrag *d_out_frags,
    uint32_t *d_frag_count,
    uint32_t max_frags,
    bool &out_overflowed
) {
    out_overflowed = false;
    CUDA_CHECK_VOID(cudaMemset(d_long_count, 0, sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMemset(d_frag_count, 0, sizeof(uint32_t)));
    if (n_hits == 0) return 0;

    extend_short_kernel<<<(n_hits + 255) / 256, 256>>>(
        d_filtered_hits, n_hits, seqX, seqY, seqY_rc,
        wl, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
        d_long_seeds, d_long_count, d_out_frags, d_frag_count, max_frags);
    CUDA_CHECK_VOID(cudaGetLastError());

    uint32_t h_nlong = 0;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nlong, d_long_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nlong > 0) {
        const uint32_t warps_per_block = 8;
        uint32_t grid = (h_nlong + warps_per_block - 1) / warps_per_block;
        extend_long_warp_kernel<<<grid, warps_per_block * 32>>>(
            d_long_seeds, h_nlong, seqX, seqY, seqY_rc,
            wl, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
            d_out_frags, d_frag_count, max_frags);
        CUDA_CHECK_VOID(cudaGetLastError());
    }

    uint32_t h_nfrags = 0;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nfrags, d_frag_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nfrags > max_frags) {
        out_overflowed = true;
        return 0;
    }
    return h_nfrags;
}
