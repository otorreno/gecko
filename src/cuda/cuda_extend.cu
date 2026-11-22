/*
 * cuda_extend.cu
 *
 * Stage 3: Two-tier GPU bidirectional ungapped extension matching GECKO's exact
 * scoring model (+4 match, -4 mismatch, -1 for N, stop on 0x00 or fscore < 0):
 *  - Tier 1 (Short, 1 thread/seed, 512 bp cap): filters out 95%+ of noise hits.
 *  - Tier 2 (Long, 1 full 32-thread warp/seed): zero-divergence warp prefix sums.
 */

#include <climits>
#include "cuda_common.cuh"

static __global__ void extend_short_kernel(
    const uint64_t * __restrict__ hits,
    uint32_t n_hits,
    const uint8_t * __restrict__ seqX,
    uint64_t lenX,
    const uint8_t * __restrict__ seqY,
    const uint8_t * __restrict__ seqY_rc,
    uint64_t lenY,
    int wl,
    uint32_t min_len,
    float sim_th,
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

    const uint8_t *dy = is_rev ? seqY_rc : seqY;

    int64_t xfil = (int64_t)posX + wl;
    int64_t ycol = (int64_t)posY + wl;
    int64_t rem_fwd = ((int64_t)lenX - xfil < (int64_t)lenY - ycol) ? ((int64_t)lenX - xfil) : ((int64_t)lenY - ycol);
    int64_t lim_fwd = (rem_fwd < GECKO_SHORT_CAP) ? rem_fwd : GECKO_SHORT_CAP;

    int64_t xfilmax = (int64_t)posX + wl - 1;
    int64_t ycolmax = (int64_t)posY + wl - 1;
    int fscore = GECKO_POINT * wl;
    int fscoreMax = fscore;
    int nIdent = wl, maxIdent = wl;

    bool fwd_terminated = false;
    for (int64_t step = 0; step < lim_fwd; step++) {
        uint8_t cx = seqX[xfil + step] & 0xDF;
        uint8_t cy = dy[ycol + step] & 0xDF;
        if (cx == 0 || cy == 0) { fwd_terminated = true; break; }
        if (cx == 'N' || cy == 'N') {
            fscore -= 1;
        } else if (cx == cy) {
            fscore += GECKO_POINT;
            nIdent++;
            if (fscoreMax < fscore) {
                fscoreMax = fscore;
                xfilmax = xfil + step;
                ycolmax = ycol + step;
                maxIdent = nIdent;
            }
        } else {
            fscore -= GECKO_POINT;
        }
        if (fscore < 0) { fwd_terminated = true; break; }
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
    for (int64_t step = 0; step < lim_bwd; step++) {
        uint8_t cx = seqX[xfil2 - step] & 0xDF;
        uint8_t cy = dy[ycol2 - step] & 0xDF;
        if (cx == 0 || cy == 0) { bwd_terminated = true; break; }
        if (cx == 'N' || cy == 'N') {
            fscore -= 1;
        } else if (cx == cy) {
            fscore += GECKO_POINT;
            nIdent++;
            if (fscoreMax < fscore) {
                fscoreMax = fscore;
                xfilmax2 = xfil2 - step;
                ycolmax2 = ycol2 - step;
                maxIdent = nIdent;
            }
        } else {
            fscore -= GECKO_POINT;
        }
        if (fscore < 0) { bwd_terminated = true; break; }
    }
    if (!bwd_terminated && rem_bwd > GECKO_SHORT_CAP) {
        uint32_t slot = atomicAdd(long_count, 1U);
        long_seeds[slot] = key;
        return;
    }

    uint32_t length = (uint32_t)(xfilmax - xfilmax2 + 1);
    float sim = (float)((double)fscoreMax * 100.0 / (double)((uint64_t)length * GECKO_POINT));
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
    const uint8_t * __restrict__ seqX,
    uint64_t lenX,
    const uint8_t * __restrict__ seqY,
    const uint8_t * __restrict__ seqY_rc,
    uint64_t lenY,
    int wl,
    uint32_t min_len,
    float sim_th,
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

    const uint8_t *dy = is_rev ? seqY_rc : seqY;

    /* Forward pass (32 bases per warp step) */
    int64_t xfil = (int64_t)posX + wl;
    int64_t ycol = (int64_t)posY + wl;
    int64_t rem_fwd = ((int64_t)lenX - xfil < (int64_t)lenY - ycol) ? ((int64_t)lenX - xfil) : ((int64_t)lenY - ycol);

    int64_t xfilmax = (int64_t)posX + wl - 1;
    int64_t ycolmax = (int64_t)posY + wl - 1;
    int fscore = GECKO_POINT * wl;
    int fscoreMax = fscore;
    int nIdent = wl, maxIdent = wl;

    for (int64_t off = 0; off < rem_fwd; off += 32) {
        bool in_bounds = (off + lane < rem_fwd);
        int delta = 0, is_match = 0;
        bool is_sep = !in_bounds;

        if (in_bounds) {
            uint8_t cx = seqX[xfil + off + lane] & 0xDF;
            uint8_t cy = dy[ycol + off + lane] & 0xDF;
            if (cx == 0 || cy == 0) is_sep = true;
            else if (cx == 'N' || cy == 'N') delta = -1;
            else if (cx == cy) { delta = GECKO_POINT; is_match = 1; }
            else delta = -GECKO_POINT;
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

    /* Backward pass (32 bases per warp step, starting with fscore = fscoreMax) */
    int64_t xfil2 = (int64_t)posX - 1;
    int64_t ycol2 = (int64_t)posY - 1;
    int64_t rem_bwd = (posX < posY) ? (int64_t)posX : (int64_t)posY;
    int64_t xfilmax2 = (int64_t)posX;
    int64_t ycolmax2 = (int64_t)posY;
    fscore = fscoreMax;
    nIdent = maxIdent;

    for (int64_t off = 0; off < rem_bwd; off += 32) {
        bool in_bounds = (off + lane < rem_bwd);
        int delta = 0, is_match = 0;
        bool is_sep = !in_bounds;

        if (in_bounds) {
            uint8_t cx = seqX[xfil2 - off - lane] & 0xDF;
            uint8_t cy = dy[ycol2 - off - lane] & 0xDF;
            if (cx == 0 || cy == 0) is_sep = true;
            else if (cx == 'N' || cy == 'N') delta = -1;
            else if (cx == cy) { delta = GECKO_POINT; is_match = 1; }
            else delta = -GECKO_POINT;
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
        float sim = (float)((double)fscoreMax * 100.0 / (double)((uint64_t)length * GECKO_POINT));
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
    const uint8_t *d_seqX,
    uint64_t lenX,
    const uint8_t *d_seqY,
    const uint8_t *d_seqY_rc,
    uint64_t lenY,
    int wl,
    uint32_t min_len,
    float sim_th,
    uint64_t *d_long_seeds,
    uint32_t *d_long_count,
    GpuFrag *d_out_frags,
    uint32_t *d_frag_count,
    uint32_t max_frags
) {
    CUDA_CHECK(cudaMemset(d_long_count, 0, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemset(d_frag_count, 0, sizeof(uint32_t)));
    if (n_hits == 0) return 0;

    extend_short_kernel<<<(n_hits + 255) / 256, 256>>>(
        d_filtered_hits, n_hits, d_seqX, lenX, d_seqY, d_seqY_rc, lenY,
        wl, min_len, sim_th,
        d_long_seeds, d_long_count, d_out_frags, d_frag_count, max_frags);
    CUDA_CHECK(cudaGetLastError());

    uint32_t h_nlong = 0;
    CUDA_CHECK(cudaMemcpy(&h_nlong, d_long_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nlong > 0) {
        const uint32_t warps_per_block = 8;
        uint32_t grid = (h_nlong + warps_per_block - 1) / warps_per_block;
        extend_long_warp_kernel<<<grid, warps_per_block * 32>>>(
            d_long_seeds, h_nlong, d_seqX, lenX, d_seqY, d_seqY_rc, lenY,
            wl, min_len, sim_th,
            d_out_frags, d_frag_count, max_frags);
        CUDA_CHECK(cudaGetLastError());
    }

    uint32_t h_nfrags = 0;
    CUDA_CHECK(cudaMemcpy(&h_nfrags, d_frag_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    return (h_nfrags > max_frags) ? max_frags : h_nfrags;
}
