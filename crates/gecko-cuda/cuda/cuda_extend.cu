/*
 * cuda_extend.cu
 *
 * Stage 3: Zero-stack two-tier 2-bit GPU bidirectional ungapped extension
 * matching GECKO's exact scoring model (+match, -mismatch, -n_penalty,
 * SEQ_SEP, rare IUPAC):
 *  - F12: Zero local stack frame (0 B stack, 0 STL/LDL across sm_75/80/86/89)
 *    via explicit register pointers + deferred cold special metadata selects,
 *    and guarded FP64 similarity divide behind (length > min_len && fscoreMax > 0).
 *  - F13: O(1) bounded diagonal segmentation (SEG_CAP=64) with provably lossless
 *    register covered_x skipping (eliminates O(L^2/d) redundant long-HSP work),
 *    Fast-path 2 popcnt decay shortcut (order-independent on backward), and
 *    1,024-bp coalesced warp stepping in extend_long_warp_kernel.
 *  - F14/F15: Tile-relative HitKeyLayout unpacking and 24-byte GpuFrag output.
 *  - W3-F13c: __launch_bounds__(256, 4) on both kernels (keeping score_special_pair
 *    __forceinline__) trims extend_long_warp_kernel from 76 -> <= 64 regs (100% T4
 *    occupancy) with 0 B stack / 0 STL / 0 LDL across sm_75/80/86/89.
 *
 * NOTE: NEVER #include <cub/cub.cuh> in this translation unit so standalone
 * static verification inspects strictly the extension kernels.
 */

#include <climits>
#include "cuda_common.cuh"

static __device__ __forceinline__ int score_special_pair(
    const uint64_t * __restrict__ sx_bits,
    const uint32_t * __restrict__ sx_sep,
    uint32_t sx_nsep,
    const uint32_t * __restrict__ sx_iupac_pos,
    const uint8_t  * __restrict__ sx_iupac_byte,
    uint32_t sx_niupac,
    uint32_t px,
    bool sp_x,
    const uint64_t * __restrict__ dy_bits,
    const uint32_t * __restrict__ dy_sep,
    uint32_t dy_nsep,
    const uint32_t * __restrict__ dy_iupac_pos,
    const uint8_t  * __restrict__ dy_iupac_byte,
    uint32_t dy_niupac,
    uint32_t py,
    bool sp_y,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    bool &is_sep,
    int &is_match
) {
    uint8_t cx = sp_x
                     ? decode_special_byte_raw(sx_sep, sx_nsep, sx_iupac_pos, sx_iupac_byte, sx_niupac, px)
                     : base_2bit_to_ascii(base_2bit_dev(sx_bits, px));
    uint8_t cy = sp_y
                     ? decode_special_byte_raw(dy_sep, dy_nsep, dy_iupac_pos, dy_iupac_byte, dy_niupac, py)
                     : base_2bit_to_ascii(base_2bit_dev(dy_bits, py));
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

static __device__ __forceinline__ int warp_inclusive_sum(int val, int lane) {
    #pragma unroll
    for (int off = 1; off < 32; off <<= 1) {
        int n = __shfl_up_sync(0xFFFFFFFFU, val, off);
        if (lane >= off) val += n;
    }
    return val;
}

static __device__ __forceinline__ int warp_sum(int val) {
    #pragma unroll
    for (int off = 16; off > 0; off >>= 1) {
        val += __shfl_xor_sync(0xFFFFFFFFU, val, off);
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

static __device__ __forceinline__ uint64_t pack_long_seg(uint32_t seg_start, uint32_t seg_len) {
    return ((uint64_t)seg_start << 32) | (uint64_t)seg_len;
}

static __global__ __launch_bounds__(256, 4) void extend_short_kernel(
    const uint64_t * __restrict__ hits,
    uint32_t n_hits,
    HitKeyLayout key_layout,
    GpuSeqView seqX,
    GpuSeqView seqY,
    GpuSeqView seqY_rc,
    int wl,
    bool is_spaced,
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

    /* F13: O(1) fixed-alignment bounded diagonal segmentation (SEG_CAP = 64) */
    uint64_t key0 = hits[idx];
    uint64_t sd0 = hit_strand_diag_tile(key_layout, key0);
    bool same_prev = (idx > 0U) && (hit_strand_diag_tile(key_layout, hits[idx - 1U]) == sd0);
    if (same_prev && ((idx & 63U) != 0U)) return;

    uint32_t max_end = (idx & ~63U) + 64U;
    if (max_end > n_hits) max_end = n_hits;
    uint32_t seg_end = idx + 1U;
    if (seg_end < max_end && hit_strand_diag_tile(key_layout, hits[seg_end]) == sd0) {
        uint32_t lo = seg_end + 1U, hi = max_end;
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (hit_strand_diag_tile(key_layout, hits[mid]) == sd0) lo = mid + 1U;
            else hi = mid;
        }
        seg_end = lo;
    }

    /* F12: Select ONLY the 3 hot sequence fields into registers (zero stack frame!) */
    uint32_t is_rev0 = (uint32_t)(key0 >> key_layout.strand_shift);
    const uint64_t * __restrict__ sx_bits = seqX.bits;
    const uint32_t * __restrict__ sx_sp   = seqX.special_bm;
    const uint64_t * __restrict__ dy_bits = is_rev0 ? seqY_rc.bits       : seqY.bits;
    const uint32_t * __restrict__ dy_sp   = is_rev0 ? seqY_rc.special_bm : seqY.special_bm;
    int32_t sx_len = (int32_t)seqX.len;
    int32_t dy_len = (int32_t)(is_rev0 ? seqY_rc.len : seqY.len);

    int32_t covered_x = -1;
    int32_t covered_min_x = 0;

    for (uint32_t k = idx; k < seg_end; k++) {
        uint64_t key = (k == idx) ? key0 : hits[k];
        uint32_t is_rev, posX, posY;
        int64_t diag_unused;
        unpack_hit_key_tile(key_layout, key, is_rev, posX, posY, diag_unused);

        if ((int32_t)posX <= covered_x && (!is_spaced || (int32_t)posX >= covered_min_x)) {
            continue;
        }

        int32_t xfilmax = (int32_t)posX + wl - 1;
        int nIdent = wl;
        int fscore = match_score * wl;
        if (is_spaced) {
            uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
            uint64_t diff = (load_32_2bit(sx_bits, posX) ^ load_32_2bit(dy_bits, posY)) & span_mask;
            int mm = __popcll((diff | (diff >> 1)) & 0x5555555555555555ULL);
            nIdent = wl - mm;
            fscore = nIdent * match_score - mm * mismatch_penalty;
        }
        const int seed_fscore = fscore;
        int fscoreMax = fscore;
        int maxIdent = nIdent;

        /* ---- Forward Extension (up to GECKO_SHORT_CAP) ---- */
        {
            int32_t xfil = (int32_t)posX + wl;
            int32_t ycol = (int32_t)posY + wl;
            int32_t rem_x = sx_len - xfil;
            int32_t rem_y = dy_len - ycol;
            int32_t rem_fwd = (rem_x < rem_y) ? rem_x : rem_y;
            int32_t lim_fwd = (rem_fwd < GECKO_SHORT_CAP) ? rem_fwd : GECKO_SHORT_CAP;

            bool fwd_terminated = false;
            for (int32_t off = 0; off < lim_fwd && !fwd_terminated; off += 32) {
                int chunk = (lim_fwd - off < 32) ? (int)(lim_fwd - off) : 32;
                uint32_t px0 = (uint32_t)(xfil + off);
                uint32_t py0 = (uint32_t)(ycol + off);
                uint64_t wx = load_32_2bit(sx_bits, px0);
                uint64_t wy = load_32_2bit(dy_bits, py0);
                uint32_t spx = load_32_1bit(sx_sp, px0);
                uint32_t spy = load_32_1bit(dy_sp, py0);
                uint32_t sp_any = spx | spy;

                /* Fast-path 1: 32-bp clean exact match */
                if (chunk == 32 && sp_any == 0U && wx == wy && fscore + match_score >= 0) {
                    fscore += 32 * match_score;
                    nIdent += 32;
                    if (fscoreMax < fscore) {
                        fscoreMax = fscore;
                        xfilmax = xfil + off + 31;
                        maxIdent = nIdent;
                    }
                    continue;
                }

                /* F13 Fast-path 2: 32-bp clean popcnt decay shortcut */
                if (chunk == 32 && sp_any == 0U) {
                    uint64_t diff = wx ^ wy;
                    int mm = __popcll((diff | (diff >> 1)) & 0x5555555555555555ULL);
                    int matches = 32 - mm;
                    if ((int64_t)fscore >= (int64_t)mm * (int64_t)mismatch_penalty &&
                        (int64_t)fscore <= (int64_t)fscoreMax - (int64_t)matches * (int64_t)match_score) {
                        fscore += matches * match_score - mm * mismatch_penalty;
                        nIdent += matches;
                        continue;
                    }
                }

                if (sp_any == 0U) {
                    for (int s = 0; s < chunk; s++) {
                        if ((wx & 3ULL) == (wy & 3ULL)) {
                            fscore += match_score;
                            nIdent++;
                            if (fscoreMax < fscore) {
                                fscoreMax = fscore;
                                xfilmax = xfil + off + s;
                                maxIdent = nIdent;
                            }
                        } else {
                            fscore -= mismatch_penalty;
                        }
                        if (fscore < 0) {
                            fwd_terminated = true;
                            break;
                        }
                        wx >>= 2;
                        wy >>= 2;
                    }
                } else {
                    for (int s = 0; s < chunk; s++) {
                        if (((sp_any >> s) & 1U) == 0U) {
                            if ((wx & 3ULL) == (wy & 3ULL)) {
                                fscore += match_score;
                                nIdent++;
                                if (fscoreMax < fscore) {
                                    fscoreMax = fscore;
                                    xfilmax = xfil + off + s;
                                    maxIdent = nIdent;
                                }
                            } else {
                                fscore -= mismatch_penalty;
                            }
                        } else {
                            bool is_sep = false;
                            int is_match = 0;
                            int delta = score_special_pair(
                                sx_bits, seqX.sep_pos, seqX.n_seps, seqX.iupac_pos, seqX.iupac_byte, seqX.n_iupac,
                                px0 + (uint32_t)s, ((spx >> s) & 1U) != 0U,
                                dy_bits,
                                is_rev0 ? seqY_rc.sep_pos    : seqY.sep_pos,
                                is_rev0 ? seqY_rc.n_seps     : seqY.n_seps,
                                is_rev0 ? seqY_rc.iupac_pos  : seqY.iupac_pos,
                                is_rev0 ? seqY_rc.iupac_byte : seqY.iupac_byte,
                                is_rev0 ? seqY_rc.n_iupac    : seqY.n_iupac,
                                py0 + (uint32_t)s, ((spy >> s) & 1U) != 0U,
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
                                    xfilmax = xfil + off + s;
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
                long_seeds[slot] = pack_long_seg(k, seg_end - k);
                return;
            }
        }

        /* ---- Backward Extension (up to GECKO_SHORT_CAP) ---- */
        int32_t xfilmax2 = (int32_t)posX;
        {
            int32_t xfil2 = (int32_t)posX - 1;
            int32_t ycol2 = (int32_t)posY - 1;
            int32_t rem_bwd = (posX < posY) ? (int32_t)posX : (int32_t)posY;
            int32_t lim_bwd = (rem_bwd < GECKO_SHORT_CAP) ? rem_bwd : GECKO_SHORT_CAP;

            fscore = fscoreMax;
            nIdent = maxIdent;

            bool bwd_terminated = false;
            for (int32_t off = 0; off < lim_bwd && !bwd_terminated; off += 32) {
                int32_t rem_chunk = lim_bwd - off;
                if (rem_chunk >= 32) {
                    uint32_t bx0 = (uint32_t)(xfil2 - off - 31);
                    uint32_t by0 = (uint32_t)(ycol2 - off - 31);
                    uint64_t wx_raw = load_32_2bit(sx_bits, bx0);
                    uint64_t wy_raw = load_32_2bit(dy_bits, by0);
                    uint32_t spx_raw = load_32_1bit(sx_sp, bx0);
                    uint32_t spy_raw = load_32_1bit(dy_sp, by0);
                    uint32_t sp_any_raw = spx_raw | spy_raw;

                    /* Order-independent Fast-path 1 (zero reversals!) */
                    if (sp_any_raw == 0U && wx_raw == wy_raw && fscore + match_score >= 0) {
                        fscore += 32 * match_score;
                        nIdent += 32;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax2 = xfil2 - off - 31;
                            maxIdent = nIdent;
                        }
                        continue;
                    }

                    if (sp_any_raw == 0U) {
                        uint64_t diff_raw = wx_raw ^ wy_raw;
                        uint64_t mm_raw = (diff_raw | (diff_raw >> 1)) & 0x5555555555555555ULL;
                        int mm = __popcll(mm_raw);
                        int matches = 32 - mm;
                        /* Order-independent Fast-path 2 popcnt decay (zero reversals!) */
                        if ((int64_t)fscore >= (int64_t)mm * (int64_t)mismatch_penalty &&
                            (int64_t)fscore <= (int64_t)fscoreMax - (int64_t)matches * (int64_t)match_score) {
                            fscore += matches * match_score - mm * mismatch_penalty;
                            nIdent += matches;
                            continue;
                        }
                        /* Single reversal for clean slow path */
                        uint64_t mm_rev = reverse_32_2bit_dev(mm_raw);
                        for (int s = 0; s < 32; s++) {
                            if ((mm_rev & 1ULL) == 0ULL) {
                                fscore += match_score;
                                nIdent++;
                                if (fscoreMax < fscore) {
                                    fscoreMax = fscore;
                                    xfilmax2 = xfil2 - off - s;
                                    maxIdent = nIdent;
                                }
                            } else {
                                fscore -= mismatch_penalty;
                            }
                            if (fscore < 0) {
                                bwd_terminated = true;
                                break;
                            }
                            mm_rev >>= 2;
                        }
                    } else {
                        uint64_t wx = reverse_32_2bit_dev(wx_raw);
                        uint64_t wy = reverse_32_2bit_dev(wy_raw);
                        uint32_t spx = __brev(spx_raw);
                        uint32_t spy = __brev(spy_raw);
                        uint32_t sp_any = spx | spy;
                        for (int s = 0; s < 32; s++) {
                            if (((sp_any >> s) & 1U) == 0U) {
                                if ((wx & 3ULL) == (wy & 3ULL)) {
                                    fscore += match_score;
                                    nIdent++;
                                    if (fscoreMax < fscore) {
                                        fscoreMax = fscore;
                                        xfilmax2 = xfil2 - off - s;
                                        maxIdent = nIdent;
                                    }
                                } else {
                                    fscore -= mismatch_penalty;
                                }
                            } else {
                                bool is_sep = false;
                                int is_match = 0;
                                int delta = score_special_pair(
                                    sx_bits, seqX.sep_pos, seqX.n_seps, seqX.iupac_pos, seqX.iupac_byte, seqX.n_iupac,
                                    (uint32_t)(xfil2 - off - s), ((spx >> s) & 1U) != 0U,
                                    dy_bits,
                                    is_rev0 ? seqY_rc.sep_pos    : seqY.sep_pos,
                                    is_rev0 ? seqY_rc.n_seps     : seqY.n_seps,
                                    is_rev0 ? seqY_rc.iupac_pos  : seqY.iupac_pos,
                                    is_rev0 ? seqY_rc.iupac_byte : seqY.iupac_byte,
                                    is_rev0 ? seqY_rc.n_iupac    : seqY.n_iupac,
                                    (uint32_t)(ycol2 - off - s), ((spy >> s) & 1U) != 0U,
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
                                        xfilmax2 = xfil2 - off - s;
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
                    for (int s = 0; s < chunk; s++) {
                        uint32_t px = (uint32_t)(xfil2 - off - s);
                        uint32_t py = (uint32_t)(ycol2 - off - s);
                        bool sp_x = is_bit_set(sx_sp, px);
                        bool sp_y = is_bit_set(dy_sp, py);
                        if (!sp_x && !sp_y) {
                            if (base_2bit_dev(sx_bits, px) == base_2bit_dev(dy_bits, py)) {
                                fscore += match_score;
                                nIdent++;
                                if (fscoreMax < fscore) {
                                    fscoreMax = fscore;
                                    xfilmax2 = (int32_t)px;
                                    maxIdent = nIdent;
                                }
                            } else {
                                fscore -= mismatch_penalty;
                            }
                        } else {
                            bool is_sep = false;
                            int is_match = 0;
                            int delta = score_special_pair(
                                sx_bits, seqX.sep_pos, seqX.n_seps, seqX.iupac_pos, seqX.iupac_byte, seqX.n_iupac,
                                px, sp_x,
                                dy_bits,
                                is_rev0 ? seqY_rc.sep_pos    : seqY.sep_pos,
                                is_rev0 ? seqY_rc.n_seps     : seqY.n_seps,
                                is_rev0 ? seqY_rc.iupac_pos  : seqY.iupac_pos,
                                is_rev0 ? seqY_rc.iupac_byte : seqY.iupac_byte,
                                is_rev0 ? seqY_rc.n_iupac    : seqY.n_iupac,
                                py, sp_y,
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
                                    xfilmax2 = (int32_t)px;
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
                long_seeds[slot] = pack_long_seg(k, seg_end - k);
                return;
            }
        }

        /* F12 + F13: Guarded FP64 divide + provably lossless covered_x update */
        uint32_t length = (uint32_t)(xfilmax - xfilmax2 + 1);
        if (length > min_len) {
            bool accepted = false;
            if (fscoreMax > 0) {
                float sim = (float)((double)fscoreMax * 100.0 /
                                    (double)((uint64_t)length * (uint64_t)match_score));
                if (sim > sim_th) {
                    accepted = true;
                    uint32_t slot = atomicAdd(frag_count, 1U);
                    if (slot < max_frags) {
                        uint32_t yStart = (uint32_t)((int32_t)posY + (xfilmax2 - (int32_t)posX));
                        out_frags[slot].xStart = (uint32_t)xfilmax2;
                        out_frags[slot].yStart = yStart;
                        out_frags[slot].length = length;
                        out_frags[slot].ident  = (uint32_t)maxIdent;
                        out_frags[slot].score  = (uint32_t)fscoreMax;
                        out_frags[slot].flags  = is_rev0 & 1U;
                    }
                }
            }
            if (accepted && seed_fscore >= 0) {
                int32_t new_cov = is_spaced ? (xfilmax - (wl - 1)) : xfilmax;
                if (new_cov > covered_x) {
                    covered_x = new_cov;
                    covered_min_x = is_spaced ? ((int32_t)posX + wl) : 0;
                }
            }
        } else if (seed_fscore >= 0) {
            int32_t new_cov = is_spaced ? (xfilmax - (wl - 1)) : xfilmax;
            int32_t new_min = is_spaced ? ((int32_t)posX + wl) : 0;
            if (new_cov > covered_x && (!is_spaced || new_cov >= new_min)) {
                covered_x = new_cov;
                covered_min_x = new_min;
            }
        }
    }
}

static __global__ __launch_bounds__(256, 4) void extend_long_warp_kernel(
    const uint64_t * __restrict__ hits,
    const uint64_t * __restrict__ long_segs,
    uint32_t n_long,
    HitKeyLayout key_layout,
    GpuSeqView seqX,
    GpuSeqView seqY,
    GpuSeqView seqY_rc,
    int wl,
    bool is_spaced,
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

    uint64_t seg = long_segs[warp_id];
    uint32_t seg_start = (uint32_t)(seg >> 32);
    uint32_t seg_len   = (uint32_t)(seg & 0xFFFFFFFFULL);
    if (seg_len == 0U) return;

    uint64_t key0 = hits[seg_start];
    uint32_t is_rev0 = (uint32_t)(key0 >> key_layout.strand_shift);

    /* F12: Explicit register pointers (0 B stack frame, 0 STL/LDL!) */
    const uint64_t * __restrict__ sx_bits = seqX.bits;
    const uint32_t * __restrict__ sx_sp   = seqX.special_bm;
    const uint64_t * __restrict__ dy_bits = is_rev0 ? seqY_rc.bits       : seqY.bits;
    const uint32_t * __restrict__ dy_sp   = is_rev0 ? seqY_rc.special_bm : seqY.special_bm;
    int32_t sx_len = (int32_t)seqX.len;
    int32_t dy_len = (int32_t)(is_rev0 ? seqY_rc.len : seqY.len);

    int32_t covered_x = -1;
    int32_t covered_min_x = 0;

    for (uint32_t m = 0U; m < seg_len; m++) {
        uint64_t key = (m == 0U) ? key0 : hits[seg_start + m];
        uint32_t is_rev, posX, posY;
        int64_t diag_unused;
        unpack_hit_key_tile(key_layout, key, is_rev, posX, posY, diag_unused);

        if ((int32_t)posX <= covered_x && (!is_spaced || (int32_t)posX >= covered_min_x)) {
            continue;
        }

        int32_t xfilmax = (int32_t)posX + wl - 1;
        int nIdent = wl;
        int fscore = match_score * wl;
        if (is_spaced) {
            uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
            uint64_t diff = (load_32_2bit(sx_bits, posX) ^ load_32_2bit(dy_bits, posY)) & span_mask;
            int mm = __popcll((diff | (diff >> 1)) & 0x5555555555555555ULL);
            nIdent = wl - mm;
            fscore = nIdent * match_score - mm * mismatch_penalty;
        }
        const int seed_fscore = fscore;
        int fscoreMax = fscore;
        int maxIdent = nIdent;

        /* ---- Forward Extension: Coalesced 1,024-bp + 32-bp Warp Stepping ---- */
        {
            int32_t xfil = (int32_t)posX + wl;
            int32_t ycol = (int32_t)posY + wl;
            int32_t rem_x = sx_len - xfil;
            int32_t rem_y = dy_len - ycol;
            int32_t rem_fwd = (rem_x < rem_y) ? rem_x : rem_y;

            int32_t off = 0;
            while (off < rem_fwd) {
                int32_t rem_chunk = rem_fwd - off;
                uint64_t wx = 0ULL, wy = 0ULL;
                uint32_t spx = 0U, spy = 0U;

                if (rem_chunk >= 1024) {
                    uint32_t l_px = (uint32_t)(xfil + off + (lane << 5));
                    uint32_t l_py = (uint32_t)(ycol + off + (lane << 5));
                    uint64_t l_wx = load_32_2bit(sx_bits, l_px);
                    uint64_t l_wy = load_32_2bit(dy_bits, l_py);
                    uint32_t l_spx = load_32_1bit(sx_sp, l_px);
                    uint32_t l_spy = load_32_1bit(dy_sp, l_py);
                    uint32_t l_sp_any = l_spx | l_spy;

                    /* Fast-path A: 1,024-bp clean exact match */
                    uint32_t non_exact = __ballot_sync(0xFFFFFFFFU, (l_sp_any != 0U) || (l_wx != l_wy));
                    if (non_exact == 0U && fscore + match_score >= 0) {
                        fscore += 1024 * match_score;
                        nIdent += 1024;
                        if (fscore > fscoreMax) {
                            fscoreMax = fscore;
                            xfilmax = xfil + off + 1023;
                            maxIdent = nIdent;
                        }
                        off += 1024;
                        continue;
                    }

                    /* Fast-path B: 1,024-bp clean coalesced decay tail */
                    uint32_t sp_mask = __ballot_sync(0xFFFFFFFFU, l_sp_any != 0U);
                    if (sp_mask == 0U &&
                        (int64_t)fscore >= 1024LL * (int64_t)mismatch_penalty &&
                        (int64_t)fscore <= (int64_t)fscoreMax - 1024LL * (int64_t)match_score) {
                        uint64_t l_diff = l_wx ^ l_wy;
                        int l_mm = __popcll((l_diff | (l_diff >> 1)) & 0x5555555555555555ULL);
                        int total_mm = warp_sum(l_mm);
                        int total_matches = 1024 - total_mm;
                        fscore += total_matches * match_score - total_mm * mismatch_penalty;
                        nIdent += total_matches;
                        off += 1024;
                        continue;
                    }

                    /* Prefix-skip exact 32-bp chunks up to first_bad, broadcast chunk first_bad */
                    int first_bad = (non_exact != 0U) ? (__ffs(non_exact) - 1) : 0;
                    if (first_bad > 0 && fscore + match_score >= 0) {
                        int exact_bp = first_bad << 5;
                        fscore += exact_bp * match_score;
                        nIdent += exact_bp;
                        if (fscore > fscoreMax) {
                            fscoreMax = fscore;
                            xfilmax = xfil + off + exact_bp - 1;
                            maxIdent = nIdent;
                        }
                        off += exact_bp;
                    } else {
                        first_bad = 0;
                    }
                    wx  = __shfl_sync(0xFFFFFFFFU, l_wx,  first_bad);
                    wy  = __shfl_sync(0xFFFFFFFFU, l_wy,  first_bad);
                    spx = __shfl_sync(0xFFFFFFFFU, l_spx, first_bad);
                    spy = __shfl_sync(0xFFFFFFFFU, l_spy, first_bad);
                    rem_chunk = rem_fwd - off;
                } else {
                    if (lane == 0) {
                        uint32_t px0 = (uint32_t)(xfil + off);
                        uint32_t py0 = (uint32_t)(ycol + off);
                        wx = load_32_2bit(sx_bits, px0);
                        wy = load_32_2bit(dy_bits, py0);
                        spx = load_32_1bit(sx_sp, px0);
                        spy = load_32_1bit(dy_sp, py0);
                    }
                    wx  = __shfl_sync(0xFFFFFFFFU, wx,  0);
                    wy  = __shfl_sync(0xFFFFFFFFU, wy,  0);
                    spx = __shfl_sync(0xFFFFFFFFU, spx, 0);
                    spy = __shfl_sync(0xFFFFFFFFU, spy, 0);
                }

                uint32_t sp_any = spx | spy;
                if (rem_chunk >= 32 && sp_any == 0U && wx == wy && fscore + match_score >= 0) {
                    fscore += 32 * match_score;
                    nIdent += 32;
                    if (fscore > fscoreMax) {
                        fscoreMax = fscore;
                        xfilmax = xfil + off + 31;
                        maxIdent = nIdent;
                    }
                    off += 32;
                    continue;
                }
                if (rem_chunk >= 32 && sp_any == 0U) {
                    uint64_t diff = wx ^ wy;
                    int mm = __popcll((diff | (diff >> 1)) & 0x5555555555555555ULL);
                    int matches = 32 - mm;
                    if ((int64_t)fscore >= (int64_t)mm * (int64_t)mismatch_penalty &&
                        (int64_t)fscore <= (int64_t)fscoreMax - (int64_t)matches * (int64_t)match_score) {
                        fscore += matches * match_score - mm * mismatch_penalty;
                        nIdent += matches;
                        off += 32;
                        continue;
                    }
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
                            sx_bits, seqX.sep_pos, seqX.n_seps, seqX.iupac_pos, seqX.iupac_byte, seqX.n_iupac,
                            (uint32_t)(xfil + off + lane), ((spx >> lane) & 1U) != 0U,
                            dy_bits,
                            is_rev0 ? seqY_rc.sep_pos    : seqY.sep_pos,
                            is_rev0 ? seqY_rc.n_seps     : seqY.n_seps,
                            is_rev0 ? seqY_rc.iupac_pos  : seqY.iupac_pos,
                            is_rev0 ? seqY_rc.iupac_byte : seqY.iupac_byte,
                            is_rev0 ? seqY_rc.n_iupac    : seqY.n_iupac,
                            (uint32_t)(ycol + off + lane), ((spy >> lane) & 1U) != 0U,
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
                int active_end = (sep_lane <= neg_lane) ? sep_lane : (neg_lane + 1);

                int cand_score = (lane < active_end) ? score_here : INT_MIN;
                int chunk_max = warp_max(cand_score);
                if (chunk_max > fscoreMax) {
                    uint32_t max_mask = __ballot_sync(0xFFFFFFFFU, (lane < active_end) && (score_here == chunk_max));
                    int best_lane = __ffs(max_mask) - 1;
                    int best_pmatch = __shfl_sync(0xFFFFFFFFU, pmatch, best_lane);
                    fscoreMax = chunk_max;
                    xfilmax = xfil + off + best_lane;
                    maxIdent = nIdent + best_pmatch;
                }

                if (sep_mask || neg_mask) break;

                fscore = __shfl_sync(0xFFFFFFFFU, score_here, 31);
                nIdent += __shfl_sync(0xFFFFFFFFU, pmatch, 31);
                off += 32;
            }
        }

        /* ---- Backward Extension: Coalesced 1,024-bp + 32-bp Warp Stepping (rlane = 31 - lane) ---- */
        int32_t xfilmax2 = (int32_t)posX;
        {
            int32_t xfil2 = (int32_t)posX - 1;
            int32_t ycol2 = (int32_t)posY - 1;
            int32_t rem_bwd = (posX < posY) ? (int32_t)posX : (int32_t)posY;
            fscore = fscoreMax;
            nIdent = maxIdent;

            int rlane = 31 - lane;
            int32_t off = 0;
            while (off < rem_bwd) {
                int32_t rem_chunk = rem_bwd - off;
                int delta = 0, is_match = 0;
                bool is_sep = false;

                if (rem_chunk >= 32) {
                    uint64_t wx = 0ULL, wy = 0ULL;
                    uint32_t spx = 0U, spy = 0U;

                    if (rem_chunk >= 1024) {
                        uint32_t bx_l = (uint32_t)(xfil2 - off - ((lane + 1) << 5) + 1);
                        uint32_t by_l = (uint32_t)(ycol2 - off - ((lane + 1) << 5) + 1);
                        uint64_t l_wx = load_32_2bit(sx_bits, bx_l);
                        uint64_t l_wy = load_32_2bit(dy_bits, by_l);
                        uint32_t l_spx = load_32_1bit(sx_sp, bx_l);
                        uint32_t l_spy = load_32_1bit(dy_sp, by_l);
                        uint32_t l_sp_any = l_spx | l_spy;

                        /* Fast-path A: 1,024-bp clean exact match */
                        uint32_t non_exact = __ballot_sync(0xFFFFFFFFU, (l_sp_any != 0U) || (l_wx != l_wy));
                        if (non_exact == 0U && fscore + match_score >= 0) {
                            fscore += 1024 * match_score;
                            nIdent += 1024;
                            if (fscore > fscoreMax) {
                                fscoreMax = fscore;
                                xfilmax2 = xfil2 - off - 1023;
                                maxIdent = nIdent;
                            }
                            off += 1024;
                            continue;
                        }

                        /* Fast-path B: 1,024-bp clean coalesced decay tail */
                        uint32_t sp_mask = __ballot_sync(0xFFFFFFFFU, l_sp_any != 0U);
                        if (sp_mask == 0U &&
                            (int64_t)fscore >= 1024LL * (int64_t)mismatch_penalty &&
                            (int64_t)fscore <= (int64_t)fscoreMax - 1024LL * (int64_t)match_score) {
                            uint64_t l_diff = l_wx ^ l_wy;
                            int l_mm = __popcll((l_diff | (l_diff >> 1)) & 0x5555555555555555ULL);
                            int total_mm = warp_sum(l_mm);
                            int total_matches = 1024 - total_mm;
                            fscore += total_matches * match_score - total_mm * mismatch_penalty;
                            nIdent += total_matches;
                            off += 1024;
                            continue;
                        }

                        int first_bad = (non_exact != 0U) ? (__ffs(non_exact) - 1) : 0;
                        if (first_bad > 0 && fscore + match_score >= 0) {
                            int exact_bp = first_bad << 5;
                            fscore += exact_bp * match_score;
                            nIdent += exact_bp;
                            if (fscore > fscoreMax) {
                                fscoreMax = fscore;
                                xfilmax2 = xfil2 - off - exact_bp + 1;
                                maxIdent = nIdent;
                            }
                            off += exact_bp;
                        } else {
                            first_bad = 0;
                        }
                        wx  = __shfl_sync(0xFFFFFFFFU, l_wx,  first_bad);
                        wy  = __shfl_sync(0xFFFFFFFFU, l_wy,  first_bad);
                        spx = __shfl_sync(0xFFFFFFFFU, l_spx, first_bad);
                        spy = __shfl_sync(0xFFFFFFFFU, l_spy, first_bad);
                    } else {
                        if (lane == 0) {
                            uint32_t bx0 = (uint32_t)(xfil2 - off - 31);
                            uint32_t by0 = (uint32_t)(ycol2 - off - 31);
                            wx = load_32_2bit(sx_bits, bx0);
                            wy = load_32_2bit(dy_bits, by0);
                            spx = load_32_1bit(sx_sp, bx0);
                            spy = load_32_1bit(dy_sp, by0);
                        }
                        wx  = __shfl_sync(0xFFFFFFFFU, wx,  0);
                        wy  = __shfl_sync(0xFFFFFFFFU, wy,  0);
                        spx = __shfl_sync(0xFFFFFFFFU, spx, 0);
                        spy = __shfl_sync(0xFFFFFFFFU, spy, 0);
                    }

                    uint32_t sp_any = spx | spy;
                    if (sp_any == 0U && wx == wy && fscore + match_score >= 0) {
                        fscore += 32 * match_score;
                        nIdent += 32;
                        if (fscore > fscoreMax) {
                            fscoreMax = fscore;
                            xfilmax2 = xfil2 - off - 31;
                            maxIdent = nIdent;
                        }
                        off += 32;
                        continue;
                    }
                    if (sp_any == 0U) {
                        uint64_t diff = wx ^ wy;
                        int mm = __popcll((diff | (diff >> 1)) & 0x5555555555555555ULL);
                        int matches = 32 - mm;
                        if ((int64_t)fscore >= (int64_t)mm * (int64_t)mismatch_penalty &&
                            (int64_t)fscore <= (int64_t)fscoreMax - (int64_t)matches * (int64_t)match_score) {
                            fscore += matches * match_score - mm * mismatch_penalty;
                            nIdent += matches;
                            off += 32;
                            continue;
                        }
                    }

                    /* Zero bit-reversal trick: forward-order words indexed via rlane = 31 - lane */
                    if (((sp_any >> rlane) & 1U) == 0U) {
                        uint32_t bx = (uint32_t)((wx >> (rlane << 1)) & 3ULL);
                        uint32_t by = (uint32_t)((wy >> (rlane << 1)) & 3ULL);
                        if (bx == by) { delta = match_score; is_match = 1; }
                        else { delta = -mismatch_penalty; }
                    } else {
                        delta = score_special_pair(
                            sx_bits, seqX.sep_pos, seqX.n_seps, seqX.iupac_pos, seqX.iupac_byte, seqX.n_iupac,
                            (uint32_t)(xfil2 - off - lane), ((spx >> rlane) & 1U) != 0U,
                            dy_bits,
                            is_rev0 ? seqY_rc.sep_pos    : seqY.sep_pos,
                            is_rev0 ? seqY_rc.n_seps     : seqY.n_seps,
                            is_rev0 ? seqY_rc.iupac_pos  : seqY.iupac_pos,
                            is_rev0 ? seqY_rc.iupac_byte : seqY.iupac_byte,
                            is_rev0 ? seqY_rc.n_iupac    : seqY.n_iupac,
                            (uint32_t)(ycol2 - off - lane), ((spy >> rlane) & 1U) != 0U,
                            match_score, mismatch_penalty, n_penalty, is_sep, is_match);
                    }
                } else {
                    bool in_bounds = (lane < rem_chunk);
                    is_sep = !in_bounds;
                    if (in_bounds) {
                        uint32_t px = (uint32_t)(xfil2 - off - lane);
                        uint32_t py = (uint32_t)(ycol2 - off - lane);
                        bool sp_x = is_bit_set(sx_sp, px);
                        bool sp_y = is_bit_set(dy_sp, py);
                        if (!sp_x && !sp_y) {
                            if (base_2bit_dev(sx_bits, px) == base_2bit_dev(dy_bits, py)) {
                                delta = match_score;
                                is_match = 1;
                            } else {
                                delta = -mismatch_penalty;
                            }
                        } else {
                            delta = score_special_pair(
                                sx_bits, seqX.sep_pos, seqX.n_seps, seqX.iupac_pos, seqX.iupac_byte, seqX.n_iupac,
                                px, sp_x,
                                dy_bits,
                                is_rev0 ? seqY_rc.sep_pos    : seqY.sep_pos,
                                is_rev0 ? seqY_rc.n_seps     : seqY.n_seps,
                                is_rev0 ? seqY_rc.iupac_pos  : seqY.iupac_pos,
                                is_rev0 ? seqY_rc.iupac_byte : seqY.iupac_byte,
                                is_rev0 ? seqY_rc.n_iupac    : seqY.n_iupac,
                                py, sp_y,
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
                int active_end = (sep_lane <= neg_lane) ? sep_lane : (neg_lane + 1);

                int cand_score = (lane < active_end) ? score_here : INT_MIN;
                int chunk_max = warp_max(cand_score);
                if (chunk_max > fscoreMax) {
                    uint32_t max_mask = __ballot_sync(0xFFFFFFFFU, (lane < active_end) && (score_here == chunk_max));
                    int best_lane = __ffs(max_mask) - 1;
                    int best_pmatch = __shfl_sync(0xFFFFFFFFU, pmatch, best_lane);
                    fscoreMax = chunk_max;
                    xfilmax2 = xfil2 - off - best_lane;
                    maxIdent = nIdent + best_pmatch;
                }

                if (sep_mask || neg_mask) break;

                fscore = __shfl_sync(0xFFFFFFFFU, score_here, 31);
                nIdent += __shfl_sync(0xFFFFFFFFU, pmatch, 31);
                off += 32;
            }
        }

        /* F12 + F13: Lane-0-only guarded FP64 divide + broadcast accepted + lossless covered_x */
        uint32_t length = (uint32_t)(xfilmax - xfilmax2 + 1);
        if (length > min_len) {
            int accepted_i = 0;
            if (lane == 0 && fscoreMax > 0) {
                float sim = (float)((double)fscoreMax * 100.0 /
                                    (double)((uint64_t)length * (uint64_t)match_score));
                if (sim > sim_th) {
                    accepted_i = 1;
                    uint32_t slot = atomicAdd(frag_count, 1U);
                    if (slot < max_frags) {
                        uint32_t yStart = (uint32_t)((int32_t)posY + (xfilmax2 - (int32_t)posX));
                        out_frags[slot].xStart = (uint32_t)xfilmax2;
                        out_frags[slot].yStart = yStart;
                        out_frags[slot].length = length;
                        out_frags[slot].ident  = (uint32_t)maxIdent;
                        out_frags[slot].score  = (uint32_t)fscoreMax;
                        out_frags[slot].flags  = is_rev0 & 1U;
                    }
                }
            }
            accepted_i = __shfl_sync(0xFFFFFFFFU, accepted_i, 0);
            if (accepted_i && seed_fscore >= 0) {
                int32_t new_cov = is_spaced ? (xfilmax - (wl - 1)) : xfilmax;
                if (new_cov > covered_x) {
                    covered_x = new_cov;
                    covered_min_x = is_spaced ? ((int32_t)posX + wl) : 0;
                }
            }
        } else if (seed_fscore >= 0) {
            int32_t new_cov = is_spaced ? (xfilmax - (wl - 1)) : xfilmax;
            int32_t new_min = is_spaced ? ((int32_t)posX + wl) : 0;
            if (new_cov > covered_x && (!is_spaced || new_cov >= new_min)) {
                covered_x = new_cov;
                covered_min_x = new_min;
            }
        }
    }
}

uint32_t gpu_extend_tile_two_tier(
    const uint64_t *d_filtered_hits,
    uint32_t n_hits,
    HitKeyLayout key_layout,
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuSeqView &seqY_rc,
    int wl,
    bool is_spaced,
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
    bool &out_overflowed,
    GeckoCudaProfile *prof
) {
    out_overflowed = false;
    CUDA_CHECK_VOID(cudaMemset(d_long_count, 0, sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMemset(d_frag_count, 0, sizeof(uint32_t)));
    if (n_hits == 0U) return 0U;

    {
        CudaEventScope ev(prof, prof ? &prof->extend_short_ms : nullptr);
        extend_short_kernel<<<(n_hits + 255U) / 256U, 256>>>(
            d_filtered_hits, n_hits, key_layout, seqX, seqY, seqY_rc,
            wl, is_spaced, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
            d_long_seeds, d_long_count, d_out_frags, d_frag_count, max_frags);
        CUDA_CHECK_VOID(cudaGetLastError());
    }

    uint32_t h_nlong = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nlong, d_long_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (prof && prof->enabled) {
        prof->long_seeds += h_nlong;
    }
    if (h_nlong > 0U) {
        const uint32_t warps_per_block = 8U;
        uint32_t grid = (h_nlong + warps_per_block - 1U) / warps_per_block;
        CudaEventScope ev(prof, prof ? &prof->extend_long_ms : nullptr);
        extend_long_warp_kernel<<<grid, warps_per_block * 32U>>>(
            d_filtered_hits, d_long_seeds, h_nlong, key_layout, seqX, seqY, seqY_rc,
            wl, is_spaced, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
            d_out_frags, d_frag_count, max_frags);
        CUDA_CHECK_VOID(cudaGetLastError());
    }

    uint32_t h_nfrags = 0U;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nfrags, d_frag_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nfrags > max_frags) {
        out_overflowed = true;
        return 0U;
    }
    return h_nfrags;
}
