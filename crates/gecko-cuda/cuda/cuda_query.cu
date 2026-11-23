/*
 * cuda_query.cu
 *
 * Stage 2: Tiled 2-bit GPU Query + Lock-Free Left/Right-Maximal Seed Filtering
 * (with 1-bit freq_limit dropped-position bitmap) + CUB RadixSort & Diagonal
 * Proximity Compaction.
 */

#include <cub/cub.cuh>
#include "cuda_common.cuh"

static __global__ void query_left_maximal_kernel(
    GpuSeqView seqX,
    GpuSeqView seqY,
    const uint32_t * __restrict__ bucket_off,
    const uint32_t * __restrict__ idx_tags,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int shift_bits,
    int tag_shift,
    int wl,
    bool compare_reverse,
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits
) {
    const uint64_t STRIDE = 16;
    uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t t_start = tile_start + tid * STRIDE;
    if (t_start >= tile_end) return;
    uint64_t t_end = t_start + STRIDE;
    if (t_end > tile_end) t_end = tile_end;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rc_shift = 2 * (wl - 1);
    int rev_rshift = 64 - 2 * wl;
    uint64_t scan_from = (t_start >= (uint64_t)(wl - 1)) ? (t_start - (uint64_t)(wl - 1)) : 0ULL;

    uint64_t kf = 0, kr = 0;
    int valid = 0;

    uint64_t i = scan_from;
    while (i < t_end) {
        uint32_t w_idx = (uint32_t)(i >> 5);
        uint32_t lane = (uint32_t)(i & 31ULL);
        uint32_t avail = (uint32_t)((32U - lane < t_end - i) ? (32U - lane) : (t_end - i));
        uint32_t inv_word = seqX.invalid_bm ? (__ldg(&seqX.invalid_bm[w_idx]) >> lane) : 0U;
        uint32_t avail_mask = (avail == 32U) ? 0xFFFFFFFFU : ((1U << avail) - 1U);
        if ((inv_word & avail_mask) == avail_mask) {
            valid = 0;
            kf = 0;
            kr = 0;
            i += avail;
            continue;
        }
        uint64_t bword = __ldg(&seqX.bits[w_idx]) >> (lane << 1);
        for (uint32_t k = 0; k < avail; k++) {
            uint64_t b = bword & 3ULL;
            bword >>= 2;
            if (((inv_word >> k) & 1U) == 0U) {
                kf = ((kf << 2) | b) & mask;
                kr = (kr >> 2) | ((b ^ 3ULL) << rc_shift);
                uint64_t cur_i = i + k;
                if (++valid >= wl && cur_i >= t_start) {
                    uint32_t posX = (uint32_t)(cur_i + 1ULL - (uint64_t)wl);
                    uint32_t x_strand = (kf <= kr) ? 0U : 1U;
                    uint64_t canon = x_strand ? kr : kf;
                    bool is_pal = (kf == kr);

                    uint64_t h = fmix64(canon);
                    uint32_t bkt = (uint32_t)(h >> shift_bits);
                    uint32_t st = __ldg(&bucket_off[bkt]);
                    uint32_t en = __ldg(&bucket_off[bkt + 1]);
                    if (st == en) continue;

                    uint32_t tag = (uint32_t)(h >> tag_shift);
                    uint32_t j = st;
                    if (en - st <= 16U) {
                        while (j < en && __ldg(&idx_tags[j]) < tag) j++;
                    } else {
                        uint32_t lo = st, hi = en;
                        while (lo < hi) {
                            uint32_t mid = lo + ((hi - lo) >> 1);
                            if (__ldg(&idx_tags[mid]) < tag) lo = mid + 1;
                            else hi = mid;
                        }
                        j = lo;
                    }
                    if (j >= en || __ldg(&idx_tags[j]) != tag) continue;

                    uint64_t x_seed_lsb = load_32_2bit(seqX.bits, posX) & mask;
                    bool prev_x_valid = (posX > 0U) && !is_bit_set(seqX.invalid_bm, posX - 1U);
                    uint32_t prev_x_2bit = prev_x_valid ? base_2bit_dev(seqX.bits, posX - 1U) : 4U;
                    bool next_x_valid = ((uint64_t)posX + (uint64_t)wl < seqX.len) &&
                                        !is_bit_set(seqX.invalid_bm, posX + (uint32_t)wl);
                    uint32_t next_x_2bit = next_x_valid ? base_2bit_dev(seqX.bits, posX + (uint32_t)wl) : 4U;

                    while (j < en && __ldg(&idx_tags[j]) == tag) {
                        uint32_t ps = __ldg(&idx_vals[j]);
                        uint32_t posY_orig = ps >> 1;
                        uint32_t y_strand = ps & 1U;
                        uint64_t y_word = load_32_2bit(seqY.bits, posY_orig);

                        if (is_pal || (x_strand == y_strand)) {
                            if (((y_word ^ x_seed_lsb) & mask) == 0ULL) {
                                bool left_ext = prev_x_valid && (posY_orig > 0U) &&
                                                !is_bit_set(seqY.invalid_bm, posY_orig - 1U) &&
                                                (base_2bit_dev(seqY.bits, posY_orig - 1U) == prev_x_2bit) &&
                                                !is_bit_set(dropped_bitmap, posY_orig - 1U);
                                bool right_ext = next_x_valid && ((uint64_t)posY_orig + (uint64_t)wl < seqY.len) &&
                                                 !is_bit_set(seqY.invalid_bm, posY_orig + (uint32_t)wl) &&
                                                 (base_2bit_dev(seqY.bits, posY_orig + (uint32_t)wl) == next_x_2bit) &&
                                                 !is_bit_set(dropped_bitmap, posY_orig + 1U);
                                if (!left_ext || !right_ext) {
                                    if (*hit_count <= max_hits) {
                                        uint32_t slot = atomicAdd(hit_count, 1U);
                                        if (slot < max_hits) {
                                            out_hits[slot] = pack_hit_key(0U, posX, posY_orig);
                                        }
                                    }
                                }
                            }
                        }
                        if (compare_reverse && (is_pal || (x_strand != y_strand))) {
                            uint64_t y_rc_lsb = (~reverse_32_2bit_dev(y_word)) >> rev_rshift;
                            if (((y_rc_lsb ^ x_seed_lsb) & mask) == 0ULL) {
                                /* Check mirrored neighbors directly in seqY (avoids storing/reading seqY_rc.invalid_bm):
                                 * Y_rc[posY_rc - 1]  == comp(Y[posY_orig + wl])
                                 * Y_rc[posY_rc + wl] == comp(Y[posY_orig - 1]) */
                                uint32_t posY_rc = (uint32_t)(seqY.len - (uint64_t)posY_orig - (uint64_t)wl);
                                bool left_ext = prev_x_valid && (posY_rc > 0U) &&
                                                !is_bit_set(seqY.invalid_bm, posY_orig + (uint32_t)wl) &&
                                                ((base_2bit_dev(seqY.bits, posY_orig + (uint32_t)wl) ^ 3U) == prev_x_2bit) &&
                                                !is_bit_set(dropped_bitmap, posY_orig + 1U);
                                bool right_ext = next_x_valid && (posY_orig > 0U) &&
                                                 !is_bit_set(seqY.invalid_bm, posY_orig - 1U) &&
                                                 ((base_2bit_dev(seqY.bits, posY_orig - 1U) ^ 3U) == next_x_2bit) &&
                                                 !is_bit_set(dropped_bitmap, posY_orig - 1U);
                                if (!left_ext || !right_ext) {
                                    if (*hit_count <= max_hits) {
                                        uint32_t slot = atomicAdd(hit_count, 1U);
                                        if (slot < max_hits) {
                                            out_hits[slot] = pack_hit_key(1U, posX, posY_rc);
                                        }
                                    }
                                }
                            }
                        }
                        j++;
                    }
                }
            } else {
                valid = 0;
                kf = 0;
                kr = 0;
            }
        }
        i += avail;
    }
}

static __global__ void filter_sorted_hits_kernel(
    const uint64_t * __restrict__ sorted_hits,
    uint32_t n,
    uint32_t min_gap,
    uint8_t * __restrict__ keep_flags
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (i == 0) {
        keep_flags[0] = 1;
        return;
    }
    uint64_t cur = sorted_hits[i];
    uint64_t prev = sorted_hits[i - 1];
    if ((cur >> 31) != (prev >> 31)) {
        keep_flags[i] = 1;
    } else {
        uint32_t cur_x = (uint32_t)(cur & 0x7FFFFFFFU);
        uint32_t prev_x = (uint32_t)(prev & 0x7FFFFFFFU);
        keep_flags[i] = (cur_x > prev_x + min_gap) ? 1 : 0;
    }
}

static void ensure_cub_temp(void *&d_temp, size_t &cap, size_t needed) {
    if (needed <= cap) return;
    if (d_temp) {
        CUDA_CHECK_VOID(cudaFree(d_temp));
        d_temp = nullptr;
        cap = 0;
    }
    /* Round up to 4 MiB boundary to avoid repeated reallocations as hit counts grow */
    size_t alloc_bytes = (needed + (4ULL << 20) - 1ULL) & ~((4ULL << 20) - 1ULL);
    CUDA_CHECK_VOID(cudaMalloc(&d_temp, alloc_bytes));
    cap = alloc_bytes;
}

uint32_t gpu_query_tile_left_maximal(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuRefIndex &idx,
    int wl,
    bool compare_reverse,
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint8_t *d_hit_flags,
    uint32_t *d_hit_count,
    uint32_t *d_filt_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed
) {
    out_overflowed = false;
    out_raw_leftmax = 0;
    CUDA_CHECK_VOID(cudaMemset(d_hit_count, 0, sizeof(uint32_t)));

    uint64_t n_threads = ((tile_end - tile_start) + 15) / 16;
    uint32_t grid = (uint32_t)((n_threads + 255) / 256);
    query_left_maximal_kernel<<<grid, 256>>>(
        seqX, seqY,
        idx.d_bucket_off, idx.d_tags, idx.d_vals, idx.d_dropped_bitmap, idx.shift_bits, idx.tag_shift, wl,
        compare_reverse,
        tile_start, tile_end, d_hits, d_hit_count, max_hits);
    CUDA_CHECK_VOID(cudaGetLastError());

    uint32_t h_nhits = 0;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nhits, d_hit_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nhits > max_hits) {
        out_overflowed = true;
        return 0;
    }
    out_raw_leftmax = h_nhits;
    if (h_nhits <= 1) return h_nhits;

    int end_bit = compare_reverse ? 64 : 63;
    size_t sort_bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, sort_bytes, d_hits, d_hits_sorted, (int)h_nhits, 0, end_bit);
    size_t sel_bytes = 0;
    cub::DeviceSelect::Flagged(nullptr, sel_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);
    size_t need_bytes = (sort_bytes > sel_bytes) ? sort_bytes : sel_bytes;
    ensure_cub_temp(d_cub_temp, cub_temp_cap, need_bytes);

    cub::DeviceRadixSort::SortKeys(d_cub_temp, sort_bytes, d_hits, d_hits_sorted, (int)h_nhits, 0, end_bit);

    filter_sorted_hits_kernel<<<(h_nhits + 255) / 256, 256>>>(
        d_hits_sorted, h_nhits, (uint32_t)(wl / 2), d_hit_flags);

    cub::DeviceSelect::Flagged(d_cub_temp, sel_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);

    uint32_t h_nfilt = 0;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nfilt, d_filt_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    return h_nfilt;
}
