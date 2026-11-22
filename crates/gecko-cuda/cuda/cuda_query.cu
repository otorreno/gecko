/*
 * cuda_query.cu
 *
 * Stage 2: Tiled GPU Query + Lock-Free Left/Right-Maximal Seed Filtering +
 * CUB RadixSort & Diagonal Proximity Compaction.
 */

#include <cub/cub.cuh>
#include "cuda_common.cuh"

static __global__ void query_left_maximal_kernel(
    const uint8_t * __restrict__ seqX,
    uint64_t lenX,
    const uint8_t * __restrict__ seqY,
    const uint8_t * __restrict__ seqY_rc,
    uint64_t lenY,
    const uint32_t * __restrict__ bucket_off,
    const uint64_t * __restrict__ idx_keys,
    const uint32_t * __restrict__ idx_vals,
    const uint32_t * __restrict__ dropped_bitmap,
    int shift_bits,
    int wl,
    bool soft_mask,
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
    uint64_t scan_from = (t_start >= (uint64_t)(wl - 1)) ? (t_start - (uint64_t)(wl - 1)) : 0;

    uint64_t kf = 0, kr = 0;
    int valid = 0;

    for (uint64_t i = scan_from; i < t_end; i++) {
        uint8_t b = base_to_2bit(seqX[i], soft_mask);
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl && i >= t_start) {
                uint32_t posX = (uint32_t)(i + 1 - (uint64_t)wl);
                uint32_t x_strand = (kf <= kr) ? 0U : 1U;
                uint64_t canon = x_strand ? kr : kf;
                bool is_pal = (kf == kr);

                uint64_t h = fmix64(canon);
                uint32_t bkt = (uint32_t)(h >> shift_bits);
                uint32_t st = __ldg(&bucket_off[bkt]);
                uint32_t en = __ldg(&bucket_off[bkt + 1]);
                if (st == en) continue;

                uint32_t j = st;
                if (en - st <= 16U) {
                    while (j < en && __ldg(&idx_keys[j]) < h) j++;
                } else {
                    uint32_t lo = st, hi = en;
                    while (lo < hi) {
                        uint32_t mid = lo + ((hi - lo) >> 1);
                        if (__ldg(&idx_keys[mid]) < h) lo = mid + 1;
                        else hi = mid;
                    }
                    j = lo;
                }

                uint8_t prev_x = (posX > 0) ? norm_seed_byte(seqX[posX - 1], soft_mask) : 0;
                bool prev_x_valid = (posX > 0) && (base_to_2bit(prev_x, true) < 4);
                uint8_t next_x = ((uint64_t)posX + (uint64_t)wl < lenX)
                                     ? norm_seed_byte(seqX[(uint64_t)posX + (uint64_t)wl], soft_mask)
                                     : 0;
                bool next_x_valid = ((uint64_t)posX + (uint64_t)wl < lenX) && (base_to_2bit(next_x, true) < 4);

                while (j < en && __ldg(&idx_keys[j]) == h) {
                    uint32_t ps = __ldg(&idx_vals[j]);
                    uint32_t posY_orig = ps >> 1;
                    uint32_t y_strand = ps & 1U;

                    if (is_pal || (x_strand == y_strand)) {
                        bool left_ext = prev_x_valid && (posY_orig > 0) &&
                                        (norm_seed_byte(seqY[posY_orig - 1], soft_mask) == prev_x) &&
                                        !is_dropped_pos(dropped_bitmap, posY_orig - 1);
                        bool right_ext = next_x_valid && ((uint64_t)posY_orig + (uint64_t)wl < lenY) &&
                                         (norm_seed_byte(seqY[(uint64_t)posY_orig + (uint64_t)wl], soft_mask) == next_x) &&
                                         !is_dropped_pos(dropped_bitmap, posY_orig + 1);
                        if (!left_ext || !right_ext) {
                            if (*hit_count <= max_hits) {
                                uint32_t slot = atomicAdd(hit_count, 1U);
                                if (slot < max_hits) {
                                    out_hits[slot] = pack_hit_key(0U, posX, posY_orig);
                                }
                            }
                        }
                    }
                    if (compare_reverse && (is_pal || (x_strand != y_strand))) {
                        uint32_t posY_rc = (uint32_t)(lenY - (uint64_t)posY_orig - (uint64_t)wl);
                        bool left_ext = prev_x_valid && (posY_rc > 0) &&
                                        (norm_seed_byte(seqY_rc[posY_rc - 1], soft_mask) == prev_x) &&
                                        !is_dropped_pos(dropped_bitmap, posY_orig + 1);
                        bool right_ext = next_x_valid && ((uint64_t)posY_rc + (uint64_t)wl < lenY) &&
                                         (norm_seed_byte(seqY_rc[(uint64_t)posY_rc + (uint64_t)wl], soft_mask) == next_x) &&
                                         !is_dropped_pos(dropped_bitmap, posY_orig - 1);
                        if (!left_ext || !right_ext) {
                            if (*hit_count <= max_hits) {
                                uint32_t slot = atomicAdd(hit_count, 1U);
                                if (slot < max_hits) {
                                    out_hits[slot] = pack_hit_key(1U, posX, posY_rc);
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

uint32_t gpu_query_tile_left_maximal(
    const uint8_t *d_seqX,
    uint64_t lenX,
    const uint8_t *d_seqY,
    const uint8_t *d_seqY_rc,
    uint64_t lenY,
    const GpuRefIndex &idx,
    int wl,
    bool soft_mask,
    bool compare_reverse,
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint8_t *d_hit_flags,
    uint32_t *d_hit_count,
    uint32_t *d_filt_count,
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
        d_seqX, lenX, d_seqY, d_seqY_rc, lenY,
        idx.d_bucket_off, idx.d_keys, idx.d_vals, idx.d_dropped_bitmap, idx.shift_bits, wl,
        soft_mask, compare_reverse,
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

    void *d_temp = nullptr;
    size_t temp_bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, temp_bytes, d_hits, d_hits_sorted, (int)h_nhits);
    CUDA_CHECK_VOID(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceRadixSort::SortKeys(d_temp, temp_bytes, d_hits, d_hits_sorted, (int)h_nhits);
    CUDA_CHECK_VOID(cudaFree(d_temp));

    filter_sorted_hits_kernel<<<(h_nhits + 255) / 256, 256>>>(
        d_hits_sorted, h_nhits, (uint32_t)(wl / 2), d_hit_flags);

    d_temp = nullptr;
    temp_bytes = 0;
    cub::DeviceSelect::Flagged(nullptr, temp_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);
    CUDA_CHECK_VOID(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceSelect::Flagged(d_temp, temp_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);
    CUDA_CHECK_VOID(cudaFree(d_temp));

    uint32_t h_nfilt = 0;
    CUDA_CHECK_VOID(cudaMemcpy(&h_nfilt, d_filt_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    return h_nfilt;
}
