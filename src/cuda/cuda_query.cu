/*
 * cuda_query.cu
 *
 * Stage 2: Tiled GPU query scan with lock-free Left-Maximal Seed Filtering,
 * 64-bit CUB radix sort by (strand, diag, posX), and diagonal gap compaction.
 */

#include <cstdio>
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
    int shift_bits,
    int wl,
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t * __restrict__ out_hits,
    uint32_t * __restrict__ hit_count,
    uint32_t max_hits
) {
    const uint32_t STRIDE = 16;
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
        uint8_t b = base_to_2bit(seqX[i]);
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl && i >= t_start) {
                uint32_t posX = (uint32_t)(i + 1 - (uint64_t)wl);
                uint32_t x_strand = (kf <= kr) ? 0U : 1U;
                uint64_t key = x_strand ? kr : kf;
                uint32_t is_pal = (kf == kr) ? 1U : 0U;

                uint64_t h = fmix64(key);
                uint32_t bkt = (uint32_t)(h >> shift_bits);
                uint32_t st = bucket_off[bkt];
                uint32_t en = bucket_off[bkt + 1];
                if (st == en) continue;

                uint32_t j = st;
                if (en - st <= 16) {
                    while (j < en && idx_keys[j] < h) j++;
                } else {
                    uint32_t lo = st, hi = en;
                    while (lo < hi) {
                        uint32_t mid = lo + ((hi - lo) >> 1);
                        if (idx_keys[mid] < h) lo = mid + 1;
                        else hi = mid;
                    }
                    j = lo;
                }

                uint8_t prev_x = (posX > 0) ? seqX[posX - 1] : 0;
                bool prev_x_valid = (base_to_2bit(prev_x) < 4);
                uint8_t next_x = ((uint64_t)posX + (uint64_t)wl < lenX) ? seqX[(uint64_t)posX + (uint64_t)wl] : 0;
                bool next_x_valid = (base_to_2bit(next_x) < 4);

                while (j < en && idx_keys[j] == h) {
                    uint32_t ps = idx_vals[j];
                    uint32_t posY_orig = ps >> 1;
                    uint32_t y_strand = ps & 1U;

                    if (is_pal || (x_strand == y_strand)) {
                        uint32_t posY = posY_orig;
                        bool left_ext = prev_x_valid && (posY > 0) && (seqY[posY - 1] == prev_x);
                        bool right_ext = next_x_valid && ((uint64_t)posY + (uint64_t)wl < lenY) && (seqY[(uint64_t)posY + (uint64_t)wl] == next_x);
                        if (!left_ext || !right_ext) {
                            uint32_t slot = atomicAdd(hit_count, 1U);
                            if (slot < max_hits) {
                                out_hits[slot] = pack_hit_key(0U, posX, posY);
                            }
                        }
                    }
                    if (is_pal || (x_strand != y_strand)) {
                        uint32_t posY_rc = (uint32_t)(lenY - (uint64_t)posY_orig - (uint64_t)wl);
                        bool left_ext = prev_x_valid && (posY_rc > 0) && (seqY_rc[posY_rc - 1] == prev_x);
                        bool right_ext = next_x_valid && ((uint64_t)posY_rc + (uint64_t)wl < lenY) && (seqY_rc[(uint64_t)posY_rc + (uint64_t)wl] == next_x);
                        if (!left_ext || !right_ext) {
                            uint32_t slot = atomicAdd(hit_count, 1U);
                            if (slot < max_hits) {
                                out_hits[slot] = pack_hit_key(1U, posX, posY_rc);
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
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint8_t *d_hit_flags,
    uint32_t *d_hit_count,
    uint32_t *d_filt_count,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax
) {
    CUDA_CHECK(cudaMemset(d_hit_count, 0, sizeof(uint32_t)));

    uint64_t n_threads = ((tile_end - tile_start) + 15) / 16;
    uint32_t grid = (uint32_t)((n_threads + 255) / 256);
    query_left_maximal_kernel<<<grid, 256>>>(
        d_seqX, lenX, d_seqY, d_seqY_rc, lenY,
        idx.d_bucket_off, idx.d_keys, idx.d_vals, idx.shift_bits, wl,
        tile_start, tile_end, d_hits, d_hit_count, max_hits);
    CUDA_CHECK(cudaGetLastError());

    uint32_t h_nhits = 0;
    CUDA_CHECK(cudaMemcpy(&h_nhits, d_hit_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    if (h_nhits > max_hits) {
        fprintf(stderr, "Warning: tile hit buffer capped (%u > %u).\n", h_nhits, max_hits);
        h_nhits = max_hits;
    }
    out_raw_leftmax = h_nhits;
    if (h_nhits <= 1) return h_nhits;

    void *d_temp = nullptr;
    size_t temp_bytes = 0;
    cub::DeviceRadixSort::SortKeys(nullptr, temp_bytes, d_hits, d_hits_sorted, (int)h_nhits);
    CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceRadixSort::SortKeys(d_temp, temp_bytes, d_hits, d_hits_sorted, (int)h_nhits);
    CUDA_CHECK(cudaFree(d_temp));

    filter_sorted_hits_kernel<<<(h_nhits + 255) / 256, 256>>>(
        d_hits_sorted, h_nhits, (uint32_t)(wl / 2), d_hit_flags);

    d_temp = nullptr;
    temp_bytes = 0;
    cub::DeviceSelect::Flagged(nullptr, temp_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);
    CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
    cub::DeviceSelect::Flagged(d_temp, temp_bytes, d_hits_sorted, d_hit_flags, d_hits, d_filt_count, (int)h_nhits);
    CUDA_CHECK(cudaFree(d_temp));

    uint32_t h_nfilt = 0;
    CUDA_CHECK(cudaMemcpy(&h_nfilt, d_filt_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    return h_nfilt;
}
