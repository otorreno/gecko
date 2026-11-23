/*
 * cuda_index.cu
 *
 * Stage 1: Word-level 2-bit GPU reverse-complement generation, 2-bit canonical
 * k-mer extraction, CUB DeviceRadixSort::SortPairs, exact per-strand freq_limit
 * compaction + 1-bit dropped-position bitmap, and 2-pass suffix-min bucket table.
 */

#include <utility>
#include <vector>
#include <cub/cub.cuh>

#include "cuda_common.cuh"

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
        /* Fast path: every base is valid */
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
            /* Entire slice is invalid (N gap or soft-masked repeat) */
            valid = 0;
            i += avail;
            continue;
        }
        if (inv_masked == 0U) {
            /* Entire slice is valid */
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

static __global__ void extract_kmers_kernel(
    GpuSeqView seq,
    int wl,
    uint32_t num_chunks,
    const uint32_t * __restrict__ chunk_offsets,
    uint64_t * __restrict__ out_keys,
    uint32_t * __restrict__ out_vals
) {
    uint32_t cid = blockIdx.x * blockDim.x + threadIdx.x;
    if (cid >= num_chunks) return;

    uint64_t c_start = (uint64_t)cid * GECKO_EXTRACT_CHUNK_BP;
    uint64_t c_end = c_start + GECKO_EXTRACT_CHUNK_BP;
    if (c_end > seq.len) c_end = seq.len;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rc_shift = 2 * (wl - 1);
    uint64_t scan_from = (c_start >= (uint64_t)(wl - 1)) ? (c_start - (uint64_t)(wl - 1)) : 0ULL;

    uint64_t kf = 0, kr = 0;
    int valid = 0;
    uint32_t dst = chunk_offsets[cid];

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
                    uint32_t strand = (kf <= kr) ? 0U : 1U;
                    uint64_t key = strand ? kr : kf;
                    uint32_t start_pos = (uint32_t)((i + k) + 1ULL - (uint64_t)wl);
                    out_keys[dst] = fmix64(key);
                    out_vals[dst] = (start_pos << 1) | strand;
                    dst++;
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

static __global__ void mark_freq_filter_kernel(
    const uint64_t * __restrict__ keys,
    const uint32_t * __restrict__ vals,
    uint32_t n,
    uint32_t freq_limit,
    uint8_t * __restrict__ keep_flags
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    uint64_t k = keys[i];
    if (i > 0 && keys[i - 1] == k) return;

    if (i + 1 >= n || keys[i + 1] != k) {
        keep_flags[i] = 1;
        return;
    }

    uint32_t step = 1;
    uint32_t hi = i + 1;
    while (hi < n && keys[hi] == k) {
        step <<= 1;
        uint32_t next = i + step;
        hi = (next < n && next > hi) ? next : n;
    }
    uint32_t lo = i + (step >> 1);
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1);
        if (keys[mid] == k) lo = mid + 1;
        else hi = mid;
    }
    uint32_t run_end = lo;
    uint32_t run_len = run_end - i;

    if (run_len <= freq_limit) {
        for (uint32_t m = i; m < run_end; m++) keep_flags[m] = 1;
        return;
    }

    uint32_t c0 = 0, c1 = 0;
    for (uint32_t m = i; m < run_end; m++) {
        if (vals[m] & 1U) c1++;
        else c0++;
        if (c0 > freq_limit && c1 > freq_limit) return;
    }
    uint8_t k0 = (c0 <= freq_limit) ? 1 : 0;
    uint8_t k1 = (c1 <= freq_limit) ? 1 : 0;
    for (uint32_t m = i; m < run_end; m++) {
        keep_flags[m] = (vals[m] & 1U) ? k1 : k0;
    }
}

static __global__ void mark_dropped_bitmap_kernel(
    const uint32_t * __restrict__ vals,
    const uint8_t * __restrict__ keep_flags,
    uint32_t n_valid,
    uint32_t * __restrict__ dropped_bitmap
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_valid) return;
    if (keep_flags[i] == 0) {
        uint32_t pos = vals[i] >> 1;
        atomicOr(&dropped_bitmap[pos >> 5], 1U << (pos & 31U));
    }
}

static __global__ void mark_bucket_starts_and_tags_kernel(
    const uint64_t * __restrict__ keys,
    uint32_t n,
    int shift_bits,
    int tag_shift,
    uint32_t * __restrict__ bucket_off,
    uint32_t * __restrict__ out_tags
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint64_t k = keys[i];
    out_tags[i] = (uint32_t)(k >> tag_shift);
    uint32_t b = (uint32_t)(k >> shift_bits);
    if (i == 0 || b != (uint32_t)(keys[i - 1] >> shift_bits)) {
        bucket_off[b] = i;
    }
}

static __global__ void suffix_min_block_kernel(
    uint32_t * __restrict__ arr,
    uint32_t n,
    uint32_t * __restrict__ block_mins
) {
    const uint32_t BS = 256;
    __shared__ uint32_t s[BS];
    uint32_t idx = blockIdx.x * BS + threadIdx.x;
    s[threadIdx.x] = (idx < n) ? arr[idx] : 0xFFFFFFFFU;
    __syncthreads();

    for (uint32_t off = 1; off < BS; off <<= 1) {
        uint32_t other = (threadIdx.x + off < BS) ? s[threadIdx.x + off] : 0xFFFFFFFFU;
        __syncthreads();
        if (other < s[threadIdx.x]) s[threadIdx.x] = other;
        __syncthreads();
    }
    if (idx < n) arr[idx] = s[threadIdx.x];
    if (threadIdx.x == 0) block_mins[blockIdx.x] = s[0];
}

static __global__ void suffix_min_apply_kernel(
    uint32_t * __restrict__ arr,
    uint32_t n,
    const uint32_t * __restrict__ block_suf_mins,
    uint32_t num_blocks
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n || blockIdx.x + 1 >= num_blocks) return;
    uint32_t next_min = block_suf_mins[blockIdx.x + 1];
    if (next_min < arr[idx]) arr[idx] = next_min;
}

void gpu_build_revcomp_2bit(
    const uint64_t *d_src_bits,
    uint64_t *d_dst_bits,
    const uint32_t *d_src_special_bm,
    uint32_t *d_dst_special_bm,
    uint64_t len
) {
    if (len == 0) return;
    uint32_t n_words = (uint32_t)((len + 31ULL) / 32ULL);
    uint32_t grid = (n_words + 255U) / 256U;
    revcomp_2bit_kernel<<<grid, 256>>>(d_src_bits, d_dst_bits, len, n_words);
    CUDA_CHECK_VOID(cudaGetLastError());
    if (d_src_special_bm && d_dst_special_bm) {
        rev_bitmap_kernel<<<grid, 256>>>(d_src_special_bm, d_dst_special_bm, len, n_words);
        CUDA_CHECK_VOID(cudaGetLastError());
    }
}

GpuRefIndex gpu_build_ref_index(
    const GpuSeqView &seqY,
    int wl,
    uint32_t freq_limit
) {
    GpuRefIndex idx{};
    uint64_t lenY = seqY.len;
    uint32_t num_chunks = (uint32_t)((lenY + GECKO_EXTRACT_CHUNK_BP - 1) / GECKO_EXTRACT_CHUNK_BP);
    uint32_t *d_chunk_counts = nullptr, *d_chunk_offsets = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_chunk_counts, (size_t)(num_chunks > 0 ? num_chunks : 1) * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_chunk_offsets, (size_t)(num_chunks + 1) * sizeof(uint32_t)));

    void *d_cub_temp = nullptr;
    size_t cub_temp_cap = 0;
    auto ensure_idx_temp = [&](size_t needed) {
        if (needed <= cub_temp_cap) return;
        if (d_cub_temp) {
            CUDA_CHECK_VOID(cudaFree(d_cub_temp));
            d_cub_temp = nullptr;
            cub_temp_cap = 0;
        }
        CUDA_CHECK_VOID(cudaMalloc(&d_cub_temp, needed));
        cub_temp_cap = needed;
    };

    if (num_chunks > 0) {
        count_valid_kmers_kernel<<<(num_chunks + 255) / 256, 256>>>(
            seqY, wl, num_chunks, d_chunk_counts);
        CUDA_CHECK_VOID(cudaGetLastError());

        size_t temp_bytes = 0;
        cub::DeviceScan::ExclusiveSum(nullptr, temp_bytes, d_chunk_counts, d_chunk_offsets, num_chunks);
        ensure_idx_temp(temp_bytes);
        cub::DeviceScan::ExclusiveSum(d_cub_temp, temp_bytes, d_chunk_counts, d_chunk_offsets, num_chunks);

        uint32_t last_off = 0, last_cnt = 0;
        CUDA_CHECK_VOID(cudaMemcpy(&last_off, d_chunk_offsets + (num_chunks - 1), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        CUDA_CHECK_VOID(cudaMemcpy(&last_cnt, d_chunk_counts + (num_chunks - 1), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        idx.n_valid = last_off + last_cnt;
    }
    CUDA_CHECK_VOID(cudaFree(d_chunk_counts));

    size_t n_alloc = (idx.n_valid > 0) ? (size_t)idx.n_valid : 1;
    uint64_t *d_keys_a = nullptr, *d_keys_b = nullptr;
    uint32_t *d_vals_a = nullptr, *d_vals_b = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_keys_a, n_alloc * sizeof(uint64_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_vals_a, n_alloc * sizeof(uint32_t)));

    if (idx.n_valid > 0) {
        extract_kmers_kernel<<<(num_chunks + 255) / 256, 256>>>(
            seqY, wl, num_chunks, d_chunk_offsets, d_keys_a, d_vals_a);
        CUDA_CHECK_VOID(cudaGetLastError());
    }
    CUDA_CHECK_VOID(cudaFree(d_chunk_offsets));

    CUDA_CHECK_VOID(cudaMalloc(&d_keys_b, n_alloc * sizeof(uint64_t)));
    CUDA_CHECK_VOID(cudaMalloc(&d_vals_b, n_alloc * sizeof(uint32_t)));
    cub::DoubleBuffer<uint64_t> db_keys(d_keys_a, d_keys_b);
    cub::DoubleBuffer<uint32_t> db_vals(d_vals_a, d_vals_b);

    if (idx.n_valid > 1) {
        size_t temp_bytes = 0;
        cub::DeviceRadixSort::SortPairs(nullptr, temp_bytes, db_keys, db_vals, (int)idx.n_valid);
        ensure_idx_temp(temp_bytes);
        cub::DeviceRadixSort::SortPairs(d_cub_temp, temp_bytes, db_keys, db_vals, (int)idx.n_valid);
    }

    uint64_t *d_sorted_keys = db_keys.Current();
    uint32_t *d_sorted_vals = db_vals.Current();
    uint64_t *d_alt_keys = db_keys.Alternate();
    uint32_t *d_alt_vals = db_vals.Alternate();

    idx.n_kept = idx.n_valid;
    if (freq_limit > 0 && idx.n_valid > 0) {
        uint8_t *d_keep = nullptr;
        uint32_t *d_num_out = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_keep, n_alloc * sizeof(uint8_t)));
        CUDA_CHECK_VOID(cudaMemset(d_keep, 0, n_alloc * sizeof(uint8_t)));
        CUDA_CHECK_VOID(cudaMalloc(&d_num_out, sizeof(uint32_t)));

        mark_freq_filter_kernel<<<(idx.n_valid + 255) / 256, 256>>>(
            d_sorted_keys, d_sorted_vals, idx.n_valid, freq_limit, d_keep);

        size_t temp_bytes = 0, temp_bytes2 = 0;
        cub::DeviceSelect::Flagged(nullptr, temp_bytes, d_sorted_keys, d_keep, d_alt_keys, d_num_out, (int)idx.n_valid);
        cub::DeviceSelect::Flagged(nullptr, temp_bytes2, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
        size_t need_sel = (temp_bytes2 > temp_bytes) ? temp_bytes2 : temp_bytes;
        ensure_idx_temp(need_sel);
        cub::DeviceSelect::Flagged(d_cub_temp, temp_bytes, d_sorted_keys, d_keep, d_alt_keys, d_num_out, (int)idx.n_valid);
        cub::DeviceSelect::Flagged(d_cub_temp, temp_bytes2, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
        CUDA_CHECK_VOID(cudaMemcpy(&idx.n_kept, d_num_out, sizeof(uint32_t), cudaMemcpyDeviceToHost));

        if (idx.n_kept < idx.n_valid && lenY > 0) {
            size_t bitmap_words = (size_t)((lenY + 31ULL) / 32ULL);
            CUDA_CHECK_VOID(cudaMalloc(&idx.d_dropped_bitmap, bitmap_words * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMemset(idx.d_dropped_bitmap, 0, bitmap_words * sizeof(uint32_t)));
            mark_dropped_bitmap_kernel<<<(idx.n_valid + 255) / 256, 256>>>(
                d_sorted_vals, d_keep, idx.n_valid, idx.d_dropped_bitmap);
            CUDA_CHECK_VOID(cudaGetLastError());
        }

        CUDA_CHECK_VOID(cudaFree(d_keep));
        CUDA_CHECK_VOID(cudaFree(d_num_out));
        std::swap(d_sorted_keys, d_alt_keys);
        std::swap(d_sorted_vals, d_alt_vals);
    }
    if (d_cub_temp) CUDA_CHECK_VOID(cudaFree(d_cub_temp));
    CUDA_CHECK_VOID(cudaFree(d_alt_keys));
    CUDA_CHECK_VOID(cudaFree(d_alt_vals));

    int bucket_bits = 14;
    while (bucket_bits < 25 && (1ULL << bucket_bits) * 4ULL < idx.n_kept) bucket_bits++;
    uint32_t num_buckets = 1U << bucket_bits;
    int shift_bits = 64 - bucket_bits;
    int tag_shift = shift_bits - 32; /* always >= 7 since bucket_bits <= 25 */

    uint32_t *d_bucket_off = nullptr;
    CUDA_CHECK_VOID(cudaMalloc(&d_bucket_off, (size_t)(num_buckets + 1) * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMemset(d_bucket_off, 0xFF, (size_t)num_buckets * sizeof(uint32_t)));
    CUDA_CHECK_VOID(cudaMemcpy(d_bucket_off + num_buckets, &idx.n_kept, sizeof(uint32_t), cudaMemcpyHostToDevice));

    uint32_t *d_tags = nullptr;
    size_t tag_alloc = (idx.n_kept > 0) ? (size_t)idx.n_kept : 1;
    CUDA_CHECK_VOID(cudaMalloc(&d_tags, tag_alloc * sizeof(uint32_t)));

    if (idx.n_kept > 0) {
        mark_bucket_starts_and_tags_kernel<<<(idx.n_kept + 255) / 256, 256>>>(
            d_sorted_keys, idx.n_kept, shift_bits, tag_shift, d_bucket_off, d_tags);
        CUDA_CHECK_VOID(cudaGetLastError());

        uint32_t n_b_blocks = (num_buckets + 1 + 255) / 256;
        uint32_t *d_block_mins = nullptr;
        CUDA_CHECK_VOID(cudaMalloc(&d_block_mins, (size_t)n_b_blocks * sizeof(uint32_t)));
        suffix_min_block_kernel<<<n_b_blocks, 256>>>(d_bucket_off, num_buckets + 1, d_block_mins);

        std::vector<uint32_t> h_bmins(n_b_blocks);
        CUDA_CHECK_VOID(cudaMemcpy(h_bmins.data(), d_block_mins, (size_t)n_b_blocks * sizeof(uint32_t), cudaMemcpyDeviceToHost));
        for (int64_t bi = (int64_t)n_b_blocks - 2; bi >= 0; bi--) {
            if (h_bmins[bi + 1] < h_bmins[bi]) h_bmins[bi] = h_bmins[bi + 1];
        }
        CUDA_CHECK_VOID(cudaMemcpy(d_block_mins, h_bmins.data(), (size_t)n_b_blocks * sizeof(uint32_t), cudaMemcpyHostToDevice));
        suffix_min_apply_kernel<<<n_b_blocks, 256>>>(d_bucket_off, num_buckets + 1, d_block_mins, n_b_blocks);
        CUDA_CHECK_VOID(cudaFree(d_block_mins));
    } else {
        CUDA_CHECK_VOID(cudaMemset(d_bucket_off, 0, (size_t)(num_buckets + 1) * sizeof(uint32_t)));
    }

    /* Free 64-bit keys immediately — steady-state VRAM stores only 32-bit tags (~50% smaller) */
    CUDA_CHECK_VOID(cudaFree(d_sorted_keys));

    idx.d_tags = d_tags;
    idx.d_vals = d_sorted_vals;
    idx.d_bucket_off = d_bucket_off;
    idx.bucket_bits = bucket_bits;
    idx.shift_bits = shift_bits;
    idx.tag_shift = tag_shift;
    return idx;
}

void gpu_free_ref_index(GpuRefIndex &idx) {
    if (idx.d_tags) CUDA_CHECK_VOID(cudaFree(idx.d_tags));
    if (idx.d_vals) CUDA_CHECK_VOID(cudaFree(idx.d_vals));
    if (idx.d_bucket_off) CUDA_CHECK_VOID(cudaFree(idx.d_bucket_off));
    if (idx.d_dropped_bitmap) CUDA_CHECK_VOID(cudaFree(idx.d_dropped_bitmap));
    idx = {};
}
