/*
 * cuda_index.cu
 *
 * Stage 1: GPU reverse-complement generation, canonical k-mer extraction,
 * CUB DeviceRadixSort::SortPairs, exact per-strand freq_limit compaction,
 * and 2-pass reverse-suffix-min bucket offset table construction.
 */

#include <utility>
#include <vector>
#include <cub/cub.cuh>

#include "cuda_common.cuh"

static __global__ void revcomp_kernel(const uint8_t * __restrict__ src, uint8_t * __restrict__ dst, uint64_t len) {
    uint64_t i = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < len) {
        dst[i] = comp_base_dev(src[len - 1 - i]);
    }
}

static __global__ void count_valid_kmers_kernel(
    const uint8_t * __restrict__ seq,
    uint64_t len,
    int wl,
    uint32_t num_chunks,
    uint32_t * __restrict__ chunk_counts
) {
    uint32_t cid = blockIdx.x * blockDim.x + threadIdx.x;
    if (cid >= num_chunks) return;

    uint64_t c_start = (uint64_t)cid * GECKO_EXTRACT_CHUNK_BP;
    uint64_t c_end = c_start + GECKO_EXTRACT_CHUNK_BP;
    if (c_end > len) c_end = len;

    uint64_t scan_from = (c_start >= (uint64_t)(wl - 1)) ? (c_start - (uint64_t)(wl - 1)) : 0;
    int valid = 0;
    uint32_t cnt = 0;
    for (uint64_t i = scan_from; i < c_end; i++) {
        if (base_to_2bit(seq[i]) < 4) {
            if (++valid >= wl && i >= c_start) cnt++;
        } else {
            valid = 0;
        }
    }
    chunk_counts[cid] = cnt;
}

static __global__ void extract_kmers_kernel(
    const uint8_t * __restrict__ seq,
    uint64_t len,
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
    if (c_end > len) c_end = len;

    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rc_shift = 2 * (wl - 1);
    uint64_t scan_from = (c_start >= (uint64_t)(wl - 1)) ? (c_start - (uint64_t)(wl - 1)) : 0;

    uint64_t kf = 0, kr = 0;
    int valid = 0;
    uint32_t dst = chunk_offsets[cid];

    for (uint64_t i = scan_from; i < c_end; i++) {
        uint8_t b = base_to_2bit(seq[i]);
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl && i >= c_start) {
                uint32_t strand = (kf <= kr) ? 0U : 1U;
                uint64_t key = strand ? kr : kf;
                uint32_t start_pos = (uint32_t)(i + 1 - (uint64_t)wl);
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

static __global__ void mark_bucket_starts_kernel(
    const uint64_t * __restrict__ keys,
    uint32_t n,
    int shift_bits,
    uint32_t * __restrict__ bucket_off
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint32_t b = (uint32_t)(keys[i] >> shift_bits);
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

void gpu_build_revcomp(const uint8_t *d_src, uint8_t *d_dst, uint64_t len) {
    if (len == 0) return;
    uint32_t grid = (uint32_t)((len + 255) / 256);
    revcomp_kernel<<<grid, 256>>>(d_src, d_dst, len);
    CUDA_CHECK(cudaGetLastError());
}

GpuRefIndex gpu_build_ref_index(const uint8_t *d_seqY, uint64_t lenY, int wl, uint32_t freq_limit) {
    GpuRefIndex idx{};
    uint32_t num_chunks = (uint32_t)((lenY + GECKO_EXTRACT_CHUNK_BP - 1) / GECKO_EXTRACT_CHUNK_BP);
    uint32_t *d_chunk_counts = nullptr, *d_chunk_offsets = nullptr;
    CUDA_CHECK(cudaMalloc(&d_chunk_counts, (size_t)(num_chunks > 0 ? num_chunks : 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_chunk_offsets, (size_t)(num_chunks + 1) * sizeof(uint32_t)));

    if (num_chunks > 0) {
        count_valid_kmers_kernel<<<(num_chunks + 255) / 256, 256>>>(d_seqY, lenY, wl, num_chunks, d_chunk_counts);
        void *d_temp = nullptr;
        size_t temp_bytes = 0;
        cub::DeviceScan::ExclusiveSum(nullptr, temp_bytes, d_chunk_counts, d_chunk_offsets, num_chunks);
        CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
        cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes, d_chunk_counts, d_chunk_offsets, num_chunks);
        CUDA_CHECK(cudaFree(d_temp));

        uint32_t last_off = 0, last_cnt = 0;
        CUDA_CHECK(cudaMemcpy(&last_off, d_chunk_offsets + (num_chunks - 1), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&last_cnt, d_chunk_counts + (num_chunks - 1), sizeof(uint32_t), cudaMemcpyDeviceToHost));
        idx.n_valid = last_off + last_cnt;
    }
    CUDA_CHECK(cudaFree(d_chunk_counts));

    size_t n_alloc = (idx.n_valid > 0) ? (size_t)idx.n_valid : 1;
    uint64_t *d_keys_a = nullptr, *d_keys_b = nullptr;
    uint32_t *d_vals_a = nullptr, *d_vals_b = nullptr;
    CUDA_CHECK(cudaMalloc(&d_keys_a, n_alloc * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_vals_a, n_alloc * sizeof(uint32_t)));

    if (idx.n_valid > 0) {
        extract_kmers_kernel<<<(num_chunks + 255) / 256, 256>>>(
            d_seqY, lenY, wl, num_chunks, d_chunk_offsets, d_keys_a, d_vals_a);
        CUDA_CHECK(cudaGetLastError());
    }
    CUDA_CHECK(cudaFree(d_chunk_offsets));

    CUDA_CHECK(cudaMalloc(&d_keys_b, n_alloc * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_vals_b, n_alloc * sizeof(uint32_t)));
    cub::DoubleBuffer<uint64_t> db_keys(d_keys_a, d_keys_b);
    cub::DoubleBuffer<uint32_t> db_vals(d_vals_a, d_vals_b);

    if (idx.n_valid > 1) {
        void *d_temp = nullptr;
        size_t temp_bytes = 0;
        cub::DeviceRadixSort::SortPairs(nullptr, temp_bytes, db_keys, db_vals, (int)idx.n_valid);
        CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
        cub::DeviceRadixSort::SortPairs(d_temp, temp_bytes, db_keys, db_vals, (int)idx.n_valid);
        CUDA_CHECK(cudaFree(d_temp));
    }

    uint64_t *d_sorted_keys = db_keys.Current();
    uint32_t *d_sorted_vals = db_vals.Current();
    uint64_t *d_alt_keys = db_keys.Alternate();
    uint32_t *d_alt_vals = db_vals.Alternate();

    idx.n_kept = idx.n_valid;
    if (freq_limit > 0 && idx.n_valid > 0) {
        uint8_t *d_keep = nullptr;
        uint32_t *d_num_out = nullptr;
        CUDA_CHECK(cudaMalloc(&d_keep, n_alloc * sizeof(uint8_t)));
        CUDA_CHECK(cudaMemset(d_keep, 0, n_alloc * sizeof(uint8_t)));
        CUDA_CHECK(cudaMalloc(&d_num_out, sizeof(uint32_t)));

        mark_freq_filter_kernel<<<(idx.n_valid + 255) / 256, 256>>>(
            d_sorted_keys, d_sorted_vals, idx.n_valid, freq_limit, d_keep);

        void *d_temp = nullptr;
        size_t temp_bytes = 0, temp_bytes2 = 0;
        cub::DeviceSelect::Flagged(nullptr, temp_bytes, d_sorted_keys, d_keep, d_alt_keys, d_num_out, (int)idx.n_valid);
        cub::DeviceSelect::Flagged(nullptr, temp_bytes2, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
        if (temp_bytes2 > temp_bytes) temp_bytes = temp_bytes2;
        CUDA_CHECK(cudaMalloc(&d_temp, temp_bytes));
        cub::DeviceSelect::Flagged(d_temp, temp_bytes, d_sorted_keys, d_keep, d_alt_keys, d_num_out, (int)idx.n_valid);
        cub::DeviceSelect::Flagged(d_temp, temp_bytes, d_sorted_vals, d_keep, d_alt_vals, d_num_out, (int)idx.n_valid);
        CUDA_CHECK(cudaMemcpy(&idx.n_kept, d_num_out, sizeof(uint32_t), cudaMemcpyDeviceToHost));

        CUDA_CHECK(cudaFree(d_temp));
        CUDA_CHECK(cudaFree(d_keep));
        CUDA_CHECK(cudaFree(d_num_out));
        std::swap(d_sorted_keys, d_alt_keys);
        std::swap(d_sorted_vals, d_alt_vals);
    }
    CUDA_CHECK(cudaFree(d_alt_keys));
    CUDA_CHECK(cudaFree(d_alt_vals));

    int bucket_bits = 14;
    while (bucket_bits < 25 && (1ULL << bucket_bits) * 4ULL < idx.n_kept) bucket_bits++;
    uint32_t num_buckets = 1U << bucket_bits;
    int shift_bits = 64 - bucket_bits;

    uint32_t *d_bucket_off = nullptr;
    CUDA_CHECK(cudaMalloc(&d_bucket_off, (size_t)(num_buckets + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemset(d_bucket_off, 0xFF, (size_t)num_buckets * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_bucket_off + num_buckets, &idx.n_kept, sizeof(uint32_t), cudaMemcpyHostToDevice));

    if (idx.n_kept > 0) {
        mark_bucket_starts_kernel<<<(idx.n_kept + 255) / 256, 256>>>(d_sorted_keys, idx.n_kept, shift_bits, d_bucket_off);
        uint32_t n_b_blocks = (num_buckets + 1 + 255) / 256;
        uint32_t *d_block_mins = nullptr;
        CUDA_CHECK(cudaMalloc(&d_block_mins, (size_t)n_b_blocks * sizeof(uint32_t)));
        suffix_min_block_kernel<<<n_b_blocks, 256>>>(d_bucket_off, num_buckets + 1, d_block_mins);

        std::vector<uint32_t> h_bmins(n_b_blocks);
        CUDA_CHECK(cudaMemcpy(h_bmins.data(), d_block_mins, (size_t)n_b_blocks * sizeof(uint32_t), cudaMemcpyDeviceToHost));
        for (int64_t bi = (int64_t)n_b_blocks - 2; bi >= 0; bi--) {
            if (h_bmins[bi + 1] < h_bmins[bi]) h_bmins[bi] = h_bmins[bi + 1];
        }
        CUDA_CHECK(cudaMemcpy(d_block_mins, h_bmins.data(), (size_t)n_b_blocks * sizeof(uint32_t), cudaMemcpyHostToDevice));
        suffix_min_apply_kernel<<<n_b_blocks, 256>>>(d_bucket_off, num_buckets + 1, d_block_mins, n_b_blocks);
        CUDA_CHECK(cudaFree(d_block_mins));
    } else {
        CUDA_CHECK(cudaMemset(d_bucket_off, 0, (size_t)(num_buckets + 1) * sizeof(uint32_t)));
    }

    idx.d_keys = d_sorted_keys;
    idx.d_vals = d_sorted_vals;
    idx.d_bucket_off = d_bucket_off;
    idx.bucket_bits = bucket_bits;
    idx.shift_bits = shift_bits;
    return idx;
}

void gpu_free_ref_index(GpuRefIndex &idx) {
    if (idx.d_keys) CUDA_CHECK(cudaFree(idx.d_keys));
    if (idx.d_vals) CUDA_CHECK(cudaFree(idx.d_vals));
    if (idx.d_bucket_off) CUDA_CHECK(cudaFree(idx.d_bucket_off));
    idx = {};
}
