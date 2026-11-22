/*
 * cuda_bridge.cu
 *
 * Clean 2-phase C-ABI bridge between Rust (gecko-cuda) and CUDA kernels:
 *   1. gecko_cuda_probe()
 *   2. gecko_cuda_prepare_ref()  — uploads seqY, builds Y_rc & GPU index (can overlap with Rust loading seqX!)
 *   3. gecko_cuda_query_and_extend() — uploads seqX, runs tiled query + two-tier warp extension
 *   4. gecko_cuda_free_ref() / gecko_cuda_free_frags()
 */

#include <algorithm>
#include <cstring>
#include <vector>

#include "cuda_common.cuh"

struct GeckoCudaRefHandle {
    uint8_t *d_seqY;
    uint8_t *d_seqY_rc;
    uint64_t lenY;
    int wl;
    bool soft_mask;
    bool compare_reverse;
    GpuRefIndex idx;
};

struct GeckoCudaResult {
    GpuFrag *frags;
    uint64_t n_frags;
    uint64_t leftmax_hits;
    uint64_t extended_hits;
    uint32_t index_kept;
    uint32_t index_valid;
    int32_t bucket_bits;
};

extern "C" int gecko_cuda_probe(
    char *name_buf,
    size_t name_cap,
    int *sm_major,
    int *sm_minor,
    uint64_t *vram_bytes
) {
    int dev_count = 0;
    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) {
        return -1;
    }
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess) {
        return -1;
    }
    if (name_buf && name_cap > 0) {
        strncpy(name_buf, prop.name, name_cap - 1);
        name_buf[name_cap - 1] = '\0';
    }
    if (sm_major) *sm_major = prop.major;
    if (sm_minor) *sm_minor = prop.minor;
    if (vram_bytes) *vram_bytes = (uint64_t)prop.totalGlobalMem;
    return 0;
}

extern "C" GeckoCudaRefHandle *gecko_cuda_prepare_ref(
    const uint8_t *seqY,
    uint64_t lenY,
    int wl,
    uint32_t freq_limit,
    int soft_mask,
    int compare_reverse
) {
    int dev_count = 0;
    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) {
        return nullptr;
    }

    auto *h = new GeckoCudaRefHandle{};
    h->lenY = lenY;
    h->wl = wl;
    h->soft_mask = (soft_mask != 0);
    h->compare_reverse = (compare_reverse != 0);

    size_t allocY = (size_t)lenY + 64;
    if (cudaMalloc(&h->d_seqY, allocY) != cudaSuccess) {
        delete h;
        return nullptr;
    }
    if (lenY > 0) {
        CUDA_CHECK_VOID(cudaMemcpy(h->d_seqY, seqY, (size_t)lenY, cudaMemcpyHostToDevice));
    }
    CUDA_CHECK_VOID(cudaMemset(h->d_seqY + lenY, 0, 64));

    if (h->compare_reverse) {
        CUDA_CHECK_VOID(cudaMalloc(&h->d_seqY_rc, allocY));
        CUDA_CHECK_VOID(cudaMemset(h->d_seqY_rc + lenY, 0, 64));
        gpu_build_revcomp(h->d_seqY, h->d_seqY_rc, lenY);
    } else {
        CUDA_CHECK_VOID(cudaMalloc(&h->d_seqY_rc, 64));
        CUDA_CHECK_VOID(cudaMemset(h->d_seqY_rc, 0, 64));
    }

    h->idx = gpu_build_ref_index(h->d_seqY, lenY, wl, freq_limit, h->soft_mask);
    return h;
}

extern "C" int gecko_cuda_query_and_extend(
    const GeckoCudaRefHandle *h,
    const uint8_t *seqX,
    uint64_t lenX,
    uint32_t min_len,
    float sim_th,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    GeckoCudaResult *out
) {
    if (!h || !out) return -1;
    memset(out, 0, sizeof(*out));
    out->index_kept = h->idx.n_kept;
    out->index_valid = h->idx.n_valid;
    out->bucket_bits = h->idx.bucket_bits;

    uint8_t *d_seqX = nullptr;
    CUDA_CHECK(cudaMalloc(&d_seqX, (size_t)lenX + 64));
    if (lenX > 0) {
        CUDA_CHECK(cudaMemcpy(d_seqX, seqX, (size_t)lenX, cudaMemcpyHostToDevice));
    }
    CUDA_CHECK(cudaMemset(d_seqX + lenX, 0, 64));

    uint64_t *d_hits = nullptr, *d_hits_sorted = nullptr, *d_long_seeds = nullptr;
    uint8_t *d_hit_flags = nullptr;
    uint32_t *d_hit_count = nullptr, *d_filt_count = nullptr, *d_long_count = nullptr, *d_frag_count = nullptr;
    GpuFrag *d_out_frags = nullptr;

    CUDA_CHECK(cudaMalloc(&d_hits, GECKO_MAX_TILE_HITS * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_hits_sorted, GECKO_MAX_TILE_HITS * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_long_seeds, GECKO_MAX_TILE_HITS * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_hit_flags, GECKO_MAX_TILE_HITS * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&d_out_frags, GECKO_MAX_TILE_FRAGS * sizeof(GpuFrag)));
    CUDA_CHECK(cudaMalloc(&d_hit_count, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_filt_count, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_long_count, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_frag_count, sizeof(uint32_t)));

    std::vector<GpuFrag> all_frags;
    std::vector<GpuFrag> tile_frags;

    for (uint64_t tile_start = 0; tile_start < lenX; tile_start += GECKO_QUERY_TILE_BP) {
        uint64_t tile_end = std::min<uint64_t>(tile_start + GECKO_QUERY_TILE_BP, lenX);
        uint32_t raw_leftmax = 0;
        uint32_t n_filt = gpu_query_tile_left_maximal(
            d_seqX, lenX, h->d_seqY, h->d_seqY_rc, h->lenY, h->idx,
            h->wl, h->soft_mask, h->compare_reverse,
            tile_start, tile_end,
            d_hits, d_hits_sorted, d_hit_flags, d_hit_count, d_filt_count,
            (uint32_t)GECKO_MAX_TILE_HITS, raw_leftmax);
        out->leftmax_hits += raw_leftmax;
        out->extended_hits += n_filt;
        if (n_filt == 0) continue;

        uint32_t n_frags = gpu_extend_tile_two_tier(
            d_hits, n_filt, d_seqX, lenX, h->d_seqY, h->d_seqY_rc, h->lenY,
            h->wl, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
            d_long_seeds, d_long_count, d_out_frags, d_frag_count, (uint32_t)GECKO_MAX_TILE_FRAGS);

        if (n_frags > 0) {
            tile_frags.resize(n_frags);
            CUDA_CHECK(cudaMemcpy(tile_frags.data(), d_out_frags, (size_t)n_frags * sizeof(GpuFrag), cudaMemcpyDeviceToHost));
            all_frags.insert(all_frags.end(), tile_frags.begin(), tile_frags.end());
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaFree(d_seqX));
    CUDA_CHECK(cudaFree(d_hits));
    CUDA_CHECK(cudaFree(d_hits_sorted));
    CUDA_CHECK(cudaFree(d_long_seeds));
    CUDA_CHECK(cudaFree(d_hit_flags));
    CUDA_CHECK(cudaFree(d_out_frags));
    CUDA_CHECK(cudaFree(d_hit_count));
    CUDA_CHECK(cudaFree(d_filt_count));
    CUDA_CHECK(cudaFree(d_long_count));
    CUDA_CHECK(cudaFree(d_frag_count));

    if (!all_frags.empty()) {
        out->frags = (GpuFrag *)malloc(all_frags.size() * sizeof(GpuFrag));
        if (!out->frags) return -2;
        memcpy(out->frags, all_frags.data(), all_frags.size() * sizeof(GpuFrag));
        out->n_frags = all_frags.size();
    }
    return 0;
}

extern "C" void gecko_cuda_free_ref(GeckoCudaRefHandle *h) {
    if (!h) return;
    gpu_free_ref_index(h->idx);
    if (h->d_seqY) cudaFree(h->d_seqY);
    if (h->d_seqY_rc) cudaFree(h->d_seqY_rc);
    delete h;
}

extern "C" void gecko_cuda_free_frags(GpuFrag *frags) {
    free(frags);
}
