/*
 * cuda_common.cuh
 *
 * Shared constants, device inline helpers, and C-ABI types for the
 * modular GECKO CUDA + CUB backend linked into gecko-cuda (Rust).
 */

#ifndef GECKO_RUST_CUDA_COMMON_CUH
#define GECKO_RUST_CUDA_COMMON_CUH

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

#define GECKO_EXTRACT_CHUNK_BP 4096U
#define GECKO_QUERY_TILE_BP    (16ULL * 1024ULL * 1024ULL)
#define GECKO_MAX_TILE_HITS    (32ULL * 1024ULL * 1024ULL)
#define GECKO_MAX_TILE_FRAGS   (4ULL * 1024ULL * 1024ULL)
#define GECKO_SHORT_CAP        512

#define CUDA_CHECK(call)                                                        \
    do {                                                                        \
        cudaError_t _err = (call);                                              \
        if (_err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n",                        \
                    __FILE__, __LINE__, cudaGetErrorString(_err));              \
            return -1;                                                          \
        }                                                                       \
    } while (0)

#define CUDA_CHECK_VOID(call)                                                   \
    do {                                                                        \
        cudaError_t _err = (call);                                              \
        if (_err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n",                        \
                    __FILE__, __LINE__, cudaGetErrorString(_err));              \
            exit(1);                                                            \
        }                                                                       \
    } while (0)

struct GpuRefIndex {
    uint64_t *d_keys;
    uint32_t *d_vals;
    uint32_t *d_bucket_off;
    uint32_t *d_dropped_bitmap;
    uint32_t n_valid;
    uint32_t n_kept;
    int bucket_bits;
    int shift_bits;
};

struct GpuFrag {
    int64_t diag;
    uint32_t xStart;
    uint32_t yStart;
    uint32_t xEnd;
    uint32_t yEnd;
    uint32_t ident;
    uint32_t score;
    uint8_t strand; /* 0 = 'f', 1 = 'r' */
    uint8_t pad[7];
};

/* ---- Device inline helpers ---- */
__device__ __forceinline__ uint8_t norm_seed_byte(uint8_t c, bool soft_mask) {
    return soft_mask ? c : (uint8_t)(c & 0xDFU);
}

__device__ __forceinline__ uint8_t base_to_2bit(uint8_t c, bool soft_mask) {
    uint8_t u = norm_seed_byte(c, soft_mask);
    if (u == 'A') return 0;
    if (u == 'C') return 1;
    if (u == 'G') return 2;
    if (u == 'T' || u == 'U') return 3;
    return 4;
}

__device__ __forceinline__ bool is_sep_byte(uint8_t c) {
    return (c == 0) || (c == '*');
}

__device__ __forceinline__ bool is_dropped_pos(
    const uint32_t * __restrict__ dropped_bitmap,
    uint32_t pos
) {
    if (!dropped_bitmap) return false;
    return ((__ldg(&dropped_bitmap[pos >> 5]) >> (pos & 31U)) & 1U) != 0U;
}

__device__ __forceinline__ uint8_t comp_base_dev(uint8_t c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'U': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'u': return 'a';
        default:  return c;
    }
}

/* Invertible 64-bit bijection (MurmurHash3 fmix64) */
__device__ __forceinline__ uint64_t fmix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

/* Pack (strand, diag, posX) into a single sortable uint64_t */
__device__ __forceinline__ uint64_t pack_hit_key(uint32_t is_rev, uint32_t posX, uint32_t posY) {
    int64_t diag = (int64_t)posX - (int64_t)posY;
    uint64_t biased_diag = (uint64_t)(uint32_t)(diag + (1LL << 31));
    return ((uint64_t)(is_rev & 1U) << 63) | (biased_diag << 31) | (uint64_t)(posX & 0x7FFFFFFFU);
}

__device__ __forceinline__ void unpack_hit_key(
    uint64_t key, uint32_t &is_rev, uint32_t &posX, uint32_t &posY, int64_t &diag
) {
    is_rev = (uint32_t)(key >> 63);
    posX = (uint32_t)(key & 0x7FFFFFFFU);
    diag = (int64_t)((key >> 31) & 0xFFFFFFFFULL) - (1LL << 31);
    posY = (uint32_t)((int64_t)posX - diag);
}

/* ---- Internal module declarations ---- */
void gpu_build_revcomp(const uint8_t *d_src, uint8_t *d_dst, uint64_t len);
GpuRefIndex gpu_build_ref_index(
    const uint8_t *d_seqY,
    uint64_t lenY,
    int wl,
    uint32_t freq_limit,
    bool soft_mask
);
void gpu_free_ref_index(GpuRefIndex &idx);

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
);

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
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    uint64_t *d_long_seeds,
    uint32_t *d_long_count,
    GpuFrag *d_out_frags,
    uint32_t *d_frag_count,
    uint32_t max_frags,
    bool &out_overflowed
);

#endif /* GECKO_RUST_CUDA_COMMON_CUH */
