/*
 * cuda_common.cuh
 *
 * Shared constants, 2-bit PackedSeq device helpers, and C-ABI types for the
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

/* Host-side zero-copy descriptor for Rust PackedSeq + sparse bitmaps */
struct GeckoPackedSeqHost {
    const uint64_t *bits;        /* (len + 31)/32 + 2 u64 words (LSB-first 2-bit/bp) */
    uint64_t bits_words;
    const uint32_t *special_bm;  /* NULL or (len + 31)/32 + 2 u32 words (N / SEP / IUPAC) */
    const uint32_t *invalid_bm;  /* NULL or (len + 31)/32 + 2 u32 words (may alias special_bm) */
    uint64_t bm_words;
    const uint32_t *sep_pos;     /* NULL or n_seps sorted SEQ_SEP positions */
    uint32_t n_seps;
    const uint32_t *iupac_pos;   /* NULL or n_iupac sorted IUPAC positions */
    const uint8_t  *iupac_byte;  /* NULL or n_iupac uppercase ASCII IUPAC bytes */
    uint32_t n_iupac;
    uint64_t len;
};

/* Device-side read-only view passed by value into CUDA kernels */
struct GpuSeqView {
    const uint64_t *bits;
    const uint32_t *special_bm;
    const uint32_t *invalid_bm;
    const uint32_t *sep_pos;
    const uint32_t *iupac_pos;
    const uint8_t  *iupac_byte;
    uint64_t len;
    uint32_t n_seps;
    uint32_t n_iupac;
};

/* Device-side owning container for an uploaded or GPU-built 2-bit sequence */
struct GpuSeqBuf {
    uint64_t *d_bits;
    uint32_t *d_special_bm;
    uint32_t *d_invalid_bm;
    bool owns_invalid_bm;
    uint32_t *d_sep_pos;
    uint32_t *d_iupac_pos;
    uint8_t  *d_iupac_byte;
    uint64_t len;
    uint32_t n_seps;
    uint32_t n_iupac;

    GpuSeqView view() const {
        GpuSeqView v{};
        v.bits = d_bits;
        v.special_bm = d_special_bm;
        v.invalid_bm = d_invalid_bm;
        v.sep_pos = d_sep_pos;
        v.iupac_pos = d_iupac_pos;
        v.iupac_byte = d_iupac_byte;
        v.len = len;
        v.n_seps = n_seps;
        v.n_iupac = n_iupac;
        return v;
    }
};

struct GpuRefIndex {
    uint32_t *d_tags;           /* 32-bit hash tags: (uint32_t)(h >> tag_shift), sorted per bucket */
    uint32_t *d_vals;
    uint32_t *d_bucket_off;
    uint32_t *d_dropped_bitmap;
    uint32_t n_valid;
    uint32_t n_kept;
    int bucket_bits;
    int shift_bits;
    int tag_shift;
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

/* ---- Device inline helpers (2-bit PackedSeq & 1-bit bitmaps) ---- */

__device__ __forceinline__ uint32_t base_2bit_dev(const uint64_t * __restrict__ bits, uint32_t pos) {
    uint64_t w = __ldg(&bits[pos >> 5]);
    return (uint32_t)((w >> ((pos & 31U) << 1)) & 3ULL);
}

__device__ __forceinline__ uint64_t load_32_2bit(const uint64_t * __restrict__ bits, uint32_t pos) {
    uint32_t w = pos >> 5;
    uint32_t shift = (pos & 31U) << 1;
    uint64_t lo = __ldg(&bits[w]);
    if (shift == 0U) return lo;
    uint64_t hi = __ldg(&bits[w + 1]);
    return (lo >> shift) | (hi << (64U - shift));
}

__device__ __forceinline__ uint64_t reverse_32_2bit_dev(uint64_t x) {
    uint32_t lo = (uint32_t)x;
    uint32_t hi = (uint32_t)(x >> 32);
    uint32_t rhi = __byte_perm(lo, 0, 0x0123U);
    uint32_t rlo = __byte_perm(hi, 0, 0x0123U);
    uint64_t y = ((uint64_t)rhi << 32) | (uint64_t)rlo;
    y = ((y & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((y >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    y = ((y & 0x3333333333333333ULL) << 2) | ((y >> 2) & 0x3333333333333333ULL);
    return y;
}

__device__ __forceinline__ bool is_bit_set(const uint32_t * __restrict__ bm, uint32_t pos) {
    if (!bm) return false;
    return ((__ldg(&bm[pos >> 5]) >> (pos & 31U)) & 1U) != 0U;
}

__device__ __forceinline__ uint32_t load_32_1bit(const uint32_t * __restrict__ bm, uint32_t pos) {
    if (!bm) return 0U;
    uint32_t w = pos >> 5;
    uint32_t shift = pos & 31U;
    uint32_t lo = __ldg(&bm[w]);
    if (shift == 0U) return lo;
    uint32_t hi = __ldg(&bm[w + 1]);
    return (lo >> shift) | (hi << (32U - shift));
}

__device__ __forceinline__ uint8_t base_2bit_to_ascii(uint32_t b) {
    return (b == 0U) ? (uint8_t)'A' : (b == 1U) ? (uint8_t)'C' : (b == 2U) ? (uint8_t)'G' : (uint8_t)'T';
}

/* Called ONLY when is_bit_set(seq.special_bm, pos) is true. Returns 0 for SEQ_SEP,
 * uppercase IUPAC byte for rare ambiguity codes, or 'N' otherwise. */
__device__ __forceinline__ uint8_t decode_special_byte(const GpuSeqView &seq, uint32_t pos) {
    if (seq.n_seps > 0U) {
        uint32_t lo = 0U, hi = seq.n_seps;
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (__ldg(&seq.sep_pos[mid]) < pos) lo = mid + 1U;
            else hi = mid;
        }
        if (lo < seq.n_seps && __ldg(&seq.sep_pos[lo]) == pos) return 0U;
    }
    if (seq.n_iupac > 0U) {
        uint32_t lo = 0U, hi = seq.n_iupac;
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (__ldg(&seq.iupac_pos[mid]) < pos) lo = mid + 1U;
            else hi = mid;
        }
        if (lo < seq.n_iupac && __ldg(&seq.iupac_pos[lo]) == pos) return __ldg(&seq.iupac_byte[lo]);
    }
    return (uint8_t)'N';
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
void gpu_build_revcomp_2bit(
    const uint64_t *d_src_bits,
    uint64_t *d_dst_bits,
    const uint32_t *d_src_special_bm,
    uint32_t *d_dst_special_bm,
    uint64_t len
);

GpuRefIndex gpu_build_ref_index(
    const GpuSeqView &seqY,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit
);

void gpu_free_ref_index(GpuRefIndex &idx);

uint32_t gpu_query_tile_left_maximal(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuRefIndex &idx,
    int wl,
    uint64_t seed_mask,
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
);

uint32_t gpu_extend_tile_two_tier(
    const uint64_t *d_filtered_hits,
    uint32_t n_hits,
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
    bool &out_overflowed
);

#endif /* GECKO_RUST_CUDA_COMMON_CUH */
