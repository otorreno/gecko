/*
 * cuda_common.cuh
 *
 * Shared constants, 2-bit PackedSeq device helpers, tile-relative key layout,
 * and C-ABI types for the modular GECKO CUDA + CUB backend linked into
 * gecko-cuda (Rust).
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
#define GECKO_HIT_STAGE_CAP    512U

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

enum GeckoCudaQueryMode : int {
    GECKO_QMODE_AUTO       = 0,
    GECKO_QMODE_MONOLITHIC = 1,
    GECKO_QMODE_L2_PART    = 2
};

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
    uint32_t *d_tags;           /* Hash-CSR: 32-bit hash tags sorted per bucket; Direct-CSR: nullptr */
    uint32_t *d_vals;           /* (posY << 1) | strand_bit */
    uint32_t *d_bucket_off;     /* num_buckets + 1 CSR offsets */
    uint32_t *d_dropped_bitmap; /* 1-bit/bp dropped bitmap (owned iff owns_dropped_bitmap) */
    uint32_t *d_ref_unique_bm;  /* F16: 1-bit/bp unique non-palindromic ref k-mer bitmap */
    uint32_t n_valid;
    uint32_t n_kept;
    uint32_t num_buckets;       /* Exact bucket count */
    int bucket_bits;            /* 2*wl for Direct-CSR; 14..26 for Hash-CSR */
    int shift_bits;             /* 0 for Direct-CSR; 64 - bucket_bits for Hash-CSR */
    int tag_shift;              /* 0 for Direct-CSR; shift_bits - 32 for Hash-CSR */
    uint32_t part_start;        /* active coarse hash partition range [part_start .. part_end) in 0..256 */
    uint32_t part_end;
    uint32_t bucket_base;       /* 0 for Direct-CSR; part_start << (bucket_bits - 8) for Hash-CSR */
    bool is_direct_csr;         /* F14: true iff contiguous wl <= 12 Direct-CSR */
    bool owns_dropped_bitmap;   /* F15: false when borrowing h->d_global_dropped_bm in multi-slab mode */
    bool has_global_dropped_bm; /* F15: true iff dropped_bitmap reflects global drops across all slabs */
};

/* Compact 24-byte HSP output struct (matches Rust #[repr(C)] GpuFrag identically) */
struct GpuFrag {
    uint32_t xStart;  /* raw data offset in seqX */
    uint32_t yStart;  /* raw data offset in seqY (f) or seqY_rc (r) */
    uint32_t length;  /* ungapped HSP length >= 1 */
    uint32_t ident;   /* exact identity count */
    uint32_t score;   /* peak HSP score */
    uint32_t flags;   /* bit 0 = strand (0 = 'f', 1 = 'r'); bits 1..31 = 0 */
};
static_assert(sizeof(GpuFrag) == 24, "GpuFrag must be exactly 24 bytes");

/* 24-byte value-passed tile-relative hit key descriptor (cmem[0], 0 bytes stack) */
struct HitKeyLayout {
    int64_t  min_diag;     /* (int64_t)min_x - (int64_t)max_y */
    uint32_t min_x;        /* minimum possible posX in tile */
    uint32_t x_mask;       /* (1U << x_bits) - 1U */
    uint32_t diag_mask;    /* (diag_bits == 32) ? 0xFFFFFFFFU : ((1U << diag_bits) - 1U) */
    uint8_t  x_bits;       /* in [1..31] */
    uint8_t  strand_shift; /* x_bits + diag_bits, in [2..63] */
    uint8_t  end_bit;      /* compare_reverse ? (strand_shift + 1) : strand_shift, in [2..64] */
    uint8_t  pad;
};
static_assert(sizeof(HitKeyLayout) == 24, "HitKeyLayout must be exactly 24 bytes");

static inline int ceil_log2_u64(uint64_t n) {
    if (n <= 1ULL) return 1;
    int bits = 0;
    uint64_t v = n - 1ULL;
    while (v > 0ULL) {
        v >>= 1;
        bits++;
    }
    return bits;
}

static inline HitKeyLayout make_hit_key_layout(
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t lenY,
    int wl,
    bool compare_reverse
) {
    HitKeyLayout lay{};
    /* cur_i in [max(tile_start, wl - 1), tile_end - 1], so posX = cur_i + 1 - wl
     * can precede tile_start by up to (wl - 1) bases. */
    uint64_t min_x = (tile_start >= (uint64_t)(wl - 1))
                         ? (tile_start - (uint64_t)(wl - 1))
                         : 0ULL;
    uint64_t max_x = (tile_end >= (uint64_t)wl) ? (tile_end - (uint64_t)wl) : 0ULL;
    uint64_t x_span = (max_x >= min_x) ? (max_x - min_x + 1ULL) : 1ULL;

    uint64_t max_y = (lenY >= (uint64_t)wl) ? (lenY - (uint64_t)wl) : 0ULL;
    int64_t min_diag = (int64_t)min_x - (int64_t)max_y;
    int64_t max_diag = (int64_t)max_x;
    uint64_t diag_span = (max_diag >= min_diag)
                             ? ((uint64_t)(max_diag - min_diag) + 1ULL)
                             : 1ULL;

    int x_bits = ceil_log2_u64(x_span);
    int diag_bits = ceil_log2_u64(diag_span);
    int strand_shift = x_bits + diag_bits;
    int total_bits = strand_shift + (compare_reverse ? 1 : 0);

    if (x_bits > 31 || diag_bits > 32 || strand_shift > 63 || total_bits > 64) {
        lay.min_diag = -(1LL << 31);
        lay.min_x = 0U;
        lay.x_mask = 0x7FFFFFFFU;
        lay.diag_mask = 0xFFFFFFFFU;
        lay.x_bits = 31U;
        lay.strand_shift = 63U;
        lay.end_bit = compare_reverse ? 64U : 63U;
        lay.pad = 0U;
        return lay;
    }

    lay.min_diag = min_diag;
    lay.min_x = (uint32_t)min_x;
    lay.x_mask = (x_bits == 32) ? 0xFFFFFFFFU : ((1U << x_bits) - 1U);
    lay.diag_mask = (diag_bits == 32) ? 0xFFFFFFFFU : ((1U << diag_bits) - 1U);
    lay.x_bits = (uint8_t)x_bits;
    lay.strand_shift = (uint8_t)strand_shift;
    lay.end_bit = (uint8_t)total_bits;
    lay.pad = 0U;
    return lay;
}

static inline void gpu_ensure_cub_temp(void *&d_temp, size_t &cap, size_t needed) {
    if (needed <= cap) return;
    if (d_temp) {
        CUDA_CHECK_VOID(cudaFree(d_temp));
        d_temp = nullptr;
        cap = 0;
    }
    size_t align = (needed < (1ULL << 20)) ? (64ULL << 10) : (4ULL << 20);
    size_t alloc_bytes = (needed + align - 1ULL) & ~(align - 1ULL);
    CUDA_CHECK_VOID(cudaMalloc(&d_temp, alloc_bytes));
    cap = alloc_bytes;
}

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

/* Raw-pointer special-byte decoder (avoids passing GpuSeqView by reference in hot kernels) */
__device__ __forceinline__ uint8_t decode_special_byte_raw(
    const uint32_t * __restrict__ sep_pos,
    uint32_t n_seps,
    const uint32_t * __restrict__ iupac_pos,
    const uint8_t  * __restrict__ iupac_byte,
    uint32_t n_iupac,
    uint32_t pos
) {
    if (n_seps > 0U) {
        uint32_t lo = 0U, hi = n_seps;
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (__ldg(&sep_pos[mid]) < pos) lo = mid + 1U;
            else hi = mid;
        }
        if (lo < n_seps && __ldg(&sep_pos[lo]) == pos) return 0U;
    }
    if (n_iupac > 0U) {
        uint32_t lo = 0U, hi = n_iupac;
        while (lo < hi) {
            uint32_t mid = lo + ((hi - lo) >> 1);
            if (__ldg(&iupac_pos[mid]) < pos) lo = mid + 1U;
            else hi = mid;
        }
        if (lo < n_iupac && __ldg(&iupac_pos[lo]) == pos) return __ldg(&iupac_byte[lo]);
    }
    return (uint8_t)'N';
}

/* Called ONLY when is_bit_set(seq.special_bm, pos) is true. Returns 0 for SEQ_SEP,
 * uppercase IUPAC byte for rare ambiguity codes, or 'N' otherwise. */
__device__ __forceinline__ uint8_t decode_special_byte(const GpuSeqView &seq, uint32_t pos) {
    return decode_special_byte_raw(
        seq.sep_pos, seq.n_seps, seq.iupac_pos, seq.iupac_byte, seq.n_iupac, pos);
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

/* Exact mathematical inverse of MurmurHash3 fmix64 (bijection on uint64_t) */
__device__ __forceinline__ uint64_t fmix64_inv(uint64_t x) {
    x ^= x >> 33;
    x *= 0x9cb4b2f8129337dbULL;
    x ^= x >> 33;
    x *= 0x4f74430c22a54005ULL;
    x ^= x >> 33;
    return x;
}

/* Exact zero-memory-read palindrome test on a canonical masked k-mer (works for
 * all contiguous wl in 4..=32 AND all palindromic spaced seeds, even or odd span). */
__device__ __forceinline__ bool is_palindromic_canon(uint64_t canon, int wl, uint64_t eff_mask) {
    uint64_t rev = ((~reverse_32_2bit_dev(canon)) >> (64 - 2 * wl)) & eff_mask;
    return rev == canon;
}

/* Pack (strand, tile-relative diag, tile-relative posX) into a bit-bounded sortable uint64_t */
__device__ __forceinline__ uint64_t pack_hit_key_tile(
    HitKeyLayout lay, uint32_t is_rev, uint32_t posX, uint32_t posY
) {
    int64_t diag = (int64_t)posX - (int64_t)posY;
    uint64_t rel_diag = (uint64_t)(diag - lay.min_diag);
    uint64_t rel_x = (uint64_t)(posX - lay.min_x);
    return ((uint64_t)(is_rev & 1U) << lay.strand_shift)
         | (rel_diag << lay.x_bits)
         | rel_x;
}

__device__ __forceinline__ void unpack_hit_key_tile(
    HitKeyLayout lay, uint64_t key,
    uint32_t &is_rev, uint32_t &posX, uint32_t &posY, int64_t &diag
) {
    is_rev = (uint32_t)(key >> lay.strand_shift);
    uint32_t rel_x = (uint32_t)(key & (uint64_t)lay.x_mask);
    posX = lay.min_x + rel_x;
    uint64_t rel_diag = (key >> lay.x_bits) & (uint64_t)lay.diag_mask;
    diag = (int64_t)rel_diag + lay.min_diag;
    posY = (uint32_t)((int64_t)posX - diag);
}

__device__ __forceinline__ uint64_t hit_strand_diag_tile(HitKeyLayout lay, uint64_t key) {
    return key >> lay.x_bits;
}

/* ---- Internal module declarations ---- */
void gpu_build_revcomp_2bit(
    const uint64_t *d_src_bits,
    uint64_t *d_dst_bits,
    const uint32_t *d_src_special_bm,
    uint32_t *d_dst_special_bm,
    uint64_t len
);

uint64_t gpu_count_valid_kmers(
    const GpuSeqView &seqY,
    int wl,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t *&out_d_chunk_offsets,
    uint32_t &out_num_chunks
);

uint32_t *gpu_build_global_dropped_bitmap(
    const GpuSeqView &seqY,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit,
    uint32_t num_slabs,
    void *&d_cub_temp,
    size_t &cub_temp_cap
);

GpuRefIndex gpu_build_ref_index(
    const GpuSeqView &seqY,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t part_start = 0U,
    uint32_t part_end = 256U,
    const uint32_t *precomp_chunk_offsets = nullptr,
    uint32_t precomp_n_valid = 0U,
    bool build_dropped_bm = true
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
    HitKeyLayout key_layout,
    int l2_cache_bytes,
    int query_mode,
    uint32_t *d_stats_buf,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint8_t *d_hit_flags,
    uint32_t *d_hit_count,
    uint32_t *d_filt_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed,
    uint64_t *out_valid_q = nullptr,
    uint64_t *out_bypassed_q = nullptr
);

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
    bool &out_overflowed
);

#endif /* GECKO_RUST_CUDA_COMMON_CUH */
