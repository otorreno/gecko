/*
 * cuda_common.cuh
 *
 * Shared constants, 2-bit PackedSeq device helpers, tile-relative key layout,
 * and C-ABI types for the modular GECKO CUDA + CUB backend linked into
 * gecko-cuda (Rust).
 */

#ifndef GECKO_RUST_CUDA_COMMON_CUH
#define GECKO_RUST_CUDA_COMMON_CUH

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <utility>
#include <vector>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define GECKO_HD __host__ __device__ __forceinline__
#else
#define GECKO_HD static inline
#endif

#define GECKO_EXTRACT_CHUNK_BP 4096U
#define GECKO_QUERY_TILE_BP    (16ULL * 1024ULL * 1024ULL)
#define GECKO_MAX_TILE_HITS    (32ULL * 1024ULL * 1024ULL)
#define GECKO_MAX_TILE_FRAGS   (4ULL * 1024ULL * 1024ULL)
#define GECKO_SHORT_CAP        512
#define GECKO_HIT_STAGE_CAP    512U

#ifdef __CUDACC__
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
#endif /* __CUDACC__ */

enum GeckoCudaQueryMode : int {
    GECKO_QMODE_AUTO       = 0,
    GECKO_QMODE_MONOLITHIC = 1,
    GECKO_QMODE_L2_PART    = 2,
    GECKO_QMODE_HYBRID     = 3
};

/* W3-F12b: 8-byte half-open [start, end) run descriptor for GPU bitmap materialization */
struct GeckoRun {
    uint32_t start; /* inclusive, half-open [start, end) */
    uint32_t end;   /* exclusive, <= len */
};
static_assert(sizeof(GeckoRun) == 8, "GeckoRun must be exactly 8 bytes");

/* Host-side zero-copy descriptor for Rust PackedSeq + run lists (replaces 60-120 MiB host bitmaps) */
struct GeckoPackedSeqHost {
    const uint64_t *bits;          /* (len + 31)/32 + 2 u64 words (LSB-first 2-bit/bp) */
    uint64_t bits_words;
    const GeckoRun *special_runs;  /* NULL or n_special_runs sorted disjoint [start, end) */
    const GeckoRun *soft_runs;     /* NULL or n_soft_runs sorted disjoint [start, end) */
    const uint32_t *sep_pos;       /* NULL or n_seps sorted SEQ_SEP positions */
    const uint32_t *iupac_pos;     /* NULL or n_iupac sorted IUPAC positions */
    const uint8_t  *iupac_byte;    /* NULL or n_iupac uppercase ASCII IUPAC bytes */
    uint64_t len;
    uint32_t n_special_runs;
    uint32_t n_soft_runs;
    uint32_t n_seps;
    uint32_t n_iupac;
    uint32_t soft_mask;            /* 1 if soft_runs apply to invalid_bm, else 0 */
    uint32_t _pad;
};
static_assert(sizeof(GeckoPackedSeqHost) == 88, "GeckoPackedSeqHost must be exactly 88 bytes");

/* Seeding mode values matching Rust gecko_core::pipeline::SeedingMode explicitly */
enum GeckoCudaSeedingMode : uint8_t {
    GECKO_SEEDING_EXACT     = 0U,
    GECKO_SEEDING_MINIMIZER = 1U,
    GECKO_SEEDING_AUTO      = 2U
};

/* Unified C-ABI config struct shared between Rust (lib.rs) and C++ (cuda_bridge.cu).
 * Strict 8-byte aligned 48-byte layout with zero implicit padding. */
struct GeckoCudaConfig {
    uint64_t seed_mask;        /* off  0: 0 for contiguous, else 2-bit spaced mask */
    uint32_t freq_limit;       /* off  8: max ref occs per strand (0 = drop all, UINT32_MAX = unlimited) */
    uint32_t min_len;          /* off 12: min HSP length in bp (used by query_and_extend) */
    uint32_t two_hit_window;   /* off 16: 0 = disabled, > 0 = inclusive two-hit window bp */
    float    sim_th;           /* off 20: min similarity % [0..100] */
    int32_t  wl;               /* off 24: seed word length / span in [8..32] */
    int32_t  match_score;      /* off 28: > 0 */
    int32_t  mismatch_penalty; /* off 32: >= 0 */
    int32_t  n_penalty;        /* off 36: >= 0 */
    uint8_t  soft_mask;        /* off 40: 0 or 1 */
    uint8_t  compare_reverse;  /* off 41: 0 or 1 */
    uint8_t  seeding_mode;     /* off 42: 0=Exact, 1=Minimizer, 2=Auto */
    uint8_t  _pad0;            /* off 43: must be 0 */
    uint32_t _pad1;            /* off 44: explicit tail padding to 48 bytes, must be 0 */
};
static_assert(sizeof(GeckoCudaConfig) == 48, "GeckoCudaConfig must be exactly 48 bytes");

/* Persistent CudaScratchPool telemetry counters (48 bytes, zero implicit padding) */
struct GeckoCudaScratchPoolStats {
    uint64_t hits;
    uint64_t misses;
    uint64_t recycles;
    uint64_t evictions;
    uint64_t idle_count;
    uint64_t cached_bytes;
};
static_assert(sizeof(GeckoCudaScratchPoolStats) == 48, "GeckoCudaScratchPoolStats must be exactly 48 bytes");

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
    GeckoRun *d_special_runs; /* transient GPU run buffer freed in free_gpu_seq (or post-sync) */
    GeckoRun *d_soft_runs;    /* transient GPU run buffer freed in free_gpu_seq (or post-sync) */
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
    uint32_t *d_tags;              /* Hash-CSR: 32-bit hash tags sorted per bucket; Direct-CSR: nullptr */
    uint32_t *d_vals;              /* (posY << 1) | strand_bit */
    uint32_t *d_bucket_off;        /* num_buckets + 1 CSR offsets */
    uint32_t *d_dropped_bitmap;    /* 1-bit/bp dropped bitmap (owned iff owns_dropped_bitmap) */
    uint32_t *d_ref_non_unique_bm; /* W3-F13a: 1-bit/bp non-unique or palindromic ref k-mer bitmap (single-slab contiguous Hash-CSR only) */
    uint32_t n_valid;
    uint32_t n_kept;
    uint32_t num_buckets;          /* Exact bucket count */
    int bucket_bits;               /* 2*wl for Direct-CSR; 14..26 for Hash-CSR */
    int shift_bits;                /* 0 for Direct-CSR; 64 - bucket_bits for Hash-CSR */
    int tag_shift;                 /* 0 for Direct-CSR; shift_bits - 32 for Hash-CSR */
    uint32_t part_start;           /* active coarse hash partition range [part_start .. part_end) in 0..256 */
    uint32_t part_end;
    uint32_t bucket_base;          /* 0 for Direct-CSR; part_start << (bucket_bits - 8) for Hash-CSR */
    bool is_direct_csr;            /* F14: true iff contiguous wl <= 12 Direct-CSR */
    bool owns_dropped_bitmap;      /* F15: false when borrowing h->d_global_dropped_bm in multi-slab mode */
    bool has_global_dropped_bm;    /* F15: true iff dropped_bitmap reflects global drops across all slabs */
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

/* Unified C-ABI query/extend result struct (56 bytes, explicit filtered_hits for --two-hit) */
struct GeckoCudaResult {
    GpuFrag *frags;         /* off  0 */
    uint64_t n_frags;       /* off  8 */
    uint64_t leftmax_hits;  /* off 16: raw L/R-maximal hits prior to --two-hit */
    uint64_t filtered_hits; /* off 24: post-two-hit surviving hits (= leftmax_hits when disabled) */
    uint64_t extended_hits; /* off 32: hits handed to ungapped extension */
    uint32_t index_kept;    /* off 40 */
    uint32_t index_valid;   /* off 44 */
    int32_t  bucket_bits;   /* off 48 */
    uint32_t _pad;          /* off 52: explicit tail padding to 56 bytes */
};
static_assert(sizeof(GeckoCudaResult) == 56, "GeckoCudaResult must be exactly 56 bytes");

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

/* ---- Pure host/device shared helpers (compilable under both nvcc and g++ -std=c++17) ---- */

/* W3-F11a: Matches gecko_core::profile_enabled() env parsing semantics verbatim */
static inline bool gecko_parse_profile_str(const char *val) {
    if (!val) return false;
    while (*val == ' ' || *val == '\t' || *val == '\r' || *val == '\n') val++;
    size_t len = 0;
    while (val[len] != '\0') len++;
    while (len > 0 && (val[len - 1] == ' ' || val[len - 1] == '\t' ||
                       val[len - 1] == '\r' || val[len - 1] == '\n')) {
        len--;
    }
    if (len == 0) return false;
    if (len == 1 && val[0] == '0') return false;
    auto eq_ci = [&](const char *lit, size_t n) -> bool {
        if (len != n) return false;
        for (size_t i = 0; i < n; i++) {
            char c = val[i];
            if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
            if (c != lit[i]) return false;
        }
        return true;
    };
    if (eq_ci("false", 5) || eq_ci("no", 2) || eq_ci("off", 3)) return false;
    return true;
}

static inline bool gecko_cuda_profile_enabled() {
    static const bool enabled = gecko_parse_profile_str(std::getenv("GECKO_PROFILE"));
    return enabled;
}

/* W3-F12b: Exact run-to-word span calculation shared by GPU fill_runs_bitmap_kernel & host tests */
GECKO_HD bool gecko_run_word_span(
    uint32_t start,
    uint32_t end,
    uint32_t len,
    uint32_t &w_s,
    uint32_t &mask_s,
    uint32_t &w_e,
    uint32_t &mask_e,
    uint32_t &inner_lo,
    uint32_t &inner_hi
) {
    uint32_t s = (start < len) ? start : len;
    uint32_t e = (end   < len) ? end   : len;
    if (s >= e) return false;
    w_s = s >> 5;
    uint32_t b_s = s & 31U;
    uint32_t last = e - 1U;
    w_e = last >> 5;
    uint32_t b_e = last & 31U;
    if (w_s == w_e) {
        mask_s = (0xFFFFFFFFU << b_s) & (0xFFFFFFFFU >> (31U - b_e));
        mask_e = 0U;
        inner_lo = 0U;
        inner_hi = 0U;
    } else {
        mask_s = 0xFFFFFFFFU << b_s;
        mask_e = 0xFFFFFFFFU >> (31U - b_e);
        inner_lo = w_s + 1U;
        inner_hi = w_e;
    }
    return true;
}

static inline void gecko_fill_runs_bitmap_host(
    const GeckoRun *runs,
    uint32_t n_runs,
    uint32_t len,
    uint32_t *bm
) {
    if (!runs || !bm || n_runs == 0U || len == 0U) return;
    for (uint32_t i = 0U; i < n_runs; i++) {
        uint32_t w_s, mask_s, w_e, mask_e, inner_lo, inner_hi;
        if (!gecko_run_word_span(runs[i].start, runs[i].end, len,
                                 w_s, mask_s, w_e, mask_e, inner_lo, inner_hi)) {
            continue;
        }
        bm[w_s] |= mask_s;
        if (w_s != w_e) {
            for (uint32_t w = inner_lo; w < inner_hi; w++) {
                bm[w] = 0xFFFFFFFFU;
            }
            bm[w_e] |= mask_e;
        }
    }
}

/* Invertible 64-bit bijection (MurmurHash3 fmix64) */
GECKO_HD uint64_t gecko_fmix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

/* Exact mathematical inverse of MurmurHash3 fmix64 (bijection on uint64_t) */
GECKO_HD uint64_t gecko_fmix64_inv(uint64_t x) {
    x ^= x >> 33;
    x *= 0x9cb4b2f8129337dbULL;
    x ^= x >> 33;
    x *= 0x4f74430c22a54005ULL;
    x ^= x >> 33;
    return x;
}

/* Reverse 32 packed 2-bit bases in a uint64_t */
GECKO_HD uint64_t gecko_reverse_32_2bit(uint64_t x) {
#ifdef __CUDA_ARCH__
    uint32_t lo = (uint32_t)x;
    uint32_t hi = (uint32_t)(x >> 32);
    uint32_t rhi = __byte_perm(lo, 0, 0x0123U);
    uint32_t rlo = __byte_perm(hi, 0, 0x0123U);
    uint64_t y = ((uint64_t)rhi << 32) | (uint64_t)rlo;
#else
    uint64_t y = __builtin_bswap64(x);
#endif
    y = ((y & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((y >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    y = ((y & 0x3333333333333333ULL) << 2) | ((y >> 2) & 0x3333333333333333ULL);
    return y;
}

/* Exact zero-memory-read palindrome test on a canonical masked k-mer (works for
 * all contiguous wl in 4..=32 AND all palindromic spaced seeds, even or odd span). */
GECKO_HD bool gecko_is_palindromic_canon(uint64_t canon, int wl, uint64_t eff_mask) {
    uint64_t rev = ((~gecko_reverse_32_2bit(canon)) >> (64 - 2 * wl)) & eff_mask;
    return rev == canon;
}

/* ---- W3-F14 / W3-F15: Dual-compiled Minimizer, mix48, 2-bit expand, two-hit & ScratchPool ---- */

static constexpr int      GECKO_MINIMIZER_MIN_K      = 18;
static constexpr uint64_t GECKO_MASK48               = 0x0000FFFFFFFFFFFFULL;
static constexpr uint64_t GECKO_MIN_SALT48           = 0x00009E3779B97F4AULL;
static constexpr uint32_t GECKO_MIN_OVF_SENTINEL     = 0x80000000U;
static constexpr uint32_t GECKO_MIN_TAG_MASK         = 0x7FFFFFFFU;
static constexpr uint32_t GECKO_MIN_OVF_VAL_SENTINEL = 0xFFFFFFFFU;

/* Matches crates/gecko-core/src/minimizer.rs::minimizer_params verbatim:
 * kp in [14..21], w = wl - kp + 1 in [5..12], 2*kp <= 42 <= 48. */
GECKO_HD bool gecko_minimizer_params(int wl, int &out_kp, int &out_w) {
    if (wl < GECKO_MINIMIZER_MIN_K || wl > 32) return false;
    int kp = ((wl + 1) / 2) + 5;
    if (kp > 21) kp = 21;
    out_kp = kp;
    out_w = wl - kp + 1;
    return true;
}

/* M_CAP = max(64, 2 * freq_limit) with uint32 saturation */
GECKO_HD uint32_t gecko_minimizer_m_cap(uint32_t freq_limit) {
    if (freq_limit >= 0x80000000U) return 0xFFFFFFFFU;
    uint32_t dbl = freq_limit * 2U;
    return (dbl > 64U) ? dbl : 64U;
}

/* Host/device activation predicate matching gecko_core::pipeline::CompareConfig::use_minimizer_seeding */
GECKO_HD bool gecko_cuda_should_use_minimizer(
    uint8_t seeding_mode,
    int wl,
    bool is_spaced,
    uint32_t two_hit_window,
    uint32_t freq_limit
) {
    if (seeding_mode != GECKO_SEEDING_MINIMIZER && seeding_mode != GECKO_SEEDING_AUTO) {
        return false;
    }
    if (is_spaced || two_hit_window > 0U || freq_limit == 0U) return false;
    int kp = 0, w = 0;
    return gecko_minimizer_params(wl, kp, w);
}

/* 48-bit bijective secondary hash decorrelated from minimum-skewed fmix64(canon_kp).
 * Both multipliers are odd (gcd(C, 2^48) == 1) with verified modular inverses mod 2^48. */
GECKO_HD uint64_t gecko_mix48(uint64_t canon_kp) {
    uint64_t x = (canon_kp ^ GECKO_MIN_SALT48) & GECKO_MASK48;
    x ^= x >> 24;
    x = (x * 0x0000BF58476D1CE5ULL) & GECKO_MASK48;
    x ^= x >> 24;
    x = (x * 0x000094D049BB1331ULL) & GECKO_MASK48;
    x ^= x >> 24;
    return x;
}

GECKO_HD uint64_t gecko_mix48_inv(uint64_t h48) {
    uint64_t x = h48 & GECKO_MASK48;
    x ^= x >> 24;
    x = (x * 0x00008EC7602365D1ULL) & GECKO_MASK48;
    x ^= x >> 24;
    x = (x * 0x000030B52C8740EDULL) & GECKO_MASK48;
    x ^= x >> 24;
    return (x ^ GECKO_MIN_SALT48) & GECKO_MASK48;
}

/* O(w) zero-stack left/right minimizer window predicate.
 * Given scores[center - (w-1) .. center + (w-1)] with ~0ULL for invalid/OOB kp-mers,
 * computes maximal left/right extents L, R in [0 .. w-1].
 * Emit iff L + R >= w - 1; covered k-mer starts are [p - L .. p + R - (w - 1)]. */
GECKO_HD bool gecko_min_window_eval(
    const uint64_t *scores,
    int center,
    int w,
    bool emit_ties,
    int &out_L,
    int &out_R
) {
    uint64_t s = scores[center];
    if (s == ~0ULL) {
        out_L = 0;
        out_R = 0;
        return false;
    }
    int max_d = w - 1;
    int L = 0;
    for (int d = 1; d <= max_d; d++) {
        uint64_t ls = scores[center - d];
        if (ls == ~0ULL || (emit_ties ? (ls < s) : (ls <= s))) break;
        L = d;
    }
    int R = 0;
    for (int d = 1; d <= max_d; d++) {
        uint64_t rs = scores[center + d];
        if (rs == ~0ULL || rs < s) break;
        R = d;
    }
    out_L = L;
    out_R = R;
    return (L + R >= max_d);
}

/* UB-safe host/device bit-scan & 2-bit/1-bit load helpers */
GECKO_HD uint32_t gecko_clz64(uint64_t x) {
#ifdef __CUDA_ARCH__
    return (uint32_t)__clzll((long long)x);
#else
    return (x == 0ULL) ? 64U : (uint32_t)__builtin_clzll(x);
#endif
}

GECKO_HD uint32_t gecko_clz32(uint32_t x) {
#ifdef __CUDA_ARCH__
    return (uint32_t)__clz((int)x);
#else
    return (x == 0U) ? 32U : (uint32_t)__builtin_clz(x);
#endif
}

GECKO_HD uint32_t gecko_ctz64(uint64_t x) {
#ifdef __CUDA_ARCH__
    return (x == 0ULL) ? 64U : (uint32_t)(__ffsll((long long)x) - 1);
#else
    return (x == 0ULL) ? 64U : (uint32_t)__builtin_ctzll(x);
#endif
}

GECKO_HD uint32_t gecko_ctz32(uint32_t x) {
#ifdef __CUDA_ARCH__
    return (x == 0U) ? 32U : (uint32_t)(__ffs((int)x) - 1);
#else
    return (x == 0U) ? 32U : (uint32_t)__builtin_ctz(x);
#endif
}

GECKO_HD uint32_t gecko_brev32(uint32_t x) {
#ifdef __CUDA_ARCH__
    return __brev(x);
#else
    x = ((x & 0x55555555U) << 1) | ((x >> 1) & 0x55555555U);
    x = ((x & 0x33333333U) << 2) | ((x >> 2) & 0x33333333U);
    x = ((x & 0x0F0F0F0FU) << 4) | ((x >> 4) & 0x0F0F0F0FU);
    return __builtin_bswap32(x);
#endif
}

GECKO_HD uint64_t gecko_load_32_2bit(const uint64_t *bits, uint32_t pos) {
    uint32_t w = pos >> 5;
    uint32_t shift = (pos & 31U) << 1;
#ifdef __CUDA_ARCH__
    uint64_t lo = __ldg(&bits[w]);
    if (shift == 0U) return lo;
    uint64_t hi = __ldg(&bits[w + 1]);
#else
    uint64_t lo = bits[w];
    if (shift == 0U) return lo;
    uint64_t hi = bits[w + 1];
#endif
    return (lo >> shift) | (hi << (64U - shift));
}

GECKO_HD uint32_t gecko_load_32_1bit(const uint32_t *bm, uint32_t pos) {
    if (!bm) return 0U;
    uint32_t w = pos >> 5;
    uint32_t shift = pos & 31U;
#ifdef __CUDA_ARCH__
    uint32_t lo = __ldg(&bm[w]);
    if (shift == 0U) return lo;
    uint32_t hi = __ldg(&bm[w + 1]);
#else
    uint32_t lo = bm[w];
    if (shift == 0U) return lo;
    uint32_t hi = bm[w + 1];
#endif
    return (lo >> shift) | (hi << (32U - shift));
}

GECKO_HD bool gecko_is_bit_set(const uint32_t *bm, uint32_t pos) {
    if (!bm) return false;
#ifdef __CUDA_ARCH__
    return ((__ldg(&bm[pos >> 5]) >> (pos & 31U)) & 1U) != 0U;
#else
    return ((bm[pos >> 5] >> (pos & 31U)) & 1U) != 0U;
#endif
}

/* Pack / unpack tile-relative hit keys (dual-compiled for GPU kernels & host parity tests) */
GECKO_HD uint64_t gecko_pack_hit_key_tile(
    HitKeyLayout lay, uint32_t is_rev, uint32_t posX, uint32_t posY
) {
    int64_t diag = (int64_t)posX - (int64_t)posY;
    uint64_t rel_diag = (uint64_t)(diag - lay.min_diag);
    uint64_t rel_x = (uint64_t)(posX - lay.min_x);
    return ((uint64_t)(is_rev & 1U) << lay.strand_shift)
         | (rel_diag << lay.x_bits)
         | rel_x;
}

GECKO_HD void gecko_unpack_hit_key_tile(
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

/* W3-F14: Dual-compiled 2-bit exact-k expansion around a matching canonical kp-mer (p_q, p_r).
 * Expands up to need_left + 1 <= w <= 12 left and need_right + 1 <= w <= 12 right so
 * left/right maximality at [r_lo .. r_hi] needs ZERO extra sequence loads. */
template <bool IS_REV>
GECKO_HD void gecko_min_expand_extents(
    const uint64_t *seqX_bits,
    const uint32_t *seqX_inv,
    uint64_t lenX,
    const uint64_t *seqY_bits,
    const uint32_t *seqY_inv,
    uint64_t lenY,
    int kp,
    uint32_t p_q,
    uint32_t p_r,
    uint32_t need_left,
    uint32_t need_right,
    uint32_t &out_ext_left,
    uint32_t &out_ext_right
) {
    uint32_t ext_left = 0U;
    uint32_t ext_right = 0U;
    uint32_t kp_u = (uint32_t)kp;
    uint32_t want_l = need_left + 1U;
    uint32_t want_r = need_right + 1U;

    if constexpr (!IS_REV) {
        uint32_t avail_l = (p_q < p_r) ? p_q : p_r;
        uint32_t max_l = (want_l < avail_l) ? want_l : avail_l;
        if (max_l > 0U) {
            uint64_t diff_l = (gecko_load_32_2bit(seqX_bits, p_q - max_l) ^
                               gecko_load_32_2bit(seqY_bits, p_r - max_l)) << (64U - 2U * max_l);
            uint32_t inv_l = (gecko_load_32_1bit(seqX_inv, p_q - max_l) |
                              gecko_load_32_1bit(seqY_inv, p_r - max_l)) << (32U - max_l);
            uint32_t ml = gecko_clz64(diff_l) >> 1;
            uint32_t vl = gecko_clz32(inv_l);
            uint32_t lim = (ml < vl) ? ml : vl;
            ext_left = (max_l < lim) ? max_l : lim;
        }

        uint32_t rem_x = (uint32_t)(lenX - (uint64_t)p_q - (uint64_t)kp_u);
        uint32_t rem_y = (uint32_t)(lenY - (uint64_t)p_r - (uint64_t)kp_u);
        uint32_t avail_r = (rem_x < rem_y) ? rem_x : rem_y;
        uint32_t max_r = (want_r < avail_r) ? want_r : avail_r;
        if (max_r > 0U) {
            uint64_t diff_r = gecko_load_32_2bit(seqX_bits, p_q + kp_u) ^
                              gecko_load_32_2bit(seqY_bits, p_r + kp_u);
            uint32_t inv_r = gecko_load_32_1bit(seqX_inv, p_q + kp_u) |
                             gecko_load_32_1bit(seqY_inv, p_r + kp_u);
            uint32_t mr = gecko_ctz64(diff_r) >> 1;
            uint32_t vr = gecko_ctz32(inv_r);
            uint32_t lim = (mr < vr) ? mr : vr;
            ext_right = (max_r < lim) ? max_r : lim;
        }
    } else {
        uint32_t rem_y = (uint32_t)(lenY - (uint64_t)p_r - (uint64_t)kp_u);
        uint32_t avail_l = (p_q < rem_y) ? p_q : rem_y;
        uint32_t max_l = (want_l < avail_l) ? want_l : avail_l;
        if (max_l > 0U) {
            uint64_t xw_l = gecko_load_32_2bit(seqX_bits, p_q - max_l) << (64U - 2U * max_l);
            uint64_t y_rc_l = ~gecko_reverse_32_2bit(gecko_load_32_2bit(seqY_bits, p_r + kp_u));
            uint64_t diff_l = xw_l ^ y_rc_l;
            uint32_t inv_l = (gecko_load_32_1bit(seqX_inv, p_q - max_l) << (32U - max_l)) |
                             gecko_brev32(gecko_load_32_1bit(seqY_inv, p_r + kp_u));
            uint32_t ml = gecko_clz64(diff_l) >> 1;
            uint32_t vl = gecko_clz32(inv_l);
            uint32_t lim = (ml < vl) ? ml : vl;
            ext_left = (max_l < lim) ? max_l : lim;
        }

        uint32_t rem_x = (uint32_t)(lenX - (uint64_t)p_q - (uint64_t)kp_u);
        uint32_t avail_r = (rem_x < p_r) ? rem_x : p_r;
        uint32_t max_r = (want_r < avail_r) ? want_r : avail_r;
        if (max_r > 0U) {
            uint64_t y_rc_r = (~gecko_reverse_32_2bit(gecko_load_32_2bit(seqY_bits, p_r - max_r)))
                              >> (64U - 2U * max_r);
            uint64_t diff_r = gecko_load_32_2bit(seqX_bits, p_q + kp_u) ^ y_rc_r;
            uint32_t inv_r = gecko_load_32_1bit(seqX_inv, p_q + kp_u) |
                             (gecko_brev32(gecko_load_32_1bit(seqY_inv, p_r - max_r)) >> (32U - max_r));
            uint32_t mr = gecko_ctz64(diff_r) >> 1;
            uint32_t vr = gecko_ctz32(inv_r);
            uint32_t lim = (mr < vr) ? mr : vr;
            ext_right = (max_r < lim) ? max_r : lim;
        }
    }
    out_ext_left = ext_left;
    out_ext_right = ext_right;
}

/* W3-F15b: Dual-compiled BLAST/LASTZ sorted-diagonal two-hit backward scan (0-stack).
 * Checks whether sorted hit `idx` has any prior non-overlapping hit `j < idx` on the same
 * (strand, diag) with wl <= cur_rx - prev_rx <= wl + two_hit_window. */
GECKO_HD uint8_t gecko_two_hit_scan_backward(
    const uint64_t *sorted_hits,
    uint32_t idx,
    HitKeyLayout lay,
    uint64_t wl_u64,
    uint64_t max_dist
) {
    uint64_t key = sorted_hits[idx];
    uint64_t sd = key >> lay.x_bits;
    uint32_t cur_rx = (uint32_t)(key & (uint64_t)lay.x_mask);
    uint32_t j = idx;
    while (j > 0U) {
        j--;
        uint64_t pkey = sorted_hits[j];
        if ((pkey >> lay.x_bits) != sd) break;
        uint32_t prev_rx = (uint32_t)(pkey & (uint64_t)lay.x_mask);
        uint64_t dist = (uint64_t)cur_rx - (uint64_t)prev_rx;
        if (dist > max_dist) break;
        if (dist >= wl_u64) return 1U;
    }
    return 0U;
}

/* W3-F15b: 4-instruction backward self-2hit check for contiguous L/R-maximal MEMs (L >= 2*wl).
 * Recovers interior-suppressed exact matches where only endpoints [a, b] were emitted. */
GECKO_HD uint8_t gecko_two_hit_self_prev(
    uint64_t key,
    HitKeyLayout lay,
    int wl,
    const uint64_t *seqX_bits,
    const uint32_t *seqX_invalid_bm,
    const uint64_t *seqY_bits,
    const uint32_t *seqY_invalid_bm,
    const uint64_t *seqY_rc_bits,
    const uint32_t *dropped_bitmap,
    uint64_t lenY
) {
    uint32_t is_rev, posX, posY_s;
    int64_t diag;
    gecko_unpack_hit_key_tile(lay, key, is_rev, posX, posY_s, diag);
    uint32_t wl_u = (uint32_t)wl;
    if (posX < wl_u || posY_s < wl_u) return 0U;

    uint32_t prev_x = posX - wl_u;
    uint32_t valid_mask_wl = (wl == 32) ? 0xFFFFFFFFU : ((1U << wl) - 1U);
    if ((gecko_load_32_1bit(seqX_invalid_bm, prev_x) & valid_mask_wl) != 0U) return 0U;

    uint32_t prev_y_s = posY_s - wl_u;
    uint32_t prev_y_orig = (is_rev == 0U)
        ? prev_y_s
        : (uint32_t)(lenY - (uint64_t)posY_s); /* (lenY - posY_s - wl) + wl */
    if ((gecko_load_32_1bit(seqY_invalid_bm, prev_y_orig) & valid_mask_wl) != 0U) return 0U;
    if (gecko_is_bit_set(dropped_bitmap, prev_y_orig)) return 0U;

    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    const uint64_t *y_bits = (is_rev == 0U) ? seqY_bits : seqY_rc_bits;
    uint64_t xw = gecko_load_32_2bit(seqX_bits, prev_x);
    uint64_t yw = gecko_load_32_2bit(y_bits, prev_y_s);
    return (((xw ^ yw) & span_mask) == 0ULL) ? 1U : 0U;
}

/* W3-F15a: Thread-safe, lock-free-teardown generic scratch pool core (zero CUDA dependency so
 * host g++ unit & concurrency tests verify the exact same state machine as cuda_bridge.cu).
 * NEVER invokes free_fn while holding mu_. */
template <typename T, typename BytesFn, typename FreeFn>
class GeckoScratchPoolCore {
public:
    GeckoScratchPoolCore(BytesFn bytes_fn, FreeFn free_fn, bool initial_enabled = true, size_t initial_max_idle = 1)
        : bytes_fn_(bytes_fn), free_fn_(free_fn), max_idle_(initial_max_idle), enabled_(initial_enabled) {}

    bool checkout(T &out) {
        std::lock_guard<std::mutex> lk(mu_);
        if (enabled_ && !idle_.empty()) {
            out = std::move(idle_.back());
            idle_.pop_back();
            stats_.hits++;
            update_idle_stats_locked();
            return true;
        }
        stats_.misses++;
        out = T{};
        return false;
    }

    bool checkin(T &&item, bool healthy) {
        bool should_free = false;
        T to_free{};
        {
            std::lock_guard<std::mutex> lk(mu_);
            if (healthy && enabled_ && max_idle_ > 0 && idle_.size() < max_idle_) {
                idle_.push_back(std::move(item));
                stats_.recycles++;
                update_idle_stats_locked();
                return true;
            }
            stats_.evictions++;
            to_free = std::move(item);
            should_free = true;
        }
        if (should_free) {
            free_fn_(to_free);
        }
        return false;
    }

    size_t release_all() {
        std::vector<T> to_free;
        {
            std::lock_guard<std::mutex> lk(mu_);
            to_free.swap(idle_);
            stats_.evictions += (uint64_t)to_free.size();
            update_idle_stats_locked();
        }
        for (auto &e : to_free) {
            free_fn_(e);
        }
        return to_free.size();
    }

    void set_enabled(bool en) {
        std::vector<T> to_free;
        {
            std::lock_guard<std::mutex> lk(mu_);
            enabled_ = en;
            if (!enabled_ && !idle_.empty()) {
                to_free.swap(idle_);
                stats_.evictions += (uint64_t)to_free.size();
                update_idle_stats_locked();
            }
        }
        for (auto &e : to_free) {
            free_fn_(e);
        }
    }

    bool is_enabled() const {
        std::lock_guard<std::mutex> lk(mu_);
        return enabled_;
    }

    void set_max_idle(size_t cap) {
        std::vector<T> to_free;
        {
            std::lock_guard<std::mutex> lk(mu_);
            max_idle_ = cap;
            while (idle_.size() > max_idle_) {
                to_free.push_back(std::move(idle_.back()));
                idle_.pop_back();
            }
            stats_.evictions += (uint64_t)to_free.size();
            update_idle_stats_locked();
        }
        for (auto &e : to_free) {
            free_fn_(e);
        }
    }

    size_t idle_count() const {
        std::lock_guard<std::mutex> lk(mu_);
        return idle_.size();
    }

    GeckoCudaScratchPoolStats stats() const {
        std::lock_guard<std::mutex> lk(mu_);
        return stats_;
    }

private:
    void update_idle_stats_locked() {
        stats_.idle_count = (uint64_t)idle_.size();
        uint64_t total = 0ULL;
        for (const auto &e : idle_) {
            total += bytes_fn_(e);
        }
        stats_.cached_bytes = total;
    }

    mutable std::mutex mu_;
    std::vector<T> idle_;
    BytesFn bytes_fn_;
    FreeFn free_fn_;
    size_t max_idle_;
    bool enabled_;
    GeckoCudaScratchPoolStats stats_{};
};

#ifdef __CUDACC__

/* W3-F14: Dedicated GPU Canonical Minimizer Hash-CSR + High-Freq Exact-k Overflow CSR.
 * Main CSR uses 48-bit bijective mix48 with bucket_bits in [17..26], leaving bit 31 of d_tags
 * strictly 0 for normal tags (0.000% kp-mer tag collisions) and reserving bit 31 as
 * GECKO_MIN_OVF_SENTINEL for buckets routed to the tiny exact-k overflow CSR. */
struct GpuMinIndex {
    uint32_t *d_tags;             /* [n_kept_main], low (48 - bucket_bits) tag | (0x80000000U if ovf) */
    uint32_t *d_vals;             /* [n_kept_main], (p_r << 1) | r_strand (0xFFFFFFFFU for sentinel) */
    uint32_t *d_bucket_off;       /* [num_buckets + 1] */
    uint32_t *d_dropped_bitmap;   /* 1-bit/bp exact k-mer freq_limit dropped bitmap (nullable) */
    uint32_t n_raw_min;
    uint32_t n_kept_main;
    uint32_t num_buckets;
    int bucket_bits;              /* in [17 .. 26], guarantees shift_bits = 48 - bucket_bits in [22 .. 31] */
    int shift_bits;               /* 48 - bucket_bits */
    uint32_t tag_mask;            /* (1U << shift_bits) - 1U */
    int kp;
    int w;
    uint32_t m_cap;               /* max(64, 2 * freq_limit) */

    /* Tiny exact-k Overflow Hash-CSR (64-bit fmix64 exact k-mer hashing; populated iff ovf_n_kept > 0) */
    uint32_t *d_ovf_tags;         /* [ovf_n_kept] */
    uint32_t *d_ovf_vals;         /* [ovf_n_kept], (r_y << 1) | k_strand */
    uint32_t *d_ovf_bucket_off;   /* [ovf_num_buckets + 1] */
    uint32_t ovf_n_kept;
    uint32_t ovf_num_buckets;
    uint32_t ovf_routed_keys;     /* count of overflowed minimizer keys (GECKO_PROFILE) */
    int ovf_bucket_bits;          /* in [10 .. 24] */
    int ovf_shift_bits;           /* 64 - ovf_bucket_bits */
    int ovf_tag_shift;            /* ovf_shift_bits - 32 */
};

/* W3-F11a: Zero-overhead-when-unset CUDA Event stage-by-stage profiler (stderr only) */
struct GeckoCudaProfile {
    bool enabled = false;
    cudaEvent_t ev_start = nullptr;
    cudaEvent_t ev_stop = nullptr;
    /* Dedicated async-stream events so profiling never forces host sync on h2d_stream */
    cudaEvent_t h2d_start = nullptr;
    cudaEvent_t h2d_stop = nullptr;
    cudaEvent_t h2d_bm_start = nullptr;
    cudaEvent_t h2d_bm_stop = nullptr;
    bool h2d_recorded = false;
    bool h2d_bm_recorded = false;

    /* Prepare / Ref index GPU event ms */
    double upload_y_ms = 0.0;
    double fill_bm_y_ms = 0.0;
    double revcomp_ms = 0.0;
    double count_valid_ms = 0.0;
    double extract_kmers_ms = 0.0;
    double sort_pairs_ms = 0.0;
    double freq_filter_ms = 0.0;
    double select_vals_ms = 0.0;
    double select_tags_ms = 0.0;
    double bucket_build_ms = 0.0;
    double global_dropped_ms = 0.0;
    /* W3-F14: Minimizer prepare timers */
    double min_extract_ms = 0.0;
    double min_hf_ms = 0.0;
    double min_ovf_ms = 0.0;

    /* Query & Extend GPU event ms */
    double ensure_scratch_ms = 0.0;
    double upload_x_ms = 0.0;
    double fill_bm_x_ms = 0.0;
    double slab_rebuild_ms = 0.0;
    double l2_pre_bypass_ms = 0.0;
    double l2_count_ms = 0.0;
    double l2_scatter_ms = 0.0;
    double l2_probe_ms = 0.0;
    double mono_probe_ms = 0.0;
    double sort_hits_ms = 0.0;
    double two_hit_ms = 0.0;
    double extend_short_ms = 0.0;
    double extend_long_ms = 0.0;
    double d2h_frags_ms = 0.0;

    /* Ref-level metadata (preserved across repeated CudaPreparedRef::query calls) */
    bool is_minimizer = false;
    int min_kp = 0;
    int min_w = 0;
    uint32_t min_ref_occs = 0U;
    uint32_t min_ovf_buckets = 0U;
    uint32_t min_ovf_kmers = 0U;
    bool scratch_pool_hit = false;

    /* Query metadata */
    uint32_t l2_part_P = 1U;
    uint32_t num_tiles = 0U;
    uint32_t two_hit_window = 0U;
    uint64_t two_hit_dropped = 0ULL;
    uint64_t tile0_valid_q = 0ULL;
    uint64_t tile0_bypassed_q = 0ULL;
    uint64_t long_seeds = 0ULL;
    bool upload_x_overlapped = false;

    void init(bool enable) {
        enabled = enable;
        if (!enabled) return;
        if (cudaEventCreate(&ev_start) != cudaSuccess ||
            cudaEventCreate(&ev_stop) != cudaSuccess ||
            cudaEventCreate(&h2d_start) != cudaSuccess ||
            cudaEventCreate(&h2d_stop) != cudaSuccess ||
            cudaEventCreate(&h2d_bm_start) != cudaSuccess ||
            cudaEventCreate(&h2d_bm_stop) != cudaSuccess) {
            destroy();
        }
    }

    void destroy() {
        if (ev_start)     { cudaEventDestroy(ev_start);     ev_start = nullptr; }
        if (ev_stop)      { cudaEventDestroy(ev_stop);      ev_stop = nullptr; }
        if (h2d_start)    { cudaEventDestroy(h2d_start);    h2d_start = nullptr; }
        if (h2d_stop)     { cudaEventDestroy(h2d_stop);     h2d_stop = nullptr; }
        if (h2d_bm_start) { cudaEventDestroy(h2d_bm_start); h2d_bm_start = nullptr; }
        if (h2d_bm_stop)  { cudaEventDestroy(h2d_bm_stop);  h2d_bm_stop = nullptr; }
        enabled = false;
    }

    void reset_query_counters() {
        ensure_scratch_ms = 0.0;
        upload_x_ms = 0.0;
        fill_bm_x_ms = 0.0;
        slab_rebuild_ms = 0.0;
        l2_pre_bypass_ms = 0.0;
        l2_count_ms = 0.0;
        l2_scatter_ms = 0.0;
        l2_probe_ms = 0.0;
        mono_probe_ms = 0.0;
        sort_hits_ms = 0.0;
        two_hit_ms = 0.0;
        extend_short_ms = 0.0;
        extend_long_ms = 0.0;
        d2h_frags_ms = 0.0;
        h2d_recorded = false;
        h2d_bm_recorded = false;
        l2_part_P = 1U;
        num_tiles = 0U;
        two_hit_window = 0U;
        two_hit_dropped = 0ULL;
        tile0_valid_q = 0ULL;
        tile0_bypassed_q = 0ULL;
        long_seeds = 0ULL;
        upload_x_overlapped = false;
    }
};

/* Non-nesting RAII CUDA Event scope timer (0 overhead when prof == nullptr or !enabled) */
struct CudaEventScope {
    GeckoCudaProfile *prof;
    double *accum_ms;
    cudaStream_t stream;

    CudaEventScope(GeckoCudaProfile *p, double *acc, cudaStream_t s = 0)
        : prof((p && p->enabled) ? p : nullptr), accum_ms(acc), stream(s) {
        if (prof) {
            (void)cudaEventRecord(prof->ev_start, stream);
        }
    }
    ~CudaEventScope() {
        if (prof) {
            float ms = 0.0f;
            if (cudaEventRecord(prof->ev_stop, stream) == cudaSuccess &&
                cudaEventSynchronize(prof->ev_stop) == cudaSuccess &&
                cudaEventElapsedTime(&ms, prof->ev_start, prof->ev_stop) == cudaSuccess) {
                *accum_ms += (double)ms;
            }
        }
    }
};

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
    return gecko_reverse_32_2bit(x);
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

__device__ __forceinline__ uint64_t fmix64(uint64_t x) {
    return gecko_fmix64(x);
}

__device__ __forceinline__ uint64_t fmix64_inv(uint64_t x) {
    return gecko_fmix64_inv(x);
}

__device__ __forceinline__ bool is_palindromic_canon(uint64_t canon, int wl, uint64_t eff_mask) {
    return gecko_is_palindromic_canon(canon, wl, eff_mask);
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
int gpu_fill_runs_bitmap(
    const GeckoRun *d_runs,
    uint32_t n_runs,
    uint64_t len,
    uint32_t *d_bm,
    cudaStream_t stream = 0
);

void gpu_build_revcomp_2bit(
    const uint64_t *d_src_bits,
    uint64_t *d_dst_bits,
    const uint32_t *d_src_special_bm,
    uint32_t *d_dst_special_bm,
    uint64_t len,
    GeckoCudaProfile *prof = nullptr
);

uint64_t gpu_count_valid_kmers(
    const GpuSeqView &seqY,
    int wl,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t *&out_d_chunk_offsets,
    uint32_t &out_num_chunks,
    GeckoCudaProfile *prof = nullptr
);

uint32_t *gpu_build_global_dropped_bitmap(
    const GpuSeqView &seqY,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit,
    uint32_t num_slabs,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    GeckoCudaProfile *prof = nullptr
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
    bool build_dropped_bm = true,
    GeckoCudaProfile *prof = nullptr
);

void gpu_free_ref_index(GpuRefIndex &idx);

uint32_t gpu_query_tile_left_maximal(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuSeqView &seqY_rc,
    const GpuRefIndex &idx,
    int wl,
    uint64_t seed_mask,
    bool compare_reverse,
    uint32_t two_hit_window,
    uint64_t tile_start,
    uint64_t tile_end,
    HitKeyLayout key_layout,
    int l2_cache_bytes,
    int query_mode,
    uint32_t *d_stats_buf,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    GpuFrag *d_out_frags,
    uint32_t *d_hit_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed,
    uint64_t *&out_sorted_hits,
    uint64_t *&out_alt_buf,
    uint64_t *out_valid_q = nullptr,
    uint64_t *out_bypassed_q = nullptr,
    uint32_t *out_two_hit_dropped = nullptr,
    GeckoCudaProfile *prof = nullptr
);

/* W3-F14: GPU Canonical Minimizer Super-Seeding (with HF Exact-k Overflow CSR) */
bool gpu_build_min_index(
    const GpuSeqView &seqY,
    int wl,
    uint32_t freq_limit,
    uint64_t owned_scratch_bytes,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    GpuMinIndex &out_idx,
    GeckoCudaProfile *prof = nullptr
);

void gpu_free_min_index(GpuMinIndex &idx);

uint32_t gpu_query_tile_minimizer(
    const GpuSeqView &seqX,
    const GpuSeqView &seqY,
    const GpuMinIndex &idx,
    int wl,
    bool compare_reverse,
    uint64_t tile_start,
    uint64_t tile_end,
    HitKeyLayout key_layout,
    int l2_cache_bytes,
    int query_mode,
    uint64_t *d_hits,
    uint64_t *d_hits_sorted,
    uint32_t *d_hit_count,
    void *&d_cub_temp,
    size_t &cub_temp_cap,
    uint32_t max_hits,
    uint32_t &out_raw_leftmax,
    bool &out_overflowed,
    uint64_t *&out_sorted_hits,
    uint64_t *&out_alt_buf,
    GeckoCudaProfile *prof = nullptr
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
    bool &out_overflowed,
    GeckoCudaProfile *prof = nullptr
);

#endif /* __CUDACC__ */

#endif /* GECKO_RUST_CUDA_COMMON_CUH */
