/*
 * cuda_bridge.cu
 *
 * Clean 2-phase C-ABI bridge between Rust (gecko-cuda) and 2-bit PackedSeq CUDA kernels:
 *   1. gecko_cuda_probe()
 *   2. gecko_cuda_prepare_ref()  — zero-copy uploads PackedSeq Y, queries L2 cache size,
 *      builds 2-bit Y_rc, single-pass chunk offsets, GPU index (or multi-slab global
 *      dropped_bitmap for F15).
 *   3. gecko_cuda_query_and_extend() — zero-copy uploads PackedSeq X, dynamically sizes
 *      aliased tile scratch (<= 640 MiB max / 20 MiB small genomes), runs adaptive
 *      monolithic-vs-L2-partitioned tiled query + zero-stack segmented warp extension.
 *   4. gecko_cuda_free_ref() / gecko_cuda_free_frags()
 */

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include "cuda_common.cuh"

static uint8_t comp_iupac_host(uint8_t c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; case 'U': return 'A';
        case 'R': return 'Y'; case 'Y': return 'R'; case 'M': return 'K'; case 'K': return 'M';
        case 'S': return 'S'; case 'W': return 'W';
        case 'B': return 'V'; case 'V': return 'B'; case 'D': return 'H'; case 'H': return 'D';
        default: return c;
    }
}

static constexpr size_t PINNED_SLOT_BYTES = 8ULL << 20;      /* 8 MiB per ping-pong slot (16 MiB total) */
static constexpr size_t PINNED_MIN_BYTES  = 64ULL << 10;     /* stage transfers >= 64 KiB */

struct GeckoCudaScratch {
    uint64_t *d_hits = nullptr;
    uint64_t *d_hits_sorted = nullptr;
    uint64_t *d_long_seeds = nullptr; /* F15: aliases d_hits_sorted — NEVER cudaFree separately! */
    uint8_t *d_hit_flags = nullptr;
    GpuFrag *d_out_frags = nullptr;
    uint32_t *d_hit_count = nullptr;
    uint32_t *d_filt_count = nullptr;
    uint32_t *d_long_count = nullptr;
    uint32_t *d_frag_count = nullptr;
    uint32_t *d_stats_buf = nullptr;  /* F16: 2 u32 telemetry counters (valid_q, bypassed_q) */
    uint32_t hits_cap = 0U;
    uint32_t frags_cap = 0U;
    void *d_cub_temp = nullptr;
    size_t cub_temp_cap = 0;
    uint8_t *h_pinned = nullptr;
    size_t h_pinned_slot_cap = 0;
    cudaEvent_t slot_done[2] = {nullptr, nullptr};
    bool pinned_failed = false;
    bool allocated = false;
};

static int ensure_pinned_host(GeckoCudaScratch *sc, size_t needed) {
    if (!sc || sc->pinned_failed) return -1;
    size_t want_slot = std::min<size_t>(needed, PINNED_SLOT_BYTES);
    if (want_slot <= sc->h_pinned_slot_cap && sc->h_pinned) return 0;

    if (sc->h_pinned) {
        cudaFreeHost(sc->h_pinned);
        sc->h_pinned = nullptr;
        sc->h_pinned_slot_cap = 0;
    }
    size_t alloc_slot = (want_slot + (1ULL << 20) - 1ULL) & ~((1ULL << 20) - 1ULL);
    if (alloc_slot > PINNED_SLOT_BYTES) alloc_slot = PINNED_SLOT_BYTES;
    if (cudaHostAlloc(&sc->h_pinned, alloc_slot * 2, cudaHostAllocDefault) == cudaSuccess) {
        if (!sc->slot_done[0]) {
            if (cudaEventCreateWithFlags(&sc->slot_done[0], cudaEventDisableTiming) != cudaSuccess ||
                cudaEventCreateWithFlags(&sc->slot_done[1], cudaEventDisableTiming) != cudaSuccess) {
                cudaFreeHost(sc->h_pinned);
                sc->h_pinned = nullptr;
                sc->pinned_failed = true;
                return -1;
            }
        }
        sc->h_pinned_slot_cap = alloc_slot;
        return 0;
    }
    sc->h_pinned = nullptr;
    sc->h_pinned_slot_cap = 0;
    sc->pinned_failed = true;
    return -1;
}

static int upload_h2d_staged(void *dst, const void *src, size_t bytes, GeckoCudaScratch *sc) {
    if (bytes == 0) return 0;
    if (bytes < PINNED_MIN_BYTES || ensure_pinned_host(sc, bytes) != 0) {
        CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice));
        return 0;
    }
    const uint8_t *s_ptr = static_cast<const uint8_t *>(src);
    uint8_t *d_ptr = static_cast<uint8_t *>(dst);
    size_t slot_cap = sc->h_pinned_slot_cap;
    size_t off = 0;
    int chunk = 0;
    while (off < bytes) {
        size_t n = std::min<size_t>(slot_cap, bytes - off);
        int s = chunk & 1;
        if (chunk >= 2) {
            CUDA_CHECK(cudaEventSynchronize(sc->slot_done[s]));
        }
        uint8_t *buf = sc->h_pinned + (size_t)s * slot_cap;
        memcpy(buf, s_ptr + off, n);
        CUDA_CHECK(cudaMemcpyAsync(d_ptr + off, buf, n, cudaMemcpyHostToDevice, 0));
        CUDA_CHECK(cudaEventRecord(sc->slot_done[s], 0));
        off += n;
        chunk++;
    }
    if (chunk >= 1) CUDA_CHECK(cudaEventSynchronize(sc->slot_done[(chunk - 1) & 1]));
    if (chunk >= 2) CUDA_CHECK(cudaEventSynchronize(sc->slot_done[(chunk - 2) & 1]));
    return 0;
}

static int download_d2h_staged(void *dst, const void *src, size_t bytes, GeckoCudaScratch *sc) {
    if (bytes == 0) return 0;
    if (bytes < PINNED_MIN_BYTES || ensure_pinned_host(sc, bytes) != 0) {
        CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost));
        return 0;
    }
    const uint8_t *s_ptr = static_cast<const uint8_t *>(src);
    uint8_t *d_ptr = static_cast<uint8_t *>(dst);
    size_t slot_cap = sc->h_pinned_slot_cap;

    size_t cur_off = 0;
    size_t cur_n = std::min<size_t>(slot_cap, bytes);
    int s = 0;
    CUDA_CHECK(cudaMemcpyAsync(sc->h_pinned, s_ptr, cur_n, cudaMemcpyDeviceToHost, 0));
    CUDA_CHECK(cudaEventRecord(sc->slot_done[0], 0));

    while (cur_off < bytes) {
        size_t next_off = cur_off + cur_n;
        size_t next_n = (next_off < bytes) ? std::min<size_t>(slot_cap, bytes - next_off) : 0;
        if (next_n > 0) {
            uint8_t *next_buf = sc->h_pinned + (size_t)(s ^ 1) * slot_cap;
            CUDA_CHECK(cudaMemcpyAsync(next_buf, s_ptr + next_off, next_n, cudaMemcpyDeviceToHost, 0));
            CUDA_CHECK(cudaEventRecord(sc->slot_done[s ^ 1], 0));
        }
        CUDA_CHECK(cudaEventSynchronize(sc->slot_done[s]));
        memcpy(d_ptr + cur_off, sc->h_pinned + (size_t)s * slot_cap, cur_n);
        cur_off = next_off;
        cur_n = next_n;
        s ^= 1;
    }
    return 0;
}

static int upload_packed_seq(const GeckoPackedSeqHost *src, GpuSeqBuf &dst, GeckoCudaScratch *sc) {
    memset(&dst, 0, sizeof(dst));
    dst.len = src->len;
    dst.n_seps = src->n_seps;
    dst.n_iupac = src->n_iupac;

    size_t bits_words = (size_t)((src->len + 31ULL) / 32ULL) + 2;
    size_t copy_words = (src->bits_words < bits_words) ? (size_t)src->bits_words : bits_words;
    CUDA_CHECK(cudaMalloc(&dst.d_bits, bits_words * sizeof(uint64_t)));
    if (copy_words > 0 && src->bits) {
        if (upload_h2d_staged(dst.d_bits, src->bits, copy_words * sizeof(uint64_t), sc) != 0) return -1;
    }
    if (bits_words > copy_words) {
        CUDA_CHECK(cudaMemset(dst.d_bits + copy_words, 0, (bits_words - copy_words) * sizeof(uint64_t)));
    }

    size_t bm_words = (size_t)((src->len + 31ULL) / 32ULL) + 2;
    size_t bm_copy = (src->bm_words < bm_words) ? (size_t)src->bm_words : bm_words;
    if (src->special_bm) {
        CUDA_CHECK(cudaMalloc(&dst.d_special_bm, bm_words * sizeof(uint32_t)));
        if (bm_copy > 0) {
            if (upload_h2d_staged(dst.d_special_bm, src->special_bm, bm_copy * sizeof(uint32_t), sc) != 0) return -1;
        }
        if (bm_words > bm_copy) {
            CUDA_CHECK(cudaMemset(dst.d_special_bm + bm_copy, 0, (bm_words - bm_copy) * sizeof(uint32_t)));
        }
    }
    if (src->invalid_bm) {
        if (src->invalid_bm == src->special_bm && dst.d_special_bm) {
            dst.d_invalid_bm = dst.d_special_bm;
            dst.owns_invalid_bm = false;
        } else {
            CUDA_CHECK(cudaMalloc(&dst.d_invalid_bm, bm_words * sizeof(uint32_t)));
            if (bm_copy > 0) {
                if (upload_h2d_staged(dst.d_invalid_bm, src->invalid_bm, bm_copy * sizeof(uint32_t), sc) != 0) return -1;
            }
            if (bm_words > bm_copy) {
                CUDA_CHECK(cudaMemset(dst.d_invalid_bm + bm_copy, 0, (bm_words - bm_copy) * sizeof(uint32_t)));
            }
            dst.owns_invalid_bm = true;
        }
    }

    if (src->n_seps > 0 && src->sep_pos) {
        CUDA_CHECK(cudaMalloc(&dst.d_sep_pos, (size_t)src->n_seps * sizeof(uint32_t)));
        if (upload_h2d_staged(dst.d_sep_pos, src->sep_pos, (size_t)src->n_seps * sizeof(uint32_t), sc) != 0) return -1;
    }
    if (src->n_iupac > 0 && src->iupac_pos && src->iupac_byte) {
        CUDA_CHECK(cudaMalloc(&dst.d_iupac_pos, (size_t)src->n_iupac * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&dst.d_iupac_byte, (size_t)src->n_iupac * sizeof(uint8_t)));
        if (upload_h2d_staged(dst.d_iupac_pos, src->iupac_pos, (size_t)src->n_iupac * sizeof(uint32_t), sc) != 0) return -1;
        if (upload_h2d_staged(dst.d_iupac_byte, src->iupac_byte, (size_t)src->n_iupac * sizeof(uint8_t), sc) != 0) return -1;
    }
    return 0;
}

static void free_gpu_seq(GpuSeqBuf &s) {
    if (s.d_bits) cudaFree(s.d_bits);
    if (s.d_special_bm) cudaFree(s.d_special_bm);
    if (s.owns_invalid_bm && s.d_invalid_bm) cudaFree(s.d_invalid_bm);
    if (s.d_sep_pos) cudaFree(s.d_sep_pos);
    if (s.d_iupac_pos) cudaFree(s.d_iupac_pos);
    if (s.d_iupac_byte) cudaFree(s.d_iupac_byte);
    memset(&s, 0, sizeof(s));
}

/* F15: Dynamic small-genome scratch sizing (20 MiB for 4 kb TYLCV, <= 640 MiB max)
 * with monotonic growth across all-vs-all queries and d_long_seeds = d_hits_sorted aliasing. */
static int ensure_scratch(GeckoCudaScratch &sc, uint64_t lenX) {
    uint32_t want_hits = (uint32_t)GECKO_MAX_TILE_HITS;   /* 32M */
    uint32_t want_frags = (uint32_t)GECKO_MAX_TILE_FRAGS; /* 4M */
    if (lenX < (4ULL << 20)) {
        uint64_t scaled = std::max<uint64_t>(1ULL << 20, lenX * 4ULL);
        uint32_t p2 = 1U << 20;
        while (p2 < scaled && p2 < (uint32_t)GECKO_MAX_TILE_HITS) p2 <<= 1;
        want_hits = p2;
        want_frags = std::max<uint32_t>(128U << 10, want_hits >> 3);
    }
    if (sc.allocated && sc.hits_cap >= want_hits && sc.frags_cap >= want_frags) return 0;

    if (sc.d_hits)        { cudaFree(sc.d_hits);        sc.d_hits = nullptr; }
    if (sc.d_hits_sorted) { cudaFree(sc.d_hits_sorted); sc.d_hits_sorted = nullptr; }
    sc.d_long_seeds = nullptr;
    if (sc.d_hit_flags)   { cudaFree(sc.d_hit_flags);   sc.d_hit_flags = nullptr; }
    if (sc.d_out_frags)   { cudaFree(sc.d_out_frags);   sc.d_out_frags = nullptr; }
    sc.hits_cap = 0U;
    sc.frags_cap = 0U;
    sc.allocated = false;

    CUDA_CHECK(cudaMalloc(&sc.d_hits, (size_t)want_hits * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_hits_sorted, (size_t)want_hits * sizeof(uint64_t)));
    sc.d_long_seeds = sc.d_hits_sorted;
    CUDA_CHECK(cudaMalloc(&sc.d_hit_flags, (size_t)want_hits * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_out_frags, (size_t)want_frags * sizeof(GpuFrag)));
    if (!sc.d_hit_count)  CUDA_CHECK(cudaMalloc(&sc.d_hit_count, sizeof(uint32_t)));
    if (!sc.d_filt_count) CUDA_CHECK(cudaMalloc(&sc.d_filt_count, sizeof(uint32_t)));
    if (!sc.d_long_count) CUDA_CHECK(cudaMalloc(&sc.d_long_count, sizeof(uint32_t)));
    if (!sc.d_frag_count) CUDA_CHECK(cudaMalloc(&sc.d_frag_count, sizeof(uint32_t)));
    if (!sc.d_stats_buf)  CUDA_CHECK(cudaMalloc(&sc.d_stats_buf, 2 * sizeof(uint32_t)));
    sc.hits_cap = want_hits;
    sc.frags_cap = want_frags;
    sc.allocated = true;
    return 0;
}

static void free_scratch(GeckoCudaScratch &sc) {
    if (sc.d_hits) cudaFree(sc.d_hits);
    if (sc.d_hits_sorted) cudaFree(sc.d_hits_sorted);
    /* sc.d_long_seeds aliases sc.d_hits_sorted — NEVER cudaFree(sc.d_long_seeds)! */
    sc.d_long_seeds = nullptr;
    if (sc.d_hit_flags) cudaFree(sc.d_hit_flags);
    if (sc.d_out_frags) cudaFree(sc.d_out_frags);
    if (sc.d_hit_count) cudaFree(sc.d_hit_count);
    if (sc.d_filt_count) cudaFree(sc.d_filt_count);
    if (sc.d_long_count) cudaFree(sc.d_long_count);
    if (sc.d_frag_count) cudaFree(sc.d_frag_count);
    if (sc.d_stats_buf) cudaFree(sc.d_stats_buf);
    if (sc.d_cub_temp) cudaFree(sc.d_cub_temp);
    if (sc.h_pinned) cudaFreeHost(sc.h_pinned);
    if (sc.slot_done[0]) cudaEventDestroy(sc.slot_done[0]);
    if (sc.slot_done[1]) cudaEventDestroy(sc.slot_done[1]);
    sc = GeckoCudaScratch{};
}

struct GeckoCudaRefHandle {
    GpuSeqBuf seqY;
    GpuSeqBuf seqY_rc;
    int wl;
    uint64_t seed_mask;
    uint32_t freq_limit;
    uint32_t num_slabs;
    int l2_cache_bytes;
    bool is_spaced;
    bool soft_mask;
    bool compare_reverse;
    uint32_t *d_global_dropped_bm; /* F15: built once for multi-slab contiguous freq_limit > 0 */
    GpuRefIndex idx;
    GeckoCudaScratch scratch;
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

static int parse_query_mode_env() {
    const char *env = std::getenv("GECKO_CUDA_QUERY_MODE");
    if (!env || env[0] == '\0') return GECKO_QMODE_AUTO;
    if (strcmp(env, "monolithic") == 0) return GECKO_QMODE_MONOLITHIC;
    if (strcmp(env, "l2part") == 0) return GECKO_QMODE_L2_PART;
    return GECKO_QMODE_AUTO;
}

static uint32_t choose_cuda_slab_count(uint64_t total_valid_kmers) {
    if (total_valid_kmers == 0ULL) return 1U;

    uint64_t max_kmers_per_slab = 0ULL;
    if (const char *env = std::getenv("GECKO_CUDA_MAX_SLAB_KMERS")) {
        unsigned long long parsed = std::strtoull(env, nullptr, 10);
        if (parsed >= 100000ULL) {
            max_kmers_per_slab = (uint64_t)parsed;
        }
    }

    if (max_kmers_per_slab == 0ULL) {
        size_t free_vram = 0, total_vram = 0;
        if (cudaMemGetInfo(&free_vram, &total_vram) == cudaSuccess) {
            /* Reserve VRAM for query sequence, aliased tile scratch (<= 640 MiB max:
             * d_hits 256 + d_hits_sorted/d_long_seeds 256 + d_hit_flags 32 + d_out_frags 96),
             * CUB temp, and bitmaps. */
            size_t reserve = std::max<size_t>(896ULL << 20, free_vram / 5ULL);
            size_t avail = (free_vram > reserve) ? (free_vram - reserve) : (256ULL << 20);
            /* Peak GPU build footprint: d_keys_a/b (16 B) + d_vals_a/b (8 B) + CUB temp + keep (~6 B) = ~30 B/k-mer */
            max_kmers_per_slab = (uint64_t)(avail / 30ULL);
            if (max_kmers_per_slab < 10000000ULL) max_kmers_per_slab = 10000000ULL;
            if (max_kmers_per_slab > 0x7FFFFFFFULL) max_kmers_per_slab = 0x7FFFFFFFULL;
        } else {
            return 1U;
        }
    }

    if (total_valid_kmers <= max_kmers_per_slab) return 1U;
    uint64_t slabs = (total_valid_kmers + max_kmers_per_slab - 1ULL) / max_kmers_per_slab;
    if (slabs < 2ULL) slabs = 2ULL;
    if (slabs > 256ULL) slabs = 256ULL;
    return (uint32_t)slabs;
}

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
    const GeckoPackedSeqHost *ref_seq,
    int wl,
    uint64_t seed_mask,
    uint32_t freq_limit,
    int soft_mask,
    int compare_reverse
) {
    if (!ref_seq) return nullptr;
    int dev_count = 0;
    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) {
        return nullptr;
    }

    auto *h = new GeckoCudaRefHandle{};
    h->wl = wl;
    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    h->seed_mask = (seed_mask != 0ULL) ? seed_mask : span_mask;
    h->freq_limit = freq_limit;
    h->num_slabs = 1U;
    h->is_spaced = (h->seed_mask != span_mask);
    h->soft_mask = (soft_mask != 0);
    h->compare_reverse = (compare_reverse != 0);
    h->d_global_dropped_bm = nullptr;

    /* F16: Query hardware L2 cache size once for L2-resident partitioned probe auto-tuning */
    int dev = 0;
    cudaGetDevice(&dev);
    int l2_bytes = 0;
    h->l2_cache_bytes = 4 << 20;
    if (cudaDeviceGetAttribute(&l2_bytes, cudaDevAttrL2CacheSize, dev) == cudaSuccess && l2_bytes > 0) {
        h->l2_cache_bytes = l2_bytes;
    }

    if (upload_packed_seq(ref_seq, h->seqY, &h->scratch) != 0) {
        free_scratch(h->scratch);
        free_gpu_seq(h->seqY);
        delete h;
        return nullptr;
    }

    if (h->compare_reverse) {
        h->seqY_rc.len = ref_seq->len;
        h->seqY_rc.n_seps = ref_seq->n_seps;
        h->seqY_rc.n_iupac = ref_seq->n_iupac;

        size_t bits_words = (size_t)((ref_seq->len + 31ULL) / 32ULL) + 2;
        CUDA_CHECK_VOID(cudaMalloc(&h->seqY_rc.d_bits, bits_words * sizeof(uint64_t)));
        CUDA_CHECK_VOID(cudaMemset(h->seqY_rc.d_bits + (bits_words - 2), 0, 2 * sizeof(uint64_t)));

        if (h->seqY.d_special_bm) {
            size_t bm_words = (size_t)((ref_seq->len + 31ULL) / 32ULL) + 2;
            CUDA_CHECK_VOID(cudaMalloc(&h->seqY_rc.d_special_bm, bm_words * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMemset(h->seqY_rc.d_special_bm + (bm_words - 2), 0, 2 * sizeof(uint32_t)));
        }

        gpu_build_revcomp_2bit(
            h->seqY.d_bits,
            h->seqY_rc.d_bits,
            h->seqY.d_special_bm,
            h->seqY_rc.d_special_bm,
            ref_seq->len);

        if (ref_seq->n_seps > 0 && ref_seq->sep_pos) {
            std::vector<uint32_t> rev_seps(ref_seq->n_seps);
            for (uint32_t k = 0; k < ref_seq->n_seps; k++) {
                rev_seps[k] = (uint32_t)(ref_seq->len - 1ULL - ref_seq->sep_pos[ref_seq->n_seps - 1U - k]);
            }
            CUDA_CHECK_VOID(cudaMalloc(&h->seqY_rc.d_sep_pos, (size_t)ref_seq->n_seps * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMemcpy(h->seqY_rc.d_sep_pos, rev_seps.data(), (size_t)ref_seq->n_seps * sizeof(uint32_t), cudaMemcpyHostToDevice));
        }
        if (ref_seq->n_iupac > 0 && ref_seq->iupac_pos && ref_seq->iupac_byte) {
            std::vector<uint32_t> rev_pos(ref_seq->n_iupac);
            std::vector<uint8_t> rev_byte(ref_seq->n_iupac);
            for (uint32_t k = 0; k < ref_seq->n_iupac; k++) {
                uint32_t src_i = ref_seq->n_iupac - 1U - k;
                rev_pos[k] = (uint32_t)(ref_seq->len - 1ULL - ref_seq->iupac_pos[src_i]);
                rev_byte[k] = comp_iupac_host(ref_seq->iupac_byte[src_i]);
            }
            CUDA_CHECK_VOID(cudaMalloc(&h->seqY_rc.d_iupac_pos, (size_t)ref_seq->n_iupac * sizeof(uint32_t)));
            CUDA_CHECK_VOID(cudaMalloc(&h->seqY_rc.d_iupac_byte, (size_t)ref_seq->n_iupac * sizeof(uint8_t)));
            CUDA_CHECK_VOID(cudaMemcpy(h->seqY_rc.d_iupac_pos, rev_pos.data(), (size_t)ref_seq->n_iupac * sizeof(uint32_t), cudaMemcpyHostToDevice));
            CUDA_CHECK_VOID(cudaMemcpy(h->seqY_rc.d_iupac_byte, rev_byte.data(), (size_t)ref_seq->n_iupac * sizeof(uint8_t), cudaMemcpyHostToDevice));
        }
    }

    /* F15: Single count + ExclusiveSum pass reused directly in single-slab gpu_build_ref_index */
    uint32_t *d_chunk_offsets = nullptr;
    uint32_t num_chunks = 0U;
    uint64_t total_valid = gpu_count_valid_kmers(
        h->seqY.view(), wl, h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
        d_chunk_offsets, num_chunks);
    h->num_slabs = choose_cuda_slab_count(total_valid);
    if (freq_limit == 0U) {
        h->num_slabs = 1U;
    }
    if (h->num_slabs == 1U) {
        h->idx = gpu_build_ref_index(
            h->seqY.view(), wl, h->seed_mask, freq_limit,
            h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
            0U, 256U, d_chunk_offsets, (uint32_t)total_valid, true);
        if (d_chunk_offsets) CUDA_CHECK_VOID(cudaFree(d_chunk_offsets));
    } else {
        if (d_chunk_offsets) CUDA_CHECK_VOID(cudaFree(d_chunk_offsets));
        /* F15: Build global 1-bit dropped bitmap once for multi-slab contiguous when drops are possible */
        if (!h->is_spaced && freq_limit > 0U && (uint64_t)freq_limit < total_valid) {
            h->d_global_dropped_bm = gpu_build_global_dropped_bitmap(
                h->seqY.view(), wl, h->seed_mask, freq_limit, h->num_slabs,
                h->scratch.d_cub_temp, h->scratch.cub_temp_cap);
        }
    }
    return h;
}

extern "C" int gecko_cuda_query_and_extend(
    GeckoCudaRefHandle *h,
    const GeckoPackedSeqHost *query_seq,
    uint32_t min_len,
    float sim_th,
    int match_score,
    int mismatch_penalty,
    int n_penalty,
    GeckoCudaResult *out
) {
    if (!h || !query_seq || !out) return -1;
    memset(out, 0, sizeof(*out));

    if (ensure_scratch(h->scratch, query_seq->len) != 0) {
        free_scratch(h->scratch);
        return -1;
    }

    GpuSeqBuf seqX{};
    if (upload_packed_seq(query_seq, seqX, &h->scratch) != 0) {
        free_gpu_seq(seqX);
        return -1;
    }

    uint64_t *d_hits = h->scratch.d_hits;
    uint64_t *d_hits_sorted = h->scratch.d_hits_sorted;
    uint64_t *d_long_seeds = h->scratch.d_long_seeds;
    uint8_t *d_hit_flags = h->scratch.d_hit_flags;
    GpuFrag *d_out_frags = h->scratch.d_out_frags;
    uint32_t *d_hit_count = h->scratch.d_hit_count;
    uint32_t *d_filt_count = h->scratch.d_filt_count;
    uint32_t *d_long_count = h->scratch.d_long_count;
    uint32_t *d_frag_count = h->scratch.d_frag_count;
    uint32_t *d_stats_buf = h->scratch.d_stats_buf;
    uint32_t hits_cap = h->scratch.hits_cap;
    uint32_t frags_cap = h->scratch.frags_cap;

    GpuSeqView viewX = seqX.view();
    GpuSeqView viewY = h->seqY.view();
    GpuSeqView viewY_rc = h->compare_reverse ? h->seqY_rc.view() : viewY;

    std::vector<GpuFrag> all_frags;

    /* Choose an initial tile size tuned to k-mer density, and adaptively halve
     * any tile that overflows hits_cap or frags_cap so zero hits and zero
     * fragments are ever dropped. */
    uint64_t init_tile_bp = GECKO_QUERY_TILE_BP;
    if (h->wl <= 14) {
        init_tile_bp = 1ULL * 1024ULL * 1024ULL;
    } else if (h->wl <= 19) {
        init_tile_bp = 4ULL * 1024ULL * 1024ULL;
    }
    if (!h->soft_mask && init_tile_bp > 1ULL * 1024ULL * 1024ULL) {
        init_tile_bp >>= 1;
    }

    uint64_t lenX = query_seq->len;
    std::vector<std::pair<uint64_t, uint64_t>> base_tiles;
    for (uint64_t ts = 0; ts < lenX; ts += init_tile_bp) {
        uint64_t te = std::min<uint64_t>(ts + init_tile_bp, lenX);
        base_tiles.push_back({ts, te});
    }
    std::reverse(base_tiles.begin(), base_tiles.end());

    int env_query_mode = parse_query_mode_env();

    for (uint32_t slab = 0; slab < h->num_slabs; slab++) {
        GpuRefIndex transient_idx{};
        const GpuRefIndex *active_idx = &h->idx;
        if (h->num_slabs > 1U) {
            uint32_t part_start = (slab * 256U) / h->num_slabs;
            uint32_t part_end = ((slab + 1U) * 256U) / h->num_slabs;
            if (part_start >= part_end) continue;
            transient_idx = gpu_build_ref_index(
                viewY, h->wl, h->seed_mask, h->freq_limit,
                h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                part_start, part_end, nullptr, 0U, false);
            transient_idx.d_dropped_bitmap = h->d_global_dropped_bm;
            transient_idx.owns_dropped_bitmap = false;
            transient_idx.has_global_dropped_bm = !h->is_spaced;
            active_idx = &transient_idx;
        }

        out->index_kept += active_idx->n_kept;
        out->index_valid += active_idx->n_valid;
        out->bucket_bits = active_idx->bucket_bits;

        /* F16: Regime-optimal adaptive query mode per slab */
        int slab_query_mode = env_query_mode;
        if (slab_query_mode == GECKO_QMODE_AUTO) {
            if (active_idx->is_direct_csr || !active_idx->d_tags) {
                slab_query_mode = GECKO_QMODE_MONOLITHIC;
            } else if (!active_idx->d_ref_unique_bm) {
                slab_query_mode = GECKO_QMODE_L2_PART;
            }
        }

        std::vector<std::pair<uint64_t, uint64_t>> tile_stack = base_tiles;
        while (!tile_stack.empty()) {
            uint64_t tile_start = tile_stack.back().first;
            uint64_t tile_end = tile_stack.back().second;
            tile_stack.pop_back();

            HitKeyLayout key_layout = make_hit_key_layout(
                tile_start, tile_end, viewY.len, h->wl, h->compare_reverse);

            uint32_t raw_leftmax = 0U;
            bool hit_overflow = false;
            uint64_t tile_valid_q = 0ULL;
            uint64_t tile_bypassed_q = 0ULL;
            bool probe_auto_sample = (slab_query_mode == GECKO_QMODE_AUTO);
            int run_mode = probe_auto_sample ? GECKO_QMODE_MONOLITHIC : slab_query_mode;

            uint32_t n_filt = gpu_query_tile_left_maximal(
                viewX, viewY, *active_idx,
                h->wl, h->seed_mask, h->compare_reverse,
                tile_start, tile_end, key_layout,
                h->l2_cache_bytes, run_mode,
                probe_auto_sample ? d_stats_buf : nullptr,
                d_hits, d_hits_sorted, d_hit_flags, d_hit_count, d_filt_count,
                h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                hits_cap, raw_leftmax, hit_overflow,
                probe_auto_sample ? &tile_valid_q : nullptr,
                probe_auto_sample ? &tile_bypassed_q : nullptr);

            if (probe_auto_sample && tile_valid_q >= 65536ULL) {
                if (tile_bypassed_q * 4ULL >= tile_valid_q) {
                    slab_query_mode = GECKO_QMODE_MONOLITHIC;
                } else {
                    slab_query_mode = GECKO_QMODE_L2_PART;
                }
            }

            if (hit_overflow && (tile_end - tile_start) > 1ULL) {
                uint64_t mid = tile_start + (tile_end - tile_start) / 2ULL;
                tile_stack.push_back({mid, tile_end});
                tile_stack.push_back({tile_start, mid});
                continue;
            }
            if (n_filt == 0U) {
                out->leftmax_hits += raw_leftmax;
                continue;
            }

            bool frag_overflow = false;
            uint32_t n_frags = gpu_extend_tile_two_tier(
                d_hits, n_filt, key_layout, viewX, viewY, viewY_rc,
                h->wl, h->is_spaced, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
                d_long_seeds, d_long_count, d_out_frags, d_frag_count,
                frags_cap, frag_overflow);
            if (frag_overflow && (tile_end - tile_start) > 1ULL) {
                uint64_t mid = tile_start + (tile_end - tile_start) / 2ULL;
                tile_stack.push_back({mid, tile_end});
                tile_stack.push_back({tile_start, mid});
                continue;
            }

            out->leftmax_hits += raw_leftmax;
            out->extended_hits += n_filt;

            if (n_frags > 0U) {
                size_t old_n = all_frags.size();
                all_frags.resize(old_n + n_frags);
                if (download_d2h_staged(
                        all_frags.data() + old_n,
                        d_out_frags,
                        (size_t)n_frags * sizeof(GpuFrag),
                        &h->scratch) != 0) {
                    if (h->num_slabs > 1U) gpu_free_ref_index(transient_idx);
                    free_gpu_seq(seqX);
                    return -1;
                }
            }
        }

        if (h->num_slabs > 1U) {
            gpu_free_ref_index(transient_idx);
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    free_gpu_seq(seqX);

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
    free_scratch(h->scratch);
    gpu_free_ref_index(h->idx);
    if (h->d_global_dropped_bm) cudaFree(h->d_global_dropped_bm);
    free_gpu_seq(h->seqY);
    free_gpu_seq(h->seqY_rc);
    delete h;
}

extern "C" void gecko_cuda_free_frags(GpuFrag *frags) {
    free(frags);
}
