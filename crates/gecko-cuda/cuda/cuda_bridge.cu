/*
 * cuda_bridge.cu
 *
 * Clean 2-phase C-ABI bridge between Rust (gecko-cuda) and 2-bit PackedSeq CUDA kernels:
 *   1. gecko_cuda_probe() / gecko_cuda_warmup()
 *   2. gecko_cuda_prepare_ref()  — zero-copy uploads PackedSeq Y + runs (GPU fill_runs_bitmap),
 *      queries L2 cache size, builds 2-bit Y_rc, single-pass chunk offsets, GPU index
 *      (or multi-slab global dropped_bitmap for F15).
 *   3. gecko_cuda_upload_query_async() / gecko_cuda_query_and_extend_uploaded() /
 *      gecko_cuda_free_uploaded_query() — overlapped non-blocking seqX H2D + run-bitmap fill.
 *   4. gecko_cuda_query_and_extend() — dynamically sizes aliased tile scratch (<= 608 MiB max /
 *      19 MiB small genomes), runs adaptive hybrid/L2-part/monolithic tiled query +
 *      DoubleBuffer SortKeys + zero-stack segmented warp extension + GECKO_PROFILE=1 events.
 *   5. gecko_cuda_free_ref() / gecko_cuda_free_frags()
 */

#include <algorithm>
#include <chrono>
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
static constexpr size_t PINNED_MIN_BYTES  = 64ULL << 10;     /* stage stream-0 transfers >= 64 KiB */

/* Reusable ping-pong pinned host staging buffer shared by GeckoCudaScratch and GeckoCudaUploadedQuery */
struct PinnedStager {
    uint8_t *h_pinned = nullptr;
    size_t slot_cap = 0;
    cudaEvent_t slot_done[2] = {nullptr, nullptr};
    bool slot_in_use[2] = {false, false};
    int next_slot = 0;
    bool failed = false;

    int ensure(size_t needed) {
        if (failed) return -1;
        size_t want_slot = std::min<size_t>(std::max<size_t>(needed, 1ULL), PINNED_SLOT_BYTES);
        if (want_slot <= slot_cap && h_pinned) return 0;

        if (sync_all() != 0) {
            failed = true;
            return -1;
        }
        if (h_pinned) {
            cudaFreeHost(h_pinned);
            h_pinned = nullptr;
            slot_cap = 0;
        }
        size_t alloc_slot = (want_slot + (1ULL << 20) - 1ULL) & ~((1ULL << 20) - 1ULL);
        if (alloc_slot > PINNED_SLOT_BYTES) alloc_slot = PINNED_SLOT_BYTES;
        if (cudaHostAlloc(&h_pinned, alloc_slot * 2, cudaHostAllocDefault) == cudaSuccess) {
            if (!slot_done[0]) {
                if (cudaEventCreateWithFlags(&slot_done[0], cudaEventDisableTiming) != cudaSuccess ||
                    cudaEventCreateWithFlags(&slot_done[1], cudaEventDisableTiming) != cudaSuccess) {
                    if (slot_done[0]) { cudaEventDestroy(slot_done[0]); slot_done[0] = nullptr; }
                    if (slot_done[1]) { cudaEventDestroy(slot_done[1]); slot_done[1] = nullptr; }
                    cudaFreeHost(h_pinned);
                    h_pinned = nullptr;
                    (void)cudaGetLastError();
                    failed = true;
                    return -1;
                }
            }
            slot_cap = alloc_slot;
            return 0;
        }
        (void)cudaGetLastError();
        h_pinned = nullptr;
        slot_cap = 0;
        failed = true;
        return -1;
    }

    int sync_all() {
        for (int s = 0; s < 2; s++) {
            if (slot_in_use[s]) {
                CUDA_CHECK(cudaEventSynchronize(slot_done[s]));
                slot_in_use[s] = false;
            }
        }
        return 0;
    }

    void destroy() {
        if (h_pinned) { cudaFreeHost(h_pinned); h_pinned = nullptr; }
        if (slot_done[0]) { cudaEventDestroy(slot_done[0]); slot_done[0] = nullptr; }
        if (slot_done[1]) { cudaEventDestroy(slot_done[1]); slot_done[1] = nullptr; }
        slot_cap = 0;
        slot_in_use[0] = false;
        slot_in_use[1] = false;
        next_slot = 0;
        failed = false;
    }
};

struct GeckoCudaScratch {
    uint64_t *d_hits = nullptr;
    uint64_t *d_hits_sorted = nullptr;
    GpuFrag *d_out_frags = nullptr;
    uint32_t *d_hit_count = nullptr;
    uint32_t *d_long_count = nullptr;
    uint32_t *d_frag_count = nullptr;
    uint32_t *d_stats_buf = nullptr;  /* F16: 2 u32 telemetry counters (valid_q, bypassed_q) */
    uint32_t hits_cap = 0U;
    uint32_t frags_cap = 0U;
    void *d_cub_temp = nullptr;
    size_t cub_temp_cap = 0;
    PinnedStager stager;
    cudaStream_t h2d_stream = nullptr;
    cudaEvent_t seqx_ready_ev = nullptr;
    bool allocated = false;
};

struct GeckoCudaUploadedQuery {
    GpuSeqBuf seqX{};
    PinnedStager stager;
    cudaStream_t stream = nullptr;
    cudaEvent_t ready_ev = nullptr;
    /* Optional profiling events (NULL when GECKO_PROFILE is unset) */
    cudaEvent_t prof_start = nullptr;
    cudaEvent_t prof_stop = nullptr;
    cudaEvent_t prof_bm_start = nullptr;
    cudaEvent_t prof_bm_stop = nullptr;
    bool prof_recorded = false;
    bool prof_bm_recorded = false;
};

/* Stream-aware staged H2D upload.
 * CRITICAL HOST-LIFETIME GUARANTEE: On async stream != 0, ALWAYS copy even < 64 KiB
 * buffers through stager->h_pinned (or sync stream if stager failed) so caller may
 * safely free/drop host `src` immediately upon return! */
static int upload_h2d_staged_stream(
    void *dst,
    const void *src,
    size_t bytes,
    PinnedStager *stager,
    cudaStream_t stream
) {
    if (bytes == 0) return 0;

    bool must_stage = (stream != nullptr) || (bytes >= PINNED_MIN_BYTES);
    if (!must_stage || !stager || stager->ensure(bytes) != 0) {
        if (stream == nullptr) {
            CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice));
        } else {
            CUDA_CHECK(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, stream));
            /* Sync so host pageable `src` is never read after caller frees it */
            CUDA_CHECK(cudaStreamSynchronize(stream));
        }
        return 0;
    }

    const uint8_t *s_ptr = static_cast<const uint8_t *>(src);
    uint8_t *d_ptr = static_cast<uint8_t *>(dst);
    size_t slot_cap = stager->slot_cap;
    size_t off = 0;
    while (off < bytes) {
        size_t n = std::min<size_t>(slot_cap, bytes - off);
        int s = stager->next_slot;
        stager->next_slot ^= 1;
        if (stager->slot_in_use[s]) {
            CUDA_CHECK(cudaEventSynchronize(stager->slot_done[s]));
            stager->slot_in_use[s] = false;
        }
        uint8_t *buf = stager->h_pinned + (size_t)s * slot_cap;
        memcpy(buf, s_ptr + off, n);
        CUDA_CHECK(cudaMemcpyAsync(d_ptr + off, buf, n, cudaMemcpyHostToDevice, stream));
        CUDA_CHECK(cudaEventRecord(stager->slot_done[s], stream));
        stager->slot_in_use[s] = true;
        off += n;
    }
    return 0;
}

static int download_d2h_staged(void *dst, const void *src, size_t bytes, PinnedStager *stager) {
    if (bytes == 0) return 0;
    if (bytes < PINNED_MIN_BYTES || !stager || stager->ensure(bytes) != 0) {
        CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost));
        return 0;
    }
    /* Ensure any previous async H2D using stager slots has finished before D2H reuses h_pinned */
    if (stager->sync_all() != 0) return -1;

    const uint8_t *s_ptr = static_cast<const uint8_t *>(src);
    uint8_t *d_ptr = static_cast<uint8_t *>(dst);
    size_t slot_cap = stager->slot_cap;

    size_t cur_off = 0;
    size_t cur_n = std::min<size_t>(slot_cap, bytes);
    int s = 0;
    CUDA_CHECK(cudaMemcpyAsync(stager->h_pinned, s_ptr, cur_n, cudaMemcpyDeviceToHost, 0));
    CUDA_CHECK(cudaEventRecord(stager->slot_done[0], 0));

    while (cur_off < bytes) {
        size_t next_off = cur_off + cur_n;
        size_t next_n = (next_off < bytes) ? std::min<size_t>(slot_cap, bytes - next_off) : 0;
        if (next_n > 0) {
            uint8_t *next_buf = stager->h_pinned + (size_t)(s ^ 1) * slot_cap;
            CUDA_CHECK(cudaMemcpyAsync(next_buf, s_ptr + next_off, next_n, cudaMemcpyDeviceToHost, 0));
            CUDA_CHECK(cudaEventRecord(stager->slot_done[s ^ 1], 0));
        }
        CUDA_CHECK(cudaEventSynchronize(stager->slot_done[s]));
        memcpy(d_ptr + cur_off, stager->h_pinned + (size_t)s * slot_cap, cur_n);
        cur_off = next_off;
        cur_n = next_n;
        s ^= 1;
    }
    return 0;
}

/* W3-F12a + W3-F12b: Stream-ordered PackedSeq upload + GPU fill_runs_bitmap materialization.
 * Uses ONLY cudaMemcpyAsync + cudaMemsetAsync + stream-ordered fill_runs_bitmap_kernel. */
static int upload_packed_seq_stream(
    const GeckoPackedSeqHost *src,
    GpuSeqBuf &dst,
    PinnedStager *stager,
    cudaStream_t stream,
    bool wait_host_done,
    cudaEvent_t bm_ev_start = nullptr,
    cudaEvent_t bm_ev_stop = nullptr,
    bool *bm_recorded = nullptr
) {
    memset(&dst, 0, sizeof(dst));
    dst.len = src->len;
    dst.n_seps = src->n_seps;
    dst.n_iupac = src->n_iupac;

    size_t bits_words = (size_t)((src->len + 31ULL) / 32ULL) + 2;
    size_t copy_words = (src->bits_words < bits_words) ? (size_t)src->bits_words : bits_words;
    CUDA_CHECK(cudaMalloc(&dst.d_bits, bits_words * sizeof(uint64_t)));
    if (copy_words > 0 && src->bits) {
        if (upload_h2d_staged_stream(dst.d_bits, src->bits, copy_words * sizeof(uint64_t), stager, stream) != 0) {
            return -1;
        }
    }
    if (bits_words > copy_words) {
        CUDA_CHECK(cudaMemsetAsync(dst.d_bits + copy_words, 0, (bits_words - copy_words) * sizeof(uint64_t), stream));
    }

    /* W3-F12b: 4-case (has_special, has_soft) aliasing matrix materialized on GPU */
    bool has_special = (src->len > 0ULL && src->n_special_runs > 0U && src->special_runs != nullptr);
    bool has_soft    = (src->len > 0ULL && src->soft_mask != 0U && src->n_soft_runs > 0U && src->soft_runs != nullptr);
    size_t bm_words  = (size_t)((src->len + 31ULL) / 32ULL) + 2;

    if (has_special) {
        CUDA_CHECK(cudaMalloc(&dst.d_special_runs, (size_t)src->n_special_runs * sizeof(GeckoRun)));
        if (upload_h2d_staged_stream(
                dst.d_special_runs, src->special_runs,
                (size_t)src->n_special_runs * sizeof(GeckoRun), stager, stream) != 0) {
            return -1;
        }
    }
    if (has_soft) {
        CUDA_CHECK(cudaMalloc(&dst.d_soft_runs, (size_t)src->n_soft_runs * sizeof(GeckoRun)));
        if (upload_h2d_staged_stream(
                dst.d_soft_runs, src->soft_runs,
                (size_t)src->n_soft_runs * sizeof(GeckoRun), stager, stream) != 0) {
            return -1;
        }
    }

    if (has_special || has_soft) {
        if (bm_ev_start) (void)cudaEventRecord(bm_ev_start, stream);
        if (has_special && !has_soft) {
            /* Case (true, false): d_invalid_bm aliases d_special_bm */
            CUDA_CHECK(cudaMalloc(&dst.d_special_bm, bm_words * sizeof(uint32_t)));
            CUDA_CHECK(cudaMemsetAsync(dst.d_special_bm, 0, bm_words * sizeof(uint32_t), stream));
            if (gpu_fill_runs_bitmap(dst.d_special_runs, src->n_special_runs, src->len, dst.d_special_bm, stream) != 0) {
                return -1;
            }
            dst.d_invalid_bm = dst.d_special_bm;
            dst.owns_invalid_bm = false;
        } else if (!has_special && has_soft) {
            /* Case (false, true): d_special_bm = nullptr; d_invalid_bm from soft_runs */
            CUDA_CHECK(cudaMalloc(&dst.d_invalid_bm, bm_words * sizeof(uint32_t)));
            dst.owns_invalid_bm = true;
            CUDA_CHECK(cudaMemsetAsync(dst.d_invalid_bm, 0, bm_words * sizeof(uint32_t), stream));
            if (gpu_fill_runs_bitmap(dst.d_soft_runs, src->n_soft_runs, src->len, dst.d_invalid_bm, stream) != 0) {
                return -1;
            }
        } else {
            /* Case (true, true): separate d_special_bm and d_invalid_bm = special | soft */
            CUDA_CHECK(cudaMalloc(&dst.d_special_bm, bm_words * sizeof(uint32_t)));
            CUDA_CHECK(cudaMemsetAsync(dst.d_special_bm, 0, bm_words * sizeof(uint32_t), stream));
            if (gpu_fill_runs_bitmap(dst.d_special_runs, src->n_special_runs, src->len, dst.d_special_bm, stream) != 0) {
                return -1;
            }

            CUDA_CHECK(cudaMalloc(&dst.d_invalid_bm, bm_words * sizeof(uint32_t)));
            dst.owns_invalid_bm = true;
            CUDA_CHECK(cudaMemsetAsync(dst.d_invalid_bm, 0, bm_words * sizeof(uint32_t), stream));
            if (gpu_fill_runs_bitmap(dst.d_special_runs, src->n_special_runs, src->len, dst.d_invalid_bm, stream) != 0) {
                return -1;
            }
            if (gpu_fill_runs_bitmap(dst.d_soft_runs, src->n_soft_runs, src->len, dst.d_invalid_bm, stream) != 0) {
                return -1;
            }
        }
        if (bm_ev_stop) {
            (void)cudaEventRecord(bm_ev_stop, stream);
            if (bm_recorded) *bm_recorded = true;
        }
    }

    if (src->n_seps > 0 && src->sep_pos) {
        CUDA_CHECK(cudaMalloc(&dst.d_sep_pos, (size_t)src->n_seps * sizeof(uint32_t)));
        if (upload_h2d_staged_stream(
                dst.d_sep_pos, src->sep_pos,
                (size_t)src->n_seps * sizeof(uint32_t), stager, stream) != 0) {
            return -1;
        }
    }
    if (src->n_iupac > 0 && src->iupac_pos && src->iupac_byte) {
        CUDA_CHECK(cudaMalloc(&dst.d_iupac_pos, (size_t)src->n_iupac * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&dst.d_iupac_byte, (size_t)src->n_iupac * sizeof(uint8_t)));
        if (upload_h2d_staged_stream(
                dst.d_iupac_pos, src->iupac_pos,
                (size_t)src->n_iupac * sizeof(uint32_t), stager, stream) != 0) {
            return -1;
        }
        if (upload_h2d_staged_stream(
                dst.d_iupac_byte, src->iupac_byte,
                (size_t)src->n_iupac * sizeof(uint8_t), stager, stream) != 0) {
            return -1;
        }
    }

    if (wait_host_done && stager) {
        if (stager->sync_all() != 0) return -1;
    }
    return 0;
}

static void free_gpu_seq(GpuSeqBuf &s) {
    if (s.d_bits) cudaFree(s.d_bits);
    if (s.d_special_bm) cudaFree(s.d_special_bm);
    if (s.owns_invalid_bm && s.d_invalid_bm) cudaFree(s.d_invalid_bm);
    if (s.d_special_runs) cudaFree(s.d_special_runs);
    if (s.d_soft_runs) cudaFree(s.d_soft_runs);
    if (s.d_sep_pos) cudaFree(s.d_sep_pos);
    if (s.d_iupac_pos) cudaFree(s.d_iupac_pos);
    if (s.d_iupac_byte) cudaFree(s.d_iupac_byte);
    memset(&s, 0, sizeof(s));
}

/* F15 + W3-F11b: Dynamic small-genome scratch sizing (19 MiB for 4 kb TYLCV, <= 608 MiB max:
 * d_hits 256 + d_hits_sorted 256 + d_out_frags 96; d_hit_flags & d_filt_count eliminated!) */
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
    if (!sc.h2d_stream) {
        CUDA_CHECK(cudaStreamCreateWithFlags(&sc.h2d_stream, cudaStreamNonBlocking));
    }
    if (!sc.seqx_ready_ev) {
        CUDA_CHECK(cudaEventCreateWithFlags(&sc.seqx_ready_ev, cudaEventDisableTiming));
    }
    if (sc.allocated && sc.hits_cap >= want_hits && sc.frags_cap >= want_frags) return 0;

    if (sc.d_hits)        { cudaFree(sc.d_hits);        sc.d_hits = nullptr; }
    if (sc.d_hits_sorted) { cudaFree(sc.d_hits_sorted); sc.d_hits_sorted = nullptr; }
    if (sc.d_out_frags)   { cudaFree(sc.d_out_frags);   sc.d_out_frags = nullptr; }
    sc.hits_cap = 0U;
    sc.frags_cap = 0U;
    sc.allocated = false;

    CUDA_CHECK(cudaMalloc(&sc.d_hits, (size_t)want_hits * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_hits_sorted, (size_t)want_hits * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_out_frags, (size_t)want_frags * sizeof(GpuFrag)));
    if (!sc.d_hit_count)  CUDA_CHECK(cudaMalloc(&sc.d_hit_count, sizeof(uint32_t)));
    if (!sc.d_long_count) CUDA_CHECK(cudaMalloc(&sc.d_long_count, sizeof(uint32_t)));
    if (!sc.d_frag_count) CUDA_CHECK(cudaMalloc(&sc.d_frag_count, sizeof(uint32_t)));
    if (!sc.d_stats_buf)  CUDA_CHECK(cudaMalloc(&sc.d_stats_buf, 2 * sizeof(uint32_t)));
    sc.hits_cap = want_hits;
    sc.frags_cap = want_frags;
    sc.allocated = true;
    return 0;
}

static void free_scratch(GeckoCudaScratch &sc) {
    if (sc.h2d_stream) {
        (void)cudaStreamSynchronize(sc.h2d_stream);
        (void)cudaStreamDestroy(sc.h2d_stream);
        sc.h2d_stream = nullptr;
    }
    if (sc.seqx_ready_ev) { (void)cudaEventDestroy(sc.seqx_ready_ev); sc.seqx_ready_ev = nullptr; }
    if (sc.d_hits) (void)cudaFree(sc.d_hits);
    if (sc.d_hits_sorted) (void)cudaFree(sc.d_hits_sorted);
    if (sc.d_out_frags) (void)cudaFree(sc.d_out_frags);
    if (sc.d_hit_count) (void)cudaFree(sc.d_hit_count);
    if (sc.d_long_count) (void)cudaFree(sc.d_long_count);
    if (sc.d_frag_count) (void)cudaFree(sc.d_frag_count);
    if (sc.d_stats_buf) (void)cudaFree(sc.d_stats_buf);
    if (sc.d_cub_temp) (void)cudaFree(sc.d_cub_temp);
    sc.stager.destroy();
    sc = GeckoCudaScratch{};
}

/* W3-F15a: Persistent CudaScratchPool (immortal construct-on-first-use singleton to avoid
 * C++ static destructor crashes after libcudart teardown at process exit). */
static inline uint64_t scratch_device_bytes(const GeckoCudaScratch &sc) {
    uint64_t bytes = (uint64_t)sc.cub_temp_cap;
    if (sc.allocated) {
        bytes += (uint64_t)sc.hits_cap * sizeof(uint64_t) * 2ULL;
        bytes += (uint64_t)sc.frags_cap * sizeof(GpuFrag);
    }
    return bytes;
}

static inline uint64_t scratch_total_cached_bytes(const GeckoCudaScratch &sc) {
    return scratch_device_bytes(sc) + (uint64_t)sc.stager.slot_cap * 2ULL;
}

static inline bool parse_scratch_pool_enabled_env() {
    const char *env = std::getenv("GECKO_CUDA_SCRATCH_POOL");
    if (!env) return true;
    return gecko_parse_profile_str(env);
}

using GpuScratchPool = GeckoScratchPoolCore<
    GeckoCudaScratch,
    uint64_t (*)(const GeckoCudaScratch &),
    void (*)(GeckoCudaScratch &)>;

static GpuScratchPool &get_scratch_pool() {
    static GpuScratchPool *pool = new GpuScratchPool(
        &scratch_total_cached_bytes,
        &free_scratch,
        parse_scratch_pool_enabled_env(),
        /*initial_max_idle=*/1);
    return *pool;
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
    bool is_minimizer;
    bool scratch_healthy;
    uint32_t *d_global_dropped_bm; /* F15: built once for multi-slab contiguous freq_limit > 0 */
    GpuRefIndex idx;
    GpuMinIndex min_idx;
    GeckoCudaScratch scratch;
    GeckoCudaProfile prof;
};

static int parse_query_mode_env() {
    const char *env = std::getenv("GECKO_CUDA_QUERY_MODE");
    if (!env || env[0] == '\0') return GECKO_QMODE_AUTO;
    if (strcmp(env, "monolithic") == 0) return GECKO_QMODE_MONOLITHIC;
    if (strcmp(env, "l2part") == 0) return GECKO_QMODE_L2_PART;
    if (strcmp(env, "hybrid") == 0) return GECKO_QMODE_HYBRID;
    return GECKO_QMODE_AUTO;
}

static const char *query_mode_name(int mode, bool is_direct) {
    if (is_direct) return "direct";
    switch (mode) {
        case GECKO_QMODE_MONOLITHIC: return "monolithic";
        case GECKO_QMODE_L2_PART:    return "l2part";
        case GECKO_QMODE_HYBRID:     return "hybrid";
        default:                     return "auto";
    }
}

static uint32_t choose_cuda_slab_count(uint64_t total_valid_kmers, uint64_t owned_scratch_bytes) {
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
            /* Reserve VRAM for query sequence, aliased tile scratch (<= 608 MiB max:
             * d_hits 256 + d_hits_sorted 256 + d_out_frags 96), CUB temp, and bitmaps,
             * avoiding double-counting already-owned pooled scratch bytes. */
            size_t base_reserve = std::max<size_t>(896ULL << 20, free_vram / 5ULL);
            size_t reserve = (base_reserve > owned_scratch_bytes + (256ULL << 20))
                ? (base_reserve - (size_t)owned_scratch_bytes)
                : (256ULL << 20);
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
        (void)cudaGetLastError();
        return -1;
    }
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess) {
        (void)cudaGetLastError();
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

/* W3-F12a: Warm up device 0 primary CUDA context on background thread during ref FASTA load */
extern "C" int gecko_cuda_warmup(void) {
    int dev_count = 0;
    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) {
        (void)cudaGetLastError();
        return -1;
    }
    if (cudaSetDevice(0) != cudaSuccess) {
        (void)cudaGetLastError();
        return -1;
    }
    if (cudaFree(0) != cudaSuccess) {
        (void)cudaGetLastError();
        return -1;
    }
    return 0;
}

static void destroy_uploaded_query(GeckoCudaUploadedQuery *uq) {
    if (!uq) return;
    if (uq->stream) {
        (void)cudaStreamSynchronize(uq->stream);
    }
    free_gpu_seq(uq->seqX);
    uq->stager.destroy();
    if (uq->ready_ev)      { cudaEventDestroy(uq->ready_ev);      uq->ready_ev = nullptr; }
    if (uq->prof_start)    { cudaEventDestroy(uq->prof_start);    uq->prof_start = nullptr; }
    if (uq->prof_stop)     { cudaEventDestroy(uq->prof_stop);     uq->prof_stop = nullptr; }
    if (uq->prof_bm_start) { cudaEventDestroy(uq->prof_bm_start); uq->prof_bm_start = nullptr; }
    if (uq->prof_bm_stop)  { cudaEventDestroy(uq->prof_bm_stop);  uq->prof_bm_stop = nullptr; }
    if (uq->stream)        { cudaStreamDestroy(uq->stream);       uq->stream = nullptr; }
    delete uq;
}

/* W3-F12a: Async seqX H2D upload + GPU fill_runs_bitmap on dedicated non-blocking stream.
 * Graceful OOM fallback: on any CUDA error, cleans up cleanly and returns nullptr so caller
 * falls back to uploading seqX after prepare_ref releases temporary SortPairs buffers. */
extern "C" GeckoCudaUploadedQuery *gecko_cuda_upload_query_async(
    const GeckoPackedSeqHost *query_seq
) {
    if (!query_seq) return nullptr;
    if (cudaSetDevice(0) != cudaSuccess) {
        (void)cudaGetLastError();
        return nullptr;
    }

    auto *uq = new GeckoCudaUploadedQuery{};
    if (cudaStreamCreateWithFlags(&uq->stream, cudaStreamNonBlocking) != cudaSuccess ||
        cudaEventCreateWithFlags(&uq->ready_ev, cudaEventDisableTiming) != cudaSuccess) {
        (void)cudaGetLastError();
        destroy_uploaded_query(uq);
        return nullptr;
    }

    if (gecko_cuda_profile_enabled()) {
        if (cudaEventCreate(&uq->prof_start) == cudaSuccess &&
            cudaEventCreate(&uq->prof_stop) == cudaSuccess &&
            cudaEventCreate(&uq->prof_bm_start) == cudaSuccess &&
            cudaEventCreate(&uq->prof_bm_stop) == cudaSuccess) {
            (void)cudaEventRecord(uq->prof_start, uq->stream);
        } else {
            (void)cudaGetLastError();
            if (uq->prof_start)    { cudaEventDestroy(uq->prof_start);    uq->prof_start = nullptr; }
            if (uq->prof_stop)     { cudaEventDestroy(uq->prof_stop);     uq->prof_stop = nullptr; }
            if (uq->prof_bm_start) { cudaEventDestroy(uq->prof_bm_start); uq->prof_bm_start = nullptr; }
            if (uq->prof_bm_stop)  { cudaEventDestroy(uq->prof_bm_stop);  uq->prof_bm_stop = nullptr; }
        }
    }

    if (upload_packed_seq_stream(
            query_seq, uq->seqX, &uq->stager, uq->stream,
            /*wait_host_done=*/false,
            uq->prof_bm_start, uq->prof_bm_stop, &uq->prof_bm_recorded) != 0) {
        (void)cudaGetLastError();
        destroy_uploaded_query(uq);
        return nullptr;
    }

    if (uq->prof_stop) {
        if (cudaEventRecord(uq->prof_stop, uq->stream) == cudaSuccess) {
            uq->prof_recorded = true;
        }
    }
    if (cudaEventRecord(uq->ready_ev, uq->stream) != cudaSuccess) {
        (void)cudaGetLastError();
        destroy_uploaded_query(uq);
        return nullptr;
    }
    return uq;
}

extern "C" void gecko_cuda_free_uploaded_query(GeckoCudaUploadedQuery *uq) {
    destroy_uploaded_query(uq);
}

extern "C" GeckoCudaRefHandle *gecko_cuda_prepare_ref(
    const GeckoPackedSeqHost *ref_seq,
    const GeckoCudaConfig *cfg
) {
    if (!ref_seq || !cfg) return nullptr;
    int dev_count = 0;
    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) {
        (void)cudaGetLastError();
        return nullptr;
    }
    if (cudaSetDevice(0) != cudaSuccess) {
        (void)cudaGetLastError();
        return nullptr;
    }

    int wl = cfg->wl;
    uint64_t seed_mask = cfg->seed_mask;
    uint32_t freq_limit = cfg->freq_limit;

    auto *h = new GeckoCudaRefHandle{};
    h->prof.init(gecko_cuda_profile_enabled());
    std::chrono::steady_clock::time_point t_prep_start{};
    if (h->prof.enabled) {
        t_prep_start = std::chrono::steady_clock::now();
    }

    /* W3-F15a: Automatically check out persistent scratch workspace from CudaScratchPool */
    bool pool_hit = get_scratch_pool().checkout(h->scratch);
    h->scratch_healthy = true;
    h->prof.scratch_pool_hit = pool_hit;

    /* Evict any OTHER idle pool entries if free VRAM is below 2 GiB */
    if (get_scratch_pool().idle_count() > 0) {
        size_t free_vram = 0, total_vram = 0;
        if (cudaMemGetInfo(&free_vram, &total_vram) == cudaSuccess && free_vram < (2ULL << 30)) {
            get_scratch_pool().release_all();
        }
    }

    h->wl = wl;
    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    h->seed_mask = (seed_mask != 0ULL) ? seed_mask : span_mask;
    h->freq_limit = freq_limit;
    h->num_slabs = 1U;
    h->is_spaced = (h->seed_mask != span_mask);
    h->soft_mask = (cfg->soft_mask != 0U);
    h->compare_reverse = (cfg->compare_reverse != 0U);
    h->is_minimizer = false;
    h->d_global_dropped_bm = nullptr;

    /* F16: Query hardware L2 cache size once for L2-resident partitioned probe auto-tuning */
    int dev = 0;
    cudaGetDevice(&dev);
    int l2_bytes = 0;
    h->l2_cache_bytes = 4 << 20;
    if (cudaDeviceGetAttribute(&l2_bytes, cudaDevAttrL2CacheSize, dev) == cudaSuccess && l2_bytes > 0) {
        h->l2_cache_bytes = l2_bytes;
    }

    bool y_bm_recorded = false;
    if (h->prof.enabled) {
        (void)cudaEventRecord(h->prof.h2d_start, 0);
    }
    if (upload_packed_seq_stream(
            ref_seq, h->seqY, &h->scratch.stager, /*stream=*/nullptr,
            /*wait_host_done=*/true,
            h->prof.enabled ? h->prof.h2d_bm_start : nullptr,
            h->prof.enabled ? h->prof.h2d_bm_stop : nullptr,
            &y_bm_recorded) != 0) {
        /* W3-M3 Fix 3: Drain stream 0 before freeing pinned stager & device buffers */
        (void)cudaStreamSynchronize(0);
        h->scratch_healthy = false;
        h->prof.destroy();
        free_scratch(h->scratch);
        free_gpu_seq(h->seqY);
        delete h;
        return nullptr;
    }
    if (h->prof.enabled) {
        float ms = 0.0f;
        if (cudaEventRecord(h->prof.h2d_stop, 0) == cudaSuccess &&
            cudaEventSynchronize(h->prof.h2d_stop) == cudaSuccess &&
            cudaEventElapsedTime(&ms, h->prof.h2d_start, h->prof.h2d_stop) == cudaSuccess) {
            h->prof.upload_y_ms = (double)ms;
        }
        if (y_bm_recorded &&
            cudaEventElapsedTime(&ms, h->prof.h2d_bm_start, h->prof.h2d_bm_stop) == cudaSuccess) {
            h->prof.fill_bm_y_ms = (double)ms;
        }
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
            ref_seq->len,
            &h->prof);

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

    /* Extreme VRAM pressure guard: release large pooled tile buffers before index build if < 512 MiB free */
    if (h->scratch.allocated) {
        size_t free_vram = 0, total_vram = 0;
        if (cudaMemGetInfo(&free_vram, &total_vram) == cudaSuccess && free_vram < (512ULL << 20)) {
            if (h->scratch.d_hits)        { (void)cudaFree(h->scratch.d_hits);        h->scratch.d_hits = nullptr; }
            if (h->scratch.d_hits_sorted) { (void)cudaFree(h->scratch.d_hits_sorted); h->scratch.d_hits_sorted = nullptr; }
            if (h->scratch.d_out_frags)   { (void)cudaFree(h->scratch.d_out_frags);   h->scratch.d_out_frags = nullptr; }
            h->scratch.hits_cap = 0U;
            h->scratch.frags_cap = 0U;
            h->scratch.allocated = false;
        }
    }
    uint64_t owned_sc_bytes = scratch_device_bytes(h->scratch);

    /* W3-F14: Try 100% lossless GPU Canonical Minimizer Super-Seeding (w/ HF exact-k overflow CSR)
     * when eligible; cleanly fall back to exact/multi-slab if VRAM or GECKO_CUDA_MAX_SLAB_KMERS rejects it. */
    uint64_t total_valid = 0ULL;
    bool want_min = gecko_cuda_should_use_minimizer(
                        cfg->seeding_mode, wl, h->is_spaced, cfg->two_hit_window, freq_limit)
                    && ref_seq->len >= (uint64_t)wl;
    if (want_min) {
        if (gpu_build_min_index(
                h->seqY.view(), wl, freq_limit, owned_sc_bytes,
                h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                h->min_idx, &h->prof)) {
            h->is_minimizer = true;
            h->num_slabs = 1U;
            total_valid = h->min_idx.n_raw_min;
            h->prof.is_minimizer = true;
            h->prof.min_kp = h->min_idx.kp;
            h->prof.min_w = h->min_idx.w;
            h->prof.min_ref_occs = h->min_idx.n_raw_min;
            h->prof.min_ovf_buckets = h->min_idx.ovf_routed_keys;
            h->prof.min_ovf_kmers = h->min_idx.ovf_n_kept;
        }
    }

    if (!h->is_minimizer) {
        /* F15: Single count + ExclusiveSum pass reused directly in single-slab gpu_build_ref_index */
        uint32_t *d_chunk_offsets = nullptr;
        uint32_t num_chunks = 0U;
        total_valid = gpu_count_valid_kmers(
            h->seqY.view(), wl, h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
            d_chunk_offsets, num_chunks, &h->prof);

        h->num_slabs = choose_cuda_slab_count(total_valid, owned_sc_bytes);
        if (freq_limit == 0U) {
            h->num_slabs = 1U;
        }
        if (h->num_slabs == 1U) {
            h->idx = gpu_build_ref_index(
                h->seqY.view(), wl, h->seed_mask, freq_limit,
                h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                0U, 256U, d_chunk_offsets, (uint32_t)total_valid, true, &h->prof);
            if (d_chunk_offsets) (void)cudaFree(d_chunk_offsets);
        } else {
            if (d_chunk_offsets) (void)cudaFree(d_chunk_offsets);
            /* F15: Build global 1-bit dropped bitmap once for multi-slab contiguous when drops are possible */
            if (!h->is_spaced && freq_limit > 0U && (uint64_t)freq_limit < total_valid) {
                h->d_global_dropped_bm = gpu_build_global_dropped_bitmap(
                    h->seqY.view(), wl, h->seed_mask, freq_limit, h->num_slabs,
                    h->scratch.d_cub_temp, h->scratch.cub_temp_cap, &h->prof);
            }
        }
    }

    /* Free transient seqY GPU run arrays now that stream 0 has synced inside index build */
    if (h->seqY.d_special_runs) { (void)cudaFree(h->seqY.d_special_runs); h->seqY.d_special_runs = nullptr; }
    if (h->seqY.d_soft_runs)    { (void)cudaFree(h->seqY.d_soft_runs);    h->seqY.d_soft_runs = nullptr; }

    if (h->prof.enabled) {
        double total_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t_prep_start).count();
        uint32_t rep_valid = h->is_minimizer
            ? h->min_idx.n_raw_min
            : ((h->num_slabs == 1U) ? h->idx.n_valid : (uint32_t)total_valid);
        uint32_t rep_kept = h->is_minimizer
            ? (h->min_idx.n_kept_main + h->min_idx.ovf_n_kept)
            : ((h->num_slabs == 1U) ? h->idx.n_kept : 0U);
        int rep_bbits = h->is_minimizer
            ? h->min_idx.bucket_bits
            : ((h->num_slabs == 1U) ? h->idx.bucket_bits : 0);
        int rep_direct = (!h->is_minimizer && h->num_slabs == 1U && h->idx.is_direct_csr) ? 1 : 0;
        fprintf(
            stderr,
            "[gecko-cuda-profile] prepare_gpu total_ms=%.2f upload_y_ms=%.2f fill_bm_y_ms=%.2f "
            "revcomp_ms=%.2f count_ms=%.2f extract_ms=%.2f sort_pairs_ms=%.2f freq_filter_ms=%.2f "
            "sel_vals_ms=%.2f sel_tags_ms=%.2f bucket_ms=%.2f global_drop_ms=%.2f "
            "slabs=%u valid=%u kept=%u bucket_bits=%d direct=%d l2_kib=%d "
            "seeding=%s min_kp=%d min_w=%d min_ref_occs=%u min_ovf_buckets=%u min_ovf_kmers=%u "
            "min_extract_ms=%.2f min_hf_ms=%.2f min_ovf_ms=%.2f scratch_pool=%s\n",
            total_ms, h->prof.upload_y_ms, h->prof.fill_bm_y_ms,
            h->prof.revcomp_ms, h->prof.count_valid_ms, h->prof.extract_kmers_ms,
            h->prof.sort_pairs_ms, h->prof.freq_filter_ms,
            h->prof.select_vals_ms, h->prof.select_tags_ms, h->prof.bucket_build_ms,
            h->prof.global_dropped_ms,
            h->num_slabs, rep_valid, rep_kept, rep_bbits, rep_direct,
            h->l2_cache_bytes / 1024,
            h->prof.is_minimizer ? "minimizer" : "exact",
            h->prof.min_kp, h->prof.min_w, h->prof.min_ref_occs,
            h->prof.min_ovf_buckets, h->prof.min_ovf_kmers,
            h->prof.min_extract_ms, h->prof.min_hf_ms, h->prof.min_ovf_ms,
            h->prof.scratch_pool_hit ? "hit" : "miss");
    }

    return h;
}

extern "C" int gecko_cuda_ref_is_minimizer(const GeckoCudaRefHandle *h) {
    return (h && h->is_minimizer) ? 1 : 0;
}

static int run_query_and_extend_impl(
    GeckoCudaRefHandle *h,
    const GeckoPackedSeqHost *query_seq,
    GeckoCudaUploadedQuery *uq,
    const GeckoCudaConfig *cfg,
    GeckoCudaResult *out
) {
    if (!h || (!query_seq && !uq) || !cfg || !out) return -1;
    memset(out, 0, sizeof(*out));

    uint32_t min_len = cfg->min_len;
    float sim_th = cfg->sim_th;
    int match_score = cfg->match_score;
    int mismatch_penalty = cfg->mismatch_penalty;
    int n_penalty = cfg->n_penalty;
    uint32_t two_hit_window = cfg->two_hit_window;

    if (h->is_minimizer && two_hit_window > 0U) {
        return -1;
    }

    if (cudaSetDevice(0) != cudaSuccess) {
        h->scratch_healthy = false;
        return -1;
    }

    h->prof.reset_query_counters();
    h->prof.two_hit_window = two_hit_window;
    std::chrono::steady_clock::time_point t_q_start{};
    if (h->prof.enabled) {
        t_q_start = std::chrono::steady_clock::now();
    }

    uint64_t lenX = uq ? uq->seqX.len : query_seq->len;

    {
        std::chrono::steady_clock::time_point t_sc{};
        if (h->prof.enabled) t_sc = std::chrono::steady_clock::now();
        if (ensure_scratch(h->scratch, lenX) != 0) {
            h->scratch_healthy = false;
            free_scratch(h->scratch);
            return -1;
        }
        if (h->prof.enabled) {
            h->prof.ensure_scratch_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t_sc).count();
        }
    }

    GpuSeqBuf local_seqX{};
    bool owns_local_seqX = (uq == nullptr);
    GpuRefIndex transient_idx{};
    bool has_transient_idx = false;
    cudaEvent_t seqx_ready_ev = nullptr;
    bool waited_seqx = false;

    /* W3-M3 Fix 4: Centralized error teardown freeing transient_idx & local_seqX and poisoning scratch */
    auto fail_cleanup = [&](const char *ctx, cudaError_t err) -> int {
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA error in %s: %s\n", ctx, cudaGetErrorString(err));
        }
        h->scratch_healthy = false;
        (void)cudaStreamSynchronize(0);
        if (h->scratch.h2d_stream) (void)cudaStreamSynchronize(h->scratch.h2d_stream);
        if (has_transient_idx) {
            gpu_free_ref_index(transient_idx);
            has_transient_idx = false;
        }
        if (owns_local_seqX) {
            free_gpu_seq(local_seqX);
            owns_local_seqX = false;
        }
        return -1;
    };

    if (owns_local_seqX) {
        if (h->prof.enabled) {
            (void)cudaEventRecord(h->prof.h2d_start, h->scratch.h2d_stream);
        }
        if (upload_packed_seq_stream(
                query_seq, local_seqX, &h->scratch.stager, h->scratch.h2d_stream,
                /*wait_host_done=*/false,
                h->prof.enabled ? h->prof.h2d_bm_start : nullptr,
                h->prof.enabled ? h->prof.h2d_bm_stop : nullptr,
                &h->prof.h2d_bm_recorded) != 0) {
            return fail_cleanup("upload_packed_seq_stream(seqX)", cudaSuccess);
        }
        if (h->prof.enabled) {
            if (cudaEventRecord(h->prof.h2d_stop, h->scratch.h2d_stream) == cudaSuccess) {
                h->prof.h2d_recorded = true;
            }
        }
        cudaError_t ev_err = cudaEventRecord(h->scratch.seqx_ready_ev, h->scratch.h2d_stream);
        if (ev_err != cudaSuccess) {
            return fail_cleanup("cudaEventRecord(seqx_ready_ev)", ev_err);
        }
        seqx_ready_ev = h->scratch.seqx_ready_ev;
    } else {
        seqx_ready_ev = uq->ready_ev;
        h->prof.upload_x_overlapped = true;
    }

    uint64_t *d_hits = h->scratch.d_hits;
    uint64_t *d_hits_sorted = h->scratch.d_hits_sorted;
    GpuFrag *d_out_frags = h->scratch.d_out_frags;
    uint32_t *d_hit_count = h->scratch.d_hit_count;
    uint32_t *d_long_count = h->scratch.d_long_count;
    uint32_t *d_frag_count = h->scratch.d_frag_count;
    uint32_t *d_stats_buf = h->scratch.d_stats_buf;
    uint32_t hits_cap = h->scratch.hits_cap;
    uint32_t frags_cap = h->scratch.frags_cap;

    GpuSeqView viewX = uq ? uq->seqX.view() : local_seqX.view();
    GpuSeqView viewY = h->seqY.view();
    GpuSeqView viewY_rc = h->compare_reverse ? h->seqY_rc.view() : viewY;

    std::vector<GpuFrag> all_frags;

    /* Choose an initial tile size tuned to seeding mode & k-mer density, and adaptively halve
     * any tile that overflows hits_cap or frags_cap so zero hits and zero fragments are dropped. */
    uint64_t init_tile_bp = GECKO_QUERY_TILE_BP;
    if (h->is_minimizer) {
        init_tile_bp = (h->wl >= 24) ? (64ULL << 20) : (32ULL << 20);
    } else {
        if (h->wl <= 14) {
            init_tile_bp = 1ULL * 1024ULL * 1024ULL;
        } else if (h->wl <= 19) {
            init_tile_bp = 4ULL * 1024ULL * 1024ULL;
        }
        if (!h->soft_mask && init_tile_bp > 1ULL * 1024ULL * 1024ULL) {
            init_tile_bp >>= 1;
        }
    }

    /* W3-F15b: Eliminate tile-boundary splits of gapped --two-hit pairs via minimal tile overlap */
    uint64_t th_overlap = 0ULL;
    if (two_hit_window > 0U) {
        th_overlap = std::min<uint64_t>((uint64_t)h->wl + (uint64_t)two_hit_window, init_tile_bp / 4ULL);
    }
    uint64_t tile_step = (init_tile_bp > th_overlap) ? (init_tile_bp - th_overlap) : init_tile_bp;

    std::vector<std::pair<uint64_t, uint64_t>> base_tiles;
    for (uint64_t ts = 0; ts < lenX; ) {
        uint64_t te = std::min<uint64_t>(ts + init_tile_bp, lenX);
        base_tiles.push_back({ts, te});
        if (te == lenX) break;
        ts += tile_step;
    }
    std::reverse(base_tiles.begin(), base_tiles.end());

    int env_query_mode = parse_query_mode_env();
    int last_mode_for_prof = env_query_mode;
    bool last_is_direct_for_prof = false;
    bool recorded_tile0_stats = false;

    if (h->is_minimizer) {
        out->index_kept = h->min_idx.n_kept_main + h->min_idx.ovf_n_kept;
        out->index_valid = h->min_idx.n_raw_min;
        out->bucket_bits = h->min_idx.bucket_bits;

        if (!waited_seqx && seqx_ready_ev) {
            cudaError_t w_err = cudaStreamWaitEvent(0, seqx_ready_ev, 0);
            if (w_err != cudaSuccess) {
                return fail_cleanup("cudaStreamWaitEvent(seqx_ready_ev)", w_err);
            }
            waited_seqx = true;
        }

        last_mode_for_prof = (env_query_mode == GECKO_QMODE_MONOLITHIC)
            ? GECKO_QMODE_MONOLITHIC
            : GECKO_QMODE_L2_PART;

        std::vector<std::pair<uint64_t, uint64_t>> tile_stack = base_tiles;
        while (!tile_stack.empty()) {
            uint64_t tile_start = tile_stack.back().first;
            uint64_t tile_end = tile_stack.back().second;
            tile_stack.pop_back();
            if (h->prof.enabled) h->prof.num_tiles++;

            HitKeyLayout key_layout = make_hit_key_layout(
                tile_start, tile_end, viewY.len, h->wl, h->compare_reverse);

            uint32_t raw_leftmax = 0U;
            bool hit_overflow = false;
            uint64_t *tile_sorted_hits = d_hits;
            uint64_t *tile_long_seeds = d_hits_sorted;

            uint32_t n_hits = gpu_query_tile_minimizer(
                viewX, viewY, h->min_idx,
                h->wl, h->compare_reverse,
                tile_start, tile_end, key_layout,
                h->l2_cache_bytes, env_query_mode,
                d_hits, d_hits_sorted, d_hit_count,
                h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                hits_cap, raw_leftmax, hit_overflow,
                tile_sorted_hits, tile_long_seeds, &h->prof);

            if (hit_overflow && (tile_end - tile_start) > 1ULL) {
                uint64_t mid = tile_start + (tile_end - tile_start) / 2ULL;
                tile_stack.push_back({mid, tile_end});
                tile_stack.push_back({tile_start, mid});
                continue;
            }
            if (n_hits == 0U) {
                out->leftmax_hits += raw_leftmax;
                continue;
            }

            bool frag_overflow = false;
            uint32_t n_frags = gpu_extend_tile_two_tier(
                tile_sorted_hits, n_hits, key_layout, viewX, viewY, viewY_rc,
                h->wl, false, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
                tile_long_seeds, d_long_count, d_out_frags, d_frag_count,
                frags_cap, frag_overflow, &h->prof);
            if (frag_overflow && (tile_end - tile_start) > 1ULL) {
                uint64_t mid = tile_start + (tile_end - tile_start) / 2ULL;
                tile_stack.push_back({mid, tile_end});
                tile_stack.push_back({tile_start, mid});
                continue;
            }

            out->leftmax_hits += raw_leftmax;
            out->filtered_hits += n_hits;
            out->extended_hits += n_hits;

            if (n_frags > 0U) {
                size_t old_n = all_frags.size();
                all_frags.resize(old_n + n_frags);
                std::chrono::steady_clock::time_point t_d2h{};
                if (h->prof.enabled) t_d2h = std::chrono::steady_clock::now();
                if (download_d2h_staged(
                        all_frags.data() + old_n,
                        d_out_frags,
                        (size_t)n_frags * sizeof(GpuFrag),
                        &h->scratch.stager) != 0) {
                    return fail_cleanup("download_d2h_staged", cudaSuccess);
                }
                if (h->prof.enabled) {
                    h->prof.d2h_frags_ms += std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t_d2h).count();
                }
            }
        }
    } else {
        for (uint32_t slab = 0; slab < h->num_slabs; slab++) {
            const GpuRefIndex *active_idx = &h->idx;
            if (h->num_slabs > 1U) {
                uint32_t part_start = (slab * 256U) / h->num_slabs;
                uint32_t part_end = ((slab + 1U) * 256U) / h->num_slabs;
                if (part_start >= part_end) continue;
                std::chrono::steady_clock::time_point t_slab{};
                if (h->prof.enabled) t_slab = std::chrono::steady_clock::now();
                transient_idx = gpu_build_ref_index(
                    viewY, h->wl, h->seed_mask, h->freq_limit,
                    h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                    part_start, part_end, nullptr, 0U, false, nullptr);
                has_transient_idx = true;
                if (h->prof.enabled) {
                    h->prof.slab_rebuild_ms += std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t_slab).count();
                }
                transient_idx.d_dropped_bitmap = h->d_global_dropped_bm;
                transient_idx.owns_dropped_bitmap = false;
                transient_idx.has_global_dropped_bm = !h->is_spaced;
                active_idx = &transient_idx;
            }

            out->index_kept += active_idx->n_kept;
            out->index_valid += active_idx->n_valid;
            out->bucket_bits = active_idx->bucket_bits;

            /* Ensure async seqX H2D + fill_runs_bitmap_kernel has completed on device before
             * the first query kernel reads viewX (placed AFTER slab-0 index build for overlap!). */
            if (!waited_seqx && seqx_ready_ev) {
                cudaError_t w_err = cudaStreamWaitEvent(0, seqx_ready_ev, 0);
                if (w_err != cudaSuccess) {
                    return fail_cleanup("cudaStreamWaitEvent(seqx_ready_ev)", w_err);
                }
                waited_seqx = true;
            }

            /* F16 + W3-F13b: Regime-optimal adaptive query mode per slab */
            int slab_query_mode = env_query_mode;
            if (slab_query_mode == GECKO_QMODE_AUTO) {
                if (active_idx->is_direct_csr || !active_idx->d_tags) {
                    slab_query_mode = GECKO_QMODE_MONOLITHIC;
                } else if (!active_idx->d_ref_non_unique_bm) {
                    slab_query_mode = GECKO_QMODE_L2_PART;
                }
            }
            last_is_direct_for_prof = active_idx->is_direct_csr;

            std::vector<std::pair<uint64_t, uint64_t>> tile_stack = base_tiles;
            while (!tile_stack.empty()) {
                uint64_t tile_start = tile_stack.back().first;
                uint64_t tile_end = tile_stack.back().second;
                tile_stack.pop_back();
                if (h->prof.enabled) h->prof.num_tiles++;

                HitKeyLayout key_layout = make_hit_key_layout(
                    tile_start, tile_end, viewY.len, h->wl, h->compare_reverse);

                uint32_t raw_leftmax = 0U;
                uint32_t tile_th_drop = 0U;
                bool hit_overflow = false;
                uint64_t tile_valid_q = 0ULL;
                uint64_t tile_bypassed_q = 0ULL;
                bool probe_auto_sample = (slab_query_mode == GECKO_QMODE_AUTO);
                int run_mode = probe_auto_sample ? GECKO_QMODE_HYBRID : slab_query_mode;
                bool want_stats = probe_auto_sample ||
                                  (h->prof.enabled && !recorded_tile0_stats && active_idx->d_ref_non_unique_bm != nullptr);

                uint64_t *tile_sorted_hits = d_hits;
                uint64_t *tile_long_seeds = d_hits_sorted;

                uint32_t n_hits = gpu_query_tile_left_maximal(
                    viewX, viewY, viewY_rc, *active_idx,
                    h->wl, h->seed_mask, h->compare_reverse, two_hit_window,
                    tile_start, tile_end, key_layout,
                    h->l2_cache_bytes, run_mode,
                    want_stats ? d_stats_buf : nullptr,
                    d_hits, d_hits_sorted, d_out_frags, d_hit_count,
                    h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
                    hits_cap, raw_leftmax, hit_overflow,
                    tile_sorted_hits, tile_long_seeds,
                    want_stats ? &tile_valid_q : nullptr,
                    want_stats ? &tile_bypassed_q : nullptr,
                    &tile_th_drop, &h->prof);

                if (want_stats && !recorded_tile0_stats && tile_valid_q > 0ULL) {
                    h->prof.tile0_valid_q = tile_valid_q;
                    h->prof.tile0_bypassed_q = tile_bypassed_q;
                    recorded_tile0_stats = true;
                }

                if (probe_auto_sample && tile_valid_q >= 65536ULL) {
                    if (tile_bypassed_q * 10ULL < tile_valid_q) {
                        slab_query_mode = GECKO_QMODE_L2_PART; /* < 10% bypassed: skip anchor probes */
                    } else {
                        slab_query_mode = GECKO_QMODE_HYBRID;  /* >= 10% bypassed: keep Hybrid */
                    }
                }
                last_mode_for_prof = (slab_query_mode == GECKO_QMODE_AUTO) ? run_mode : slab_query_mode;

                if (hit_overflow && (tile_end - tile_start) > 1ULL) {
                    uint64_t tlen = tile_end - tile_start;
                    uint64_t mid = tile_start + tlen / 2ULL;
                    uint64_t right_ts = (two_hit_window > 0U && tlen > 4ULL * th_overlap)
                        ? (mid - th_overlap)
                        : mid;
                    tile_stack.push_back({right_ts, tile_end});
                    tile_stack.push_back({tile_start, mid});
                    continue;
                }
                if (n_hits == 0U) {
                    out->leftmax_hits += raw_leftmax;
                    h->prof.two_hit_dropped += tile_th_drop;
                    continue;
                }

                bool frag_overflow = false;
                uint32_t n_frags = gpu_extend_tile_two_tier(
                    tile_sorted_hits, n_hits, key_layout, viewX, viewY, viewY_rc,
                    h->wl, h->is_spaced, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
                    tile_long_seeds, d_long_count, d_out_frags, d_frag_count,
                    frags_cap, frag_overflow, &h->prof);
                if (frag_overflow && (tile_end - tile_start) > 1ULL) {
                    uint64_t tlen = tile_end - tile_start;
                    uint64_t mid = tile_start + tlen / 2ULL;
                    uint64_t right_ts = (two_hit_window > 0U && tlen > 4ULL * th_overlap)
                        ? (mid - th_overlap)
                        : mid;
                    tile_stack.push_back({right_ts, tile_end});
                    tile_stack.push_back({tile_start, mid});
                    continue;
                }

                out->leftmax_hits += raw_leftmax;
                out->filtered_hits += n_hits;
                out->extended_hits += n_hits;
                h->prof.two_hit_dropped += tile_th_drop;

                if (n_frags > 0U) {
                    size_t old_n = all_frags.size();
                    all_frags.resize(old_n + n_frags);
                    std::chrono::steady_clock::time_point t_d2h{};
                    if (h->prof.enabled) t_d2h = std::chrono::steady_clock::now();
                    if (download_d2h_staged(
                            all_frags.data() + old_n,
                            d_out_frags,
                            (size_t)n_frags * sizeof(GpuFrag),
                            &h->scratch.stager) != 0) {
                        return fail_cleanup("download_d2h_staged", cudaSuccess);
                    }
                    if (h->prof.enabled) {
                        h->prof.d2h_frags_ms += std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - t_d2h).count();
                    }
                }
            }

            if (has_transient_idx) {
                gpu_free_ref_index(transient_idx);
                has_transient_idx = false;
            }
        }
    }

    cudaError_t sync_err = cudaDeviceSynchronize();
    if (sync_err != cudaSuccess) {
        return fail_cleanup("cudaDeviceSynchronize", sync_err);
    }

    if (h->prof.enabled) {
        float ms = 0.0f;
        if (uq) {
            if (uq->prof_recorded &&
                cudaEventElapsedTime(&ms, uq->prof_start, uq->prof_stop) == cudaSuccess) {
                h->prof.upload_x_ms = (double)ms;
            }
            if (uq->prof_bm_recorded &&
                cudaEventElapsedTime(&ms, uq->prof_bm_start, uq->prof_bm_stop) == cudaSuccess) {
                h->prof.fill_bm_x_ms = (double)ms;
            }
        } else {
            if (h->prof.h2d_recorded &&
                cudaEventElapsedTime(&ms, h->prof.h2d_start, h->prof.h2d_stop) == cudaSuccess) {
                h->prof.upload_x_ms = (double)ms;
            }
            if (h->prof.h2d_bm_recorded &&
                cudaEventElapsedTime(&ms, h->prof.h2d_bm_start, h->prof.h2d_bm_stop) == cudaSuccess) {
                h->prof.fill_bm_x_ms = (double)ms;
            }
        }
        double total_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t_q_start).count();
        double bypass_pct = (h->prof.tile0_valid_q > 0ULL)
            ? (100.0 * (double)h->prof.tile0_bypassed_q / (double)h->prof.tile0_valid_q)
            : 0.0;
        fprintf(
            stderr,
            "[gecko-cuda-profile] query_gpu total_ms=%.2f scratch_ms=%.2f upload_x_ms=%.2f "
            "fill_bm_x_ms=%.2f upload_x_ovl=%d slab_idx_ms=%.2f l2_pre_bypass_ms=%.2f "
            "l2_count_ms=%.2f l2_scatter_ms=%.2f l2_probe_ms=%.2f mono_probe_ms=%.2f "
            "sort_hits_ms=%.2f extend_short_ms=%.2f extend_long_ms=%.2f d2h_ms=%.2f "
            "mode=%s P=%u slabs=%u tiles=%u tile0_bypass_pct=%.2f valid_q0=%llu "
            "bypass_q0=%llu raw_hits=%llu ext_hits=%llu long_seeds=%llu frags=%llu "
            "seeding=%s two_hit_ms=%.2f two_hit_w=%u two_hit_drop=%llu scratch_pool=%s\n",
            total_ms, h->prof.ensure_scratch_ms, h->prof.upload_x_ms,
            h->prof.fill_bm_x_ms, h->prof.upload_x_overlapped ? 1 : 0,
            h->prof.slab_rebuild_ms, h->prof.l2_pre_bypass_ms,
            h->prof.l2_count_ms, h->prof.l2_scatter_ms, h->prof.l2_probe_ms,
            h->prof.mono_probe_ms, h->prof.sort_hits_ms,
            h->prof.extend_short_ms, h->prof.extend_long_ms, h->prof.d2h_frags_ms,
            query_mode_name(last_mode_for_prof, last_is_direct_for_prof),
            h->prof.l2_part_P, h->num_slabs, h->prof.num_tiles, bypass_pct,
            (unsigned long long)h->prof.tile0_valid_q,
            (unsigned long long)h->prof.tile0_bypassed_q,
            (unsigned long long)out->leftmax_hits,
            (unsigned long long)out->extended_hits,
            (unsigned long long)h->prof.long_seeds,
            (unsigned long long)all_frags.size(),
            h->prof.is_minimizer ? "minimizer" : "exact",
            h->prof.two_hit_ms, h->prof.two_hit_window,
            (unsigned long long)h->prof.two_hit_dropped,
            h->prof.scratch_pool_hit ? "hit" : "miss");
    }

    if (owns_local_seqX) {
        free_gpu_seq(local_seqX);
    }

    if (!all_frags.empty()) {
        out->frags = (GpuFrag *)malloc(all_frags.size() * sizeof(GpuFrag));
        if (!out->frags) return -2;
        memcpy(out->frags, all_frags.data(), all_frags.size() * sizeof(GpuFrag));
        out->n_frags = all_frags.size();
    }
    return 0;
}

extern "C" int gecko_cuda_query_and_extend(
    GeckoCudaRefHandle *h,
    const GeckoPackedSeqHost *query_seq,
    const GeckoCudaConfig *cfg,
    GeckoCudaResult *out
) {
    return run_query_and_extend_impl(h, query_seq, nullptr, cfg, out);
}

extern "C" int gecko_cuda_query_and_extend_uploaded(
    GeckoCudaRefHandle *h,
    GeckoCudaUploadedQuery *uq,
    const GeckoCudaConfig *cfg,
    GeckoCudaResult *out
) {
    return run_query_and_extend_impl(h, nullptr, uq, cfg, out);
}

extern "C" void gecko_cuda_free_ref(GeckoCudaRefHandle *h) {
    if (!h) return;
    h->prof.destroy();
    gpu_free_min_index(h->min_idx);
    gpu_free_ref_index(h->idx);
    if (h->d_global_dropped_bm) (void)cudaFree(h->d_global_dropped_bm);
    free_gpu_seq(h->seqY);
    free_gpu_seq(h->seqY_rc);

    /* W3-F15a: Recycle scratch back to CudaScratchPool ONLY after stream drain & clean health check */
    bool can_recycle = h->scratch_healthy && !h->scratch.stager.failed;
    if (can_recycle) {
        if (cudaPeekAtLastError() != cudaSuccess) {
            can_recycle = false;
        } else {
            cudaError_t e0 = cudaStreamSynchronize(0);
            cudaError_t e1 = h->scratch.h2d_stream
                ? cudaStreamSynchronize(h->scratch.h2d_stream)
                : cudaSuccess;
            if (e0 != cudaSuccess || e1 != cudaSuccess || cudaPeekAtLastError() != cudaSuccess) {
                can_recycle = false;
            } else {
                h->scratch.stager.slot_in_use[0] = false;
                h->scratch.stager.slot_in_use[1] = false;
                h->scratch.stager.next_slot = 0;
            }
        }
    }
    get_scratch_pool().checkin(std::move(h->scratch), can_recycle);
    h->scratch = GeckoCudaScratch{};
    delete h;
}

extern "C" void gecko_cuda_free_frags(GpuFrag *frags) {
    free(frags);
}

/* ---- W3-F15a: Persistent CudaScratchPool Control & Telemetry FFI ---- */
extern "C" void gecko_cuda_scratch_pool_release(void) {
    get_scratch_pool().release_all();
}

extern "C" void gecko_cuda_scratch_pool_stats(GeckoCudaScratchPoolStats *out) {
    if (!out) return;
    *out = get_scratch_pool().stats();
}

extern "C" void gecko_cuda_scratch_pool_set_enabled(int enabled) {
    get_scratch_pool().set_enabled(enabled != 0);
}

extern "C" int gecko_cuda_scratch_pool_is_enabled(void) {
    return get_scratch_pool().is_enabled() ? 1 : 0;
}

extern "C" void gecko_cuda_scratch_pool_set_max_idle(size_t max_idle) {
    get_scratch_pool().set_max_idle(max_idle);
}

/* ---- Host-only C-ABI parity test wrappers (zero GPU calls, execute exact header logic) ---- */
extern "C" void gecko_cuda_fill_runs_bitmap_host_ref(
    const GeckoRun *runs,
    uint32_t n_runs,
    uint32_t len,
    uint32_t *out_bm
) {
    gecko_fill_runs_bitmap_host(runs, n_runs, len, out_bm);
}

extern "C" int gecko_cuda_parse_profile_env_str(const char *s) {
    return gecko_parse_profile_str(s) ? 1 : 0;
}

extern "C" uint64_t gecko_cuda_mix48_host(uint64_t canon_kp) {
    return gecko_mix48(canon_kp);
}

extern "C" uint64_t gecko_cuda_mix48_inv_host(uint64_t h48) {
    return gecko_mix48_inv(h48);
}

extern "C" int gecko_cuda_minimizer_params_host(int wl, int *out_kp, int *out_w) {
    int kp = 0, w = 0;
    if (!gecko_minimizer_params(wl, kp, w)) return 0;
    if (out_kp) *out_kp = kp;
    if (out_w) *out_w = w;
    return 1;
}

extern "C" void gecko_cuda_two_hit_scan_host(
    const uint64_t *sorted_hits,
    uint32_t n_hits,
    uint64_t tile_start,
    uint64_t tile_end,
    uint64_t lenY,
    int wl,
    int compare_reverse,
    uint32_t two_hit_window,
    uint8_t *out_flags
) {
    if (!sorted_hits || !out_flags || n_hits == 0U) return;
    HitKeyLayout lay = make_hit_key_layout(tile_start, tile_end, lenY, wl, compare_reverse != 0);
    uint64_t wl_u64 = (uint64_t)wl;
    uint64_t max_dist = wl_u64 + (uint64_t)two_hit_window;
    for (uint32_t i = 0U; i < n_hits; i++) {
        out_flags[i] = gecko_two_hit_scan_backward(sorted_hits, i, lay, wl_u64, max_dist);
    }
}
