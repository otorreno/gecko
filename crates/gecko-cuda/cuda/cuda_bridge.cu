/*
 * cuda_bridge.cu
 *
 * Clean 2-phase C-ABI bridge between Rust (gecko-cuda) and 2-bit PackedSeq CUDA kernels:
 *   1. gecko_cuda_probe()
 *   2. gecko_cuda_prepare_ref()  — zero-copy uploads PackedSeq Y, builds 2-bit Y_rc & GPU index
 *   3. gecko_cuda_query_and_extend() — zero-copy uploads PackedSeq X, runs adaptive tiled query + warp extension
 *   4. gecko_cuda_free_ref() / gecko_cuda_free_frags()
 */

#include <algorithm>
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
    uint64_t *d_long_seeds = nullptr;
    uint8_t *d_hit_flags = nullptr;
    GpuFrag *d_out_frags = nullptr;
    uint32_t *d_hit_count = nullptr;
    uint32_t *d_filt_count = nullptr;
    uint32_t *d_long_count = nullptr;
    uint32_t *d_frag_count = nullptr;
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
    /* Graceful fallback to pageable host memory if pinned quota is exhausted */
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

static int ensure_scratch(GeckoCudaScratch &sc) {
    if (sc.allocated) return 0;
    CUDA_CHECK(cudaMalloc(&sc.d_hits, GECKO_MAX_TILE_HITS * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_hits_sorted, GECKO_MAX_TILE_HITS * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_long_seeds, GECKO_MAX_TILE_HITS * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_hit_flags, GECKO_MAX_TILE_HITS * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_out_frags, GECKO_MAX_TILE_FRAGS * sizeof(GpuFrag)));
    CUDA_CHECK(cudaMalloc(&sc.d_hit_count, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_filt_count, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_long_count, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&sc.d_frag_count, sizeof(uint32_t)));
    sc.allocated = true;
    return 0;
}

static void free_scratch(GeckoCudaScratch &sc) {
    if (sc.d_hits) cudaFree(sc.d_hits);
    if (sc.d_hits_sorted) cudaFree(sc.d_hits_sorted);
    if (sc.d_long_seeds) cudaFree(sc.d_long_seeds);
    if (sc.d_hit_flags) cudaFree(sc.d_hit_flags);
    if (sc.d_out_frags) cudaFree(sc.d_out_frags);
    if (sc.d_hit_count) cudaFree(sc.d_hit_count);
    if (sc.d_filt_count) cudaFree(sc.d_filt_count);
    if (sc.d_long_count) cudaFree(sc.d_long_count);
    if (sc.d_frag_count) cudaFree(sc.d_frag_count);
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
    bool is_spaced;
    bool soft_mask;
    bool compare_reverse;
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
    h->is_spaced = (h->seed_mask != span_mask);
    h->soft_mask = (soft_mask != 0);
    h->compare_reverse = (compare_reverse != 0);

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

    h->idx = gpu_build_ref_index(h->seqY.view(), wl, h->seed_mask, freq_limit);
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
    out->index_kept = h->idx.n_kept;
    out->index_valid = h->idx.n_valid;
    out->bucket_bits = h->idx.bucket_bits;

    if (ensure_scratch(h->scratch) != 0) {
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

    GpuSeqView viewX = seqX.view();
    GpuSeqView viewY = h->seqY.view();
    GpuSeqView viewY_rc = h->compare_reverse ? h->seqY_rc.view() : viewY;

    std::vector<GpuFrag> all_frags;

    /* Choose an initial tile size tuned to k-mer density, and adaptively halve
     * any tile that overflows GECKO_MAX_TILE_HITS or GECKO_MAX_TILE_FRAGS so
     * zero hits and zero fragments are ever dropped. */
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
    std::vector<std::pair<uint64_t, uint64_t>> tile_stack;
    for (uint64_t ts = 0; ts < lenX; ts += init_tile_bp) {
        uint64_t te = std::min<uint64_t>(ts + init_tile_bp, lenX);
        tile_stack.push_back({ts, te});
    }
    std::reverse(tile_stack.begin(), tile_stack.end());

    while (!tile_stack.empty()) {
        uint64_t tile_start = tile_stack.back().first;
        uint64_t tile_end = tile_stack.back().second;
        tile_stack.pop_back();

        uint32_t raw_leftmax = 0;
        bool hit_overflow = false;
        uint32_t n_filt = gpu_query_tile_left_maximal(
            viewX, viewY, h->idx,
            h->wl, h->seed_mask, h->compare_reverse,
            tile_start, tile_end,
            d_hits, d_hits_sorted, d_hit_flags, d_hit_count, d_filt_count,
            h->scratch.d_cub_temp, h->scratch.cub_temp_cap,
            (uint32_t)GECKO_MAX_TILE_HITS, raw_leftmax, hit_overflow);
        if (hit_overflow && (tile_end - tile_start) > 1) {
            uint64_t mid = tile_start + (tile_end - tile_start) / 2;
            tile_stack.push_back({mid, tile_end});
            tile_stack.push_back({tile_start, mid});
            continue;
        }
        if (n_filt == 0) {
            out->leftmax_hits += raw_leftmax;
            continue;
        }

        bool frag_overflow = false;
        uint32_t n_frags = gpu_extend_tile_two_tier(
            d_hits, n_filt, viewX, viewY, viewY_rc,
            h->wl, h->is_spaced, min_len, sim_th, match_score, mismatch_penalty, n_penalty,
            d_long_seeds, d_long_count, d_out_frags, d_frag_count,
            (uint32_t)GECKO_MAX_TILE_FRAGS, frag_overflow);
        if (frag_overflow && (tile_end - tile_start) > 1) {
            uint64_t mid = tile_start + (tile_end - tile_start) / 2;
            tile_stack.push_back({mid, tile_end});
            tile_stack.push_back({tile_start, mid});
            continue;
        }

        out->leftmax_hits += raw_leftmax;
        out->extended_hits += n_filt;

        if (n_frags > 0) {
            size_t old_n = all_frags.size();
            all_frags.resize(old_n + n_frags);
            if (download_d2h_staged(
                    all_frags.data() + old_n,
                    d_out_frags,
                    (size_t)n_frags * sizeof(GpuFrag),
                    &h->scratch) != 0) {
                free_gpu_seq(seqX);
                return -1;
            }
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
    free_gpu_seq(h->seqY);
    free_gpu_seq(h->seqY_rc);
    delete h;
}

extern "C" void gecko_cuda_free_frags(GpuFrag *frags) {
    free(frags);
}
