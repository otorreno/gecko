/*
 * gecko_cuda_main.cu
 *
 * Thin CLI orchestrator for the modular GECKO CUDA backend:
 *   cuda_io -> cuda_index -> cuda_query -> cuda_extend -> cuda_io
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <algorithm>
#include <vector>

#include "cuda_common.cuh"

int main(int argc, char **argv) {
    if (argc == 2 && strcmp(argv[1], "--check") == 0) {
        int dev_count = 0;
        if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) return 2;
        cudaDeviceProp prop;
        if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess) return 2;
        printf("CUDA device: %s (SM %d.%d, %.1f GB VRAM)\n",
               prop.name, prop.major, prop.minor, (double)prop.totalGlobalMem / (1024.0 * 1024.0 * 1024.0));
        return 0;
    }

    if (argc < 8 || argc > 10) {
        fprintf(stderr, "Usage: %s <query.fa> <ref.fa> <out_prefix> <length> <similarity> <WL> <fixedL> [unused_threads] [freqLimit]\n", argv[0]);
        return 1;
    }

    const char *seqX_path = argv[1];
    const char *seqY_path = argv[2];
    const char *out_prefix = argv[3];
    uint64_t Lmin = strtoull(argv[4], nullptr, 10);
    uint64_t SimTh = strtoull(argv[5], nullptr, 10);
    int WL = atoi(argv[6]);
    int fixedL = atoi(argv[7]);
    uint32_t freq_limit = (argc >= 10) ? (uint32_t)strtoul(argv[9], nullptr, 10) : 1000U;

    if (WL < 8 || WL > 32) {
        fprintf(stderr, "Error: WL must be between 8 and 32 (got %d)\n", WL);
        return 1;
    }

    int dev_count = 0;
    if (cudaGetDeviceCount(&dev_count) != cudaSuccess || dev_count == 0) {
        fprintf(stderr, "Error: no CUDA-capable GPU detected.\n");
        return 2;
    }
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));

    double t0 = gecko_now_sec();

    /* 1. Load FASTA sequences on host */
    FastSeq *sX = gecko_load_fasta(seqX_path);
    FastSeq *sY = gecko_load_fasta(seqY_path);
    if (sX->buf_len >= (1ULL << 31) || sY->buf_len >= (1ULL << 31)) {
        fprintf(stderr, "Error: sequences >= 2.14 Gb exceed 31-bit packed hit coordinates.\n");
        return 1;
    }
    uint64_t min_len = fixedL ? Lmin : (uint64_t)((double)std::min(sX->n_bases, sY->n_bases) * ((double)Lmin / 100.0));
    FastSeq sY_rc_meta = gecko_build_rc_meta(sY);
    double t_load = gecko_now_sec();

    /* 2. Upload sequences & build reverse complement on GPU */
    uint8_t *d_seqX = nullptr, *d_seqY = nullptr, *d_seqY_rc = nullptr;
    CUDA_CHECK(cudaMalloc(&d_seqX, sX->buf_len + 64));
    CUDA_CHECK(cudaMalloc(&d_seqY, sY->buf_len + 64));
    CUDA_CHECK(cudaMalloc(&d_seqY_rc, sY->buf_len + 64));
    CUDA_CHECK(cudaMemcpy(d_seqX, sX->data, sX->buf_len + 64, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_seqY, sY->data, sY->buf_len + 64, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_seqY_rc + sY->buf_len, 0, 64));
    gpu_build_revcomp(d_seqY, d_seqY_rc, sY->buf_len);

    /* 3. Stage 1: GPU Canonical k-mer Index Construction */
    GpuRefIndex idx = gpu_build_ref_index(d_seqY, sY->buf_len, WL, freq_limit);
    CUDA_CHECK(cudaDeviceSynchronize());
    double t_index = gecko_now_sec();

    /* 4. Stage 2 & 3: Tiled GPU Query + Left-Maximal Filter + Two-Tier Extension */
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

    std::vector<struct FragFile> frags_f, frags_r;
    std::vector<GpuFrag> h_tile_frags;
    long coverage[1000][100];
    memset(coverage, 0, sizeof(coverage));
    uint64_t tot_leftmax_hits = 0, tot_extended = 0;

    for (uint64_t tile_start = 0; tile_start < sX->buf_len; tile_start += GECKO_QUERY_TILE_BP) {
        uint64_t tile_end = std::min<uint64_t>(tile_start + GECKO_QUERY_TILE_BP, sX->buf_len);
        uint32_t raw_leftmax = 0;
        uint32_t n_filt = gpu_query_tile_left_maximal(
            d_seqX, sX->buf_len, d_seqY, d_seqY_rc, sY->buf_len, idx, WL,
            tile_start, tile_end,
            d_hits, d_hits_sorted, d_hit_flags, d_hit_count, d_filt_count,
            (uint32_t)GECKO_MAX_TILE_HITS, raw_leftmax);
        tot_leftmax_hits += raw_leftmax;
        tot_extended += n_filt;
        if (n_filt == 0) continue;

        uint32_t n_frags = gpu_extend_tile_two_tier(
            d_hits, n_filt, d_seqX, sX->buf_len, d_seqY, d_seqY_rc, sY->buf_len,
            WL, (uint32_t)min_len, (float)SimTh,
            d_long_seeds, d_long_count, d_out_frags, d_frag_count, (uint32_t)GECKO_MAX_TILE_FRAGS);

        if (n_frags > 0) {
            h_tile_frags.resize(n_frags);
            CUDA_CHECK(cudaMemcpy(h_tile_frags.data(), d_out_frags, (size_t)n_frags * sizeof(GpuFrag), cudaMemcpyDeviceToHost));
            for (uint32_t fi = 0; fi < n_frags; fi++) {
                const GpuFrag &gf = h_tile_frags[fi];
                uint64_t seqX_id = gecko_find_seq_id(sX, gf.xStart);
                uint64_t seqY_buf_id = gf.strand ? gecko_find_seq_id(&sY_rc_meta, gf.yStart) : gecko_find_seq_id(sY, gf.yStart);
                uint64_t length = (uint64_t)(gf.xEnd - gf.xStart + 1);
                float sim = (float)((double)gf.score * 100.0 / (double)(length * GECKO_POINT));

                struct FragFile ff;
                ff.diag = gf.diag;
                ff.xStart = (uint64_t)gf.xStart - seqX_id;
                ff.yStart = (uint64_t)gf.yStart - seqY_buf_id;
                ff.xEnd = (uint64_t)gf.xEnd - seqX_id;
                ff.yEnd = (uint64_t)gf.yEnd - seqY_buf_id;
                ff.length = length;
                ff.ident = gf.ident;
                ff.score = gf.score;
                ff.similarity = sim;
                ff.seqX = seqX_id;
                ff.seqY = gf.strand ? (sY->n_seqs - seqY_buf_id - 1) : seqY_buf_id;
                ff.block = 0;
                ff.strand = gf.strand ? 'r' : 'f';

                int lb = (length < 999) ? (int)length : 999;
                int sb = (int)sim;
                if (sb < 0) sb = 0;
                if (sb > 99) sb = 99;
                coverage[lb][sb]++;

                if (gf.strand) frags_r.push_back(ff);
                else frags_f.push_back(ff);
            }
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    double t_gpu = gecko_now_sec();

    /* 5. Deduplicate & write GECKO-compatible outputs */
    gecko_dedup_frags(frags_f);
    gecko_dedup_frags(frags_r);

    gecko_write_outputs(out_prefix, seqX_path, seqY_path, sX, sY, Lmin, SimTh,
                        tot_leftmax_hits, tot_extended, frags_f, frags_r, coverage);

    double t_end = gecko_now_sec();

    fprintf(stdout,
        "gecko_cuda summary (%s, SM %d.%d):\n"
        "  SeqX: %" PRIu64 " bp (%" PRIu64 " seqs) | SeqY: %" PRIu64 " bp (%" PRIu64 " seqs)\n"
        "  Index: GPU-Hash-CSR (WL=%d, FL=%u, P=%d, kept=%u/%u)\n"
        "  Left-max seeds: %" PRIu64 " | Extended: %" PRIu64 " | Fragments: %zu (F:%zu R:%zu)\n"
        "  Timing: load=%.3fs  gpu_index=%.3fs  gpu_query+ext=%.3fs  write=%.3fs  TOTAL=%.3fs\n",
        prop.name, prop.major, prop.minor,
        sX->n_bases, sX->n_seqs, sY->n_bases, sY->n_seqs,
        WL, freq_limit, idx.bucket_bits, idx.n_kept, idx.n_valid,
        tot_leftmax_hits, tot_extended, frags_f.size() + frags_r.size(), frags_f.size(), frags_r.size(),
        t_load - t0, t_index - t_load, t_gpu - t_index, t_end - t_gpu, t_end - t0);

    gpu_free_ref_index(idx);
    CUDA_CHECK(cudaFree(d_seqX));
    CUDA_CHECK(cudaFree(d_seqY));
    CUDA_CHECK(cudaFree(d_seqY_rc));
    CUDA_CHECK(cudaFree(d_hits));
    CUDA_CHECK(cudaFree(d_hits_sorted));
    CUDA_CHECK(cudaFree(d_long_seeds));
    CUDA_CHECK(cudaFree(d_hit_flags));
    CUDA_CHECK(cudaFree(d_out_frags));
    CUDA_CHECK(cudaFree(d_hit_count));
    CUDA_CHECK(cudaFree(d_filt_count));
    CUDA_CHECK(cudaFree(d_long_count));
    CUDA_CHECK(cudaFree(d_frag_count));

    free(sX->data); delete sX;
    free(sY->data); delete sY;
    return 0;
}
