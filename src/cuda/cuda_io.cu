/*
 * cuda_io.cu
 *
 * Host FASTA loading, multi-FASTA coordinate resolution, fragment deduplication,
 * and GECKO-compatible output serialization (.frags, .frags.INF, .frags.MAT, .csv).
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <algorithm>
#include <sys/time.h>

#include "cuda_common.cuh"

double gecko_now_sec() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

FastSeq *gecko_load_fasta(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "Error: cannot open FASTA file: %s\n", path);
        exit(1);
    }
    fseeko(f, 0, SEEK_END);
    off_t fsize = ftello(f);
    fseeko(f, 0, SEEK_SET);
    if (fsize <= 0) {
        fprintf(stderr, "Error: empty FASTA file: %s\n", path);
        exit(1);
    }

    FastSeq *fs = new FastSeq();
    memset(fs->ident, 0, sizeof(fs->ident));
    fs->data = (uint8_t *)calloc((size_t)fsize + 64, 1);
    if (!fs->data) {
        fprintf(stderr, "OOM sequence buffer\n");
        exit(1);
    }

    const size_t io_buf_sz = 4 * 1024 * 1024;
    uint8_t *io_buf = (uint8_t *)malloc(io_buf_sz);
    if (!io_buf) {
        fprintf(stderr, "OOM io_buf\n");
        exit(1);
    }

    uint64_t pos = 0;
    int in_header = 0, capture_first_id = 0, seen_first_gt = 0;
    size_t id_len = 0;

    size_t nread;
    while ((nread = fread(io_buf, 1, io_buf_sz, f)) > 0) {
        for (size_t i = 0; i < nread; i++) {
            uint8_t c = io_buf[i];
            if (in_header) {
                if (c == '\n') {
                    in_header = 0;
                    capture_first_id = 0;
                } else if (capture_first_id) {
                    if (id_len == 0 && c == ' ') {
                    } else if (c == ' ' || c == '\r' || c == '\t') {
                        capture_first_id = 0;
                    } else if (id_len < MAXLID) {
                        fs->ident[id_len++] = (char)c;
                    }
                }
                continue;
            }
            if (c == '>') {
                in_header = 1;
                if (!seen_first_gt) {
                    seen_first_gt = 1;
                    capture_first_id = 1;
                    fs->seq_starts.push_back(0);
                } else {
                    fs->data[pos++] = 0x00;
                    fs->seq_starts.push_back(pos);
                }
                continue;
            }
            if (!seen_first_gt) continue;
            if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
                fs->data[pos++] = c;
            } else if (c == '*') {
                fs->data[pos++] = 0x00;
            }
        }
    }

    fs->buf_len = pos;
    fs->n_seqs = fs->seq_starts.size();
    fs->n_bases = (fs->n_seqs > 0) ? (pos - (fs->n_seqs - 1)) : 0;
    fs->seq_starts.push_back(pos + 1);

    free(io_buf);
    fclose(f);
    return fs;
}

FastSeq gecko_build_rc_meta(const FastSeq *sY) {
    FastSeq rc;
    memset(rc.ident, 0, sizeof(rc.ident));
    rc.data = nullptr;
    rc.buf_len = sY->buf_len;
    rc.n_bases = sY->n_bases;
    rc.n_seqs = sY->n_seqs;
    rc.seq_starts.push_back(0);
    for (uint64_t i = 0; i < sY->buf_len; i++) {
        if (sY->data[sY->buf_len - 1 - i] == 0x00) {
            rc.seq_starts.push_back(i + 1);
        }
    }
    rc.seq_starts.push_back(sY->buf_len + 1);
    return rc;
}

uint64_t gecko_find_seq_id(const FastSeq *fs, uint64_t pos) {
    if (fs->n_seqs <= 1) return 0;
    uint64_t lo = 0, hi = fs->n_seqs - 1;
    while (lo < hi) {
        uint64_t mid = lo + ((hi - lo + 1) >> 1);
        if (fs->seq_starts[mid] <= pos) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

void gecko_dedup_frags(std::vector<struct FragFile> &v) {
    if (v.size() <= 1) return;
    std::sort(v.begin(), v.end(), [](const struct FragFile &a, const struct FragFile &b) {
        if (a.diag != b.diag) return a.diag < b.diag;
        if (a.xStart != b.xStart) return a.xStart < b.xStart;
        if (a.xEnd != b.xEnd) return a.xEnd > b.xEnd;
        return a.score > b.score;
    });
    size_t out = 0;
    for (size_t i = 0; i < v.size(); i++) {
        if (out > 0 &&
            v[i].diag == v[out - 1].diag &&
            v[i].seqX == v[out - 1].seqX &&
            v[i].seqY == v[out - 1].seqY &&
            v[i].xEnd <= v[out - 1].xEnd) {
            continue;
        }
        v[out++] = v[i];
    }
    v.resize(out);
}

void gecko_write_outputs(
    const char *out_prefix,
    const char *seqX_path,
    const char *seqY_path,
    const FastSeq *sX,
    const FastSeq *sY,
    uint64_t Lmin,
    uint64_t SimTh,
    uint64_t tot_hits,
    uint64_t tot_hits_used,
    const std::vector<struct FragFile> &frags_f,
    const std::vector<struct FragFile> &frags_r,
    long coverage[1000][100]
) {
    char path[1024];
    uint64_t n0 = sX->n_bases, n1 = sY->n_bases;

    snprintf(path, sizeof(path), "%s.frags", out_prefix);
    FILE *fOut = fopen(path, "wb");
    if (!fOut) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    writeSequenceLength(&n0, fOut);
    writeSequenceLength(&n1, fOut);
    for (auto tmp : frags_f) writeFragment(&tmp, fOut);
    for (auto tmp : frags_r) writeFragment(&tmp, fOut);
    fclose(fOut);

    snprintf(path, sizeof(path), "%s.frags.MAT", out_prefix);
    FILE *fMat = fopen(path, "wt");
    if (!fMat) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 100; j++) fprintf(fMat, "%ld\t", coverage[i][j]);
        fprintf(fMat, "\n");
    }
    fclose(fMat);

    snprintf(path, sizeof(path), "%s.frags.INF", out_prefix);
    FILE *fInf = fopen(path, "wt");
    if (!fInf) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    fprintf(fInf, "All by-Identity Ungapped Fragments (gecko_cuda GPU engine)\n");
    fprintf(fInf, "[GECKO-compatible output]\n");
    fprintf(fInf, "SeqX filename        : %s\n", seqX_path);
    fprintf(fInf, "SeqY filename        : %s\n", seqY_path);
    fprintf(fInf, "SeqX name            : %s\n", sX->ident);
    fprintf(fInf, "SeqY name            : %s\n", sY->ident);
    fprintf(fInf, "SeqX length          : %" PRIu64 "\n", n0);
    fprintf(fInf, "SeqY length          : %" PRIu64 "\n", n1);
    fprintf(fInf, "Min.fragment.length  : %" PRIu64 "\n", Lmin);
    fprintf(fInf, "Min.Identity         : %" PRIu64 "\n", SimTh);
    fprintf(fInf, "Tot Hits (seeds)     : %" PRIu64 "\n", tot_hits);
    fprintf(fInf, "Tot Hits (seeds) used: %" PRIu64 "\n", tot_hits_used);
    fprintf(fInf, "Total fragments      : %zu\n", frags_f.size() + frags_r.size());
    fprintf(fInf, "========================================================\n");
    fclose(fInf);

    snprintf(path, sizeof(path), "%s.csv", out_prefix);
    FILE *fCsv = fopen(path, "wt");
    if (!fCsv) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    fprintf(fCsv, "All by-Identity Ungapped Fragments (gecko_cuda GPU engine)\n");
    fprintf(fCsv, "[GECKO-compatible output]\n");
    fprintf(fCsv, "SeqX filename        : %s\n", seqX_path);
    fprintf(fCsv, "SeqY filename        : %s\n", seqY_path);
    fprintf(fCsv, "SeqX name            : %s\n", sX->ident);
    fprintf(fCsv, "SeqY name            : %s\n", sY->ident);
    fprintf(fCsv, "SeqX length          : %" PRIu64 "\n", n0);
    fprintf(fCsv, "SeqY length          : %" PRIu64 "\n", n1);
    fprintf(fCsv, "Min.fragment.length  : %" PRIu64 "\n", Lmin);
    fprintf(fCsv, "Min.Identity         : %" PRIu64 "\n", SimTh);
    fprintf(fCsv, "Tot Hits (seeds)     : %" PRIu64 "\n", tot_hits);
    fprintf(fCsv, "Tot Hits (seeds) used: %" PRIu64 "\n", tot_hits_used);
    fprintf(fCsv, "Total fragments      : %zu\n", frags_f.size() + frags_r.size());
    fprintf(fCsv, "========================================================\n");
    fprintf(fCsv, "Total CSB: 0\n");
    fprintf(fCsv, "========================================================\n");
    fprintf(fCsv, "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");

    for (const auto &f : frags_f) {
        double pct_ident = 100.0 * ((double)f.ident / (double)f.length);
        fprintf(fCsv, "Frag,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%c,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.2f,%.2f,%" PRIu64 ",%" PRIu64 "\n",
                f.xStart, f.yStart, f.xEnd, f.yEnd, f.strand, (uint64_t)f.block,
                f.length, f.score, f.ident, f.similarity, pct_ident, f.seqX, f.seqY);
    }
    for (const auto &f : frags_r) {
        double pct_ident = 100.0 * ((double)f.ident / (double)f.length);
        uint64_t ys = n1 - f.yStart - 1;
        uint64_t ye = n1 - f.yEnd - 1;
        fprintf(fCsv, "Frag,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%c,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.2f,%.2f,%" PRIu64 ",%" PRIu64 "\n",
                f.xStart, ys, f.xEnd, ye, f.strand, (uint64_t)f.block,
                f.length, f.score, f.ident, f.similarity, pct_ident, f.seqX, f.seqY);
    }
    fclose(fCsv);
}
