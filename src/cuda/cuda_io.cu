/*
 * cuda_io.cu
 *
 * High-speed 4 MB chunked + line-width stride-predicted FASTA loading,
 * O(n_seqs) multi-FASTA coordinate resolution, fragment deduplication, and
 * fast buffered GECKO-compatible output serialization (.frags, .frags.INF, .frags.MAT, .csv).
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <algorithm>
#include <sys/time.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "cuda_common.cuh"

double gecko_now_sec() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static inline bool line_is_pure_alpha(const uint8_t *p, size_t len) {
    uint8_t bad = 0;
    for (size_t i = 0; i < len; i++) {
        uint8_t u = (uint8_t)((p[i] & 0xDFU) - (uint8_t)'A');
        bad |= (uint8_t)(u >= 26U);
    }
    return bad == 0;
}

static inline void append_seq_line(
    const uint8_t *p,
    size_t line_len,
    uint8_t *dst,
    uint64_t &pos
) {
    if (line_is_pure_alpha(p, line_len)) {
        memcpy(dst + pos, p, line_len);
        pos += line_len;
    } else {
        for (size_t i = 0; i < line_len; i++) {
            uint8_t c = p[i];
            if ((uint8_t)((c & 0xDFU) - (uint8_t)'A') < 26U) {
                dst[pos++] = c;
            } else if (c == '*') {
                dst[pos++] = 0x00;
            }
        }
    }
}

FastSeq *gecko_load_fasta(const char *path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "Error: cannot open FASTA file: %s\n", path);
        exit(1);
    }
    struct stat st;
    if (fstat(fd, &st) != 0 || st.st_size <= 0) {
        close(fd);
        fprintf(stderr, "Error: empty or unreadable FASTA file: %s\n", path);
        exit(1);
    }
    size_t fsize = (size_t)st.st_size;

    FastSeq *fs = new FastSeq();
    memset(fs->ident, 0, sizeof(fs->ident));
    fs->data = (uint8_t *)malloc(fsize + 64);
    if (!fs->data) {
        close(fd);
        fprintf(stderr, "OOM sequence buffer\n");
        exit(1);
    }

    const size_t CHUNK_CAP = 4 * 1024 * 1024;
    uint8_t *buf = (uint8_t *)malloc(CHUNK_CAP);
    if (!buf) {
        close(fd);
        fprintf(stderr, "OOM read buffer\n");
        exit(1);
    }

    uint8_t *dst = fs->data;
    uint64_t pos = 0;
    int seen_first_gt = 0;
    size_t carry = 0;
    size_t pred_w = 0; /* predicted constant FASTA line width (e.g. 60 or 80) */

    ssize_t nr;
    while ((nr = read(fd, buf + carry, CHUNK_CAP - carry)) > 0) {
        size_t tot = carry + (size_t)nr;
        const uint8_t *c = buf;
        const uint8_t *e = buf + tot;

        while (c < e) {
            /* Ultra-fast path: constant line width (skips memchr on 99.99% of FASTA lines) */
            if (seen_first_gt && pred_w > 0 && c + pred_w < e && c[pred_w] == '\n' && line_is_pure_alpha(c, pred_w)) {
                memcpy(dst + pos, c, pred_w);
                pos += pred_w;
                c += pred_w + 1;
                continue;
            }

            const uint8_t *nl = (const uint8_t *)memchr(c, '\n', (size_t)(e - c));
            if (!nl) break;

            size_t raw_len = (size_t)(nl - c);
            size_t line_len = raw_len;
            bool has_cr = (line_len > 0 && c[line_len - 1] == '\r');
            if (has_cr) line_len--;

            if (line_len > 0 && c[0] == '>') {
                if (!seen_first_gt) {
                    seen_first_gt = 1;
                    fs->seq_starts.push_back(0);
                    size_t k = 1;
                    while (k < line_len && c[k] == ' ') k++;
                    size_t id_len = 0;
                    while (k < line_len && c[k] != ' ' && c[k] != '\t' && id_len < MAXLID) {
                        fs->ident[id_len++] = (char)c[k++];
                    }
                    fs->ident[id_len] = '\0';
                } else {
                    dst[pos++] = 0x00;
                    fs->seq_starts.push_back(pos);
                }
            } else if (seen_first_gt && line_len > 0) {
                if (!has_cr && line_len >= 16) pred_w = line_len;
                append_seq_line(c, line_len, dst, pos);
            }
            c = nl + 1;
        }

        carry = (size_t)(e - c);
        if (carry > 0) {
            if (carry == CHUNK_CAP) {
                /* Line longer than 4MB: process directly */
                if (seen_first_gt && buf[0] != '>') append_seq_line(buf, carry, dst, pos);
                carry = 0;
            } else {
                memmove(buf, c, carry);
            }
        }
    }

    /* Process trailing partial line without final newline */
    if (carry > 0) {
        size_t line_len = carry;
        if (buf[line_len - 1] == '\r') line_len--;
        if (seen_first_gt && line_len > 0 && buf[0] != '>') {
            append_seq_line(buf, line_len, dst, pos);
        }
    }

    free(buf);
    close(fd);

    memset(dst + pos, 0, 64);
    fs->buf_len = pos;
    fs->n_seqs = fs->seq_starts.size();
    fs->n_bases = (fs->n_seqs > 0) ? (pos - (fs->n_seqs - 1)) : 0;
    fs->seq_starts.push_back(pos + 1);
    return fs;
}

/* O(n_seqs) reverse-complement sequence start computation (no full-sequence scan!) */
FastSeq gecko_build_rc_meta(const FastSeq *sY) {
    FastSeq rc;
    memset(rc.ident, 0, sizeof(rc.ident));
    rc.data = nullptr;
    rc.buf_len = sY->buf_len;
    rc.n_bases = sY->n_bases;
    rc.n_seqs = sY->n_seqs;
    rc.seq_starts.reserve((size_t)sY->n_seqs + 1);
    rc.seq_starts.push_back(0);
    /* Separator k (1 <= k < n_seqs) in sY sits at sY->seq_starts[k] - 1.
     * In reversed buffer of length L, it lands at L - sY->seq_starts[k],
     * so the next RC sequence starts at L - sY->seq_starts[k] + 1. */
    for (uint64_t k = sY->n_seqs; k > 1; k--) {
        uint64_t orig_start = sY->seq_starts[(size_t)(k - 1)];
        rc.seq_starts.push_back(sY->buf_len - orig_start + 1);
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

/* ---- Fast big-endian binary & ASCII serialization helpers ---- */
static inline void put_be64(uint8_t *dst, uint64_t v) {
    uint64_t be = __builtin_bswap64(v);
    memcpy(dst, &be, 8);
}

static inline void put_be32(uint8_t *dst, uint32_t v) {
    uint32_t be = __builtin_bswap32(v);
    memcpy(dst, &be, 4);
}

static inline void pack_frag_be93(uint8_t *dst, const struct FragFile &f) {
    put_be64(dst +  0, (uint64_t)f.diag);
    put_be64(dst +  8, f.xStart);
    put_be64(dst + 16, f.yStart);
    put_be64(dst + 24, f.xEnd);
    put_be64(dst + 32, f.yEnd);
    put_be64(dst + 40, f.length);
    put_be64(dst + 48, f.ident);
    put_be64(dst + 56, f.score);
    uint32_t sim_bits;
    memcpy(&sim_bits, &f.similarity, 4);
    put_be32(dst + 64, sim_bits);
    put_be64(dst + 68, f.seqX);
    put_be64(dst + 76, f.seqY);
    put_be64(dst + 84, (uint64_t)f.block);
    dst[92] = (uint8_t)f.strand;
}

static inline char *append_u64(char *p, uint64_t v) {
    char tmp[24];
    int n = 0;
    do {
        tmp[n++] = (char)('0' + (v % 10));
        v /= 10;
    } while (v);
    while (n > 0) *p++ = tmp[--n];
    return p;
}

static inline char *append_fixed2(char *p, double val) {
    if (val < 0.0) { *p++ = '-'; val = -val; }
    uint64_t scaled = (uint64_t)(val * 100.0 + 0.5);
    p = append_u64(p, scaled / 100);
    *p++ = '.';
    uint64_t frac = scaled % 100;
    *p++ = (char)('0' + (frac / 10));
    *p++ = (char)('0' + (frac % 10));
    return p;
}

static inline char *format_csv_frag(char *p, const struct FragFile &f, uint64_t ys, uint64_t ye) {
    double pct_ident = 100.0 * ((double)f.ident / (double)f.length);
    memcpy(p, "Frag,", 5); p += 5;
    p = append_u64(p, f.xStart); *p++ = ',';
    p = append_u64(p, ys); *p++ = ',';
    p = append_u64(p, f.xEnd); *p++ = ',';
    p = append_u64(p, ye); *p++ = ',';
    *p++ = f.strand; *p++ = ',';
    p = append_u64(p, (uint64_t)f.block); *p++ = ',';
    p = append_u64(p, f.length); *p++ = ',';
    p = append_u64(p, f.score); *p++ = ',';
    p = append_u64(p, f.ident); *p++ = ',';
    p = append_fixed2(p, (double)f.similarity); *p++ = ',';
    p = append_fixed2(p, pct_ident); *p++ = ',';
    p = append_u64(p, f.seqX); *p++ = ',';
    p = append_u64(p, f.seqY); *p++ = '\n';
    return p;
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
    const size_t BUF_CAP = 4 * 1024 * 1024;
    std::vector<uint8_t> io_buf(BUF_CAP);

    /* 1. Fast buffered binary .frags (exact 93-byte big-endian records) */
    snprintf(path, sizeof(path), "%s.frags", out_prefix);
    FILE *fOut = fopen(path, "wb");
    if (!fOut) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    put_be64(io_buf.data(), n0);
    put_be64(io_buf.data() + 8, n1);
    size_t bpos = 16;
    auto push_frag = [&](const struct FragFile &f) {
        if (bpos + 93 > BUF_CAP) {
            fwrite(io_buf.data(), 1, bpos, fOut);
            bpos = 0;
        }
        pack_frag_be93(io_buf.data() + bpos, f);
        bpos += 93;
    };
    for (const auto &f : frags_f) push_frag(f);
    for (const auto &f : frags_r) push_frag(f);
    if (bpos > 0) fwrite(io_buf.data(), 1, bpos, fOut);
    fclose(fOut);

    /* 2. Fast buffered .frags.MAT */
    snprintf(path, sizeof(path), "%s.frags.MAT", out_prefix);
    FILE *fMat = fopen(path, "wb");
    if (!fMat) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    char *cp = (char *)io_buf.data();
    char *cend = cp + BUF_CAP - 64;
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 100; j++) {
            if (cp >= cend) {
                fwrite(io_buf.data(), 1, (size_t)(cp - (char *)io_buf.data()), fMat);
                cp = (char *)io_buf.data();
            }
            cp = append_u64(cp, (uint64_t)coverage[i][j]);
            *cp++ = '\t';
        }
        *cp++ = '\n';
    }
    if (cp > (char *)io_buf.data()) {
        fwrite(io_buf.data(), 1, (size_t)(cp - (char *)io_buf.data()), fMat);
    }
    fclose(fMat);

    /* 3. .frags.INF */
    char inf_hdr[2048];
    int hdr_len = snprintf(inf_hdr, sizeof(inf_hdr),
        "All by-Identity Ungapped Fragments (gecko_cuda GPU engine)\n"
        "[GECKO-compatible output]\n"
        "SeqX filename        : %s\n"
        "SeqY filename        : %s\n"
        "SeqX name            : %s\n"
        "SeqY name            : %s\n"
        "SeqX length          : %" PRIu64 "\n"
        "SeqY length          : %" PRIu64 "\n"
        "Min.fragment.length  : %" PRIu64 "\n"
        "Min.Identity         : %" PRIu64 "\n"
        "Tot Hits (seeds)     : %" PRIu64 "\n"
        "Tot Hits (seeds) used: %" PRIu64 "\n"
        "Total fragments      : %zu\n"
        "========================================================\n",
        seqX_path, seqY_path, sX->ident, sY->ident, n0, n1, Lmin, SimTh,
        tot_hits, tot_hits_used, frags_f.size() + frags_r.size());

    snprintf(path, sizeof(path), "%s.frags.INF", out_prefix);
    FILE *fInf = fopen(path, "wb");
    if (!fInf) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    fwrite(inf_hdr, 1, (size_t)hdr_len, fInf);
    fclose(fInf);

    /* 4. Fast buffered .csv */
    snprintf(path, sizeof(path), "%s.csv", out_prefix);
    FILE *fCsv = fopen(path, "wb");
    if (!fCsv) { fprintf(stderr, "Error opening %s\n", path); exit(1); }
    fwrite(inf_hdr, 1, (size_t)hdr_len, fCsv);
    static const char csv_cols[] =
        "Total CSB: 0\n"
        "========================================================\n"
        "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY\n";
    fwrite(csv_cols, 1, sizeof(csv_cols) - 1, fCsv);

    cp = (char *)io_buf.data();
    cend = cp + BUF_CAP - 512;
    for (const auto &f : frags_f) {
        if (cp >= cend) {
            fwrite(io_buf.data(), 1, (size_t)(cp - (char *)io_buf.data()), fCsv);
            cp = (char *)io_buf.data();
        }
        cp = format_csv_frag(cp, f, f.yStart, f.yEnd);
    }
    for (const auto &f : frags_r) {
        if (cp >= cend) {
            fwrite(io_buf.data(), 1, (size_t)(cp - (char *)io_buf.data()), fCsv);
            cp = (char *)io_buf.data();
        }
        cp = format_csv_frag(cp, f, n1 - f.yStart - 1, n1 - f.yEnd - 1);
    }
    if (cp > (char *)io_buf.data()) {
        fwrite(io_buf.data(), 1, (size_t)(cp - (char *)io_buf.data()), fCsv);
    }
    fclose(fCsv);
}
