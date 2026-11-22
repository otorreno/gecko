/*
 * gecko_fast.c
 *
 * High-performance in-memory pairwise genome comparison engine for GECKO.
 * Produces GECKO-compatible ungapped High-scoring Segment Pairs (.frags, .INF, .MAT, .csv)
 * without materializing dictionaries or hit files to disk.
 *
 * Key algorithmic improvements over the out-of-core GECKO workflow:
 *  1. Single canonical k-mer index of Reference Y only (no X dictionary, no Y-revcomp dictionary).
 *  2. Direct CSR counting-sort index for WL <= 12; 2-level cache-friendly parallel radix-bucketed
 *     hash-CSR index for 13 <= WL <= 32 (O(N) build, ~12 bytes/bp, ~4 entries/bucket).
 *  3. Streaming multi-threaded query scan over X with 64-wide software prefetching, eliminating
 *     hit materialization, external sorting, and disk I/O completely.
 *  4. Lock-free 2-tier per-thread diagonal cache (separate HSP and noise caches) so hits on
 *     covered diagonals are skipped in O(1) without noise hits evicting active HSPs.
 *  5. Bit-exact AVX2-accelerated bidirectional ungapped extension matching GECKO's scoring model.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define POINT 4
#define PREFETCH_BATCH 64
#define HSP_CACHE_BITS 17     /* 131,072 slots * 16B = 2 MB per strand per thread */
#define NOISE_CACHE_BITS 16   /*  65,536 slots * 16B = 1 MB per strand per thread */
#define HSP_CACHE_SIZE (1U << HSP_CACHE_BITS)
#define NOISE_CACHE_SIZE (1U << NOISE_CACHE_BITS)
#define COARSE_BITS 8
#define NUM_PARTITIONS (1U << COARSE_BITS)
#define SCATTER_BUF_LEN 64

/* ============================================================================
 * Timing utilities
 * ============================================================================ */
static double now_sec(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

/* ============================================================================
 * Sequence representation
 * ============================================================================
 * Sequences in multi-FASTA are concatenated into a single byte array separated
 * by 0x00 bytes (matching GECKO's '*' separator offsets).
 * Original ASCII case is preserved so seeding skips soft-masked lowercase bases
 * (exact words.c behavior) while extension uppercases via (c & ~0x20) (exact
 * LeeSeqDB behavior, with 0x00 & ~0x20 == 0x00 for sequence separators).
 */
typedef struct {
    char ident[MAXLID + 1];
    unsigned char *data;     /* length = buf_len + 64 (zero-padded for safe AVX2 loads) */
    uint64_t buf_len;        /* total bytes including (n_seqs - 1) separator 0x00 bytes */
    uint64_t n_bases;        /* total nucleotide bases excluding separators */
    uint64_t n_seqs;         /* number of sequences in FASTA */
    uint64_t *seq_starts;    /* starting offset in data[] for each sequence, size n_seqs + 1 */
} FastSeq;

static uint8_t BASE_TO_2BIT[256];

static void init_base_table(void) {
    memset(BASE_TO_2BIT, 4, sizeof(BASE_TO_2BIT));
    BASE_TO_2BIT['A'] = 0;
    BASE_TO_2BIT['C'] = 1;
    BASE_TO_2BIT['G'] = 2;
    BASE_TO_2BIT['T'] = 3;
}

static unsigned char comp_base(unsigned char c) {
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

static inline int line_is_pure_alpha(const unsigned char *p, size_t len) {
    uint8_t bad = 0;
    for (size_t i = 0; i < len; i++) {
        uint8_t u = (uint8_t)((p[i] & 0xDFU) - (uint8_t)'A');
        bad |= (uint8_t)(u >= 26U);
    }
    return bad == 0;
}

static inline void append_seq_line(const unsigned char *p, size_t line_len, unsigned char *dst, uint64_t *pos) {
    if (line_is_pure_alpha(p, line_len)) {
        memcpy(dst + *pos, p, line_len);
        *pos += line_len;
    } else {
        for (size_t i = 0; i < line_len; i++) {
            unsigned char c = p[i];
            if ((uint8_t)((c & 0xDFU) - (uint8_t)'A') < 26U) {
                dst[(*pos)++] = c;
            } else if (c == '*') {
                dst[(*pos)++] = 0x00;
            }
        }
    }
}

static FastSeq *load_fasta(const char *path) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "Error: cannot open FASTA file: %s\n", path);
        exit(1);
    }
    struct stat st;
    if (fstat(fd, &st) != 0 || st.st_size <= 0) {
        close(fd);
        fprintf(stderr, "Error: empty FASTA file: %s\n", path);
        exit(1);
    }
    size_t fsize = (size_t)st.st_size;

    FastSeq *fs = (FastSeq *)calloc(1, sizeof(FastSeq));
    if (!fs) terror("OOM FastSeq");

    unsigned char *raw_alloc = (unsigned char *)malloc(fsize + 128);
    if (!raw_alloc) terror("OOM sequence buffer");
    memset(raw_alloc, 0, 64);
    fs->data = raw_alloc + 64;

    size_t starts_cap = 1024;
    fs->seq_starts = (uint64_t *)malloc(starts_cap * sizeof(uint64_t));
    if (!fs->seq_starts) terror("OOM seq_starts");

    const size_t CHUNK_CAP = 4 * 1024 * 1024;
    unsigned char *buf = (unsigned char *)malloc(CHUNK_CAP);
    if (!buf) terror("OOM read buffer");

    unsigned char *dst = fs->data;
    uint64_t pos = 0;
    uint64_t n_seqs = 0;
    int seen_first_gt = 0;
    size_t carry = 0;
    size_t pred_w = 0;

    ssize_t nr;
    while ((nr = read(fd, buf + carry, CHUNK_CAP - carry)) > 0) {
        size_t tot = carry + (size_t)nr;
        const unsigned char *c = buf;
        const unsigned char *e = buf + tot;

        while (c < e) {
            if (seen_first_gt && pred_w > 0 && c + pred_w < e && c[pred_w] == '\n' && line_is_pure_alpha(c, pred_w)) {
                memcpy(dst + pos, c, pred_w);
                pos += pred_w;
                c += pred_w + 1;
                continue;
            }
            const unsigned char *nl = (const unsigned char *)memchr(c, '\n', (size_t)(e - c));
            if (!nl) break;

            size_t line_len = (size_t)(nl - c);
            int has_cr = (line_len > 0 && c[line_len - 1] == '\r');
            if (has_cr) line_len--;

            if (line_len > 0 && c[0] == '>') {
                if (!seen_first_gt) {
                    seen_first_gt = 1;
                    fs->seq_starts[0] = 0;
                    n_seqs = 1;
                    size_t k = 1;
                    while (k < line_len && c[k] == ' ') k++;
                    size_t id_len = 0;
                    while (k < line_len && c[k] != ' ' && c[k] != '\t' && id_len < MAXLID) {
                        fs->ident[id_len++] = (char)c[k++];
                    }
                    fs->ident[id_len] = '\0';
                } else {
                    dst[pos++] = 0x00;
                    if (n_seqs + 1 >= starts_cap) {
                        starts_cap *= 2;
                        fs->seq_starts = (uint64_t *)realloc(fs->seq_starts, starts_cap * sizeof(uint64_t));
                        if (!fs->seq_starts) terror("OOM seq_starts realloc");
                    }
                    fs->seq_starts[n_seqs++] = pos;
                }
            } else if (seen_first_gt && line_len > 0) {
                if (!has_cr && line_len >= 16) pred_w = line_len;
                append_seq_line(c, line_len, dst, &pos);
            }
            c = nl + 1;
        }

        carry = (size_t)(e - c);
        if (carry > 0) {
            if (carry == CHUNK_CAP) {
                if (seen_first_gt && buf[0] != '>') append_seq_line(buf, carry, dst, &pos);
                carry = 0;
            } else {
                memmove(buf, c, carry);
            }
        }
    }

    if (carry > 0) {
        size_t line_len = carry;
        if (buf[line_len - 1] == '\r') line_len--;
        if (seen_first_gt && line_len > 0 && buf[0] != '>') {
            append_seq_line(buf, line_len, dst, &pos);
        }
    }

    free(buf);
    close(fd);

    memset(dst + pos, 0, 64);
    fs->buf_len = pos;
    fs->n_seqs = n_seqs;
    fs->n_bases = (n_seqs > 0) ? (pos - (n_seqs - 1)) : 0;
    fs->seq_starts[n_seqs] = pos + 1;
    return fs;
}

static FastSeq *build_revcomp_seq(const FastSeq *src) {
    FastSeq *rc = (FastSeq *)calloc(1, sizeof(FastSeq));
    if (!rc) terror("OOM FastSeq rc");

    memcpy(rc->ident, src->ident, sizeof(rc->ident));
    rc->buf_len = src->buf_len;
    rc->n_bases = src->n_bases;
    rc->n_seqs = src->n_seqs;

    unsigned char *raw_alloc = (unsigned char *)malloc(src->buf_len + 128);
    if (!raw_alloc) terror("OOM rc data");
    memset(raw_alloc, 0, 64);
    rc->data = raw_alloc + 64;
    memset(rc->data + src->buf_len, 0, 64);

    unsigned char comp_lut[256];
    for (int i = 0; i < 256; i++) comp_lut[i] = comp_base((unsigned char)i);

    uint64_t L = src->buf_len;
    const unsigned char *sdata = src->data;
    unsigned char *rdata = rc->data;
    for (uint64_t i = 0; i < L; i++) {
        rdata[i] = comp_lut[sdata[L - 1 - i]];
    }

    rc->seq_starts = (uint64_t *)malloc((src->n_seqs + 1) * sizeof(uint64_t));
    if (!rc->seq_starts) terror("OOM rc seq_starts");
    rc->seq_starts[0] = 0;
    uint64_t s = 1;
    for (uint64_t k = src->n_seqs; k > 1; k--) {
        rc->seq_starts[s++] = L - src->seq_starts[k - 1] + 1;
    }
    rc->seq_starts[rc->n_seqs] = L + 1;
    return rc;
}

static void free_fast_seq(FastSeq *fs) {
    if (!fs) return;
    if (fs->data) free(fs->data - 64);
    free(fs->seq_starts);
    free(fs);
}

static inline uint64_t find_seq_id(const FastSeq *fs, uint64_t pos) {
    if (fs->n_seqs <= 1) return 0;
    uint64_t lo = 0, hi = fs->n_seqs - 1;
    while (lo < hi) {
        uint64_t mid = lo + ((hi - lo + 1) >> 1);
        if (fs->seq_starts[mid] <= pos) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

/* ============================================================================
 * Invertible 64-bit hash bijection (MurmurHash3 fmix64)
 * ============================================================================ */
static inline uint64_t hash64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

/* ============================================================================
 * Reference Index (Dual-mode: Direct CSR for WL <= 12, Partitioned Hash-CSR for WL > 12)
 * ============================================================================ */
typedef struct {
    int wl;
    int direct_mode;         /* 1 if WL <= 12 (direct 4^WL table), 0 if bucketed hash-CSR */
    int bucket_bits;         /* P: total bucket bits (hash-CSR mode) */
    int sub_bits;            /* P - COARSE_BITS */
    uint32_t sub_buckets;    /* 1 << sub_bits */
    uint64_t n_indexed;      /* total indexed positions before FL filtering */
    uint32_t *offsets;       /* direct: size (4^WL + 1); hash-CSR: size NUM_PARTITIONS * (sub_buckets + 1) */
    uint64_t *keys;          /* NULL in direct mode; canonical k-mer keys in hash-CSR mode */
    uint32_t *pos_strand;    /* (pos << 1) | strand_bit */
} RefIndex;

static void free_ref_index(RefIndex *idx) {
    if (!idx) return;
    free(idx->offsets);
    free(idx->keys);
    free(idx->pos_strand);
    free(idx);
}

/* ---- Direct CSR Index Builder (WL <= 12) ---- */
static RefIndex *build_direct_index(const FastSeq *seqY, int wl, uint64_t freq_limit) {
    RefIndex *idx = (RefIndex *)calloc(1, sizeof(RefIndex));
    if (!idx) terror("OOM RefIndex");
    idx->wl = wl;
    idx->direct_mode = 1;

    uint64_t n_kmers = 1ULL << (2 * wl);
    uint64_t mask = n_kmers - 1;
    int rc_shift = 2 * (wl - 1);

    /* Count occurrences per strand separately for exact per-strand FL filtering */
    uint32_t *cnt0 = (uint32_t *)calloc(n_kmers, sizeof(uint32_t));
    uint32_t *cnt1 = (uint32_t *)calloc(n_kmers, sizeof(uint32_t));
    if (!cnt0 || !cnt1) terror("OOM direct counts");

    const unsigned char *data = seqY->data;
    uint64_t L = seqY->buf_len;
    uint64_t kf = 0, kr = 0;
    int valid = 0;

    for (uint64_t i = 0; i < L; i++) {
        uint8_t b = BASE_TO_2BIT[data[i]];
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl) {
                if (kf <= kr) cnt0[kf]++;
                else cnt1[kr]++;
            }
        } else {
            valid = 0;
            kf = 0;
            kr = 0;
        }
    }

    uint32_t *offsets = (uint32_t *)malloc((n_kmers + 1) * sizeof(uint32_t));
    if (!offsets) terror("OOM direct offsets");

    uint64_t total = 0;
    for (uint64_t k = 0; k < n_kmers; k++) {
        offsets[k] = (uint32_t)total;
        uint32_t c0 = cnt0[k];
        uint32_t c1 = cnt1[k];
        if (freq_limit > 0) {
            if (c0 > freq_limit) { cnt0[k] = 0; c0 = 0; }
            if (c1 > freq_limit) { cnt1[k] = 0; c1 = 0; }
        }
        total += (uint64_t)c0 + (uint64_t)c1;
    }
    offsets[n_kmers] = (uint32_t)total;

    uint32_t *pos_strand = (uint32_t *)malloc((total > 0 ? total : 1) * sizeof(uint32_t));
    if (!pos_strand) terror("OOM direct pos_strand");

    /* Write cursor reuses offsets[] temporarily, then restores it */
    uint32_t *cursor = (uint32_t *)malloc(n_kmers * sizeof(uint32_t));
    if (!cursor) terror("OOM direct cursor");
    memcpy(cursor, offsets, n_kmers * sizeof(uint32_t));

    kf = 0; kr = 0; valid = 0;
    for (uint64_t i = 0; i < L; i++) {
        uint8_t b = BASE_TO_2BIT[data[i]];
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl) {
                uint32_t start_pos = (uint32_t)(i + 1 - (uint64_t)wl);
                if (kf <= kr) {
                    if (cnt0[kf]) pos_strand[cursor[kf]++] = (start_pos << 1);
                } else {
                    if (cnt1[kr]) pos_strand[cursor[kr]++] = (start_pos << 1) | 1U;
                }
            }
        } else {
            valid = 0;
            kf = 0;
            kr = 0;
        }
    }

    free(cnt0);
    free(cnt1);
    free(cursor);

    idx->offsets = offsets;
    idx->pos_strand = pos_strand;
    idx->n_indexed = total;
    return idx;
}

/* ---- 2-Level Parallel Radix-Bucketed Hash-CSR Builder (13 <= WL <= 32) ---- */
typedef struct {
    const FastSeq *seqY;
    int wl;
    uint64_t start;
    uint64_t end;
    uint64_t local_counts[NUM_PARTITIONS];
    uint64_t write_offsets[NUM_PARTITIONS];
    uint64_t *keys;
    uint32_t *pos_strand;
} BuildPart1Task;

static void *build_pass1_count(void *arg) {
    BuildPart1Task *t = (BuildPart1Task *)arg;
    const unsigned char *data = t->seqY->data;
    int wl = t->wl;
    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rc_shift = 2 * (wl - 1);

    memset(t->local_counts, 0, sizeof(t->local_counts));

    /* Start (wl - 1) bases earlier to warm up rolling k-mer at t->start */
    uint64_t scan_from = (t->start >= (uint64_t)(wl - 1)) ? (t->start - (uint64_t)(wl - 1)) : 0;
    uint64_t kf = 0, kr = 0;
    int valid = 0;

    for (uint64_t i = scan_from; i < t->end; i++) {
        uint8_t b = BASE_TO_2BIT[data[i]];
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl && i >= t->start) {
                uint64_t key = (kf <= kr) ? kf : kr;
                uint32_t part = (uint32_t)(hash64(key) >> (64 - COARSE_BITS));
                t->local_counts[part]++;
            }
        } else {
            valid = 0;
            kf = 0;
            kr = 0;
        }
    }
    return NULL;
}

static void *build_pass1_scatter(void *arg) {
    BuildPart1Task *t = (BuildPart1Task *)arg;
    const unsigned char *data = t->seqY->data;
    int wl = t->wl;
    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rc_shift = 2 * (wl - 1);

    /* Small L1-resident scatter buffers per coarse partition to avoid TLB thrashing */
    uint64_t buf_k[NUM_PARTITIONS][SCATTER_BUF_LEN];
    uint32_t buf_p[NUM_PARTITIONS][SCATTER_BUF_LEN];
    uint8_t  buf_n[NUM_PARTITIONS];
    memset(buf_n, 0, sizeof(buf_n));

    uint64_t scan_from = (t->start >= (uint64_t)(wl - 1)) ? (t->start - (uint64_t)(wl - 1)) : 0;
    uint64_t kf = 0, kr = 0;
    int valid = 0;

    for (uint64_t i = scan_from; i < t->end; i++) {
        uint8_t b = BASE_TO_2BIT[data[i]];
        if (b < 4) {
            kf = ((kf << 2) | b) & mask;
            kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
            if (++valid >= wl && i >= t->start) {
                uint32_t strand = (kf <= kr) ? 0U : 1U;
                uint64_t key = strand ? kr : kf;
                uint32_t start_pos = (uint32_t)(i + 1 - (uint64_t)wl);
                uint32_t part = (uint32_t)(hash64(key) >> (64 - COARSE_BITS));

                uint8_t n = buf_n[part];
                buf_k[part][n] = key;
                buf_p[part][n] = (start_pos << 1) | strand;
                if (++n == SCATTER_BUF_LEN) {
                    uint64_t dst = t->write_offsets[part];
                    memcpy(&t->keys[dst], buf_k[part], SCATTER_BUF_LEN * sizeof(uint64_t));
                    memcpy(&t->pos_strand[dst], buf_p[part], SCATTER_BUF_LEN * sizeof(uint32_t));
                    t->write_offsets[part] = dst + SCATTER_BUF_LEN;
                    n = 0;
                }
                buf_n[part] = n;
            }
        } else {
            valid = 0;
            kf = 0;
            kr = 0;
        }
    }

    /* Flush remaining entries */
    for (uint32_t p = 0; p < NUM_PARTITIONS; p++) {
        uint8_t n = buf_n[p];
        if (n > 0) {
            uint64_t dst = t->write_offsets[p];
            memcpy(&t->keys[dst], buf_k[p], (size_t)n * sizeof(uint64_t));
            memcpy(&t->pos_strand[dst], buf_p[p], (size_t)n * sizeof(uint32_t));
            t->write_offsets[p] = dst + n;
        }
    }
    return NULL;
}

typedef struct {
    RefIndex *idx;
    const uint64_t *part_starts;
    const uint64_t *part_ends;
    uint64_t freq_limit;
    uint32_t *next_part;
} BuildPart2Shared;

static void sort_bucket_small(uint64_t *k, uint32_t *p, uint32_t n) {
    /* Insertion sort for tiny buckets (avg size ~4) */
    for (uint32_t i = 1; i < n; i++) {
        uint64_t kv = k[i];
        uint32_t pv = p[i];
        uint32_t j = i;
        while (j > 0 && (k[j - 1] > kv || (k[j - 1] == kv && p[j - 1] > pv))) {
            k[j] = k[j - 1];
            p[j] = p[j - 1];
            j--;
        }
        k[j] = kv;
        p[j] = pv;
    }
}

static void sort_bucket_qsort(uint64_t *k, uint32_t *p, int64_t lo, int64_t hi) {
    while (hi - lo > 24) {
        int64_t mid = lo + ((hi - lo) >> 1);
        /* Median-of-three pivot */
        if (k[mid] < k[lo]) {
            uint64_t tk = k[lo]; k[lo] = k[mid]; k[mid] = tk;
            uint32_t tp = p[lo]; p[lo] = p[mid]; p[mid] = tp;
        }
        if (k[hi] < k[lo]) {
            uint64_t tk = k[lo]; k[lo] = k[hi]; k[hi] = tk;
            uint32_t tp = p[lo]; p[lo] = p[hi]; p[hi] = tp;
        }
        if (k[hi] < k[mid]) {
            uint64_t tk = k[mid]; k[mid] = k[hi]; k[hi] = tk;
            uint32_t tp = p[mid]; p[mid] = p[hi]; p[hi] = tp;
        }
        uint64_t pk = k[mid];
        uint32_t pp = p[mid];
        int64_t i = lo, j = hi;
        while (i <= j) {
            while (k[i] < pk || (k[i] == pk && p[i] < pp)) i++;
            while (k[j] > pk || (k[j] == pk && p[j] > pp)) j--;
            if (i <= j) {
                uint64_t tk = k[i]; k[i] = k[j]; k[j] = tk;
                uint32_t tp = p[i]; p[i] = p[j]; p[j] = tp;
                i++; j--;
            }
        }
        if (j - lo < hi - i) {
            if (lo < j) sort_bucket_qsort(k, p, lo, j);
            lo = i;
        } else {
            if (i < hi) sort_bucket_qsort(k, p, i, hi);
            hi = j;
        }
    }
    if (hi > lo) {
        sort_bucket_small(k + lo, p + lo, (uint32_t)(hi - lo + 1));
    }
}

static void *build_pass2_worker(void *arg) {
    BuildPart2Shared *sh = (BuildPart2Shared *)arg;
    RefIndex *idx = sh->idx;
    uint32_t sub_buckets = idx->sub_buckets;
    int shift_sub = 64 - idx->bucket_bits;
    uint32_t sub_mask = sub_buckets - 1;
    uint64_t fl = sh->freq_limit;

    uint32_t *sub_cnt = (uint32_t *)malloc(sub_buckets * sizeof(uint32_t));
    uint32_t *sub_cur = (uint32_t *)malloc(sub_buckets * sizeof(uint32_t));
    if (!sub_cnt || !sub_cur) terror("OOM pass2 sub_cnt");

    size_t scratch_cap = 0;
    uint64_t *tmp_k = NULL;
    uint32_t *tmp_p = NULL;

    while (1) {
        uint32_t part = __atomic_fetch_add(sh->next_part, 1U, __ATOMIC_RELAXED);
        if (part >= NUM_PARTITIONS) break;

        uint64_t p_start = sh->part_starts[part];
        uint64_t p_end = sh->part_ends[part];
        uint32_t p_len = (uint32_t)(p_end - p_start);
        uint32_t *off_base = &idx->offsets[(uint64_t)part * (sub_buckets + 1)];

        if (p_len == 0) {
            for (uint32_t s = 0; s <= sub_buckets; s++) off_base[s] = (uint32_t)p_start;
            continue;
        }

        if (p_len > scratch_cap) {
            scratch_cap = (size_t)p_len + (size_t)(p_len >> 3) + 1024;
            free(tmp_k);
            free(tmp_p);
            tmp_k = (uint64_t *)malloc(scratch_cap * sizeof(uint64_t));
            tmp_p = (uint32_t *)malloc(scratch_cap * sizeof(uint32_t));
            if (!tmp_k || !tmp_p) terror("OOM pass2 scratch");
        }

        uint64_t *gk = &idx->keys[p_start];
        uint32_t *gp = &idx->pos_strand[p_start];

        /* Count fine sub-buckets */
        memset(sub_cnt, 0, sub_buckets * sizeof(uint32_t));
        for (uint32_t i = 0; i < p_len; i++) {
            uint32_t sb = (uint32_t)((hash64(gk[i]) >> shift_sub) & sub_mask);
            sub_cnt[sb]++;
        }

        /* Prefix sum into sub_cur */
        uint32_t acc = 0;
        for (uint32_t s = 0; s < sub_buckets; s++) {
            sub_cur[s] = acc;
            acc += sub_cnt[s];
        }

        /* Scatter into local scratch */
        for (uint32_t i = 0; i < p_len; i++) {
            uint64_t kv = gk[i];
            uint32_t sb = (uint32_t)((hash64(kv) >> shift_sub) & sub_mask);
            uint32_t dst = sub_cur[sb]++;
            tmp_k[dst] = kv;
            tmp_p[dst] = gp[i];
        }

        /* Sort each tiny sub-bucket, apply per-strand FL filtering, and write back */
        uint32_t out_pos = 0;
        uint32_t in_pos = 0;
        for (uint32_t s = 0; s < sub_buckets; s++) {
            off_base[s] = (uint32_t)(p_start + out_pos);
            uint32_t n = sub_cnt[s];
            if (n == 0) continue;

            uint64_t *bk = &tmp_k[in_pos];
            uint32_t *bp = &tmp_p[in_pos];
            in_pos += n;

            if (n > 1) {
                if (n <= 24) sort_bucket_small(bk, bp, n);
                else sort_bucket_qsort(bk, bp, 0, (int64_t)n - 1);
            }

            if (fl == 0 || n <= fl) {
                /* Fast path: impossible to exceed freq_limit */
                memcpy(&gk[out_pos], bk, (size_t)n * sizeof(uint64_t));
                memcpy(&gp[out_pos], bp, (size_t)n * sizeof(uint32_t));
                out_pos += n;
            } else {
                /* Per-strand frequency filtering for each distinct canonical key */
                uint32_t i = 0;
                while (i < n) {
                    uint64_t kv = bk[i];
                    uint32_t j = i + 1;
                    while (j < n && bk[j] == kv) j++;
                    uint32_t c0 = 0, c1 = 0;
                    for (uint32_t m = i; m < j; m++) {
                        if (bp[m] & 1U) c1++;
                        else c0++;
                    }
                    int keep0 = (c0 <= fl);
                    int keep1 = (c1 <= fl);
                    if (keep0 && keep1) {
                        uint32_t run = j - i;
                        memcpy(&gk[out_pos], &bk[i], (size_t)run * sizeof(uint64_t));
                        memcpy(&gp[out_pos], &bp[i], (size_t)run * sizeof(uint32_t));
                        out_pos += run;
                    } else if (keep0 || keep1) {
                        uint32_t want = keep1 ? 1U : 0U;
                        for (uint32_t m = i; m < j; m++) {
                            if ((bp[m] & 1U) == want) {
                                gk[out_pos] = kv;
                                gp[out_pos] = bp[m];
                                out_pos++;
                            }
                        }
                    }
                    i = j;
                }
            }
        }
        off_base[sub_buckets] = (uint32_t)(p_start + out_pos);
    }

    free(tmp_k);
    free(tmp_p);
    free(sub_cnt);
    free(sub_cur);
    return NULL;
}

static RefIndex *build_hash_csr_index(const FastSeq *seqY, int wl, uint64_t freq_limit, int n_threads) {
    RefIndex *idx = (RefIndex *)calloc(1, sizeof(RefIndex));
    if (!idx) terror("OOM RefIndex");
    idx->wl = wl;
    idx->direct_mode = 0;

    /* Choose P so average bucket size is ~3-4 entries, clamped to [14, 26] */
    int p = 14;
    while (p < 26 && (1ULL << p) * 4ULL < seqY->n_bases) p++;
    idx->bucket_bits = p;
    idx->sub_bits = p - COARSE_BITS;
    idx->sub_buckets = 1U << idx->sub_bits;

    int nt = n_threads;
    if (nt < 1) nt = 1;
    if ((uint64_t)nt > seqY->buf_len / 65536ULL + 1ULL) {
        nt = (int)(seqY->buf_len / 65536ULL + 1ULL);
    }

    BuildPart1Task *tasks = (BuildPart1Task *)calloc((size_t)nt, sizeof(BuildPart1Task));
    pthread_t *ths = (pthread_t *)malloc((size_t)nt * sizeof(pthread_t));
    if (!tasks || !ths) terror("OOM pass1 tasks");

    uint64_t chunk = (seqY->buf_len + (uint64_t)nt - 1ULL) / (uint64_t)nt;
    for (int t = 0; t < nt; t++) {
        tasks[t].seqY = seqY;
        tasks[t].wl = wl;
        tasks[t].start = (uint64_t)t * chunk;
        tasks[t].end = tasks[t].start + chunk;
        if (tasks[t].end > seqY->buf_len) tasks[t].end = seqY->buf_len;
        if (tasks[t].start > seqY->buf_len) tasks[t].start = seqY->buf_len;
        pthread_create(&ths[t], NULL, build_pass1_count, &tasks[t]);
    }
    for (int t = 0; t < nt; t++) pthread_join(ths[t], NULL);

    /* Prefix sum coarse partitions across threads */
    uint64_t part_starts[NUM_PARTITIONS];
    uint64_t part_ends[NUM_PARTITIONS];
    uint64_t total = 0;
    for (uint32_t part = 0; part < NUM_PARTITIONS; part++) {
        part_starts[part] = total;
        for (int t = 0; t < nt; t++) {
            tasks[t].write_offsets[part] = total;
            total += tasks[t].local_counts[part];
        }
        part_ends[part] = total;
    }

    idx->n_indexed = total;
    idx->keys = (uint64_t *)malloc((total > 0 ? total : 1) * sizeof(uint64_t));
    idx->pos_strand = (uint32_t *)malloc((total > 0 ? total : 1) * sizeof(uint32_t));
    uint64_t off_count = (uint64_t)NUM_PARTITIONS * ((uint64_t)idx->sub_buckets + 1ULL);
    idx->offsets = (uint32_t *)malloc(off_count * sizeof(uint32_t));
    if (!idx->keys || !idx->pos_strand || !idx->offsets) terror("OOM hash-CSR arrays");

    for (int t = 0; t < nt; t++) {
        tasks[t].keys = idx->keys;
        tasks[t].pos_strand = idx->pos_strand;
        pthread_create(&ths[t], NULL, build_pass1_scatter, &tasks[t]);
    }
    for (int t = 0; t < nt; t++) pthread_join(ths[t], NULL);

    /* Pass 2: fine bucketing + sorting + FL filtering in parallel across 256 partitions */
    uint32_t next_part = 0;
    BuildPart2Shared sh = {
        .idx = idx,
        .part_starts = part_starts,
        .part_ends = part_ends,
        .freq_limit = freq_limit,
        .next_part = &next_part
    };
    for (int t = 0; t < nt; t++) {
        pthread_create(&ths[t], NULL, build_pass2_worker, &sh);
    }
    for (int t = 0; t < nt; t++) pthread_join(ths[t], NULL);

    free(tasks);
    free(ths);
    return idx;
}

static RefIndex *build_ref_index(const FastSeq *seqY, int wl, uint64_t freq_limit, int n_threads) {
    if (wl <= 10 || (wl <= 12 && (1ULL << (2 * wl)) <= 2ULL * seqY->n_bases)) {
        return build_direct_index(seqY, wl, freq_limit);
    }
    return build_hash_csr_index(seqY, wl, freq_limit, n_threads);
}

/* ============================================================================
 * 2-Tier Lock-Free Diagonal Cache (per thread, per strand)
 * ============================================================================ */
typedef struct {
    int64_t diag;
    int64_t last_x;
} DiagSlot;

typedef struct {
    DiagSlot hsp[HSP_CACHE_SIZE];
    DiagSlot noise[NOISE_CACHE_SIZE];
} StrandDiagCache;

static void init_diag_cache(StrandDiagCache *c) {
    for (uint32_t i = 0; i < HSP_CACHE_SIZE; i++) {
        c->hsp[i].diag = INT64_MIN;
        c->hsp[i].last_x = -1;
    }
    for (uint32_t i = 0; i < NOISE_CACHE_SIZE; i++) {
        c->noise[i].diag = INT64_MIN;
        c->noise[i].last_x = -1;
    }
}

static inline uint32_t diag_hash_hsp(int64_t d) {
    return (uint32_t)(((uint64_t)d * 11400714819323198485ULL) >> (64 - HSP_CACHE_BITS));
}

static inline uint32_t diag_hash_noise(int64_t d) {
    return (uint32_t)(((uint64_t)d * 11400714819323198485ULL) >> (64 - NOISE_CACHE_BITS));
}

static inline int diag_is_covered(const StrandDiagCache *c, int64_t d, int64_t posX) {
    const DiagSlot *hs = &c->hsp[diag_hash_hsp(d)];
    if (hs->diag == d && posX <= hs->last_x) return 1;
    const DiagSlot *ns = &c->noise[diag_hash_noise(d)];
    if (ns->diag == d && posX <= ns->last_x) return 1;
    return 0;
}

static inline void diag_record_hsp(StrandDiagCache *c, int64_t d, int64_t end_x) {
    DiagSlot *hs = &c->hsp[diag_hash_hsp(d)];
    if (hs->diag != d || end_x > hs->last_x) {
        hs->diag = d;
        hs->last_x = end_x;
    }
}

static inline void diag_record_noise(StrandDiagCache *c, int64_t d, int64_t end_x) {
    DiagSlot *ns = &c->noise[diag_hash_noise(d)];
    ns->diag = d;
    ns->last_x = end_x;
}

/* ============================================================================
 * Bit-Exact AVX2-Accelerated Ungapped Bidirectional Extension
 * ============================================================================
 * Matches GECKO's FragFromHit() scoring exactly:
 *   match:    +4, nIdentities++
 *   mismatch: -4
 *   N in either: -1
 *   separator (0x00): stop
 *   fscore < 0: stop
 *   Backward pass starts with fscore = fscoreMax (from forward pass).
 */
typedef struct {
    struct FragFile *items;
    size_t count;
    size_t cap;
} FragVec;

static void frag_vec_push(FragVec *v, const struct FragFile *f) {
    if (v->count == v->cap) {
        v->cap = (v->cap == 0) ? 1024 : (v->cap * 2);
        v->items = (struct FragFile *)realloc(v->items, v->cap * sizeof(struct FragFile));
        if (!v->items) terror("OOM FragVec");
    }
    v->items[v->count++] = *f;
}

static int extend_hit(
    const FastSeq *sX,
    const FastSeq *sY,
    uint64_t posX,
    uint64_t posY,
    uint64_t seqX_id,
    int wl,
    uint64_t min_len,
    uint64_t sim_th,
    char strand,
    struct FragFile *out_frag,
    int64_t *out_xfilmax,
    long coverage[1000][100]
) {
    const unsigned char *dx = sX->data;
    const unsigned char *dy = sY->data;
    int64_t n0 = (int64_t)sX->buf_len;
    int64_t n1 = (int64_t)sY->buf_len;

    int64_t xfil = (int64_t)posX + wl;
    int64_t ycol = (int64_t)posY + wl;
    int64_t rem_fwd = (n0 - xfil < n1 - ycol) ? (n0 - xfil) : (n1 - ycol);

    int64_t xfilmax = (int64_t)posX + wl - 1;
    int64_t ycolmax = (int64_t)posY + wl - 1;
    int fscore = POINT * wl;
    int fscoreMax = fscore;
    int nIdentities = wl;
    int maxIdentities = wl;

#if defined(__AVX2__)
    const __m256i v_case_mask = _mm256_set1_epi8((char)0xDF);
    const __m256i v_zero = _mm256_setzero_si256();
    const __m256i v_N = _mm256_set1_epi8('N');
    const __m256i v_two = _mm256_set1_epi8(2);
    const __m256i v_neg1 = _mm256_set1_epi8(-1);
    const __m256i v_ff = _mm256_set1_epi8((char)0xFF);
    const __m128i v_bcast15 = _mm_set1_epi8(15);
    const __m256i v_rev16 = _mm256_setr_epi8(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
    );

    /* Forward pass: pure AVX2 32-byte prefix-sum / prefix-max (bit-exact with scalar) */
    while (rem_fwd >= 32) {
        __m256i vx = _mm256_and_si256(_mm256_loadu_si256((const __m256i *)(dx + xfil)), v_case_mask);
        __m256i vy = _mm256_and_si256(_mm256_loadu_si256((const __m256i *)(dy + ycol)), v_case_mask);

        __m256i sep = _mm256_or_si256(_mm256_cmpeq_epi8(vx, v_zero), _mm256_cmpeq_epi8(vy, v_zero));
        __m256i hasN = _mm256_or_si256(_mm256_cmpeq_epi8(vx, v_N), _mm256_cmpeq_epi8(vy, v_N));
        if (_mm256_testz_si256(_mm256_or_si256(sep, hasN), v_ff)) {
            __m256i eq = _mm256_cmpeq_epi8(vx, vy);
            uint32_t eq_mask = (uint32_t)_mm256_movemask_epi8(eq);
            if (eq_mask == 0xFFFFFFFFU) {
                fscore += 32 * POINT;
                nIdentities += 32;
                xfil += 32;
                ycol += 32;
                rem_fwd -= 32;
                if (fscore > fscoreMax) {
                    fscoreMax = fscore;
                    xfilmax = xfil - 1;
                    ycolmax = ycol - 1;
                    maxIdentities = nIdentities;
                }
                continue;
            }
            if (fscore >= 32 * POINT && fscore + 32 * POINT <= fscoreMax) {
                int m = __builtin_popcount(eq_mask);
                fscore += (2 * m - 32) * POINT;
                nIdentities += m;
                xfil += 32;
                ycol += 32;
                rem_fwd -= 32;
                continue;
            }
            /* Pure AVX2 32-byte prefix-sum / prefix-max for near-peak or near-zero chunks */
            __m256i s = _mm256_add_epi8(_mm256_and_si256(eq, v_two), v_neg1);
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 1));
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 2));
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 4));
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 8));
            __m128i lo_lane = _mm256_castsi256_si128(s);
            __m128i lo_sum = _mm_shuffle_epi8(lo_lane, v_bcast15);
            __m256i psum = _mm256_add_epi8(s, _mm256_inserti128_si256(v_zero, lo_sum, 1));

            __m128i p_lo = _mm256_castsi256_si128(psum);
            __m128i p_hi = _mm256_extracti128_si256(psum, 1);
            __m128i mn = _mm_min_epi8(p_lo, p_hi);
            __m128i mx = _mm_max_epi8(p_lo, p_hi);
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 8));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 8));
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 4));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 4));
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 2));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 2));
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 1));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 1));
            int chunk_min = (int)(int8_t)_mm_extract_epi8(mn, 0);
            int chunk_max = (int)(int8_t)_mm_extract_epi8(mx, 0);

            if (__builtin_expect(fscore + POINT * chunk_min < 0, 0)) {
                for (int k = 0; k < 32; k++) {
                    if ((eq_mask >> k) & 1U) {
                        fscore += POINT;
                        nIdentities++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax = xfil + k;
                            ycolmax = ycol + k;
                            maxIdentities = nIdentities;
                        }
                    } else {
                        fscore -= POINT;
                        if (fscore < 0) goto fwd_done;
                    }
                }
            }

            if (fscore + POINT * chunk_max > fscoreMax) {
                uint32_t max_mask = (uint32_t)_mm256_movemask_epi8(
                    _mm256_cmpeq_epi8(psum, _mm256_set1_epi8((char)chunk_max)));
                int k = __builtin_ctz(max_mask);
                fscoreMax = fscore + POINT * chunk_max;
                xfilmax = xfil + k;
                ycolmax = ycol + k;
                maxIdentities = nIdentities + ((k + 1 + chunk_max) >> 1);
            }

            int total_delta = (int)(int8_t)_mm_extract_epi8(p_hi, 15);
            fscore += POINT * total_delta;
            nIdentities += (32 + total_delta) >> 1;
            xfil += 32;
            ycol += 32;
            rem_fwd -= 32;
            continue;
        }
        /* Chunk contains N or separator: handle scalar up to 32 bytes */
        for (int k = 0; k < 32; k++) {
            unsigned char vx_c = dx[xfil] & 0xDF;
            unsigned char vy_c = dy[ycol] & 0xDF;
            if (vx_c == 0 || vy_c == 0) goto fwd_done;
            if (vx_c == 'N' || vy_c == 'N') {
                fscore -= 1;
            } else if (vx_c == vy_c) {
                fscore += POINT;
                nIdentities++;
                if (fscoreMax < fscore) {
                    fscoreMax = fscore;
                    xfilmax = xfil;
                    ycolmax = ycol;
                    maxIdentities = nIdentities;
                }
            } else {
                fscore -= POINT;
            }
            xfil++;
            ycol++;
            rem_fwd--;
            if (fscore < 0) goto fwd_done;
        }
    }
#endif

    /* Forward scalar tail */
    while (rem_fwd > 0) {
        unsigned char vx_c = dx[xfil] & 0xDF;
        unsigned char vy_c = dy[ycol] & 0xDF;
        if (vx_c == 0 || vy_c == 0) break;
        if (vx_c == 'N' || vy_c == 'N') {
            fscore -= 1;
        } else if (vx_c == vy_c) {
            fscore += POINT;
            nIdentities++;
            if (fscoreMax < fscore) {
                fscoreMax = fscore;
                xfilmax = xfil;
                ycolmax = ycol;
                maxIdentities = nIdentities;
            }
        } else {
            fscore -= POINT;
        }
        xfil++;
        ycol++;
        rem_fwd--;
        if (fscore < 0) break;
    }

fwd_done:;
    /* Backward pass: starts with fscore = fscoreMax (exact FragHits.c semantics) */
    int64_t xfil2 = (int64_t)posX - 1;
    int64_t ycol2 = (int64_t)posY - 1;
    int64_t rem_bwd = (posX < posY) ? (int64_t)posX : (int64_t)posY;
    int64_t xfilmax2 = (int64_t)posX;
    int64_t ycolmax2 = (int64_t)posY;
    fscore = fscoreMax;
    nIdentities = maxIdentities;

#if defined(__AVX2__)
    while (rem_bwd >= 32) {
        int64_t bx = xfil2 - 31;
        int64_t by = ycol2 - 31;
        __m256i vx = _mm256_and_si256(_mm256_loadu_si256((const __m256i *)(dx + bx)), v_case_mask);
        __m256i vy = _mm256_and_si256(_mm256_loadu_si256((const __m256i *)(dy + by)), v_case_mask);

        __m256i sep = _mm256_or_si256(_mm256_cmpeq_epi8(vx, v_zero), _mm256_cmpeq_epi8(vy, v_zero));
        __m256i hasN = _mm256_or_si256(_mm256_cmpeq_epi8(vx, v_N), _mm256_cmpeq_epi8(vy, v_N));
        if (_mm256_testz_si256(_mm256_or_si256(sep, hasN), v_ff)) {
            __m256i eq = _mm256_cmpeq_epi8(vx, vy);
            uint32_t eq_mask_fwd = (uint32_t)_mm256_movemask_epi8(eq);
            if (eq_mask_fwd == 0xFFFFFFFFU) {
                fscore += 32 * POINT;
                nIdentities += 32;
                xfil2 -= 32;
                ycol2 -= 32;
                rem_bwd -= 32;
                if (fscore > fscoreMax) {
                    fscoreMax = fscore;
                    xfilmax2 = xfil2 + 1;
                    ycolmax2 = ycol2 + 1;
                    maxIdentities = nIdentities;
                }
                continue;
            }
            if (fscore >= 32 * POINT && fscore + 32 * POINT <= fscoreMax) {
                int m = __builtin_popcount(eq_mask_fwd);
                fscore += (2 * m - 32) * POINT;
                nIdentities += m;
                xfil2 -= 32;
                ycol2 -= 32;
                rem_bwd -= 32;
                continue;
            }
            /* Reverse 32 bytes so index 0 corresponds to xfil2 and index 31 to xfil2 - 31 */
            __m256i eq_rev = _mm256_permute2x128_si256(eq, eq, 0x01);
            eq_rev = _mm256_shuffle_epi8(eq_rev, v_rev16);

            __m256i s = _mm256_add_epi8(_mm256_and_si256(eq_rev, v_two), v_neg1);
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 1));
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 2));
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 4));
            s = _mm256_add_epi8(s, _mm256_slli_si256(s, 8));
            __m128i lo_lane = _mm256_castsi256_si128(s);
            __m128i lo_sum = _mm_shuffle_epi8(lo_lane, v_bcast15);
            __m256i psum = _mm256_add_epi8(s, _mm256_inserti128_si256(v_zero, lo_sum, 1));

            __m128i p_lo = _mm256_castsi256_si128(psum);
            __m128i p_hi = _mm256_extracti128_si256(psum, 1);
            __m128i mn = _mm_min_epi8(p_lo, p_hi);
            __m128i mx = _mm_max_epi8(p_lo, p_hi);
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 8));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 8));
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 4));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 4));
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 2));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 2));
            mn = _mm_min_epi8(mn, _mm_srli_si128(mn, 1));
            mx = _mm_max_epi8(mx, _mm_srli_si128(mx, 1));
            int chunk_min = (int)(int8_t)_mm_extract_epi8(mn, 0);
            int chunk_max = (int)(int8_t)_mm_extract_epi8(mx, 0);

            if (__builtin_expect(fscore + POINT * chunk_min < 0, 0)) {
                uint32_t eq_mask = (uint32_t)_mm256_movemask_epi8(eq_rev);
                for (int k = 0; k < 32; k++) {
                    if ((eq_mask >> k) & 1U) {
                        fscore += POINT;
                        nIdentities++;
                        if (fscoreMax < fscore) {
                            fscoreMax = fscore;
                            xfilmax2 = xfil2 - k;
                            ycolmax2 = ycol2 - k;
                            maxIdentities = nIdentities;
                        }
                    } else {
                        fscore -= POINT;
                        if (fscore < 0) goto bwd_done;
                    }
                }
            }

            if (fscore + POINT * chunk_max > fscoreMax) {
                uint32_t max_mask = (uint32_t)_mm256_movemask_epi8(
                    _mm256_cmpeq_epi8(psum, _mm256_set1_epi8((char)chunk_max)));
                int k = __builtin_ctz(max_mask);
                fscoreMax = fscore + POINT * chunk_max;
                xfilmax2 = xfil2 - k;
                ycolmax2 = ycol2 - k;
                maxIdentities = nIdentities + ((k + 1 + chunk_max) >> 1);
            }

            int total_delta = (int)(int8_t)_mm_extract_epi8(p_hi, 15);
            fscore += POINT * total_delta;
            nIdentities += (32 + total_delta) >> 1;
            xfil2 -= 32;
            ycol2 -= 32;
            rem_bwd -= 32;
            continue;
        }
        for (int k = 0; k < 32; k++) {
            unsigned char vx_c = dx[xfil2] & 0xDF;
            unsigned char vy_c = dy[ycol2] & 0xDF;
            if (vx_c == 0 || vy_c == 0) goto bwd_done;
            if (vx_c == 'N' || vy_c == 'N') {
                fscore -= 1;
            } else if (vx_c == vy_c) {
                fscore += POINT;
                nIdentities++;
                if (fscoreMax < fscore) {
                    fscoreMax = fscore;
                    xfilmax2 = xfil2;
                    ycolmax2 = ycol2;
                    maxIdentities = nIdentities;
                }
            } else {
                fscore -= POINT;
            }
            xfil2--;
            ycol2--;
            rem_bwd--;
            if (fscore < 0) goto bwd_done;
        }
    }
#endif

    /* Backward scalar tail */
    while (rem_bwd > 0) {
        unsigned char vx_c = dx[xfil2] & 0xDF;
        unsigned char vy_c = dy[ycol2] & 0xDF;
        if (vx_c == 0 || vy_c == 0) break;
        if (vx_c == 'N' || vy_c == 'N') {
            fscore -= 1;
        } else if (vx_c == vy_c) {
            fscore += POINT;
            nIdentities++;
            if (fscoreMax < fscore) {
                fscoreMax = fscore;
                xfilmax2 = xfil2;
                ycolmax2 = ycol2;
                maxIdentities = nIdentities;
            }
        } else {
            fscore -= POINT;
        }
        xfil2--;
        ycol2--;
        rem_bwd--;
        if (fscore < 0) break;
    }

bwd_done:;
    uint64_t length = (uint64_t)(xfilmax - xfilmax2 + 1);
    float similarity = (float)((double)fscoreMax * 100.0 / (double)(length * POINT));

    int len_bin = (length < 999) ? (int)length : 999;
    int sim_bin = (int)similarity;
    if (sim_bin < 0) sim_bin = 0;
    if (sim_bin > 99) sim_bin = 99;
    coverage[len_bin][sim_bin]++;

    *out_xfilmax = xfilmax;

    if (length > min_len && similarity > (float)sim_th) {
        uint64_t seqY_id_in_buf = find_seq_id(sY, posY);
        out_frag->diag = (int64_t)posX - (int64_t)posY;
        out_frag->xStart = (uint64_t)xfilmax2 - seqX_id;
        out_frag->yStart = (uint64_t)ycolmax2 - seqY_id_in_buf;
        out_frag->xEnd = (uint64_t)xfilmax - seqX_id;
        out_frag->yEnd = (uint64_t)ycolmax - seqY_id_in_buf;
        out_frag->length = length;
        out_frag->ident = (uint64_t)maxIdentities;
        out_frag->score = (uint64_t)fscoreMax;
        out_frag->similarity = similarity;
        out_frag->seqX = seqX_id;
        out_frag->seqY = (strand == 'f') ? seqY_id_in_buf : (sY->n_seqs - seqY_id_in_buf - 1);
        out_frag->block = 0;
        out_frag->strand = strand;
        return 1;
    }
    return 0;
}

/* ============================================================================
 * Multi-Threaded Streaming Query Worker with 64-Wide Software Prefetching
 * ============================================================================ */
typedef struct {
    uint64_t posX;
    uint64_t seqX;
    uint64_t key;
    uint32_t b_idx;
    uint32_t start;
    uint32_t end;
    uint8_t  x_strand;
    uint8_t  is_pal;
} QueryItem;

typedef struct {
    const FastSeq *sX;
    const FastSeq *sY;
    const FastSeq *sY_rc;
    const RefIndex *idx;
    int wl;
    uint64_t min_len;
    uint64_t sim_th;
    uint64_t chunk_size;
    uint64_t num_chunks;
    uint64_t *next_chunk;
} QueryShared;

typedef struct {
    const QueryShared *sh;
    FragVec frags_f;
    FragVec frags_r;
    uint64_t hits_seen;
    uint64_t hits_extended;
    long coverage[1000][100];
} QueryThreadResult;

/*
 * Exact rejected-span skip bound:
 * Given a rejected extension from posX reaching forward peak xfilmax, partition
 * [posX .. xfilmax] into blocks of size wl. For any later seed p in block b,
 * p - 1 lies in block <= b and every t >= p + wl lies in block >= b + 1.
 * Whenever max_{u <= b} S(u) <= min_{v >= b + 1} S(v), every seed p in block b
 * is mathematically guaranteed to reach the exact same forward peak xfilmax and
 * identical backward peak xfilmax2, producing the exact same rejected span.
 */
static int64_t compute_exact_reject_skip(
    const FastSeq *sX,
    const FastSeq *sY,
    uint64_t posX,
    uint64_t posY,
    int64_t xfilmax,
    int wl
) {
    int64_t span = xfilmax - (int64_t)posX + 1;
    size_t n_blk = (size_t)((span + wl - 1) / wl);
    if (n_blk < 3) return (int64_t)posX + (int64_t)(2 * wl - 1);

    int32_t *blk_min = (int32_t *)malloc(n_blk * sizeof(int32_t));
    int32_t *blk_max = (int32_t *)malloc(n_blk * sizeof(int32_t));
    if (!blk_min || !blk_max) {
        free(blk_min);
        free(blk_max);
        return (int64_t)posX + (int64_t)(2 * wl - 1);
    }

    const unsigned char *dx = sX->data;
    const unsigned char *dy = sY->data;
    int32_t sc = 0;

    for (size_t b = 0; b < n_blk; b++) {
        int64_t b_start = (int64_t)b * wl;
        int64_t b_end = b_start + wl;
        if (b_end > span) b_end = span;
        int32_t bmin = INT32_MAX;
        int32_t bmax = INT32_MIN;
        for (int64_t off = b_start; off < b_end; off++) {
            unsigned char cx = dx[posX + (uint64_t)off] & 0xDF;
            unsigned char cy = dy[posY + (uint64_t)off] & 0xDF;
            if (cx == 'N' || cy == 'N') sc -= 1;
            else if (cx == cy) sc += POINT;
            else sc -= POINT;
            if (sc < bmin) bmin = sc;
            if (sc > bmax) bmax = sc;
        }
        blk_min[b] = bmin;
        blk_max[b] = bmax;
    }

    /* Suffix minimum over blk_min[b .. n_blk - 1] in-place */
    for (size_t b = n_blk - 1; b > 0; b--) {
        if (blk_min[b] < blk_min[b - 1]) {
            blk_min[b - 1] = blk_min[b];
        }
    }

    /* Find largest contiguous prefix of blocks satisfying running_max <= suf_min[b + 2].
     * For any exact wl-mer starting at p in block b, [p .. p+wl-1] adds +4*wl and there
     * are at most (wl - 1) bases remaining in block b + 1, so S(t) > S(p - 1) automatically
     * holds throughout block b + 1; checking suf_min[b + 2] covers all t in blocks >= b + 2. */
    int32_t run_max = 0;
    size_t safe_b = 0;
    for (size_t b = 0; b + 2 < n_blk; b++) {
        if (blk_max[b] > run_max) run_max = blk_max[b];
        if (run_max <= blk_min[b + 2]) {
            safe_b = b + 1; /* all positions in block b are provably safe to skip */
        } else {
            break;
        }
    }

    free(blk_min);
    free(blk_max);

    int64_t skip_end = (int64_t)posX + (int64_t)safe_b * wl - 1;
    int64_t min_skip = (int64_t)posX + (int64_t)(wl / 2);
    return (skip_end > min_skip) ? skip_end : min_skip;
}

static inline void handle_match_hit(
    const QueryShared *sh,
    QueryThreadResult *res,
    StrandDiagCache *cache_f,
    StrandDiagCache *cache_r,
    uint64_t posX,
    uint64_t seqX,
    uint32_t posY_orig,
    int is_rev
) {
    res->hits_seen++;
    int wl = sh->wl;
    int64_t noise_end = (int64_t)posX + (int64_t)(2 * wl - 1);
    struct FragFile frag;
    int64_t xfilmax;

    if (!is_rev) {
        uint64_t posY = (uint64_t)posY_orig;
        int64_t diag = (int64_t)posX - (int64_t)posY;
        if (diag_is_covered(cache_f, diag, (int64_t)posX)) return;

        res->hits_extended++;
        int accepted = extend_hit(sh->sX, sh->sY, posX, posY, seqX, wl,
                                  sh->min_len, sh->sim_th, 'f', &frag, &xfilmax, res->coverage);
        if (accepted) {
            frag_vec_push(&res->frags_f, &frag);
            diag_record_hsp(cache_f, diag, xfilmax);
        } else if (xfilmax - (int64_t)posX > 256) {
            int64_t skip_x = compute_exact_reject_skip(sh->sX, sh->sY, posX, posY, xfilmax, wl);
            diag_record_hsp(cache_f, diag, skip_x);
        } else {
            diag_record_noise(cache_f, diag, noise_end);
        }
    } else {
        uint64_t posY_rc = sh->sY->buf_len - (uint64_t)posY_orig - (uint64_t)wl;
        int64_t diag = (int64_t)posX - (int64_t)posY_rc;
        if (diag_is_covered(cache_r, diag, (int64_t)posX)) return;

        res->hits_extended++;
        int accepted = extend_hit(sh->sX, sh->sY_rc, posX, posY_rc, seqX, wl,
                                  sh->min_len, sh->sim_th, 'r', &frag, &xfilmax, res->coverage);
        if (accepted) {
            frag_vec_push(&res->frags_r, &frag);
            diag_record_hsp(cache_r, diag, xfilmax);
        } else if (xfilmax - (int64_t)posX > 256) {
            int64_t skip_x = compute_exact_reject_skip(sh->sX, sh->sY_rc, posX, posY_rc, xfilmax, wl);
            diag_record_hsp(cache_r, diag, skip_x);
        } else {
            diag_record_noise(cache_r, diag, noise_end);
        }
    }
}

static void process_batch_direct(
    const QueryShared *sh,
    QueryThreadResult *res,
    StrandDiagCache *cache_f,
    StrandDiagCache *cache_r,
    QueryItem *batch,
    int count
) {
    const uint32_t *offsets = sh->idx->offsets;
    const uint32_t *pos_strand = sh->idx->pos_strand;

    /* Stage 1: prefetch offset table entries */
    for (int i = 0; i < count; i++) {
        __builtin_prefetch(&offsets[batch[i].key], 0, 1);
    }

    /* Stage 2: fetch ranges and prefetch non-empty posting lists */
    for (int i = 0; i < count; i++) {
        uint32_t st = offsets[batch[i].key];
        uint32_t en = offsets[batch[i].key + 1];
        batch[i].start = st;
        batch[i].end = en;
        if (st < en) {
            __builtin_prefetch(&pos_strand[st], 0, 1);
        }
    }

    /* Stage 3: process hits */
    for (int i = 0; i < count; i++) {
        uint32_t st = batch[i].start;
        uint32_t en = batch[i].end;
        if (st == en) continue;

        uint64_t posX = batch[i].posX;
        uint64_t seqX = batch[i].seqX;
        uint8_t x_strand = batch[i].x_strand;
        uint8_t is_pal = batch[i].is_pal;

        for (uint32_t j = st; j < en; j++) {
            uint32_t ps = pos_strand[j];
            uint32_t posY = ps >> 1;
            uint8_t y_strand = (uint8_t)(ps & 1U);
            if (__builtin_expect(is_pal, 0)) {
                handle_match_hit(sh, res, cache_f, cache_r, posX, seqX, posY, 0);
                handle_match_hit(sh, res, cache_f, cache_r, posX, seqX, posY, 1);
            } else {
                handle_match_hit(sh, res, cache_f, cache_r, posX, seqX, posY, (x_strand != y_strand));
            }
        }
    }
}

static void process_batch_hash_csr(
    const QueryShared *sh,
    QueryThreadResult *res,
    StrandDiagCache *cache_f,
    StrandDiagCache *cache_r,
    QueryItem *batch,
    int count
) {
    const RefIndex *idx = sh->idx;
    const uint32_t *offsets = idx->offsets;
    const uint64_t *keys = idx->keys;
    const uint32_t *pos_strand = idx->pos_strand;

    /* Stage 1: prefetch bucket offset table */
    for (int i = 0; i < count; i++) {
        __builtin_prefetch(&offsets[batch[i].b_idx], 0, 1);
    }

    /* Stage 2: fetch bucket ranges and prefetch key cache lines */
    for (int i = 0; i < count; i++) {
        uint32_t b_idx = batch[i].b_idx;
        uint32_t st = offsets[b_idx];
        uint32_t en = offsets[b_idx + 1];
        batch[i].start = st;
        batch[i].end = en;
        if (st < en) {
            __builtin_prefetch(&keys[st], 0, 1);
        }
    }

    /* Stage 3: probe buckets (linear/AVX2 for small buckets, binary search for large) */
    for (int i = 0; i < count; i++) {
        uint32_t st = batch[i].start;
        uint32_t en = batch[i].end;
        if (st == en) continue;

        uint64_t key = batch[i].key;
        uint32_t len = en - st;
        uint32_t j = st;

        if (len <= 16) {
#if defined(__AVX2__)
            __m256i vkey = _mm256_set1_epi64x((long long)key);
            while (j + 4 <= en) {
                __m256i vcand = _mm256_loadu_si256((const __m256i *)&keys[j]);
                int mask = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(vcand, vkey)));
                if (mask) {
                    j += (uint32_t)__builtin_ctz(mask);
                    goto found_key;
                }
                if (keys[j + 3] > key) goto next_item;
                j += 4;
            }
#endif
            while (j < en && keys[j] < key) j++;
            if (j >= en || keys[j] != key) continue;
        } else {
            /* Binary search lower_bound */
            uint32_t lo = st, hi = en;
            while (lo < hi) {
                uint32_t mid = lo + ((hi - lo) >> 1);
                if (keys[mid] < key) lo = mid + 1;
                else hi = mid;
            }
            j = lo;
            if (j >= en || keys[j] != key) continue;
        }

#if defined(__AVX2__)
found_key:;
#endif
        uint64_t posX = batch[i].posX;
        uint64_t seqX = batch[i].seqX;
        uint8_t x_strand = batch[i].x_strand;
        uint8_t is_pal = batch[i].is_pal;

        do {
            uint32_t ps = pos_strand[j];
            uint32_t posY = ps >> 1;
            uint8_t y_strand = (uint8_t)(ps & 1U);
            if (__builtin_expect(is_pal, 0)) {
                handle_match_hit(sh, res, cache_f, cache_r, posX, seqX, posY, 0);
                handle_match_hit(sh, res, cache_f, cache_r, posX, seqX, posY, 1);
            } else {
                handle_match_hit(sh, res, cache_f, cache_r, posX, seqX, posY, (x_strand != y_strand));
            }
            j++;
        } while (j < en && keys[j] == key);

#if defined(__AVX2__)
next_item:;
#endif
    }
}

static void *query_worker(void *arg) {
    QueryThreadResult *res = (QueryThreadResult *)arg;
    const QueryShared *sh = res->sh;
    const FastSeq *sX = sh->sX;
    const RefIndex *idx = sh->idx;
    const unsigned char *dataX = sX->data;

    int wl = sh->wl;
    uint64_t mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rc_shift = 2 * (wl - 1);
    int direct = idx->direct_mode;
    int shift_bucket = 64 - idx->bucket_bits;
    int sub_bits = idx->sub_bits;

    StrandDiagCache *cache_f = (StrandDiagCache *)malloc(sizeof(StrandDiagCache));
    StrandDiagCache *cache_r = (StrandDiagCache *)malloc(sizeof(StrandDiagCache));
    if (!cache_f || !cache_r) terror("OOM diag cache");
    init_diag_cache(cache_f);
    init_diag_cache(cache_r);

    QueryItem batch[PREFETCH_BATCH];

    while (1) {
        uint64_t cid = __atomic_fetch_add(sh->next_chunk, 1ULL, __ATOMIC_RELAXED);
        if (cid >= sh->num_chunks) break;

        uint64_t c_start = cid * sh->chunk_size;
        uint64_t c_end = c_start + sh->chunk_size;
        if (c_end > sX->buf_len) c_end = sX->buf_len;

        uint64_t scan_from = (c_start >= (uint64_t)(wl - 1)) ? (c_start - (uint64_t)(wl - 1)) : 0;
        uint64_t cur_seqX = find_seq_id(sX, scan_from);

        uint64_t kf = 0, kr = 0;
        int valid = 0;
        int bcount = 0;

        for (uint64_t i = scan_from; i < c_end; i++) {
            unsigned char ch = dataX[i];
            if (ch == 0x00) {
                cur_seqX++;
                valid = 0;
                kf = 0;
                kr = 0;
                continue;
            }
            uint8_t b = BASE_TO_2BIT[ch];
            if (b < 4) {
                kf = ((kf << 2) | b) & mask;
                kr = (kr >> 2) | ((uint64_t)(b ^ 3) << rc_shift);
                if (++valid >= wl && i >= c_start) {
                    uint8_t x_strand = (kf <= kr) ? 0U : 1U;
                    uint64_t key = x_strand ? kr : kf;
                    QueryItem *qi = &batch[bcount++];
                    qi->posX = i + 1 - (uint64_t)wl;
                    qi->seqX = cur_seqX;
                    qi->key = key;
                    qi->x_strand = x_strand;
                    qi->is_pal = (uint8_t)(kf == kr);
                    if (!direct) {
                        uint32_t b_num = (uint32_t)(hash64(key) >> shift_bucket);
                        qi->b_idx = b_num + (b_num >> sub_bits);
                    }
                    if (bcount == PREFETCH_BATCH) {
                        if (direct) process_batch_direct(sh, res, cache_f, cache_r, batch, bcount);
                        else process_batch_hash_csr(sh, res, cache_f, cache_r, batch, bcount);
                        bcount = 0;
                    }
                }
            } else {
                valid = 0;
                kf = 0;
                kr = 0;
            }
        }
        if (bcount > 0) {
            if (direct) process_batch_direct(sh, res, cache_f, cache_r, batch, bcount);
            else process_batch_hash_csr(sh, res, cache_f, cache_r, batch, bcount);
        }
    }

    free(cache_f);
    free(cache_r);
    return NULL;
}

/* ============================================================================
 * Fragment Sorting & Exact Deduplication (handles chunk boundaries & collisions)
 * ============================================================================ */
static int cmp_frags_dedup(const void *a, const void *b) {
    const struct FragFile *fa = (const struct FragFile *)a;
    const struct FragFile *fb = (const struct FragFile *)b;
    if (fa->diag != fb->diag) return (fa->diag < fb->diag) ? -1 : 1;
    if (fa->xStart != fb->xStart) return (fa->xStart < fb->xStart) ? -1 : 1;
    if (fa->xEnd != fb->xEnd) return (fa->xEnd > fb->xEnd) ? -1 : 1; /* longer first */
    if (fa->score != fb->score) return (fa->score > fb->score) ? -1 : 1;
    return 0;
}

static size_t dedup_frags(struct FragFile *arr, size_t n) {
    if (n <= 1) return n;
    qsort(arr, n, sizeof(struct FragFile), cmp_frags_dedup);
    size_t out = 0;
    for (size_t i = 0; i < n; i++) {
        if (out > 0 &&
            arr[i].diag == arr[out - 1].diag &&
            arr[i].seqX == arr[out - 1].seqX &&
            arr[i].seqY == arr[out - 1].seqY &&
            arr[i].xEnd <= arr[out - 1].xEnd) {
            /* Fully contained within previous fragment on the exact same diagonal */
            continue;
        }
        arr[out++] = arr[i];
    }
    return out;
}

/* ============================================================================
 * Fast Buffered Output Writers (100% GECKO-compatible .frags, .INF, .MAT, .csv)
 * ============================================================================ */
static inline void put_be64(uint8_t *dst, uint64_t v) {
    uint64_t be = __builtin_bswap64(v);
    memcpy(dst, &be, 8);
}

static inline void put_be32(uint8_t *dst, uint32_t v) {
    uint32_t be = __builtin_bswap32(v);
    memcpy(dst, &be, 4);
}

static inline void pack_frag_be93(uint8_t *dst, const struct FragFile *f) {
    put_be64(dst +  0, (uint64_t)f->diag);
    put_be64(dst +  8, f->xStart);
    put_be64(dst + 16, f->yStart);
    put_be64(dst + 24, f->xEnd);
    put_be64(dst + 32, f->yEnd);
    put_be64(dst + 40, f->length);
    put_be64(dst + 48, f->ident);
    put_be64(dst + 56, f->score);
    uint32_t sim_bits;
    memcpy(&sim_bits, &f->similarity, 4);
    put_be32(dst + 64, sim_bits);
    put_be64(dst + 68, f->seqX);
    put_be64(dst + 76, f->seqY);
    put_be64(dst + 84, (uint64_t)f->block);
    dst[92] = (uint8_t)f->strand;
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

static inline char *format_csv_frag(char *p, const struct FragFile *f, uint64_t ys, uint64_t ye) {
    double pct_ident = 100.0 * ((double)f->ident / (double)f->length);
    memcpy(p, "Frag,", 5); p += 5;
    p = append_u64(p, f->xStart); *p++ = ',';
    p = append_u64(p, ys); *p++ = ',';
    p = append_u64(p, f->xEnd); *p++ = ',';
    p = append_u64(p, ye); *p++ = ',';
    *p++ = f->strand; *p++ = ',';
    p = append_u64(p, (uint64_t)f->block); *p++ = ',';
    p = append_u64(p, f->length); *p++ = ',';
    p = append_u64(p, f->score); *p++ = ',';
    p = append_u64(p, f->ident); *p++ = ',';
    p = append_fixed2(p, (double)f->similarity); *p++ = ',';
    p = append_fixed2(p, pct_ident); *p++ = ',';
    p = append_u64(p, f->seqX); *p++ = ',';
    p = append_u64(p, f->seqY); *p++ = '\n';
    return p;
}

static void write_outputs(
    const char *out_prefix,
    const char *seqX_path,
    const char *seqY_path,
    const FastSeq *sX,
    const FastSeq *sY,
    uint64_t Lmin,
    uint64_t SimTh,
    uint64_t tot_hits,
    uint64_t tot_hits_used,
    const struct FragFile *frags_f,
    size_t nf_f,
    const struct FragFile *frags_r,
    size_t nf_r,
    long coverage[1000][100]
) {
    char path[1024];
    uint64_t n0 = sX->n_bases;
    uint64_t n1 = sY->n_bases;
    const size_t BUF_CAP = 4 * 1024 * 1024;
    uint8_t *io_buf = (uint8_t *)malloc(BUF_CAP);
    if (!io_buf) terror("OOM output io_buf");

    /* 1. Fast buffered binary .frags file (GECKO big-endian 93-byte records) */
    snprintf(path, sizeof(path), "%s.frags", out_prefix);
    FILE *fOut = fopen(path, "wb");
    if (!fOut) terror("opening output .frags file");
    put_be64(io_buf, n0);
    put_be64(io_buf + 8, n1);
    size_t bpos = 16;
    for (size_t i = 0; i < nf_f; i++) {
        if (bpos + 93 > BUF_CAP) { fwrite(io_buf, 1, bpos, fOut); bpos = 0; }
        pack_frag_be93(io_buf + bpos, &frags_f[i]);
        bpos += 93;
    }
    for (size_t i = 0; i < nf_r; i++) {
        if (bpos + 93 > BUF_CAP) { fwrite(io_buf, 1, bpos, fOut); bpos = 0; }
        pack_frag_be93(io_buf + bpos, &frags_r[i]);
        bpos += 93;
    }
    if (bpos > 0) fwrite(io_buf, 1, bpos, fOut);
    fclose(fOut);

    /* 2. Fast buffered .frags.MAT file */
    snprintf(path, sizeof(path), "%s.frags.MAT", out_prefix);
    FILE *fMat = fopen(path, "wb");
    if (!fMat) terror("opening output .frags.MAT file");
    char *cp = (char *)io_buf;
    char *cend = cp + BUF_CAP - 64;
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 100; j++) {
            if (cp >= cend) { fwrite(io_buf, 1, (size_t)(cp - (char *)io_buf), fMat); cp = (char *)io_buf; }
            cp = append_u64(cp, (uint64_t)coverage[i][j]);
            *cp++ = '\t';
        }
        *cp++ = '\n';
    }
    if (cp > (char *)io_buf) fwrite(io_buf, 1, (size_t)(cp - (char *)io_buf), fMat);
    fclose(fMat);

    /* 3. .frags.INF file */
    char inf_hdr[2048];
    int hdr_len = snprintf(inf_hdr, sizeof(inf_hdr),
        "All by-Identity Ungapped Fragments (gecko_fast in-memory engine)\n"
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
        "Total fragments      : %" PRIu64 "\n"
        "========================================================\n",
        seqX_path, seqY_path, sX->ident, sY->ident, n0, n1, Lmin, SimTh,
        tot_hits, tot_hits_used, (uint64_t)(nf_f + nf_r));

    snprintf(path, sizeof(path), "%s.frags.INF", out_prefix);
    FILE *fInf = fopen(path, "wb");
    if (!fInf) terror("opening output .frags.INF file");
    fwrite(inf_hdr, 1, (size_t)hdr_len, fInf);
    fclose(fInf);

    /* 4. Fast buffered .csv file */
    snprintf(path, sizeof(path), "%s.csv", out_prefix);
    FILE *fCsv = fopen(path, "wb");
    if (!fCsv) terror("opening output .csv file");
    fwrite(inf_hdr, 1, (size_t)hdr_len, fCsv);
    static const char csv_cols[] =
        "Total CSB: 0\n"
        "========================================================\n"
        "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY\n";
    fwrite(csv_cols, 1, sizeof(csv_cols) - 1, fCsv);

    cp = (char *)io_buf;
    cend = cp + BUF_CAP - 512;
    for (size_t i = 0; i < nf_f; i++) {
        if (cp >= cend) { fwrite(io_buf, 1, (size_t)(cp - (char *)io_buf), fCsv); cp = (char *)io_buf; }
        cp = format_csv_frag(cp, &frags_f[i], frags_f[i].yStart, frags_f[i].yEnd);
    }
    for (size_t i = 0; i < nf_r; i++) {
        if (cp >= cend) { fwrite(io_buf, 1, (size_t)(cp - (char *)io_buf), fCsv); cp = (char *)io_buf; }
        cp = format_csv_frag(cp, &frags_r[i], n1 - frags_r[i].yStart - 1, n1 - frags_r[i].yEnd - 1);
    }
    if (cp > (char *)io_buf) fwrite(io_buf, 1, (size_t)(cp - (char *)io_buf), fCsv);
    fclose(fCsv);
    free(io_buf);
}

typedef struct {
    const char *path;
    FastSeq *out;
} LoadTask;

static void *load_fasta_thread(void *arg) {
    LoadTask *t = (LoadTask *)arg;
    t->out = load_fasta(t->path);
    return NULL;
}

/* ============================================================================
 * Main
 * ============================================================================ */
static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <query.fasta> <ref.fasta> <out_prefix> <length> <similarity> <WL> <fixedL> [threads] [freqLimit]\n"
        "\n"
        "  query.fasta  Query FASTA file\n"
        "  ref.fasta    Reference FASTA file\n"
        "  out_prefix   Output path prefix (creates <prefix>.frags, .frags.INF, .frags.MAT, .csv)\n"
        "  length       Minimum HSP length (bp if fixedL=1, percentage if fixedL=0)\n"
        "  similarity   Minimum similarity threshold (0-100)\n"
        "  WL           Seed word length (8..32)\n"
        "  fixedL       1 for fixed bp length (standard), 0 for percentage\n"
        "  threads      Number of worker threads (default: all available cores)\n"
        "  freqLimit    Per-strand k-mer frequency cap in reference (default: 1000, 0 = unlimited)\n",
        prog);
}

int main(int argc, char **argv) {
    if (argc < 8 || argc > 10) {
        print_usage(argv[0]);
        return 1;
    }

    const char *seqX_path = argv[1];
    const char *seqY_path = argv[2];
    const char *out_prefix = argv[3];
    uint64_t Lmin = (uint64_t)strtoull(argv[4], NULL, 10);
    uint64_t SimTh = (uint64_t)strtoull(argv[5], NULL, 10);
    int WL = atoi(argv[6]);
    int fixedL = atoi(argv[7]);

    int n_threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (n_threads < 1) n_threads = 1;
    if (argc >= 9) {
        int req = atoi(argv[8]);
        if (req > 0) n_threads = req;
    }

    uint64_t freq_limit = 1000;
    if (argc >= 10) {
        freq_limit = (uint64_t)strtoull(argv[9], NULL, 10);
    }

    if (WL < 8 || WL > 32) {
        fprintf(stderr, "Error: WL must be between 8 and 32 (got %d)\n", WL);
        return 1;
    }

    init_base_table();
    double t0 = now_sec();

    /* 1. Load seqX and seqY concurrently on 2 threads, then build in-memory revcomp of Y */
    LoadTask load_x_task = { .path = seqX_path, .out = NULL };
    pthread_t load_x_th;
    pthread_create(&load_x_th, NULL, load_fasta_thread, &load_x_task);

    FastSeq *sY = load_fasta(seqY_path);
    if (sY->buf_len >= (1ULL << 31)) {
        pthread_join(load_x_th, NULL);
        fprintf(stderr, "Error: reference sequence exceeds 2.14 Gb limit for 32-bit position packing.\n");
        return 1;
    }
    FastSeq *sY_rc = build_revcomp_seq(sY);
    pthread_join(load_x_th, NULL);
    FastSeq *sX = load_x_task.out;

    double t_load = now_sec();

    /* Compute effective minimum length */
    uint64_t min_len = Lmin;
    if (!fixedL) {
        uint64_t min_seq = (sX->n_bases < sY->n_bases) ? sX->n_bases : sY->n_bases;
        min_len = (uint64_t)((double)min_seq * ((double)Lmin / 100.0));
    }

    /* 2. Build in-memory canonical k-mer index for Reference Y */
    RefIndex *idx = build_ref_index(sY, WL, freq_limit, n_threads);
    double t_index = now_sec();

    /* 3. Parallel streaming query scan + on-the-fly diagonal extension */
    uint64_t chunk_size = 1024 * 1024; /* 1 Mb chunks for fine-grained load balancing */
    if (sX->buf_len / (uint64_t)(n_threads * 8) < chunk_size) {
        chunk_size = sX->buf_len / (uint64_t)(n_threads * 8);
        if (chunk_size < 65536) chunk_size = 65536;
    }
    uint64_t num_chunks = (sX->buf_len + chunk_size - 1) / chunk_size;
    uint64_t next_chunk = 0;

    QueryShared qsh = {
        .sX = sX,
        .sY = sY,
        .sY_rc = sY_rc,
        .idx = idx,
        .wl = WL,
        .min_len = min_len,
        .sim_th = SimTh,
        .chunk_size = chunk_size,
        .num_chunks = num_chunks,
        .next_chunk = &next_chunk
    };

    QueryThreadResult *t_res = (QueryThreadResult *)calloc((size_t)n_threads, sizeof(QueryThreadResult));
    pthread_t *ths = (pthread_t *)malloc((size_t)n_threads * sizeof(pthread_t));
    if (!t_res || !ths) terror("OOM query threads");

    for (int t = 0; t < n_threads; t++) {
        t_res[t].sh = &qsh;
        pthread_create(&ths[t], NULL, query_worker, &t_res[t]);
    }
    for (int t = 0; t < n_threads; t++) {
        pthread_join(ths[t], NULL);
    }
    double t_query = now_sec();

    /* 4. Merge thread results and deduplicate */
    size_t tot_f = 0, tot_r = 0;
    uint64_t tot_hits = 0, tot_hits_used = 0;
    long coverage[1000][100];
    memset(coverage, 0, sizeof(coverage));

    for (int t = 0; t < n_threads; t++) {
        tot_f += t_res[t].frags_f.count;
        tot_r += t_res[t].frags_r.count;
        tot_hits += t_res[t].hits_seen;
        tot_hits_used += t_res[t].hits_extended;
        for (int i = 0; i < 1000; i++) {
            for (int j = 0; j < 100; j++) {
                coverage[i][j] += t_res[t].coverage[i][j];
            }
        }
    }

    struct FragFile *all_f = (struct FragFile *)malloc((tot_f > 0 ? tot_f : 1) * sizeof(struct FragFile));
    struct FragFile *all_r = (struct FragFile *)malloc((tot_r > 0 ? tot_r : 1) * sizeof(struct FragFile));
    if (!all_f || !all_r) terror("OOM merge frags");

    size_t off_f = 0, off_r = 0;
    for (int t = 0; t < n_threads; t++) {
        if (t_res[t].frags_f.count) {
            memcpy(&all_f[off_f], t_res[t].frags_f.items, t_res[t].frags_f.count * sizeof(struct FragFile));
            off_f += t_res[t].frags_f.count;
        }
        if (t_res[t].frags_r.count) {
            memcpy(&all_r[off_r], t_res[t].frags_r.items, t_res[t].frags_r.count * sizeof(struct FragFile));
            off_r += t_res[t].frags_r.count;
        }
        free(t_res[t].frags_f.items);
        free(t_res[t].frags_r.items);
    }

    size_t nf_f = dedup_frags(all_f, tot_f);
    size_t nf_r = dedup_frags(all_r, tot_r);

    /* 5. Write GECKO-compatible outputs */
    write_outputs(out_prefix, seqX_path, seqY_path, sX, sY, Lmin, SimTh,
                  tot_hits, tot_hits_used, all_f, nf_f, all_r, nf_r, coverage);

    double t_end = now_sec();

    fprintf(stdout,
        "gecko_fast summary:\n"
        "  SeqX: %" PRIu64 " bp (%" PRIu64 " seqs) | SeqY: %" PRIu64 " bp (%" PRIu64 " seqs)\n"
        "  Index: %s (WL=%d, FL=%" PRIu64 ", threads=%d)\n"
        "  Hits matched: %" PRIu64 " | Extended: %" PRIu64 " | Fragments: %zu (F:%zu R:%zu)\n"
        "  Timing: load=%.3fs  index=%.3fs  query+ext=%.3fs  write=%.3fs  TOTAL=%.3fs\n",
        sX->n_bases, sX->n_seqs, sY->n_bases, sY->n_seqs,
        idx->direct_mode ? "direct-CSR" : "hash-CSR", WL, freq_limit, n_threads,
        tot_hits, tot_hits_used, nf_f + nf_r, nf_f, nf_r,
        t_load - t0, t_index - t_load, t_query - t_index, t_end - t_query, t_end - t0);

    free(all_f);
    free(all_r);
    free(t_res);
    free(ths);
    free_ref_index(idx);
    free_fast_seq(sX);
    free_fast_seq(sY);
    free_fast_seq(sY_rc);

    return 0;
}
