#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "quicksortWordForward.h"
#include "quicksortWordReverse.h"

int wordcmp(unsigned char *w1, unsigned char *w2, int n) {
    int i = 0, limit;

    if (n % 4 != 0) {
        w1[n / 4] = w1[n / 4] >> (2 * (3 - ((n - 1) % 4)));
        w1[n / 4] = w1[n / 4] << (2 * (3 - ((n - 1) % 4)));
        w2[n / 4] = w2[n / 4] >> (2 * (3 - ((n - 1) % 4)));
        w2[n / 4] = w2[n / 4] << (2 * (3 - ((n - 1) % 4)));
        limit = (n / 4) + 1;
    } else {
        limit = n / 4;
    }

    for (i = 0; i < limit; i++) {
        if (w1[i] < w2[i]) return -1;
        if (w1[i] > w2[i]) return +1;
    }
    return 0;
}

void showWord(word *w, char *ws) {
    char Alf[] = {'A', 'C', 'G', 'T'};
    int i;
    int wsize = 8;
    unsigned char c;
    for (i = 0; i < wsize; i++) {
        c = w->b[i];
        c = c >> 6;
        ws[4 * i] = Alf[(int) c];
        c = w->b[i];
        c = c << 2;
        c = c >> 6;
        ws[4 * i + 1] = Alf[(int) c];
        c = w->b[i];
        c = c << 4;
        c = c >> 6;
        ws[4 * i + 2] = Alf[(int) c];
        c = w->b[i];
        c = c << 6;
        c = c >> 6;
        ws[4 * i + 3] = Alf[(int) c];
    }
}

int GTWF(wentryF a1, wentryF a2) {
    if (a1.w.b[0] > a2.w.b[0]) return 1;
    else if (a1.w.b[0] < a2.w.b[0]) return 0;

    if (a1.w.b[1] > a2.w.b[1]) return 1;
    else if (a1.w.b[1] < a2.w.b[1]) return 0;

    if (a1.w.b[2] > a2.w.b[2]) return 1;
    else if (a1.w.b[2] < a2.w.b[2]) return 0;

    if (a1.w.b[3] > a2.w.b[3]) return 1;
    else if (a1.w.b[3] < a2.w.b[3]) return 0;

    if (a1.w.b[4] > a2.w.b[4]) return 1;
    else if (a1.w.b[4] < a2.w.b[4]) return 0;

    if (a1.w.b[5] > a2.w.b[5]) return 1;
    else if (a1.w.b[5] < a2.w.b[5]) return 0;

    if (a1.w.b[6] > a2.w.b[6]) return 1;
    else if (a1.w.b[6] < a2.w.b[6]) return 0;

    if (a1.w.b[7] > a2.w.b[7]) return 1;
    else if (a1.w.b[7] < a2.w.b[7]) return 0;

    if (a1.seq > a2.seq) return 1;
    else if (a1.seq < a2.seq) return 0;

    if (a1.pos > a2.pos) return 1;
    return 0;
}

int GTWR(wentryR a1, wentryR a2) {
    if (a1.w.b[0] > a2.w.b[0]) return 1;
    else if (a1.w.b[0] < a2.w.b[0]) return 0;

    if (a1.w.b[1] > a2.w.b[1]) return 1;
    else if (a1.w.b[1] < a2.w.b[1]) return 0;

    if (a1.w.b[2] > a2.w.b[2]) return 1;
    else if (a1.w.b[2] < a2.w.b[2]) return 0;

    if (a1.w.b[3] > a2.w.b[3]) return 1;
    else if (a1.w.b[3] < a2.w.b[3]) return 0;

    if (a1.w.b[4] > a2.w.b[4]) return 1;
    else if (a1.w.b[4] < a2.w.b[4]) return 0;

    if (a1.w.b[5] > a2.w.b[5]) return 1;
    else if (a1.w.b[5] < a2.w.b[5]) return 0;

    if (a1.w.b[6] > a2.w.b[6]) return 1;
    else if (a1.w.b[6] < a2.w.b[6]) return 0;

    if (a1.w.b[7] > a2.w.b[7]) return 1;
    else if (a1.w.b[7] < a2.w.b[7]) return 0;

    if (a1.seq > a2.seq) return 1;
    else if (a1.seq < a2.seq) return 0;

    if (a1.pos > a2.pos) return 1;
    return 0;
}

void shift_word(word *w) {
    int i;
    for (i = 0; i < BYTES_IN_WORD - 1; i++) {
        w->b[i] <<= 2;
        w->b[i] |= (w->b[i + 1] >> 6);
    }
    w->b[BYTES_IN_WORD - 1] <<= 2;
}

void *dictionary(void *a) {
    DictionaryArgs *args = (DictionaryArgs *) a;
    FILE *f;
    char c;
    uint64_t SIZE = 0;
    wentryF *words = NULL;
    wentryF *more_words = NULL;
    char *seq = NULL;
    uint64_t i = 0, r = 0;

    if ((f = fopen(args->seqFile, "rt")) == NULL) {
        fprintf(stdout, "opening sequence file: %s\n", args->seqFile);
        terror("opening sequence file");
    }

    fseek(f, 0, SEEK_END);
    SIZE = ftell(f);
    fseek(f, 0, SEEK_SET);

    if (SIZE == 0)
        terror("Empty sequence X file");

    if ((words = calloc(SIZE, sizeof(wentryF))) == NULL) {
        terror("not enough memory for words array");
    }

    if ((seq = calloc(READBUF, sizeof(char))) == NULL) {
        terror("not enough memory for read buffer");
    }

#ifdef VERBOSE
    fprintf(stdout, "Memory allocated: %" PRIu64 "\n", SIZE * sizeof(wentryF));
    fflush(stdout);
#endif

    //To force the read
    i = READBUF + 1;

    c = buffered_fgetc(seq, &i, &r, f);
    while (c != '\n') {
        c = buffered_fgetc(seq, &i, &r, f);
    }

    wentryF WE;
    WE.seq = 0;
    uint64_t index = 0;
    uint64_t inEntry = 0;
    uint64_t NW = 0;
    uint64_t Tot = 0;
    uint64_t NoACGT = 0;
    uint64_t NoC = 0;
    uint64_t nLocs = 0;
    uint64_t loc_size = 0;
    c = buffered_fgetc(seq, &i, &r, f);
    //TODO check the loop conditions
    while (!feof(f)) {
        if (!isupper(toupper(c))) {
            if (c == '>') {
                c = buffered_fgetc(seq, &i, &r, f);
                while (c != '\n')
                    c = buffered_fgetc(seq, &i, &r, f);
                WE.seq++;
                inEntry = 0;
                index++;
            }
            NoC++;
            c = buffered_fgetc(seq, &i, &r, f);
            continue;
        }
        shift_word(&WE.w);
        switch (c) {
            case 'A':
                inEntry++;
                break;
            case 'C':
                WE.w.b[BYTES_IN_WORD - 1] |= 1;
                inEntry++;
                break;
            case 'G':
                WE.w.b[BYTES_IN_WORD - 1] |= 2;
                inEntry++;
                break;
            case 'T':
                WE.w.b[BYTES_IN_WORD - 1] |= 3;
                inEntry++;
                break;
            default :
                inEntry = 0;
                NoACGT++;
                break;
        }
        index++;
        Tot++;
        if (inEntry >= (uint64_t) WORD_SIZE) {
            WE.pos = index - WORD_SIZE;
            memcpy(&words[NW], &WE, sizeof(wentryF));
            NW++;
        }
        c = buffered_fgetc(seq, &i, &r, f);

    }
    //printf("FILE: Create %d Words --(seqLen=%d NoACGT=%d noChar=%d\n",NW,Tot,NoACGT, NoC);
    free(seq);
    fclose(f);

    if (NW == 0)
        terror("Words array empty");

#ifdef VERBOSE
    fprintf(stdout, "Reallocating words of SeqX to: %" PRIu64 "\n", NW);
    fflush(stdout);
#endif

    words = realloc(words, NW * sizeof(wentryF));
    if (words == NULL)
        terror("Error reallocating words of seqX array");

#ifdef VERBOSE
    fprintf(stdout, "Before sorting words\n");
    fflush(stdout);
#endif

    psortWF(32, words, NW);

#ifdef VERBOSE
    fprintf(stdout, "After sorting words\n");
    fflush(stdout);
#endif

#ifdef VERBOSE
    fprintf(stdout, "Before w2hd\n");
#endif
    hashentryF *entries = NULL;
    if (NW == 0)
        terror("Words array empty");
    if ((entries = calloc(REALLOC_FREQ, sizeof(hashentryF))) == NULL) {
        terror("not enough memory for hashentry array");
    }

    memcpy(&entries[0].w.b[0], &words[0].w.b[0], 8);
    entries[0].num = 0;
    entries[0].locs = NULL;
    if ((entries[0].locs = calloc(SIZE_LOC, sizeof(locationF))) == NULL) {
        terror("not enough memory for locs array");
    }
    loc_size = SIZE_LOC;

#ifdef VERBOSE
    fprintf(stdout, "memory allocated w2hd\n");
#endif

    i = 0;
    uint64_t j = 0;
    uint64_t k = 0;
    uint64_t l = REALLOC_FREQ;
    locationF loc;
    while (i < NW) {
        loc.pos = words[k].pos;
        loc.seq = words[k].seq;
        if (wordcmp(&entries[j].w.b[0], &words[k].w.b[0], 32) != 0) {
            if (nLocs == 0)
                terror("nLocs is 0");
            entries[j].locs = realloc(entries[j].locs, nLocs * sizeof(locationF));
            if (entries[j].locs == NULL)
                terror("Error reallocating location array");
            j++;
            if (j >= l) {
                l += REALLOC_FREQ;
                entries = realloc(entries, l * sizeof(hashentryF));
                if (entries == NULL) {
                    terror("Error reallocating entries of seqX array");
                }
            }
            memcpy(&entries[j].w.b[0], &words[k].w.b[0], 8);
            entries[j].num = 0;
            nLocs = 0;
            if ((entries[j].locs = calloc(SIZE_LOC, sizeof(locationF))) == NULL) {
                terror("not enough memory for locs array");
            }
            loc_size = SIZE_LOC;
        }

        if (nLocs >= loc_size) {
            loc_size += SIZE_LOC;
//            fprintf(stdout, "Reallocating memory from: %lu to: %lu\n",loc_size-SIZE_LOC,loc_size);
            entries[j].locs = realloc(entries[j].locs, loc_size * sizeof(locationF));
            if (entries[j].locs == NULL) {
                terror("Error re-allocation location array");
            }

        }
        memcpy(&entries[j].locs[nLocs++], &loc, sizeof(locationF));
        entries[j].num++;
        i++;
        k++;
        if (i % REALLOC_FREQ == 0) {
            memmove(words, words + REALLOC_FREQ, NW - i);
            more_words = realloc(words, (NW - i) * sizeof(wentryF));
            if (more_words == NULL) {
                free(words);
                terror("Error reallocating words of seqX array");
            } else {
                words = more_words;
            }
            k -= REALLOC_FREQ;
        }
    }
#ifdef VERBOSE
    fprintf(stdout, "After w2hd\n");
#endif

    if (j == 0)
        terror("Number of entries in the hash table is 0");
    free(entries[j].locs);
    entries = realloc(entries, j * sizeof(hashentryF));
    if (entries == NULL)
        terror("Error reallocating entries of seqX array");
    free(words);
    *(args->nEntries) = j;
    *(args->seqLen) = Tot;
    return entries;
}

void shift_word_right(word *w) {
    int i;
    for (i = BYTES_IN_WORD - 1; i > 0; i--) {
        w->b[i] >>= 2;
        w->b[i] |= (w->b[i - 1] << 6);
    }
    w->b[i] >>= 2;
}

void *dictionaryWithReverse(void *a) {
    DictionaryArgs *args = (DictionaryArgs *) a;
    FILE *f;
    char c;
    wentryR *words = NULL;
    wentryR *more_words = NULL;
    char *seq = NULL;
    uint64_t Tot = 0;
    uint64_t i = 0, r = 0;
    uint64_t j = 0;
    uint64_t k = 0;
    uint64_t l = REALLOC_FREQ;

    if ((f = fopen(args->seqFile, "rt")) == NULL) {
        fprintf(stdout, "opening sequence file: %s\n", args->seqFile);
        terror("opening sequence file");
    }

    if ((seq = calloc(READBUF, sizeof(char))) == NULL) {
        terror("not enough memory for read buffersequence");
    }

    //To force read
    i = READBUF + 1;

    c = buffered_fgetc(seq, &i, &r, f);
    while (c != '\n') {
        c = buffered_fgetc(seq, &i, &r, f);
    }

#ifdef VERBOSE
    fprintf(stdout, "Calculating sequence length of SeqY\n");
    fflush(stdout);
#endif

    c = buffered_fgetc(seq, &i, &r, f);
    while (!feof(f)) {
        if (!isupper(toupper(c))) {
            if (c == '>') {
                c = buffered_fgetc(seq, &i, &r, f);
                while (c != '\n')
                    c = buffered_fgetc(seq, &i, &r, f);
                c = buffered_fgetc(seq, &i, &r, f);
                continue;
            }
        }
        Tot++;
        c = buffered_fgetc(seq, &i, &r, f);
    }

#ifdef VERBOSE
    fprintf(stdout, "SeqY Len: %" PRIu64 "\n", Tot);
    fflush(stdout);
#endif

    fseek(f, 0, SEEK_SET);
    //To force read
    i = READBUF + 1;

    c = buffered_fgetc(seq, &i, &r, f);
    while (c != '\n') {
        c = buffered_fgetc(seq, &i, &r, f);
    }

#ifdef VERBOSE
    fprintf(stdout, "Allocating memory for words of seqY (including reverse). Size: %" PRIu64 "\n", 2 * Tot);
    fflush(stdout);
#endif
    if (Tot == 0)
        terror("Empty sequence Y file");

    if ((words = calloc(2 * Tot, sizeof(wentryR))) == NULL) {
        terror("not enough memory for words array");
    }
#ifdef VERBOSE
    fprintf(stdout, "Allocated memory for words of seqY (including reverse).\n");
    fflush(stdout);
#endif


    wentryR WE, WER;
    WE.strand = 'f';
    WER.strand = 'r';
    WE.seq = 0;
    WER.seq = 0;
    uint64_t index = 0;
    uint64_t inEntry = 0;
    uint64_t NW = 0;
    uint64_t NoACGT = 0;
    uint64_t NoC = 0;
    uint64_t nLocs = 0;
    uint64_t loc_size = 0;
    c = buffered_fgetc(seq, &i, &r, f);
    while (!feof(f)) {
        if (!isupper(toupper(c))) {
            if (c == '>') {
                c = buffered_fgetc(seq, &i, &r, f);
                while (c != '\n')
                    c = buffered_fgetc(seq, &i, &r, f);
                WE.seq++;
                WER.seq++;
                inEntry = 0;
                index++;
            }
            NoC++;
            c = buffered_fgetc(seq, &i, &r, f);
            continue;
        }
        shift_word(&WE.w);
        shift_word_right(&WER.w);
        switch (c) {
            case 'A':
                WER.w.b[0] |= 192;
                inEntry++;
                break;
            case 'C':
                WER.w.b[0] |= 128;
                WE.w.b[BYTES_IN_WORD - 1] |= 1;
                inEntry++;
                break;
            case 'G':
                WER.w.b[0] |= 64;
                WE.w.b[BYTES_IN_WORD - 1] |= 2;
                inEntry++;
                break;
            case 'T':
                WE.w.b[BYTES_IN_WORD - 1] |= 3;
                inEntry++;
                break;
            default :
                inEntry = 0;
                NoACGT++;
                break;
        }
        index++;
        if (inEntry >= (uint64_t) WORD_SIZE) {
            WE.pos = index - WORD_SIZE;
            if (WE.pos > Tot)
                terror("position of forward word out of the sequence");
            WER.pos = index - 1;
            if (WE.pos > Tot)
                terror("position of reverse word out of the sequence");
            memcpy(&words[NW++], &WE, sizeof(wentryR));
            memcpy(&words[NW++], &WER, sizeof(wentryR));
        }
        c = buffered_fgetc(seq, &i, &r, f);
    }

    fclose(f);
    free(seq);

    if (NW == 0)
        terror("Words array empty");

#ifdef VERBOSE
    fprintf(stdout, "Reallocating words of SeqY to: %" PRIu64 "\n", NW);
    fflush(stdout);
#endif

    words = realloc(words, NW * sizeof(wentryR));
    if (words == NULL)
        terror("Error reallocating words of seqY array");

#ifdef VERBOSE
    fprintf(stdout, "Before sorting\n");
    fflush(stdout);
#endif

    psortWR(32, words, NW);

#ifdef VERBOSE
    fprintf(stdout, "After sorting\n");
    fflush(stdout);
#endif

#ifdef VERBOSE
    fprintf(stdout, "Before w2hd\n");
#endif
    hashentryR *entries = NULL;
    if (NW == 0)
        terror("Words array empty");
    if ((entries = calloc(REALLOC_FREQ, sizeof(hashentryR))) == NULL) {
        terror("not enough memory for hashentry array");
    }

    memcpy(&entries[0].w.b[0], &words[0].w.b[0], 8);
    entries[0].num = 0;
    entries[0].locs = NULL;
    if ((entries[0].locs = calloc(SIZE_LOC, sizeof(locationR))) == NULL) {
        terror("not enough memory for locs array");
    }
    loc_size = SIZE_LOC;

#ifdef VERBOSE
    fprintf(stdout, "memory allocated w2hd\n");
#endif

    locationR loc;
    i = 0;
    k = 0;
    while (i < NW) {
        loc.pos = words[k].pos;
        loc.seq = words[k].seq;
        loc.strand = words[k].strand;
        if (wordcmp(&entries[j].w.b[0], &words[k].w.b[0], 32) != 0) {
            if (nLocs == 0)
                terror("nLocs is 0");
            entries[j].locs = realloc(entries[j].locs, nLocs * sizeof(locationR));
            if (entries[j].locs == NULL)
                terror("Error reallocating location array");
            j++;
            if (j >= l) {
                l += REALLOC_FREQ;
                entries = realloc(entries, l * sizeof(hashentryR));
                if (entries == NULL) {
                    terror("Error reallocating entries of seqX array");
                }
            }
            memcpy(&entries[j].w.b[0], &words[k].w.b[0], 8);
            entries[j].num = 0;
            nLocs = 0;
            if ((entries[j].locs = calloc(SIZE_LOC, sizeof(locationR))) == NULL) {
                terror("not enough memory for locs array");
            }
            loc_size = SIZE_LOC;
        }

        if (nLocs >= loc_size) {
            loc_size += SIZE_LOC;
//            fprintf(stdout, "Reallocating memory from: %lu to: %lu\n",loc_size-SIZE_LOC,loc_size);
            entries[j].locs = realloc(entries[j].locs, loc_size * sizeof(locationR));
            if (entries[j].locs == NULL) {
                terror("Error re-allocation location array");
            }

        }
        memcpy(&entries[j].locs[nLocs++], &loc, sizeof(locationR));
        entries[j].num++;
        i++;
        k++;
        if (i % REALLOC_FREQ == 0) {
            memmove(words, words + REALLOC_FREQ, NW - i);
            more_words = realloc(words, (NW - i) * sizeof(wentryR));
            if (more_words == NULL) {
                free(words);
                terror("Error reallocating words of seqY array");
            } else {
                words = more_words;
            }
            k -= REALLOC_FREQ;
        }
    }
#ifdef VERBOSE
    fprintf(stdout, "After w2hd\n");
#endif

    if (j == 0)
        terror("Number of entries in the hash table is 0");
    free(entries[j].locs);
    entries = realloc(entries, j * sizeof(hashentryR));
    if (entries == NULL)
        terror("Error reallocating entries of seqY array");
    free(words);
    *(args->nEntries) = j;
    *(args->seqLen) = Tot;
    return entries;
}