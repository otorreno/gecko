#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "quicksortWordForward.h"
#include "quicksortWordReverse.h"
#include "karlin.h"

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
    uint64_t Tot = 0;

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

    #ifdef VERBOSE
        fprintf(stdout, "Calculating sequence length of SeqX\n");
        fflush(stdout);
    #endif

    args->nSeqs = 0;
    c = buffered_fgetc(seq, &i, &r, f);
    while ((!feof(f) || (feof(f) && i < r))) {
        if (!isupper(toupper(c))) {
            if (c == '>') {
                args->nSeqs++;
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

    fseek(f, 0, SEEK_SET);

    //To force the read
    i = READBUF + 1;

    c = buffered_fgetc(seq, &i, &r, f);
    while (c != '\n') {
        c = buffered_fgetc(seq, &i, &r, f);
    }

    #ifdef VERBOSE
        fprintf(stdout, "SeqX Len: %" PRIu64 " composed of %" PRIu64" sequence(s)\n", Tot, args->nSeqs);
        fflush(stdout);
    #endif

    args->seqStats = (struct statsHSP *) malloc(args->nSeqs * sizeof(struct statsHSP));
    if(args->seqStats == NULL) terror("Could not allocate memory for stats array");

    wentryF WE;
    WE.seq = 0;
    uint64_t index = 0;
    uint64_t inEntry = 0;
    uint64_t NW = 0;
    Tot = 0;
    uint64_t NoACGT = 0;
    uint64_t NoC = 0;
    uint64_t nLocs = 0;
    uint64_t loc_size = 0;
    args->seqStats[0].tf.total = 0;
    c = buffered_fgetc(seq, &i, &r, f);
    //TODO check the loop conditions
    while (!feof(f) || (feof(f) && i < r)){
        if (!isupper(toupper(c))) {
            if (c == '>') {
                c = buffered_fgetc(seq, &i, &r, f);
                while (c != '\n'){
                    c = buffered_fgetc(seq, &i, &r, f);
                }
                WE.seq++;
                args->seqStats[WE.seq].tf.total = 0;
                inEntry = 0;
                index++;
            }
            NoC++;
            c = buffered_fgetc(seq, &i, &r, f);
            continue;
        }
        args->seqStats[WE.seq].tf.total++;
        shift_word(&WE.w);
        switch (c) {
            case 'A':
                inEntry++;
                args->seqStats[WE.seq].tf.A++;
                break;
            case 'C':
                WE.w.b[BYTES_IN_WORD - 1] |= 1;
                inEntry++;
                args->seqStats[WE.seq].tf.C++;
                break;
            case 'G':
                WE.w.b[BYTES_IN_WORD - 1] |= 2;
                inEntry++;
                args->seqStats[WE.seq].tf.G++;
                break;
            case 'T':
                WE.w.b[BYTES_IN_WORD - 1] |= 3;
                args->seqStats[WE.seq].tf.T++;
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

    //Compute karlin and lambda parameters. Not needed since PAM matrix uniquely determines a pair of Karlin,Lambda
    //computeKarlinLambda(args->seqStats, args->nSeqs);

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
    fprintf(stdout, "Before sorting words of seqX\n");
    fflush(stdout);
#endif

    psortWF(32, words, NW);

#ifdef VERBOSE
    fprintf(stdout, "After sorting words of seqX\n");
    fflush(stdout);
#endif

#ifdef VERBOSE
    fprintf(stdout, "Before w2hd of seqX\n");
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
    fprintf(stdout, "memory allocated w2hd for seqX\n");
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
            memmove(words, words + REALLOC_FREQ, sizeof(wentryF)*(NW - i));
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
    fprintf(stdout, "After w2hd of seqX\n");
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
    /*
    c = buffered_fgetc(seq, &i, &r, f);
    while (c != '\n') {
        c = buffered_fgetc(seq, &i, &r, f);
    }
    */

#ifdef VERBOSE
    fprintf(stdout, "Calculating sequence length of SeqY\n");
    fflush(stdout);
#endif

    args->nSeqs = 0;
    c = buffered_fgetc(seq, &i, &r, f);
    while ((!feof(f) || (feof(f) && i < r))) {
        if (!isupper(toupper(c))) {
            if (c == '>') {
                args->nSeqs++;
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
    fprintf(stdout, "SeqY Len: %" PRIu64 " composed of %" PRIu64" sequence(s)\n", Tot, args->nSeqs);
    fflush(stdout);
#endif

    args->seqStats = (struct statsHSP *) malloc(args->nSeqs * sizeof(struct statsHSP));
    if(args->seqStats == NULL) terror("Could not allocate memory for stats array");

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
    args->seqStats[0].tf.total = 0;
    uint64_t index = 0;
    uint64_t inEntry = 0;
    uint64_t NW = 0;
    uint64_t NoACGT = 0;
    uint64_t NoC = 0;
    uint64_t nLocs = 0;
    uint64_t loc_size = 0;
    c = buffered_fgetc(seq, &i, &r, f);
    while (!feof(f) || (feof(f) && i < r)){
	if (!isupper(toupper(c))) {
            if (c == '>') {
                c = buffered_fgetc(seq, &i, &r, f);
                while (c != '\n'){
                    c = buffered_fgetc(seq, &i, &r, f);
                }
                WE.seq++;
                WER.seq++;
                inEntry = 0;
                index++;
            }
            NoC++;
            c = buffered_fgetc(seq, &i, &r, f);
            continue;
        }
        args->seqStats[WE.seq].tf.total++;
        shift_word(&WE.w);
        shift_word_right(&WER.w);
        switch (c) {
            case 'A':
                WER.w.b[0] |= 192;
                inEntry++;
                args->seqStats[WE.seq].tf.A++;
                break;
            case 'C':
                WER.w.b[0] |= 128;
                WE.w.b[BYTES_IN_WORD - 1] |= 1;
                inEntry++;
                args->seqStats[WE.seq].tf.C++;
                break;
            case 'G':
                WER.w.b[0] |= 64;
                WE.w.b[BYTES_IN_WORD - 1] |= 2;
                inEntry++;
                args->seqStats[WE.seq].tf.G++;
                break;
            case 'T':
                WE.w.b[BYTES_IN_WORD - 1] |= 3;
                inEntry++;
                args->seqStats[WE.seq].tf.T++;
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

    //Not needed since the karlin and lambda parameters do not change if PAM matrix is the same
    //computeKarlinLambda(args->seqStats, args->nSeqs);

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
    fprintf(stdout, "Before sorting seqY\n");
    fflush(stdout);
#endif

    psortWR(32, words, NW);

#ifdef VERBOSE
    fprintf(stdout, "After sorting seqY\n");
    fflush(stdout);
#endif

#ifdef VERBOSE
    fprintf(stdout, "Before w2hd seqY\n");
#endif
    hashentryR *entries = NULL;
    if (NW == 0)
        terror("Words array empty");
    if ((entries = calloc(REALLOC_FREQ, sizeof(hashentryR))) == NULL) {
        terror("not enough memory for hashentry array");
    }

    //Copy first entry
    memcpy(&entries[0].w.b[0], &words[0].w.b[0], 8);
    entries[0].num = 0;
    entries[0].locs = NULL;
    if ((entries[0].locs = calloc(SIZE_LOC, sizeof(locationR))) == NULL) {
        terror("not enough memory for locs array");
    }
    loc_size = SIZE_LOC;

#ifdef VERBOSE
    fprintf(stdout, "memory allocated w2hd seqY\n");
#endif

    locationR loc;
    i = 0;
    k = 0;
    //While there are still words in the dictionary
    while (i < NW) {
	//Copy current location to aux var loc
        loc.pos = words[k].pos;
        loc.seq = words[k].seq;
        loc.strand = words[k].strand;
	//If the last added entry is not equal to the current word in the dictionary (new entry)
        if (wordcmp(&entries[j].w.b[0], &words[k].w.b[0], 32) != 0) {
            if (nLocs == 0)
                terror("nLocs is 0");
	    //Realloc the locations of the last entry to not consume that much memory
            entries[j].locs = realloc(entries[j].locs, nLocs * sizeof(locationR));
            if (entries[j].locs == NULL)
                terror("Error reallocating location array");
	    //Since its a new entry, increase j to work with this new one 
            j++;
	    //If we have no more free allocated entries, realloc everything to get more
            if (j >= l) {
                l += REALLOC_FREQ;
                entries = realloc(entries, l * sizeof(hashentryR));
                if (entries == NULL) {
                    terror("Error reallocating entries of seqX array");
                }
            }
	    //Copy the actual word to the new entry
            memcpy(&entries[j].w.b[0], &words[k].w.b[0], 8);
            entries[j].num = 0;
            nLocs = 0;
	    //Allocate an initial SIZE_LOC of locations for the new entry
            if ((entries[j].locs = calloc(SIZE_LOC, sizeof(locationR))) == NULL) {
                terror("not enough memory for locs array");
            }
            loc_size = SIZE_LOC;
        }
	//If it was not a new entry, but it has too many repetitions, reallocate the locations with an extra SIZE_LOC
        if (nLocs >= loc_size) {
            loc_size += SIZE_LOC;
//            fprintf(stdout, "Reallocating memory from: %lu to: %lu\n",loc_size-SIZE_LOC,loc_size);
            entries[j].locs = realloc(entries[j].locs, loc_size * sizeof(locationR));
            if (entries[j].locs == NULL) {
                terror("Error re-allocation location array");
            }

        }
	//Copy the new location (repetition) to the current entry
        memcpy(&entries[j].locs[nLocs++], &loc, sizeof(locationR));
        entries[j].num++;
        i++;
        k++;
	//This part frees words while on runtime [Currently there is a bug]
	
        if (i % REALLOC_FREQ == 0) {
            memmove(words, words + REALLOC_FREQ, sizeof(wentryR)*(NW - i));
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
    fprintf(stdout, "After w2hd seqY\n");
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

void computeKarlinLambda(struct statsHSP * seqStats, uint64_t nSeqs){
       

    int max_value = POINT;
    int min_value = -POINT;
    //Range from -POINT to POINT
    double * prob_array_range = (double *) malloc((2*POINT + 1) * sizeof(double));

    uint64_t i, j;

    double prob_match = 0, prob_miss = 0;
    double pA, pC, pG, pT;
    double H; //Used in Karlin function
    //double total;
    //Compute K and lambda for each sequence
    /*
    for(i=0;i<nSeqs;i++){
        total = seqStats[i].tf.A + seqStats[i].tf.C + seqStats[i].tf.G + seqStats[i].tf.T;
        pA = ((double) seqStats[i].tf.A/total);
        pC = ((double) seqStats[i].tf.C/total);
        pG = ((double) seqStats[i].tf.G/total);
        pT = ((double) seqStats[i].tf.T/total);
        prob_match =  pA * pA + pC * pC + pG * pG + pT * pT;
        prob_miss = pA * (pC + pG + pT) + pC * (pA + pG + pT) + pG * (pA + pC + pT) + pT * (pA + pC + pG);
        prob_array_range[0] = prob_miss;
        prob_array_range[2*POINT] = prob_match;
        for(j=1;j<2*POINT;j++){
            prob_array_range[j] = 0.0;
        }
        karlin(min_value, max_value, prob_array_range, &seqStats[i].lambda, &seqStats[i].karlin, &H);

        fprintf(stdout, "[INFO] K: %e L:%e\n", seqStats[i].lambda, seqStats[i].karlin);
        getchar();
    }
    */
    
    uint64_t tA=0, tC=0, tG=0, tT=0, all=0;
    for(i=0;i<nSeqs;i++){
        tA += seqStats[i].tf.A;
        tC += seqStats[i].tf.C;
        tG += seqStats[i].tf.G;
        tT += seqStats[i].tf.T;
        all = tA + tC + tG + tT;
    }
    pA = ((double) tA/all);
    pC = ((double) tC/all);
    pG = ((double) tG/all);
    pT = ((double) tT/all);

    prob_match =  pA * pA + pC * pC + pG * pG + pT * pT;
    prob_miss = pA * (pC + pG + pT) + pC * (pA + pG + pT) + pG * (pA + pC + pT) + pT * (pA + pC + pG);
    prob_array_range[0] = prob_miss;
    prob_array_range[2*POINT] = prob_match;
    for(j=1;j<2*POINT;j++){
        prob_array_range[j] = 0.0;
    }
    karlin(min_value, max_value, prob_array_range, &seqStats[0].lambda, &seqStats[0].karlin, &H);
    fprintf(stdout, "[INFO] K: %e L:%e\n", seqStats[0].lambda, seqStats[0].karlin);
    
    free(prob_array_range);
}


boolean karlin(int low,int high,double *pr,double *lambda,double *K,double *H)
{
    int i,j,range,lo,hi,first,last;
    double up,new,sum,Sum,av,beta,oldsum,ratio,ftemp;
    double *p,*P,*ptrP,*ptr1,*ptr2;

    /* Check that scores and their associated probabilities are valid     */

    if (low>=0) {
       fprintf(stderr,"Lowest score must be negative.\n");
       return FALSE;
    }
    for (i=range=high-low;i> -low && !pr[i];--i);
    if (i<= -low) {
       fprintf(stderr,"A positive score must be possible.\n");
       return FALSE;
    }
    for (sum=i=0;i<=range;sum+=pr[i++]) if (pr[i]<0) {
       fprintf(stderr,"Negative probabilities not allowed.\n");
       return FALSE;
    }
    if (sum<0.99995 || sum>1.00005)
       fprintf(stderr,"Probabilities sum to %.4f.  Normalizing.\n",sum);
    //NEW(p,range+1,double);
    p=(double *) calloc(range+1,sizeof(double));
    if(p==NULL) fprintf(stderr,"Out of Memory.");

    for (Sum=low,i=0;i<=range;++i) Sum+=i*(p[i]=pr[i]/sum);
    if (Sum>=0) {
       fprintf(stderr,"Invalid (non-negative) expected score:  %.3f\n",Sum);
       return FALSE;
    }

    /* Calculate the parameter lambda */

    up=0.5;
    do {
        up*=2;
        ptr1=p;
        beta=exp(up);
        ftemp=exp(up*(low-1));
        for (sum=i=0;i<=range;++i) sum+= *ptr1++ * (ftemp*=beta);
    }
    while (sum<1.0);
    for (*lambda=j=0;j<25;++j) {
        new=(*lambda+up)/2.0;
        beta=exp(new);
        ftemp=exp(new*(low-1));
        ptr1=p;
        for (sum=i=0;i<=range;++i) sum+= *ptr1++ * (ftemp*=beta);
        if (sum>1.0) up=new;
        else *lambda=new;
    }

    /* Calculate the pamameter K */

    ptr1=p;
    ftemp=exp(*lambda*(low-1));
    for (av=0,i=low;i<=high;++i) av+= *ptr1++ *i*(ftemp*=beta);
        *H= *lambda*av/log(2.0);
    Sum=lo=hi=0;
    
    //NEW(P,KARLINMAXIT*range+1,double);
    P=(double *) calloc(KARLINMAXIT*range+1,sizeof(double));


    for (*P=sum=oldsum=j=1;j<=KARLINMAXIT && sum>0.001;Sum+=sum/=j++) {
        first=last=range;
        for (ptrP=P+(hi+=high)-(lo+=low);ptrP>=P;*ptrP-- =sum) {
            ptr1=ptrP-first;
            ptr2=p+first;
            for (sum=0,i=first;i<=last;++i) sum+= *ptr1-- * *ptr2++;
            if (first) --first;
            if (ptrP-P<=range) --last;
        }
        ftemp=exp(*lambda*(lo-1));
        for (sum=0,i=lo;i;++i) sum+= *++ptrP * (ftemp*=beta);
        for (;i<=hi;++i) sum+= *++ptrP;
        ratio=sum/oldsum;
        oldsum=sum;
    }
    for (;j<=200;Sum+=oldsum/j++) oldsum*=ratio;
    for (i=low;!p[i-low];++i);
    for (j= -i;i<high && j>1;) if (p[++i-low]) j=karlin_gcd(j,i);
    *K = (j*exp(-2*Sum))/(av*(1.0-exp(- *lambda*j)));
    free(p);
    free(P);
    return TRUE;        /* Parameters calculated successfully */
}

double ExpectedInformation(a_type A, double lambda, double *freq)
{
        int i,j;
        double sum,tot,fij,eij,mij;

        sum = tot = 0.0;

        for (i=1; i<=20; i++){
           for (j=1; j<=20; j++) {
                mij = valAlphaR(i,j,A);
                fij = freq[i]*freq[j];
                tot += fij;
                eij = mij*fij*exp(lambda*mij);
                sum += eij;
           }
        }
        return(lambda*sum/tot);
}

int karlin_gcd(int a,int b)
{
    int c;

    if (b<0) b= -b;
    if (b>a) { c=a; a=b; b=c; }
    for (;b;b=c) { c=a%b; a=b; }
    return a;
}