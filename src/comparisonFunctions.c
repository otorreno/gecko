#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "quicksortHit.h"
#include "comparisonFunctions.h"

#define MAXBUF 1000000
#define MaxREP 100

int FragFromHit(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
                uint64_t n0, struct Sequence *sY,
                uint64_t n1, uint64_t nSeqs1, uint64_t LenOrPer, uint64_t SimTh, int WL);

int FragFromHitReverse(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
                       uint64_t n0, struct Sequence *sY,
                       uint64_t n1, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL);

struct FragFile *frags(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL,
                       uint64_t *nF, uint64_t *nHU, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY, long double e_value);

struct FragFile *fragsReverse(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL,
                              uint64_t *nF, uint64_t *nHU, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY, long double e_value);

long double computeExpectedValueOfFrag(struct FragFile * frag, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY);

inline char complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return c;
    }
}

uint64_t scoreMax(char *seq, char *seq2, uint64_t len, int point) {
    //CHANGE WHEN USING PAM MATRIX
    //Scor=0; for (i=0;ii<len;i++) if (Seq[i]==Seq2[i])scor+=Ptos
    /* TO DELETE WARNINGS*/
    if (seq + 1) {
    }

    if (seq2 + 1) {
    }
    /**********************/

    return len * point;
}

void endianessConversion(char *source, char *target, int numberOfBytes) {
    int i, j;
    for (i = numberOfBytes - 1; i >= 0; i--) {
        j = numberOfBytes - 1 - i;
        target[j] = source[i];
    }
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];

    if (htons(1) == 1) {
        //big endian
        if (fread(&frag->diag, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return;
            terror("Error reading the HSP diagonal");
        }
        if (fread(&frag->xStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        if (fread(&frag->yStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        if (fread(&frag->xEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        if (fread(&frag->yEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        if (fread(&frag->length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        if (fread(&frag->ident, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identities");
        }
        if (fread(&frag->score, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        if (fread(&frag->similarity, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP similarity");
        }
        if (fread(&frag->seqX, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        if (fread(&frag->seqY, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        if (fread(&frag->block, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        frag->strand = fgetc(f);
        if (fread(&frag->evalue, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return;
            terror("Error reading the HSP diagonal");
        }
        endianessConversion(tmpArray, (char *) (&frag->diag), sizeof(int64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->xStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->yStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->xEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->yEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        endianessConversion(tmpArray, (char *) (&frag->length), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identity");
        }
        endianessConversion(tmpArray, (char *) (&frag->ident), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        endianessConversion(tmpArray, (char *) (&frag->score), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP float");
        }
        endianessConversion(tmpArray, (char *) (&frag->similarity), sizeof(float));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqX), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqY), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        endianessConversion(tmpArray, (char *) (&frag->block), sizeof(int64_t));
        frag->strand = fgetc(f);
        if (fread(tmpArray, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
        endianessConversion(tmpArray, (char *) (&frag->evalue), sizeof(long double));
    }
}

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];
    if (htons(1) == 1) {
        //Big endian
        fwrite(&frag->diag, sizeof(int64_t), 1, f);
        fwrite(&frag->xStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->yStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->xEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->yEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->length, sizeof(uint64_t), 1, f);
        fwrite(&frag->ident, sizeof(uint64_t), 1, f);
        fwrite(&frag->score, sizeof(uint64_t), 1, f);
        fwrite(&frag->similarity, sizeof(float), 1, f);
        fwrite(&frag->seqX, sizeof(uint64_t), 1, f);
        fwrite(&frag->seqY, sizeof(uint64_t), 1, f);
        fwrite(&frag->block, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        fwrite(&frag->evalue, sizeof(long double), 1, f);
    } else {
        //Little endian
        endianessConversion((char *) (&frag->diag), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        endianessConversion((char *) (&frag->xStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->xEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->length), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->ident), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->score), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->similarity), tmpArray, sizeof(float));
        fwrite(tmpArray, sizeof(float), 1, f);
        endianessConversion((char *) (&frag->seqX), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->seqY), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->block), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        endianessConversion((char *) (&frag->evalue), tmpArray, sizeof(long double));
        fwrite(tmpArray, sizeof(long double), 1, f);
    }
}

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f) {
    char tmpArray[8];
    if (htons(1) == 1) {
        //big endian
        if (fread(length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading sequence length");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading sequence length");
        }
        endianessConversion(tmpArray, (char *) length, sizeof(uint64_t));
    }
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f) {
    char tmpArray[8];
    if (htons(1) == 1) {
        //big endian
        fwrite(length, sizeof(uint64_t), 1, f);
    } else {
        //little endian
        endianessConversion((char *) length, tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
    }
}

/**
 * Function to return the sizeof a fragment.
 * Due to architeture issues the value of sizeof built-in
 * function could be different
 */
long int sizeofFragment() {
    return 9 * sizeof(uint64_t) + 2 * sizeof(int64_t) + 1 * sizeof(float) + 1 * sizeof(char) + 1*sizeof(long double);
}

int differentSequences(hit h1, hit h2) {
    return h1.seqX != h2.seqX || h1.seqY != h2.seqY;
}

uint64_t filterHits(hit *hBuf, uint64_t hitsInBuf, int wSize) {
    if (hitsInBuf == 0 || hBuf == NULL) {
        return 0;
    }

    int64_t diagonal;
    uint64_t lastPosition;
    uint64_t originalNumberOfHits = 0, finalNumberOfHits = 0;

    lastPosition = 0;

    originalNumberOfHits = hitsInBuf;
    uint64_t i = 0;
    hit lastHitNotFiltered;

    memcpy(&lastHitNotFiltered, &hBuf[0], sizeof(hit));
    finalNumberOfHits++;
    diagonal = lastHitNotFiltered.posX - lastHitNotFiltered.posY; 
    while (i < hitsInBuf-1) {
        if (differentSequences(lastHitNotFiltered, hBuf[i + 1])) {
            lastPosition = hBuf[i+1].posX + (2*wSize-1);
            diagonal = hBuf[i+1].posX - hBuf[i+1].posY;
	        memcpy(&lastHitNotFiltered, &hBuf[i+1], sizeof(hit));
	        memmove(&hBuf[finalNumberOfHits++], &hBuf[i+1], sizeof(hit));
            i++;
            continue;
        }

        if (diagonal != (hBuf[i + 1].posX - hBuf[i + 1].posY) || hBuf[i + 1].posX > lastPosition) {
            lastPosition = hBuf[i+1].posX + (2 * wSize - 1);
            diagonal = hBuf[i+1].posX - hBuf[i+1].posY;
	        memcpy(&lastHitNotFiltered, &hBuf[i+1], sizeof(hit));
            memmove(&hBuf[finalNumberOfHits++], &hBuf[i+1], sizeof(hit));
        }
        i++;
    }

    if (diagonal != (hBuf[i].posX - hBuf[i].posY) || hBuf[i].posX > (lastPosition)) {
        memmove(&hBuf[finalNumberOfHits++], &hBuf[i], sizeof(hit));
    }

    if (finalNumberOfHits == 0) {
        perror("0 Hits. At least one hit should be kept in the filtering");
    }
    hBuf = realloc(hBuf, finalNumberOfHits * sizeof(hit));
    if (hBuf == NULL)
        perror("Error reallocating filtered hits array");

    fprintf(stdout,
            "filterHits original number of Hits=%"
    PRIu64
    "  Final number of hits=%"
    PRIu64
    "\n",
            originalNumberOfHits, finalNumberOfHits);

    return finalNumberOfHits;
    //End of filterhits
}

uint64_t filterHitsReverse(hit *hBuf, uint64_t hitsInBuf, int wSize, uint64_t minSeqXLenSeqYLen) {
    if (hitsInBuf == 0 || hBuf == NULL) {
        return 0;
    }

    int64_t diagonal;
    uint64_t lastPosition;
    uint64_t originalNumberOfHits = 0, finalNumberOfHits = 0;

    lastPosition = 0;
    originalNumberOfHits = hitsInBuf;
    uint64_t i = 0;

    hit lastHitNotFiltered;

    memcpy(&lastHitNotFiltered, &hBuf[0], sizeof(hit));
    finalNumberOfHits++;
    diagonal = lastHitNotFiltered.posX - lastHitNotFiltered.posY; 
    while (i < hitsInBuf-1) {
        if (differentSequences(lastHitNotFiltered, hBuf[i + 1])) {
            lastPosition = hBuf[i+1].posX + (2*wSize-1);
            diagonal = hBuf[i+1].posX + hBuf[i+1].posY - minSeqXLenSeqYLen;
	        memcpy(&lastHitNotFiltered, &hBuf[i+1], sizeof(hit));
            memmove(&hBuf[finalNumberOfHits++], &hBuf[i+1], sizeof(hit));
            i++;
            continue;
        }

        if (diagonal != (hBuf[i + 1].posX + hBuf[i + 1].posY - minSeqXLenSeqYLen) || hBuf[i + 1].posX > lastPosition) {
            lastPosition = hBuf[i+1].posX + (2 * wSize - 1);
            diagonal = hBuf[i+1].posX + hBuf[i+1].posY - minSeqXLenSeqYLen;
	        memcpy(&lastHitNotFiltered, &hBuf[i+1], sizeof(hit));
            memmove(&hBuf[finalNumberOfHits++], &hBuf[i+1], sizeof(hit));
        }
        i++;
    }

    if (diagonal != (hBuf[i].posX + hBuf[i].posY - minSeqXLenSeqYLen) || hBuf[i].posX > (lastPosition)) {
        memmove(&hBuf[finalNumberOfHits++], &hBuf[i], sizeof(hit));
    }

    if (finalNumberOfHits == 0) {
        perror("0 Hits. At least one hit should be kept in the filtering");
    }
    hBuf = realloc(hBuf, finalNumberOfHits * sizeof(hit));
    if (hBuf == NULL)
        perror("Error reallocating filtered hits array");

    fprintf(stdout,
            "filterHits original number of Hits=%"
    PRIu64
    "  Final number of hits=%"
    PRIu64
    "\n",
            originalNumberOfHits, finalNumberOfHits);

    return finalNumberOfHits;
    //End of filterhits
}

void *sortHitsFilterHitsFragHitsTh(void *a) {
    ComparisonArgs *args = (ComparisonArgs *) a;
    uint64_t HIB;

#ifdef ELAPSEDTIME
    //time variables
    clock_t begin, end;
    double elapsed_secs;

    begin = clock();
#endif

    if (args->nHits > 0) {
#ifdef VERBOSE
        fprintf(stdout, "Sorting Hits 1\n");
        fflush(stdout);
#endif
        psortHF(32, args->hits, args->nHits);

#ifdef ELAPSEDTIME
        end = clock();
        elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
        fprintf(stdout, "[SORTHITS] Elapsed time: %lf\n", elapsed_secs);
        begin = clock();
#endif

#ifdef VERBOSE
        fprintf(stdout, "End of Sorting Hits 1\n");
        fflush(stdout);
#endif

#ifdef VERBOSE
        fprintf(stdout, "Filtering Hits 1\n");
        fflush(stdout);
#endif

        HIB = filterHits(args->hits, args->nHits, args->wSize);

#ifdef ELAPSEDTIME
        end = clock();
        elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
        fprintf(stdout, "[FILTERHITS] Elapsed time: %lf\n", elapsed_secs);
        begin = clock();
#endif

#ifdef VERBOSE
        fprintf(stdout, "End of filtering Hits 1\n");
        fflush(stdout);
#endif


#ifdef VERBOSE
        fprintf(stdout, "Frag Hits 1\n");
        fflush(stdout);
#endif

        struct FragFile *fragsBuf = NULL;
        if (HIB > 0) {
            fragsBuf = frags(args->seqX, args->seqY, args->hits, HIB, args->Lmin, args->SimTh,
                             args->wSize, args->nFrags, args->nHitsUsed, args->seqStatsX, args->seqStatsY, args->e_value);
#ifdef ELAPSEDTIME
            end = clock();
            elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
            fprintf(stdout, "[FRAGHITSREV] Elapsed time: %lf\n", elapsed_secs);
#endif
        }

#ifdef VERBOSE
        fprintf(stdout, "End of Frag Hits 1\n");
        fflush(stdout);
#endif

        if (args->hits != NULL)
            free(args->hits);
        return fragsBuf;
    }

    return NULL;
}

void *sortHitsFilterHitsFragHitsReverseTh(void *a) {
    ComparisonArgs *args = (ComparisonArgs *) a;
    uint64_t HIB;

#ifdef ELAPSEDTIME
    //time variables
    clock_t begin, end;
    double elapsed_secs;

    begin = clock();
#endif

    if (args->nHits > 0) {
        psortHR(32, args->hits, args->nHits, args->minSeqLen);
#ifdef ELAPSEDTIME
        end = clock();
        elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
        fprintf(stdout, "[SORTHITSREV] Elapsed time: %lf\n", elapsed_secs);
        begin = clock();
#endif

        HIB = filterHitsReverse(args->hits, args->nHits, args->wSize, args->minSeqLen);

#ifdef ELAPSEDTIME
        end = clock();
        elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
        fprintf(stdout, "[FILTERHITSREV] Elapsed time: %lf\n", elapsed_secs);
        begin = clock();
#endif

#ifdef VERBOSE
        fprintf(stdout, "Frag Hits 2\n");
        fflush(stdout);
#endif


        struct FragFile *fragsBuf = NULL;
        if (HIB > 0) {
            fragsBuf = fragsReverse(args->seqX, args->seqY, args->hits, HIB, args->Lmin, args->SimTh,
                                    args->wSize, args->nFrags, args->nHitsUsed, args->seqStatsX, args->seqStatsY, args->e_value);
#ifdef ELAPSEDTIME
            end = clock();
            elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
            fprintf(stdout, "[FRAGHITSREV] Elapsed time: %lf\n", elapsed_secs);
#endif
        }
#ifdef VERBOSE
        fprintf(stdout, "End of Frag Hits 2\n");
        fflush(stdout);
#endif


        if (args->hits != NULL)
            free(args->hits);

        return fragsBuf;
    }
    return NULL;
}

struct FragFile *hitsAndFrags(char *seqX, char *seqY, char *out, uint64_t seqXLen, uint64_t seqYLen,
                              hashentryF *entriesX, uint64_t nEntriesX, hashentryR *entriesY, uint64_t nEntriesY,
                              int wSize, uint64_t Lmin, uint64_t SimTh, uint64_t *nFrags, struct statsHSP * seqStatsX,
                              struct statsHSP * seqStatsY, long double e_value) {
    int bufSizeForward = MAXBUF, bufSizeReverse = MAXBUF;
    uint64_t hitsInBufForward = 0, hitsInBufReverse = 0, hitsInBufForwardUsed = 0, hitsInBufReverseUsed = 0,
            nFragsForward = 0, nFragsReverse = 0;
    uint64_t i, j, k = 0, l = 0;
    int comp;
    int firstMatch = 0, endMatch = 0;
    uint64_t nHits = 0, wordMatches = 0;
    int stepX, stepY;
    long offsetWordBefore = 0, offsetWordAfter = 0;
    hit *hBufForward, *hBufReverse;
    uint64_t minSeqXLenSeqYLen = min(seqXLen, seqYLen);
    FILE *fOut, *fInf;
    char infoFileName[1024];

#ifdef ELAPSEDTIME
    //time variables
    clock_t begin, end;
    double elapsed_secs;

    begin = clock();
#endif

    if ((hBufForward = (hit *) calloc(MAXBUF, sizeof(hit))) == NULL)
        terror("HITS: memory for I-O buffer");

    if ((hBufReverse = (hit *) calloc(MAXBUF, sizeof(hit))) == NULL)
        terror("HITS: memory for I-O buffer");

#ifdef VERBOSE
    fprintf(stdout, "Memory allocated for the hits buffer\n");
    fflush(stdout);
#endif

    while (k < nEntriesX && l < nEntriesY) {
        comp = wordcmp(&entriesX[k].w.b[0], &entriesY[l].w.b[0], wSize);
//		fprintf(stdout, "k: %" PRIu64 " l:%" PRIu64 "\n", k, l);
        if (comp < 0) {
            k++;
            //Save position of first missmatch after matches and rewind
            if (firstMatch) {
                offsetWordAfter = l - 1;
                l = offsetWordBefore + 1;
                firstMatch = 0;
                endMatch = 1;
            }
            continue;
        }
        if (comp > 0) {
            //No more matches, go to the next word
            if (endMatch) {
                l = offsetWordAfter;
                endMatch = 0;
            }
            l++;
            continue;
        }

        wordMatches++;
//		fprintf(stdout, "Hay match\n");

        // Saving the offset of the first match
        if (wSize < 32 && !firstMatch) {
            offsetWordBefore = l - 1;
            firstMatch = 1;
        }

        // Hits-----------------------
        if (entriesX[k].num > MaxREP)
            stepX = entriesX[k].num / MaxREP;
        else
            stepX = 1;
        if (entriesY[l].num > MaxREP)
            stepY = entriesY[l].num / MaxREP;
        else
            stepY = 1;

        for (i = 0; i < entriesX[k].num; i += stepX)
            for (j = 0; j < entriesY[l].num; j += stepY) {
                if (entriesY[l].locs[j].strand == 'f') {
                    hBufForward[hitsInBufForward].posX = entriesX[k].locs[i].pos;
                    hBufForward[hitsInBufForward].seqX = entriesX[k].locs[i].seq;
                    hBufForward[hitsInBufForward].posY = entriesY[l].locs[j].pos;
                    hBufForward[hitsInBufForward].seqY = entriesY[l].locs[j].seq;

                    hitsInBufForward++;
                    if (hitsInBufForward >= bufSizeForward) {
                        hBufForward = realloc(hBufForward, (bufSizeForward + MAXBUF) * sizeof(hit));
                        if (hBufForward == NULL)
                            perror("Error reallocating forward hits buffer");
                        bufSizeForward += MAXBUF;
                    }
                } else {
                    hBufReverse[hitsInBufReverse].posX = entriesX[k].locs[i].pos;
                    hBufReverse[hitsInBufReverse].seqX = entriesX[k].locs[i].seq;
                    hBufReverse[hitsInBufReverse].posY = entriesY[l].locs[j].pos;
                    hBufReverse[hitsInBufReverse].seqY = entriesY[l].locs[j].seq;

                    hitsInBufReverse++;
                    if (hitsInBufReverse >= bufSizeReverse) {
//                    fprintf(stdout, "reallocating\n");
                        hBufReverse = realloc(hBufReverse, (bufSizeReverse + MAXBUF) * sizeof(hit));
                        if (hBufReverse == NULL)
                            perror("Error reallocating reverse complement hits buffer");
                        bufSizeReverse += MAXBUF;
                    }
                }
            }

        nHits += ((entriesX[k].num / stepX) * (entriesY[l].num / stepY));

        if (!firstMatch)k++;
        l++;
    }

    if (hitsInBufForward > 0)
        hBufForward = realloc(hBufForward, hitsInBufForward * sizeof(hit));
    if (hBufForward == NULL)
        perror("Error reallocating forward hits buffer");

#ifdef VERBOSE
    fprintf(stdout, "hitsInBufForward: %"
    PRIu64
    "\n", hitsInBufForward);
    fflush(stdout);
#endif

    if (hitsInBufReverse > 0)
        hBufReverse = realloc(hBufReverse, hitsInBufReverse * sizeof(hit));
    if (hBufReverse == NULL)
        perror("Error reallocating reverse complement hits buffer");

#ifdef VERBOSE
    fprintf(stdout, "hitsInBufReverse: %"
    PRIu64
    "\n", hitsInBufReverse);
    fflush(stdout);
#endif

#ifdef ELAPSEDTIME
    end = clock();
    elapsed_secs = (double)(end-begin)/CLOCKS_PER_SEC;
    fprintf(stdout, "[HITS] Elapsed time: %lf\n", elapsed_secs);
#endif



    ComparisonArgs argsForward, argsReverse;
    pthread_t thF, thR;

    argsForward.seqX = seqX;
    argsForward.seqY = seqY;
    argsForward.nHits = hitsInBufForward;
    argsForward.nHitsUsed = &hitsInBufForwardUsed;
    argsForward.hits = hBufForward;
    argsForward.nFrags = &nFragsForward;
    argsForward.Lmin = Lmin;
    argsForward.SimTh = SimTh;
    argsForward.wSize = wSize;
    argsForward.minSeqLen = minSeqXLenSeqYLen;
    argsForward.seqStatsX = seqStatsX;
    argsForward.seqStatsY = seqStatsY;
    argsForward.e_value = e_value;

    argsReverse.seqX = seqX;
    argsReverse.seqY = seqY;
    argsReverse.nHits = hitsInBufReverse;
    argsReverse.nHitsUsed = &hitsInBufReverseUsed;
    argsReverse.hits = hBufReverse;
    argsReverse.nFrags = &nFragsReverse;
    argsReverse.Lmin = Lmin;
    argsReverse.SimTh = SimTh;
    argsReverse.wSize = wSize;
    argsReverse.minSeqLen = minSeqXLenSeqYLen;
    argsReverse.seqStatsX = seqStatsX;
    argsReverse.seqStatsY = seqStatsY;
    argsReverse.e_value = e_value;

    /*
    FILE * fileout = fopen("db_nosort_f.hits", "wb");
    fwrite(hBufForward, sizeof(hit), hitsInBufForward, fileout);
    fclose(fileout);
    fileout = fopen("db_nosort_r.hits", "wb");
    fwrite(hBufReverse, sizeof(hit), hitsInBufReverse, fileout);
    fclose(fileout);
    */

    pthread_create(&thF, NULL, sortHitsFilterHitsFragHitsTh, (void *) (&argsForward));
    pthread_create(&thR, NULL, sortHitsFilterHitsFragHitsReverseTh, (void *) (&argsReverse));

    struct FragFile *fragsBufForward;
    struct FragFile *fragsBufReverse;
    pthread_join(thF, (void **) &fragsBufForward);
    pthread_join(thR, (void **) &fragsBufReverse);

    fprintf(stdout, "Frags in forward: %"
    PRIu64
    "; Frags in reverse: %"
    PRIu64
    "\n", nFragsForward, nFragsReverse);
    fflush(stdout);

    *nFrags = nFragsForward + nFragsReverse;

#ifdef VERBOSE
    fprintf(stdout, "Writing fragments to disk\n");
    fflush(stdout);
#endif

    if ((fOut = fopen(out, "wb")) == NULL)
        terror("opening fragments output file");

#ifdef VERBOSE
    fprintf(stdout, "SeqXLen: %"
    PRIu64
    " SeqYLen: %"
    PRIu64
    "\n", seqXLen, seqYLen);
    fflush(stdout);
#endif

    writeSequenceLength(&seqXLen, fOut);
    writeSequenceLength(&seqYLen, fOut);

#ifdef VERBOSE
    fprintf(stdout, "Writting forward frags\n");
    fflush(stdout);
#endif

    if (fragsBufForward != NULL) {
        for (i = 0; i < nFragsForward; i++) {
            writeFragment(&fragsBufForward[i], fOut);
        }
        free(fragsBufForward);
    }

#ifdef VERBOSE
    fprintf(stdout, "Forward frags written\n");
    fflush(stdout);
#endif

    if (fragsBufReverse != NULL) {
        for (i = 0; i < nFragsReverse; i++) {
            //this was causing errors
            //uint64_t tmp;
            //tmp = seqYLen - fragsBufReverse[i].yEnd;
            //fragsBufReverse[i].yStart = seqYLen - fragsBufReverse[i].yStart;
            //fragsBufReverse[i].yEnd = tmp;
            writeFragment(&fragsBufReverse[i], fOut);
        }
        free(fragsBufReverse);
    }

#ifdef VERBOSE
    fprintf(stdout, "Reverse frags written\n");
    fflush(stdout);
#endif

    fclose(fOut);

    // metadata info
    sprintf(infoFileName, "%s.INF", out);
    if ((fInf = fopen(infoFileName, "wt")) == NULL)
        terror("opening INFO file");

    fprintf(fInf, "GECKO in memory\n");
    fprintf(fInf, "Jun. 2016 -- <oscart@uma.es>\n");
    fprintf(fInf, "SeqX filename        : %s\n", seqX);
    fprintf(fInf, "SeqY filename        : %s\n", seqY);
    fprintf(fInf, "SeqX name            : %s\n", seqX);
    fprintf(fInf, "SeqY name            : %s\n", seqY);
    fprintf(fInf, "SeqX length          : %"
    PRIu64
    "\n", seqXLen);
    fprintf(fInf, "SeqY length          : %"
    PRIu64
    "\n", seqYLen);
    fprintf(fInf, "Min.fragment.length  : %"
    PRIu64
    "\n", Lmin);
    fprintf(fInf, "Min.Identity         : %"
    PRIu64
    "\n", SimTh);
    fprintf(fInf, "Tot Hits (seeds)     : %"
    PRIu64
    "\n", hitsInBufForward + hitsInBufReverse);
    fprintf(fInf, "Tot Hits (seeds) used: %"
    PRIu64
    "\n", hitsInBufForwardUsed + hitsInBufReverseUsed);
    fprintf(fInf, "Total fragments      : %"
    PRIu64
    "\n", *nFrags);
    fprintf(fInf, "========================================================\n");
    fclose(fInf);

#ifdef VERBOSE
    fprintf(stdout, "End of hits and frags function\n");
    fflush(stdout);
#endif

    //TODO change if we want to have the fragments outside this function
    return NULL;
}

struct FragFile *frags(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL,
                       uint64_t *nF, uint64_t *nHU, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY, long double e_value) {
    struct Sequence *sX, *sY;
    uint64_t n0, n1, nSeqs0 = 0, nSeqs1 = 0;
    uint64_t i, j;
    int newFrag;
    uint64_t e_val_filtered = 0;
    uint64_t nFrags = 0, nHitsUsed = 0;
    int64_t lastOffset, lastDiagonal;
    struct FragFile myF;
    FILE *f;

    //MAT
    long coverage[1000][100];
    for (i = 0; i < 1000; i++) {
        for (j = 0; j < 100; j++) {
            coverage[i][j] = 0;
        }
    }
    //---

    if (!nHits) {
        *nHU = 0;
        *nF = 0;
        return NULL;
    }

    struct FragFile *fragsBuf;
    fragsBuf = (struct FragFile *) calloc(nHits, sizeof(struct FragFile));
    if (fragsBuf == NULL)
        terror("Not enough memory for frags array");

    //Open files

    if ((f = fopen(seqY, "rt")) == NULL)
        terror("opening seqY file");
    sY = LeeSeqDB(f, &n1, &nSeqs1, 0);
    fclose(f);

   

    if ((f = fopen(seqX, "rt")) == NULL)
        terror("opening seqX file");
    sX = LeeSeqDB(f, &n0, &nSeqs0, 0);
    fclose(f);
    

    n0 += min(nSeqs0 - 1, 0);
    n1 += min(nSeqs1 - 1, 0);

    lastDiagonal = hits[0].posX - hits[0].posY;
    lastOffset = hits[0].posX - 1;

    i = 0;
    while (i < nHits) {
        if (lastDiagonal != (hits[i].posX - hits[i].posY)) {
            //Different diagonal update the variables
            lastDiagonal = hits[i].posX - hits[i].posY;
            lastOffset = hits[i].posX - 1;
        }
        //We have a casting here because of a funny error
        //where the program was saying that -1 > 0
        if (lastOffset > ((int64_t) hits[i].posX)) {
            //The hit is covered by a previous frag in the same diagonal
            i++;
            continue;
        }

        nHitsUsed++;
        newFrag = FragFromHit(coverage, &myF, &hits[i], sX, n0, sY, n1, nSeqs1, Lmin, SimTh,
                              WL);

        if(newFrag){
            if(computeExpectedValueOfFrag(&myF, seqStatsX, seqStatsY) <= e_value){
                newFrag = 1;
            }else{
                newFrag = 0;
                e_val_filtered++;
            }
        }
	
        if (newFrag) {
            memcpy(&fragsBuf[nFrags], &myF, sizeof(struct FragFile));
            lastOffset = hits[i].posX + myF.length;
            nFrags++;
        }
        i++;
    }

    *nF = nFrags;
    *nHU = nHitsUsed;

#ifdef VERBOSE
    fprintf(stdout, "Forward strand frags calculated\n");
    fprintf(stdout, "%"PRIu64" out of %"PRIu64" frags filtered due to e-value\n", e_val_filtered, nFrags);
    fflush(stdout);
#endif

    if (nFrags > 0)
        fragsBuf = realloc(fragsBuf, nFrags * sizeof(struct FragFile));
    if (fragsBuf == NULL)
        perror("Error reallocating memory of forward frags array");

    if (sX != NULL) {
        free(sX->datos);
        free(sX);
    }

    if (sY != NULL) {
        free(sY->datos);
        free(sY);
    }

    return fragsBuf;
}

struct FragFile *fragsReverse(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL,
                              uint64_t *nF, uint64_t *nHU, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY, long double e_value) {
    struct Sequence *sX, *sY;
    uint64_t n0, n1, nSeqs0 = 0, nSeqs1 = 0;
    uint64_t i, j;
    int newFrag;
    uint64_t e_val_filtered = 0;
    uint64_t nFrags = 0, nHitsUsed = 0;
    int64_t lastOffset, lastDiagonal;
    struct FragFile myF;
    FILE *f;

    //MAT
    long coverage[1000][100];
    for (i = 0; i < 1000; i++) {
        for (j = 0; j < 100; j++) {
            coverage[i][j] = 0;
        }
    }
    //---

    if (!nHits) {
        *nHU = 0;
        *nF = 0;
        return NULL;
    }

    struct FragFile *fragsBuf;
    fragsBuf = (struct FragFile *) calloc(nHits, sizeof(struct FragFile));
    if (fragsBuf == NULL)
        terror("Not enough memory for frags array");

    //Open files
    if ((f = fopen(seqX, "rt")) == NULL)
        terror("opening seqX file");
    sX = LeeSeqDB(f, &n0, &nSeqs0, 0);
    fclose(f);

    if ((f = fopen(seqY, "rt")) == NULL)
        terror("opening seqY file");
    sY = LeeSeqDB(f, &n1, &nSeqs1, 0);
    fclose(f);

    n0 += min(nSeqs0 - 1, 0);
    n1 += min(nSeqs1 - 1, 0);

    // read Hits
    if (nHits == 0) {
        terror("Empty hits file");
    }
    lastDiagonal = hits[0].posX + hits[0].posY - min(n0, n1);
    lastOffset = hits[0].posX - 1;

    i = 0;
    while (i < nHits) {
        if (lastDiagonal != (hits[i].posX + hits[i].posY - min(n0, n1))) {
            //Different diagonal update the variables
            lastDiagonal = hits[i].posX + hits[i].posY - min(n0, n1);
            lastOffset = hits[i].posX - 1;
        }
        //We have a casting here because of a funny error
        //where the program was saying that -1 > 0
        if (lastOffset > ((int64_t) hits[i].posX)) {
            //The hit is covered by a previous frag in the same diagonal
            i++;
            continue;
        }

        nHitsUsed++;
        newFrag = FragFromHitReverse(coverage, &myF, &hits[i], sX, n0, sY, n1, nSeqs1, Lmin, SimTh,
                                     WL);

        if(newFrag){
            if(computeExpectedValueOfFrag(&myF, seqStatsX, seqStatsY) <= e_value){
                newFrag = 1;
            }else{
                newFrag = 0;
                e_val_filtered++;
            }
        }

        if (newFrag) {
            memcpy(&fragsBuf[nFrags], &myF, sizeof(struct FragFile));
            lastOffset = hits[i].posX + myF.length;
            nFrags++;
        }
        i++;
    }

    *nF = nFrags;
    *nHU = nHitsUsed;

#ifdef VERBOSE
    fprintf(stdout, "Reverse strand frags calculated\n");
    fprintf(stdout, "%"PRIu64" out of %"PRIu64" frags filtered due to e-value\n", e_val_filtered, nFrags);
    fflush(stdout);
#endif

    if (nFrags > 0)
        fragsBuf = realloc(fragsBuf, nFrags * sizeof(struct FragFile));
    if (fragsBuf == NULL)
        perror("Error reallocating memory of reverse frags array");

    if (sX != NULL) {
        free(sX->datos);
        free(sX);
    }

    if (sY != NULL) {
        free(sY->datos);
        free(sY);
    }

    return fragsBuf;
}


/**
 * Compute a fragments from one seed point
 * Similarirty thershold and length > mimL
 */
int FragFromHit(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
                uint64_t n0, struct Sequence *sY,
                uint64_t n1, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL) {
    int64_t ldiag, ldiag2;
    int64_t xfil, ycol;
    /* for version with backward search */
    int64_t xfil2, ycol2;
    int fragmentLength = 0;
    /* for version Maximum global---*/
    int64_t xfilmax, ycolmax;
    /* for version with backward search */
    int64_t xfilmax2, ycolmax2;
    int nIdentities = 0, maxIdentities = 0;
    char valueX, valueY;
    int fscore = 0, fscoreMax = 0; // full score

    uint64_t minLength = Lm;

    // Initialize values
    ldiag = min(n0 - H->posX, n1 - H->posY);
    //var to know how much we have until we reach the origin of coordinates
    ldiag2 = min(H->posX, H->posY);
    xfil = H->posX + WL;
    xfil2 = H->posX - 1;
    ycol = H->posY + WL;
    ycol2 = H->posY - 1;
    fragmentLength += WL;
    xfilmax = xfil;
    xfilmax2 = xfil2;
    ycolmax = ycol;
    ycolmax2 = ycol2;
    nIdentities = maxIdentities = WL;
    fscore = POINT * WL;
    fscoreMax = fscore;

    // now, looking for end_frag---
    while (fragmentLength < ldiag) {
        valueX = getValue(sX, xfil);
        valueY = getValue(sY, ycol);
        if (valueX == '*' || valueY == '*') {
            //separator between sequences
            break;
        }

        if (valueX == 'N' || valueY == 'N') {
            fscore -= 1;
        } else {
            if (valueX == valueY) {
                // match
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
        }

        xfil++;
        ycol++;
        fragmentLength++;
        if (fscore < 0)
            break;
    }

    /**
     * Backward search --- Oscar (Sept.2013)
     **/
    fragmentLength = 0;
    fscore = fscoreMax;
    xfilmax2 = H->posX;
    ycolmax2 = H->posY;
    nIdentities = maxIdentities;
    if (xfil2 >= 0 && ycol2 >= 0)
        while (fragmentLength < ldiag2) {
            valueX = getValue(sX, xfil2);
            valueY = getValue(sY, ycol2);
            if (valueX == '*' || valueY == '*') {
                //separator between sequences
                break;
            }

            if (valueX == 'N' || valueY == 'N') {
                fscore -= 1;
            } else {
                if (valueX == valueY) {
                    // matches----
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
            }

            xfil2--;
            ycol2--;
            fragmentLength++;
            if (fscore < 0)
                break;
        }

    // Set the values of the FragFile
    myF->diag = H->posX - H->posY;
    myF->xStart = (uint64_t) xfilmax2 - H->seqX;
    myF->yStart = (uint64_t) ycolmax2 - H->seqY;
    myF->xEnd = (uint64_t) xfilmax - H->seqX;
    myF->yEnd = (uint64_t) ycolmax - H->seqY;;
    myF->length = myF->xEnd - myF->xStart + 1;
    myF->score = fscoreMax;
    myF->ident = maxIdentities;
    myF->similarity = myF->score * 100.0
                      / scoreMax(&sX->datos[myF->xStart], &sY->datos[myF->yStart],
                                 myF->length, POINT);
    myF->seqX = H->seqX;
    myF->seqY = H->seqY;
    myF->block = 0;
    myF->strand = 'f';

    M[min(myF->length, 999)][(int) myF->similarity]++;

    if (myF->length > minLength && myF->similarity > SimTh)
        return 1;
    else
        return 0;
}

/**
 * Compute a fragments from one seed point
 * Similarirty thershold and length > mimL
 */
int FragFromHitReverse(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
                       uint64_t n0, struct Sequence *sY,
                       uint64_t n1, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL) {
    int64_t ldiag, ldiag2;
    int64_t xfil, ycol;
    /* for version with backward search */
    int64_t xfil2, ycol2;
    int fragmentLength = 0;
    /* for version Maximum global---*/
    int64_t xfilmax, ycolmax;
    /* for version with backward search */
    int64_t xfilmax2, ycolmax2;
    int nIdentities = 0, maxIdentities = 0;
    char valueX, valueY;
    int fscore = 0, fscoreMax = 0; // full score

    uint64_t minLength = Lm;

    // Initialize values
    ldiag = min(H->posY, n0 - H->posX);
    //var to know how much we have until we reach the origin of coordinates
    ldiag2 = min(n1 - H->posY, H->posX);
    xfil = H->posX + WL;
    xfil2 = H->posX - 1;
    ycol = H->posY - WL;
    ycol2 = H->posY + 1;
    fragmentLength += WL;
    xfilmax = xfil;
    xfilmax2 = xfil2;
    ycolmax = ycol;
    ycolmax2 = ycol2;
    nIdentities = maxIdentities = WL;
    fscore = POINT * WL;
    fscoreMax = fscore;

    // now, looking for end_frag---
    while (fragmentLength < ldiag) {
        valueX = getValue(sX, xfil);
        valueY = complement(getValue(sY, ycol));
        if (valueX == '*' || valueY == '*') {
            //separator between sequences
            break;
        }

        if (valueX == 'N' || valueY == 'N') {
            fscore -= 1;
        } else {
            if (valueX == valueY) {
                // match
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
        }

        xfil++;
        ycol--;
        fragmentLength++;
        if (fscore < 0)
            break;
    }

    /**
     * Backward search --- Oscar (Sept.2013)
     **/
    fragmentLength = 0;
    fscore = fscoreMax;
    xfilmax2 = H->posX;
    ycolmax2 = H->posY;
    nIdentities = maxIdentities;
    if (xfil2 >= 0 && ycol2 < n1)
        while (fragmentLength < ldiag2) {
            valueX = getValue(sX, xfil2);
            valueY = complement(getValue(sY, ycol2));
            if (valueX == '*' || valueY == '*') {
                //separator between sequences
                break;
            }

            if (valueX == 'N' || valueY == 'N') {
                fscore -= 1;
            } else {
                if (valueX == valueY) {
                    // matches----
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
            }

            xfil2--;
            ycol2++;
            fragmentLength++;
            if (fscore < 0)
                break;
        }

    // Set the values of the FragFile
    myF->diag = H->posX + H->posY - min(n0, n1);
    myF->xStart = (uint64_t) xfilmax2 - H->seqX;
    myF->yStart = (uint64_t) ycolmax2 - H->seqY;
    myF->xEnd = (uint64_t) xfilmax - H->seqX;
    myF->yEnd = (uint64_t) ycolmax - H->seqY;;
    myF->length = myF->xEnd - myF->xStart + 1;
    myF->score = fscoreMax;
    myF->ident = maxIdentities;
    myF->similarity = myF->score * 100.0
                      / scoreMax(&sX->datos[myF->xStart], &sY->datos[myF->yStart],
                                 myF->length, POINT);
    myF->seqX = H->seqX;
    myF->seqY = H->seqY;
    myF->block = 0;
    myF->strand = 'r';

    M[min(myF->length, 999)][(int) myF->similarity]++;

    if (myF->length > minLength && myF->similarity > SimTh)
        return 1;
    else
        return 0;
}

int64_t GTHF(hit a1, hit a2) {
    int64_t a1Diag;
    int64_t a2Diag;

    a1Diag = a1.posX - a1.posY;
    a2Diag = a2.posX - a2.posY;

    if (a1Diag > a2Diag)
        return 1;
    else if (a1Diag < a2Diag)
        return 0;
    if (a1.posX > a2.posX)
        return 1;

    return 0;
}

int64_t GTHR(hit a1, hit a2, uint64_t minSeqXLenSeqYLen) {
    int64_t a1Diag;
    int64_t a2Diag;

    a1Diag = a1.posX + a1.posY - minSeqXLenSeqYLen;
    a2Diag = a2.posX + a2.posY - minSeqXLenSeqYLen;

    if (a1Diag > a2Diag)
        return 1;
    else if (a1Diag < a2Diag)
        return 0;
    if (a1.posX > a2.posX)
        return 1;

    return 0;
}

long double computeExpectedValueOfFrag(struct FragFile * frag, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY){
    
    long double lambda = 0.275;
    long double karlin = 0.333;

    long double rawscore = (frag->ident*POINT) - (frag->length - frag->ident)*(POINT);
    

    //printf("Seqx: %"PRIu64", Seqy: %"PRIu64"\n", frag->seqX, frag->seqY);


    long double t_len_x = seqStatsX[frag->seqX].tf.total;
    long double t_len_y = seqStatsY[frag->seqY].tf.total;

    //printf("Real length seqX: %Le, Real length seqY: %Le\n", t_len_x, t_len_y);

    //I verified that it is correct.
    long double expected = karlin*t_len_y*t_len_x*expl(-lambda * rawscore);


    //printf("Len: %"PRIu64" Ident: %"PRIu64"\n Raw: %Le Expected: %Le\n", frag->length, frag->ident, rawscore, expected);
    //getchar();
    frag->evalue = expected;
    return expected;
}




/************************/

