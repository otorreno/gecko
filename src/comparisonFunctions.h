#ifndef COMPARISON_FUNCTIONS_H
#define COMPARISON_FUNCTIONS_H

#define min(x, y)    (((x) < (y)) ? (x) : (y))

typedef struct {
    char *seqX;
    char *seqY;
    hit *hits;
    hit **hits_realloc;
    uint64_t nHits;
    uint64_t *nHitsUsed;
    uint64_t Lmin;
    uint64_t SimTh;
    int wSize;
    uint64_t *nFrags;
    uint64_t minSeqLen;
    struct statsHSP * seqStatsX;
    struct statsHSP * seqStatsY;
    long double e_value;
} ComparisonArgs;

/**
 * Get the maximum score possible between the two sequences
 */
uint64_t scoreMax(char *seq, char *seq2, uint64_t len, int point);

/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f);

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f);

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f);

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f);

/**
 * Function to return the sizeof a fragment.
 * Due to architeture issues the value of sizeof built-in
 * function could be different
 */
long int sizeofFragment();

/**
 * Function to calculate the hits between two dictionaries
 * in memory
 */
struct FragFile *hitsAndFrags(char *seqX, char *seqY, char *out, uint64_t seqXLen, uint64_t seqYLen,
                              hashentryF *entriesX,
                              uint64_t nEntriesX, hashentryR *entriesY, uint64_t nEntriesY, int wSize, uint64_t Lmin,
                              uint64_t SimTh, uint64_t *nFrags, struct statsHSP * seqStatsX, struct statsHSP * seqStatsY,
                              long double e_value);

#endif /* COMPARISON_FUNCTIONS_H */
