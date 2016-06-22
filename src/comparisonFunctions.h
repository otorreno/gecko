#ifndef COMPARISON_FUNCTIONS_H
#define COMPARISON_FUNCTIONS_H

#define min(x,y)    (((x) < (y)) ? (x) : (y))

typedef struct {
    char *seqX;
    char *seqY;
    hit *hits;
    uint64_t nHits;
    uint64_t *nHitsUsed;
    uint64_t Lmin;
    uint64_t SimTh;
    int wSize;
    uint64_t* nFrags;
} ComparisonArgs;

/**
 * Read a new hash entry from 'f' with no more occurences
 * than 'freqThr'
 */
int readHashEntry(hashentry *h, FILE *f, uint64_t freqThr);

/**
 * Read the ocurrences of the given hash entry from the
 * ocurrences file 'f' and stored them in 'pos'
 */
void loadWordOcurrences(hashentry he, location** pos, FILE** f);

/**
 * Get the maximum score possible between the two sequences
 */
unsigned long scoreMax(char *seq, char *seq2, uint64_t len, int point);

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
struct FragFile *hitsAndFrags(char *seqX, char *seqY, char *out, uint64_t seqXLen, uint64_t seqYLen, hashentry *entriesX,
                              uint64_t nEntriesX, hashentry *entriesY, uint64_t nEntriesY, int wSize, uint64_t Lmin,
                              uint64_t SimTh, uint64_t *nFrags);

#endif /* COMPARISON_FUNCTIONS_H */
