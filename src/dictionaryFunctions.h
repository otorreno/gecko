#ifndef DICTIONARY_FUNCTIONS_H
#define DICTIONARY_FUNCTIONS_H

#define WORD_SIZE 32
//static int BITS_PER_BASE=2;
#define BYTES_IN_WORD 8//(int)ceil(WORD_SIZE/8.*BITS_PER_BASE);
#define SIZE_LOC 20

typedef struct {
    char* seqFile;
    uint64_t* nEntries;
} DictionaryArgs;

/**
 * Compress the word stored in 'buf' using 2 bits per letter
 * The result will be at word 'w'
 */
int seq2word(char* buf, int wsize, word* w);

/**
 * Function to skip the identification line of a fasta sequence
 */
void skipIDLine(FILE *fIn);

/**
 * Function to convert the alphabet letter
 * to an index.
 * coding (a=0,c=1,g=2,t=3,'>'=4 others=9)
 */
int letterToIndex(char c);

/**
 * Function to print in stdout the given compressed word
 */
void showWord(word* w, char *ws);

/**
 * Function to compare two k-mers
 * The function returns 0 if equal, 1 if greater
 * than and -1 otherwise
 */
int wordcmp(unsigned char *w1, unsigned char*w2, int n);

/**
 * Function to be executed by a pthread to calculate
 * a sequence dictionary
 */
void *dictionary(void *a);

/**
 * Function to be executed by a pthread to calculate
 * a sequence dictionary
 */
void *dictionaryWithReverse(void *a);

#endif /* DICTIONARY_FUNCTIONS_H */
