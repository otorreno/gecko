#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

/**
 * Print the error message 's' and exit(-1)
 */
void terror(char *s);

unsigned long timestart();

unsigned long timestop(unsigned long start);

/**
 * Read the sequence from disk and store it in a list of Sequence struct
 * n: sequence length
 * ns: number of nodes in the list
 */
struct Sequence*LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nStruct, uint64_t *nSeqs, int fAst);

/**
 * Get the value of the sequence in a given position of the list node ns
 */
char getValue(struct Sequence *s, uint64_t pos, int ns);

/**
 * Get the length of the sequence 's'
 */
long getSeqLength(struct Sequence *s, uint64_t start, int ns);

#endif /* COMMON_FUNCTIONS_H */
