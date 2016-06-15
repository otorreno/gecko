#ifndef DICTIONARY_FUNCTIONS_H
#define DICTIONARY_FUNCTIONS_H

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
 * Function to load the sequence.
 * The code was present inside words.c program
 */
void loadSequence(char *fileName, char *seq, uint64_t *Tot);

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
 * Function to destroy the binary tree
 */
void destroy_tree(node **leaf);

/**
 * Insert the given word 'key' in the binary tree
 */
int insert(char *key, int wsize, location loc, node **leaf, node *nodePool, int *actualNodes, list_node *locationPool, int *actualLocations);

/**
 * Write the binary tree in order (left,
 * root, right)
 */
void write(node *leaf, FILE* fOut);

#endif /* DICTIONARY_FUNCTIONS_H */
