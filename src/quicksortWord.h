#ifndef QUICKSORTWORD_H
#define QUICKSORTWORD_H

#include <inttypes.h>

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements
 */
int64_t psortW(int64_t nproc, wentry *array, int64_t n);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int GTW(wentry a1, wentry a2);

#endif /* QUICKSORTWORD_H */
