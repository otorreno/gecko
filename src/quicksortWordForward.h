#ifndef QUICKSORTWORDFORWARD_H
#define QUICKSORTWORDFORWARD_H

#include <inttypes.h>

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements
 */
int64_t psortWF(int64_t nproc, wentryF *array, int64_t n);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int GTWF(wentryF a1, wentryF a2);

#endif /* QUICKSORTWORDFORWARD_H */
