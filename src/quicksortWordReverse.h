#ifndef QUICKSORTWORDREVERSE_H
#define QUICKSORTWORDREVERSE_H

#include <inttypes.h>

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements
 */
int64_t psortWR(int64_t nproc, wentryR *array, int64_t n);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int GTWR(wentryR a1, wentryR a2);

#endif /* QUICKSORTWORDREVERSE_H */
