#ifndef QUICKSORT_H
#define QUICKSORT_H

#include <inttypes.h>

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements
 */
int psort(int nproc, BaseType* array, uint64_t n);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int GT(BaseType a1, BaseType a2);

#endif /* QUICKSORT_H */
