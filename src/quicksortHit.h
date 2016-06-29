//
// Created by otorreno on 17/06/16.
//

#ifndef QUICKSORTHIT_H
#define QUICKSORTHIT_H

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements
 */
int psortH(int nproc, hit *array, uint64_t n);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int GTH(hit a1, hit a2);

#endif //QUICKSORTHIT_H