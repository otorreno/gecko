//
// Created by otorreno on 17/06/16.
//

#ifndef QUICKSORTHIT_H
#define QUICKSORTHIT_H

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements.
 * For hits with the FORWARD strand
 */
int64_t psortHF(int64_t nproc, hit *array, int64_t n);

/**
 * Function to order 'array' in memory with
 * a thread number of 'nproc'. 'array' has 'n' elements.
 * For hits with the REVERSE strand
 */
int64_t psortHR(int64_t nproc, hit *array, int64_t n, uint64_t minSeqXLenSeqYLen);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * For hits with the FORWARD strand
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int64_t GTHF(hit a1, hit a2);

/**
 * Function to determine if object 1 is strictly greater than 2.
 * For hits with the REVERSE strand
 * Returns 1 if true, 0 otherwise
 * Has to be defined in the main procedure
 */
int64_t GTHR(hit a1, hit a2, uint64_t minSeqXLenSeqYLen);

#endif //QUICKSORTHIT_H
