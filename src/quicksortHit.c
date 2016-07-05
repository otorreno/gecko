//
// Created by otorreno on 17/06/16.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksortHit.h"

#define SWAP(a, b, t) t=a; a=b; b=t;
#define STRING_SIZE 1024

typedef struct {
    hit *a;
    int64_t l;
    int64_t r;
    int64_t nth;
    uint64_t minSeqLen;
} PQSortArgsH;

// For hits with the FORWARD strand
int64_t partitionHF(hit *a, int64_t l, int64_t r) {
    int64_t i = l;
    int64_t j = r + 1;
    hit t;

    // l sera el pivote
    // y contendra la mediana de l, r y (l+r)/2
    int64_t mid = (l + r) / 2;

    if (GTHF(a[mid], a[r])) {
        SWAP(a[mid], a[r], t);
    }

    if (GTHF(a[mid], a[l])) {
        SWAP(a[mid], a[l], t);
    }

    if (GTHF(a[l], a[r])) {
        SWAP(a[l], a[r], t);
    }

    while (1) {
        do {
            ++i;
        } while (!GTHF(a[i], a[l]) && i <= r);

        do {
            --j;
        } while (GTHF(a[j], a[l]) && j >= l);

        if (i >= j) break;

        SWAP(a[i], a[j], t)
    }

    SWAP(a[l], a[j], t)

    return j;
}

int64_t QsortCHF(hit *a, int64_t l, int64_t r) {
    int64_t j;

    if (l < r) {
        // divide and conquer
        j = partitionHF(a, l, r);
        //  j=(l+r)/2;
        QsortCHF(a, l, j - 1);
        QsortCHF(a, j + 1, r);
    }
    return 0;
}

void *PQSortF(void *a) {
    PQSortArgsH *args = (PQSortArgsH *) a;
    if (args->nth > 1) {
        int64_t j = partitionHF(args->a, args->l, args->r);
        int64_t np = 1;
        if (args->r - args->l > 0)
            np = (args->nth * (j - args->l)) / (args->r - args->l);
        if (np < 1) np = 1;
        if (np >= args->nth) np = args->nth - 1;
        //printf("%d\t%d (%d)\t%d (%d)\n",args->r-args->l,j-args->l,np,args->r-j,args->nth-np);

        pthread_t *th = (pthread_t *) calloc(2, sizeof(pthread_t));

        PQSortArgsH *nargs = (PQSortArgsH *) calloc(2, sizeof(PQSortArgsH));

        nargs[0].a = args->a;
        nargs[0].l = args->l;
        nargs[0].r = j - 1;
        nargs[0].nth = np;
        pthread_create(th, NULL, PQSortF, (void *) (nargs));

        nargs[1].a = args->a;
        nargs[1].l = j + 1;
        nargs[1].r = args->r;
        nargs[1].nth = args->nth - np;
        pthread_create(th + 1, NULL, PQSortF, (void *) (nargs + 1));

        pthread_join(th[0], NULL);
        pthread_join(th[1], NULL);

        free(th);
        free(nargs);
    } else {
        QsortCHF(args->a, args->l, args->r);
    }
    pthread_exit(NULL);
}

int64_t psortHF(int64_t nproc, hit *a, int64_t n) {
    int64_t np = nproc;
    if (np < 1) np = 1;

//    printf("Stage1: Quicksorts\n");
//    unsigned long t = timestart();

    //Quicksort:
    // printf("Quicksort %d\n",n);
    pthread_t th;
    PQSortArgsH args;

    args.a = a;
    args.l = 0;
    args.r = n - 1;
    args.nth = np;
    pthread_create(&th, NULL, PQSortF, (void *) (&args));

    //Wait:
    pthread_join(th, NULL);

//    printf("Stage1: %lu\n\n",timestop(t));

    return 0;

}

/////////////////////////////////////////////////////////////////////
//For hits with the REVERSE strand
int64_t partitionHR(hit *a, int64_t l, int64_t r, uint64_t minSeqXLenSeqYLen) {
    int64_t i = l;
    int64_t j = r + 1;
    hit t;

    // l sera el pivote
    // y contendra la mediana de l, r y (l+r)/2
    int64_t mid = (l + r) / 2;

    if (GTHR(a[mid], a[r], minSeqXLenSeqYLen)) {
        SWAP(a[mid], a[r], t);
    }

    if (GTHR(a[mid], a[l], minSeqXLenSeqYLen)) {
        SWAP(a[mid], a[l], t);
    }

    if (GTHR(a[l], a[r], minSeqXLenSeqYLen)) {
        SWAP(a[l], a[r], t);
    }

    while (1) {
        do {
            ++i;
        } while (!GTHR(a[i], a[l], minSeqXLenSeqYLen) && i <= r);

        do {
            --j;
        } while (GTHR(a[j], a[l], minSeqXLenSeqYLen) && j >= l);

        if (i >= j) break;

        SWAP(a[i], a[j], t)
    }

    SWAP(a[l], a[j], t)

    return j;
}

int64_t QsortCHR(hit *a, int64_t l, int64_t r, uint64_t minSeqXLenSeqYLen) {
    int64_t j;

    if (l < r) {
        // divide and conquer
        j = partitionHR(a, l, r, minSeqXLenSeqYLen);
        //  j=(l+r)/2;
        QsortCHR(a, l, j - 1, minSeqXLenSeqYLen);
        QsortCHR(a, j + 1, r, minSeqXLenSeqYLen);
    }
    return 0;
}

void *PQSortR(void *a) {
    PQSortArgsH *args = (PQSortArgsH *) a;
    if (args->nth > 1) {
        int64_t j = partitionHR(args->a, args->l, args->r, args->minSeqLen);
        int64_t np = 1;
        if (args->r - args->l > 0)
            np = (args->nth * (j - args->l)) / (args->r - args->l);
        if (np < 1) np = 1;
        if (np >= args->nth) np = args->nth - 1;
//        printf("%"PRId64"\t%"PRId64" %"PRId64"\t(%"PRId64") (%"PRId64")\n",args->l,args->r,j,np,args->nth-np);

        pthread_t *th = (pthread_t *) calloc(2, sizeof(pthread_t));

        PQSortArgsH *nargs = (PQSortArgsH *) calloc(2, sizeof(PQSortArgsH));

        nargs[0].a = args->a;
        nargs[0].l = args->l;
        nargs[0].r = j - 1;
        nargs[0].nth = np;
        nargs[0].minSeqLen = args->minSeqLen;
        pthread_create(th, NULL, PQSortR, (void *) (nargs));

        nargs[1].a = args->a;
        nargs[1].l = j + 1;
        nargs[1].r = args->r;
        nargs[1].nth = args->nth - np;
        nargs[1].minSeqLen = args->minSeqLen;
        pthread_create(th + 1, NULL, PQSortR, (void *) (nargs + 1));

        pthread_join(th[0], NULL);
        pthread_join(th[1], NULL);

        free(th);
        free(nargs);
    } else {
        QsortCHR(args->a, args->l, args->r, args->minSeqLen);
    }
    pthread_exit(NULL);
}

int64_t psortHR(int64_t nproc, hit *a, int64_t n, uint64_t minSeqXLenSeqYLen) {
    int64_t np = nproc;
    if (np < 1) np = 1;

//    printf("Stage1: Quicksorts\n");
//    unsigned long t = timestart();

    //Quicksort:
    // printf("Quicksort %d\n",n);
    pthread_t th;
    PQSortArgsH args;

    args.a = a;
    args.l = 0;
    args.r = n - 1;
    args.nth = np;
    args.minSeqLen = minSeqXLenSeqYLen;
    pthread_create(&th, NULL, PQSortR, (void *) (&args));

    //Wait:
    pthread_join(th, NULL);

//    printf("Stage1: %lu\n\n",timestop(t));

    return 0;

}
