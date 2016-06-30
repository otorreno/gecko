#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksortWord.h"

#define SWAP(a, b, t) t=a; a=b; b=t;
#define STRING_SIZE 1024

typedef struct {
    wentry *a;
    int64_t l;
    int64_t r;
    int64_t nth;
} PQSortArgsW;

int64_t partitionW(wentry *a, int64_t l, int64_t r) {
    int64_t i = l;
    int64_t j = r + 1;
    wentry t;

    // l sera el pivote
    // y contendra la mediana de l, r y (l+r)/2
    int64_t mid = (l + r) / 2;

    if (GTW(a[mid], a[r])) {
        SWAP(a[mid], a[r], t);
    }

    if (GTW(a[mid], a[l])) {
        SWAP(a[mid], a[l], t);
    }

    if (GTW(a[l], a[r])) {
        SWAP(a[l], a[r], t);
    }

    while (1) {
        do {
            ++i;
        } while (!GTW(a[i], a[l]) && i <= r);

        do {
            --j;
        } while (GTW(a[j], a[l]) && j >= l);

        if (i >= j) break;

        SWAP(a[i], a[j], t)
    }

    SWAP(a[l], a[j], t)

    return j;
}

int64_t QsortCW(wentry *a, int64_t l, int64_t r) {
    int64_t j;

    if (l < r) {
        // divide and conquer
        j = partitionW(a, l, r);
        //  j=(l+r)/2;
        QsortCW(a, l, j - 1);
        QsortCW(a, j + 1, r);
    }
    return 0;
}

void *PQSortW(void *a) {
    PQSortArgsW *args = (PQSortArgsW *) a;
    if (args->nth > 1) {
        int64_t j = partitionW(args->a, args->l, args->r);
        int64_t np = 1;
        if (args->r - args->l > 0)
            np = (args->nth * (j - args->l)) / (args->r - args->l);
        if (np < 1) np = 1;
        if (np >= args->nth) np = args->nth - 1;
        //printf("%d\t%d (%d)\t%d (%d)\n",args->r-args->l,j-args->l,np,args->r-j,args->nth-np);

        pthread_t *th = (pthread_t *) calloc(2, sizeof(pthread_t));

        PQSortArgsW *nargs = (PQSortArgsW *) calloc(2, sizeof(PQSortArgsW));

        nargs[0].a = args->a;
        nargs[0].l = args->l;
        nargs[0].r = j - 1;
        nargs[0].nth = np;
        pthread_create(th, NULL, PQSortW, (void *) (nargs));

        nargs[1].a = args->a;
        nargs[1].l = j + 1;
        nargs[1].r = args->r;
        nargs[1].nth = args->nth - np;
        pthread_create(th + 1, NULL, PQSortW, (void *) (nargs + 1));

        pthread_join(th[0], NULL);
        pthread_join(th[1], NULL);

        free(th);
        free(nargs);
    } else {
        QsortCW(args->a, args->l, args->r);
    }
    pthread_exit(NULL);
}

int64_t psortW(int64_t nproc, wentry *a, int64_t n) {
    int64_t np = nproc;
    if (np < 1) np = 1;

#ifdef VERBOSE
    printf("Quicksorts\n");
    unsigned long t = timestart();
#endif

    pthread_t th;
    PQSortArgsW args;

    args.a = a;
    args.l = 0;
    args.r = n - 1;
    args.nth = np;
    pthread_create(&th, NULL, PQSortW, (void *) (&args));

    //Wait:
    pthread_join(th, NULL);

#ifdef VERBOSE
    printf("Quicksorts: %lu\n\n",timestop(t));
#endif

    return 0;

}
