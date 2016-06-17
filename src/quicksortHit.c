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

#define SWAP(a,b,t) t=a; a=b; b=t;
#define STRING_SIZE 1024

typedef struct {
    hit* a;
    int l;
    int r;
    int nth;
} PQSortArgsH;

int partitionH(hit* a, int l, int r) {
    int i=l;
    int j=r+1;
    hit t;

    // l sera el pivote
    // y contendra la mediana de l, r y (l+r)/2
    int mid = (l+r)/2;

    if(GTH(a[mid],a[r])) {
        SWAP(a[mid],a[r],t);
    }

    if(GTH(a[mid],a[l])) {
        SWAP(a[mid],a[l],t);
    }

    if(GTH(a[l],a[r])) {
        SWAP(a[l],a[r],t);
    }

    while (1) {
        do{
            ++i;
        }while( !GTH(a[i],a[l]) && i <= r );

        do{
            --j;
        }while( GTH(a[j],a[l]) && j >= l);

        if( i >= j ) break;

        SWAP(a[i],a[j],t)
    }

    SWAP(a[l],a[j],t)

    return j;
}

int QsortCH(hit* a, int l,int r) {
    int j;

    if( l < r ) {
        // divide and conquer
        j = partitionH( a, l, r);
        //  j=(l+r)/2;
        QsortCH( a, l, j-1);
        QsortCH( a, j+1, r);
    }
    return 0;
}

void *PQSort(void* a){
    PQSortArgsH *args=(PQSortArgsH*)a;
    if(args->nth>1){
        int j = partitionH(args->a,args->l,args->r);
        int np=1;
        if(args->r - args->l > 0)
            np = (args->nth*(j-args->l))/(args->r-args->l);
        if(np<1) np=1;
        if(np>=args->nth) np=args->nth-1;
        //printf("%d\t%d (%d)\t%d (%d)\n",args->r-args->l,j-args->l,np,args->r-j,args->nth-np);

        pthread_t* th = (pthread_t*)calloc(2,sizeof(pthread_t));

        PQSortArgsH* nargs = (PQSortArgsH*)calloc(2,sizeof(PQSortArgsH));

        nargs[0].a=args->a;
        nargs[0].l=args->l;
        nargs[0].r=j-1;
        nargs[0].nth=np;
        pthread_create(th,NULL,PQSort,(void*)(nargs));

        nargs[1].a=args->a;
        nargs[1].l=j+1;
        nargs[1].r=args->r;
        nargs[1].nth=args->nth-np;
        pthread_create(th+1,NULL,PQSort,(void*)(nargs+1));

        pthread_join(th[0],NULL);
        pthread_join(th[1],NULL);

        free(th);
        free(nargs);
    }else{
        QsortCH(args->a,args->l,args->r);
    }
    pthread_exit(NULL);
}

int psortH(int nproc, hit* a, uint64_t n){
    int np=nproc;
    if(np<1) np=1;

    printf("Stage1: Quicksorts\n");
    unsigned long t = timestart();

    //Quicksort:
    // printf("Quicksort %d\n",n);
    pthread_t th;
    PQSortArgsH args;

    args.a=a;
    args.l=0;
    args.r=n-1;
    args.nth=np;
    pthread_create(&th,NULL,PQSort,(void*)(&args));

    //Wait:
    pthread_join(th,NULL);

    printf("Stage1: %lu\n\n",timestop(t));

    return 0;

}
