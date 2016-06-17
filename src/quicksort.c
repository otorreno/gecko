#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "quicksort.h"

#define SWAP(a,b,t) t=a; a=b; b=t;
#define STRING_SIZE 1024

typedef struct {
	BaseType* a;
	int l;
	int r;
	int nth;
} PQSortArgs;

void assertNotNull(void* p, char* msg){
	if(p==NULL){
		terror(msg);
	}
}

void assertIntEQ(int a, int b, char* msg){
	if(a!=b){
		terror(msg);
	}
}

void assertIntGE(int a, int b, char* msg){
	if(a<b){
		terror(msg);
	}
}

int partition(BaseType* a, int l, int r) {
   int i=l;
   int j=r+1;
   BaseType t;

   // l sera el pivote
   // y contendra la mediana de l, r y (l+r)/2
   int mid = (l+r)/2;

   if(GT(a[mid],a[r])) {
		 SWAP(a[mid],a[r],t);
   }

   if(GT(a[mid],a[l])) {
		 SWAP(a[mid],a[l],t);
   }

   if(GT(a[l],a[r])) {
		 SWAP(a[l],a[r],t);
	 }

	while (1) {
		do{
			++i;
		}while( !GT(a[i],a[l]) && i <= r );

		do{
			--j;
		}while( GT(a[j],a[l]) && j >= l);

		if( i >= j ) break;

		SWAP(a[i],a[j],t)
	}

	SWAP(a[l],a[j],t)

	return j;
}

int QsortC(BaseType* a, int l,int r) {
   int j;

	if( l < r ) {
 	// divide and conquer
       j = partition( a, l, r);
       //  j=(l+r)/2;
       QsortC( a, l, j-1);
       QsortC( a, j+1, r);
   }
   return 0;
}

void *PQSort(void* a){
	PQSortArgs *args=(PQSortArgs*)a;
	if(args->nth>1){
		int tmp;
		int j = partition(args->a,args->l,args->r);
		int np=1;
		if(args->r - args->l > 0)
			np = (args->nth*(j-args->l))/(args->r-args->l);
		if(np<1) np=1;
		if(np>=args->nth) np=args->nth-1;
		//printf("%d\t%d (%d)\t%d (%d)\n",args->r-args->l,j-args->l,np,args->r-j,args->nth-np);

		pthread_t* th = (pthread_t*)calloc(2,sizeof(pthread_t));
		assertNotNull(th,"calloc");

		PQSortArgs* nargs = (PQSortArgs*)calloc(2,sizeof(PQSortArgs));
		assertNotNull(args,"calloc");

		nargs[0].a=args->a;
		nargs[0].l=args->l;
		nargs[0].r=j-1;
		nargs[0].nth=np;
		tmp=pthread_create(th,NULL,PQSort,(void*)(nargs));
		assertIntEQ(tmp,0,"pthread_create");

		nargs[1].a=args->a;
		nargs[1].l=j+1;
		nargs[1].r=args->r;
		nargs[1].nth=args->nth-np;
		tmp=pthread_create(th+1,NULL,PQSort,(void*)(nargs+1));
		assertIntEQ(tmp,0,"pthread_create");

		tmp=pthread_join(th[0],NULL);
		assertIntEQ(tmp,0,"pthread_join");
		tmp=pthread_join(th[1],NULL);
		assertIntEQ(tmp,0,"pthread_join");

		free(th);
		free(nargs);
	}else{
		QsortC(args->a,args->l,args->r);
	}
	pthread_exit(NULL);
}

unsigned long timestart(){
	struct timeval tv;

	gettimeofday(&tv,NULL);

	return (tv.tv_usec/1000) + (tv.tv_sec*1000);
}

unsigned long timestop(unsigned long start){
	struct timeval tv;

	gettimeofday(&tv,NULL);

	return (tv.tv_usec/1000) + (tv.tv_sec*1000) - start;
}

int psort(int nproc, BaseType* a, uint64_t n){
    fprintf(stdout, "en psort\n");
	int tmp;
	int np=nproc;
	if(np<1) np=1;

	printf("Stage1: Quicksorts\n");
	unsigned long t = timestart();

    //Quicksort:
    // printf("Quicksort %d\n",n);
    pthread_t th;
    PQSortArgs args;

    args.a=a;
    args.l=0;
    args.r=n-1;
    args.nth=np;
    tmp=pthread_create(&th,NULL,PQSort,(void*)(&args));
    assertIntEQ(tmp,0,"pthread_create");

    //Wait:
    tmp=pthread_join(th,NULL);
    assertIntEQ(tmp,0,"pthread_join");

	printf("Stage1: %lu\n\n",timestop(t));

	return 0;

}
