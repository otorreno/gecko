#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "comparisonFunctions.h"

int main(int ac, char** av){
    if(ac!=7){
        terror("USE: gecko seqFileX.IN seqFileX.IN fragsFile.OUT Lmin SimTh WL");
    }

    char *seqX = av[1];
    char *seqY = av[2];
    char *out = av[3];

    uint64_t Lmin = atoi(av[4]);
    uint64_t SimTh = atoi(av[5]);
    int wSize = atoi(av[6]);

    uint64_t nEntriesX = 0;
    hashentry *entriesX;
    uint64_t nEntriesY = 0;
    hashentry *entriesY;

    uint64_t hitsInBuf;

    pthread_t thX,thY;
    DictionaryArgs argsX,argsY;

    argsX.seqFile=strdup(seqX);
    argsX.nEntries=&nEntriesX;
    pthread_create(&thX,NULL,dictionary,(void*)(&argsX));

    argsY.seqFile=strdup(seqY);
    argsY.nEntries=&nEntriesY;
    pthread_create(&thY,NULL,dictionaryWithReverse,(void*)(&argsY));

    //Wait:
    pthread_join(thX,(void **)&entriesX);
    pthread_join(thY,(void **)&entriesY);

    free(argsX.seqFile);
    free(argsY.seqFile);

//    char wordString[33];
//	wordString[32] = '\0';
//    uint64_t i,j;
//	for(i=0;i<nEntriesX;i++){
//		showWord(&entriesX[i].w, wordString);
//		fprintf(stdout, "Words(%" PRIu64 "):%s Repetitions: %" PRIu64 " Positions: ", i, wordString, entriesX[i].num);
//        for(j=0;j<entriesX[i].num;j++){
//            fprintf(stdout,"(%" PRIu64 ",%" PRIu64 ") ",entriesX[i].locs[j].pos,entriesX[i].locs[j].seq);
//        }
//        fprintf(stdout, "\n");
//	}

    //Hits, SortHits and filterHits
    struct FragFile *fragsBuf = hitsAndFrags(seqX, seqY, out, entriesX, nEntriesX, entriesY, nEntriesY, wSize, &hitsInBuf, Lmin, SimTh);
    uint64_t i;

//    for(i=0;i<hitsInBuf; i++) {
//        fprintf(stdout, "d=%-7" PRId64 " pX=%-7" PRIu64 " pY=%-7" PRIu64 " seqX=%-7" PRIu64 " seqY=%-7" PRIu64 "\n", hBuf[i].diag, hBuf[i].posX, hBuf[i].posY, hBuf[i].seqX, hBuf[i].seqY);
//    }
	for(i=0;i<nEntriesX;i++){
        free(entriesX[i].locs);
    }
    free(entriesX);
    for(i=0;i<nEntriesY;i++){
        free(entriesY[i].locs);
    }
    free(entriesY);



    free(fragsBuf);

    return 0;
}