#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "comparisonFunctions.h"

int main(int ac, char **av) {
    if (ac != 7) {
        terror("USE: gecko seqFileX.IN seqFileX.IN fragsFile.OUT Lmin SimTh WL");
    }

    char *seqX = av[1];
    char *seqY = av[2];
    char *out = av[3];

    uint64_t Lmin = atoi(av[4]);
    uint64_t SimTh = atoi(av[5]);
    int wSize = atoi(av[6]);

    uint64_t nEntriesX = 0;
    uint64_t seqXLen = 0;
    hashentry *entriesX;
    uint64_t nEntriesY = 0;
    uint64_t seqYLen = 0;
    hashentry *entriesY;

    uint64_t nFrags;

    pthread_t thX, thY;
    DictionaryArgs argsX, argsY;

    argsX.seqFile = strdup(seqX);
    argsX.nEntries = &nEntriesX;
    argsX.seqLen = &seqXLen;
    pthread_create(&thX, NULL, dictionary, (void *) (&argsX));

    argsY.seqFile = strdup(seqY);
    argsY.nEntries = &nEntriesY;
    argsY.seqLen = &seqYLen;
    pthread_create(&thY, NULL, dictionaryWithReverse, (void *) (&argsY));

    //Wait:
    pthread_join(thX, (void **) &entriesX);
    pthread_join(thY, (void **) &entriesY);

    free(argsX.seqFile);
    free(argsY.seqFile);

    //Hits, SortHits and filterHits
    hitsAndFrags(seqX, seqY, out, seqXLen, seqYLen, entriesX, nEntriesX, entriesY, nEntriesY, wSize, Lmin, SimTh,
                 &nFrags);
    uint64_t i;

    for (i = 0; i < nEntriesX; i++) {
        free(entriesX[i].locs);
    }
    free(entriesX);
    for (i = 0; i < nEntriesY; i++) {
        free(entriesY[i].locs);
    }
    free(entriesY);

    return 0;
}