#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "comparisonFunctions.h"



int main(int ac, char **av) {
    if (ac < 6) {
        terror("USE: gecko <path_list> <fragsFile.OUT> <Lmin> <SimTh> <WL> <*min_e_value>");
    }

    //Load files to compare
    uint64_t i, n_files, t_alloc, j, curr_comp = 0, k;

    char ** paths_to_files = read_all_vs_all_files(av[1], &n_files, &t_alloc);
    if(n_files < 2) terror("At least two files need to be compared");



    uint64_t n_comps = (n_files * (n_files-1)) / 2;
    fprintf(stdout, "[INFO] N-comps: %"PRIu64"\n", n_comps);

    char *out = av[2];

    uint64_t Lmin = atoi(av[3]);
    uint64_t SimTh = atoi(av[4]);
    int wSize = atoi(av[5]);

    long double e_value = LDBL_MAX; //If no expected value was provided, use Long Double Max which will never happend
    if(ac == 7) e_value = (long double) atof(av[6]);

    uint64_t nEntriesX = 0;
    uint64_t seqXLen = 0;
    
    uint64_t nEntriesY = 0;
    uint64_t seqYLen = 0;

    hashentryF *entriesX;
    hashentryR *entriesY;

    uint64_t nFrags[n_comps];
    FragsFandR frags_list[n_comps];

    FILE * fOut = fopen64(out, "wb");
    if(fOut == NULL) terror("Could not open output frags file");
    
    uint64_t unused_len = 0;
    writeSequenceLength(&unused_len, fOut);
    writeSequenceLength(&unused_len, fOut);


    for(i=0;i<n_files;i++){

        for(j=i+1;j<n_files;j++){

            pthread_t thX, thY;

            DictionaryArgs argsX, argsY;

            nEntriesX = 0;
            seqXLen = 0;
            
            nEntriesY = 0;
            seqYLen = 0;

            argsX.seqFile = strdup(paths_to_files[i]);
            argsX.nEntries = &nEntriesX;
            argsX.seqLen = &seqXLen;

            argsY.seqFile = strdup(paths_to_files[j]);
            argsY.nEntries = &nEntriesY;
            argsY.seqLen = &seqYLen;

            pthread_create(&thX, NULL, dictionary, (void *) (&argsX));
            pthread_create(&thY, NULL, dictionaryWithReverse, (void *) (&argsY));

            pthread_join(thX, (void **) &entriesX);
            pthread_join(thY, (void **) &entriesY);

            free(argsX.seqFile);
            free(argsY.seqFile);

            frags_list[curr_comp] = hitsAndFrags(paths_to_files[i], paths_to_files[j], out, seqXLen, seqYLen, entriesX, nEntriesX, entriesY, nEntriesY, wSize, Lmin, SimTh,
                 &nFrags[curr_comp], argsX.seqStats, argsY.seqStats, e_value);
            
            for(k=0;k<frags_list[curr_comp].t_forward_fragments;k++){
                
                frags_list[curr_comp].fforward[k].seqX = i;                
                frags_list[curr_comp].fforward[k].seqY = j;

                //printf("this frag: %"PRIu64"; %"PRIu64"; %"PRIu64"; \n", frags_list[curr_comp].fforward[k].seqX, frags_list[curr_comp].fforward[k].seqY, frags_list[curr_comp].fforward[k].ident);

                writeFragment(&frags_list[curr_comp].fforward[k], fOut);
                //fwrite(&frags_list[curr_comp].fforward[k], sizeof(struct FragFile), 1, fOut);
                
            }

            for(k=0;k<frags_list[curr_comp].t_reverse_fragments;k++){
                
                
                frags_list[curr_comp].freverse[k].seqX = i;
                frags_list[curr_comp].freverse[k].seqY = j;
                
                writeFragment(&frags_list[curr_comp].freverse[k], fOut);
                //fwrite(&frags_list[curr_comp].freverse[k], sizeof(struct FragFile), 1, fOut);
            }


            curr_comp++;

            for (k = 0; k < nEntriesX; k++) {
                free(entriesX[k].locs);
            }
            free(entriesX);

            for (k = 0; k < nEntriesY; k++) {
                free(entriesY[k].locs);
            }
            free(entriesY);

            free(argsX.seqStats);
            free(argsY.seqStats);

        }
    }

    for(i=0;i<n_comps;i++){
        fprintf(stdout, "[INFO] Fragments in comparison %"PRIu64": %"PRIu64"\n", i, nFrags[i]);
    }

    for(i=0;i<t_alloc;i++){
        free(paths_to_files[i]);
    }
    free(paths_to_files);

    fclose(fOut);
    return 0;
}