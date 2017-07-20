#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "comparisonFunctionsAll.h"

void write_header(FILE * f, uint64_t sx_len, uint64_t sy_len);
void write_csvs(char * generic_path, uint64_t n_comparison, struct FragFile * f, uint64_t n_frags_f, struct FragFile * r, uint64_t n_frags_r, uint64_t xlen, uint64_t ylen, uint64_t lmin, uint64_t simmin, int kmer);

int main(int ac, char **av) {
    if (ac < 6) {
        terror("USE: gecko <path_list> <fragsFile.OUT> <Lmin> <SimTh> <WL> <*min_e_value>");
    }

    //Load files to compare
    uint64_t i, n_files, j, curr_comp = 0, k, bijectve_i_j = 0;

    char ** paths_to_files = read_all_vs_all_files(av[1], &n_files);
    if(n_files < 2) terror("At least two files need to be compared");
    unsigned char wrote_lengths[n_files]; //To keep track of which lengths were written or not
    for(i=0;i<n_files;i++) { wrote_lengths[i] = 0; }



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
    char path_for_lengths[READLINE];
    strcpy(path_for_lengths, out);
    strcat(path_for_lengths, ".lengths");
    FILE * fOutLengths = fopen64(path_for_lengths, "wb");
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
                 
            write_csvs("multicomp", bijectve_i_j, frags_list[curr_comp].fforward, frags_list[curr_comp].t_forward_fragments, frags_list[curr_comp].freverse, frags_list[curr_comp].t_reverse_fragments, *argsX.seqLen, *argsY.seqLen, Lmin, SimTh, wSize);
            bijectve_i_j++;
            
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


            //This will write the first sequence's length
            if(wrote_lengths[i] == 0){
                fwrite(&argsX.seqStats->tf.total, sizeof(uint64_t), 1, fOutLengths);
                wrote_lengths[i] = 1;
            }
            //This will write the rest throughout the second loop
            if(wrote_lengths[j] == 0){
                fwrite(&argsY.seqStats->tf.total, sizeof(uint64_t), 1, fOutLengths);
                wrote_lengths[j] = 1;
            }

            free(argsY.seqStats);

        }
    }

    for(i=0;i<n_comps;i++){
        fprintf(stdout, "[INFO] Fragments in comparison %"PRIu64": %"PRIu64"\n", i, nFrags[i]);
    }

    for(i=0;i<n_files;i++){
        free(paths_to_files[i]);
    }
    free(paths_to_files);

    fclose(fOut);
    fclose(fOutLengths);
    return 0;
}

void write_header(FILE * f, uint64_t sx_len, uint64_t sy_len){
    fprintf(f, "CSV file\n");
    fprintf(f, "[Jul.17 -- < estebanpw@uma.es >\n");
    fprintf(f, "SeqX filename	: DATA1.dat\n");
    fprintf(f, "SeqY filename	: DATA2.dat\n");
    fprintf(f, "SeqX name	: S1\n");
    fprintf(f, "SeqY name	: S2\n");
    fprintf(f, "SeqX length	: %"PRIu64"\n", sx_len);
    fprintf(f, "SeqY length	: %"PRIu64"\n", sy_len);
    fprintf(f, "Min.fragment.length	: 0\n");
    fprintf(f, "Min.Identity	: 0.0\n");
    fprintf(f, "Total hits	: 0\n");
    fprintf(f, "Total hits (used)	: 0\n");
    fprintf(f, "Total fragments	: 0\n");
    fprintf(f, "Total CSBs:	: 0\n");
    fprintf(f, "frag/CSB,xStart,yStart,xEnd,YEnd,strand,block,length,score,ident,similarity,identity,geneX,geneY\n");
    fprintf(f, "========================================================\n");
}

void write_csvs(char * generic_path, uint64_t n_comparison, struct FragFile * f, uint64_t n_frags_f, struct FragFile * r, uint64_t n_frags_r, uint64_t xlen, uint64_t ylen, uint64_t lmin, uint64_t simmin, int kmer){

    FILE * out;
    char buffer[2048]; buffer[0]='\0';
    sprintf(buffer, "%s_%"PRIu64"_%"PRIu64"_%"PRIu64"_%d", generic_path, n_comparison, lmin, simmin, kmer);
    sprintf(buffer, "%s.csv", buffer);
    out = fopen64(buffer, "wt");
    if(out == NULL) terror("Could not open output csv file");


    write_header(out, xlen, ylen);


	double similarity, likeness;
	uint64_t i;
	for(i=0;i<n_frags_f;i++){
		similarity=(((double)f[i].score)/((double)f[i].length*4.0));
		likeness=(((double)f[i].ident)/((double)f[i].length));
		

		fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRId64",%"PRIu64",%.2f,%.2f,0,0,%"PRIu64",%"PRIu64"\n", f[i].xStart, f[i].yStart, f[i].xEnd, f[i].yEnd,f[i].strand, f[i].block, f[i].length, f[i].score, f[i].ident, similarity, likeness, (uint64_t)0, (uint64_t)1);
	}
	
	for(i=0;i<n_frags_r;i++){
		similarity=(((double)r[i].score)/((double)r[i].length*4.0));
		likeness=(((double)r[i].ident)/((double)r[i].length));
		


		fprintf(out, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRId64",%"PRIu64",%.2f,%.2f,0,0,%"PRIu64",%"PRIu64"\n", r[i].xStart, r[i].yStart, r[i].xEnd, r[i].yEnd,r[i].strand, r[i].block, r[i].length, r[i].score, r[i].ident, similarity, likeness, (uint64_t)0, (uint64_t)1);
	}
	fclose(out);
}
