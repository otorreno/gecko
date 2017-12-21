/*

filterFrags.c

Exports frags to CSV filtering by thresholds

Use: filterFrags file.frags length similarity

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#define  MAX_LEVELS  900000


int main(int ac,char** av){

	
	if(ac<4){
		printf("Use: filterFrags <file.frags> <length> <similarity>\n");
		exit(-1);
	}
	
	// Read fragments
	struct FragFile frag;
	uint64_t xtotal,ytotal;
	//f=readFragmentsv2(av[1],&nf,&xtotal,&ytotal);

	FILE * fFrags = fopen(av[1], "rb");
	if(fFrags == NULL){ fprintf(stderr, "Could not open input file\n"); exit(-1); }

	readSequenceLength(&xtotal, fFrags);
	readSequenceLength(&ytotal, fFrags);
	

	

	uint64_t min_l = (uint64_t) atoi(av[2]);
	double min_sim = (double) atof(av[3])/100;
	
	/******************************************/
	
	// Print header. Frags file info
	printf("Total CSB: 0\n");
	printf("========================================================\n");
	
	printf("Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");
	
	double similarity,likeness;
				
	while(!feof(fFrags)){
		readFragment(&frag, fFrags);

		similarity=(((double)frag.score)/((double)frag.length*4.0));
		likeness=(((double)frag.ident)/((double)frag.length));
		
		if(frag.strand=='r'){
			frag.yStart = ytotal - frag.yStart - 1;
			frag.yEnd = ytotal - frag.yEnd - 1;
		}
		
		
		if(similarity >= min_sim && (uint64_t)frag.length >= min_l){
			printf("Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,%"PRIu64",%"PRIu64"\n",frag.xStart,frag.yStart,frag.xEnd,frag.yEnd,frag.strand,frag.block,frag.length,frag.score,frag.ident,similarity,likeness,frag.seqX,frag.seqY);
		}
				
	}

	fclose(fFrags);
	
	return 0;
}
