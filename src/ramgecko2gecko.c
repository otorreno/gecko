/*********

File		ramgecko2gecko.c
Author		EPW <estebanpw@uma.es>
Description	Parses the .frags results of a GECKO-RAM comparison to the old GECKO coordinates version
	

USAGE		<comparison.frags> 	A comparison file produced by GECKO
			<comparison_fixed.frags>	Output file name
			<non_sorted_index.IND2>	Non sorted index for the query produced by mgReadsIndex
			<non_sorted_index.IND2> Non sorted index for the genomes
			<non_sorted_index.IND2> Non sorted index for the reversed genomes

**********/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>

#include "commonFunctions.h"
#include "comparisonFunctions.h"
#include "structs.h"

#define MAXFRAGS 8000
#define MAXPATH 1000



void fixCoordinates(uint64_t * yStart, uint64_t * yEnd, int genLength);
struct rIndex2* loadReadsIndex(char *filename, uint64_t *nReads);



int main(int argc, char **av){
	FILE * fragfile, * out;
	if(argc != 5)terror("USE: ramgecko2gecko <comparison.frags> <comparison_fixed.frags> <query_non_sorted_index.IND2> <db_non_sorted_index.IND2>\n");
	
	//Open frag file
	fragfile = fopen64(av[1], "rb");
	if(fragfile==NULL)terror("Could not open fragfile 1.\n");
	
	//Open output file
	out = fopen64(av[2], "wb");
	if(out==NULL)terror("Could not open output file.\n");
	 
	uint64_t lengthX, lengthY;
	readSequenceLength(&lengthX, fragfile);
	readSequenceLength(&lengthY, fragfile);

	writeSequenceLength(&lengthX, out);
	writeSequenceLength(&lengthY, out);
	
	struct rIndex2 *RI;
	uint64_t nReads;
	RI = loadReadsIndex(av[3], &nReads);
	fprintf(stdout, "Total of %"PRIu64" reads.\n", nReads);
	
	uint64_t nGenomes = 0;

	struct rIndex2 *RI2;
	RI2 = loadReadsIndex(av[4], &nGenomes);

	uint64_t * RI3 = (uint64_t *) malloc(nGenomes * sizeof(uint64_t));
	if(RI3 == NULL) terror("Could not allocate memory for the third index");

	//Calculate reverse LAC
	uint64_t i;
	RI3[0] = 0;
	for(i=1;i<nGenomes;i++){
		RI3[i] = RI3[i-1] + RI2[nGenomes-i].rLen;
	}

	//Read first frag
	struct FragFile countFrag;
	readFragment(&countFrag, fragfile);

	uint64_t old_seq_y, new_seq_y;
	uint64_t old_yStart, old_yEnd, new_yStart, new_yEnd;


	while(!feof(fragfile)){
		//read fragment		
		readFragment(&countFrag, fragfile);

		//printf("%"PRIu64" -> %"PRIu64", %"PRIu64"\n", countFrag.seqY, countFrag.yStart, countFrag.yEnd);

		//Add the "*"
		countFrag.xStart += countFrag.seqX;
		countFrag.xEnd += countFrag.seqX;
	
		if(countFrag.strand == 'r'){
	
			old_seq_y = countFrag.seqY;
			old_yStart = countFrag.yStart;
			old_yEnd = countFrag.yEnd;

			//Convert to local coordinates
			new_yStart = old_yStart - RI2[old_seq_y].Lac + 1;
			new_yEnd = old_yEnd - RI2[old_seq_y].Lac + 1;

			//Invert
			fixCoordinates(&new_yStart, &new_yEnd, RI2[old_seq_y].rLen);

			//Calculate new sequence position
			new_seq_y = nGenomes - 1 - old_seq_y;

			new_yStart += new_seq_y;
			new_yEnd += new_seq_y;

			//Add the accumulated length in reverse
			new_yStart += RI3[new_seq_y] - 1;
			new_yEnd += RI3[new_seq_y] - 1;
		
			//And write to file
			countFrag.seqY = new_seq_y;
			countFrag.yStart = new_yStart;
			countFrag.yEnd = new_yEnd;
		}else{
			countFrag.yStart += countFrag.seqY;
	                countFrag.yEnd += countFrag.seqY;
		}

		//Recalculate diagonal
		countFrag.diag = countFrag.xStart - countFrag.yStart;
		writeFragment(&countFrag, out);
		
	}
	
	
	fclose(fragfile);
	fclose(out);
	free(RI);
	free(RI2);
	free(RI3);
	return 0;
}




/*
Gecko out

Posx-Posy:	-529966
xStart:		192
yStart:		530158
xEnd:		244
yEnd:		530210
length:		53
ident:		46
score:		156
similarity:	74
seqX:		32
seqY:		35
Synteny:	0
strand:		r


*/





//invert frags
void fixCoordinates(uint64_t * yStart, uint64_t * yEnd, int genLength){

	*yStart = (genLength +1) - *yStart;
	*yEnd = (genLength +1) - *yEnd;
	
}
struct rIndex2* loadReadsIndex(char *filename, uint64_t *nReads){
	struct rIndex2 *RR;
	uint64_t nR=0,i;
	FILE *f;
	uint64_t fsize;

    if ((f=fopen(filename,"rb"))==NULL) terror("Could not open index input file");

	fseeko(f, 0, SEEK_END);
	fsize = ftello(f);
	rewind(f);
	nR = fsize/sizeof(struct rIndex2);

    if ((RR =(struct rIndex2*) calloc(nR,sizeof(struct rIndex2)))==NULL) terror("Could not allocate index");

	for (i=0; i<nR; i++){
		if(0 == fread(&RR[i],sizeof(struct rIndex2),1,f)) break;
	}
	fclose(f);
	(*nReads) = nR;
	return RR;
}
/*
Printing properties of a frag
		printf("Posx-Posy:	%"PRId64"\n",countFrag.diag);
		printf("xStart:		%"PRIu64"\n",countFrag.xStart);
		printf("yStart:		%"PRIu64"\n",countFrag.yStart);
		printf("xEnd:		%"PRIu64"\n",countFrag.xEnd);
		printf("yEnd:		%"PRIu64"\n",countFrag.yEnd);
		printf("length:		%"PRIu64"\n",countFrag.length);
		printf("ident:		%"PRIu64"\n",countFrag.ident);
		printf("score:		%"PRIu64"\n",countFrag.score);
		printf("similarity:	%.f\n",countFrag.similarity);
		printf("seqX:		%"PRIu64" which is %s\n",countFrag.seqX, RI[(int)countFrag.seqX].id);
		printf("seqY:		%"PRIu64" which is %s\n",countFrag.seqY, gL[(int)countFrag.seqY].name);
		printf("Synteny:	%"PRId64"\n",countFrag.block);
		printf("strand:		%c\n",countFrag.strand);	
		*/

