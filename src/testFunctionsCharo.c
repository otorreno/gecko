#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <math.h>


#include "fragmentv3.h"
#include "fragmentv2.h"
#include "lista.h"
#include "comparisonFunctions.h"
#define  MAX_LEVELS  900000

void terror(char*);
void readFragmentRaw(struct FragFile *, FILE *);
void readSequenceLengthRaw(uint64_t *, FILE *);
void writeSequenceLengthRaw(uint64_t *, FILE *);

int main(int ac,char** av){
	
	// Read fragments
	FILE *fichero;
	struct FragFile *f, frag;
	int nf = 0; // Number of fragments
	int n = 0;
	uint64_t xtotal,ytotal;
	
	f = readFragmentsv2(av[1], &nf,&xtotal,&ytotal); 
	fichero = fopen(av[2], "wb");
	if(fichero == NULL){
		terror("Abriendo fichero de salida");
	}

	fwrite(&xtotal,sizeof(uint64_t), 1, fichero);
	fwrite(&ytotal,sizeof(uint64_t), 1, fichero);
	for(n ;n < nf; n++){
		writeFragmentRaw(&f[n], fichero);	
	}

	fclose(fichero);
	fichero = fopen(av[2], "rb"); //Se vuelve a abrir el fichero, pero esta vez de lectura
	if(fichero == NULL){
		terror("Abriendo fichero de salida");
	}

	fread(&xtotal,sizeof(uint64_t), 1, fichero);
	fread(&ytotal,sizeof(uint64_t), 1, fichero);
	//fprintf(stdout,"%" PRIu64"\n", xtotal);
	//fprintf(stdout,"%" PRIu64"\n", ytotal);

	while(!feof(fichero)){
        readFragmentRaw(&frag, fichero); 
        /*printf( "diag: %" PRId64   " xStart, ystart: %" PRIu64" %"PRIu64" xEnd yEnd: %" PRIu64" %"PRIu64     "  length: %" PRId64    " ident: %" PRId64    
		" score: %" PRId64   "   similarity: %f"    "  seqX, seqY: %" PRIu64" %"PRIu64   "   strand: %c \n", 
		frag.diag, frag.xStart, frag.yStart, frag.xEnd, frag.yEnd, frag.length, frag.ident, frag.score, frag.similarity, frag.seqX, frag.seqY, frag.strand);
		*/
    }

	fclose(fichero);
	return 0;
	
}
