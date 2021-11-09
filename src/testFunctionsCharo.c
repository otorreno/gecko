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

int main(int ac,char** av){
	
	// Read fragments
	struct FragFile *f;
	int nf; // Number of fragments
	int n = 0;
	uint64_t xtotal,ytotal;
	nf = 0;
    FILE *fichero;

	f = readFragmentsv2(av[1],&nf,&xtotal,&ytotal);
    fprintf(stdout,"%" PRIu64"\n", xtotal);
	fprintf(stdout,"%" PRIu64"\n", ytotal);
	fprintf(stdout,"%" PRIu64"\n", nf);
	fichero = fopen(av[2], "wb");
	if(fichero == NULL){
		terror("Abriendo fichero de salida");
	}
	fwrite(&xtotal,sizeof(uint64_t), 1, fichero);
	fwrite(&ytotal,sizeof(uint64_t), 1, fichero);
	for(n;n < nf; n++){
		writeFragmentRaw(&f[n], fichero);
		
	}

	fclose(fichero);
	return 0;
	
}