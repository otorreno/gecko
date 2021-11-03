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

int main(int ac,char** av){
	
	// Read fragments
	struct FragFile* f;
	int nf; // Number of fragments
	uint64_t xtotal,ytotal;
	nf = 0;
    FILE *fichero;
	f = readFragmentsv2(av[1],&nf,&xtotal,&ytotal);
    fichero = fopen("ficheroCharo.txt", "rb");
    writeFragmentRaw(f, fichero);
	return 0;
}