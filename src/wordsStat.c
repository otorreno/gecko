/* wordsStat.c 

 This program shows the position of each word and the
 sequence it corresponds to
 -----------------------------------------------------12.Nov.2012
 oscart @ uma.es
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"

int main(int ac, char** av) {
	FILE *f;
	uint64_t i = 0;
	hashentryR wordEntry;
	char wordString[33];
	wordString[32] = '\0';

	if (ac != 2)
		terror("USE: wordsStat words");

	if ((f = fopen(av[1], "rt")) == NULL)
		terror("opening words file");

	if(fread(&wordEntry, sizeof(hashentryR), 1, f)!=1){
		terror("Empty words file");
	}
	while (!feof(f)) {
		showWord(&wordEntry.w, wordString);
		fprintf(stdout, "%s\n", wordString);
		if(fread(&wordEntry, sizeof(hashentryR), 1, f)!=1){
			if(ferror(f))terror("Error reading words file");
		}
		i++;
	}
	showWord(&wordEntry.w, wordString);
        fprintf(stdout, "%s\n", wordString);


	fclose(f);

	return 0;
}

