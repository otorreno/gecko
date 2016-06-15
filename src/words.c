/* words.c
 new version of "words" program. Instead of working with the bin file
 this version works over the plain-text sequence (and do not uses the
 seqio library -for managing big data sets-)
 Problems have been detected in the previous version.

 Using this option makes unnecesary to "format" the sequence (moving
 to binary code). Thus, after masking the Low-Complex regions this
 program can be used as next step

 This new uses "./words seq.IN words.OUT
 where seq.IN is a plain-text sequence
 words.OUT is a binary file of "wentry" structures

 NEXT: this program load the full seq into memory. Need to be modified
 to load partial chunks of sequence (not difficult)
 -----------------------------------------------------4.Feb.2012
 ortrelles @ uma.es
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"

static int WORD_SIZE=32;
//static int BITS_PER_BASE=2;
static int BYTES_IN_WORD=8;//(int)ceil(WORD_SIZE/8.*BITS_PER_BASE);

void shift_word(word * w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}

void main_FILE(char * inFile, char * outFile){
	FILE *f;
	FILE *f2;
	char c;

	if ((f=fopen(inFile,"rt"))==NULL){
	perror("opening sequence file");
	}
	if ((f2=fopen(outFile,"wb"))==NULL) {
		terror("opening OUT sequence Words file");
	}
	
	c=fgetc(f);
	while(c!='\n'){
		c=fgetc(f);
	}

	wentry WE;
	WE.seq=0;
	unsigned long index=0;
	unsigned long inEntry=0;
	unsigned long NW=0;
	unsigned long Tot=0;
	unsigned long NoACGT=0;
	unsigned long NoC=0;
	c=fgetc(f);
	while(!feof(f)){
		if (!isupper(toupper(c))){
			if(c=='>'){
				c = fgetc(f);
				while (c != '\n')
					c = fgetc(f);
				WE.seq++;
				inEntry=0;
				index++;
			}
			NoC++;
			c=fgetc(f);
			continue;
		}
		shift_word(&WE.w);
		switch (c) {
			case 'A': inEntry++; break;
			case 'C':
				WE.w.b[BYTES_IN_WORD-1]|=1;
				inEntry++;
				break;
			case 'G':
				WE.w.b[BYTES_IN_WORD-1]|=2;
				inEntry++;
				break;
			case 'T':
				WE.w.b[BYTES_IN_WORD-1]|=3;
				inEntry++;
				break;
			default :
				inEntry=0; NoACGT++; break;
		}
		index++;
		Tot++;
		if(inEntry>=(unsigned long)WORD_SIZE){
			WE.pos=index-WORD_SIZE;
			NW++;
			fwrite(&WE,sizeof(wentry),1,f2);
		}
		c=fgetc(f);
	}
	//printf("FILE: Create %d Words --(seqLen=%d NoACGT=%d noChar=%d\n",NW,Tot,NoACGT, NoC);
	fclose(f);

}

int main(int ac, char** av){
	if(ac!=3){
		terror("USE: words seqFile.IN words.OUT");
	}
	main_FILE(av[1], av[2]);
	return 0;
}

