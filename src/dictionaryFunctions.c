#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "quicksortWord.h"

int seq2word(char* buf, int wsize, word* w) {
	int i;
	int b = 6;
	memset(w, 0, sizeof(word));

	for (i = 0; i < wsize; i++) {
		if (buf[i] >= 4)
			return -1;
		w->b[i / 4] |= buf[i] << b;
		b -= 2;
		if (b < 0)
			b = 6;
	}
	return 0;
}

void skipIDLine(FILE *fIn) {
	char c;
	// first line (skip ID fasta Line)
	c = fgetc(fIn);
	while (c != '\n')
		c = fgetc(fIn);
}

int letterToIndex(char c) {
	// coding (a=0,c=1,g=2,t=3,'>'=4 others=9 )
	switch (c) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	case '>':
		return 4;
	default:
		return 9;
	}
}

int loadSequence(char *fileName, char *seq, uint64_t *Tot) {
	FILE *fIn;
	char c;

	//Opening input and output files
	if ((fIn = fopen(fileName, "rt")) == NULL)
		terror("opening sequence file");

	//Skip the identification of the sequence
	skipIDLine(fIn);

	// Load Sequence into memory
	c = fgetc(fIn);
	while (!feof(fIn)) {
		//Check if is a letter
		if (!isupper(c)) {
			/*
			 * If not is a start of a sequence,
			 * then read a new char and continue
			 */
			if (c != '>') {
				c = fgetc(fIn);
				continue;
			}
		}

		//Get the index of the letter
		seq[*Tot] = letterToIndex(c);

		//Check if is a multi-sequence file
		if (c == '>') {
			skipIDLine(fIn);
		}

		(*Tot)++;
		c = fgetc(fIn);
	}

	fclose(fIn);
	return 0;
}

int wordcmp(unsigned char *w1, unsigned char *w2, int n) {

	int i = 0, limit;

	if(n%4 != 0){
		w1[n/4] = w1[n/4] >> (2*(3-((n-1)%4)));
		w1[n/4] = w1[n/4] << (2*(3-((n-1)%4)));
		w2[n/4] = w2[n/4] >> (2*(3-((n-1)%4)));
		w2[n/4] = w2[n/4] << (2*(3-((n-1)%4)));
		limit=(n/4)+1;
	} else {
		limit = n/4;
	}

	for (i=0;i<limit;i++) {
		if (w1[i]<w2[i]) return -1;
		if (w1[i]>w2[i]) return +1;
	}
	return 0;
}

void showWord(word* w, char *ws) {
	char Alf[] = { 'A', 'C', 'G', 'T' };
	int i;
	int wsize = 8;
	unsigned char c;
	for (i = 0; i < wsize; i++) {
		c = w->b[i];
		c = c >> 6;
		ws[4*i] = Alf[(int) c];
		c = w->b[i];
		c = c << 2;
		c = c >> 6;
		ws[4*i+1] = Alf[(int) c];
		c = w->b[i];
		c = c << 4;
		c = c >> 6;
		ws[4*i+2] = Alf[(int) c];
		c = w->b[i];
		c = c << 6;
		c = c >> 6;
		ws[4*i+3] = Alf[(int) c];
	}
}

int GTW(wentry a1, wentry a2){
	if(a1.w.b[0] > a2.w.b[0]) return 1;
	else if(a1.w.b[0] < a2.w.b[0]) return 0;

	if(a1.w.b[1] > a2.w.b[1]) return 1;
	else if(a1.w.b[1] < a2.w.b[1]) return 0;

	if(a1.w.b[2] > a2.w.b[2]) return 1;
	else if(a1.w.b[2] < a2.w.b[2]) return 0;

	if(a1.w.b[3] > a2.w.b[3]) return 1;
	else if(a1.w.b[3] < a2.w.b[3]) return 0;

	if(a1.w.b[4] > a2.w.b[4]) return 1;
	else if(a1.w.b[4] < a2.w.b[4]) return 0;

	if(a1.w.b[5] > a2.w.b[5]) return 1;
	else if(a1.w.b[5] < a2.w.b[5]) return 0;

	if(a1.w.b[6] > a2.w.b[6]) return 1;
	else if(a1.w.b[6] < a2.w.b[6]) return 0;

	if(a1.w.b[7] > a2.w.b[7]) return 1;
	else if(a1.w.b[7] < a2.w.b[7]) return 0;

	if(a1.seq > a2.seq) return 1;
	else if(a1.seq < a2.seq) return 0;

	if(a1.pos > a2.pos) return 1;
	return 0;
}

void shift_word(word * w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}

void *dictionary(void *a){
	DictionaryArgs *args=(DictionaryArgs*)a;
	FILE *f;
	char c;
	uint64_t SIZE = 0;
	wentry *words = NULL;

	if ((f=fopen(args->seqFile,"rt"))==NULL){
		fprintf(stdout, "opening sequence file: %s\n", args->seqFile);
		perror("opening sequence file");
	}

	fseek(f,0,SEEK_END);
	SIZE=ftell(f);
	fseek(f,0,SEEK_SET);

	if((words = calloc(SIZE, sizeof(wentry)))==NULL){
		perror("not enough memory for words array");
	}

//    fprintf(stdout, "Memoria reservada: %" PRIu64 "\n", SIZE);

	c=fgetc(f);
	while(c!='\n'){
		c=fgetc(f);
	}

	wentry WE;
	WE.seq=0;
	unsigned long index=0;
	unsigned long inEntry=0;
	uint64_t NW=0;
	unsigned long Tot=0;
	unsigned long NoACGT=0;
	unsigned long NoC=0;
	unsigned long nLocs = 0;
	unsigned long loc_size = 0;
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
			memcpy(&words[NW],&WE,sizeof(wentry));
		}
		c=fgetc(f);
	}
	//printf("FILE: Create %d Words --(seqLen=%d NoACGT=%d noChar=%d\n",NW,Tot,NoACGT, NoC);
	fclose(f);

	words = (wentry *)realloc(words,NW*sizeof(wentry));

//    fprintf(stdout, "Antes del sort\n");

	psortW(32,words,NW);

//    fprintf(stdout, "Despues del sort\n");

//    fprintf(stdout, "Antes del w2hd\n");
	hashentry* entries = NULL;
	if((entries = calloc(NW, sizeof(hashentry)))==NULL){
		perror("not enough memory for hashentry array");
	}

	memcpy(&entries[0].w.b[0],&words[0].w.b[0],8);
	entries[0].num=0;
	entries[0].locs = NULL;
	if((entries[0].locs = calloc(SIZE_LOC, sizeof(location)))==NULL){
		perror("not enough memory for locs array");
	}
	loc_size=SIZE_LOC;

//    fprintf(stdout, "memoria reservada w2hd\n");

	uint64_t i=0;
	uint64_t j=0;
	location loc;
	while (i<NW){
		loc.pos=words[i].pos;
		loc.seq=words[i].seq;
		if (wordcmp(&entries[j].w.b[0],&words[i].w.b[0],32)!=0) {
			entries[j].locs = realloc(entries[j].locs,nLocs);
			j++;
			memcpy(&entries[j].w.b[0],&words[i].w.b[0],8);
			entries[j].num=0;
			nLocs=0;
			if((entries[j].locs = calloc(SIZE_LOC, sizeof(location)))==NULL){
				perror("not enough memory for locs array");
			}
			loc_size=SIZE_LOC;
		}

		if(nLocs >= loc_size){
			loc_size += SIZE_LOC;
//            fprintf(stdout, "Reallocating memory from: %lu to: %lu\n",loc_size-SIZE_LOC,loc_size);
			entries[j].locs = realloc(entries[j].locs,loc_size*sizeof(location));
			if(entries[j].locs == NULL){
				perror("Error re-allocation location array");
			}

		}
		memcpy(&entries[j].locs[nLocs++],&loc,sizeof(location));
		entries[j].num++;
		i++;
	}
//    fprintf(stdout, "Despues del w2hd\n");

//    char wordString[33];
//	wordString[32] = '\0';
//    uint64_t nEntries = j;
//	for(i=0;i<nEntries;i++){
//		showWord(&entries[i].w, wordString);
//		fprintf(stdout, "Words(%" PRIu64 "):%s Repetitions: %" PRIu64 "Positions: ", i, wordString, entries[i].num);
//        for(j=0;j<entries[i].num;j++){
//            fprintf(stdout,"(%" PRIu64 ",%" PRIu64 ") ",entries[i].locs[j].pos,entries[i].locs[j].seq);
//        }
//        fprintf(stdout, "\n");
//	}

	free(words);
	*(args->nEntries)=j;
	return entries;
}
