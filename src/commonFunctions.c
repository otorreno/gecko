#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include "structs.h"

void terror(char *s) {
	printf("ERR**** %s ****\n", s);
	exit(-1);
}

unsigned long timestart(){
	struct timeval tv;

	gettimeofday(&tv,NULL);

	return (tv.tv_usec/1000) + (tv.tv_sec*1000);
}

unsigned long timestop(unsigned long start){
	struct timeval tv;

	gettimeofday(&tv,NULL);

	return (tv.tv_usec/1000) + (tv.tv_sec*1000) - start;
}

struct Sequence* LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nStruct, uint64_t *nSeqs,
						  int fAst) {
	char c;
	uint64_t lon = 0, k = 0, ns, seqs = 0;
	uint64_t lonFinal = 0;
	struct Sequence *sX, *sX2; //sX will be the first elem. sX2 will generate all the structure

	//Initialize
	*n = 0;
	*nStruct = 0;

	//Memory
	ns = 1;
	if ((sX = (struct Sequence*) malloc(sizeof(struct Sequence))) == NULL)
		terror("Memory...");

	if (!fAst)
		while ((c = getc(f)) != '>' && !feof(f))
			; //start seq
	if (feof(f))
		return 0;
	while ((c = getc(f)) == ' ')
		;

	while (k < MAXLID && c != '\n' && c != ' ') {
		if (feof(f))
			return 0;

		sX->ident[k++] = c;
		c = getc(f);
	}

	sX->ident[k] = 0; //end of data.
	while (c != '\n')
		c = getc(f);
	c = getc(f);

	//start list with sX2
	sX2 = sX;
	while (/*c!='*'&&*/!feof(f)) {
		c = toupper(c);
		if (c == '>') {
			fAst = 1;
			seqs++;
			sX2->datos[lon++] = '*';
			while (c != '\n') {
				if (feof(f))
					return 0;
				c = getc(f);
			}
			//break;
		}
		if (isupper(c))
			sX2->datos[lon++] = c;
		if (c == '*') {
			sX2->datos[lon++] = c;
		}
		c = getc(f);

		//Check if the length is the end of this struct
		if (lon >= MAXLS) {
			lonFinal += lon;
			lon = 0;
			ns++;
			if ((sX = (struct Sequence*) realloc(sX,
												 ns * sizeof(struct Sequence))) == NULL)
				terror("Memory...");
			sX2 = sX + ns - 1;
		}
	}

	if (lon < MAXLS)
		sX2->datos[lon] = 0x00;

	lonFinal += lon;
	*nStruct = ns;
	*nSeqs = seqs + 1;
	*n = lonFinal - seqs;
	return sX;
}

char getValue(struct Sequence *s, uint64_t pos, int ns) {
	struct Sequence *aux = s;
	int nActual = 1;

	while (pos >= MAXLS) {
		aux++;
		pos -= MAXLS;
		nActual++;
		if (nActual > ns){
			terror("out of sequence.");
		}
	}

	return aux->datos[pos];
}

long getSeqLength(struct Sequence *s, uint64_t start, int ns) {
	int nActual = 1;
	struct Sequence *aux = s;
	while (start >= MAXLS) {
		aux++;
		start -= MAXLS;
		nActual++;
		if (nActual > ns)
			terror("out of sequence.");
	}
	uint64_t s1 = start;
	while (s1 > 0 && aux->datos[s1] != '*') {
		s1--;
	}
	s1++;
	char *tmp = strchr(aux->datos + s1, '*');
	if (tmp == NULL) {
		return strlen(aux->datos) - s1 + 1;
	}
	return tmp - (aux->datos + s1) + 1;
}
