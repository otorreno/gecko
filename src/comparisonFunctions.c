#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "quicksortHit.h"

#define MAXBUF 10000000
#define MaxREP 10000

#define BaseType hit

unsigned long scoreMax(char *seq, char *seq2, uint64_t len, int point) {
	//CHANGE WHEN USING PAM MATRIX
	//Scor=0; for (i=0;ii<len;i++) if (Seq[i]==Seq2[i])scor+=Ptos
	/* TO DELETE WARNINGS*/
	if(seq+1){
	}

	if(seq2+1){
	}
	/**********************/

	return len * point;
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

void endianessConversion(char *source, char *target, int numberOfBytes){
	int i,j;
	for(i=numberOfBytes-1;i>=0;i--){
		j=numberOfBytes-1-i;
		target[j]=source[i];
	}
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f){
	char tmpArray[8];

	if(htons(1)==1){
		//big endian
		if(fread(&frag->diag, sizeof(int64_t), 1, f)!=1){
			if(feof(f))return;
			terror("Error reading the HSP diagonal");
		}
		if(fread(&frag->xStart, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xStart");
		}
		if(fread(&frag->yStart, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yStart");
		}
		if(fread(&frag->xEnd, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xEnd");
		}
		if(fread(&frag->yEnd, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yEnd");
		}
		if(fread(&frag->length, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP length");
		}
		if(fread(&frag->ident, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP identities");
		}
		if(fread(&frag->score, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP score");
		}
		if(fread(&frag->similarity, sizeof(float), 1, f)!=1){
			terror("Error reading the HSP similarity");
		}
		if(fread(&frag->seqX, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqX");
		}
		if(fread(&frag->seqY, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqY");
		}
		if(fread(&frag->block, sizeof(int64_t), 1, f)!=1){
			terror("Error reading the HSP block");
		}
		frag->strand = fgetc(f);
	} else {
		//little endian
		if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
			if(feof(f))return;
			terror("Error reading the HSP diagonal");
		}
		endianessConversion(tmpArray, (char *)(&frag->diag), sizeof(int64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xStart");
		}
		endianessConversion(tmpArray, (char *)(&frag->xStart), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yStart");
		}
		endianessConversion(tmpArray, (char *)(&frag->yStart), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP xEnd");
		}
		endianessConversion(tmpArray, (char *)(&frag->xEnd), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP yEnd");
		}
		endianessConversion(tmpArray, (char *)(&frag->yEnd), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP length");
		}
		endianessConversion(tmpArray, (char *)(&frag->length), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP identity");
		}
		endianessConversion(tmpArray, (char *)(&frag->ident), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP score");
		}
		endianessConversion(tmpArray, (char *)(&frag->score), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(float), 1, f)!=1){
			terror("Error reading the HSP float");
		}
		endianessConversion(tmpArray, (char *)(&frag->similarity), sizeof(float)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqX");
		}
		endianessConversion(tmpArray, (char *)(&frag->seqX), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading the HSP seqY");
		}
		endianessConversion(tmpArray, (char *)(&frag->seqY), sizeof(uint64_t)); 
		if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
			terror("Error reading the HSP block");
		}
		endianessConversion(tmpArray, (char *)(&frag->block), sizeof(int64_t)); 
		frag->strand = fgetc(f);
	}
}

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//Big endian
		fwrite(&frag->diag, sizeof(int64_t), 1, f);
		fwrite(&frag->xStart, sizeof(uint64_t), 1, f);
		fwrite(&frag->yStart, sizeof(uint64_t), 1, f);
		fwrite(&frag->xEnd, sizeof(uint64_t), 1, f);
		fwrite(&frag->yEnd, sizeof(uint64_t), 1, f);
		fwrite(&frag->length, sizeof(uint64_t), 1, f);
		fwrite(&frag->ident, sizeof(uint64_t), 1, f);
		fwrite(&frag->score, sizeof(uint64_t), 1, f);
		fwrite(&frag->similarity, sizeof(float), 1, f);
		fwrite(&frag->seqX, sizeof(uint64_t), 1, f);
		fwrite(&frag->seqY, sizeof(uint64_t), 1, f);
		fwrite(&frag->block, sizeof(int64_t), 1, f);
		fputc(frag->strand, f);
	} else {
		//Little endian
		endianessConversion((char *)(&frag->diag), tmpArray, sizeof(int64_t));
		fwrite(tmpArray, sizeof(int64_t), 1, f);
		endianessConversion((char *)(&frag->xStart), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->yStart), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->xEnd), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->yEnd), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->length), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->ident), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->score), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->similarity), tmpArray, sizeof(float));
		fwrite(tmpArray, sizeof(float), 1, f);
		endianessConversion((char *)(&frag->seqX), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->seqY), tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
		endianessConversion((char *)(&frag->block), tmpArray, sizeof(int64_t));
		fwrite(tmpArray, sizeof(int64_t), 1, f);
		fputc(frag->strand, f);
	}
}

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		if(fread(length, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading sequence length");
		}
	} else {
		//little endian
		if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
			terror("Error reading sequence length");
		}
		endianessConversion(tmpArray, (char *)length, sizeof(uint64_t));
	}
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f){
	char tmpArray[8];
	if(htons(1)==1){
		//big endian
		fwrite(length, sizeof(uint64_t), 1, f);
	} else {
		//little endian
		endianessConversion((char *)length, tmpArray, sizeof(uint64_t));
		fwrite(tmpArray, sizeof(uint64_t), 1, f);
	}
}

/**
 * Function to return the sizeof a fragment.
 * Due to architeture issues the value of sizeof built-in
 * function could be different
 */
long int sizeofFragment(){
	return 9*sizeof(uint64_t)+2*sizeof(int64_t)+1*sizeof(float)+1*sizeof(char);
}

/**************** ARJONA *******************/
void cpyFrag2(struct FragFile *f, struct FragFile g){

	f->diag = g.diag;
	f->xStart = g.xStart;
	f->xEnd = g.xEnd;
	f->yStart = g.yStart;
	f->yEnd = g.yEnd;
	f->length = g.length;
	f->score = g.score;
	f->similarity = g.similarity;
	f->seqX = g.seqX;
	f->seqY = g.seqY;
	f->ident = g.ident;
	f->block = g.block;
	f->strand = g.strand;
}

/******/
struct FragFile* readFragments(char* s,int* nf,uint64_t *xtotal,uint64_t *ytotal){
//Fragment* readFragments(char* s,int* nf,int *xtotal,int *ytotal){

	FILE* fe;

	struct FragFile* fs,f;
	int n;

	if((fe=fopen(s,"rb"))==NULL){
		printf("***ERROR Opening input file");
		exit(-1);
	}
	n=0;
	readSequenceLength(xtotal, fe);
	readSequenceLength(ytotal, fe);
//	long int longFile;
//	fseek(fe,0,SEEK_END);
//	longFile=ftell(fe);
//	n=(int)(longFile-2*sizeof(uint64_t))/sizeof(struct FragFile);
//	printf("\n ReadFragments Complete\nnum: %d\n",n);

	// Alternativa +++++++++++++++++
	n=0;
	readFragment(&f, fe);
	while(!feof(fe)){
		readFragment(&f, fe);
		n++;
	}
//	printf("\n ReadFragments Complete\nnum: %d\n",n);
	//+++++++++++++
	rewind(fe);
	fs=(struct FragFile*)malloc(sizeof(struct FragFile)*(n));
	if(fs==NULL){printf("****ERROR: Out of memory\n");exit(-1);}

	*nf=n;
	readSequenceLength(xtotal, fe);
	readSequenceLength(ytotal, fe);
	n=0;
	readFragment(&f, fe);
	while(!feof(fe)){



		if(f.length>0){
			cpyFrag2(&fs[n],f);
		}
		n++;
		readFragment(&f, fe);
		//fprintf(stdout,"%d\t%" PRId64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%c" "\t%" PRIu64 "\n",n,fs[n].xStart, fs[n].yStart, fs[n].xEnd, fs[n].yEnd, fs[n].length, fs[n].strand, fs[n].ident);

	}
	*nf=n;


	fclose(fe);
	return fs;
}

int differentSequences(hit h1, hit h2){
    return h1.seqX != h2.seqX || h1.seqY != h2.seqY;
}

hit *hits(hashentry *entriesX, uint64_t nEntriesX, hashentry *entriesY, uint64_t nEntriesY, int wSize, uint64_t *HIB){
    int bufSize=MAXBUF;
    uint64_t hitsInBuf = 0;
    uint64_t i, j, k=0, l=0;
    int comp;
    int firstMatch = 0, endMatch = 0;
    uint64_t nHits = 0, wordMatches = 0;
    int stepX, stepY;
    long offsetWordBefore, offsetWordAfter;
    hit *hBuf;

    if ((hBuf = (hit*) calloc(sizeof(hit), MAXBUF)) == NULL)
        terror("HITS: memory for I-O buffer");

//	fprintf(stdout, "Memoria reservada para el buffer de hits\n");

    while (k < nEntriesX && l < nEntriesY){
		comp = wordcmp(&entriesX[k].w.b[0], &entriesY[l].w.b[0], wSize);
//		fprintf(stdout, "k: %" PRIu64 " l:%" PRIu64 "\n", k, l);
		if (comp < 0) {
			k++;
			//Save position of first missmatch after matches and rewind
			if (firstMatch) {
				offsetWordAfter = l-1;
				l = offsetWordBefore;
				firstMatch = 0;
				endMatch = 1;
			}
			continue;
		}
		if (comp > 0) {
			//No more matches, go to the next word
			if (endMatch) {
				l = offsetWordAfter;
				endMatch = 0;
			}
            l++;
			continue;
		}

		wordMatches++;
//		fprintf(stdout, "Hay match\n");

		// Saving the offset of the first match
		if (wSize < 32 && !firstMatch) {
			offsetWordBefore = l-1;
			firstMatch = 1;
		}

		// Hits-----------------------
		if (entriesX[k].num > MaxREP)
			stepX = entriesX[k].num / MaxREP;
		else
			stepX = 1;
		if (entriesY[l].num > MaxREP)
			stepY = entriesY[l].num / MaxREP;
		else
			stepY = 1;

//		fprintf(stdout, "antes del for\n");

		for (i = 0; i < entriesX[k].num; i += stepX)
			for (j = 0; j < entriesY[l].num; j += stepY) {
                hBuf[hitsInBuf].diag = entriesX[k].locs[i].pos - entriesY[l].locs[j].pos;
                hBuf[hitsInBuf].posX = entriesX[k].locs[i].pos;
                hBuf[hitsInBuf].seqX = entriesX[k].locs[i].seq;
                hBuf[hitsInBuf].posY = entriesY[l].locs[j].pos;
                hBuf[hitsInBuf].seqY = entriesY[l].locs[j].seq;

                hitsInBuf++;
				if (hitsInBuf >= bufSize ) {
//                    fprintf(stdout, "reallocating\n");
                    hBuf = realloc(hBuf,(hitsInBuf+MAXBUF)*sizeof(hit));
					bufSize += MAXBUF;
				}
			}

//        fprintf(stdout, "despues del for\n");

        nHits += ((entriesX[k].num / stepX) * (entriesY[l].num / stepY));

		if (!firstMatch)k++;
		l++;
	}

    hBuf = realloc(hBuf,hitsInBuf*sizeof(hit));
    fprintf(stdout, "hitsInBuf: %" PRIu64 "\n", hitsInBuf);

    //SortHits
    psortH(32,hBuf,hitsInBuf);

    //FilterHits
    hit *hBuf2 = (hit *)calloc(hitsInBuf, sizeof(hit));
    if(hBuf2==NULL)
        perror("Not enough memory for filtered hits buffer");

    int64_t diagonal;
    uint64_t lastPosition;
    uint64_t originalNumberOfHits = 0, finalNumberOfHits = 0;

    lastPosition = 0;

    originalNumberOfHits = hitsInBuf;
    i=0;
    diagonal = hBuf[i].diag;
    while (i < (hitsInBuf - 1)) {
        if(differentSequences(hBuf[i], hBuf[i+1])){
            lastPosition=0;
            diagonal = hBuf[i].diag;
            i++;
            continue;
        }

        if (diagonal != hBuf[i+1].diag || hBuf[i+1].posX > lastPosition) {
            memcpy(&hBuf2[finalNumberOfHits++],&hBuf[i], sizeof(hit));
            lastPosition = hBuf[i].posX + (2 * wSize - 1);
            diagonal = hBuf[i].diag;
        }
        i++;
    }

    if(diagonal != hBuf[i].diag || hBuf[i].posX > (lastPosition - wSize)){
        memcpy(&hBuf2[finalNumberOfHits++],&hBuf[i], sizeof(hit));
    }

    hBuf2 = realloc(hBuf2,finalNumberOfHits*sizeof(hit));
    free(hBuf);

    fprintf(stdout,
            "\nfilterHits\noriginal number of Hits=%" PRIu64 "  Final number of hits=%" PRIu64 "\n",
            originalNumberOfHits, finalNumberOfHits);
    //End of filterhits

    *HIB = finalNumberOfHits;

	return hBuf2;
}

int GTH(hit a1, hit a2) {
    if(a1.diag > a2.diag)
        return 1;
    else if (a1.diag < a2.diag)
        return 0;
    if (a1.posX > a2.posX)
        return 1;

    return 0;
}
/************************/

