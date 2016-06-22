#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include <pthread.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "quicksortHit.h"
#include "comparisonFunctions.h"

#define MAXBUF 10000000
#define MaxREP 10000
#define POINT 4

int FragFromHit(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
                uint64_t n0, uint64_t nsx, struct Sequence *sY,
                uint64_t n1, uint64_t nsy, uint64_t nSeqs1, uint64_t LenOrPer, uint64_t SimTh, int WL);
int FragFromHitReverse(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
					   uint64_t n0, uint64_t nsx, struct Sequence *sY,
					   uint64_t n1, uint64_t nsy, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL);
struct FragFile *frags(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL, uint64_t *nF, uint64_t *nHU);
struct FragFile *fragsReverse(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL, uint64_t *nF, uint64_t *nHU);

inline char complement(char c){
    switch(c){
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return c;
    }
}

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

uint64_t filterHits(hit *hBuf, uint64_t hitsInBuf, int wSize, hit **output){
    hit *hBuf2 = (hit *)calloc(hitsInBuf, sizeof(hit));
    if(hBuf2==NULL)
        perror("Not enough memory for filtered hits buffer");

    int64_t diagonal;
    uint64_t lastPosition;
    uint64_t originalNumberOfHits = 0, finalNumberOfHits = 0;

    lastPosition = 0;

    originalNumberOfHits = hitsInBuf;
    uint64_t i=0;
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

    *output = hBuf2 = realloc(hBuf2,finalNumberOfHits*sizeof(hit));
    if(hBuf2 == NULL)
        perror("Error reallocating filtered hits array");

    fprintf(stdout,
            "\nfilterHits\noriginal number of Hits=%" PRIu64 "  Final number of hits=%" PRIu64 "\n",
            originalNumberOfHits, finalNumberOfHits);

    return finalNumberOfHits;
    //End of filterhits
}

void *sortHitsFilterHitsFragHitsTh(void *a){
    ComparisonArgs *args = (ComparisonArgs *)a;
    uint64_t HIB;
    hit *hBuf2;

    if(args->nHits > 0) {
        psortH(32, args->hits, args->nHits);

        HIB = filterHits(args->hits, args->nHits, args->wSize, &hBuf2);
        free(args->hits);

        struct FragFile *fragsBuf = frags(args->seqX, args->seqY, hBuf2, HIB, args->Lmin, args->SimTh,
                                          args->wSize, args->nFrags, args->nHitsUsed);
        free(hBuf2);
        return fragsBuf;
    }

    return NULL;
}

void *sortHitsFilterHitsFragHitsReverseTh(void *a){
	ComparisonArgs *args = (ComparisonArgs *)a;
	uint64_t HIB;
	hit *hBuf2;

	if(args->nHits > 0) {
		psortH(32, args->hits, args->nHits);

		HIB = filterHits(args->hits, args->nHits, args->wSize, &hBuf2);
		free(args->hits);

		struct FragFile *fragsBuf = fragsReverse(args->seqX, args->seqY, hBuf2, HIB, args->Lmin, args->SimTh,
										  args->wSize, args->nFrags, args->nHitsUsed);
		free(hBuf2);
		return fragsBuf;
	}

	return NULL;
}

struct FragFile *hitsAndFrags(char *seqX, char *seqY, char *out, uint64_t seqXLen, uint64_t seqYLen,
                              hashentry *entriesX, uint64_t nEntriesX, hashentry *entriesY, uint64_t nEntriesY,
                              int wSize, uint64_t Lmin, uint64_t SimTh, uint64_t *nFrags){
    int bufSizeForward=MAXBUF, bufSizeReverse=MAXBUF;
    uint64_t hitsInBufForward = 0, hitsInBufReverse = 0, hitsInBufForwardUsed = 0, hitsInBufReverseUsed = 0,
            nFragsForward, nFragsReverse;
    uint64_t i, j, k=0, l=0;
    int comp;
    int firstMatch = 0, endMatch = 0;
    uint64_t nHits = 0, wordMatches = 0;
    int stepX, stepY;
    long offsetWordBefore = 0, offsetWordAfter = 0;
    hit *hBufForward, *hBufReverse;
    uint64_t minSeqXLenSeqYLen = min(seqXLen, seqYLen);
    FILE *fOut, *fInf;
    char infoFileName[1024];

    if ((hBufForward = (hit*) calloc(MAXBUF, sizeof(hit))) == NULL)
        terror("HITS: memory for I-O buffer");

    if ((hBufReverse = (hit*) calloc(MAXBUF, sizeof(hit))) == NULL)
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
				l = offsetWordBefore + 1;
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
				if(entriesY[l].locs[j].strand == 'f') {
                    hBufForward[hitsInBufForward].diag = entriesX[k].locs[i].pos - entriesY[l].locs[j].pos;
                    hBufForward[hitsInBufForward].posX = entriesX[k].locs[i].pos;
                    hBufForward[hitsInBufForward].seqX = entriesX[k].locs[i].seq;
                    hBufForward[hitsInBufForward].posY = entriesY[l].locs[j].pos;
                    hBufForward[hitsInBufForward].seqY = entriesY[l].locs[j].seq;

                    hitsInBufForward++;
                    if (hitsInBufForward >= bufSizeForward) {
//                    fprintf(stdout, "reallocating\n");
                        hBufForward = realloc(hBufForward, (hitsInBufForward + MAXBUF) * sizeof(hit));
                        if(hBufForward == NULL)
                            perror("Error reallocating forward hits buffer");
                        bufSizeForward += MAXBUF;
                    }
                } else {
                    hBufReverse[hitsInBufReverse].diag = entriesX[k].locs[i].pos + entriesY[l].locs[j].pos - minSeqXLenSeqYLen;
                    hBufReverse[hitsInBufReverse].posX = entriesX[k].locs[i].pos;
                    hBufReverse[hitsInBufReverse].seqX = entriesX[k].locs[i].seq;
                    hBufReverse[hitsInBufReverse].posY = entriesY[l].locs[j].pos;
                    hBufReverse[hitsInBufReverse].seqY = entriesY[l].locs[j].seq;

                    hitsInBufReverse++;
                    if (hitsInBufReverse >= bufSizeReverse) {
//                    fprintf(stdout, "reallocating\n");
                        hBufReverse = realloc(hBufReverse, (hitsInBufReverse + MAXBUF) * sizeof(hit));
                        if(hBufReverse == NULL)
                            perror("Error reallocating reverse complement hits buffer");
                        bufSizeReverse += MAXBUF;
                    }
                }
			}

//        fprintf(stdout, "despues del for\n");

        nHits += ((entriesX[k].num / stepX) * (entriesY[l].num / stepY));

		if (!firstMatch)k++;
		l++;
	}

    hBufForward = realloc(hBufForward,hitsInBufForward*sizeof(hit));
    if(hBufForward == NULL)
        perror("Error reallocating forward hits buffer");
    fprintf(stdout, "hitsInBufForward: %" PRIu64 "\n", hitsInBufForward);

    hBufReverse = realloc(hBufReverse,hitsInBufReverse*sizeof(hit));
    if(hBufReverse == NULL)
        perror("Error reallocating reverse complement hits buffer");
    fprintf(stdout, "hitsInBufReverse: %" PRIu64 "\n", hitsInBufReverse);

    ComparisonArgs argsForward, argsReverse;
    pthread_t thF, thR;

    argsForward.seqX = seqX;
    argsForward.seqY = seqY;
    argsForward.nHits = hitsInBufForward;
    argsForward.nHitsUsed = &hitsInBufForwardUsed;
    argsForward.hits = hBufForward;
    argsForward.nFrags = &nFragsForward;
    argsForward.Lmin = Lmin;
    argsForward.SimTh = SimTh;
    argsForward.wSize = wSize;

    argsReverse.seqX = seqX;
    argsReverse.seqY = seqY;
    argsReverse.nHits = hitsInBufReverse;
    argsReverse.nHitsUsed = &hitsInBufReverseUsed;
    argsReverse.hits = hBufReverse;
    argsReverse.nFrags = &nFragsReverse;
    argsReverse.Lmin = Lmin;
    argsReverse.SimTh = SimTh;
    argsReverse.wSize = wSize;

    pthread_create(&thF,NULL,sortHitsFilterHitsFragHitsTh,(void*)(&argsForward));
    pthread_create(&thR,NULL,sortHitsFilterHitsFragHitsReverseTh,(void*)(&argsReverse));

    struct FragFile *fragsBufForward;
    struct FragFile *fragsBufReverse;
    pthread_join(thF,(void **)&fragsBufForward);
    pthread_join(thR,(void **)&fragsBufReverse);

	fprintf(stdout, "Frags in forward: %" PRIu64 "; Frags in reverse: %" PRIu64 "\n", nFragsForward, nFragsReverse);

    *nFrags=nFragsForward + nFragsReverse;

    fprintf(stdout, "Writing fragments to disk\n");

    if ((fOut = fopen(out, "wb")) == NULL)
        terror("opening fragments output file");

    fprintf(stdout, "SeqXLen: %" PRIu64 " SeqYLen: %" PRIu64 "\n", seqXLen, seqYLen);

    writeSequenceLength(&seqXLen, fOut);
    writeSequenceLength(&seqYLen, fOut);

    for(i=0;i<nFragsForward;i++){
        writeFragment(&fragsBufForward[i],fOut);
    }
    for(i=0;i<nFragsReverse;i++){
        uint64_t tmp;
        tmp = seqYLen - fragsBufReverse[i].yEnd;
        fragsBufReverse[i].yStart = seqYLen - fragsBufReverse[i].yStart;
        fragsBufReverse[i].yEnd = tmp;
        writeFragment(&fragsBufReverse[i],fOut);
    }

    fclose(fOut);

    // metadata info
    sprintf(infoFileName, "%s.INF", out);
    if ((fInf = fopen(infoFileName, "wt")) == NULL)
        terror("opening INFO file");

    fprintf(fInf, "GECKO in memory\n");
    fprintf(fInf, "Jun. 2016 -- <oscart@uma.es>\n");
    fprintf(fInf, "SeqX filename        : %s\n", seqX);
    fprintf(fInf, "SeqY filename        : %s\n", seqY);
    fprintf(fInf, "SeqX name            : %s\n", seqX);
    fprintf(fInf, "SeqY name            : %s\n", seqY);
    fprintf(fInf, "SeqX length          : %" PRIu64 "\n", seqXLen);
    fprintf(fInf, "SeqY length          : %" PRIu64 "\n", seqYLen);
    fprintf(fInf, "Min.fragment.length  : %" PRIu64 "\n", Lmin);
    fprintf(fInf, "Min.Identity         : %" PRIu64 "\n", SimTh);
    fprintf(fInf, "Tot Hits (seeds)     : %" PRIu64 "\n", hitsInBufForward+hitsInBufReverse);
    fprintf(fInf, "Tot Hits (seeds) used: %" PRIu64 "\n", hitsInBufForwardUsed+hitsInBufReverseUsed);
    fprintf(fInf, "Total fragments      : %" PRIu64 "\n", *nFrags);
    fprintf(fInf, "========================================================\n");
    fclose(fInf);

    fprintf(stdout, "End of hits and frags function\n");

    //TODO change if we want to have the fragments outside this function
	return NULL;
}

struct FragFile *frags(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL, uint64_t *nF, uint64_t *nHU) {
    struct Sequence *sX, *sY;
    uint64_t n0, n1, ns0, ns1, nSeqs0, nSeqs1;
    uint64_t i, j;
    int newFrag;
    uint64_t nFrags = 0, nHitsUsed = 0;
    int64_t lastOffset, lastDiagonal;
    struct FragFile myF;
    FILE *f;

    //MAT
    long coverage[1000][100];
    for(i=0; i<1000; i++){
        for(j=0; j<100; j++){
            coverage[i][j]=0;
        }
    }
    //---

    struct FragFile *fragsBuf;
    fragsBuf = (struct FragFile *)calloc(nHits, sizeof(struct FragFile));
    if(fragsBuf == NULL)
        terror("Not enough memory for frags array");

    //Open files
    if ((f = fopen(seqX, "rt")) == NULL)
        terror("opening seqX file");
    sX = LeeSeqDB(f, &n0, &ns0, &nSeqs0, 0);
    fclose(f);

    if ((f = fopen(seqY, "rt")) == NULL)
        terror("opening seqY file");
    sY = LeeSeqDB(f, &n1, &ns1, &nSeqs1, 0);
    fclose(f);

    n0 += nSeqs0 - 1;
    n1 += nSeqs1 - 1;

    // read Hits
    if(nHits == 0){
        terror("Empty hits file");
    }
    lastDiagonal = hits[0].diag;
    lastOffset = hits[0].posX - 1;

    i=0;
    while (i < nHits) {
        if (lastDiagonal != hits[i].diag) {
            //Different diagonal update the variables
            lastDiagonal = hits[i].diag;
            lastOffset = hits[i].posX - 1;
        }
        //We have a casting here because of a funny error
        //where the program was saying that -1 > 0
        if (lastOffset > ((int64_t) hits[i].posX)) {
            //The hit is covered by a previous frag in the same diagonal
            i++;
            continue;
        }

        nHitsUsed++;
        newFrag = FragFromHit(coverage, &myF, &hits[i], sX, n0, ns0, sY, n1, ns1, nSeqs1, Lmin, SimTh,
                              WL);
        if (newFrag) {
            memcpy(&fragsBuf[nFrags],&myF,sizeof(struct FragFile));
            lastOffset = hits[i].posX + myF.length;
            nFrags++;
        }
        i++;
    }

    *nF = nFrags;
    *nHU = nHitsUsed;

    fprintf(stdout, "Forward strand frags calculated\n");

    fragsBuf = realloc(fragsBuf,nFrags * sizeof(struct FragFile));
    if(fragsBuf == NULL)
        perror("Error reallocating memory of forward frags array");

    return fragsBuf;
}

struct FragFile *fragsReverse(char *seqX, char *seqY, hit *hits, uint64_t nHits, uint64_t Lmin, uint64_t SimTh, int WL, uint64_t *nF, uint64_t *nHU) {
    struct Sequence *sX, *sY;
	uint64_t n0, n1, ns0, ns1, nSeqs0, nSeqs1;
	uint64_t i, j;
	int newFrag;
	uint64_t nFrags = 0, nHitsUsed = 0;
	int64_t lastOffset, lastDiagonal;
	struct FragFile myF;
	FILE *f;

	//MAT
	long coverage[1000][100];
	for(i=0; i<1000; i++){
		for(j=0; j<100; j++){
			coverage[i][j]=0;
		}
	}
	//---

	struct FragFile *fragsBuf;
	fragsBuf = (struct FragFile *)calloc(nHits, sizeof(struct FragFile));
	if(fragsBuf == NULL)
		terror("Not enough memory for frags array");

	//Open files
	if ((f = fopen(seqX, "rt")) == NULL)
		terror("opening seqX file");
	sX = LeeSeqDB(f, &n0, &ns0, &nSeqs0, 0);
	fclose(f);

	if ((f = fopen(seqY, "rt")) == NULL)
		terror("opening seqY file");
	sY = LeeSeqDB(f, &n1, &ns1, &nSeqs1, 0);
	fclose(f);

	n0 += nSeqs0 - 1;
	n1 += nSeqs1 - 1;

	// read Hits
	if(nHits == 0){
		terror("Empty hits file");
	}
	lastDiagonal = hits[0].diag;
	lastOffset = hits[0].posX - 1;

    i=0;
	while (i < nHits) {
		if (lastDiagonal != hits[i].diag) {
			//Different diagonal update the variables
			lastDiagonal = hits[i].diag;
			lastOffset = hits[i].posX - 1;
		}
		//We have a casting here because of a funny error
		//where the program was saying that -1 > 0
		if (lastOffset > ((int64_t) hits[i].posX)) {
			//The hit is covered by a previous frag in the same diagonal
			i++;
			continue;
		}

		nHitsUsed++;
		newFrag = FragFromHitReverse(coverage, &myF, &hits[i], sX, n0, ns0, sY, n1, ns1, nSeqs1, Lmin, SimTh,
							  WL);
		if (newFrag) {
			memcpy(&fragsBuf[nFrags],&myF,sizeof(struct FragFile));
			lastOffset = hits[i].posX + myF.length;
			nFrags++;
		}
		i++;
	}

	*nF = nFrags;
    *nHU = nHitsUsed;

    fprintf(stdout, "Reverse strand frags calculated\n");

    fragsBuf = realloc(fragsBuf,nFrags * sizeof(struct FragFile));
    if(fragsBuf == NULL)
        perror("Error reallocating memory of reverse frags array");

    return fragsBuf;
}


/**
 * Compute a fragments from one seed point
 * Similarirty thershold and length > mimL
 */
int FragFromHit(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
                uint64_t n0, uint64_t nsx, struct Sequence *sY,
                uint64_t n1, uint64_t nsy, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL) {
    int64_t ldiag, ldiag2;
    int64_t xfil, ycol;
    /* for version with backward search */
    int64_t xfil2, ycol2;
    int fragmentLength = 0;
    /* for version Maximum global---*/
    int64_t xfilmax, ycolmax;
    /* for version with backward search */
    int64_t xfilmax2, ycolmax2;
    int nIdentities = 0, maxIdentities = 0;
    char valueX, valueY;
    int fscore = 0, fscoreMax = 0; // full score

    uint64_t minLength = Lm;

    // Initialize values
    ldiag = min(n0 - H->posX, n1 - H->posY);
    //var to know how much we have until we reach the origin of coordinates
    ldiag2 = min(H->posX, H->posY);
    xfil = H->posX + WL;
    xfil2 = H->posX - 1;
    ycol = H->posY + WL;
    ycol2 = H->posY - 1;
    fragmentLength += WL;
    xfilmax = xfil;
    xfilmax2 = xfil2;
    ycolmax = ycol;
    ycolmax2 = ycol2;
    nIdentities = maxIdentities = WL;
    fscore = POINT * WL;
    fscoreMax = fscore;

    // now, looking for end_frag---
    while (fragmentLength < ldiag) {
        valueX = getValue(sX, xfil, nsx);
        valueY = getValue(sY, ycol, nsy);
        if (valueX == '*' || valueY == '*') {
            //separator between sequences
            break;
        }

        if (valueX == 'N' || valueY == 'N') {
            fscore -= 1;
        } else {
            if (valueX == valueY) {
                // match
                fscore += POINT;
                nIdentities++;
                if (fscoreMax < fscore) {
                    fscoreMax = fscore;
                    xfilmax = xfil;
                    ycolmax = ycol;
                    maxIdentities = nIdentities;
                }
            } else {
                fscore -= POINT;
            }
        }

        xfil++;
        ycol++;
        fragmentLength++;
        if (fscore < 0)
            break;
    }

    /**
     * Backward search --- Oscar (Sept.2013)
     **/
    fragmentLength = 0;
    fscore = fscoreMax;
    xfilmax2 = H->posX;
    ycolmax2 = H->posY;
    nIdentities = maxIdentities;
    if (xfil2 >= 0 && ycol2 >= 0)
        while (fragmentLength < ldiag2) {
            valueX = getValue(sX, xfil2, nsx);
            valueY = getValue(sY, ycol2, nsy);
            if (valueX == '*' || valueY == '*') {
                //separator between sequences
                break;
            }

            if (valueX == 'N' || valueY == 'N') {
                fscore -= 1;
            } else {
                if (valueX == valueY) {
                    // matches----
                    fscore += POINT;
                    nIdentities++;
                    if (fscoreMax < fscore) {
                        fscoreMax = fscore;
                        xfilmax2 = xfil2;
                        ycolmax2 = ycol2;
                        maxIdentities = nIdentities;
                    }
                } else {
                    fscore -= POINT;
                }
            }

            xfil2--;
            ycol2--;
            fragmentLength++;
            if (fscore < 0)
                break;
        }

    // Set the values of the FragFile
    myF->diag = H->diag;
    myF->xStart = (uint64_t) xfilmax2 - H->seqX;
    myF->yStart = (uint64_t) ycolmax2 - H->seqY;
    myF->xEnd = (uint64_t) xfilmax - H->seqX;
    myF->yEnd = (uint64_t) ycolmax - H->seqY;;
    myF->length = myF->xEnd - myF->xStart + 1;
    myF->score = fscoreMax;
    myF->ident = maxIdentities;
    myF->similarity = myF->score * 100.0
                      / scoreMax(&sX->datos[myF->xStart], &sY->datos[myF->yStart],
                                 myF->length, POINT);
    myF->seqX = H->seqX;
    myF->seqY = H->seqY;
    myF->block = 0;
    myF->strand = 'f';

    M[min(myF->length, 999)][(int)myF->similarity]++;

    if (myF->length > minLength && myF->similarity > SimTh)
        return 1;
    else
        return 0;
}

/**
 * Compute a fragments from one seed point
 * Similarirty thershold and length > mimL
 */
int FragFromHitReverse(long M[1000][100], struct FragFile *myF, hit *H, struct Sequence *sX,
				uint64_t n0, uint64_t nsx, struct Sequence *sY,
				uint64_t n1, uint64_t nsy, uint64_t nSeqs1, uint64_t Lm, uint64_t SimTh, int WL) {
	int64_t ldiag, ldiag2;
	int64_t xfil, ycol;
	/* for version with backward search */
	int64_t xfil2, ycol2;
	int fragmentLength = 0;
	/* for version Maximum global---*/
	int64_t xfilmax, ycolmax;
	/* for version with backward search */
	int64_t xfilmax2, ycolmax2;
	int nIdentities = 0, maxIdentities = 0;
	char valueX, valueY;
	int fscore = 0, fscoreMax = 0; // full score

	uint64_t minLength = Lm;

	// Initialize values
	ldiag = min(n0 - H->posY, n1 - H->posX);
	//var to know how much we have until we reach the origin of coordinates
	ldiag2 = min(H->posX, n1 - H->posY);;
	xfil = H->posX + WL;
	xfil2 = H->posX - 1;
	ycol = H->posY - WL;
	ycol2 = H->posY + 1;
	fragmentLength += WL;
	xfilmax = xfil;
	xfilmax2 = xfil2;
	ycolmax = ycol;
	ycolmax2 = ycol2;
	nIdentities = maxIdentities = WL;
	fscore = POINT * WL;
	fscoreMax = fscore;

	// now, looking for end_frag---
	while (fragmentLength < ldiag) {
		valueX = getValue(sX, xfil, nsx);
		valueY = complement(getValue(sY, ycol, nsy));
//        fprintf(stdout, "%c--%c\n",valueX,valueY);
		if (valueX == '*' || valueY == '*') {
			//separator between sequences
			break;
		}

		if (valueX == 'N' || valueY == 'N') {
			fscore -= 1;
		} else {
			if (valueX == valueY) {
				// match
				fscore += POINT;
				nIdentities++;
				if (fscoreMax < fscore) {
					fscoreMax = fscore;
					xfilmax = xfil;
					ycolmax = ycol;
					maxIdentities = nIdentities;
				}
			} else {
				fscore -= POINT;
			}
		}

		xfil++;
		ycol--;
		fragmentLength++;
		if (fscore < 0)
			break;
	}

	/**
     * Backward search --- Oscar (Sept.2013)
     **/
	fragmentLength = 0;
	fscore = fscoreMax;
	xfilmax2 = H->posX;
	ycolmax2 = H->posY;
	nIdentities = maxIdentities;
	if (xfil2 >= 0 && ycol2 < n1)
		while (fragmentLength < ldiag2) {
			valueX = getValue(sX, xfil2, nsx);
			valueY = complement(getValue(sY, ycol2, nsy));
			if (valueX == '*' || valueY == '*') {
				//separator between sequences
				break;
			}

			if (valueX == 'N' || valueY == 'N') {
				fscore -= 1;
			} else {
				if (valueX == valueY) {
					// matches----
					fscore += POINT;
					nIdentities++;
					if (fscoreMax < fscore) {
						fscoreMax = fscore;
						xfilmax2 = xfil2;
						ycolmax2 = ycol2;
						maxIdentities = nIdentities;
					}
				} else {
					fscore -= POINT;
				}
			}

			xfil2--;
			ycol2++;
			fragmentLength++;
			if (fscore < 0)
				break;
		}

	// Set the values of the FragFile
	myF->diag = H->diag;
	myF->xStart = (uint64_t) xfilmax2 - H->seqX;
	myF->yStart = (uint64_t) ycolmax2 - H->seqY;
	myF->xEnd = (uint64_t) xfilmax - H->seqX;
	myF->yEnd = (uint64_t) ycolmax - H->seqY;;
	myF->length = myF->xEnd - myF->xStart + 1;
	myF->score = fscoreMax;
	myF->ident = maxIdentities;
	myF->similarity = myF->score * 100.0
					  / scoreMax(&sX->datos[myF->xStart], &sY->datos[myF->yStart],
								 myF->length, POINT);
	myF->seqX = H->seqX;
	myF->seqY = H->seqY;
	myF->block = 0;
	myF->strand = 'r';

	M[min(myF->length, 999)][(int)myF->similarity]++;

	if (myF->length > minLength && myF->similarity > SimTh)
		return 1;
	else
		return 0;
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

