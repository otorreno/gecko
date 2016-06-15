/* hits : determine HITS between two sequences based on their 
 dictionaries from disk
 Syntax: hits prefixNameX prefixNameY Outfile

 prefixNameX & Y :refers to *.d2hW : index of words-Pos-Ocurrences
 and *.d2hP : words positions
 Outfile is in the form of [Diag][posX][posY]

 Jan.2012: define an I-O buffer to reduce disk activity
 Feb.2012: - define a function to compare words instead of chars
 - use static buffer for positions
 - New parameter: FreqThre (frequency threshold) to avoid high
 frequency words (word.freq> FreqThr are skipped)

 Feb.6   : new parameter (same meaning as in w2hd): PrefixSize

 May 22  : for long repetitions some of them will be skipped (the step is
 computed as NREP / MaxREP

 ortrelles@uma.es / Dic.2011
 ---------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "structs.h"
#include "commonFunctions.h"
#include "dictionaryFunctions.h"
#include "comparisonFunctions.h"

#define MAXBUF 10000000
#define MaxREP 10000

int main(int ac, char** av) {

	char fname[1024];
	int hitsInBuf = 0;
	int i, j;
	int comp;
	int freqThr;
	uint64_t nHits = 0, wordMatches = 0;
	int wSize;
	int stepX, stepY;
	FILE *fX, *fY, *fOut;

	location *posX = NULL, *posY = NULL;
	hashentry heX, heY;
	hit *hBuf;

	if (ac != 6)
		terror(
				"USE: hits prefixNameX prefixNameY  Outfile FreqThre PrefixSize");
	freqThr = atoi(av[4]);
	wSize = atoi(av[5]);

	// I-O buffer
	if ((hBuf = (hit*) calloc(sizeof(hit), MAXBUF)) == NULL)
		terror("HITS: memory for I-O buffer");
	// word positions buffers
	if ((posX = (location*) calloc(MAXBUF, sizeof(location))) == NULL)
		terror("memory for posX buffer.. using MAXBUF=10MB");
	if ((posY = (location*) calloc(MAXBUF, sizeof(location))) == NULL)
		terror("memory for posY buffer.. using MAXBUF=10MB");

	// Sequence X file
	sprintf(fname, "%s.dict", av[1]);
	if ((fX = fopen(fname, "rb")) == NULL) {
		terror("opening seqX.dict file");
	}

	// Sequence Y file
	sprintf(fname, "%s.dict", av[2]);
	if ((fY = fopen(fname, "rb")) == NULL) {
		terror("opening seqX.dict file");
	}

	// OUT file
	if ((fOut = fopen(av[3], "wb")) == NULL)
		terror("opening OUT file");

	// kick-off
	if (readHashEntry(&heX, fX, freqThr) == -1)
		terror("no hash (1)");
	if (readHashEntry(&heY, fY, freqThr) == -1)
		terror("no hash (2)");

	while (!feof(fX) && !feof(fY)) {

		comp = wordcmp(&heX.w.b[0], &heY.w.b[0], wSize);
		if (comp < 0) {
			if (fread(posX, sizeof(location), heX.num, fX) != heX.num)
				terror("Error reading word occurrences in seqX");
			readHashEntry(&heX, fX, freqThr);
			continue;
		}
		if (comp > 0) {
			if (fread(posY, sizeof(location), heY.num, fY) != heY.num)
				terror("Error reading word occurrences in seqY");
			readHashEntry(&heY, fY, freqThr);
			continue;
		}
		wordMatches++;

		// Load word position for seqX
		loadWordOcurrences(heX, &posX, &fX);
		// Load word position for seqY
		loadWordOcurrences(heY, &posY, &fY);

		// Hits-----------------------
		if (heX.num > MaxREP)
			stepX = heX.num / MaxREP;
		else
			stepX = 1;
		if (heY.num > MaxREP)
			stepY = heY.num / MaxREP;
		else
			stepY = 1;

		for (i = 0; i < heX.num; i += stepX)
			for (j = 0; j < heY.num; j += stepY) {
				hBuf[hitsInBuf].diag = posX[i].pos - posY[j].pos;
				hBuf[hitsInBuf].posX = posX[i].pos;
				hBuf[hitsInBuf].seqX = posX[i].seq;
				hBuf[hitsInBuf].posY = posY[j].pos;
				hBuf[hitsInBuf].seqY = posY[j].seq;

				hitsInBuf++;
				if (hitsInBuf == MAXBUF - 1) {
					fwrite(hBuf, sizeof(hit), hitsInBuf, fOut);
					hitsInBuf = 0;
				}
			}

		nHits += ((heX.num / stepX) * (heY.num / stepY));

		readHashEntry(&heX, fX, freqThr);
		readHashEntry(&heY, fY, freqThr);

	}

	//Closing dictionary files
	fclose(fX);
	fclose(fY);

	//Checking if there is something still at the buffer
	if (hitsInBuf != 0) {
		fwrite(hBuf, sizeof(hit), hitsInBuf, fOut);
	}
	fclose(fOut);

	fprintf(stdout, "HITS: matches=%" PRIu64 " Tot HITS=%" PRIu64 "\n", wordMatches, nHits);

	return 0;
}

