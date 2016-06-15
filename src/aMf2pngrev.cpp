/* aMf2pngrev: plots a png graph from MULTIPLE fragment files
 *             (based on af2pngrev.cpp - M stands for Multiple genomes)
 * syntax:
 *   af2Mpngrev FragsFileLIST.txt outFile.png nameX nameY
 *
 *   where FragsFileLIST.txt contains a filename in each line and the first one
 *         must be the "reference" genome (X axes)
 *   ortrelles at uma.es ----Aug.2014
 *
 * ------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <inttypes.h>
#include "pngwriter.h"
#include "structs.h"
#include "comparisonFunctions.h"

using namespace std;

#define min(x,y)    (((x) < (y)) ? (x) : (y))
#define fontpath "/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansMono.ttf"
#define NEGRO  0
#define ROJO   1
#define MAXFRAME  800
#define MAXLINE  1000

struct Files { 
   char *fname;
   uint64_t nx, ny;;
};




void drawMat(char*, char*, char*, string);
int transform(float, uint64_t);
int howManyFiles(char *);

void terror(string s) {
	cout << endl << "***ERR***" << s << endl;
	exit(-1);
}


int howManyFiles(char *fname){
	FILE *f;
	char line[MAXLINE];
	int nL=0;

	if ((f=fopen(fname,"rt"))==NULL)
	    terror((string)"Opening frags file LIST");


        if(fgets(line,MAXLINE,f)==NULL){
		terror("Not able to read a line from 'list of files' file");
	}
        while(!feof(f)) {
 		if (line[0]!='#' && (int)strlen(line)!=0) nL++;
        	if(fgets(line,MAXLINE,f)==NULL){
			terror("Not able to read a line from 'list of files' file");
		}
	}
	fclose(f);
	return nL;
}



// Load to memory a list of files --------------------------------------------------
// datafile format: fileName[tab]nSeq[tab][format][newLINE]

struct Files* LoadListOfFiles(char*fListName, int *nF) {

        FILE *f, *ff;
        struct Files*L=NULL;
        char line[MAXLINE];
        int N=0,nFiles;
		uint64_t xnx,xny;

	nFiles = howManyFiles(fListName);

	if ((L=(struct Files*) calloc(nFiles,sizeof(struct Files)))==NULL) 
           terror((string)"memory for list of files");

        if ((f=fopen(fListName,"rt"))==NULL) 
		terror((string)"opening List of Files");

       	if(fgets(line,MAXLINE,f)==NULL){
		terror("Not able to read a line from 'list of files' file");
	}
        while(!feof(f)) {

 		if (line[strlen(line)-1]=='\n') line[strlen(line)-1]=0x00;
		 
 		if (line[0]!='#' && (int)strlen(line)!=0) {
			L[N].fname = strdup(line);         

			if ((ff=fopen(L[N].fname,"rt"))==NULL)
				terror((string)"Opening frags file from LIST X X XX");

    		readSequenceLength(&xnx, ff);
		    readSequenceLength(&xny, ff);

			L[N].nx = xnx;
			L[N].ny = xny;
			if (N!=0 && L[N].nx !=L[N-1].nx)
				terror((string)"uncoherent NX values --reference genome is not the same");
			N++;
			fclose(ff);
		}
       		if(fgets(line,MAXLINE,f)==NULL){
			terror("Not able to read a line from 'list of files' file");
		}
	}
	fclose(f);
	(*nF) = nFiles;
	return L;
}


void drawMat(char *fname, char *namex, char* namey, string outputFile) {
    // black, blue, red, green, yellow,magenta, cyan, navy, darkGreen , purple
	double R[] = { 0., 0., 1., 0., 1., 1., 0., 0.0, 0.0, 0.5 };
	double G[] = { 0., 0., 0., 1., 1., 0., 1., 0.0, 0.5, 0.0 };
	double B[] = { 0., 1., 0., 0., 0., 1., 1., 0.5, 0.0, 0.5 };

	int sizeH, sizeV;
	int fil, col;
	int border = 5;
	int txtEdge = 100;
	char stmp[256];
	float scale;
	uint64_t i;
	int ij;
	FILE *fDP;
	uint64_t xnx, xny;
	uint64_t nFrags = 0, nFragsRev = 0;
	int elemSize;
	struct FragFile frag;
	struct Files *L;
	int nFiles, NX;
	unsigned int maxNY=0;
	int COLOR;

	L = LoadListOfFiles(fname, &nFiles);
	for (ij=0;ij<nFiles; ij++) 
		if (L[ij].ny > maxNY) maxNY=L[ij].ny; 
	NX = L[0].nx;

	// global scale (ONE for all) to compute framework-size & elemSize
	scale = min((float)MAXFRAME/(float)NX, (float)MAXFRAME/(float)maxNY); 
	if (scale > 1)
		scale = 1;
	if (scale < 1)
		elemSize = 1;
	else
		elemSize = (int) scale;

	sizeH = (int) (NX * scale * elemSize + 2 * border);
	sizeV = (int) (maxNY * scale * elemSize + txtEdge);

	cout << "nx=" << NX << " ny=" << maxNY << " Scale=" << scale << " sizeH=" << sizeH << " sizeV=" << sizeV << " ElemSize=" << elemSize << endl;

	pngwriter png(sizeH, sizeV, 65500, outputFile.c_str());



	for (ij=0;ij<nFiles; ij++) {

		if ((fDP = fopen(L[ij].fname, "rb")) == NULL)
			terror((string)"Opening fragments binary file");

    	readSequenceLength(&xnx, fDP); // skip
	    readSequenceLength(&xny, fDP); // skip
		// individual scale
		scale = min((float)MAXFRAME/(float)xnx, (float)MAXFRAME/(float)xny); 
 		COLOR = ij%10;
  		readFragment(&frag, fDP);

		//DRAW HSPs
		while (!feof(fDP)) {
			if(frag.strand=='f'){
				nFrags++;
			}else{
				nFragsRev++;
			}
			for (i = 0; i < frag.length; i++) {
				col = border + transform(scale, frag.xStart + i);
				if (frag.strand == 'f') {
					fil = border + txtEdge + transform(scale, frag.yStart + i);
					png.filledsquare(col, fil, col + elemSize - 1,
						fil + elemSize - 1, R[COLOR], G[COLOR], B[COLOR]);
						//fil + elemSize - 1, R[ROJO], G[NEGRO], B[NEGRO]);
				} else {
					fil = txtEdge + transform(scale, (maxNY - frag.yStart) - i);
					png.filledsquare(col, fil, col + elemSize - 1,
						fil + elemSize - 1, R[COLOR], G[COLOR], B[COLOR]);
						//fil + elemSize - 1, R[NEGRO], G[NEGRO], B[NEGRO]);
				}
			}

  			readFragment(&frag, fDP);
		}

		fclose(fDP);
	}

	cout << "nFrags: " << nFrags << " nFragsRev: " << nFragsRev << endl;

	//X-axis
	png.line(border - 2, txtEdge + border - 2, sizeH - border,
			txtEdge + border - 2, 0.0, 0.0, 0.0);
	//Y-axis
	png.line(border - 2, txtEdge + border - 2, border - 2, sizeV - border, 0.0,
			0.0, 0.0);

	//DRAW TITLE & RANGES
	std::ostringstream out;
	out << "nx(Xaxis)=" << NX << " ny(Yaxis)=" << maxNY;
	png.plot_text((char*) fontpath, 12, 20, 80, 0.0, strdup(out.str().c_str()), R[NEGRO], G[NEGRO],
			B[NEGRO]);

	out.str("");
	png.filledsquare(21, 60, 42, 70, R[ROJO], G[ROJO], B[ROJO]);
	out << "Forward. (" << nFrags << ")";
	png.plot_text((char*) fontpath, 10, 50, 62, 0.0, strdup(out.str().c_str()), R[NEGRO], G[NEGRO],
			B[NEGRO]);

	out.str("");
	png.filledsquare(21, 40, 42, 50, R[NEGRO], G[NEGRO], B[NEGRO]);
	out << "Reverse. (" << nFragsRev << ")";
	png.plot_text((char*) fontpath, 10, 50, 42, 0.0, strdup(out.str().c_str()), R[NEGRO], G[NEGRO],
			B[NEGRO]);

	out.str("");
	out << namex << " VS " << namey;
	cout << "sizeH=" << sizeH << ", strlen(nombres)=" << strlen(strdup(out.str().c_str())) << endl;
	int posXnombre = (sizeH - border - strlen(stmp)) / 2;
	png.plot_text((char*) fontpath, 10, posXnombre, 22, 0.0, strdup(out.str().c_str()), R[NEGRO],
			G[NEGRO], B[NEGRO]);

	png.close();
}

int transform(float Scale, uint64_t v) {
	return (int) (Scale * v);
}



int main(int ac, char **av) {
	if (ac != 5)
		terror((string)"Use aMf2pngrev FragsFileLIST.txt outFile.png nameX nameY");

	drawMat(av[1], av[3], av[4], string(av[2]));

	return 1;

}

