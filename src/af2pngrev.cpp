/* Laboratory : af2pngrev plots a png graph from a fragment file
 *
 * syntax:
 *   af2pngrev FragsFile.BIN outFile.png nameX nameY
 * ------------------------------------------*/

#include <stdio.h>
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
#define ROJO   2
#define MAXFRAME 800

void drawMat(char*, char*, char*, string);
int transform(float, uint64_t);

void terror(string s) {
	cout << endl << "***ERR***" << s << endl;
	exit(-1);
}

int main(int ac, char **av) {
	if (ac != 5)
		terror("Use af2pngrev FragsFile.BIN outFile.png nameX nameY");

	drawMat(av[1], av[3], av[4], string(av[2]));

	return 1;

}

void drawMat(char *fname, char *namex, char* namey, string outputFile) {

	double R[] = { 0., 0., 1., 1., 0., 0.6, 1.0, 1.0 };
	double G[] = { 0., 0., 0., 0.55, 1., 0.4, 1.0, 1.0 };
	double B[] = { 0., 1., 0., 0., 0., 1.0, 0.0, 1.0 };

	int sizeH, sizeV;
	int fil, col;
	int border = 5;
	int txtEdge = 100;
	char stmp[256];
	float scale;
	uint64_t i;
	FILE *fDP;
	uint64_t nx, ny;
	uint64_t nFrags = 0, nFragsRev = 0;
	int elemSize;
	struct FragFile frag;

	if ((fDP = fopen(fname, "rb")) == NULL)
		terror("Opening DP matrix binary file");

	readSequenceLength(&nx, fDP);
	readSequenceLength(&ny, fDP);

	scale = min((float)MAXFRAME/(float)nx, (float)MAXFRAME/(float)ny);
	if (scale > 1)
		scale = 1;
	if (scale < 1)
		elemSize = 1;
	else
		elemSize = (int) scale;

	sizeH = (int) (nx * scale * elemSize + 2 * border);
	sizeV = (int) (ny * scale * elemSize + txtEdge);
	cout << "nx=" << nx << " ny=" << ny << " Scale=" << scale << " sizeH=" << sizeH << " sizeV=" << sizeV << " ElemSize=" << elemSize << endl;

	pngwriter png(sizeH, sizeV, 65500, outputFile.c_str());

	readFragment(&frag, fDP);
	//DRAW SPOTS
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
						fil + elemSize - 1, R[ROJO], G[NEGRO], B[NEGRO]);
			} else {
				fil = txtEdge + transform(scale, (ny - frag.yStart) - i);
				png.filledsquare(col, fil, col + elemSize - 1,
						fil + elemSize - 1, R[NEGRO], G[NEGRO], B[NEGRO]);
			}
		}

		readFragment(&frag, fDP);
	}
	cout << "nFrags: " << nFrags << " nFragsRev: " << nFragsRev << endl;
	fclose(fDP);

	//X-axis
	png.line(border - 2, txtEdge + border - 2, sizeH - border,
			txtEdge + border - 2, 0.0, 0.0, 0.0);
	//Y-axis
	png.line(border - 2, txtEdge + border - 2, border - 2, sizeV - border, 0.0,
			0.0, 0.0);

	//DRAW TITLE & RANGES
	std::ostringstream out;
	out << "nx(Xaxis)=" << nx << " ny(Yaxis)=" << ny;
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

