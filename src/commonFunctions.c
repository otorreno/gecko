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

unsigned long timestart() {
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return (tv.tv_usec / 1000) + (tv.tv_sec * 1000);
}

unsigned long timestop(unsigned long start) {
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return (tv.tv_usec / 1000) + (tv.tv_sec * 1000) - start;
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

struct Sequence *LeeSeqDB(FILE *f, uint64_t *n, uint64_t *nSeqs,
                          int fAst) {
    char c;
    uint64_t lon = 0, k = 0, seqs = 0;
    uint64_t lonFinal = 0;
    uint64_t SIZE = 0;
    uint64_t i = 0, r = 0;
    struct Sequence *sX;
    char *seq = NULL;

    //Initialize
    *n = 0;

    //Memory
    if ((sX = calloc(1, sizeof(struct Sequence))) == NULL)
        terror("Memory...");

    fseek(f, 0, SEEK_END);
    SIZE = ftell(f);
    fseek(f, 0, SEEK_SET);

    if ((sX->datos = calloc(SIZE, sizeof(char))) == NULL) {
        terror("Memory for sequence...");
    }

    if ((seq = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Memory for sequence...");
    }

    i = READBUF + 1;

    if (!fAst)
        while ((c = buffered_fgetc(seq, &i, &r, f)) != '>' && r > 0); //start seq
    if (feof(f))
        return 0;
    while ((c = buffered_fgetc(seq, &i, &r, f)) == ' ');

    while (k < MAXLID && c != '\n' && c != ' ') {
        if (r <= 0)
            return 0;

        sX->ident[k++] = c;
        c = buffered_fgetc(seq, &i, &r, f);
    }

    sX->ident[k] = 0; //end of data.
    while (c != '\n')
        c = buffered_fgetc(seq, &i, &r, f);
    c = buffered_fgetc(seq, &i, &r, f);

    //start list with sX2
    while (r > 0) {
        c = toupper(c);
        if (c == '>') {
            fAst = 1;
            seqs++;
            sX->datos[lon++] = '*';
            while (c != '\n') {
                if (r <= 0)
                    return 0;
                c = buffered_fgetc(seq, &i, &r, f);
            }
            //break;
        }
        if (isupper(c))
            sX->datos[lon++] = c;
        if (c == '*') {
            sX->datos[lon++] = c;
        }
        c = buffered_fgetc(seq, &i, &r, f);
    }

    free(seq);

    sX->datos[lon] = 0x00;

    lonFinal += lon;
    *nSeqs = seqs + 1;
    *n = lonFinal - seqs;
    return sX;
}

char getValue(struct Sequence *s, uint64_t pos) {
    return s->datos[pos];
}

long getSeqLength(struct Sequence *s, uint64_t start, int ns) {
    uint64_t s1 = 0;
    while (s1 > 0 && s->datos[s1] != '*') {
        s1--;
    }
    s1++;
    char *tmp = strchr(s->datos + s1, '*');
    if (tmp == NULL) {
        return strlen(s->datos) - s1 + 1;
    }
    return tmp - (s->datos + s1) + 1;
}
