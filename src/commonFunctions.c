#include <stdio.h>
#include <stdlib.h>

void terror(char *s) {
	printf("ERR**** %s ****\n", s);
	exit(-1);
}
