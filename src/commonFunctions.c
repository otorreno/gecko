#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

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