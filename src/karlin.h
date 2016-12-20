#if !defined(KARLIN)
#define KARLIN
#include <math.h>
#include "stdinc.h"
#include "alphabet.h"

#define KARLINMAXIT 50  /* Max. # iterations used in calculating K */
boolean	karlin(int low,int high,double *pr,double *lambda,double *K,double *H);
double ExpectedInformation(a_type A, double lambda, double *freq);
int 	karlin_gcd(int a,int b);

#endif

