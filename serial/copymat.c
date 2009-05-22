#include "petu.h"


void copymat(float *A, float *B) 
{
	int i;
	for(i=0;i<16;i++)
	A[i] = B[i];
}
