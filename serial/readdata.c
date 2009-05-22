#include "petu.h"


void readdata(int ndat, FILE *filer, float *cosfi, float *sinfi, float *cossi, float *sinsi) 
{
	int i;
	float fi,si;
	for (i = 1; i <= ndat; i++) 
	{
		fscanf(filer,"%f %f",&fi,&si);
		cosfi[i] = cos(fi*0.017453);
		sinfi[i] = sin(fi*0.017453);    
		cossi[i] = cos(si*0.017453);
		sinsi[i] = sin(si*0.017453);    
	}
}
