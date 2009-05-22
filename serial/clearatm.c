#include "petu.h"

void clearatm( ATOM *patm, int nres)
{
	int i; 
	for (i=1; i <= 3*nres; i++) 
	{
		patm[i].x = 0.0;
		patm[i].y = 0.0;
		patm[i].z = 0.0;    
	}
}
