#include "petu.h"

//This function puts 0 in all the atoms coordinates

void clearatm( ATOM *patm, int nres)
{
	int i; 
	for (i=0; i < 3*nres; i++) 
	{
		patm[i].x = 0.0;
		patm[i].y = 0.0;
		patm[i].z = 0.0;    
	}
}
