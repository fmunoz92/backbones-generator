#include "petu.h"


int poneres(float *pR, float cossi, float sinsi, float cosfi, float sinfi, ATOM *patm, int resN)
{
	int at;
	float T[16];

	/*Guardo la anterior*/
	copymat(T,pR);
 
	at = 3*(resN-1); 
	at++;    
	int2car(T, b_C_N, cos_a_CA_C_N, sin_a_CA_C_N, cossi, sinsi, patm, at, N);
	if(isclash(patm,at,resN) == MAL)
	{  
		return MAL;
	}

    
	at++;    
	int2car(T, b_N_CA, cos_a_C_N_CA, sin_a_C_N_CA, cos_OMEGA ,sin_OMEGA, patm, at, CA);
	if(isclash(patm,at,resN) == MAL)
	{  
		return MAL;
	}
	if(islong(patm,at,dmax2) == MAL)
	{ 
		return MAL;
	}


	at++;    
	int2car(T, b_CA_C, cos_a_N_CA_C, sin_a_N_CA_C, cosfi, sinfi, patm, at, C);
	if(isclash(patm,at,resN) == MAL)
	{  
		return MAL;
	}

 
	copymat(pR,T);
	return BIEN;      
}

