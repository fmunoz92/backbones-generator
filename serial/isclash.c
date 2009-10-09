#include "petu.h"
#include <math.h>

float r[3][3][3];

//This function detect colitions between atoms
//The matrix r have the reference radius (actually the sum of the squares of the raduis)
//The last index of r is 0 for 1-4 clashes, 1 for 1-5 and 2 for the rest

int isclash(ATOM *patm, int at, int resN)
{
	int i;
	float d2,dx2,dy2,dz2;

        //This is to check for the so-called 1-4 clashes, 
        //i.e. a clash between atom at position i with atom at position i+3 
	i  = at-3;
	dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
	dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
	dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
	d2  = dx2+dy2+dz2;
	if (d2< r[patm[at].vdw][patm[i].vdw][0])
	{ 
#ifdef VERBOSE
	        printf("Clash 1-4 between atmom=%i and atom=%i distancia=%2.3f\n",at,i,sqrt(d2));
#endif	        
		return MAL;
	}
        
        //If this is atom #3 and passes the check for 1-4 clashes, then there is nothing else to check.
	if(at ==  3)
	{ 
		return BIEN;
	}

        //This is to check for the so-called 1-5 clashes; 
        //i.e. a clash between atom at position i with atom at position i+4 
	i = at-4;
	dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
	dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
	dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
	d2  = dx2+dy2+dz2;
	if (d2< r[patm[at].vdw][patm[i].vdw][1]) 
	{    
#ifdef VERBOSE
	        printf("Clash 1-5 between atmom=%i and atom=%i distancia=%2.3f\n",at,i,sqrt(d2));
#endif
		return MAL;
	}

        //If this is atom #4 and passes check for 1-4 and 1-5 clashes, then there is nothing else to check.
	if(at == 4) 
	{
		return BIEN;
	}
	
	//The rest of the clashes until the end of the chain.
	for (i = at-5; i >= 0; i--) 
	{
		dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
		dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
		dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
		d2  = dx2+dy2+dz2;
		if (d2< r[patm[at].vdw][patm[i].vdw][2])
		{ 
#ifdef VERBOSE
		        printf("Clash > 1-5 between atmom=%i and atom=%i distancia=%2.3f\n",at,i,sqrt(d2));
#endif
			return MAL;
		}
	}

	return BIEN;
}         	

