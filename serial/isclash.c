#include "petu.h"

float r[3][3][3];
float dmax2;


int isclash(ATOM *patm, int at, int resN)
{
	int i,to;
	float d2,dx2,dy2,dz2;

	i  = at-3;
	dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
	dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
	dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
	d2  = dx2+dy2+dz2;
	if (d2< r[patm[at].vdw][patm[i].vdw][0])
	{ 
		return MAL;
	}
 
	if(at <=  4)
	{ 
		return BIEN;
	}
	i = at-4;
	dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
	dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
	dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
	d2  = dx2+dy2+dz2;
	if (d2< r[patm[at].vdw][patm[i].vdw][1]) 
	{    
		return MAL;
	}

	if(at <= 5) 
	{
		return BIEN;
	}
	for (i = at-5; i >= 1; i--) 
	{
		dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
		dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
		dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
		d2  = dx2+dy2+dz2;
		if (d2< r[patm[at].vdw][patm[i].vdw][2])
		{ 
			return MAL;
		}
	}

	return BIEN;
}         	

