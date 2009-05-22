#include "petu.h"


int islong(ATOM *patm, int at, int dmax2)
{ 
	int i;
	float dx2,dy2,dz2,d2;
   
	if(at<9)
	{
		return BIEN;
	}
   
	for (i=at-9;i>=2;i-=3) 
	{     
		dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
		dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
		dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
		d2  = dx2+dy2+dz2;
		if (d2>dmax2)
		{ 
			return MAL;
		}
	}
	return BIEN;
}         	

