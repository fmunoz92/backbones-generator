#include "petu.h"
#include <math.h>

//This function calculates the gyration radius of the CA atoms
//It was meant to filter long chains. Now it is replaced by the volume filter

int calcRdG(ATOM *patm, int nres, float rgmax)
{  
	int i;
	float Rcm,xcm,ycm,zcm,dx2,dy2,dz2;
  
	Rcm=xcm=ycm=zcm=dx2=dy2=dz2=0;
  
	for (i=1;i<=(nres*3)-2;i+=3) 
	{
		xcm += patm[i].x;
		ycm += patm[i].y;
		zcm += patm[i].z;
	}    
	xcm /= nres;
	ycm /= nres;
	zcm /= nres;
	for (i=1;i<=(nres*3)-2;i+=3) 
	{
		dx2= (patm[i].x-xcm) * (patm[i].x-xcm);
		dy2= (patm[i].y-ycm) * (patm[i].y-ycm);
		dz2= (patm[i].z-zcm) * (patm[i].z-zcm);
		Rcm += (dx2+dy2+dz2);
	}
  
        // printf("radio de giro %f  %f\n",sqrt(Rcm/nres),rgmax);
  
	if(sqrt(Rcm/nres) > rgmax)
	{
		return MAL ;
	}
	else 
	{
		return BIEN;
	}    
}  

