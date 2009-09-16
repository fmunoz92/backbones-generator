#include "petu.h"

//This function filter long chains.
//Calculates the CA-CA distances of the recently added CA with the rest.
//Then compares this distances with dmax2, wich cames from the analisys of protein 
//structures experimetally determined

int islong(ATOM *patm, int at, float dmax2)
{ 
	int i;
	float dx2,dy2,dz2,d2;
   
        // Until we reach residue #5, the check for long chain is meaningless
        // at=13 is the CA of residue #5
	if(at<13)
	{
		return BIEN;
	}
   
        // at-12 is the CA atom that is four residues down the chain
        // at=1 is the first CA in the chain
	for (i=at-12;i>=1;i-=3) 
	{     
		dx2 = (patm[at].x - patm[i].x) * (patm[at].x - patm[i].x);
		dy2 = (patm[at].y - patm[i].y) * (patm[at].y - patm[i].y);
		dz2 = (patm[at].z - patm[i].z) * (patm[at].z - patm[i].z);
		d2  = dx2+dy2+dz2;
		if (d2>dmax2)
		{ 
		        printf("Largo=%f\n",d2);
			return MAL;
		}
	}
	return BIEN;
}         	

