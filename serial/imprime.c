#include "petu.h"


void imprime(FILE *filew, ATOM *patm, int Nres, int Mod)
{   
	int i,j;
	j=0;  
	fprintf(filew,"MODEL    %5i\n",Mod);
	for (i=1;i<=Nres;i++) 
	{
		j++;
		fprintf(filew,"ATOM  %5i  N   ALA   %3i    %8.3f%8.3f%8.3f\n",j,i,patm[j].x,patm[j].y,patm[j].z);
		j++;
		fprintf(filew,"ATOM  %5i  CA  ALA   %3i    %8.3f%8.3f%8.3f\n",j,i,patm[j].x,patm[j].y,patm[j].z);
		j++;
		fprintf(filew,"ATOM  %5i  C   ALA   %3i    %8.3f%8.3f%8.3f\n",j,i,patm[j].x,patm[j].y,patm[j].z);
	}
	fprintf(filew,"ENDMDL\n");
}
