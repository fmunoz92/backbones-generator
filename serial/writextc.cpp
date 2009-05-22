#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif

#include "petu.h"

enum Ejes
{
	EjeX,
	EjeY,
	EjeZ,
	// -------------------
	CantidadDeEjes  // dejar este valor a lo ultimo
};

typedef float EjesCartesianos[CantidadDeEjes];

void writextc (int xfp, int nres, int n, ATOM *patm)
{

	const int natom=3*nres;
	int i;  
	EjesCartesianos box[3];
	EjesCartesianos* cxtc = new EjesCartesianos[natom]; /* cosas para escribir el xtc*/
  
	box[0][EjeX]= 9.0;  box[0][EjeY]= 0.0;  box[0][EjeZ]= 0.0;
	box[1][EjeX]= 0.0;  box[1][EjeY]= 9.0;  box[1][EjeZ]= 0.0;
	box[2][EjeX]= 0.0;  box[2][EjeY]= 0.0;  box[2][EjeZ]= 9.0;

	for (i=0;i<natom;i++) 
	{
		cxtc[i][EjeX]=patm[i+1].x/10.0;
		cxtc[i][EjeY]=patm[i+1].y/10.0;
		cxtc[i][EjeZ]=patm[i+1].z/10.0;	
	}
	write_xtc(xfp,natom,n,(float) n*1.0,box,cxtc,1000.0);

	delete[] cxtc;
}
