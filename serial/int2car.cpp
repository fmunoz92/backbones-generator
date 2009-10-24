#include "int2car.h"
#include "copymat.h"
//This performs the actual convertin of internal coordenates to cartesians

void int2car( float *pT, float bon, float cosang, float sinang, 
              float costor, float sintor, ATOM *patm, int at, AtomType tipo )
{
	float A[16],rot[11];

	rot[0]= -cosang;
	rot[1]= -sinang;
	rot[2]=  costor; 
	rot[3]= -sintor;           
	rot[4]=  bon*rot[0];
	rot[5]= -rot[1]*rot[2];
	rot[6]=  rot[0]*rot[2];
	rot[7]=  bon*rot[5];
	rot[8]=  rot[1]*rot[3];
	rot[9]= -rot[0]*rot[3];
	rot[10]= bon*rot[8];

	A[0]  = pT[0] *rot[0] + pT[1] *rot[5] + pT[2] *rot[8];
	A[4]  = pT[4] *rot[0] + pT[5] *rot[5] + pT[6] *rot[8]; 
	A[8]  = pT[8] *rot[0] + pT[9] *rot[5] + pT[10]*rot[8]; 
	A[12] = pT[12]*rot[0] + pT[13]*rot[5] + pT[14]*rot[8]; 
 
	A[1]  = pT[0] *rot[1] + pT[1] *rot[6] + pT[2] *rot[9];
	A[5]  = pT[4] *rot[1] + pT[5] *rot[6] + pT[6] *rot[9]; 
	A[9]  = pT[8] *rot[1] + pT[9] *rot[6] + pT[10]*rot[9]; 
	A[13] = pT[12]*rot[1] + pT[13]*rot[6] + pT[14]*rot[9];

	A[2]  = pT[1] *rot[3] + pT[2] *rot[2];
	A[6]  = pT[5] *rot[3] + pT[6] *rot[2];
	A[10] = pT[9] *rot[3] + pT[10]*rot[2];
	A[14] = pT[13]*rot[3] + pT[14]*rot[2];

	A[3]  = pT[0] *rot[4] + pT[1] *rot[7] + pT[2] *rot[10] + pT[3];	    
	A[7]  = pT[4] *rot[4] + pT[5] *rot[7] + pT[6] *rot[10] + pT[7];	    
	A[11] = pT[8] *rot[4] + pT[9] *rot[7] + pT[10]*rot[10] + pT[11];	    
	A[15] = pT[12]*rot[4] + pT[13]*rot[7] + pT[14]*rot[10] + pT[15];	    
 
	patm[at].x   = A[3];
	patm[at].y   = A[7];
	patm[at].z   = A[11];
	patm[at].vdw = tipo;
 
	copymat(pT,A);
}
