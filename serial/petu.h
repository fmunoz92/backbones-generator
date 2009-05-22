#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <xdrfile/xdrfile_xtc.h>

#define N  0
#define CA 1
#define C  2

#define b_C_N  1.330
#define b_N_CA 1.460
#define b_CA_C 1.525

#define a_CA_C_N 2.059454 
#define cos_a_CA_C_N -0.469441085 
#define sin_a_CA_C_N  0.882963797

#define a_C_N_CA 2.111813
#define cos_a_C_N_CA -0.515007735 
#define sin_a_C_N_CA  0.85718553

#define a_N_CA_C 2.024582156
#define cos_a_N_CA_C -0.438371348
#define sin_a_N_CA_C  0.898793948

#define OMEGA 3.124087
#define cos_OMEGA -0.999847695
#define sin_OMEGA 0.017452406

#define MAL  0
#define BIEN 1

#define NRESMAX 50

/*#define MAX(A,B)  ((A)>(B)?(A):(B))*/
/*#define DEBUG*/

typedef struct {
  int   vdw;
  float x;
  float y;
  float z;
} ATOM;

#ifdef CPLUSPLUS
extern "C" {
#endif


/* a esto se le pone extern, para decir que quien importe el .h va a ser usuario de estas variables */
extern float r[3][3][3];
extern float dmax2;

#ifdef DEBUG
int   poneres( float *, float, float, float, float, ATOM *, int, FILE *);
int   isclash(ATOM *, int, int , FILE *);
#else
int   poneres( float *, float, float, float, float, ATOM *, int);
int   isclash(ATOM *, int, int);
#endif
int   islong(ATOM *, int , int );
int   calcRdG(ATOM *, int , float);
void  setr(float , float, float, float, float);
void  int2car(float *, float , float , float, float, float ,ATOM *, int , int );
void  copymat(float *, float *); 
void  semilla(ATOM *, float *);
void  imprime(FILE *, ATOM *, int , int );
void  readdata(int , FILE *, float *, float *, float *, float *); 
void  clearatm(ATOM *, int);
void  writextc (XDRFILE* xfp, int nres, int n, ATOM *patm);

#ifdef CPLUSPLUS
};
#endif
