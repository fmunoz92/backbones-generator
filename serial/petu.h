#ifndef PETU_H
#define PETU_H
#include <vector>

#include "grillado.h"

#define CPLUSPLUS
#include <xdrfile/xdrfile_xtc.h>
#undef CPLUSPLUS


enum AtomType {N,CA,C};

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



enum FilterResultType  {FILTER_FAIL, FILTER_OK};

//#define DEBUG

typedef struct {
  AtomType   vdw;
  float x;
  float y;
  float z;
} ATOM;


// Datos a compartir por todos los niveles:



			 
			 
			 
struct ArbolData
{
	float rgmax, dmax2;
	std::vector<float> cosfi, cossi, sinfi, sinsi;    // constantes
	ATOM  * atm;       // estructura parcial   
	unsigned int nres, ndat;    // constantes
	long int cont;              // cantidad de estructuras exitosas hasta el momento
	XDRFILE* xfp;                    // file handler (constante)
	bool hubo_algun_exito;      // si encendido, dice que hubo al menos una rama que llego al final
	Grillado *grilla; 		// Utilizamos el grillado para aproximar el volumen parcial
};
struct Residuo 
{		
		Residuo(const esferaId &param_at1, const esferaId &param_at2, const esferaId &param_at3) : at1(param_at1), at2(param_at2), at3(param_at3){};
		Residuo() {};
		esferaId at1;
		esferaId at2;
		esferaId at3;
}; 

/* a esto se le pone extern, para decir que quien importe el .h va a ser usuario de estas variables */
extern float r[3][3][3];


#endif
