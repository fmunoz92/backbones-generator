#ifndef PETU_H
#define PETU_H
#include <vector>

#include "grillado.h"
#include "prot-filer/definitions.h"

enum FilterResultType  {FILTER_FAIL, FILTER_OK};


#include "prot-filer/format_filer.h"
#include "prot-filer/angles.h"

// Datos a compartir por todos los niveles:

class FormatFiler;
			 
struct ArbolData
{
	float rgmax, dmax2;
	std::vector<float> cosfi, cossi, sinfi, sinsi;    // constantes
	ATOM  * atm;       // estructura parcial   
	unsigned int nres, ndat;    // constantes
	long int cont;              // cantidad de estructuras exitosas hasta el momento
	bool hubo_algun_exito;      // si encendido, dice que hubo al menos una rama que llego al final
	Grillado *grilla; 		// Utilizamos el grillado para aproximar el volumen parcial
	AnglesData *angles_data; // Used only when writing compressed data.
	FormatFiler *filer;
	AnglesMapping *angles_mapping;
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
