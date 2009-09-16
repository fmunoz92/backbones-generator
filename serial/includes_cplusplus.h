#include <fstream>
#include <vector>
#include "grillado.h"
#include <cmath>
			 
			 
			 
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
			 
			 
int   poneres( float *, float, float, float, float, ATOM *, int, ArbolData *, float dmax2, Residuo &);
void  readdata( std::ifstream &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &); 
void  semilla(ArbolData *, float *, Residuo &);


inline float cubic_root(float a) {
	return (pow(a, 1.0f/3.0f));
}


// Devuelve 1 si value esta entre min y max.
// Se separo en una funcion distinta para mejorar el rendimiento de volumen_en_rango.
inline bool in_range(float value, float min, float max) {
	return value < max && value > min;
}

inline void sacar_residuo(ArbolData * arbol_data, Residuo & residuo) {
	//arbol_data->grilla->sacar_esfera(residuo.at1);
	arbol_data->grilla->sacar_esfera(residuo.at2);
	//arbol_data->grilla->sacar_esfera(residuo.at3);
}


