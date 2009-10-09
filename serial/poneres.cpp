
#include "petu.h"

// En el modo DEBUG se deshabilitan los chequeos por lo que
// siempre devuelve BIEN.
#ifdef COMBINATIONS_DEBUG

int poneres(float *pR, float cossi, float sinsi, float cosfi, float sinfi, ATOM *patm, int resN, ArbolData *arbol_data, float dmax2, Residuo  &residuo)
{
	int at;
	float T[16];

	/* Varios de los parametros pueden ser eliminados, ya que son todos campos de arbol_data,
         pero no entiendo bien que rol juega indice_nivel_anterior en generar_nivel_intermedio().
	*/
	
	/*Guardo la anterior*/
	copymat(T,pR);
 
	at = 3*(resN-1)-1;
	at++;    
	int2car(T, b_C_N, cos_a_CA_C_N, sin_a_CA_C_N, cossi, sinsi, patm, at, N);

    
	at++;    
	int2car(T, b_N_CA, cos_a_C_N_CA, sin_a_C_N_CA, cos_OMEGA ,sin_OMEGA, patm, at, CA);


	at++;    
	int2car(T, b_CA_C, cos_a_N_CA_C, sin_a_N_CA_C, cosfi, sinfi, patm, at, C);

	//residuo = Residuo ( arbol_data->grilla->agregar_esfera(patm[at-2].x,patm[at-2].y,patm[at-2].z), arbol_data->grilla->agregar_esfera(patm[at-1].x,patm[at-1].y,patm[at-1].z), arbol_data->grilla->agregar_esfera(patm[at].x,patm[at].y,patm[at].z));
	//
residuo.at2 = arbol_data->grilla->agregar_esfera(patm[at-1].x,patm[at-1].y,patm[at-1].z);
	
	copymat(pR,T);
	return BIEN;      
}
#else
int poneres(float *pR, float cossi, float sinsi, float cosfi, float sinfi, ATOM *patm, int resN, ArbolData *arbol_data, float dmax2, Residuo  &residuo)
{
	int at;
	float T[16];

	/* Varios de los parametros pueden ser eliminados, ya que son todos campos de arbol_data,
         pero no entiendo bien que rol juega indice_nivel_anterior en generar_nivel_intermedio().
	*/

	
	/*Guardo la anterior*/
	copymat(T,pR);
 
	at = 3*(resN-1)-1;
	at++;    
	int2car(T, b_C_N, cos_a_CA_C_N, sin_a_CA_C_N, cossi, sinsi, patm, at, N);
	if(isclash(patm,at,resN) == MAL)
	{  
		return MAL;
	}

    
	at++;    
	int2car(T, b_N_CA, cos_a_C_N_CA, sin_a_C_N_CA, cos_OMEGA ,sin_OMEGA, patm, at, CA);
	if(isclash(patm,at,resN) == MAL)
	{  
		return MAL;
	}
	if(islong(patm,at,dmax2) == MAL)
	{ 
		return MAL;
	}


	at++;    
	int2car(T, b_CA_C, cos_a_N_CA_C, sin_a_N_CA_C, cosfi, sinfi, patm, at, C);
	if(isclash(patm,at,resN) == MAL)
	{  
		return MAL;
	}
	/*
	residuo = Residuo ( arbol_data->grilla->agregar_esfera(patm[at-2].x,patm[at-2].y,patm[at-2].z), arbol_data->grilla->agregar_esfera(patm[at-1].x,patm[at-1].y,patm[at-1].z), arbol_data->grilla->agregar_esfera(patm[at].x,patm[at].y,patm[at].z));
	*/
	residuo.at2 = arbol_data->grilla->agregar_esfera(patm[at-1].x,patm[at-1].y,patm[at-1].z);

	copymat(pR,T);
	return BIEN;      
}

#endif

