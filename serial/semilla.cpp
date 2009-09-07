#include "petu.h"


void semilla(ArbolData *arbol_data, float *R, Residuo & residuo)
{
  /* Los valores semilla */
	arbol_data->atm[1].x=0.000; arbol_data->atm[1].y=0.000; arbol_data->atm[1].z=0.000;
	arbol_data->atm[2].x=1.460; arbol_data->atm[2].y=0.000; arbol_data->atm[2].z=0.000;
	arbol_data->atm[3].x=2.011; arbol_data->atm[3].y=1.422; arbol_data->atm[3].z=0.000;
	arbol_data->atm[1].vdw=0  ; arbol_data->atm[2].vdw=1  ; arbol_data->atm[3].vdw=2  ;
	  
	R[0]  = 0.361594; R[1]  =-0.932336;
	R[2]  = 0.000000; R[3]  = 2.011431;
	R[4]  = 0.932336; R[5]  = 0.361594;
	R[6]  = 0.000000; R[7]  = 1.421812;
	R[8]  = 0.000000; R[9]  = 0.000000;
	R[10] = 1.000000; R[11] = 0.000000;
	R[12] = 0.000000; R[13] = 0.000000;
	R[14] = 0.000000; R[15] = 1.000000;
	residuo = Residuo(arbol_data->grilla->agregar_esfera(arbol_data->atm[1].x,arbol_data->atm[1].y,arbol_data->atm[1].z),
	arbol_data->grilla->agregar_esfera(arbol_data->atm[2].x,arbol_data->atm[2].y,arbol_data->atm[2].z),
	arbol_data->grilla->agregar_esfera(arbol_data->atm[3].x,arbol_data->atm[3].y,arbol_data->atm[3].z));
	
}
