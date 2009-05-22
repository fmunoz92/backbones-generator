#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif

#include "petu.h"
#include <assert.h>
#include <popt.h>


// Datos a compartir por todos los niveles:
struct ArbolData
{
	float rgmax, cosfi[100], cossi[100], sinfi[100], sinsi[100];    // constantes
	ATOM  atm[3*NRESMAX];       // estructura parcial
	unsigned int nres, ndat;    // constantes
	long int cont;              // cantidad de estructuras exitosas hasta el momento
	XDRFILE* xfp;                    // file handler (constante)
	bool hubo_algun_exito;      // si encendido, dice que hubo al menos una rama que llego al final
};


// La estructura del algoritmo esta dividida en tres funciones, que se llaman en este orden:
static void generar_arbol(ArbolData* arbol_data);

static void generar_nivel_intermedio(unsigned int nivel, 
                                     float R_inicial[16], 
                                     unsigned int indice_nivel_anterior, 
                                     ArbolData* arbol_data);

bool procesar_ultimo_nivel(ArbolData* arbol_data); //bool void...?




int main (int argc , char **argv) 
{

// the following code is used to parse arguments with popt...

	 int    NivMax, 		//cantidad maxima de niveles 	
		Nres,		//cantidad de aminoacidos
		Ndat;		//cantidad de datos de entrada en el archivo data		

		
	float 	RN,			// Radio de nitrogeno
		RCa,			//Radio de carbono alfa
		RC,			//Radio de carbono
		Scal_1_4,		//
		Scal_1_5,
		DMax,   		//max distance
		RgMax;			//Radio de giro maximo

	int  argument;	    	/* argument passed by  line command*/
	char *File_Name;		/*name of vectors file*/
	struct poptOption options_table[] = {
		{"nres",(char)'r',POPT_ARG_INT,(void*)&Nres,0,"total structures","Nres"},
		{"ndat",(char)'t',POPT_ARG_INT,(void*)&Ndat,0,"total data in data file","Ndat"},
		{"rgmax",(char)'R',POPT_ARG_FLOAT,(void*)&RgMax,0,"max turn radius","RgMax"},
		{"dmax",(char)'d',POPT_ARG_FLOAT,(void*)&DMax,0,"Max distance","DMax"},
		{"rn",(char)'n',POPT_ARG_FLOAT,(void*)&RN,0,"nitrogen radius","RN"},
		{"rca",(char)'a',POPT_ARG_FLOAT,(void*)&RCa,0,"alpha carbon radius","RCa"},
		{"rc",(char)'c',POPT_ARG_FLOAT,(void*)&RC,0,"carbon radius","RC"},
		{"scal_1_4",(char)'s',POPT_ARG_FLOAT,(void*)&Scal_1_4,0,"follow with scal_1_4 value","Scal_1_4"},
		{"scal_1_5",(char)'l',POPT_ARG_FLOAT,(void*)&Scal_1_5,0,"follow with scal_1_5 value","Scal_1_5"},
		POPT_AUTOHELP{ NULL, 0, 0, NULL, 0 }    };	/*options table for line commands*/

 	poptContext opt_con;   

	opt_con = poptGetContext(NULL, argc, (const char **)argv,(struct poptOption* ) &options_table, 0);
	argument = poptGetNextOpt(opt_con);
	
	
	poptFreeContext(opt_con); /*frees memory destinated to context*/

// till here.



	FILE *filer,*filew,*fileras;
 
	ArbolData arbol_data;
  

	arbol_data.nres = Nres;
	arbol_data.ndat = Ndat;
 
	arbol_data.cont = 0;
	arbol_data.hubo_algun_exito = false;

	arbol_data.rgmax= RgMax;
	dmax2= DMax; dmax2 *= dmax2;
	setr(RN,RCa,RC,Scal_1_4,Scal_1_5);

	arbol_data.xfp = xdrfile_open("traj.xtc","w");
	filer=fopen("data","r");
	readdata(arbol_data.ndat, filer, arbol_data.cosfi, arbol_data.sinfi, arbol_data.cossi, arbol_data.sinsi);
	generar_arbol(&arbol_data);
	printf("%li\n",arbol_data.cont);        
	return(0);   
}


// =============================================== ALGORITMO ================================

void generar_arbol(ArbolData* arbol_data)
{
	float R_inicial[16];
	unsigned int i;

    //inicializar_arbol(arbol_data);  
	i = 1;
	while ( i <= arbol_data->ndat && !arbol_data->hubo_algun_exito )
	{
		clearatm( arbol_data->atm, arbol_data->nres);
		semilla(arbol_data->atm,R_inicial);
		generar_nivel_intermedio(2, R_inicial, i, arbol_data);
		i++;
	}
}

// se asume que nivel > 1
void generar_nivel_intermedio(unsigned int nivel, float R_inicial[16], unsigned int indice_nivel_anterior, ArbolData* arbol_data)
{
	float R_local[16];
	int resultado;
	bool exito = false; // solo aplicable si somos anteultimo nivel
	unsigned int i;

	assert(nivel > 1);  // pre condicion

	i = 1;
	while (i <= arbol_data->ndat && !exito)
	{
		
		copymat(R_local, R_inicial);

		resultado = poneres(R_local,
                        arbol_data->cossi[indice_nivel_anterior],   // OJO!!!!!!!!!!! i empieza de 1, esta bien eso??
                        arbol_data->sinsi[indice_nivel_anterior],   // idem
                        arbol_data->cosfi[i],                       // idem
                        arbol_data->sinfi[i],                       // idem
                        arbol_data->atm,
                        nivel);

		if ( resultado == BIEN )
		{
			if (nivel < arbol_data->nres)
			{              
				generar_nivel_intermedio(nivel+1, R_local, i, arbol_data);
			}
			else
			{
				exito = procesar_ultimo_nivel(arbol_data);
			}
		}
        
	i++;
	}
}

bool procesar_ultimo_nivel(ArbolData* arbol_data)
{
	bool exito = false;
    
	if ( calcRdG(arbol_data->atm, arbol_data->nres, arbol_data->rgmax) == BIEN )        
	{
		writextc(arbol_data->xfp, 
                 arbol_data->nres, 
                 arbol_data->cont, 
                 arbol_data->atm);
        /*imprime(filew,atm,nres,cont);*/       
		arbol_data->cont++;
		arbol_data->hubo_algun_exito = exito = true;
	}
	return exito;
}
