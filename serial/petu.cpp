#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif

#include "petu.h"
#include <assert.h>
#include <iostream>
#include <stdlib.h>

#include "getopt_pp_standalone.h"

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


static void show_usage();

int main (int argc , char **argv) 
{
    using namespace GetOpt;

	 int    NivMax, //cantidad maxima de niveles 	
            Nres,   //cantidad de aminoacidos
            Ndat;   //cantidad de datos de entrada en el archivo data		


	float   RN,         // Radio de nitrogeno
	        RCa,        // Radio de carbono alfa
            RC,         // Radio de carbono
            Scal_1_4,   // ??
            Scal_1_5,
            RgMax;      // Radio de giro maximo

#if 0
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
#endif

    GetOpt_pp ops(argc, argv);

    if (ops >> Option('r', "nres", Nres))
    {
        ops
            >> Option('t', "ndat", Ndat)    /* FIXME! Should be automatic */
            >> Option('R', "rgmax", RgMax)
            >> Option('n', "rn", RN)
            >> Option('a', "rca", RCa)
            >> Option('c', "rc", RC)
            >> Option('c', "rc", RC)
            >> Option('s', "scal_1_4", Scal_1_4)
            >> Option('l', "scal_1_5", Scal_1_5);

	        FILE *filer;
         
	        ArbolData arbol_data;
          

	        arbol_data.nres = Nres;
	        arbol_data.ndat = Ndat;
         
	        arbol_data.cont = 0;
	        arbol_data.hubo_algun_exito = false;

	        arbol_data.rgmax= RgMax;
	        setr(RN,RCa,RC,Scal_1_4,Scal_1_5);

	        arbol_data.xfp = xdrfile_open("traj.xtc","w");
	        filer=fopen("data","r");
	        readdata(arbol_data.ndat, filer, arbol_data.cosfi, arbol_data.sinfi, arbol_data.cossi, arbol_data.sinsi);
	        generar_arbol(&arbol_data);
	        printf("%li\n",arbol_data.cont);        

	        return EXIT_SUCCESS;
    }
    else
    {
        show_usage();
        return EXIT_FAILURE;
    }
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

void show_usage()
{
    std::cerr << "Invalid arguments." << std::endl;
}

