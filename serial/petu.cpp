




#include <cassert>
#include <iostream>
#include <fstream>
#include "getopt_pp_standalone.h"

#include <mili/mili.h>


#include "poneres.h"
#include "readdata.h"
#include "semilla.h"
#include "copymat.h"
#include "calcrdg.h"
#include "setr.h"


#include "imprime.h"

#include "clearatm.h"
#include "writextc.h"



inline float cubic_root(float a) {
	return (pow(a, 1.0f/3.0f));
}

// Devuelve 1 si value esta entre min y max.
// Se separo en una funcion distinta para mejorar el rendimiento de volumen_en_rango.
inline bool in_range(float value, float min, float max) {
	return min < bchain(value) < max;
}

inline void sacar_residuo(ArbolData * arbol_data, Residuo & residuo) {
	//arbol_data->grilla->sacar_esfera(residuo.at1);
	arbol_data->grilla->sacar_esfera(residuo.at2);
	//arbol_data->grilla->sacar_esfera(residuo.at3);
}


// La estructura del algoritmo esta dividida en tres funciones, que se llaman en este orden:
static void generar_arbol(ArbolData* arbol_data);

static void generar_nivel_intermedio(unsigned int nivel, 
                                     float R_inicial[16], 
                                     unsigned int indice_nivel_anterior, 
                                     ArbolData* arbol_data);

static bool procesar_ultimo_nivel(ArbolData* arbol_data); //bool void...?


static FilterResultType volumen_en_rango(const ArbolData * arbol_data);
static FilterResultType filtros_ultimo_nivel(const ArbolData * arbol_data) ;
static void show_usage();



const float cota_maxima_volumen = 177.65f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.
const float pendiente_empirica = -0.0882f ; // Pendiente obtenida a partir de las pruebas de un set de datos en Grillado.
const float volumen_min_aa = 110.0f;
//const float volumen_min_aa = 0.0f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.
	

int main (int argc , char **argv) 
{
    using namespace GetOpt;

    int     Nres;   // Number of amino acids in the chains to build.

    float   RN,         // Radius of the Nitrogen atom.
            RCa,        // Radius of the Carbon atom.
            RC,         // Radius of the Carbon atom.
            Scal_1_4,   // Scaling factor for the radii. Used to check for 1-4 clashes.
            Scal_1_5;   // Scaling factor for the radii. Used to check for 1-5 clashes.

    std::string data;   // Name of the input file.
	
    size_t m; // They indicate the size of each dimention of the grid.
    size_t n;
    size_t z;



    GetOpt_pp ops(argc, argv);
		 
    if (ops >> Option('r', "Nres", Nres))
    {
	
        ops
            >> Option('n', "Rn", RN, 1.5f)
            >> Option('a', "Rca", RCa, 1.7f)
            >> Option('c', "Rc", RC, 1.6f)
            >> Option('s', "Scal_1_4", Scal_1_4, 0.85f)
            >> Option('l', "Scal_1_5", Scal_1_5, 1.0f)
	    >> Option('i', "input_file", data, "ramachandran.dat")
	    >> Option('N', "rows", n, static_cast<size_t>(100))
            >> Option('M', "cols", m, static_cast<size_t>(100))
	    >> Option('Z', "depth", z, static_cast<size_t>(100));
		// Hay que decidir si efectivamente estos valores queremos que se puedan configurar desde
		// la linea de comandos.	
	

	        std::ifstream filer;
		filer.open(data.c_str(), std::ifstream::in);
         
	        ArbolData arbol_data;

		float radius = 4.0f;
		float dist = 5.7f;
		arbol_data.grilla = new Grillado(m, n , z, radius, dist);
                // Maximun gyration radius and maximun CA-CA distance. 
                // Both equations constructed from database analisys.
		arbol_data.rgmax =  2.72 * cubic_root(float(Nres)) + 5.0;
	        arbol_data.dmax2 =  (8.0 * cubic_root(float(Nres))+25.0) * (8.0 * cubic_root(float(Nres))+25.0);

                // Number of residues
		arbol_data.nres = Nres;
		// Initialize the atm matrix.
		arbol_data.atm = new ATOM[(arbol_data.nres)*3];
	        // Initialize the chain counter
	        arbol_data.cont = 0;
	        arbol_data.hubo_algun_exito = false;

                // Fill r[][][] with the minimun squared distance between atoms
	        setr(RN,RCa,RC,Scal_1_4,Scal_1_5);
                //Default compressed output
	        arbol_data.xfp = xdrfile_open("traj.xtc","w");
	        
	        readdata(filer, arbol_data.cosfi, arbol_data.sinfi, arbol_data.cossi, arbol_data.sinsi);

		arbol_data.ndat = arbol_data.cossi.size(); // Se podria eliminar ndat por completo. Tener en cuenta a futuro.
                printf("Number of fi-si combinations in file=%i\n",arbol_data.ndat);

	        generar_arbol(&arbol_data);
	        printf("Number of chains generated=%li\n",arbol_data.cont);        
		
		// Se libera la memoria de la matriz atm.
		delete [] arbol_data.atm;
		delete arbol_data.grilla;
		filer.close();
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
	Residuo residuo;// Va a ser el residuo que agregue semilla en cada iteracion y al terminar del ciclo
			// se usa para sacar el residuo del grillado.
    //inicializar_arbol(arbol_data);  


	i = 0;
	while ( i < arbol_data->ndat && !arbol_data->hubo_algun_exito )
	{
		//printf("cleardata\n");
		clearatm( arbol_data->atm, arbol_data->nres);
		semilla(arbol_data,R_inicial, residuo);
		
		generar_nivel_intermedio(2, R_inicial, i, arbol_data);
		sacar_residuo(arbol_data, residuo);
		i++;
	}
}

// se asume que nivel > 1
void generar_nivel_intermedio(unsigned int nivel, float R_inicial[16], unsigned int indice_nivel_anterior, ArbolData* arbol_data)
{
	float R_local[16];
	FilterResultType resultado;
	bool exito = false; // solo aplicable si somos anteultimo nivel
	unsigned int i;
	assert(nivel > 1);  // pre condicion
	Residuo residuo;

	i = 0;
	while (i < arbol_data->ndat && !exito)
	{
		
		copymat(R_local, R_inicial);

		resultado = poneres(R_local,
                        arbol_data->cossi[indice_nivel_anterior],   
                        arbol_data->sinsi[indice_nivel_anterior],   
                        arbol_data->cosfi[i],                       
                        arbol_data->sinfi[i],                       
                        arbol_data->atm,
                        nivel,
			arbol_data,
			arbol_data->dmax2,
			residuo);

		
		if ( resultado == FILTER_OK )
		{
			if (nivel < arbol_data->nres)
			{              
				generar_nivel_intermedio(nivel+1, R_local, i, arbol_data);
			}
			else
			{
				exito = procesar_ultimo_nivel(arbol_data);
			}
			sacar_residuo(arbol_data, residuo);
		}
		
		i++;
	}
}


#ifdef COMBINATIONS_DEBUG
// En el modo DEBUG se deshabilitan los chequeos.
static bool procesar_ultimo_nivel(ArbolData* arbol_data) { 
	writextc(arbol_data->xfp, 
                 arbol_data->nres, 
                 arbol_data->cont, 
                 arbol_data->atm);

	arbol_data->cont++;
	return false;
}
#else
static bool procesar_ultimo_nivel(ArbolData* arbol_data)
{
	bool exito = false;
    
	if ( filtros_ultimo_nivel(arbol_data) == FILTER_OK )        
	{
		writextc(arbol_data->xfp, 
                 arbol_data->nres, 
                 arbol_data->cont, 
                 arbol_data->atm);

		arbol_data->cont++;
		arbol_data->hubo_algun_exito = exito = true;
	}
	
	return exito;
}
#endif
// Devuelve 1 si el volumen indicado por el grillado se encuentra en el rango aceptable, 0 si no.


static FilterResultType volumen_en_rango(const ArbolData * arbol_data) {
	const float volumen_max_aa = pendiente_empirica * float(arbol_data->nres) + cota_maxima_volumen;
#ifdef VERBOSE
        printf("Maximun Volume allowed per a.a  =%f.   Volumen in this chain=%f\n",volumen_max_aa,float(arbol_data->grilla->obtener_vol_parcial())/float(arbol_data->nres) );
#endif
  FilterResultType res = FILTER_FAIL;
  if(in_range(float(arbol_data->grilla->obtener_vol_parcial())/float(arbol_data->nres), volumen_min_aa , volumen_max_aa)) {
    res = FILTER_OK;
  }
	return res;
}



// Devuelve 1 si arbol_data pasa todos los filtros.
static FilterResultType filtros_ultimo_nivel(const ArbolData * arbol_data) {
  FilterResultType res = FILTER_FAIL;
	 	if(calcRdG(arbol_data->atm, arbol_data->nres, arbol_data->rgmax) == FILTER_OK 
		&& volumen_en_rango(arbol_data) == FILTER_OK) {
		  res = FILTER_OK;
    }
		return res;
}

void show_usage()
{
    std::cerr << "Invalid arguments." << std::endl;
}

