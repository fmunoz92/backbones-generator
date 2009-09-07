

#include "petu.h"
#include <assert.h>
#include <iostream>
#include <stdlib.h>

#include <string>




#include "getopt_pp_standalone.h"




// La estructura del algoritmo esta dividida en tres funciones, que se llaman en este orden:
static void generar_arbol(ArbolData* arbol_data);

static void generar_nivel_intermedio(unsigned int nivel, 
                                     float R_inicial[16], 
                                     unsigned int indice_nivel_anterior, 
                                     ArbolData* arbol_data);

static bool procesar_ultimo_nivel(ArbolData* arbol_data); //bool void...?
static int volumen_en_rango(ArbolData * arbol_data);
static int filtros_ultimo_nivel(ArbolData * arbol_data) ;
static void show_usage();

const float cota_maxima_volumen = 177.65f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.
const float pendiente_empirica = -0.0882f ; // Pendiente obtenida a partir de las pruebas de un set de datos en Grillado.
const float volumen_min_aa = 110.0f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.
	

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
            RgMax,      // Radio de giro maximo
            DMax;       // Distancia maxima entre átomos
	std::string data; // Nombre de archivo de input.
	
	size_t m; // Indican el tamano de cada dimension del grillado.
	size_t n;
	size_t z;


		 GetOpt_pp ops(argc, argv);
		 
		 
		 
		 
    if (ops >> Option('r', "Nres", Nres))
    {
	float default_DMax = 50.0 + 0.1 * float(Nres);
        ops
            //>> Option('R', "RgMax", RgMax,2.72 * (pow(Nres,0.33333)) + 5)
            >> Option('n', "Rn", RN, 1.5f)
            >> Option('a', "Rca", RCa, 1.7f)
            >> Option('c', "Rc", RC, 1.6f)
	    >> Option('d', "DMax", DMax, default_DMax ) 
            >> Option('s', "Scal_1_4", Scal_1_4, 0.85f)
            >> Option('l', "Scal_1_5", Scal_1_5, 1.0f)
	    >> Option('i', "input_file", data, "ramachandran.dat")
	    >> Option('N', "rows", n, static_cast<size_t>(100))
            >> Option('M', "cols", m, static_cast<size_t>(100))
	    >> Option('Z', "depth", z, static_cast<size_t>(100));
		// Hay que decidir si efectivamente estos valores queremos que se puedan configurar desde
		// la linea de comandos.	
	
		RgMax = 2.72 * cubic_root(float(Nres)) + 5;
	        std::ifstream filer;
		filer.open(data.c_str(), std::ifstream::in);
         
	        ArbolData arbol_data;
		float radius = 1.0f;
		float dist = 1.5f;
		arbol_data.grilla = new Grillado(m, n , z, radius, dist);
		arbol_data.nres = Nres;
		// Se pide inicializa la matriz atm.
		arbol_data.atm = new ATOM[(arbol_data.nres)*3+1];
		// ATENCION: En clearatm se indexa la matriz a partir del 1 y la condicion del for es <=
		// con lo cual el indice alcanza el valor nres*3. Nota a futuro: revisar todos los indices.

	        
	        arbol_data.cont = 0;
	        arbol_data.hubo_algun_exito = false;

	        arbol_data.rgmax= RgMax;
	        arbol_data.dmax2= DMax*DMax;
	        setr(RN,RCa,RC,Scal_1_4,Scal_1_5);

	        arbol_data.xfp = xdrfile_open("traj.xtc","w");
	        
	        readdata(arbol_data.ndat, filer, arbol_data.cosfi, arbol_data.sinfi, arbol_data.cossi, arbol_data.sinsi);

		arbol_data.ndat = arbol_data.cossi.size()-1; // Se podria eliminar ndat por completo. Tener en cuenta a futuro.
		// Se usa size()-1 porque en readdata se agrega un elemento inicial para respetar la convencion de usar los arrays desde 1.
	        generar_arbol(&arbol_data);
	        printf("%li\n",arbol_data.cont);        
		
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
	i = 1;
	while ( i <= arbol_data->ndat && !arbol_data->hubo_algun_exito )
	{
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
	int resultado;
	bool exito = false; // solo aplicable si somos anteultimo nivel
	unsigned int i;
	assert(nivel > 1);  // pre condicion
	Residuo residuo;
	
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
                        nivel,
			arbol_data,
			arbol_data->dmax2,
			residuo);

		
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
	//arbol_data->hubo_algun_exito  = true;
	return true;
}
#else
static bool procesar_ultimo_nivel(ArbolData* arbol_data)
{
	bool exito = false;
    
	if ( filtros_ultimo_nivel(arbol_data) )        
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
// Se usa int porque el resto de las funciones booleanas estan implementadas en C y por lo tanto devuelven int.


static int volumen_en_rango(ArbolData * arbol_data) {
	float volumen_max_aa = pendiente_empirica * float(arbol_data->nres) + cota_maxima_volumen;
	return in_range(float(arbol_data->grilla->obtener_vol_parcial()), volumen_min_aa , volumen_max_aa);
}



// Devuelve 1 si arbol_data pasa todos los filtros.
static int filtros_ultimo_nivel(ArbolData * arbol_data) {
	return 	calcRdG(arbol_data->atm, arbol_data->nres, arbol_data->rgmax) == BIEN 
		&& volumen_en_rango(arbol_data);
}

void show_usage()
{
    std::cerr << "Invalid arguments." << std::endl;
}

