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

const float cota_maxima_volumen = 177.65f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.
const float pendiente_empirica = -0.0882f ; // Pendiente obtenida a partir de las pruebas de un set de datos en Grillado.
const float volumen_min_aa = 110.0f;
//const float volumen_min_aa = 0.0f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.

class TreeGenerator
{
    private:
        TreeData* tree_data;
    public:
        void generate();
        TreeGenerator(TreeData*); 
    private:
        void generar_nivel_intermedio(unsigned int nivel, float R_inicial[16], unsigned int indice_nivel_anterior);
        bool procesar_ultimo_nivel();
        FilterResultType filtros_ultimo_nivel();
        FilterResultType volumen_en_rango();
        void sacar_residuo(Residuo& residuo) 
        {
        	//tree_data->grilla->sacar_esfera(residuo.at1);
        	tree_data->grilla->sacar_esfera(residuo.at2);
        	//tree_data->grilla->sacar_esfera(residuo.at3);
        }
        // Devuelve 1 si value esta entre min y max.
        // Se separo en una funcion distinta para mejorar el rendimiento de volumen_en_rango.
        bool in_range(float value, float min, float max) 
        {
            return min < bchain(value) < max;
        }
};

inline float cubic_root(float a) {
	return (pow(a, 1.0f/3.0f));
}

static void show_usage();

int main(int argc, char **argv) 
{
    using namespace GetOpt;

    int     Nres;   // Number of amino acids in the chains to build.

    float   RN,         // Radius of the Nitrogen atom.
            RCa,        // Radius of the Carbon atom.
            RC,         // Radius of the Carbon atom.
            Scal_1_4,   // Scaling factor for the radii. Used to check for 1-4 clashes.
            Scal_1_5;   // Scaling factor for the radii. Used to check for 1-5 clashes.

    std::string data;   // Name of the input file.
    std::string write_format;
   
    size_t m; // They indicate the size of each dimention of the grid.
    size_t n;
    size_t z;
    
    FilerFactory * filer_factory = FilerFactory::get_instance();

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
	    >> Option('Z', "depth", z, static_cast<size_t>(100))
	    >> Option('w', "write_format", write_format, string("xtc"));
		// Hay que decidir si efectivamente estos valores queremos que se puedan configurar desde
		// la linea de comandos.	
	
      std::ifstream filer;
		filer.open(data.c_str(), std::ifstream::in);            
                                 
	   TreeData tree_data;
	   tree_data.filer = filer_factory->create(write_format);


		float radius = 4.0f;
		float dist = 5.7f;
		tree_data.grilla = new Grillado(m, n , z, radius, dist);
                // Maximun gyration radius and maximun CA-CA distance. 
                // Both equations constructed from database analisys.
		tree_data.rgmax =  2.72 * cubic_root(float(Nres)) + 5.0;
	    tree_data.dmax2 =  (8.0 * cubic_root(float(Nres))+25.0) * (8.0 * cubic_root(float(Nres))+25.0);
      
        // Number of residues
		tree_data.nres = Nres;
		// Initialize the atm matrix.
		tree_data.atm = new ATOM[(tree_data.nres)*3];
	    // Initialize the chain counter
	    tree_data.cont = 0;
	    tree_data.hubo_algun_exito = false;

            // Fill r[][][] with the minimun squared distance between atoms
	    setr(RN,RCa,RC,Scal_1_4,Scal_1_5);
            //Default compressed output
	    tree_data.filer->open_write("traj.xtc");
	    //tree_data.filer->open_read("traj_compressed.xtc"); delete [] tree_data.filer->read();
	    tree_data.angles_mapping = new AnglesMapping( tree_data.nres);
	    
	    readdata(filer, tree_data.cosfi, tree_data.sinfi, tree_data.cossi, tree_data.sinsi, tree_data.angles_mapping);
	    tree_data.angles_data = new AnglesData(tree_data.nres, *tree_data.angles_mapping);

		tree_data.ndat = tree_data.cossi.size(); // Se podria eliminar ndat por completo. Tener en cuenta a futuro.
        printf("Number of fi-si combinations in file=%i\n",tree_data.ndat);

	    TreeGenerator(&tree_data).generate();
	    printf("Number of chains generated=%li\n", tree_data.cont);        
		
		// Se libera la memoria de la matriz atm.
		tree_data.filer->close();
		delete [] tree_data.atm;
		delete tree_data.grilla;
		delete tree_data.filer;
		delete tree_data.angles_data;
		delete tree_data.angles_mapping;
		filer.close();
		delete filer_factory;
	    return EXIT_SUCCESS;
    }
    else
    {
        show_usage();
        return EXIT_FAILURE;
    }
}

void show_usage()
{
    std::cerr << "Invalid arguments." << std::endl;
}

TreeGenerator::TreeGenerator(TreeData* data) : tree_data(data) {};

void TreeGenerator::generate()
{
	float R_inicial[16];
	unsigned int i;
	Residuo residuo;// Va a ser el residuo que agregue semilla en cada iteracion y al terminar del ciclo
			// se usa para sacar el residuo del grillado.
    //inicializar_arbol(tree_data);  

	i = 0;
	while ( i < tree_data->ndat && !tree_data->hubo_algun_exito )
	{
		//printf("cleardata\n");
		clearatm(tree_data->atm, tree_data->nres);
		semilla(tree_data,R_inicial, residuo);
		
		generar_nivel_intermedio(2, R_inicial, i);
		sacar_residuo(residuo);
		i++;
	}
}

// se asume que nivel > 1
void TreeGenerator::generar_nivel_intermedio(unsigned int nivel, float R_inicial[16], unsigned int indice_nivel_anterior)
{
	float R_local[16];
	FilterResultType resultado;
	bool exito = false; // solo aplicable si somos anteultimo nivel
	unsigned int i;
	assert(nivel > 1);  // pre condicion
	Residuo residuo;

	i = 0;
	while (i < tree_data->ndat && !exito)
	{
		
		copymat(R_local, R_inicial);

		resultado = poneres(R_local,
                        tree_data->cossi[indice_nivel_anterior],   
                        tree_data->sinsi[indice_nivel_anterior],   
                        tree_data->cosfi[i],                       
                        tree_data->sinfi[i],                       
                        tree_data->atm,
                        nivel,
                        tree_data,
                        tree_data->dmax2,
                        residuo,
                        indice_nivel_anterior,
                        i);

		
		if ( resultado == FILTER_OK )
		{
			if (nivel < tree_data->nres)
			{              
				generar_nivel_intermedio(nivel+1, R_local, i);
			}
			else
			{
				exito = procesar_ultimo_nivel();
			}
			sacar_residuo(residuo);
		}
		
		i++;
	}
}

// True si el volumen indicado por el grillado se encuentra en el rango aceptable
bool TreeGenerator::procesar_ultimo_nivel() 
{ 
#ifdef COMBINATIONS_DEBUG 
// En el modo DEBUG se deshabilitan los chequeos.
	tree_data->filer->write(tree_data->atm, *tree_data->angles_data);
	tree_data->cont++;
	return false;
#else
	bool exito = false;
    
	if (filtros_ultimo_nivel() == FILTER_OK)
    {
        tree_data->filer->write(tree_data->atm, *tree_data->angles_data);
		tree_data->cont++;
		tree_data->hubo_algun_exito = exito = true;
	}
	return exito;
#endif
}

// Devuelve 1 si tree_data pasa todos los filtros.
FilterResultType TreeGenerator::filtros_ultimo_nivel() 
{
    FilterResultType res = FILTER_FAIL;
	if(calcRdG(tree_data->atm, tree_data->nres, tree_data->rgmax) == FILTER_OK 
		&& volumen_en_rango() == FILTER_OK) {
        res = FILTER_OK;
    }
	return res;
}

FilterResultType TreeGenerator::volumen_en_rango() 
{
	const float volumen_max_aa = pendiente_empirica * float(tree_data->nres) + cota_maxima_volumen;
#ifdef VERBOSE
    printf("Maximun Volume allowed per a.a  =%f.   Volumen in this chain=%f\n",volumen_max_aa,float(tree_data->grilla->obtener_vol_parcial())/float(tree_data->nres) );
#endif
    FilterResultType res = FILTER_FAIL;
    if(in_range(float(tree_data->grilla->obtener_vol_parcial())/float(tree_data->nres), volumen_min_aa , volumen_max_aa)) {
        res = FILTER_OK;
    }
	return res;
}
