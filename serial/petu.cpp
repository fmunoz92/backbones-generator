#include <cassert>
#include <iostream>
#include <fstream>
#include "getopt_pp.h"

#include <mili/mili.h>

#include "poneres.h"
#include "readdata.h"
#include "calcrdg.h"
#include "setr.h"
#include "imprime.h"
#include "clearatm.h"
#include "copyatm.h"

const float cota_maxima_volumen = 177.65f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.
const float pendiente_empirica = -0.0882f ; // Pendiente obtenida a partir de las pruebas de un set de datos en Grillado.
const float volumen_min_aa = 110.0f;
//const float volumen_min_aa = 0.0f; // Volumen obtenido tambien a partir de las pruebas de un set de datos en Grillado.

TreeData::TreeData(int nRes, FormatFiler* filer, Grillado* grillado) :
    nres(nRes),
    // Maximun gyration radius and maximun CA-CA distance.
    // Both equations constructed from database analisys.
    rgmax(2.72 * mili::cubic_root(nres) + 5.0),
    dmax2(mili::square(8.0 * mili::cubic_root(float(nres)) + 25.0)),
    atm(new ATOM[nres * 3]),
    cont(0),
    hubo_algun_exito(false),
    grilla(grillado),
    filer(filer),
    angles_mapping(new AnglesMapping(nres))
{
    //Default compressed output
    filer->open_write("traj.xtc");
}

TreeData::~TreeData()
{
    filer->close();
    delete [] atm;
    delete grilla;
    delete filer;
    delete angles_data;
    delete angles_mapping;
}

class TreeGenerator
{
private:
    TreeData* tree_data;
    AnglesDatabase* residue_chains;
    const unsigned int size;
    const bool with_residues_input;
public:
    void generate();
    TreeGenerator(TreeData*, AnglesDatabase* residue_chain_database);
    ~TreeGenerator()
    {
        delete residue_chains;
    }

private:
    void semilla(float* R, Residuo& residuo);
    void generar_nivel_intermedio(unsigned int nivel, float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
    FilterResultType filtros_ultimo_nivel();
    FilterResultType volumen_en_rango();
    void sacar_residuo(const Residuo& residuo)
    {
        //tree_data->grilla->sacar_esfera(residuo.at1);
        tree_data->grilla->sacar_esfera(residuo.at2);
        //tree_data->grilla->sacar_esfera(residuo.at3);
    }
    void sacar_residuos(const vector<Residuo>& residuos)
    {
        for (unsigned int i = 0; i < residuos.size(); ++i)
        {
            sacar_residuo(residuos[i]);
        }
    }
    // Devuelve 1 si value esta entre min y max.
    // Se separo en una funcion distinta para mejorar el rendimiento de volumen_en_rango.
    bool in_range(float value, float min, float max)
    {
        return min < mili::bchain(value) < max;
    }
};

static void show_usage();

int main(int argc, char** argv)
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
    std::string residues_input; //Name of the residue chains file.

    size_t m; // They indicate the size of each dimention of the grid.
    size_t n;
    size_t z;

    FilerFactory* filer_factory = FilerFactory::get_instance();

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
                >> Option('w', "write_format", write_format, string("xtc"))
                >> Option("chains_input", residues_input);
        //TODO:
        // Hay que decidir si efectivamente estos valores queremos que se puedan configurar desde
        // la linea de comandos.
        std::ifstream filer;
        filer.open(data.c_str());
        FormatFiler* outFiler = filer_factory->create(write_format);
        const float radius = 4.0f;
        const float dist = 5.7f;
        Grillado* grilla = new Grillado(m, n , z, radius, dist);
        TreeData tree_data(Nres, outFiler, grilla);

        // Fill r[][][] with the minimun squared distance between atoms
        setr(RN, RCa, RC, Scal_1_4, Scal_1_5);

        readdata(filer, tree_data.cosfi, tree_data.sinfi, tree_data.cossi, tree_data.sinsi, tree_data.angles_mapping);
        tree_data.angles_data = new AnglesData(tree_data.nres, *tree_data.angles_mapping);

        cout << "Number of fi-si combinations in file=" << tree_data.cossi.size() << endl;

        AnglesDatabase* residue_chain_database = NULL;

        if (!residues_input.empty())
        {
            residue_chain_database = read_chains(residues_input, "compressed", "full_cache");
        }

        TreeGenerator(&tree_data, residue_chain_database).generate();
        cout << "Number of chains generated=" << tree_data.cont << endl;

        FilerFactory::destroy_instance();
        ProteinCacheFactory::destroy_instance();
        AnglesCacheFactory::destroy_instance();
        //residue_chain_database is deleted by TreeGenerator
        filer.close();
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

TreeGenerator::TreeGenerator(TreeData* data, AnglesDatabase* residue_chain_database) :
    tree_data(data),
    residue_chains(residue_chain_database),
    size(residue_chain_database != NULL ? (residue_chain_database->size() - 1) : tree_data->cossi.size()),
    with_residues_input(residue_chain_database != NULL)
{};

void TreeGenerator::semilla(float* R, Residuo& residuo)
{
    ATOM* atm = tree_data->atm;
    backbones_utils::semilla(atm, R);
    ATOM* seed = tree_data->angles_data->seed;
    
    copyatm(atm[0], seed[0]);
    copyatm(atm[1], seed[1]);
    copyatm(atm[2], seed[2]);

    residuo.at2 = tree_data->grilla->agregar_esfera(tree_data->atm[1].x, tree_data->atm[1].y, tree_data->atm[1].z);
}

void TreeGenerator::generate()
{
    float R_inicial[16];
    unsigned int i;
    Residuo residuo;// Va a ser el residuo que agregue semilla en cada iteracion y al terminar del ciclo
    // se usa para sacar el residuo del grillado.

    i = 0;
    while (i < size && !tree_data->hubo_algun_exito)
    {
        clearatm(tree_data->atm, tree_data->nres);
        semilla(R_inicial, residuo);

        generar_nivel_intermedio(2, R_inicial, i);
        sacar_residuo(residuo);
        i++;
    }
}

// se asume que nivel > 1
void TreeGenerator::generar_nivel_intermedio(const unsigned int nivel, float R_inicial[16], const unsigned int indice_nivel_anterior)
{
    float R_local[16];
    FilterResultType resultado;
    bool exito = false; // solo aplicable si somos anteultimo nivel
    assert(nivel > 1);  // pre condicion
    vector<Residuo> residuos;
    IncompleteAnglesData data(1);

    unsigned int i = 0;
    while (i < size && !exito)
    {
        backbones_utils::copymat(R_local, R_inicial);

        if (!with_residues_input)
        {
            data.angles[0].si = indice_nivel_anterior;
            data.angles[0].fi = i;
        }

        resultado = addNRes(R_local,
                            nivel,
                            tree_data,
                            residuos,
                            with_residues_input ? (*residue_chains)[i] : data);

        if (resultado == FILTER_OK)
        {
            if (nivel < tree_data->nres)
            {
                generar_nivel_intermedio(nivel + residuos.size(), R_local, i);
            }
            else
            {
                exito = procesar_ultimo_nivel();
            }
            sacar_residuos(residuos);
        }
        residuos.clear();
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
    if (calcRdG(tree_data->atm, tree_data->nres, tree_data->rgmax) == FILTER_OK
            && volumen_en_rango() == FILTER_OK)
    {
        res = FILTER_OK;
    }
    return res;
}

FilterResultType TreeGenerator::volumen_en_rango()
{
    const float volumen_max_aa = pendiente_empirica * float(tree_data->nres) + cota_maxima_volumen;
#ifdef VERBOSE
    printf("Maximun Volume allowed per a.a  =%f.   Volumen in this chain=%f\n", volumen_max_aa, float(tree_data->grilla->obtener_vol_parcial()) / float(tree_data->nres));
#endif
    FilterResultType res = FILTER_FAIL;
    if (in_range(float(tree_data->grilla->obtener_vol_parcial()) / float(tree_data->nres), volumen_min_aa , volumen_max_aa))
    {
        res = FILTER_OK;
    }
    return res;
}
