#ifndef PETU_H
#define PETU_H
#include <vector>

#include "grillado.h"
#include "prot-filer/definitions.h"

enum FilterResultType  {FILTER_FAIL, FILTER_OK};

#include "prot-filer/format_filer.h"
#include "prot-filer/fragments_filer.h"
#include "prot-filer/cached_reader.h"
#include "prot-filer/angles.h"
#include "prot-filer/backbones_utils.h"

// Datos a compartir por todos los niveles:
using namespace prot_filer;
typedef CachedReader<FullCache, SimpleAnglesReader, AnglesData> FullCachedAnglesSeqReader;
typedef BasicProtein Atoms;

struct TreeData
{
    unsigned int nres;
    float rgmax, dmax2;
    std::vector<float> cosfi, cossi, sinfi, sinsi;    // constantes
    Atoms atm;       // estructura parcial
    long int cont;              // cantidad de estructuras exitosas hasta el momento
    bool hubo_algun_exito;      // si encendido, dice que hubo al menos una rama que llego al final
    Grillado* grilla;       // Utilizamos el grillado para aproximar el volumen parcial
    AnglesData* angles_data; // Used only when writing compressed data.
    AnglesMapping* angles_mapping;
    FragmentIds fragment_ids;
    const std::string output_file;

    TreeData(int nRes, Grillado* grillado) :
        nres(nRes),
        // Maximun gyration radius and maximun CA-CA distance.
        // Both equations constructed from database analisys.
        rgmax(2.72 * mili::cubic_root(nres) + 5.0),
        dmax2(mili::square(8.0 * mili::cubic_root(nres) + 25.0)),
        atm(Atoms(nres * 3)),
        cont(0),
        hubo_algun_exito(false),
        grilla(grillado),
        angles_mapping(new AnglesMapping(nres)),
        output_file("traj.xtc")
    {}
    ~TreeData()
    {
        delete grilla;
        delete angles_data;
        delete angles_mapping;
    }
};

struct Residuo
{
    Residuo(const esferaId& param_at1, const esferaId& param_at2, const esferaId& param_at3)
        : at1(param_at1),
          at2(param_at2),
          at3(param_at3)
    {};
    Residuo() {};
    esferaId at1;
    esferaId at2;
    esferaId at3;
};

/* a esto se le pone extern, para decir que quien importe el .h va a ser usuario de estas variables */
extern float r[3][3][3];

#endif
