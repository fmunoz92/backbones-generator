#ifndef TREE_DATA_H
#define TREE_DATA_H

#include <string>
#include <memory>

#include "prot-filer/format_filer.h"
#include "prot-filer/cached_reader.h"

#include "grillado.h"

enum FilterResultType //move to TreeFilters
{
    FILTER_FAIL,
    FILTER_OK
};

typedef prot_filer::CachedReader<prot_filer::FullCache, prot_filer::SimpleAnglesReader, prot_filer::AnglesData> FullCachedAnglesSeqReader;
typedef prot_filer::BasicProtein Atoms;

// Datos a compartir por todos los niveles:
class TreeData
{
public:
    //tomar referencias a prot_filer::AnglesMapping, prot_filer::AnglesData y Grillado en vez de crearlos aca
    TreeData(int nRes, size_t cols, size_t rows, size_t depth, std::istream& input_file, std::string& output_file);

    const unsigned int nres;

    const float rgmax;
    const float dmax2;

    std::vector<float> cosfi;
    std::vector<float> cossi;
    std::vector<float> sinfi;
    std::vector<float> sinsi;

    Atoms atm; // estructura parcial

    long int cont;         // cantidad de estructuras exitosas hasta el momento
    bool hubo_algun_exito; // si encendido, dice que hubo al menos una rama que llego al final

    const std::auto_ptr<Grillado>                  grilla;       // Utilizamos el grillado para aproximar el volumen parcial
    const std::auto_ptr<prot_filer::AnglesMapping> angles_mapping;
    const std::auto_ptr<prot_filer::AnglesData>    angles_data; // Used only when writing compressed data.

    prot_filer::FragmentIds fragment_ids;

    const std::string output_file;

private:
    void readdata(std::istream& filer);
};

struct Residuo
{
    Residuo(const Residuo& r)
        : at1(r.at1),
          at2(r.at2),
          at3(r.at3)
    {}

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

#endif
