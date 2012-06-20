#ifndef TREE_DATA_H
#define TREE_DATA_H

#include <string>
#include <memory>

#include "prot-filer/format_filer.h"
#include "prot-filer/cached_reader.h"

#include "grillado.h"

enum FilterResultType  {FILTER_FAIL, FILTER_OK};

typedef prot_filer::CachedReader<prot_filer::FullCache, prot_filer::SimpleAnglesReader, prot_filer::AnglesData> FullCachedAnglesSeqReader;
typedef prot_filer::BasicProtein Atoms;

// Datos a compartir por todos los niveles:
struct TreeData
{
    unsigned int nres;
    float rgmax, dmax2;
    std::vector<float> cosfi, cossi, sinfi, sinsi;    // constantes
    Atoms atm;       // estructura parcial
    long int cont;              // cantidad de estructuras exitosas hasta el momento
    bool hubo_algun_exito;      // si encendido, dice que hubo al menos una rama que llego al final
    const std::auto_ptr<Grillado>& grilla;       // Utilizamos el grillado para aproximar el volumen parcial
    const std::auto_ptr<prot_filer::AnglesMapping> angles_mapping;
    const std::auto_ptr<prot_filer::AnglesData> angles_data; // Used only when writing compressed data.
    prot_filer::FragmentIds fragment_ids;
    const std::string output_file;

    TreeData(int nRes, std::auto_ptr<Grillado>& grillado, std::istream& input_file);

private:
    void readdata(std::istream& filer);
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

#endif
