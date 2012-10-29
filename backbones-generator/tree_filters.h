#ifndef TREE_FILTERS_H
#define TREE_FILTERS_H

#include "backbones-generator/tree_data.h"

class TreeFilters
{
public:

    enum FilterResultType
    {
        FILTER_FAIL,
        FILTER_OK
    };

    void setr(float rn, float rca, float rc, float scal_1_4, float scal_1_5);

    //This function calculates the gyration radius of the CA atoms
    //It was meant to filter long chains. Now it is replaced by the volume filter
    FilterResultType calcRdG(const Atoms& patm, unsigned int nres, float rgmax) const;

    //This function filter long chains.
    //Calculates the CA-CA distances of the recently added CA with the rest.
    //Then compares this distances with dmax2, wich cames from the analisys of protein
    //structures experimetally determined
    FilterResultType islong(const Atoms& patm, unsigned int at, float dmax2) const;

    FilterResultType isclash(const Atoms& patm, unsigned int at) const;

    static inline FilterResultType volumeInRange(unsigned int nres, Volume partialVolume);
private:
    bool isParticularClash(const Atoms& patm, unsigned int x, unsigned int y, unsigned int z) const;
    //This function detect colitions between atoms
    //The matrix r have the reference radius (actually the sum of the squares of the raduis)
    //The last index of r is 0 for 1-4 clashes, 1 for 1-5 and 2 for the rest
    static inline float distance(const Atoms& patm, unsigned int at, unsigned int i);

    float r[3][3][3];
};

class ClashFilter
{
public:
    ClashFilter(const TreeData& tree_data, const TreeFilters& tree_filters);
    bool operator()(unsigned int index, const Atoms& patm, unsigned int at) const;
private:
    const TreeData& tree_data;
    const TreeFilters& tree_filters;
};

inline float TreeFilters::distance(const Atoms& patm, const unsigned int at, const unsigned int i)
{
    const float dx2 = square(patm[at].x - patm[i].x);
    const float dy2 = square(patm[at].y - patm[i].y);
    const float dz2 = square(patm[at].z - patm[i].z);

    return dx2 + dy2 + dz2;
}

inline TreeFilters::FilterResultType TreeFilters::volumeInRange(unsigned int nres, Volume partialVolume)
{
    // Valores obtenidos a partir de pruebas de un set de datos en Grillado.
    static const float cota_maxima_volumen = 177.65f;
    static const float pendiente_empirica = -0.0882f;
    static const float volumen_min_aa = 110.0f;

    const float volumen_max_aa = pendiente_empirica * float(nres) + cota_maxima_volumen;
    const float chain_volumen = float(partialVolume) / float(nres);

#ifdef VERBOSE
    std::cout << "Maximun Volume allowed per a.a  =" << volumen_max_aa << ". Volumen in this chain=" << chain_volumen << endl;
#endif

    return mili::in_range(chain_volumen, volumen_min_aa, volumen_max_aa) ? FILTER_OK : FILTER_FAIL;
}

#endif
