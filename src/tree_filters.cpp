#include <cmath>
#include <iostream>
#include "backbones-generator/tree_filters.h"

using mili::square;
using mili::in_range;

void TreeFilters::setr(float rn, float rca, float rc, float scal_1_4, float scal_1_5)
{
    /* La matriz de distancias minimas al cuadrado*/
    float radio[3];

    radio[0] = rn;
    radio[1] = rca;
    radio[2] = rc;

    for (unsigned int i = 0; i <= 2; i++)
        for (unsigned int j = 0; j <= 2; j++)
        {
            r[i][j][0] = square(radio[i] + radio[j]) * scal_1_4 * scal_1_4;
            r[i][j][1] = square(radio[i] + radio[j]) * scal_1_5 * scal_1_5;
            r[i][j][2] = square(radio[i] + radio[j]);
        }
}

TreeFilters::FilterResultType TreeFilters::calcRdG(const Atoms& patm, unsigned int nres, float rgmax) const
{
    float Rcm = 0;
    float xcm = 0;
    float ycm = 0;
    float zcm = 0;

    for (unsigned int i = 1; i <= (nres * 3) - 2; i += 3)
    {
        xcm += patm[i].x;
        ycm += patm[i].y;
        zcm += patm[i].z;
    }

    xcm /= nres;
    ycm /= nres;
    zcm /= nres;

    for (unsigned int i = 1; i <= (nres * 3) - 2; i += 3)
    {
        const float dx2 = square(patm[i].x - xcm);
        const float dy2 = square(patm[i].y - ycm);
        const float dz2 = square(patm[i].z - zcm);

        Rcm += (dx2 + dy2 + dz2);
    }

#ifdef VERBOSE
    std::cout << "Raduis of gyration = " << sqrt(Rcm / nres) << " Maximun allowed = " << rgmax << std::endl;
#endif

    return (sqrt(Rcm / nres) > rgmax) ? FILTER_FAIL : FILTER_OK;
}

TreeFilters::FilterResultType TreeFilters::islong(const Atoms& patm, unsigned int at, float dmax2) const
{
    // Until we reach residue #5, the check for long chain is meaningless
    // at=13 is the CA of residue #5
    if (at < 13)
        return FILTER_OK;

    // at-12 is the CA atom that is four residues down the chain
    // at=1 is the first CA in the chain
    for (int i = at - 12; i >= 1; i -= 3)
    {
        const float dx2 = square(patm[at].x - patm[i].x);
        const float dy2 = square(patm[at].y - patm[i].y);
        const float dz2 = square(patm[at].z - patm[i].z);
        const float d2  = dx2 + dy2 + dz2;

        if (d2 > dmax2)
        {
#ifdef VERBOSE
            std::cout << "Chain length = " << sqrt(d2) << ", while Dmax is = " << sqrt(dmax2) << std::endl;
#endif
            return FILTER_FAIL;
        }
    }
    return FILTER_OK;
}

TreeFilters::FilterResultType TreeFilters::isclash(const Atoms& patm, unsigned int at) const
{
    //This is to check for the so-called 1-4 clashes,
    //i.e. a clash between atom at position i with atom at position i+3
    int i  = at - 3;
    float d2 = distance(patm, at, i);

    if (d2 < r[patm[at].vdw][patm[i].vdw][0])
    {
#ifdef VERBOSE
        std::cout << "Clash 1-4 between atmom = " << at << " and atom = " << i << " distancia= " <<  sqrt(d2) << std::endl;
#endif
        return FILTER_FAIL;
    }

    //If this is atom #3 and passes the check for 1-4 clashes, then there is nothing else to check.
    if (at == 3)
        return FILTER_OK;

    //This is to check for the so-called 1-5 clashes;
    //i.e. a clash between atom at position i with atom at position i+4
    i = at - 4;
    d2 = distance(patm, at, i);
    if (d2 < r[patm[at].vdw][patm[i].vdw][1])
    {
#ifdef VERBOSE
        std::cout << "Clash 1-5 between atmom = " << at << " and atom = " << i << " distancia= " <<  sqrt(d2) << std::endl;
#endif
        return FILTER_FAIL;
    }

    //If this is atom #4 and passes check for 1-4 and 1-5 clashes, then there is nothing else to check.
    if (at == 4)
        return FILTER_OK;

    //The rest of the clashes until the end of the chain.
    for (i = at - 5; i >= 0; i--)
    {
        d2 = distance(patm, at, i);
        if (d2 < r[patm[at].vdw][patm[i].vdw][2])
        {
#ifdef VERBOSE
            std::cout << "Clash 1-5 between atmom = " << at << " and atom = " << i << " distancia= " <<  sqrt(d2) << std::endl;
#endif
            return FILTER_FAIL;
        }
    }

    return FILTER_OK;
}

TreeFilters::FilterResultType TreeFilters::volumeInRange(unsigned int nres, Volume vol_parcial) const
{
    // Valores obtenidos a partir de pruebas de un set de datos en Grillado.
    static const float cota_maxima_volumen = 177.65f;
    static const float pendiente_empirica = -0.0882f;
    static const float volumen_min_aa = 110.0f;

    const float volumen_max_aa = pendiente_empirica * float(nres) + cota_maxima_volumen;
    const float chain_volumen = float(vol_parcial) / float(nres);

#ifdef VERBOSE
    std::cout << "Maximun Volume allowed per a.a  =" << volumen_max_aa << ". Volumen in this chain=" << chain_volumen << endl;
#endif

    return in_range(chain_volumen, volumen_min_aa, volumen_max_aa) ? FILTER_OK : FILTER_FAIL;
}


ClashFilter::ClashFilter(const TreeData& tree_data, const TreeFilters& tree_filters) :
    tree_data(tree_data),
    tree_filters(tree_filters)
{}

bool ClashFilter::operator()(unsigned int index, const Atoms& patm, unsigned int at) const
{
    const bool clash = tree_filters.isclash(patm, at) == TreeFilters::FILTER_FAIL;

    return !clash && (index != 2 || (tree_filters.islong(patm, at, tree_data.dmax2) == TreeFilters::FILTER_OK));
}
