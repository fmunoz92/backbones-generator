#include <cmath>
#include <iostream>
#include "backbones-generator/tree_filters.h"

using mili::square;

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
    const bool isNesessaryCheck = (at >= 13);

    // at-12 is the CA atom that is four residues down the chain
    // at=1 is the first CA in the chain
    int i = at - 12;
    bool result = isNesessaryCheck;
    while (i >= 1 && result)
    {
        const float dx2 = square(patm[at].x - patm[i].x);
        const float dy2 = square(patm[at].y - patm[i].y);
        const float dz2 = square(patm[at].z - patm[i].z);
        const float d2  = dx2 + dy2 + dz2;

        result = (d2 <= dmax2);
        i -= 3;
#ifdef VERBOSE
        if (result)
            std::cout << "Chain length = " << sqrt(d2) << ", while Dmax is = " << sqrt(dmax2) << std::endl;
#endif
    }

    return (result || !isNesessaryCheck) ? FILTER_OK : FILTER_FAIL;
}

bool TreeFilters::isParticularClash(const Atoms& patm, unsigned int x, unsigned int y, unsigned int z) const
{
    float d2 = distance(patm, x, y);
    const bool result = (d2 < r[patm[x].vdw][patm[y].vdw][z]);
#ifdef VERBOSE
    if (result)
        std::cout << "Clash between atmom = " << x << " and atom = " << z << " distancia= " <<  sqrt(d2) << std::endl;
#endif

    return result;
}

TreeFilters::FilterResultType TreeFilters::isclash(const Atoms& patm, unsigned int at) const
{
    //This is to check for the so-called 1-4 clashes,
    //i.e. a clash between atom at position i with atom at position i+3

    bool result = true;
    unsigned int z = 0;
    int i  = at - 3;

    while (i >= 0 && result)
    {
        //If this is atom #3 and passes the check for 1-4 clashes or this is atom #4 and passes
        //check for 1-4 and 1-5 clashes, then there is nothing else to check.
        const bool nothingElseToCheck = (i == 0) || (i == 1);

        result = isParticularClash(patm, at, i, z) && !nothingElseToCheck;

        --i;
        if (z < 2)
            ++z;
    }

    return (result) ? FILTER_OK : FILTER_FAIL;
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
