#ifndef UTILS_H
#define UTILS_H

#include "tree_data.h"// for global enum FilterResultType

class TreeFilters
{
public:
    void setr(float rn, float rca, float rc, float scal_1_4, float scal_1_5);

    //This function calculates the gyration radius of the CA atoms
    //It was meant to filter long chains. Now it is replaced by the volume filter
    FilterResultType calcRdG(const Atoms& patm, int nres, float rgmax) const;

    //This function filter long chains.
    //Calculates the CA-CA distances of the recently added CA with the rest.
    //Then compares this distances with dmax2, wich cames from the analisys of protein
    //structures experimetally determined
    FilterResultType islong(const Atoms& patm, int at, float dmax2) const;

    FilterResultType isclash(const Atoms& patm, int at) const;

    FilterResultType volumen_en_rango(int nres, Volume vol_parcial) const;
private:

    //This function detect colitions between atoms
    //The matrix r have the reference radius (actually the sum of the squares of the raduis)
    //The last index of r is 0 for 1-4 clashes, 1 for 1-5 and 2 for the rest
    static inline float distance(const Atoms& patm, const int at, const int i)
    {
        const float dx2 = square(patm[at].x - patm[i].x);
        const float dy2 = square(patm[at].y - patm[i].y);
        const float dz2 = square(patm[at].z - patm[i].z);
        return dx2 + dy2 + dz2;
    }

    float r[3][3][3];
};


#endif
