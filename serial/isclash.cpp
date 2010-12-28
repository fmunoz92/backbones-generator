#include "isclash.h"
#include <cmath>

float r[3][3][3];

//This function detect colitions between atoms
//The matrix r have the reference radius (actually the sum of the squares of the raduis)
//The last index of r is 0 for 1-4 clashes, 1 for 1-5 and 2 for the rest

inline float distance(const ATOM* patm, const int at, const int i)
{
    const float dx2 = mili::square(patm[at].x - patm[i].x);
    const float dy2 = mili::square(patm[at].y - patm[i].y);
    const float dz2 = mili::square(patm[at].z - patm[i].z);
    return dx2 + dy2 + dz2;
}

FilterResultType isclash(const ATOM* patm, int at)
{
    //This is to check for the so-called 1-4 clashes,
    //i.e. a clash between atom at position i with atom at position i+3
    int i  = at - 3;
    float d2 = distance(patm, at, i);

    if (d2 < r[patm[at].vdw][patm[i].vdw][0])
    {
#ifdef VERBOSE
        printf("Clash 1-4 between atmom=%i and atom=%i distancia=%2.3f\n", at, i, sqrt(d2));
#endif
        return FILTER_FAIL;
    }

    //If this is atom #3 and passes the check for 1-4 clashes, then there is nothing else to check.
    if (at == 3)
    {
        return FILTER_OK;
    }

    //This is to check for the so-called 1-5 clashes;
    //i.e. a clash between atom at position i with atom at position i+4
    i = at - 4;
    d2 = distance(patm, at, i);
    if (d2 < r[patm[at].vdw][patm[i].vdw][1])
    {
#ifdef VERBOSE
        printf("Clash 1-5 between atmom=%i and atom=%i distancia=%2.3f\n", at, i, sqrt(d2));
#endif
        return FILTER_FAIL;
    }

    //If this is atom #4 and passes check for 1-4 and 1-5 clashes, then there is nothing else to check.
    if (at == 4)
    {
        return FILTER_OK;
    }

    //The rest of the clashes until the end of the chain.
    for (i = at - 5; i >= 0; i--)
    {
        d2 = distance(patm, at, i);
        if (d2 < r[patm[at].vdw][patm[i].vdw][2])
        {
#ifdef VERBOSE
            printf("Clash > 1-5 between atmom=%i and atom=%i distancia=%2.3f\n", at, i, sqrt(d2));
#endif
            return FILTER_FAIL;
        }
    }

    return FILTER_OK;
}

