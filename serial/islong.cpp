#include "islong.h"
#include <cmath>
//This function filter long chains.
//Calculates the CA-CA distances of the recently added CA with the rest.
//Then compares this distances with dmax2, wich cames from the analisys of protein
//structures experimetally determined

FilterResultType islong(const ATOM* patm, int at, float dmax2)
{
    // Until we reach residue #5, the check for long chain is meaningless
    // at=13 is the CA of residue #5
    if (at < 13)
    {
        return FILTER_OK;
    }

    // at-12 is the CA atom that is four residues down the chain
    // at=1 is the first CA in the chain
    for (int i = at - 12; i >= 1; i -= 3)
    {
        const float dx2 = mili::square(patm[at].x - patm[i].x);
        const float dy2 = mili::square(patm[at].y - patm[i].y);
        const float dz2 = mili::square(patm[at].z - patm[i].z);
        const float d2  = dx2 + dy2 + dz2;
        if (d2 > dmax2)
        {
#ifdef VERBOSE
            printf("Chain length = %f, while Dmax is= %f \n", sqrt(d2), sqrt(dmax2));
#endif
            return FILTER_FAIL;
        }
    }
    return FILTER_OK;
}

