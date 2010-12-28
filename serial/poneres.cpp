#include "poneres.h"

#include "copymat.h"
#include "isclash.h"
#include "int2car.h"
#include "islong.h"

// En el modo DEBUG se deshabilitan los chequeos por lo que
// siempre devuelve FILTER_OK.

FilterResultType poneres(float* pR, const unsigned int resN, TreeData* tree_data, Residuo& residuo, unsigned int si_index, unsigned int fi_index)
{
    float cossi = tree_data->cossi[si_index];
    float sinsi = tree_data->sinsi[si_index];
    float cosfi = tree_data->cosfi[fi_index];
    float sinfi = tree_data->sinfi[fi_index];
    ATOM* patm = tree_data->atm;
    float T[16];

    const unsigned int i = resN - 2;
    tree_data->angles_data->angles[i].si = si_index;
    tree_data->angles_data->angles[i].fi = fi_index;

    /*Guardo la anterior*/
    copymat(T, pR);
    unsigned int at = 3 * (resN - 1);

#ifdef COMBINATIONS_DEBUG
    int2car(T, b_C_N, cos_a_CA_C_N, sin_a_CA_C_N, cossi, sinsi, patm, at, N);

    at++;
    int2car(T, b_N_CA, cos_a_C_N_CA, sin_a_C_N_CA, cos_OMEGA, sin_OMEGA, patm, at, CA);

    at++;
    int2car(T, b_CA_C, cos_a_N_CA_C, sin_a_N_CA_C, cosfi, sinfi, patm, at, C);
#else
    int2car(T, b_C_N, cos_a_CA_C_N, sin_a_CA_C_N, cossi, sinsi, patm, at, N);
    if (isclash(patm, at) == FILTER_FAIL)
    {
        return FILTER_FAIL;
    }

    at++;
    int2car(T, b_N_CA, cos_a_C_N_CA, sin_a_C_N_CA, cos_OMEGA, sin_OMEGA, patm, at, CA);
    if (isclash(patm, at) == FILTER_FAIL)
    {
        return FILTER_FAIL;
    }
    if (islong(patm, at, tree_data->dmax2) == FILTER_FAIL)
    {
        return FILTER_FAIL;
    }

    at++;
    int2car(T, b_CA_C, cos_a_N_CA_C, sin_a_N_CA_C, cosfi, sinfi, patm, at, C);
    if (isclash(patm, at) == FILTER_FAIL)
    {
        return FILTER_FAIL;
    }
#endif
    /*
    residuo = Residuo(tree_data->grilla->agregar_esfera(patm[at-2].x,patm[at-2].y,patm[at-2].z), tree_data->grilla->agregar_esfera(patm[at-1].x,patm[at-1].y,patm[at-1].z), tree_data->grilla->agregar_esfera(patm[at].x,patm[at].y,patm[at].z));
    */
    const ATOM atm = patm[at - 1];
    residuo.at2 = tree_data->grilla->agregar_esfera(atm.x, atm.y, atm.z);

    copymat(pR, T);
    return FILTER_OK;
}

FilterResultType addNRes(float* pR, unsigned int resN, TreeData* tree_data, vector<Residuo> &residuos, const IncompleteAnglesData& residue_chain)
{
    FilterResultType result = FILTER_OK;
    unsigned int i = 0;
    while (result == FILTER_OK  && i < residue_chain.nres && resN <= tree_data->nres)
    {
        Residuo residuo;
        result = poneres(pR, resN, tree_data, residuo, residue_chain.angles[i].si, residue_chain.angles[i].fi);
        if (result == FILTER_OK)
        {
            mili::insert_into(residuos, residuo);
        }
        ++i;
        ++resN;
    }
    return result;
}
