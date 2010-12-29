#include "poneres.h"

#include "isclash.h"
#include "islong.h"

class ClashFilter
{
public:
    ClashFilter(TreeData* tree_data)
        : tree_data(tree_data)
    {}
    bool operator()(unsigned int index, const ATOM* patm, int at) const
    {
        const bool clash = isclash(patm, at) == FILTER_FAIL;
        return !clash && (index != 2 || (islong(patm, at, tree_data->dmax2) != FILTER_FAIL));
    }
private:
    const TreeData* tree_data;
};

//TODO: temp, deberia utilizarse backbones_utils::poneres
template <class F>
bool poneres2(float* pR, float cossi, float sinsi, float cosfi, float sinfi, ATOM* patm, int resN, const F& filter)
{
    using namespace backbones_utils;
    float T[16];

    /*Guardo la anterior*/
    copymat(T, pR);

    int at = 3 * (resN - 1) - 1;

    at++;
    int2car(T, b_C_N, cos_a_CA_C_N, sin_a_CA_C_N, cossi, sinsi, patm, at, N);
    if (!filter(1, patm, at))
    {
        return false;
    }

    at++;
    int2car(T, b_N_CA, cos_a_C_N_CA, sin_a_C_N_CA, cos_OMEGA , sin_OMEGA, patm, at, CA);
    if (!filter(2, patm, at))
    {
        return false;
    }

    at++;
    int2car(T, b_CA_C, cos_a_N_CA_C, sin_a_N_CA_C, cosfi, sinfi, patm, at, C);
    if (!filter(3, patm, at))
    {
        return false;
    }

    copymat(pR, T);
    return true;
}

FilterResultType poneres(float* pR, const unsigned int resN, TreeData* tree_data, Residuo& residuo, unsigned int si_index, unsigned int fi_index)
{
    float cossi = tree_data->cossi[si_index];
    float sinsi = tree_data->sinsi[si_index];
    float cosfi = tree_data->cosfi[fi_index];
    float sinfi = tree_data->sinfi[fi_index];
    ATOM* patm = tree_data->atm;

    const unsigned int i = resN - 2;
    tree_data->angles_data->angles[i].si = si_index;
    tree_data->angles_data->angles[i].fi = fi_index;

#ifdef COMBINATIONS_DEBUG
// En el modo DEBUG se deshabilitan los chequeos por lo que
// siempre devuelve FILTER_OK.
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(tree_data);
#endif
    //bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, patm, resN, filter);
    bool success = poneres2(pR, cossi, sinsi, cosfi, sinfi, patm, resN, filter);

    if (!success)
    {
        return FILTER_FAIL;
    }

    const ATOM atm = patm[3 * (resN - 1) + 1];
    residuo.at2 = tree_data->grilla->agregar_esfera(atm.x, atm.y, atm.z);
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
