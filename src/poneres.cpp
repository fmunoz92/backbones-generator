#include "utils.h"
#include "poneres.h"

class ClashFilter
{
public:
    ClashFilter(TreeData& tree_data)
        : tree_data(tree_data)
    {}
    bool operator()(unsigned int index, const ATOM* patm, int at) const
    {
        const bool clash = isclash(patm, at) == FILTER_FAIL;
        return !clash && (index != 2 || (islong(patm, at, tree_data.dmax2) != FILTER_FAIL));
    }
private:
    const TreeData& tree_data;
};

FilterResultType poneres(float* pR, const unsigned int resN, TreeData& tree_data, Residuo& residuo, unsigned int si_index, unsigned int fi_index)
{
    float cossi = tree_data.cossi[si_index];
    float sinsi = tree_data.sinsi[si_index];
    float cosfi = tree_data.cosfi[fi_index];
    float sinfi = tree_data.sinfi[fi_index];
    ATOM* patm = tree_data.atm;

    const unsigned int i = resN - 2;
    tree_data.angles_data->angles[i].si = si_index;
    tree_data.angles_data->angles[i].fi = fi_index;

#ifdef COMBINATIONS_DEBUG
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(tree_data);
#endif
    bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, patm, resN, filter);

    if (!success)
    {
        return FILTER_FAIL;
    }

    const ATOM atm = patm[3 * (resN - 1) + 1];
    residuo.at2 = tree_data.grilla->agregar_esfera(atm.x, atm.y, atm.z);
    return FILTER_OK;
}

FilterResultType addChain(float* pR, unsigned int resN, TreeData& tree_data, vector<Residuo> &residuos, const IncompleteAnglesData& chain, unsigned int chain_index)
{
    FilterResultType result = FILTER_OK;
    unsigned int i = 0;
    while (result == FILTER_OK  && i < (chain.nres - 1) && ((resN + i) <= tree_data.nres))
    {
        Residuo residuo;
        const unsigned int si = chain.angles[i].si;
        const unsigned int fi = chain.angles[i].fi;
        result = poneres(pR, resN + i, tree_data, residuo, si, fi);
        if (result == FILTER_OK)
        {
            mili::insert_into(residuos, residuo);
        }
        ++i;
    }
    if (result == FILTER_OK)
    {
        tree_data.chain_indexs.push_back(chain_index);
    }
    return result;
}
