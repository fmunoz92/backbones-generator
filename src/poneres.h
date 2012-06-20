/*
 * 
 * rename this file
 * 
 */ 
#ifndef PONERES_H
#define PONERES_H

#include <vector>
#include "tree_data.h"
#include "utils.h"

struct TreeHelper
{
    static inline void semilla(TreeData& tree_data, float* R, Residuo& residuo);
    static inline FilterResultType filtros_ultimo_nivel(TreeData& tree_data);
    static inline void sacar_residuo(TreeData& tree_data, const Residuo& residuo);
    static inline void sacar_residuos(TreeData& tree_data, const vector<Residuo>& residuos); 
    static inline FilterResultType poner_residuo(float* pR, const unsigned int resN, TreeData& tree_data, Residuo& residuo, unsigned int si_index, unsigned int fi_index);
    static inline FilterResultType addChain(float* pR, unsigned int resN, TreeData& tree_data, vector<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chain_index);
};

class ClashFilter
{
public:
    ClashFilter(const TreeData& tree_data);
    bool operator()(unsigned int index, const Atoms& patm, int at) const;
private:
    const TreeData& tree_data;
};


inline void TreeHelper::semilla(TreeData& tree_data, float* R, Residuo& residuo)
{
    Atoms& atm = tree_data.atm;
    backbones_utils::semilla(atm, R);
    prot_filer::ATOM* seed = tree_data.angles_data->seed;

    seed[0] = atm[0];
    seed[1] = atm[1];
    seed[2] = atm[2];

    residuo.at2 = tree_data.grilla->agregar_esfera(tree_data.atm[1].x, tree_data.atm[1].y, tree_data.atm[1].z);
}

inline FilterResultType TreeHelper::filtros_ultimo_nivel(TreeData& tree_data)
{
    bool ok = calcRdG(tree_data.atm, tree_data.nres, tree_data.rgmax) == FILTER_OK
              && volumen_en_rango(tree_data.nres, tree_data.grilla->obtener_vol_parcial()) == FILTER_OK;
    return ok ? FILTER_OK : FILTER_FAIL;
}

inline void TreeHelper::sacar_residuo(TreeData& tree_data, const Residuo& residuo)
{
    tree_data.grilla->sacar_esfera(residuo.at2);
}

inline void TreeHelper::sacar_residuos(TreeData& tree_data, const vector<Residuo>& residuos)
{
    for (unsigned int i = 0; i < residuos.size(); ++i)
    {
        sacar_residuo(tree_data, residuos[i]);
    }
}

inline FilterResultType TreeHelper::poner_residuo(float* pR, const unsigned int resN, TreeData& tree_data, Residuo& residuo, unsigned int si_index, unsigned int fi_index)
{
    float cossi = tree_data.cossi[si_index];
    float sinsi = tree_data.sinsi[si_index];
    float cosfi = tree_data.cosfi[fi_index];
    float sinfi = tree_data.sinfi[fi_index];
    Atoms& patm = tree_data.atm;

    const unsigned int i = resN - 2;
    tree_data.angles_data->angles[i].si = si_index;
    tree_data.angles_data->angles[i].fi = fi_index;

#ifdef COMBINATIONS_DEBUG
	#include "prot-filer/backbones_utils.h"
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(tree_data);
#endif
    bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, patm, resN, filter);

    if (!success)
    {
        return FILTER_FAIL;
    }

    const prot_filer::ATOM& atm = patm[3 * (resN - 1) + 1];
    residuo.at2 = tree_data.grilla->agregar_esfera(atm.x, atm.y, atm.z);
    return FILTER_OK;
}

inline FilterResultType TreeHelper::addChain(float* pR, unsigned int resN, TreeData& tree_data, vector<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chain_index)
{
    FilterResultType result = FILTER_OK;
    unsigned int i = 0;
    while (result == FILTER_OK  && i < (chain.nres - 1) && ((resN + i - 1) < (tree_data.nres - 1)))
    {
        Residuo residuo;
        const unsigned int si = chain.angles[i].si;
        const unsigned int fi = chain.angles[i].fi;
        result = poner_residuo(pR, resN + i, tree_data, residuo, si, fi);
        if (result == FILTER_OK)
        {
            mili::insert_into(residuos, residuo);
        }
        ++i;
    }
    if (result == FILTER_OK)
    {
        tree_data.fragment_ids.push_back(chain_index);
    }
    return result;
}


#endif
