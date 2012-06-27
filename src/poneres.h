/*
 *
 * rename this file
 *
 */
#ifndef PONERES_H
#define PONERES_H

#include <list>
#include "tree_data.h"
#include "utils.h"

class TreeHelper
{
public:
    TreeHelper(TreeData& tree_data);

    void putSeed(float* R, Residuo& residuo);
    FilterResultType putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int si_index, unsigned int fi_index);
    FilterResultType putChain(float* pR, unsigned int resN, list<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chain_index);

    void deleteRes(const Residuo& residuo);
    void deleteRes(const list<Residuo>& residuos);

    FilterResultType filtros_ultimo_nivel();
private:
    TreeData& tree_data;
};

class ClashFilter
{
public:
    ClashFilter(const TreeData& tree_data);
    bool operator()(unsigned int index, const Atoms& patm, int at) const;
private:
    const TreeData& tree_data;
};


#endif
