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
    TreeHelper(TreeData& tree_data, TreeFilters& tree_filters);

    void             putSeed(float* R, Residuo& residuo);
    FilterResultType putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int si_index, unsigned int fi_index);
    FilterResultType putChain(float* pR, unsigned int resN, std::list<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chain_index);

    void deleteRes(const Residuo& residuo);
    void deleteRes(const std::list<Residuo>& residuos);

    FilterResultType filtros_ultimo_nivel();

    void clearatm();//This function puts 0 in all the atoms coordinates
    void reportSuccess();
    void deleteLastFragmentId();

    prot_filer::AnglesData&  getAnglesData()  const;
    Atoms&                   getAtm()         const;
    prot_filer::FragmentIds& getFragmentIds() const;
    const std::string&       getOutputFile()  const;
    unsigned int             getNRes()        const;
    unsigned int             getNAngles()     const;

private:
    TreeData& tree_data;
    const TreeFilters& tree_filters;
};

class ClashFilter
{
public:
    ClashFilter(const TreeData& tree_data, const TreeFilters& tree_filters);
    bool operator()(unsigned int index, const Atoms& patm, int at) const;
private:
    const TreeData& tree_data;
    const TreeFilters& tree_filters;
};


#endif
