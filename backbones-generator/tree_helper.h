#ifndef TREE_HELPER_H
#define TREE_HELPER_H

#include <list>
#include "backbones-generator/tree_data.h"
#include "backbones-generator/tree_filters.h"

class TreeHelper
{
public:

    TreeHelper(TreeData& tree_data, TreeFilters& tree_filters);

    void                          putSeed(float* R, Residuo& residuo);
    TreeFilters::FilterResultType putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int si_index, unsigned int fi_index);
    TreeFilters::FilterResultType putChain(float* pR, unsigned int resN, std::list<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chain_index);

    void deleteRes(const Residuo& residuo);
    void deleteRes(const std::list<Residuo>& residuos);

    bool filterLastLevelOk();

    void clearatm();//This function puts 0 in all the atoms coordinates
    void reportSuccess();
    void deleteLastFragmentId();

    prot_filer::AnglesData&  getAnglesData()  const;
    Atoms&                   getAtm()         const;
    prot_filer::FragmentIds& getFragmentIds() const;
    const std::string&       getOutputFile()  const;
    unsigned int             getNRes()        const;
    unsigned int             getNAngles()     const;
    bool                     success()        const;

private:
    TreeData& tree_data;
    const TreeFilters& tree_filters;
};

#endif
