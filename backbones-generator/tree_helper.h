#ifndef TREE_HELPER_H
#define TREE_HELPER_H

#include <list>

#include "backbones-generator/tree_data.h"
#include "backbones-generator/tree_filters.h"

class TreeHelper
{
public:

    TreeHelper(TreeData& treeData, TreeFilters& treeFilters);

    void reportSuccess();

    IncrementalBackbone& getAtm();

    const prot_filer::AnglesData&  getAnglesData()  const;
    const prot_filer::FragmentIds& getFragmentIds() const;

    unsigned int getNRes()        const;
    unsigned int getNAngles()     const;
    bool         success()        const;

private:
    TreeData& treeData;
    const TreeFilters& treeFilters;
};

#endif
