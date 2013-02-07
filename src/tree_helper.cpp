#include "backbones-generator/tree_helper.h"

#include <fstream>

TreeHelper::TreeHelper(TreeData& treeData, TreeFilters& treeFilters)
    : treeData(treeData),
      treeFilters(treeFilters)
{}

bool TreeHelper::success() const
{
    return treeData.hubo_algun_exito;
}

void TreeHelper::pushChainIndex(unsigned int index)
{
    treeData.fragmentIds.push_back(index);
}

void TreeHelper::reportSuccess()
{
    ++treeData.cont;
#ifndef COMBINATIONS_DEBUG
    treeData.hubo_algun_exito = true;
#endif
}

void TreeHelper::deleteLastFragmentId()
{
    treeData.fragmentIds.pop_back();
}

const prot_filer::AnglesData&  TreeHelper::getAnglesData() const
{
    return treeData.incrementalBackbone.getAnglesData();
}


IncrementalBackbone& TreeHelper::getAtm()
{
    return treeData.incrementalBackbone;
}

const prot_filer::FragmentIds& TreeHelper::getFragmentIds() const
{
    return treeData.fragmentIds;
}

unsigned int TreeHelper::getNRes() const
{
    return treeData.nres;
}

unsigned int TreeHelper::getNAngles() const
{
    return treeData.cossi.size();
}
