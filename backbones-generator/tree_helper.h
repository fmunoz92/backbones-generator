#ifndef TREE_HELPER_H
#define TREE_HELPER_H

#include <list>

#include "backbones-generator/tree_data.h"
#include "backbones-generator/tree_filters.h"

class TreeHelper
{
public:
    TreeHelper(TreeData& treeData, IncrementalBackbone& incrementalBackbone);

    IncrementalBackbone& getAtm();
    TreeData& getData();
private:
    TreeData& treeData; // static data
    IncrementalBackbone& incrementalBackbone; // partial structure
};

#endif
