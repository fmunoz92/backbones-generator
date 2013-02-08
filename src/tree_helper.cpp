#include "backbones-generator/tree_helper.h"

#include <fstream>

TreeHelper::TreeHelper(TreeData& treeData, IncrementalBackbone& incrementalBackbone)
    : treeData(treeData),
      incrementalBackbone(incrementalBackbone)
{}

IncrementalBackbone& TreeHelper::getAtm()
{
    return incrementalBackbone;
}

TreeData& TreeHelper::getData()
{
    return treeData;
}
