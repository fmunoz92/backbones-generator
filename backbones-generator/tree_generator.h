#ifndef TREE_GENERATOR_H
#define TREE_GENERATOR_H

#include <string>
#include <list>
#include <mili/mili.h>

#include "backbones-generator/tree_helper.h"


template <class TOperator>
class TreeGenerator
{
public:
    inline TreeGenerator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader);

    inline void generate();

private:
    inline void expandTree(unsigned int nivel, unsigned int indice_nivel_anterior);

    inline bool appendElements(unsigned int nivel, unsigned int index);

    inline bool processLeaf();

    TOperator treeOperator;
    const unsigned int CANT_RES;
    const unsigned int CANT_ANGLES;
};


#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
