#ifndef TREE_GENERATOR_H
#define TREE_GENERATOR_H

#include <string>
#include <list>
#include <mili/mili.h>
#include "tree_helper.h"

typedef float RMatrix[16]; //TODO: move this global definition to TreeOperator class

template <class TOperator>
class TreeGenerator
{
public:
    inline TreeGenerator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader);

    inline void generate();

private:
    inline void expandTree(unsigned int nivel, const RMatrix R_inicial, unsigned int indice_nivel_anterior);
    inline bool processLeaf();
    inline bool appendElements(unsigned int nivel, unsigned int index);

    TOperator treeOperator;
    const unsigned int CANT_RES;
    const unsigned int CANT_ANGLES;
};

template <class WriterHelper>
class TreeOperator
{
public:
    enum KeepRecursion {DoRecursion, StopRecursion};

    inline TreeOperator(TreeHelper& tree_helper);
    inline TreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline void initMatrix(const RMatrix rMatrix);
    inline void copyMatrix(RMatrix R_inicial);

    inline bool putFirst(unsigned int& nivel, unsigned int fi_index, unsigned int si_index);
    inline void removeFirst(unsigned int& nivel);

    inline bool write();

protected:
    RMatrix R;
    std::list<Residuo> residuos;
    TreeHelper& tree_helper;

private:
    inline bool lastLevelOk();
    WriterHelper writer_helper;
};

template <class WriterHelper>
class SimpleTreeOperator : public TreeOperator<WriterHelper>
{
public:
    inline SimpleTreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& nivel, unsigned int index_seed);
    inline bool putNext(unsigned int& nivel, unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion);

    inline void remove(unsigned int& nivel);
    inline void removeSeed(unsigned int& nivel);
};

template <class WriterHelper>
class ChainsTreeOperator : public TreeOperator<WriterHelper>
{
public:
    inline ChainsTreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& nivel, unsigned int index_seed);
    inline bool putNext(unsigned int& nivel, unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion);

    inline void remove(unsigned int& nivel);
    inline void removeSeed(unsigned int& nivel);

private:
    inline void putChain(prot_filer::AnglesData& chain, unsigned int& nivel, unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& recursion);
    inline void putSeed(prot_filer::AnglesData& chain, unsigned int& nivel, unsigned int index_seed);

    std::list<std::list<Residuo> > vectoresParaBorrar;
    FullCachedAnglesSeqReader* const reader;
};

#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
