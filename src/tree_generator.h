#ifndef TREE_GENERATOR_H
#define TREE_GENERATOR_H

#include <string>
#include <list>
#include <mili/mili.h>
#include "tree_helper.h"

template <class TOperator>
class TreeGenerator
{
public:
    inline TreeGenerator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader);

    inline void generate();
private:
    inline void expandTree(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    inline bool processLeaf();
    inline bool appendElements(unsigned int nivel, unsigned int indice_nivel_anterior, unsigned int index, const float R_local[16]);

    TOperator treeOperator;
    const unsigned int CANT_RES;
    const unsigned int CANT_ANGLES;
};

struct TreeOperator
{
    enum KeepRecursion {DoRecursion, StopRecursion};
};

template <class WriterHelper>
class SimpleTreeOperator : public TreeOperator
{
public:
    inline SimpleTreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& nivel, unsigned int index_seed);
    inline void initMatrix(float R[16]);
    inline bool putNext(unsigned int& nivel, unsigned int index_res, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion);
    inline void remove(unsigned int& nivel);
    inline void write();
    inline bool lastLevelOk();
private:
    TreeHelper& tree_helper;
    mili::FirstTimeFlag firstTime;
    float* R;
    std::list<Residuo> paraBorrar;
    WriterHelper writer_helper;
};

template <class WriterHelper>
class ChainsTreeOperator : public TreeOperator
{
public:
    inline ChainsTreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& nivel, unsigned int index_seed);
    inline void initMatrix(float R[16]);
    inline bool putNext(unsigned int& nivel, unsigned int index_res, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion);
    inline void remove(unsigned int& nivel);
    inline void write();
    inline bool lastLevelOk();
private:
    TreeHelper& tree_helper;
    mili::FirstTimeFlag firstTime;
    float* R;
    std::list<Residuo> residuosParaBorrar;
    std::list<std::list<Residuo> > vectoresParaBorrar;
    FullCachedAnglesSeqReader* const reader;
    WriterHelper writer_helper;
};

#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
