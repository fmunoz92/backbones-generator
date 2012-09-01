#ifndef TREE_OPERATOR_H
#define TREE_OPERATOR_H

#include "backbones-generator/tree_helper.h"

template <class WriterHelper>
class TreeOperator
{
public:

    typedef float RMatrix[16];

    enum KeepRecursion
    {
        DoRecursion,
        StopRecursion
    };

    inline TreeOperator(TreeHelper& treeHelper);
    inline TreeOperator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* reader);

    inline void initMatrix(const RMatrix rMatrix);
    inline void copyMatrix(RMatrix rInicial) const;

    inline bool putFirst(unsigned int& level, unsigned int fiIndex, unsigned int siIndex);
    inline void removeFirst(unsigned int& levell);

    inline bool write();

protected:
    RMatrix R;
    std::list<Residuo> residuos;
    TreeHelper& treeHelper;

private:
    inline bool lastLevelOk() const;
    WriterHelper writerHelper;
};

/**********************************************************************/

template <class WriterHelper>
class SimpleTreeOperator : public TreeOperator<WriterHelper>
{
public:
    inline SimpleTreeOperator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& level, unsigned int indexSeed);
    inline bool putNext(unsigned int& level, unsigned int indexRes, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion);

    inline void remove(unsigned int& level);
    inline void removeSeed(unsigned int& level);
};

/**********************************************************************/

template <class WriterHelper>
class ChainsTreeOperator : public TreeOperator<WriterHelper>
{
public:
    inline ChainsTreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& level, unsigned int indexSeed);
    inline bool putNext(unsigned int& level, unsigned int indexRes, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion);

    inline void remove(unsigned int& level);
    inline void removeSeed(unsigned int& level);

private:
    inline void putChain(prot_filer::AnglesData& chain, unsigned int& level, unsigned int indexRes, typename TreeOperator<WriterHelper>::KeepRecursion& recursion);
    inline void putSeed(prot_filer::AnglesData& chain, unsigned int& level, unsigned int indexSeed);

    std::list<std::list<Residuo> > vectoresParaBorrar;
    FullCachedAnglesSeqReader* const reader;
};

#define TREE_OPERATOR_INLINE_H
#include "tree_operator_inline.h"
#undef TREE_OPERATOR_INLINE_H

#endif
