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

    inline void putSeed();
    inline void removeSeed();
    inline bool write();

protected:
    Residuo semilla;
    RMatrix R;
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

    inline bool putNext(unsigned int& level, unsigned int /*indexRes*/, unsigned int indexAngles, unsigned int previousLevelIndex, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion);

    inline void remove(unsigned int& level);
private:
    std::list<Residuo> residuos;
};

/**********************************************************************/

template <class WriterHelper>
class ChainsTreeOperator : public TreeOperator<WriterHelper>
{
public:
    inline ChainsTreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);
    inline bool putNext(unsigned int& level, unsigned int /*indexRes*/, unsigned int indexAngles, unsigned int previousLevelIndex, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion);
    inline void remove(unsigned int& level);

private:
    inline void putChain(prot_filer::AnglesData& chain, unsigned int& level, const unsigned int indexRes, unsigned int firstSi, unsigned int firstFi, typename TreeOperator<WriterHelper>::KeepRecursion& recursion);

    typedef std::list<Residuo> ChainsRes;
    typedef std::list<ChainsRes > StackChainRes;

    FullCachedAnglesSeqReader* const reader;
    StackChainRes stackChainRes;
};

#define TREE_OPERATOR_INLINE_H
#include "tree_operator_inline.h"
#undef TREE_OPERATOR_INLINE_H

#endif
