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

    inline TreeOperator(TreeHelper& tree_helper);
    inline TreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);

    inline void initMatrix(const RMatrix rMatrix);
    inline void copyMatrix(RMatrix R_inicial) const;

    inline bool putFirst(unsigned int& nivel, unsigned int fi_index, unsigned int si_index);
    inline void removeFirst(unsigned int& nivel);

    inline bool write();

protected:
    RMatrix R;
    std::list<Residuo> residuos;
    TreeHelper& tree_helper;

private:
    inline bool lastLevelOk() const;
    WriterHelper writer_helper;
};

/**********************************************************************/

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

/**********************************************************************/

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

#define TREE_OPERATOR_INLINE_H
#include "tree_operator_inline.h"
#undef TREE_OPERATOR_INLINE_H

#endif
