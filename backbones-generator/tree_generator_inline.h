#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "backbones-generator/tree_operator.h" //for RMatrix

template <class TOperator>
inline TreeGenerator<TOperator>::TreeGenerator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* const reader)
    : treeOperator(treeHelper, reader),
      CANT_RES(treeHelper.getNRes()),
      CANT_ANGLES(treeHelper.getNAngles())
{
}

template <class TOperator>
inline void TreeGenerator<TOperator>::generate()
{
    typename TOperator::RMatrix rInicial;
    unsigned int level = 1;
    unsigned int indexSeed = 0;

    treeOperator.initMatrix(rInicial);

    while (treeOperator.putNextSeed(level, indexSeed))
    {
        if (level < CANT_RES)
            expandTree(level, 0);
        else
            processLeaf();

        treeOperator.removeSeed(level);

        indexSeed++;
    }
}

/*
    given a chain C of elements
    for every angles pair P
        for every element E in the collection
            append E in C oriented to P
            recurse
            remove E
*/
/*
 * indice_nivel_anterior es siempre igual al index_angles de la llamada anterior
 */
template <class TOperator>
inline void TreeGenerator<TOperator>::expandTree(unsigned int level, unsigned int previousLevelIndex)
{
    bool lastLevelSuccess = false;//solo interesa si somos el anteultimo nivel
    unsigned int indexAngles = 0;
    typename TOperator::RMatrix rInicial;

    treeOperator.copyMatrix(rInicial);

    while (indexAngles < CANT_ANGLES && !lastLevelSuccess)
    {
        treeOperator.initMatrix(rInicial);

        if (treeOperator.putFirst(level, indexAngles, previousLevelIndex))
        {
            if (level < CANT_RES + 1)
                lastLevelSuccess = appendElements(level, indexAngles);
            else
                lastLevelSuccess = processLeaf();

            treeOperator.removeFirst(level);
        }

        indexAngles++;
    }
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::appendElements(unsigned int level, unsigned int indexAngles)
{
    typename TOperator::KeepRecursion resultRecursion;
    bool result = false;
    unsigned int indexRes = 0;

    while (treeOperator.putNext(level, indexRes, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (level < CANT_RES + 1)
                expandTree(level, indexAngles);
            else
                result = processLeaf();

            treeOperator.remove(level);
        }

        indexRes++;
    }

    return result;
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::processLeaf()
{
    return treeOperator.write();
}
