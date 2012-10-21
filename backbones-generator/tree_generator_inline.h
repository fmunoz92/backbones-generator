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
        if (level < CANT_RES)//if (level < CANT_RES + 1) para que ande igual en singlemode
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
        put single element oriented to P
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

    treeOperator.copyMatrix(rInicial);//is a get method

    while (indexAngles < CANT_ANGLES && !lastLevelSuccess)
    {
        if (treeOperator.putFirst(level, indexAngles, previousLevelIndex))
        {
            if (level < CANT_RES)//if (level < CANT_RES + 1) para que ande igual en singlemode
                lastLevelSuccess = appendElements(level, indexAngles);
            else
                lastLevelSuccess = processLeaf();

            treeOperator.removeFirst(level);
        }
        
        treeOperator.initMatrix(rInicial);//restore the original matrix
        indexAngles++;
    }
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::appendElements(unsigned int level_arg, unsigned int indexAngles)
{
	unsigned int level = level_arg;
    typename TOperator::KeepRecursion resultRecursion;
    bool result = false;
    unsigned int indexRes = 0;

    while (treeOperator.putNext(level, indexRes, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (level < CANT_RES)//if (level < CANT_RES + 1) para que ande igual en singlemode
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
