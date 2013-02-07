#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

template <class TOperator>
inline TreeGenerator<TOperator>::TreeGenerator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* const reader)
    : treeOperator(treeHelper, reader),
      CANT_RES(treeHelper.getData().nRes),
      CANT_ANGLES(treeHelper.getData().nAngles),
      count(0)
{
}

template <class TOperator>
inline unsigned int TreeGenerator<TOperator>::generate()
{
    typename TOperator::RMatrix rInicial;
    const unsigned int level = 1;//seed is level 0
    const unsigned int indexAngles = 0;
    const unsigned int previousLevelIndex = 0;

    treeOperator.initMatrix(rInicial);
    treeOperator.putSeed();

    if (level < CANT_RES)
        appendElements(level, indexAngles, previousLevelIndex);
    else
        processLeaf();

    treeOperator.removeSeed();

    return count;
}

/*
    given a chain C of elements
    for every angles pair P
        for every element E in the collection
            append E in C oriented to P
            recurse
            remove E
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
        lastLevelSuccess = appendElements(level, indexAngles, previousLevelIndex);

        treeOperator.initMatrix(rInicial);//restore the original matrix
        ++indexAngles;
    }
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::appendElements(unsigned int level, unsigned int indexAngles, unsigned int previousLevelIndex)
{
    typename TOperator::KeepRecursion resultRecursion;
    bool result = false;
    unsigned int indexRes = 0;

    while (treeOperator.putNext(level, indexRes, indexAngles, previousLevelIndex, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (level < CANT_RES)
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
    const bool success = treeOperator.write();

#ifdef COMBINATIONS_DEBUG
    ++count;
#else
    if (success)
        ++count;
#endif


    return success;
}
