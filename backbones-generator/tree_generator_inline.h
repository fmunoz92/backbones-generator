#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif
// modo cadenas:
// fragmentos = 4, angulos = 2, chain_size = 3, nres = 9
// cadenas = ceil(nres / chain_size)
// soluciones = fragmentos^cadenas * angulos^(cadenas - 1)

        
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
    const unsigned int level = 2;
    const unsigned int indexAngles = 0;
    const unsigned int previousLevelIndex = 0;

    treeOperator.initMatrix(rInicial);
    treeOperator.putSeed();

    if (level < CANT_RES + 1)
        appendElements(level, indexAngles, previousLevelIndex);
    else
        processLeaf();

    treeOperator.removeSeed();
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
