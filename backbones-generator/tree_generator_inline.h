#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "backbones-generator/tree_operator.h" //for RMatrix

template <class TOperator>
inline TreeGenerator<TOperator>::TreeGenerator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader)
    : treeOperator(tree_helper, reader),
      CANT_RES(tree_helper.getNRes()),
      CANT_ANGLES(tree_helper.getNAngles())
{
}

template <class TOperator>
inline void TreeGenerator<TOperator>::generate()
{
    typename TOperator::RMatrix R_inicial;
    unsigned int nivel = 1;
    unsigned int index_seed = 0;

    treeOperator.initMatrix(R_inicial);

    while (treeOperator.putNextSeed(nivel, index_seed))
    {
        if (nivel < CANT_RES)
            expandTree(nivel, 0);
        else
            processLeaf();

        treeOperator.removeSeed(nivel);

        index_seed++;
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
inline void TreeGenerator<TOperator>::expandTree(unsigned int nivel, unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    unsigned int index_angles = 0;
    typename TOperator::RMatrix R_inicial;

    treeOperator.copyMatrix(R_inicial);

    while (index_angles < CANT_ANGLES && !ultimo_nivel_exitoso)
    {
        treeOperator.initMatrix(R_inicial);

        if (treeOperator.putFirst(nivel, index_angles, indice_nivel_anterior))
        {
            if (nivel < CANT_RES + 1)
                ultimo_nivel_exitoso = appendElements(nivel, index_angles);
            else
                ultimo_nivel_exitoso = processLeaf();

            treeOperator.removeFirst(nivel);
        }

        index_angles++;
    }
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::appendElements(unsigned int nivel, unsigned int index)
{
    typename TOperator::KeepRecursion resultRecursion;
    bool result = false;
    unsigned int index_res = 0;

    while (treeOperator.putNext(nivel, index_res, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (nivel < CANT_RES + 1)
                expandTree(nivel, index);
            else
                result = processLeaf();

            treeOperator.remove(nivel);
        }

        index_res++;
    }

    return result;
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::processLeaf()
{
    return treeOperator.write();
}

