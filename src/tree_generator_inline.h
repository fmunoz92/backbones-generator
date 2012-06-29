#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "tree_filters.h"
#include "filer.h"


template <class TOperator>
inline TreeGenerator<TOperator>::TreeGenerator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader) :
    treeOperator(tree_helper, reader),
    CANT_RES(tree_helper.getNRes()),
    CANT_ANGLES(tree_helper.getNAngles())
{}

template <class TOperator>
inline void TreeGenerator<TOperator>::generate()
{
    float R_inicial[16];
    unsigned int nivel = 1;
    unsigned int index_seed = 0;
    treeOperator.initMatrix(R_inicial);
    while (treeOperator.putNextSeed(nivel, index_seed))
    {
        if (nivel < CANT_RES)
            expandTree(nivel, R_inicial, 0);
        else
            processLeaf();
        treeOperator.remove(nivel);
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
template <class TOperator>
inline void TreeGenerator<TOperator>::expandTree(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    float R_local[16];

    unsigned int index_angles = 0;

    while (index_angles < CANT_ANGLES && !ultimo_nivel_exitoso)
    {
        backbones_utils::copymat(R_local, R_inicial);
        treeOperator.initMatrix(R_local);

        ultimo_nivel_exitoso = appendElements(nivel, indice_nivel_anterior, index_angles, R_local);

        index_angles++;
    }
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::appendElements(unsigned int nivel, unsigned int indice_nivel_anterior, unsigned int index, const float R_local[16])
{
    typename TOperator::KeepRecursion resultRecursion;
    bool result = false;
    unsigned int index_res = 0;
    while (treeOperator.putNext(nivel, index_res, index, indice_nivel_anterior, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (nivel < CANT_RES + 1)//with ChainsTreeOperator the condition is (nivel < CANT_RES)
                expandTree(nivel, R_local, index);
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
    bool exito = false;
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    treeOperator.write();
#else
    if (treeOperator.lastLevelOk())
    {
        treeOperator.write();
        exito = true;
    }
#endif

    return exito;
}

template <class WriterHelper>
inline SimpleTreeOperator<WriterHelper>::SimpleTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader*) :
    tree_helper(t),
    writer_helper(t)
{}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel, unsigned int)
{
    bool result = false;

    if (firstTime)
    {
        Residuo residuo;
        tree_helper.clearatm();
        tree_helper.putSeed(R, residuo);
        paraBorrar.push_back(residuo);
        nivel = 2; //semilla is level 1 then next level is 2
        result = true;
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    tree_helper.deleteRes(paraBorrar.back());
    nivel--;
    paraBorrar.pop_back();
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::initMatrix(float newR[16])
{
    firstTime.reset();
    R = newR;
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::lastLevelOk()
{
    return (tree_helper.filtros_ultimo_nivel() == FILTER_OK);
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion)
{
    bool result = false;
    resultRecursion = StopRecursion;

    if (firstTime)
    {
        Residuo residuo;
        FilterResultType filerResult = tree_helper.putRes(R, nivel, residuo, si_index, fi_index);

        if (filerResult == FILTER_OK)
        {
            result = true;
            paraBorrar.push_back(residuo);
            resultRecursion = DoRecursion;
            nivel++;
        }
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::write()
{
    writer_helper.write();
    tree_helper.reportSuccess();
}

template <class WriterHelper>
inline ChainsTreeOperator<WriterHelper>::ChainsTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader* const reader) :
    tree_helper(t),
    reader(reader),
    writer_helper(t)
{}

/*Adapter*/
template <>
inline ChainsTreeOperator<FragmentsWriterHelper>::ChainsTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader* const reader) :
    tree_helper(t),
    reader(reader),
    writer_helper(t, reader)//call constructor adapter
{}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::initMatrix(float newR[16])
{
    R = newR;
    firstTime.reset();
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::lastLevelOk()
{
    return (tree_helper.filtros_ultimo_nivel() == FILTER_OK);
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel, unsigned int index_seed)
{
    bool result = true;
    prot_filer::AnglesData* chain;
    Residuo residuo;
    std::list<Residuo> residuos;

    chain = reader->read(index_seed);

    if (chain != NULL)
    {
        const unsigned int nextLevel = 2; //semilla is "level 1"
        tree_helper.clearatm();
        tree_helper.putSeed(R, residuo);
        tree_helper.putChain(R, nextLevel, residuos, *chain, index_seed);
        nivel = residuos.size() + nextLevel;
        residuosParaBorrar.push_back(residuo);
        vectoresParaBorrar.push_back(residuos);
    }
    else
        result = false;

    return result;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int index_res, unsigned int  i, unsigned int  indice_nivel_anterior,  KeepRecursion& recursion)
{
    bool result = true;
    FilterResultType filterResult;
    prot_filer::AnglesData* chain;
    recursion = StopRecursion;
    Residuo residuo;
    std::list<Residuo> residuos;

    if (firstTime)//or currentPosInChain == 0
    {
        if (tree_helper.putRes(R, nivel, residuo, indice_nivel_anterior, i) == FILTER_OK)
        {
            residuosParaBorrar.push_back(residuo);
            nivel++;
        }
        else
            result = false;
    }

    if (result && !(nivel < tree_helper.getNRes()))
    {
        recursion = DoRecursion;
        return true;
    }

    chain = reader->read(index_res);
    if (result && (chain != NULL))
    {
        filterResult = tree_helper.putChain(R, nivel, residuos, *chain, index_res);

        if (filterResult == FILTER_OK)
        {
            nivel += residuos.size();
            vectoresParaBorrar.push_back(residuos);
            recursion = DoRecursion;
        }
        else
            //saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
            tree_helper.deleteRes(residuos);
    }
    else
        result = false;

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    const unsigned int nivelesRetrocedidos = vectoresParaBorrar.back().size() + 1;

    tree_helper.deleteRes(residuosParaBorrar.back());
    residuosParaBorrar.pop_back();

    tree_helper.deleteRes(vectoresParaBorrar.back());
    vectoresParaBorrar.pop_back();

    tree_helper.deleteLastFragmentId();

    nivel -= nivelesRetrocedidos;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::write()
{
    writer_helper.write();
    tree_helper.reportSuccess();
}
