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
            //with ChainsTreeOperator the condition is (nivel < CANT_RES) but
            // with Simple is (nivel < CANT_RES + 1)
            if (nivel < CANT_RES)
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
    bool exito;
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    treeOperator.write();
    exito = true;
#else
    if (treeOperator.lastLevelOk())
    {
        treeOperator.write();
        exito = true;
    }
    else
        exito = false;
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
    bool result;

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
        else
            result = false;
    }
    else
        result = false;

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
    bool result;
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
        result = true;
    }
    else
        result = false;

    return result;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putRes(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior,  KeepRecursion& recursion)
{
    Residuo residuo;
    bool result;

    if (tree_helper.putRes(R, nivel, residuo, indice_nivel_anterior, i) == FILTER_OK)
    {
        residuosParaBorrar.push_back(residuo);
        nivel++;
        if (nivel < tree_helper.getNRes())
        {
            recursion = StopRecursion; //recursionamos cuando metemos alguna cadenita recien
            result = true;
        }
        else
        {
            recursion = DoRecursion;//recursionamos para que en realidad vaya a procesar ultimo nivel
            result = true;
        }
    }
    else
        result = false;

    return result;
}
template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putChain(unsigned int& nivel, unsigned int index_res, KeepRecursion& recursion)
{
    std::list<Residuo> residuos;
    prot_filer::AnglesData* chain;
    FilterResultType filterResult;
    bool result;

    chain = reader->read(index_res);
    if (chain != NULL)
    {
        filterResult = tree_helper.putChain(R, nivel, residuos, *chain, index_res);
        if (filterResult == FILTER_OK)
        {
            nivel += residuos.size();
            vectoresParaBorrar.push_back(residuos);
            recursion = DoRecursion;
        }
        else
        {
            //saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
            tree_helper.deleteRes(residuos);
            recursion = StopRecursion;
        }

        result = true;//vamos a ciclar mientras tengamos chains para leer
    }
    else
        result = false;

    return result;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int index_res, unsigned int  i, unsigned int  indice_nivel_anterior,  KeepRecursion& recursion)
{
    bool result;

    if (firstTime)//or index_res == 0
        result = putRes(nivel, i, indice_nivel_anterior, recursion);
    else
    {
        const unsigned int index_chain = index_res - 1;// 0 is used by putRes
        result = putChain(nivel, index_chain, recursion);
    }

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    unsigned int nivelesRetrocedidos =  0;

    if (!residuosParaBorrar.empty())
    {
        tree_helper.deleteRes(residuosParaBorrar.back());
        residuosParaBorrar.pop_back();
        nivelesRetrocedidos++;
    }

    if (!vectoresParaBorrar.empty())
    {
        tree_helper.deleteRes(vectoresParaBorrar.back());
        nivelesRetrocedidos += vectoresParaBorrar.back().size();
        vectoresParaBorrar.pop_back();
    }

    tree_helper.deleteLastFragmentId();

    nivel -= nivelesRetrocedidos;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::write()
{
    writer_helper.write();
    tree_helper.reportSuccess();
}
