#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "poneres.h"
#include "utils.h"
#include "filer.h"

template <class TOperator>
inline TreeGenerator<TOperator>::TreeGenerator(TreeData& tree_data, FullCachedAnglesSeqReader* const reader) :
    tree_data(tree_data),
    treeOperator(tree_data, reader)
{}

template <class TOperator>
inline void TreeGenerator<TOperator>::generate()
{
    float R_inicial[16];
    unsigned int nivel = 1;

    treeOperator.initMatrix(R_inicial);
    while (treeOperator.putNextSeed(nivel))
    {
        if (nivel < tree_data.nres)
            expand_tree(nivel, R_inicial, 0);
        else
            process_leaf();
        treeOperator.remove(nivel);
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
inline void TreeGenerator<TOperator>::expand_tree(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    float R_local[16];

    unsigned int index_angles = 0;
    const unsigned int angles = tree_data.cossi.size();

    while (index_angles < angles && !ultimo_nivel_exitoso)
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

    while (treeOperator.putNext(nivel, index, indice_nivel_anterior, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (nivel - 1 < tree_data.nres)
                expand_tree(nivel, R_local, index);
            else
                result = process_leaf();
            treeOperator.remove(nivel);
        }
    }

    return result;
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::process_leaf()
{
    bool exito = false;
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    treeOperator.write();
    tree_data.cont++;
#else
    if (TreeHelper::filtros_ultimo_nivel(tree_data) == FILTER_OK)
    {
        treeOperator.write();
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
#endif

    return exito;
}

template <class WriterHelper>
inline SimpleTreeOperator<WriterHelper>::SimpleTreeOperator(TreeData& t, FullCachedAnglesSeqReader*) :
    tree_data(t),
    writer_helper(t)
{}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel)
{
    bool result = false;

    if (firstTime)
    {
        Residuo residuo;
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R, residuo);
        paraBorrar.push_back(residuo);
        nivel = 2; //semilla is level 1 then next level is 2
        result = true;
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    TreeHelper::sacar_residuo(tree_data, paraBorrar.back());
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
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion)
{
    bool result = false;
    resultRecursion = StopRecursion;

    if (firstTime)
    {
        Residuo residuo;
        FilterResultType filerResult = TreeHelper::poner_residuo(R, nivel, tree_data, residuo, si_index, fi_index);

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
}

template <class WriterHelper>
inline ChainsTreeOperator<WriterHelper>::ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const reader) :
    tree_data(t),
    reader(reader),
    writer_helper(t)
{}

/*Adapter*/
template <>
inline ChainsTreeOperator<FragmentsWriterHelper>::ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const reader) :
    tree_data(t),
    reader(reader),
    writer_helper(t, reader)//call constructor adapter
{}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::initMatrix(float newR[16])
{
    R = newR;
    firstTime.reset();
    currentPosInChain = 0;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel)
{
    bool result = true;
    prot_filer::AnglesData* chain;
    Residuo residuo;
    list<Residuo> residuos;
    
    chain = reader->read(currentPosInChain);
    
    if (chain != NULL)
    {
        const unsigned int nextLevel = 2; //semilla is "level 1"
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R, residuo);
        TreeHelper::addChain(R, nextLevel, tree_data, residuos, *chain, currentPosInChain);

        nivel = residuos.size() + nextLevel;
        currentPosInChain++;
        residuosParaBorrar.push_back(residuo);
        vectoresParaBorrar.push_back(residuos);
    }
    else
        result = false;

    return result;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior,  KeepRecursion& recursion)
{
    bool result = true;
    FilterResultType filterResult;
    prot_filer::AnglesData* chain;
    recursion = StopRecursion;
    Residuo residuo;
    list<Residuo> residuos;

    if (firstTime)//or currentPosInChain == 0
    {
        if (TreeHelper::poner_residuo(R, nivel, tree_data, residuo, indice_nivel_anterior, i) == FILTER_OK)
        {
            residuosParaBorrar.push_back(residuo);
            nivel++;
        }
        else
            result = false;
    }

    chain = reader->read(currentPosInChain);

    if (result && (chain != NULL))
    {
        filterResult = TreeHelper::addChain(R, nivel, tree_data, residuos, *chain, currentPosInChain);
        nivel += residuos.size();
        if (filterResult == FILTER_OK)
        {
            vectoresParaBorrar.push_back(residuos);
            recursion = DoRecursion;
        }
        else
            //saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
	            TreeHelper::sacar_residuos(tree_data, residuos);
        currentPosInChain++;
    }

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    const unsigned int nivelesRetrocedidos = vectoresParaBorrar.back().size() + 1;// +1 for residuosParaBorrar
    TreeHelper::sacar_residuo(tree_data, residuosParaBorrar.back());
    TreeHelper::sacar_residuos(tree_data, vectoresParaBorrar.back());

    residuosParaBorrar.pop_back();
    vectoresParaBorrar.pop_back();
    tree_data.fragment_ids.pop_back();

    nivel -= nivelesRetrocedidos;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::write()
{
    writer_helper.write();
}
