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
    RMatrix R_inicial;
    unsigned int nivel = 1;
    unsigned int index_seed = 0;

    treeOperator.initMatrix(R_inicial);

    while (treeOperator.putNextSeed(nivel, index_seed))
    {
        if (nivel < CANT_RES + 1)
        {
            treeOperator.copyMatrix(R_inicial);//for changes generated in putNextSeed
            expandTree(nivel, R_inicial, 0);
        }
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
template <class TOperator>
inline void TreeGenerator<TOperator>::expandTree(unsigned int nivel, const RMatrix R_inicial, unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    unsigned int index_angles = 0;

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
    RMatrix rLocal;
    bool result = false;
    unsigned int index_res = 0;

    while (treeOperator.putNext(nivel, index_res, resultRecursion))
    {
        if (resultRecursion == TOperator::DoRecursion)
        {
            if (nivel < CANT_RES + 1)
            {
                treeOperator.copyMatrix(rLocal);//for changes generated in putNext and putFirst
                expandTree(nivel, rLocal, index);
            }
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

template <class WriterHelper>
inline TreeOperator<WriterHelper>::TreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader) :
    tree_helper(tree_helper),
    writer_helper(tree_helper, reader)//call constructor adapter
{}

template <class WriterHelper>
inline TreeOperator<WriterHelper>::TreeOperator(TreeHelper& tree_helper) :
    tree_helper(tree_helper),
    writer_helper(tree_helper)
{}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::copyMatrix(RMatrix R_inicial)
{
    backbones_utils::copymat(R_inicial, R);
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::initMatrix(const RMatrix newR)
{
    backbones_utils::copymat(R, newR);
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::putFirst(unsigned int& nivel, unsigned int fi_index, unsigned int si_index)
{

    Residuo residuo;

    bool result = tree_helper.putRes(R, nivel, residuo, si_index, fi_index) == FILTER_OK;

    if (result)
    {
        residuos.push_back(residuo);
        nivel++;
    }

    return result;
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::removeFirst(unsigned int& nivel)
{
    tree_helper.deleteRes(residuos.back());
    nivel--;
    residuos.pop_back();
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::lastLevelOk()
{
    return (tree_helper.filtros_ultimo_nivel() == FILTER_OK);
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::write()
{
    bool exito;

#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.

    writer_helper.write();
    tree_helper.reportSuccess();
    exito = true;

#else

    exito = lastLevelOk();
    if (exito)
    {
        writer_helper.write();
        tree_helper.reportSuccess();
    }

#endif

    return exito;
}

template <class WriterHelper>
inline SimpleTreeOperator<WriterHelper>::SimpleTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader*) :
    TreeOperator<WriterHelper>(t),
    tree_helper(t)
{}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel, unsigned int)
{
    bool result = !tree_helper.success();

    if (result)
    {
        Residuo residuo;
        tree_helper.clearatm();
        tree_helper.putSeed(this->R, residuo);
        this->residuos.push_back(residuo);
        nivel = 2; //semilla is level 1 then next level is 2
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::remove(unsigned int&)
{
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::removeSeed(unsigned int& nivel)
{
    this->removeFirst(nivel);
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& /*nivel*/, unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion)
{
    bool result = index_res == 0;

    if (result)
        resultRecursion = TreeOperator<WriterHelper>::DoRecursion;

    return result;
}

template <class WriterHelper>
inline ChainsTreeOperator<WriterHelper>::ChainsTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader* const reader) :
    TreeOperator<WriterHelper>(t),
    tree_helper(t),
    reader(reader)
{}

/*Adapter*/
template <>
inline ChainsTreeOperator<FragmentsWriterHelper>::ChainsTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader* const reader) :
    TreeOperator<FragmentsWriterHelper>(t, reader),
    tree_helper(t),
    reader(reader)
{}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel, unsigned int index_seed)
{
    Residuo residuo;
    std::list<Residuo> residuos;

    prot_filer::AnglesData* chain = reader->read(index_seed);

    bool result = chain != NULL;

    if (result)
    {
        const unsigned int nextLevel = 2; //semilla is "level 1"
        tree_helper.clearatm();

        tree_helper.putSeed(this->R, residuo);
        tree_helper.putChain(this->R, nextLevel, residuos, *chain, index_seed);

        nivel = residuos.size() + nextLevel;

        this->residuos.push_back(residuo);
        vectoresParaBorrar.push_back(residuos);
    }

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::putChain(prot_filer::AnglesData& chain, unsigned int& nivel, unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    std::list<Residuo> residuos;

    bool isOk = tree_helper.putChain(this->R, nivel, residuos, chain, index_res) == FILTER_OK;

    if (isOk)
    {
        nivel += residuos.size();
        vectoresParaBorrar.push_back(residuos);
        recursion = TreeOperator<WriterHelper>::DoRecursion;
    }
    else
    {
        tree_helper.deleteRes(residuos);//saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
        recursion = TreeOperator<WriterHelper>::StopRecursion;
    }
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    prot_filer::AnglesData* chain = reader->read(index_res);

    bool result = chain != NULL;//vamos a ciclar mientras tengamos chains para leer

    if (result)
        putChain(*chain, nivel, index_res, recursion);

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    tree_helper.deleteRes(vectoresParaBorrar.back());
    tree_helper.deleteLastFragmentId();
    nivel -= vectoresParaBorrar.back().size();

    vectoresParaBorrar.pop_back();
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::removeSeed(unsigned int& nivel)
{
    this->removeFirst(nivel);
    remove(nivel);
}
