#ifndef TREE_OPERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "backbones-generator/filer.h"
#include "backbones-generator/tree_filters.h"

template <class WriterHelper>
inline TreeOperator<WriterHelper>::TreeOperator(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader)
    : tree_helper(tree_helper),
      writer_helper(tree_helper, reader)//call constructor adapter
{
}

template <class WriterHelper>
inline TreeOperator<WriterHelper>::TreeOperator(TreeHelper& tree_helper)
    : tree_helper(tree_helper),
      writer_helper(tree_helper)
{
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::copyMatrix(RMatrix R_inicial) const
{
    backbones_utils::copymat(R_inicial, R);
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::initMatrix(const RMatrix newR)
{
    backbones_utils::copymat(R, newR);
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::putFirst(unsigned int& nivel, const unsigned int fi_index, const unsigned int si_index)
{
    Residuo residuo;

    const bool result = tree_helper.putRes(R, nivel, residuo, si_index, fi_index) == TreeFilters::FILTER_OK;

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
inline bool TreeOperator<WriterHelper>::write()
{
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    const bool exito = true;
#else
    const bool exito = tree_helper.filterLastLevelOk();
#endif
    if (exito)
    {
        writer_helper.write();
        tree_helper.reportSuccess();
    }

    return exito;
}

/**********************************************************************/

template <class WriterHelper>
inline SimpleTreeOperator<WriterHelper>::SimpleTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader*)
    : TreeOperator<WriterHelper>(t)
{
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel, const unsigned int index_seed)
{
    const bool result = !this->tree_helper.success() && index_seed == 0;

    if (result)
    {
        Residuo residuo;
        this->tree_helper.clearatm();
        this->tree_helper.putSeed(this->R, residuo);
        this->residuos.push_back(residuo);
        nivel = 2; //semilla is level 1 then next level is 2
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::remove(unsigned int&)
{}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::removeSeed(unsigned int& nivel)
{
    this->removeFirst(nivel);
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& /*nivel*/, const unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion)
{
    const bool result = index_res == 0;

    if (result)
        resultRecursion = TreeOperator<WriterHelper>::DoRecursion;

    return result;
}

/**********************************************************************/

template <class WriterHelper>
inline ChainsTreeOperator<WriterHelper>::ChainsTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader* const reader)
    : TreeOperator<WriterHelper>(t),
      reader(reader)
{
}

/*Adapter*/
template <>
inline ChainsTreeOperator<FragmentsWriterHelper>::ChainsTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader* const reader)
    : TreeOperator<FragmentsWriterHelper>(t, reader),
      reader(reader)
{
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::putSeed(prot_filer::AnglesData& chain, unsigned int& nivel, const unsigned int index_seed)
{
    Residuo residuo;
    std::list<Residuo> residuos;

    const unsigned int nextLevel = 2; //semilla is "level 1"

    this->tree_helper.clearatm();

    this->tree_helper.putSeed(this->R, residuo);
    this->tree_helper.putChain(this->R, nextLevel, residuos, chain, index_seed);

    nivel = residuos.size() + nextLevel;

    this->residuos.push_back(residuo);
    vectoresParaBorrar.push_back(residuos);
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel, const unsigned int index_seed)
{
    prot_filer::AnglesData* chain = reader->read(index_seed);

    const bool result = chain != NULL;

    if (result)
        putSeed(*chain, nivel, index_seed);

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::putChain(prot_filer::AnglesData& chain, unsigned int& nivel, const unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    std::list<Residuo> residuos;

    const bool isOk = this->tree_helper.putChain(this->R, nivel, residuos, chain, index_res) == TreeFilters::FILTER_OK;

    if (isOk)
    {
        nivel += residuos.size();
        vectoresParaBorrar.push_back(residuos);
        recursion = TreeOperator<WriterHelper>::DoRecursion;
    }
    else
    {
        this->tree_helper.deleteRes(residuos);//saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
        recursion = TreeOperator<WriterHelper>::StopRecursion;
    }
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& nivel, const unsigned int index_res, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    prot_filer::AnglesData* chain = reader->read(index_res);

    const bool result = chain != NULL;//vamos a ciclar mientras tengamos chains para leer

    if (result)
        putChain(*chain, nivel, index_res, recursion);

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& nivel)
{
    const unsigned int nivelesRetrocedidos = vectoresParaBorrar.back().size();

    this->tree_helper.deleteRes(vectoresParaBorrar.back());
    this->tree_helper.deleteLastFragmentId();

    nivel -= nivelesRetrocedidos;
    vectoresParaBorrar.pop_back();
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::removeSeed(unsigned int& nivel)
{
    this->removeFirst(nivel);
    remove(nivel);
}

