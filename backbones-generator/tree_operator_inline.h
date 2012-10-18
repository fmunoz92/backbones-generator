#ifndef TREE_OPERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "backbones-generator/filer.h"
#include "backbones-generator/tree_filters.h"

template <class WriterHelper>
inline TreeOperator<WriterHelper>::TreeOperator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* reader)
    : treeHelper(treeHelper),
      writerHelper(treeHelper, reader)//call constructor adapter
{
}

template <class WriterHelper>
inline TreeOperator<WriterHelper>::TreeOperator(TreeHelper& treeHelper)
    : treeHelper(treeHelper),
      writerHelper(treeHelper)
{
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::copyMatrix(RMatrix rInicial) const
{
    backbones_utils::copymat(rInicial, R);
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::initMatrix(const RMatrix newR)
{
    backbones_utils::copymat(R, newR);
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::putFirst(unsigned int& level, const unsigned int fiIndex, const unsigned int siIndex)
{
    Residuo residuo;

    const bool result = treeHelper.putRes(R, level, residuo, siIndex, fiIndex) == TreeFilters::FILTER_OK;

    if (result)
    {
        residuos.push_back(residuo);
        level++;
    }

    return result;
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::removeFirst(unsigned int& level)
{
    treeHelper.deleteRes(residuos.back());
    level--;
    residuos.pop_back();
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::write()
{
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    const bool success = true;
#else
    const bool success = treeHelper.filterLastLevelOk();
#endif
    if (success)
    {
        writerHelper.write();
        treeHelper.reportSuccess();
    }

    return success;
}

/**********************************************************************/

template <class WriterHelper>
inline SimpleTreeOperator<WriterHelper>::SimpleTreeOperator(TreeHelper& t, FullCachedAnglesSeqReader*)
    : TreeOperator<WriterHelper>(t)
{
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNextSeed(unsigned int& level, const unsigned int indexSeed)
{
    const bool result = !this->treeHelper.success() && indexSeed == 0;

    if (result)
    {
        Residuo residuo;
        this->treeHelper.clearatm();
        this->treeHelper.putSeed(this->R, residuo);
        this->residuos.push_back(residuo);
        level = 2; //semilla is level 1 then next level is 2
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::remove(unsigned int&)
{}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::removeSeed(unsigned int& level)
{
    this->removeFirst(level);
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& /*level*/, const unsigned int indexRes, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion)
{
    const bool result = indexRes == 0;

    resultRecursion = (result) ? TreeOperator<WriterHelper>::DoRecursion : TreeOperator<WriterHelper>::StopRecursion;

    return result;
}

/**********************************************************************/

template <class WriterHelper>
inline ChainsTreeOperator<WriterHelper>::ChainsTreeOperator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* const reader)
    : TreeOperator<WriterHelper>(treeHelper),
      reader(reader)
{
}

/*Adapter*/
template <>
inline ChainsTreeOperator<FragmentsWriterHelper>::ChainsTreeOperator(TreeHelper& treeHelper, FullCachedAnglesSeqReader* const reader)
    : TreeOperator<FragmentsWriterHelper>(treeHelper, reader),
      reader(reader)
{
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::putSeed(prot_filer::AnglesData& chain, unsigned int& level, const unsigned int indexSeed)
{
    Residuo residuo;
    std::list<Residuo> residuos;

    const unsigned int nextLevel = 2; //semilla is "level 1"

    this->treeHelper.clearatm();

    this->treeHelper.putSeed(this->R, residuo);
    this->treeHelper.putChain(this->R, nextLevel, residuos, chain, indexSeed);

    level = residuos.size() + nextLevel;

    this->residuos.push_back(residuo);
    vectoresParaBorrar.push_back(residuos);
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNextSeed(unsigned int& level, const unsigned int indexSeed)
{
    prot_filer::AnglesData* chain = reader->read(indexSeed);

    const bool result = chain != NULL;

    if (result)
        putSeed(*chain, level, indexSeed);

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::putChain(prot_filer::AnglesData& chain, unsigned int& level, const unsigned int indexRes, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    std::list<Residuo> residuos;

    const bool isOk = this->treeHelper.putChain(this->R, level, residuos, chain, indexRes) == TreeFilters::FILTER_OK;

    if (isOk)
    {
        level += residuos.size();
        vectoresParaBorrar.push_back(residuos);
        recursion = TreeOperator<WriterHelper>::DoRecursion;
    }
    else
    {
        this->treeHelper.deleteRes(residuos);//saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
        recursion = TreeOperator<WriterHelper>::StopRecursion;
    }
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& level, const unsigned int indexRes, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    prot_filer::AnglesData* chain = reader->read(indexRes);

    const bool result = chain != NULL;//vamos a ciclar mientras tengamos chains para leer

    if (result)
        putChain(*chain, level, indexRes, recursion);

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& level)
{
    const unsigned int nivelesRetrocedidos = vectoresParaBorrar.back().size();

    this->treeHelper.deleteRes(vectoresParaBorrar.back());
    this->treeHelper.deleteLastFragmentId();

    level -= nivelesRetrocedidos;
    vectoresParaBorrar.pop_back();
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::removeSeed(unsigned int& level)
{
    remove(level);
    this->removeFirst(level);
}
