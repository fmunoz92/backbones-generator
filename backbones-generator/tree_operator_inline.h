#ifndef TREE_OPERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

#include "backbones-generator/filer.h"
#include "backbones-generator/tree_filters.h"


#define COMBINATIONS_DEBUG 1

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
inline void TreeOperator<WriterHelper>::putSeed()
{
    this->treeHelper.clearatm();
    this->treeHelper.putSeed(this->R, semilla);
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::removeSeed()
{
    treeHelper.deleteRes(semilla);
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::write()
{   
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    const bool success = false;
    writerHelper.write();
    treeHelper.reportSuccess();
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
inline void SimpleTreeOperator<WriterHelper>::remove(unsigned int& level)
{
    this->treeHelper.deleteRes(residuos.back());
    level--;
    residuos.pop_back();
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& level, unsigned int indexRes, unsigned int indexAngles, unsigned int previousLevelIndex, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion)
{
    bool result = indexRes == 0;
    if (result)
    {
        Residuo residuo;

        result = this->treeHelper.putRes(this->R, level, residuo, indexAngles, previousLevelIndex) == TreeFilters::FILTER_OK;

        if (result)
        {
            resultRecursion = TreeOperator<WriterHelper>::DoRecursion;
            residuos.push_back(residuo);
            ++level;
        }
    }

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
inline void ChainsTreeOperator<WriterHelper>::putChain(prot_filer::AnglesData& chain, unsigned int& level, const unsigned int indexRes, unsigned int firstSi, unsigned int firstFi, typename TreeOperator<WriterHelper>::KeepRecursion& recursion)
{
    ChainsRes residuos;

    const bool isOk = this->treeHelper.putChain(this->R, level, residuos, chain, indexRes, firstSi, firstFi) == TreeFilters::FILTER_OK;

    if (isOk)
    {
        level += residuos.size();
        stackChainRes.push_back(residuos);
        recursion = TreeOperator<WriterHelper>::DoRecursion;
    }
    else
    {
        this->treeHelper.deleteRes(residuos);//saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
        recursion = TreeOperator<WriterHelper>::StopRecursion;
    }
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& level, unsigned int indexRes, unsigned int indexAngles, unsigned int previousLevelIndex, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion)
{
    prot_filer::AnglesData* chain = reader->read(indexRes);

    const bool result = chain != NULL;//vamos a ciclar mientras tengamos chains para leer

    if (result)
        putChain(*chain, level, indexRes, indexAngles, previousLevelIndex, resultRecursion);

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& level)
{
    const unsigned int nivelesRetrocedidos = stackChainRes.back().size();

    this->treeHelper.deleteRes(stackChainRes.back());
    this->treeHelper.deleteLastFragmentId();

    level -= nivelesRetrocedidos;
    stackChainRes.pop_back();
}
