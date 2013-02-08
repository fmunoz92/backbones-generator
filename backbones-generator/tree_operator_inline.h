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
inline void TreeOperator<WriterHelper>::putSeed()
{
    this->treeHelper.getAtm().clear();
    this->treeHelper.getAtm().putSeed(this->R, semilla);
}

template <class WriterHelper>
inline void TreeOperator<WriterHelper>::removeSeed()
{
    treeHelper.getAtm().deleteRes(semilla);
}

template <class WriterHelper>
inline bool TreeOperator<WriterHelper>::write()
{
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    const bool success = false;
    writerHelper.write();
#else
    const bool success = treeHelper.getAtm().filterLastLevelOk();
#endif
    if (success)
        writerHelper.write();

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
    this->treeHelper.getAtm().deleteRes(residuos.back());
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

        result = this->treeHelper.getAtm().putRes(this->R, level, residuo, indexAngles, previousLevelIndex);

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
inline bool ChainsTreeOperator<WriterHelper>::putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int siIndex, unsigned int fiIndex, std::list<Residuo>& residuos)
{
    const bool result = this->treeHelper.getAtm().putRes(pR, resN, residuo, siIndex, fiIndex);
    if (result)
        mili::insert_into(residuos, residuo);

    return result;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putChain(float* pR,
        unsigned int resN,
        std::list<Residuo>& residuos,
        const prot_filer::AnglesData& chain,
        unsigned int chainIndex,
        unsigned int firstSi,
        unsigned int firstFi)
{
    const unsigned int LENGTH_OF_CHAIN = chain.nres - 1;
    bool result;
    Residuo residuo;
    unsigned int fi;
    unsigned int si;

    //the first residue replaces the fragment's seed
    result = putRes(pR, resN, residuo, firstSi, firstFi, residuos);
    ++resN;

    // iterate on the angles of the fragment (angles between residues).
    unsigned int angle = 1;//we start from the 2nd pair of angles, since the
    // first pair is the angle between the seed and the 1st residue. We don't
    // consider the seed from the fragments.
    while (result  && (angle < LENGTH_OF_CHAIN) && (resN < this->treeHelper.getData().nRes))
    {
        fi = chain.angles[angle].fi;
        si = chain.angles[angle].si;

        result = putRes(pR, resN, residuo, si, fi, residuos);

        ++resN;
        ++angle;
    }

    if (result)
        this->treeHelper.getAtm().pushChainIndex(chainIndex);

    return result;
}


template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& level, unsigned int indexRes, unsigned int indexAngles, unsigned int previousLevelIndex, typename TreeOperator<WriterHelper>::KeepRecursion& resultRecursion)
{
    prot_filer::AnglesData* chain = reader->read(indexRes);

    const bool result = chain != NULL;//vamos a ciclar mientras tengamos chains para leer

    if (result)
    {
        ChainsRes residuos;

        const bool isOk = putChain(this->R, level, residuos, *chain, indexRes, indexAngles, previousLevelIndex);

        if (isOk)
        {
            level += residuos.size();
            stackChainRes.push_back(residuos);
            resultRecursion = TreeOperator<WriterHelper>::DoRecursion;
        }
        else
        {
            this->treeHelper.getAtm().deleteRes(residuos);//saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
            resultRecursion = TreeOperator<WriterHelper>::StopRecursion;
        }
    }

    return result;
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove(unsigned int& level)
{
    const unsigned int nivelesRetrocedidos = stackChainRes.back().size();

    this->treeHelper.getAtm().deleteRes(stackChainRes.back());
    this->treeHelper.getAtm().deleteLastFragmentId();

    level -= nivelesRetrocedidos;
    stackChainRes.pop_back();
}
