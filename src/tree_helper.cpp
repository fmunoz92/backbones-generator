#include "backbones-generator/tree_helper.h"

#include <fstream>

TreeHelper::TreeHelper(TreeData& treeData, TreeFilters& treeFilters, const std::string&  outputFile)
    : treeData(treeData),
      treeFilters(treeFilters),
      outputFile(outputFile)
{}

void TreeHelper::putSeed(float* R, Residuo& residuo)
{
    backbones_utils::semilla(treeData.atm, R);

    treeData.anglesData.seed[0] = treeData.atm[0];
    treeData.anglesData.seed[1] = treeData.atm[1];
    treeData.anglesData.seed[2] = treeData.atm[2];

    residuo.at2 = treeData.grilla.agregar_esfera(treeData.atm[1].x, treeData.atm[1].y, treeData.atm[1].z);
}

bool TreeHelper::filterLastLevelOk()
{
    const bool ok = treeFilters.calcRdG(treeData.atm, treeData.nres, treeData.rgmax) == TreeFilters::FILTER_OK &&
                    treeFilters.volumeInRange(treeData.nres, treeData.grilla.obtener_vol_parcial()) == TreeFilters::FILTER_OK;

    return ok;
}

void TreeHelper::deleteRes(const Residuo& residuo)
{
    treeData.grilla.sacar_esfera(residuo.at2);
}

void TreeHelper::deleteRes(const std::list<Residuo>& residuos)
{
    for (std::list<Residuo>::const_iterator it = residuos.begin(); it != residuos.end(); it++)
        deleteRes(*it);
}

bool TreeHelper::putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int siIndex, unsigned int fiIndex)
{
    const float cossi = treeData.cossi[siIndex];
    const float sinsi = treeData.sinsi[siIndex];
    const float cosfi = treeData.cosfi[fiIndex];
    const float sinfi = treeData.sinfi[fiIndex];

    const unsigned int angle = resN - 1;//la semilla no se considera

    treeData.anglesData.angles[angle].si = siIndex;
    treeData.anglesData.angles[angle].fi = fiIndex;

#ifdef COMBINATIONS_DEBUG
#include "prot-filer/backbones_utils.h"
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(treeData, treeFilters);
#endif

    const unsigned int RES_N = resN + 1;//prot_filer comienza desde 1
    const bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, treeData.atm, RES_N, filter);

    if (success)
    {
        const prot_filer::ATOM& atm = treeData.atm[3 * (RES_N - 1) + 1];
        residuo.at2 = treeData.grilla.agregar_esfera(atm.x, atm.y, atm.z);
    }

    return success;
}

bool TreeHelper::putResOfChain(float* pR, const unsigned int resN, Residuo& residuo, unsigned int siIndex, unsigned int fiIndex, std::list<Residuo>& residuos)
{
    const bool result = putRes(pR, resN, residuo, siIndex, fiIndex);
    if (result)
        mili::insert_into(residuos, residuo);

    return result;
}

bool TreeHelper::putChain(float* pR, unsigned int resN, std::list<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chainIndex,
                          unsigned int firstSi, unsigned int firstFi)
{
    const unsigned int LENGTH_OF_CHAIN = chain.nres - 1;
    bool result;
    Residuo residuo;
    unsigned int fi;
    unsigned int si;

    //the first residue replaces the fragment's seed
    result = putResOfChain(pR, resN, residuo, firstSi, firstFi, residuos);
    ++resN;

    // iterate on the angles of the fragment (angles between residues).
    unsigned int angle = 1;//we start from the 2nd pair of angles, since the
    // first pair is the angle between the seed and the 1st residue. We don't
    // consider the seed from the fragments.
    while (result  && (angle < LENGTH_OF_CHAIN) && (resN < treeData.nres))
    {
        fi = chain.angles[angle].fi;
        si = chain.angles[angle].si;

        result = putResOfChain(pR, resN, residuo, si, fi, residuos);

        ++resN;
        ++angle;
    }

    if (result)
        treeData.fragmentIds.push_back(chainIndex);

    return result;
}

void TreeHelper::clearatm()
{
    for (unsigned int i = 0; i < 3 * treeData.nres; i++)
    {
        treeData.atm[i].x = 0.0;
        treeData.atm[i].y = 0.0;
        treeData.atm[i].z = 0.0;
    }
}

bool TreeHelper::success() const
{
    return treeData.hubo_algun_exito;
}

void TreeHelper::reportSuccess()
{
    ++treeData.cont;
#ifndef COMBINATIONS_DEBUG
    treeData.hubo_algun_exito = true;
#endif
}

void TreeHelper::deleteLastFragmentId()
{
    treeData.fragmentIds.pop_back();
}

const prot_filer::AnglesData&  TreeHelper::getAnglesData() const
{
    return treeData.anglesData;
}

Atoms& TreeHelper::getAtm()
{
    return treeData.atm;
}

const prot_filer::FragmentIds& TreeHelper::getFragmentIds() const
{
    return treeData.fragmentIds;
}

const std::string& TreeHelper::getOutputFile() const
{
    return outputFile;
}

unsigned int TreeHelper::getNRes() const
{
    return treeData.nres;
}

unsigned int TreeHelper::getNAngles() const
{
    return treeData.cossi.size();
}
