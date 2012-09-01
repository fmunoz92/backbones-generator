#include "backbones-generator/tree_helper.h"

TreeHelper::TreeHelper(TreeData& treeData, TreeFilters& treeFilters)
    : treeData(treeData),
      treeFilters(treeFilters)
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
    bool ok = treeFilters.calcRdG(treeData.atm, treeData.nres, treeData.rgmax) == TreeFilters::FILTER_OK &&
              treeFilters.volumen_en_rango(treeData.nres, treeData.grilla.obtener_vol_parcial()) == TreeFilters::FILTER_OK;

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

TreeFilters::FilterResultType TreeHelper::putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int siIndex, unsigned int fiIndex)
{
    const float cossi = treeData.cossi[siIndex];
    const float sinsi = treeData.sinsi[siIndex];
    const float cosfi = treeData.cosfi[fiIndex];
    const float sinfi = treeData.sinfi[fiIndex];

    const unsigned int i = resN - 2;

    treeData.anglesData.angles[i].si = siIndex;
    treeData.anglesData.angles[i].fi = fiIndex;

#ifdef COMBINATIONS_DEBUG
#include "prot-filer/backbones_utils.h"
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(treeData, treeFilters);
#endif

    bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, treeData.atm, resN, filter);

    if (success)
    {
        const prot_filer::ATOM& atm = treeData.atm[3 * (resN - 1) + 1];
        residuo.at2 = treeData.grilla.agregar_esfera(atm.x, atm.y, atm.z);
    }

    return success ? TreeFilters::FILTER_OK : TreeFilters::FILTER_FAIL;
}

TreeFilters::FilterResultType TreeHelper::putChain(float* pR, unsigned int resN, std::list<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chainIndex)
{
    TreeFilters::FilterResultType result = TreeFilters::FILTER_OK;
    unsigned int i = 0;

    while (result == TreeFilters::FILTER_OK  && ((i + 1) < chain.nres) && ((resN + i) < treeData.nres))
    {
        Residuo residuo;
        result = putRes(pR, resN + i, residuo, chain.angles[i].si, chain.angles[i].fi);

        if (result == TreeFilters::FILTER_OK)
            mili::insert_into(residuos, residuo);

        ++i;
    }

    if (result == TreeFilters::FILTER_OK)
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
    treeData.cont++;
    treeData.hubo_algun_exito = true;
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
    return treeData.outputFile;
}

unsigned int TreeHelper::getNRes() const
{
    return treeData.nres;
}

unsigned int TreeHelper::getNAngles() const
{
    return treeData.cossi.size();
}
