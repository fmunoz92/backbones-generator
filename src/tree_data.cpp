#include <cmath>

#include "backbones-generator/tree_data.h"


TreeData::TreeData(unsigned int nRes)
    : nRes(nRes),
      rgmax(2.72 * mili::cubic_root(nRes) + 5.0),              // Maximun gyration radius and maximun CA-CA distance.
      dmax2(mili::square(8.0 * mili::cubic_root(nRes) + 25.0))// Both equations constructed from database analisys.
{}

BareBackbone::BareBackbone(unsigned int nRes, Grillado& grilla, prot_filer::AnglesData& anglesData, prot_filer::AnglesMapping& anglesMapping, TreeFilters& treeFilters)
    : prot_filer::BasicProtein(nRes * 3),
      treeFilters(treeFilters),
      grilla(grilla),
      anglesMapping(anglesMapping),
      anglesData(anglesData)
{
}

TreeData* BareBackbone::treeData = NULL;

void BareBackbone::pushChainIndex(unsigned int index)
{
    fragmentIds.push_back(index);
}

void BareBackbone::deleteLastFragmentId()
{
    fragmentIds.pop_back();
}

const prot_filer::FragmentIds& BareBackbone::getFragmentIds() const
{
    return fragmentIds;
}

prot_filer::AnglesData& BareBackbone::getAnglesData()
{
    return anglesData;
}

void BareBackbone::deleteRes(const Residuo& residuo)
{
    grilla.sacar_esfera(residuo.at2);
}

void BareBackbone::deleteRes(const std::list<Residuo>& residuos)
{
    for (std::list<Residuo>::const_iterator it = residuos.begin(); it != residuos.end(); it++)
        deleteRes(*it);
}

bool BareBackbone::putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int siIndex, unsigned int fiIndex)
{
    const unsigned int angle = resN - 1;//la semilla no se considera

    anglesData.angles[angle].si = siIndex;
    anglesData.angles[angle].fi = fiIndex;

#ifdef COMBINATIONS_DEBUG
#include "prot-filer/backbones_utils.h"
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(*treeData, treeFilters);
#endif

    const float cossi = treeData->cossi[siIndex];
    const float sinsi = treeData->sinsi[siIndex];
    const float cosfi = treeData->cosfi[fiIndex];
    const float sinfi = treeData->sinfi[fiIndex];
    const unsigned int RES_N = resN + 1;//prot_filer comienza desde 1

    const bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, *this, RES_N, filter);

    if (success)
    {
        const unsigned int RES_N = resN + 1;//prot_filer comienza desde 1
        const prot_filer::ATOM& atm = (*this)[3 * (RES_N - 1) + 1];
        residuo.at2 = grilla.agregar_esfera(atm.x, atm.y, atm.z);
    }

    return success;
}

void BareBackbone::putSeed(float* R, Residuo& residuo)
{
    backbones_utils::semilla(*this, R);

    anglesData.seed[0] = (*this)[0];
    anglesData.seed[1] = (*this)[1];
    anglesData.seed[2] = (*this)[2];
    residuo.at2 = grilla.agregar_esfera((*this)[1].x, (*this)[1].y, (*this)[1].z);
}


void BareBackbone::clear()
{
    for (unsigned int i = 0; i < this->size(); ++i)//TODO: verificar si size == nres * 3de treeData
    {
        (*this)[i].x = 0.0;
        (*this)[i].y = 0.0;
        (*this)[i].z = 0.0;
    }
}


IncrementalBackbone::IncrementalBackbone(unsigned int nRes, Grillado& grilla, prot_filer::AnglesData& anglesData, prot_filer::AnglesMapping& anglesMapping, TreeFilters& treeFilters)
    : BareBackbone(nRes, grilla, anglesData, anglesMapping, treeFilters)
{
}


bool IncrementalBackbone::filterLastLevelOk()
{
    const bool ok = this->treeFilters.calcRdG(*this, this->treeData->nRes, this->treeData->rgmax) == TreeFilters::FILTER_OK &&
                    this->treeFilters.volumeInRange(this->treeData->nRes, this->grilla.obtener_vol_parcial()) == TreeFilters::FILTER_OK;

    return ok;
}

void readData(std::istream& filer, TreeData& treeData, prot_filer::AnglesMapping& anglesMapping)
{
    float fi = 0.0f;
    float si = 0.0f;
    unsigned int i = 0;

    while (filer.good())
    {
        if (filer >> fi && filer >> si)
        {
            treeData.cosfi.push_back(std::cos(mili::deg2rad(fi)));
            treeData.sinfi.push_back(std::sin(mili::deg2rad(fi)));
            treeData.cossi.push_back(std::cos(mili::deg2rad(si)));
            treeData.sinsi.push_back(std::sin(mili::deg2rad(si)));

            anglesMapping.set_mapping(fi, si);
        }
        ++i;
    }
    treeData.nAngles = treeData.cosfi.size();
    // Nota a futuro: se deberia lanzar una Excepcion si el formato del archivo fuera equivocado.
}
