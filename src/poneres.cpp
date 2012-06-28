#include "poneres.h"

TreeHelper::TreeHelper(TreeData& tree_data, TreeFilters& tree_filters) :
    tree_data(tree_data),
    tree_filters(tree_filters)
{}

void TreeHelper::putSeed(float* R, Residuo& residuo)
{
    backbones_utils::semilla(tree_data.atm, R);

    tree_data.angles_data->seed[0] = tree_data.atm[0];
    tree_data.angles_data->seed[1] = tree_data.atm[1];
    tree_data.angles_data->seed[2] = tree_data.atm[2];

    residuo.at2 = tree_data.grilla->agregar_esfera(tree_data.atm[1].x, tree_data.atm[1].y, tree_data.atm[1].z);
}

FilterResultType TreeHelper::filtros_ultimo_nivel()
{
    bool ok = tree_filters.calcRdG(tree_data.atm, tree_data.nres, tree_data.rgmax) == FILTER_OK &&
              tree_filters.volumen_en_rango(tree_data.nres, tree_data.grilla->obtener_vol_parcial()) == FILTER_OK;
    return ok ? FILTER_OK : FILTER_FAIL;
}

void TreeHelper::deleteRes(const Residuo& residuo)
{
    tree_data.grilla->sacar_esfera(residuo.at2);
}

void TreeHelper::deleteRes(const std::list<Residuo>& residuos)
{
    for (std::list<Residuo>::const_iterator it = residuos.begin(); it != residuos.end(); it++)
        deleteRes(*it);
}

FilterResultType TreeHelper::putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int si_index, unsigned int fi_index)
{
    float cossi = tree_data.cossi[si_index];
    float sinsi = tree_data.sinsi[si_index];
    float cosfi = tree_data.cosfi[fi_index];
    float sinfi = tree_data.sinfi[fi_index];
    Atoms& patm = tree_data.atm;

    const unsigned int i = resN - 2;

    tree_data.angles_data->angles[i].si = si_index;
    tree_data.angles_data->angles[i].fi = fi_index;

#ifdef COMBINATIONS_DEBUG
#include "prot-filer/backbones_utils.h"
    backbones_utils::DummyFilter filter;
#else
    ClashFilter filter(tree_data, tree_filters);
#endif
    bool success = backbones_utils::poneres(pR, cossi, sinsi, cosfi, sinfi, patm, resN, filter);

    if (!success)
    {
        return FILTER_FAIL;
    }

    const prot_filer::ATOM& atm = patm[3 * (resN - 1) + 1];
    residuo.at2 = tree_data.grilla->agregar_esfera(atm.x, atm.y, atm.z);

    return FILTER_OK;
}

FilterResultType TreeHelper::putChain(float* pR, unsigned int resN, std::list<Residuo>& residuos, const prot_filer::AnglesData& chain, unsigned int chain_index)
{
    FilterResultType result = FILTER_OK;
    unsigned int i = 0;

    while (result == FILTER_OK  && i < (chain.nres - 1) && ((resN + i) < tree_data.nres))
    {
        Residuo residuo;
        result = putRes(pR, resN + i, residuo, chain.angles[i].si, chain.angles[i].fi);

        if (result == FILTER_OK)
            mili::insert_into(residuos, residuo);

        ++i;
    }

    if (result == FILTER_OK)
        tree_data.fragment_ids.push_back(chain_index);

    return result;
}

void TreeHelper::clearatm()
{
    for (unsigned int i = 0; i < 3 * tree_data.nres; i++)
    {
        tree_data.atm[i].x = 0.0;
        tree_data.atm[i].y = 0.0;
        tree_data.atm[i].z = 0.0;
    }
}

void TreeHelper::reportSuccess()
{
    tree_data.cont++;
    tree_data.hubo_algun_exito = true;
}

void TreeHelper::deleteLastFragmentId()
{
    tree_data.fragment_ids.pop_back();
}

prot_filer::AnglesData&  TreeHelper::getAnglesData() const
{
    return *tree_data.angles_data;
}

Atoms& TreeHelper::getAtm() const
{
    return tree_data.atm;
}

prot_filer::FragmentIds& TreeHelper::getFragmentIds() const
{
    return tree_data.fragment_ids;
}

const std::string& TreeHelper::getOutputFile() const
{
    return tree_data.output_file;
}

unsigned int TreeHelper::getNRes() const
{
    return tree_data.nres;
}

unsigned int TreeHelper::getNAngles() const
{
    return tree_data.cossi.size();
}

ClashFilter::ClashFilter(const TreeData& tree_data, const TreeFilters& tree_filters) :
    tree_data(tree_data),
    tree_filters(tree_filters)
{}

bool ClashFilter::operator()(unsigned int index, const Atoms& patm, int at) const
{
    const bool clash = tree_filters.isclash(patm, at) == FILTER_FAIL;
    return !clash && (index != 2 || (tree_filters.islong(patm, at, tree_data.dmax2) != FILTER_FAIL));
}

