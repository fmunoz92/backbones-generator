#include "tree_generator.h"

SimpleTreeOperator::SimpleTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const) :
    tree_data(t)
{}

bool SimpleTreeOperator::putNextSeed(unsigned int& nivel)
{
    bool result = false;

    if (firstTime)
    {
        Residuo residuo;
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R, residuo);
        mili::insert_into(paraBorrar, residuo);
        nivel = 2; //semilla is level 1 then next level is 2
        result = true;
    }

    return result;
}

void SimpleTreeOperator::remove()
{
    TreeHelper::sacar_residuo(tree_data, paraBorrar.back());
    paraBorrar.pop_back();
}

void SimpleTreeOperator::initMatrix(float newR[16])
{
    firstTime.reset();
    R = newR;
}

bool SimpleTreeOperator::putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, Result& resultRecursion)
{
    bool result = false;
    resultRecursion = stopRecursion;

    if (firstTime)
    {
        Residuo residuo;
        FilterResultType filerResult = poneres(R, nivel, tree_data, residuo, si_index, fi_index);

        if (filerResult == FILTER_OK)
        {
            result = true;
            mili::insert_into(paraBorrar, residuo);
            resultRecursion = doRecursion;
            nivel++;
        }
    }

    return result;
}

ChainsTreeOperator::ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const reader) :
    tree_data(t),
    reader(reader)
{}

void ChainsTreeOperator::initMatrix(float newR[16])
{
    R = newR;
    firstTime.reset();
    currentPosInChain = 0;
}

bool ChainsTreeOperator::putNextSeed(unsigned int& nivel)
{
    bool result = true;
    AnglesData* chain;
    Residuo residuo;
    vector<Residuo> residuos;
    if ((chain = reader->read(currentPosInChain)) != NULL)
    {
        const unsigned int nextLevel = 2; //semilla is "level 1"
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R, residuo);
        addChain(R, nextLevel, tree_data, residuos, *chain, currentPosInChain);

        nivel = residuos.size() + nextLevel;
        currentPosInChain++;
        residuosParaBorrar.push_back(residuo);
        vectoresParaBorrar.push_back(residuos);
    }
    else
        result = false;

    return result;
}

bool ChainsTreeOperator::putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& recursion)
{
    bool result = true;
    FilterResultType filterResult;
    AnglesData* chain;
    recursion = stopRecursion;
    Residuo residuo;
    vector<Residuo> residuos;

    if (firstTime)//or currentPosInChain == 0
    {
        if (poneres(R, nivel, tree_data, residuo, indice_nivel_anterior, i) == FILTER_OK)
        {
            residuosParaBorrar.push_back(residuo);
            nivel++;
        }
        else
            result = false;
    }

    if (result && (chain = reader->read(currentPosInChain)) != NULL)
    {
        filterResult = addChain(R, nivel, tree_data, residuos, *chain, currentPosInChain);
        nivel += residuos.size();
        if (filterResult == FILTER_OK)
        {
            vectoresParaBorrar.push_back(residuos);
            recursion = doRecursion;
        }
        else
            //saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
            TreeHelper::sacar_residuos(tree_data, residuos);
        currentPosInChain++;
    }
    else
        result = false;

    return result;
}

void ChainsTreeOperator::remove()
{
    TreeHelper::sacar_residuo(tree_data, residuosParaBorrar.back());
    TreeHelper::sacar_residuos(tree_data, vectoresParaBorrar.back());

    residuosParaBorrar.pop_back();
    vectoresParaBorrar.pop_back();
    tree_data.fragment_ids.pop_back();
}
