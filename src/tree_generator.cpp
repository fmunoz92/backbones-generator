#include "tree_generator.h"

SimpleTreeOperator::SimpleTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const) :
    tree_data(t)
{}

bool SimpleTreeOperator::putFirstWithSeed(unsigned int& nivel, unsigned int c)
{
    bool result = false;

    if (c == 0) //TODO: use mili::FirstTimeFlag
    {
        Residuo residuo;
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R, residuo);
        mili::insert_into(paraBorrar, residuo);
        nivel = 2;
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
    R = newR;
}

bool SimpleTreeOperator::putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, Result& resultRecursion, unsigned int c)
{
    bool result = false;
    resultRecursion = stopRecursion;

    if (c == 0) //TODO: use mili::FirstTimeFlag
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
}

bool ChainsTreeOperator::putFirstWithSeed(unsigned int& nivel, unsigned int c)
{
    bool result = true;
    AnglesData* chain;
    Residuo residuo;
    vector<Residuo> residuos;
    if ((chain = reader->read(c)) != NULL)
    {
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R, residuo);
        addChain(R, 2, tree_data, residuos, *chain, c);
        nivel = residuos.size() + 2;
        residuosParaBorrar.push_back(residuo);
        vectoresParaBorrar.push_back(residuos);
    }
    else
        result = false;

    return result;
}

bool ChainsTreeOperator::putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& recursion, unsigned int c)
{
    bool result = true;
    FilterResultType filterResult;
    AnglesData* chain;
    recursion = stopRecursion;
    Residuo residuo;
    vector<Residuo> residuos;

    if (c == 0)//TODO: use mili::FirstTimeFlag and modularize
    {
        if (poneres(R, nivel, tree_data, residuo, indice_nivel_anterior, i) == FILTER_OK)
        {
            residuosParaBorrar.push_back(residuo);
            nivel++;
        }
        else
            result = false;
    }

    if (result && (chain = reader->read(c)) != NULL) //TODO: instead c using class attribute
    {
        filterResult = addChain(R, nivel, tree_data, residuos, *chain, c);
        nivel += residuos.size();
        if (filterResult == FILTER_OK)
        {
            vectoresParaBorrar.push_back(residuos);
            recursion = doRecursion;
        }
        else
            //saco residuos apendeados antes del primer residuo que genero FILTER_FAIL
            TreeHelper::sacar_residuos(tree_data, residuos);
    }
    else
        result = false;

    return result;
}

void ChainsTreeOperator::remove()
{
    if (!residuosParaBorrar.empty())
    {
        TreeHelper::sacar_residuo(tree_data, residuosParaBorrar.back());
        residuosParaBorrar.pop_back();
    }
    if (!vectoresParaBorrar.empty())
    {
        TreeHelper::sacar_residuos(tree_data, vectoresParaBorrar.back());
        vectoresParaBorrar.pop_back();
    }
    tree_data.fragment_ids.pop_back();
}
