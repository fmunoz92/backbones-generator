#ifndef TREE_GENERATOR_CPP
#define TREE_GENERATOR_CPP

#include "tree_generator.h"

SimpleTreeOperator::SimpleTreeOperator(TreeData& t) :
    tree_data(t),
    yaPuseUnResiduo(false)
{}

SimpleTreeOperator::~SimpleTreeOperator()
{}

void SimpleTreeOperator::putFirstWithSeed(float R_inicial[16])
{
    Residuo residuo;
    clearatm(tree_data.atm, tree_data.nres);
    TreeHelper::semilla(tree_data, R_inicial, residuo);
    mili::insert_into(paraBorrar, residuo);
}

void SimpleTreeOperator::remove()
{
    TreeHelper::sacar_residuo(tree_data, paraBorrar.back());
    paraBorrar.pop_back();
}

void SimpleTreeOperator::initMatrix(float newR[16])
{
    R = newR;
    yaPuseUnResiduo = false;
}

bool SimpleTreeOperator::putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, Result& resultRecursion)
{
    bool result = false;
    resultRecursion = stopRecursion;

    if (!yaPuseUnResiduo)
    {
        yaPuseUnResiduo = true;
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

#endif
