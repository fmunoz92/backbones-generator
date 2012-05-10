#include "tree_generator.h"
#include "poneres.h"
#include "utils.h"

class TreeHelper
{
public:
    static void semilla(TreeData& tree_data, float* R, Residuo& residuo);
    static inline FilterResultType filtros_ultimo_nivel(TreeData& tree_data);
    static inline void sacar_residuo(TreeData& tree_data, const Residuo& residuo);
    static inline void sacar_residuos(TreeData& tree_data, const vector<Residuo>& residuos);
};

void TreeHelper::sacar_residuo(TreeData& tree_data, const Residuo& residuo)
{
    tree_data.grilla->sacar_esfera(residuo.at2);
}

void TreeHelper::sacar_residuos(TreeData& tree_data, const vector<Residuo>& residuos)
{
    for (unsigned int i = 0; i < residuos.size(); ++i)
    {
        sacar_residuo(tree_data, residuos[i]);
    }
}

void TreeHelper::semilla(TreeData& tree_data, float* R, Residuo& residuo)
{
    Atoms& atm = tree_data.atm;
    backbones_utils::semilla(atm, R);
    ATOM* seed = tree_data.angles_data->seed;

    seed[0] = atm[0];
    seed[1] = atm[1];
    seed[2] = atm[2];

    residuo.at2 = tree_data.grilla->agregar_esfera(tree_data.atm[1].x, tree_data.atm[1].y, tree_data.atm[1].z);
}

FilterResultType TreeHelper::filtros_ultimo_nivel(TreeData& tree_data)
{
    bool ok = calcRdG(tree_data.atm, tree_data.nres, tree_data.rgmax) == FILTER_OK
              && volumen_en_rango(tree_data.nres, tree_data.grilla->obtener_vol_parcial()) == FILTER_OK;
    return ok ? FILTER_OK : FILTER_FAIL;
}

TreeGenerator::TreeGenerator(TreeData& tree_data, WriterHelper& helper, TreeOperator& g) :
    tree_data(tree_data),
    writer_helper(helper),
    g(g)
{}

TreeGenerator::~TreeGenerator()
{}

void TreeGenerator::generate()
{
    float R_inicial[16];
    g.putFirstWithSeed(R_inicial);
    if (tree_data.nres > 1)
    {
        generar_nivel_intermedio(2, R_inicial, 0);
        g.remove();
    }

}

/*
    given a chain C of elements
    for every angles pair P
        for every element E in the collection
            append E in C oriented to P
            recurse
            remove E
*/
void TreeGenerator::generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    float R_local[16];
    TreeOperator::Result result;
    unsigned int nivelAux = nivel;
    unsigned int i = 0;

    while (i < tree_data.cossi.size() && !ultimo_nivel_exitoso)
    {
        nivelAux = nivel;
        backbones_utils::copymat(R_local, R_inicial);
        g.initMatrix(R_local);

        while (g.putNext(nivelAux, i, indice_nivel_anterior, result))
        {
            if (result == TreeOperator::doRecursion)
            {
                if (nivel < tree_data.nres)
                    generar_nivel_intermedio(nivelAux, R_local, i);
                else
                    ultimo_nivel_exitoso = procesar_ultimo_nivel();
                g.remove();
            }
        }
        i++;
    }
}

bool TreeGenerator::procesar_ultimo_nivel()
{
    bool exito = false;
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    writer_helper.write(*this);
    tree_data.cont++;
#else
    if (TreeHelper::filtros_ultimo_nivel(tree_data) == FILTER_OK)
    {
        writer_helper.write(tree_data);
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
#endif

    return exito;
}


SimpleTreeOperator::SimpleTreeOperator(TreeData& t) :
    tree_data(t),
    yaPuseUnResiduo(false)
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

bool SimpleTreeOperator::putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& resultRecursion)
{
    bool result = false;
    resultRecursion = stopRecursion;

    if (!yaPuseUnResiduo)
    {
        yaPuseUnResiduo = true;
        Residuo residuo;
        FilterResultType filerResult = poneres(R, nivel, tree_data, residuo, indice_nivel_anterior, i);

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
