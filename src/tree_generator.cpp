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
    ATOM* atm = tree_data.atm;
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

template<class Writer>
void SimpleTreeGenerator<Writer>::generate()
{
    float R_inicial[16];
    // Va a ser el residuo que agregue semilla en cada iteracion y al terminar
    // del ciclo se usa para sacar el residuo del grillado.
    Residuo residuo;
    unsigned int i = 0;
    while (i < tree_data.cossi.size() && !tree_data.hubo_algun_exito)
    {
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R_inicial, residuo);

        generar_nivel_intermedio(2, R_inicial, i);
        TreeHelper::sacar_residuo(tree_data, residuo);
        ++i;
    }
}

template<class Writer>
void SimpleTreeGenerator<Writer>::generar_nivel_intermedio(const unsigned int nivel, const float R_inicial[16], const unsigned int indice_nivel_anterior)
{
    float R_local[16];
    FilterResultType resultado;
    bool exito = false; // solo aplicable si somos anteultimo nivel
    assert(nivel > 1);  // pre condicion
    Residuo residuo;

    unsigned int i = 0;
    while (i < tree_data.cossi.size() && !exito)
    {
        backbones_utils::copymat(R_local, R_inicial);

        resultado = poneres(R_local, nivel, tree_data, residuo, indice_nivel_anterior, i);

        if (resultado == FILTER_OK)
        {
            if (nivel < tree_data.nres)
            {
                generar_nivel_intermedio(nivel + 1, R_local, i);
            }
            else
            {
                exito = procesar_ultimo_nivel();
            }
            TreeHelper::sacar_residuo(tree_data, residuo);
        }
        ++i;
    }
}

// True si el volumen indicado por el grillado se encuentra en el rango aceptable
template<class Writer>
bool SimpleTreeGenerator<Writer>::procesar_ultimo_nivel()
{
#ifdef COMBINATIONS_DEBUG
    // En el modo DEBUG se deshabilitan los chequeos.
    writer_helper.write(*this);
    tree_data.cont++;
    return false;
#else
    bool exito = false;

    if (TreeHelper::filtros_ultimo_nivel(tree_data) == FILTER_OK)
    {
        writer_helper.write(*this);
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
    return exito;
#endif
}

template<class Writer>
void ChainsTreeGenerator<Writer>::generate()
{
    float R_inicial[16];
    // Va a ser el residuo que agregue semilla en cada iteracion y al terminar del ciclo
    // se usa para sacar el residuo del grillado.
    Residuo residuo;
    vector<Residuo> residuos;

    unsigned int i = 0;
    IncompleteAnglesData* chain;
    while ((chain = reader->read(i)) != NULL)
    {
        clearatm(tree_data.atm, tree_data.nres);
        TreeHelper::semilla(tree_data, R_inicial, residuo);
        addChain(R_inicial, 2, tree_data, residuos, *chain, i);
        const unsigned int nivel = residuos.size() + 2;
        if ((nivel - 1) < tree_data.nres - 1)
        {
            generar_nivel_intermedio(nivel, R_inicial, 0);
        }
        else
        {
            procesar_ultimo_nivel();
        }
        TreeHelper::sacar_residuos(tree_data, residuos);
        residuos.clear();
        TreeHelper::sacar_residuo(tree_data, residuo);
        tree_data.fragment_ids.pop_back();
        ++i;
    }
}

template<class Writer>
void ChainsTreeGenerator<Writer>::generar_nivel_intermedio(const unsigned int nivel, const float R_inicial[16], const unsigned int indice_nivel_anterior)
{
    float R_local[16];
    FilterResultType resultado;
    bool exito = false; // solo aplicable si somos anteultimo nivel
    assert(nivel > 1);  // pre condicion
    Residuo residuo;

    unsigned int i = 0;
    while (i < tree_data.cossi.size() && !exito)
    {
        backbones_utils::copymat(R_local, R_inicial);

        resultado = poneres(R_local, nivel, tree_data, residuo, indice_nivel_anterior, i);

        if (resultado == FILTER_OK)
        {
            //comparo nivel y no nivel - 1, por el poneres
            if (nivel < (tree_data.nres - 1))
            {
                unsigned int c = 0;
                IncompleteAnglesData* chain;
                while ((chain = reader->read(c)) != NULL)
                {
                    vector<Residuo> residuos;
                    unsigned int next_lvl = nivel + 1;
                    FilterResultType r = addChain(R_local, next_lvl, tree_data, residuos, *chain, c);
                    next_lvl += residuos.size();
                    if (r == FILTER_OK)
                    {
                        if ((next_lvl - 1) < (tree_data.nres - 1))
                        {
                            generar_nivel_intermedio(next_lvl, R_local, i);
                        }
                        else
                        {
                            exito = procesar_ultimo_nivel();
                        }
                        tree_data.fragment_ids.pop_back();
                    }
                    TreeHelper::sacar_residuos(tree_data, residuos);
                    ++c;
                }
            }
            else
            {
                exito = procesar_ultimo_nivel();
            }
            TreeHelper::sacar_residuo(tree_data, residuo);
        }
        ++i;
    }
}

// True si el volumen indicado por el grillado se encuentra en el rango aceptable
template<class Writer>
bool ChainsTreeGenerator<Writer>::procesar_ultimo_nivel()
{
#ifdef COMBINATIONS_DEBUG
// En el modo DEBUG se deshabilitan los chequeos.
    writer_helper.write(*this);
    tree_data.cont++;
    return false;
#else
    bool exito = false;

    if (TreeHelper::filtros_ultimo_nivel(tree_data) == FILTER_OK)
    {
        writer_helper.write(*this);
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
    return exito;
#endif
}

template class SimpleTreeGenerator<XtcWriter>;
template class SimpleTreeGenerator<CompressedWriter>;
template class ChainsTreeGenerator<XtcWriter>;
template class ChainsTreeGenerator<CompressedWriter>;
template class ChainsTreeGenerator<FragmentsWriter>;
