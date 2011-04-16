#include "tree_generator.h"
#include "poneres.h"
#include "utils.h"

WriterAdapter* TreeGenerator::createWriter(const string& write_format)
{
    if (write_format == "xtc")
    {
        return new XtcWriterAdapter();
    }
    else if (write_format == "compressed")
    {
        return new CompressedWriterAdapter();
    }
    else
    {
        throw runtime_error("Invalid format");
    }
}

void TreeGenerator::semilla(float* R, Residuo& residuo)
{
    ATOM* atm = tree_data.atm;
    backbones_utils::semilla(atm, R);
    ATOM* seed = tree_data.angles_data->seed;

    seed[0] = atm[0];
    seed[1] = atm[1];
    seed[2] = atm[2];

    residuo.at2 = tree_data.grilla->agregar_esfera(tree_data.atm[1].x, tree_data.atm[1].y, tree_data.atm[1].z);
}

FilterResultType TreeGenerator::filtros_ultimo_nivel()
{
    bool ok = calcRdG(tree_data.atm, tree_data.nres, tree_data.rgmax) == FILTER_OK
              && volumen_en_rango(tree_data.nres, tree_data.grilla->obtener_vol_parcial()) == FILTER_OK;
    return ok ? FILTER_OK : FILTER_FAIL;
}

void SimpleTreeGenerator::generate()
{
    float R_inicial[16];
    // Va a ser el residuo que agregue semilla en cada iteracion y al terminar
    // del ciclo se usa para sacar el residuo del grillado.
    Residuo residuo;
    unsigned int i = 0;
    while (i < tree_data.cossi.size() && !tree_data.hubo_algun_exito)
    {
        clearatm(tree_data.atm, tree_data.nres);
        semilla(R_inicial, residuo);

        generar_nivel_intermedio(2, R_inicial, i);
        sacar_residuo(residuo);
        ++i;
    }
}

void SimpleTreeGenerator::generar_nivel_intermedio(const unsigned int nivel, const float R_inicial[16], const unsigned int indice_nivel_anterior)
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
            sacar_residuo(residuo);
        }
        ++i;
    }
}

// True si el volumen indicado por el grillado se encuentra en el rango aceptable
bool SimpleTreeGenerator::procesar_ultimo_nivel()
{
#ifdef COMBINATIONS_DEBUG
    // En el modo DEBUG se deshabilitan los chequeos.
    writer->write(tree_data.atm, *tree_data.angles_data);
    tree_data.cont++;
    return false;
#else
    bool exito = false;

    if (filtros_ultimo_nivel() == FILTER_OK)
    {
        writer->write(tree_data.atm, *tree_data.angles_data);
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
    return exito;
#endif
}

void ChainsTreeGenerator::generate()
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
        semilla(R_inicial, residuo);
        addChain(R_inicial, 2, tree_data, residuos, *chain, i);
        const unsigned int nivel = residuos.size() + 2;
        if (nivel < tree_data.nres)
        {
            generar_nivel_intermedio(nivel, R_inicial, i);
        }
        else
        {
            procesar_ultimo_nivel();
        }
        sacar_residuos(residuos);
        residuos.clear();
        sacar_residuo(residuo);
        tree_data.chain_indexs.pop_back();
        ++i;
    }
}

void ChainsTreeGenerator::generar_nivel_intermedio(const unsigned int nivel, const float R_inicial[16], const unsigned int indice_nivel_anterior)
{
    float R_local[16];
    FilterResultType resultado;
    bool exito = false; // solo aplicable si somos anteultimo nivel
    assert(nivel > 1);  // pre condicion
    Residuo residuo;
    vector<Residuo> residuos;

    unsigned int i = 0;
    IncompleteAnglesData* chain;
    while (i < tree_data.cossi.size())
    {
        backbones_utils::copymat(R_local, R_inicial);

        resultado = poneres(R_local, nivel, tree_data, residuo, indice_nivel_anterior, i);

        if (resultado == FILTER_OK)
        {
            unsigned int c = 0;
            while ((chain = reader->read(c)) != NULL)
            {
                addChain(R_local, nivel, tree_data, residuos, *chain, c);
                unsigned int next_lvl = nivel + residuos.size() + 1;
                if (next_lvl < tree_data.nres)
                {
                    generar_nivel_intermedio(nivel + residuos.size() + 1, R_local, i);
                }
                else
                {
                    exito = procesar_ultimo_nivel();
                }
                sacar_residuos(residuos);
                tree_data.chain_indexs.pop_back();
                residuos.clear();
                ++c;
            }
            sacar_residuo(residuo);
        }
        ++i;
    }
}

// True si el volumen indicado por el grillado se encuentra en el rango aceptable
bool ChainsTreeGenerator::procesar_ultimo_nivel()
{
    //TODO: remover
    //cout << tree_data.chain_indexs.size() << "\n";
    //cout << tree_data.chain_indexs[0] << "\n";
    //cout << "***********************************    \n";

#ifdef COMBINATIONS_DEBUG
// En el modo DEBUG se deshabilitan los chequeos.
    writer->write(tree_data.atm, *tree_data.angles_data);
    tree_data.cont++;
    return false;
#else
    bool exito = false;

    if (filtros_ultimo_nivel() == FILTER_OK)
    {
        writer->write(tree_data.atm, *tree_data.angles_data);
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
    return exito;
#endif
}

