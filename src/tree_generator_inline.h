#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

template <class TreeOperator, class WriterHelper>
TreeGenerator<TreeOperator, WriterHelper>::TreeGenerator(TreeData& tree_data) :
    tree_data(tree_data),
    g(tree_data)
{}

template <class TreeOperator, class WriterHelper>
void TreeGenerator<TreeOperator, WriterHelper>::generate()
{
    float R_inicial[16];
    g.putFirstWithSeed(R_inicial);
    if (tree_data.nres > 1)
    {
        generar_nivel_intermedio(2, R_inicial, 0);
    }
    g.remove();
}

/*
    given a chain C of elements
    for every angles pair P
        for every element E in the collection
            append E in C oriented to P
            recurse
            remove E
*/
template <class TreeOperator, class WriterHelper>
void TreeGenerator<TreeOperator, WriterHelper>::generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    float R_local[16];
    Result result;
    unsigned int nivelAux = nivel;
    unsigned int i = 0;

    while (i < tree_data.cossi.size() && !ultimo_nivel_exitoso)
    {
        nivelAux = nivel;
        backbones_utils::copymat(R_local, R_inicial);
        g.initMatrix(R_local);

        while (g.putNext(nivelAux, i, indice_nivel_anterior, result))
        {
            if (result == doRecursion)
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

template <class TreeOperator, class WriterHelper>
bool TreeGenerator<TreeOperator, WriterHelper>::procesar_ultimo_nivel()
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

inline XtcWriterHelper::XtcWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

inline XtcWriterHelper::~XtcWriterHelper()
{
    writer.close();
}

inline void XtcWriterHelper::write(TreeData& tree_data)
{
    writer.write(tree_data.atm, *tree_data.angles_data);
}

inline CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

inline CompressedWriterHelper::CompressedWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

inline void CompressedWriterHelper::write(TreeData& tree_data)
{
    writer.write(*tree_data.angles_data);
}

inline void TreeHelper::semilla(TreeData& tree_data, float* R, Residuo& residuo)
{
    Atoms& atm = tree_data.atm;
    backbones_utils::semilla(atm, R);
    ATOM* seed = tree_data.angles_data->seed;

    seed[0] = atm[0];
    seed[1] = atm[1];
    seed[2] = atm[2];

    residuo.at2 = tree_data.grilla->agregar_esfera(tree_data.atm[1].x, tree_data.atm[1].y, tree_data.atm[1].z);
}

inline FilterResultType TreeHelper::filtros_ultimo_nivel(TreeData& tree_data)
{
    bool ok = calcRdG(tree_data.atm, tree_data.nres, tree_data.rgmax) == FILTER_OK
              && volumen_en_rango(tree_data.nres, tree_data.grilla->obtener_vol_parcial()) == FILTER_OK;
    return ok ? FILTER_OK : FILTER_FAIL;
}

inline void TreeHelper::sacar_residuo(TreeData& tree_data, const Residuo& residuo)
{
    tree_data.grilla->sacar_esfera(residuo.at2);
}

inline void TreeHelper::sacar_residuos(TreeData& tree_data, const vector<Residuo>& residuos)
{
    for (unsigned int i = 0; i < residuos.size(); ++i)
    {
        sacar_residuo(tree_data, residuos[i]);
    }
}
