#ifndef TREE_GENERATOR_INLINE_H
#error Internal header file, DO NOT include this.
#endif

template <class Writer>
inline void GeneratorSimple<Writer>::generate(TreeData& tree_data)
{
    TreeGenerator<SimpleTreeOperator<Writer> > generator(tree_data, NULL);
    generator.generate();
}

template <class Writer>
inline void GeneratorChains<Writer>::generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader)
{
    TreeGenerator<ChainsTreeOperator<Writer> > generator(tree_data, reader);
    generator.generate();
}

template <class TOperator>
inline TreeGenerator<TOperator>::TreeGenerator(TreeData& tree_data, FullCachedAnglesSeqReader* const reader) :
    tree_data(tree_data),
    treeOperator(tree_data, reader)
{}

template <class TOperator>
inline void TreeGenerator<TOperator>::generate()
{
    float R_inicial[16];
    unsigned int nivel = 1;

    treeOperator.initMatrix(R_inicial);
    while (treeOperator.putNextSeed(nivel))
    {
        if (nivel < tree_data.nres)
            generar_nivel_intermedio(nivel, R_inicial, 0);
        else
            procesar_ultimo_nivel();
        treeOperator.remove();
        nivel = 0;
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
template <class TOperator>
inline void TreeGenerator<TOperator>::generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior)
{
    bool ultimo_nivel_exitoso = false;//solo interesa si somos el anteultimo nivel
    float R_local[16];
    typename TOperator::KeepRecursion result;
    unsigned int i = 0;
    unsigned int nivelAux = nivel;
    const unsigned int angles = tree_data.cossi.size();

    while (i < angles && !ultimo_nivel_exitoso)
    {
        backbones_utils::copymat(R_local, R_inicial);
        treeOperator.initMatrix(R_local);

        while (treeOperator.putNext(nivelAux, i, indice_nivel_anterior, result))
        {
            if (result == TOperator::DoRecursion)
            {
                // hay inconsistencia en implementacion vieja, el simple le resta uno
                // a lo que seria nivel aux y el de chains lo deja como esta
                if (nivelAux - 1 < tree_data.nres)
                    generar_nivel_intermedio(nivelAux, R_local, i);
                else
                    ultimo_nivel_exitoso = procesar_ultimo_nivel();
                treeOperator.remove();
            }
        }
        i++;
        nivelAux = nivel;
    }
}

template <class TOperator>
inline bool TreeGenerator<TOperator>::procesar_ultimo_nivel()
{
    bool exito = false;
#ifdef COMBINATIONS_DEBUG // En el modo DEBUG se deshabilitan los chequeos.
    treeOperator.write();
    tree_data.cont++;
#else
    if (TreeHelper::filtros_ultimo_nivel(tree_data) == FILTER_OK)
    {
        treeOperator.write();
        tree_data.cont++;
        tree_data.hubo_algun_exito = exito = true;
    }
#endif

    return exito;
}

template <class WriterHelper>
inline SimpleTreeOperator<WriterHelper>::SimpleTreeOperator(TreeData& t, FullCachedAnglesSeqReader*) :
    tree_data(t),
    writer_helper(t)
{}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel)
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

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::remove()
{
    TreeHelper::sacar_residuo(tree_data, paraBorrar.back());
    paraBorrar.pop_back();
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::initMatrix(float newR[16])
{
    firstTime.reset();
    R = newR;
}

template <class WriterHelper>
inline bool SimpleTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion)
{
    bool result = false;
    resultRecursion = StopRecursion;

    if (firstTime)
    {
        Residuo residuo;
        FilterResultType filerResult = poneres(R, nivel, tree_data, residuo, si_index, fi_index);

        if (filerResult == FILTER_OK)
        {
            result = true;
            mili::insert_into(paraBorrar, residuo);
            resultRecursion = DoRecursion;
            nivel++;
        }
    }

    return result;
}

template <class WriterHelper>
inline void SimpleTreeOperator<WriterHelper>::write()
{
    writer_helper.write();
}

template <class WriterHelper>
inline ChainsTreeOperator<WriterHelper>::ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const reader) :
    tree_data(t),
    reader(reader),
    writer_helper(t)
{}

/*Adapter*/
template <>
inline ChainsTreeOperator<FragmentsWriterHelper>::ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const reader) :
    tree_data(t),
    reader(reader),
    writer_helper(t, reader)//call constructor adapter
{}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::initMatrix(float newR[16])
{
    R = newR;
    firstTime.reset();
    currentPosInChain = 0;
}

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNextSeed(unsigned int& nivel)
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

template <class WriterHelper>
inline bool ChainsTreeOperator<WriterHelper>::putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior,  KeepRecursion& recursion)
{
    bool result = true;
    FilterResultType filterResult;
    AnglesData* chain;
    recursion = StopRecursion;
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
            recursion = DoRecursion;
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

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::remove()
{
    TreeHelper::sacar_residuo(tree_data, residuosParaBorrar.back());
    TreeHelper::sacar_residuos(tree_data, vectoresParaBorrar.back());

    residuosParaBorrar.pop_back();
    vectoresParaBorrar.pop_back();
    tree_data.fragment_ids.pop_back();
}

template <class WriterHelper>
inline void ChainsTreeOperator<WriterHelper>::write()
{
    writer_helper.write();
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

