#ifndef TREE_GENERATOR_H
#define TREE_GENERATOR_H

#include <string>
#include <vector>
#include <mili/mili.h>
#include "poneres.h"
#include "utils.h"

template <class TOperator>
class TreeGenerator
{
public:
    inline TreeGenerator(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
    inline void generate();
private:
    inline void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    inline bool procesar_ultimo_nivel();
    TreeData& tree_data;
    TOperator treeOperator;
};

struct TreeOperator
{
    enum KeepRecursion {DoRecursion, StopRecursion};
    /*
    interfaz no utilizada, por razones de eficacia?
    virtual bool putNextSeed(unsigned int& nivel) = 0;
    virtual void initMatrix(float R[16]) = 0;
    virtual bool putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion) = 0;
    virtual void remove() = 0;
    virtual void write() = 0;
    */
};

template <class WriterHelper>
class SimpleTreeOperator : public TreeOperator
{
public:
    inline SimpleTreeOperator(TreeData& t, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& nivel);
    inline void initMatrix(float R[16]);
    inline bool putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion);
    inline void remove();
    inline void write();
private:
    mili::FirstTimeFlag firstTime;
    float* R;
    TreeData& tree_data;
    vector<Residuo> paraBorrar;
    WriterHelper writer_helper;
};

template <class WriterHelper>
class ChainsTreeOperator : public TreeOperator
{
public:
    inline ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* reader);

    inline bool putNextSeed(unsigned int& nivel);
    inline void initMatrix(float R[16]);
    inline bool putNext(unsigned int& nivel, unsigned int fi_index, unsigned int si_index, KeepRecursion& resultRecursion);
    inline void remove();
    inline void write();
private:
    mili::FirstTimeFlag firstTime;
    unsigned int currentPosInChain;
    float* R;
    TreeData& tree_data;
    vector<Residuo> residuosParaBorrar;
    vector<vector<Residuo> > vectoresParaBorrar;
    FullCachedAnglesSeqReader* const reader;
    WriterHelper writer_helper;
};

struct TreeHelper
{
    static inline void semilla(TreeData& tree_data, float* R, Residuo& residuo);
    static inline FilterResultType filtros_ultimo_nivel(TreeData& tree_data);
    static inline void sacar_residuo(TreeData& tree_data, const Residuo& residuo);
    static inline void sacar_residuos(TreeData& tree_data, const vector<Residuo>& residuos);
};

class XtcWriterHelper
{
public:
    XtcWriterHelper(TreeData& tree_data);
    ~XtcWriterHelper();
    void write();
private:
    TreeData& tree_data;
    XtcWriter writer;
};

class CompressedWriterHelper
{
public:
    CompressedWriterHelper(TreeData& tree_data);
    ~CompressedWriterHelper();
    void write();
private:
    TreeData& tree_data;
    CompressedWriter writer;
};

class FragmentsWriterHelper
{
public:
    FragmentsWriterHelper(TreeData& tree_data, FullCachedAnglesSeqReader* reader);//Adapter
    ~FragmentsWriterHelper();
    void write();
private:
    TreeData& tree_data;
    const FullCachedAnglesSeqReader* reader;
    FragmentsWriter writer;
};

struct IGeneratorSimple
{
    virtual ~IGeneratorSimple() {}
    virtual void generate(TreeData& tree_data) = 0;
};

struct IGeneratorChains
{
    virtual ~IGeneratorChains() {}
    virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader) = 0;
};

template <class Writer>
class GeneratorChains : public IGeneratorChains
{
public:
    inline virtual ~GeneratorChains() {}
private:
    inline virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
};

template <class Writer>
class GeneratorSimple : public IGeneratorSimple
{
public:
    virtual ~GeneratorSimple() {}
private:
    inline virtual void generate(TreeData& tree_data);
};


#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
