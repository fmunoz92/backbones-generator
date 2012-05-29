#ifndef TREE_GENERATOR_H
#define TREE_GENERATOR_H

#include <string>
#include <vector>
#include "poneres.h"
#include "utils.h"

const string output_f = "traj.xtc";

typedef enum {doRecursion, stopRecursion} Result;

template <class TreeOperator, class WriterHelper>
class TreeGenerator
{
public:
    TreeGenerator(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
    ~TreeGenerator() {}
    void generate();
private:
    void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
    TreeData& tree_data;
    WriterHelper writer_helper;
    TreeOperator g;
};

class SimpleTreeOperator
{
public:
    SimpleTreeOperator(TreeData& t, FullCachedAnglesSeqReader* const);
    ~SimpleTreeOperator()
    {}

    bool putFirstWithSeed(unsigned int& nivel);
    void initMatrix(float R[16]);
    bool putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& doRecursion);
    void remove();
private:
    mili::FirstTimeFlag firstTime;
    float* R;
    TreeData& tree_data;
    vector<Residuo> paraBorrar;
};

class ChainsTreeOperator
{
public:

    ChainsTreeOperator(TreeData& t, FullCachedAnglesSeqReader* reader);
    ~ChainsTreeOperator()
    {}

    bool putFirstWithSeed(unsigned int& nivel);
    void initMatrix(float R[16]);
    bool putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& doRecursion);
    void remove();
private:
    mili::FirstTimeFlag firstTime;
    unsigned int currentPosInChain;
    float* R;
    TreeData& tree_data;
    vector<Residuo> residuosParaBorrar;
    vector<vector<Residuo> > vectoresParaBorrar;
    FullCachedAnglesSeqReader* const reader;
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
    inline XtcWriterHelper();
    inline void write(TreeData& tree_data);
private:
    const string output_file;
    XtcWriter writer;
};

class CompressedWriterHelper
{
public:
    inline CompressedWriterHelper();
    inline void write(TreeData& tree_data);
private:
    const string output_file;
    CompressedWriter writer;
};

class FragmentsWriterHelper
{
public:
    inline FragmentsWriterHelper();
    inline void write(TreeData& tree_data);
private:
    const string output_file;
    FragmentsWriter writer;
};

template<class T>
struct Generate
{};

template<>
struct Generate<SimpleTreeOperator>
{
    void operator()(const string& format, TreeData& tree_data, FullCachedAnglesSeqReader* const)
    {
        if (format == "xtc")
        {
            TreeGenerator<SimpleTreeOperator, XtcWriterHelper> g(tree_data, NULL);
            g.generate();
        }
        else if (format == "compressed")
        {
            TreeGenerator<SimpleTreeOperator, CompressedWriterHelper> g(tree_data, NULL);
            g.generate();
        }
        else
            throw runtime_error("wrong format");
    }
};

template<>
struct Generate<ChainsTreeOperator>
{
    void operator()(const string& format, TreeData& tree_data, FullCachedAnglesSeqReader* const reader)
    {
        if (format == "xtc")
        {
            TreeGenerator<ChainsTreeOperator, XtcWriterHelper> g(tree_data, reader);
            g.generate();
        }
        else if (format == "compressed")
        {
            TreeGenerator<ChainsTreeOperator, CompressedWriterHelper> g(tree_data, reader);
            g.generate();
        }
        else if (format == "fragments")
        {
            TreeGenerator<ChainsTreeOperator, FragmentsWriterHelper> g(tree_data, reader);
            g.generate();
        }
        else
            throw runtime_error("wrong format");
    }
};

#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
