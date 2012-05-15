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
    TreeGenerator(TreeData& tree_data);
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
    SimpleTreeOperator(TreeData& t);
    ~SimpleTreeOperator();
    void putFirstWithSeed(float R[16]);
    void initMatrix(float R[16]);
    bool putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& doRecursion);
    void remove();
private:
    float* R;
    TreeData& tree_data;
    vector<Residuo> paraBorrar;
    bool yaPuseUnResiduo;
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
    inline ~XtcWriterHelper();
    inline void write(TreeData& tree_data);
private:
    const string output_file;
    XtcWriter writer;
};

class CompressedWriterHelper
{
public:
    inline CompressedWriterHelper();
    inline ~CompressedWriterHelper();
    inline void write(TreeData& tree_data);
private:
    const string output_file;
    CompressedWriter writer;
};

struct Generate
{
    void operator()(const string& format, TreeData& tree_data)
    {
        if (format == "xtc")
        {
            TreeGenerator<SimpleTreeOperator, XtcWriterHelper> g(tree_data);
            g.generate();
        }
        else if (format == "compressed")
        {
            TreeGenerator<SimpleTreeOperator, CompressedWriterHelper> g(tree_data);
            g.generate();
        }
    }
};

#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
