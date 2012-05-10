#include <string>
#include <vector>
#include "petu.h"

struct WriterHelper
{
    virtual void write(TreeData& tree_data) = 0;
};

struct TreeOperator
{
    typedef enum {doRecursion, stopRecursion} Result;

    virtual void putFirstWithSeed(float R[16]) = 0;
    virtual void initMatrix(float R[16])       = 0;
    virtual void remove()                      = 0;
    virtual bool putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& doRecursion) = 0;
};

class TreeGenerator
{
public:
    TreeGenerator(TreeData& tree_data, WriterHelper& helper, TreeOperator& g);
    ~TreeGenerator();
    void generate();
private:
    void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
    TreeData& tree_data;
    WriterHelper& writer_helper;
    TreeOperator& g;
};

class XtcWriterHelper : public WriterHelper
{
public:
    inline XtcWriterHelper();
    virtual ~XtcWriterHelper();
private:
    inline virtual void write(TreeData& tree_data);

    const string output_file;
    XtcWriter writer;
};


class CompressedWriterHelper : public WriterHelper
{
public:
    inline CompressedWriterHelper();
    virtual ~CompressedWriterHelper();
private:
    inline virtual void write(TreeData& tree_data);

    const string output_file;
    CompressedWriter writer;
};

class SimpleTreeOperator : public TreeOperator
{
public:
    SimpleTreeOperator(TreeData& t);
private:
    virtual void putFirstWithSeed(float R[16]);
    virtual void initMatrix(float R[16]);
    virtual bool putNext(unsigned int& nivel, unsigned int  i, unsigned int  indice_nivel_anterior, Result& doRecursion);
    virtual void remove();

    float* R;
    TreeData& tree_data;
    vector<Residuo> paraBorrar;
    bool yaPuseUnResiduo;
};


struct Generate
{
    void operator()(const string& format, TreeData& tree_data)
    {
        WriterHelper* wh;

        if (format == "xtc")
            wh = new XtcWriterHelper();
        else if (format == "compressed")
            wh = new CompressedWriterHelper();

        SimpleTreeOperator st(tree_data);
        TreeGenerator g(tree_data, *wh, st);
        g.generate();

        delete wh;
    }
};

#define INTERNAL_FILER_H
#include "filer.h"
#undef INTERNAL_FILER_H
