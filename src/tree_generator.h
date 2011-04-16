#include "petu.h"
#include "filer.h"

class TreeGenerator
{
protected:
    TreeData& tree_data;
    WriterAdapter* const writer;
public:
    virtual void generate() = 0;
    TreeGenerator(TreeData& tree_data, WriterAdapter* writer) :
        tree_data(tree_data),
        writer(writer)
    {
        writer->open("traj.xtc");
    };
    virtual ~TreeGenerator()
    {
        writer->close();
        delete writer;
    };
protected:
    void semilla(float* R, Residuo& residuo);
    inline FilterResultType filtros_ultimo_nivel();
    void sacar_residuo(const Residuo& residuo)
    {
        tree_data.grilla->sacar_esfera(residuo.at2);
    }
    void sacar_residuos(const vector<Residuo>& residuos)
    {
        for (unsigned int i = 0; i < residuos.size(); ++i)
        {
            sacar_residuo(residuos[i]);
        }
    }
    static WriterAdapter* createWriter(const string& write_format);
};

class SimpleTreeGenerator : public TreeGenerator
{
public:
    SimpleTreeGenerator(TreeData& tree_data, string write_format) :
        TreeGenerator(tree_data, createWriter(write_format))
    {};
    //Just for debugging purposes:
    SimpleTreeGenerator(TreeData& tree_data, WriterAdapter* w) :
        TreeGenerator(tree_data, w)
    {};
    virtual void generate();
private:
    void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
};

class ChainsTreeGenerator : public TreeGenerator
{
private:
    FullCachedAnglesSeqReader* const reader;
public:
    ChainsTreeGenerator(TreeData& tree_data, string write_format, FullCachedAnglesSeqReader* reader_) :
        TreeGenerator(tree_data, createWriter(write_format)),
        reader(reader_)
    {};

    ~ChainsTreeGenerator()
    {
        reader->close();
        delete reader;
    }
    virtual void generate();
private:
    void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
};

