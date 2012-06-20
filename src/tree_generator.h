#ifndef TREE_GENERATOR_H
#define TREE_GENERATOR_H

#include <string>
#include <vector>
#include "poneres.h"
#include "utils.h"
#include <mili/mili.h>

const string output_f = "traj.xtc"; //TODO: llevar a TreeData
enum KeepRecursion {DoRecursion, StopRecursion}; //TODO: llevar a TreeOperator o TreeGenerator

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

template <class WriterHelper>
class SimpleTreeOperator
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
class ChainsTreeOperator
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

class FragmentsWriterHelper
{
public:
    inline FragmentsWriterHelper(FullCachedAnglesSeqReader* reader);//Adapter
    inline ~FragmentsWriterHelper();
    inline void write(TreeData& tree_data);
private:
    const string output_file;
    FragmentsWriter writer;
    const FullCachedAnglesSeqReader* reader;
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

class ChainsFormatGeneratorFragmentsWriter : public IGeneratorChains
{
public:
    virtual ~ChainsFormatGeneratorFragmentsWriter() {}
private:
    virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
};

class ChainsFormatGeneratorXtcWriter : public IGeneratorChains
{
public:
    virtual ~ChainsFormatGeneratorXtcWriter() {}
private:
    virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
};

class ChainsFormatGeneratorCompressedWriter : public IGeneratorChains
{
public:
    virtual ~ChainsFormatGeneratorCompressedWriter() {}
private:
    virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
};

class SimpleFormatGeneratorCompressedWriter : public IGeneratorSimple
{
public:
    virtual ~SimpleFormatGeneratorCompressedWriter() {}
private:
    virtual void generate(TreeData& tree_data);
};

class SimpleFormatGeneratorXtcWriter : public IGeneratorSimple
{
public:
    virtual ~SimpleFormatGeneratorXtcWriter() {}
private:
    virtual void generate(TreeData& tree_data);
};

#define TREE_GENERATOR_INLINE_H
#include "tree_generator_inline.h"
#undef INTERNAL_FILER_H

#endif
