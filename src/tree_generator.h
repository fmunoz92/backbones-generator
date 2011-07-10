#include "petu.h"

template<template <class> class Generator, class Writer>
class WriterHelper;

template<class Writer>
class SimpleTreeGenerator
{
    typedef WriterHelper< ::SimpleTreeGenerator, Writer> Helper;
public:
    SimpleTreeGenerator(TreeData& tree_data, Helper& helper) :
        tree_data(tree_data),
        writer_helper(helper)
    {
        writer_helper.open();
    };
    ~SimpleTreeGenerator()
    {
        writer_helper.close();
    }
    void generate();
    TreeData& get_tree_data()
    {
        return tree_data;
    }
private:
    void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
    TreeData& tree_data;
    Helper& writer_helper;
};

template<class Writer>
class ChainsTreeGenerator
{
    typedef WriterHelper< ::ChainsTreeGenerator, Writer> Helper;
public:
    ChainsTreeGenerator(TreeData& tree_data, FullCachedAnglesSeqReader* reader_,
                        Helper& helper) :
        tree_data(tree_data),
        reader(reader_),
        writer_helper(helper)
    {
        writer_helper.open();
    };
    ~ChainsTreeGenerator()
    {
        writer_helper.close();
        reader->close();
        delete reader;
    }
    void generate();
    TreeData& get_tree_data()
    {
        return tree_data;
    }
    size_t get_fragment_nres() const
    {
        return reader->get_reader().get_atom_number() / 3;
    }
private:
    void generar_nivel_intermedio(unsigned int nivel, const float R_inicial[16], unsigned int indice_nivel_anterior);
    bool procesar_ultimo_nivel();
    TreeData& tree_data;
    FullCachedAnglesSeqReader* const reader;
    Helper& writer_helper;
};

template <template <class> class Generator, class Writer>
struct GenerateW;

template<class W>
struct GenerateW<SimpleTreeGenerator, W>
{
    void operator()(TreeData& tree_data, FullCachedAnglesSeqReader* /*reader*/)
    {
        WriterHelper<SimpleTreeGenerator, W> helper;
        SimpleTreeGenerator<W> g(tree_data, helper);
        g.generate();
    }
};

template<class W>
struct GenerateW<ChainsTreeGenerator, W>
{
    void operator()(TreeData& tree_data, FullCachedAnglesSeqReader* reader)
    {
        WriterHelper<ChainsTreeGenerator, W> helper;
        ChainsTreeGenerator<W> g(tree_data, reader, helper);
        g.generate();
    }
};

template <template <class> class Generator>
struct Generate
{
    void operator()(const string& format, TreeData& tree_data, FullCachedAnglesSeqReader* reader)
    {
        if (format == "xtc")
        {
            GenerateW<Generator, XtcWriter> g;
            g(tree_data, reader);
        }
        else if (format == "compressed")
        {
            GenerateW<Generator, CompressedWriter> g;
            g(tree_data, reader);
        }
        else
        {
            throw runtime_error("wrong format");
        }
    }
};

template <>
struct Generate<ChainsTreeGenerator>
{
    void operator()(const string& format, TreeData& tree_data, FullCachedAnglesSeqReader* reader)
    {
        if (format == "xtc")
        {
            GenerateW<ChainsTreeGenerator, XtcWriter> g;
            g(tree_data, reader);
        }
        else if (format == "compressed")
        {
            GenerateW<ChainsTreeGenerator, CompressedWriter> g;
            g(tree_data, reader);
        }
        else if (format == "fragments")
        {
            GenerateW<ChainsTreeGenerator, FragmentsWriter> g;
            g(tree_data, reader);
        }
        else
        {
            throw runtime_error("wrong format");
        }
    }
};

#define INTERNAL_FILER_H
#include "filer.h"
#undef INTERNAL_FILER_H
