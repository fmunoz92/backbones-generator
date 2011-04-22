#include "petu.h"
#include "prot-filer/compressed_filer.h"
#include "prot-filer/xtc_filer.h"

static const string file = "traj.xtc";

template<template <class> class Generator, class Writer>
class WriterHelper
{};

template<template <class> class Generator>
class WriterHelper<Generator, XtcWriter>
{
public:
    virtual void open()
    {
        writer.open(file);
    }
    virtual void write(Generator<XtcWriter>& g)
    {
        TreeData& tree_data = g.get_tree_data();
        writer.write(tree_data.atm, *tree_data.angles_data);
    }
    virtual void close()
    {
        writer.close();
    }
private:
    XtcWriter writer;
};

template<template <class> class Generator>
class WriterHelper<Generator, CompressedWriter>
{
public:
    virtual void open()
    {
        writer.open(file);
    }
    virtual void write(Generator<CompressedWriter>& g)
    {
        TreeData& tree_data = g.get_tree_data();
        writer.write(*tree_data.angles_data);
    }
    virtual void close()
    {
        writer.close();
    }
private:
    CompressedWriter writer;
};

template<>
class WriterHelper<ChainsTreeGenerator, FragmentsWriter>
{
public:
    virtual void open()
    {
        writer.open(file);
    }
    virtual void write(ChainsTreeGenerator<FragmentsWriter>& g)
    {
        TreeData& tree_data = g.get_tree_data();
        writer.write(g.get_fragment_nres(), tree_data.fragment_ids, *tree_data.angles_data);
    }
    virtual void close()
    {
        writer.close();
    }
private:
    FragmentsWriter writer;
};
