#ifndef GENERATOR_H
#define GENERATOR_H

#include <mili/mili.h>

#include "tree_generator.h"
#include "tree_data.h"
#include "filer.h"

struct IGeneratorSimple
{
    typedef mili::FactoryRegistry<IGeneratorSimple, std::string> Factory;

    virtual void generate(TreeData& tree_data) = 0;
    virtual ~IGeneratorSimple() {}
};

struct IGeneratorChains
{
    typedef mili::FactoryRegistry<IGeneratorChains, std::string> Factory;

    virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader) = 0;
    virtual ~IGeneratorChains() {}
};

template <class Writer>
class GeneratorSimple : public IGeneratorSimple
{
public:
    virtual ~GeneratorSimple() {}
private:
    inline virtual void generate(TreeData& tree_data);
};

template <class Writer>
class GeneratorChains : public IGeneratorChains
{
public:
    virtual ~GeneratorChains() {}
private:
    inline virtual void generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader);
};

//move to _inline
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

#endif
