#ifndef GENERATOR_H
#define GENERATOR_H

#include <mili/mili.h>

#include "tree_generator.h"
#include "tree_helper.h"
#include "filer.h"

struct IGeneratorSimple
{
    typedef mili::FactoryRegistry<IGeneratorSimple, std::string> Factory;

    virtual void generate(TreeHelper& tree_helper) = 0;
    virtual ~IGeneratorSimple() {}
};

struct IGeneratorChains
{
    typedef mili::FactoryRegistry<IGeneratorChains, std::string> Factory;

    virtual void generate(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader) = 0;
    virtual ~IGeneratorChains() {}
};

template <class Writer>
class GeneratorSimple : public IGeneratorSimple
{
public:
    virtual ~GeneratorSimple() {}
private:
    inline virtual void generate(TreeHelper& tree_helper);
};

template <class Writer>
class GeneratorChains : public IGeneratorChains
{
public:
    virtual ~GeneratorChains() {}
private:
    inline virtual void generate(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader);
};

//move to _inline
template <class Writer>
inline void GeneratorSimple<Writer>::generate(TreeHelper& tree_helper)
{
    TreeGenerator<SimpleTreeOperator<Writer> > generator(tree_helper, NULL);
    generator.generate();
}

template <class Writer>
inline void GeneratorChains<Writer>::generate(TreeHelper& tree_helper, FullCachedAnglesSeqReader* const reader)
{
    TreeGenerator<ChainsTreeOperator<Writer> > generator(tree_helper, reader);
    generator.generate();
}

#endif
