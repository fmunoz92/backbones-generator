#ifndef GENERATOR_H
#define GENERATOR_H

#include "tree_data.h"

struct IGeneratorSimple
{
    virtual ~IGeneratorSimple() {}
    virtual void generate(TreeData& tree_data) = 0;
};

template <class Writer>
class GeneratorSimple : public IGeneratorSimple
{
public:
    virtual ~GeneratorSimple() {}
private:
    inline virtual void generate(TreeData& tree_data);
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

typedef mili::FactoryRegistry<IGeneratorChains, std::string> FactoryGeneratorChains;
typedef mili::FactoryRegistry<IGeneratorSimple, std::string> FactoryGeneratorSimple;

#endif
