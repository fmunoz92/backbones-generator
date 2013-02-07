#ifndef TREE_DATA_H
#define TREE_DATA_H

#include <string>

#include "prot-filer/format_filer.h"
#include "prot-filer/cached_reader.h"

//TODO: move to "definitions.h"
typedef prot_filer::BasicProtein Atoms;
typedef prot_filer::CachedReader<prot_filer::FullCache, prot_filer::SimpleAnglesReader, prot_filer::AnglesData> FullCachedAnglesSeqReader;

#include "backbones-generator/tree_filters.h"
#include "backbones-generator/grillado.h"

struct TreeData;
struct TreeFilters;

struct Residuo
{
    Residuo(const Residuo& r)
        : at1(r.at1),
          at2(r.at2),
          at3(r.at3)
    {}

    Residuo(const esferaId& param_at1, const esferaId& param_at2, const esferaId& param_at3)
        : at1(param_at1),
          at2(param_at2),
          at3(param_at3)
    {};

    Residuo() {};

    esferaId at1;
    esferaId at2;
    esferaId at3;
};


struct BareBackbone : public Atoms
{
    BareBackbone(unsigned int nRes, Grillado& grilla, prot_filer::AnglesData& anglesData, prot_filer::AnglesMapping& anglesMapping, TreeFilters& treeFilters);

    void putSeed(float* R, Residuo& residuo);

    bool putRes(float* pR, const unsigned int resN, Residuo& residuo, unsigned int siIndex, unsigned int fiIndex);

    void deleteRes(const Residuo& residuo);

    void deleteRes(const std::list<Residuo>& residuos);

    void clear(); //This function puts 0 in all the atoms coordinates

    void deleteLastFragmentId();

    void pushChainIndex(unsigned int index);

    prot_filer::AnglesData& getAnglesData();
    const prot_filer::FragmentIds& getFragmentIds() const;

    static TreeData* treeData;
protected:
    TreeFilters& treeFilters;

    Grillado& grilla; // Utilizamos el grillado para aproximar el volumen parcial
    prot_filer::AnglesMapping& anglesMapping;
    prot_filer::AnglesData&    anglesData; // Used only when writing compressed data.
    prot_filer::FragmentIds fragmentIds;
};

struct IncrementalBackbone : public BareBackbone
{
    IncrementalBackbone(unsigned int nRes, Grillado& grilla, prot_filer::AnglesData& anglesData, prot_filer::AnglesMapping& anglesMapping, TreeFilters& treeFilters);

    bool filterLastLevelOk();
};

struct TreeData
{
    TreeData(unsigned int nRes);

    const unsigned int nRes;
    unsigned int nAngles; // assigned in readData(..)

    const float rgmax;
    const float dmax2;

    std::vector<float> cosfi;
    std::vector<float> cossi;
    std::vector<float> sinfi;
    std::vector<float> sinsi;
};


void readData(std::istream& filer, TreeData& treeData, prot_filer::AnglesMapping& anglesMapping);

#endif
