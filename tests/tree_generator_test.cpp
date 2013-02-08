#include <math.h>
#include <memory>
#include <fstream>

#include <gtest/gtest.h>
#include <mili/mili.h>

#include "backbones-generator/generator.h"
#include "backbones-generator/factory_reader_chains.h"
#include "backbones-generator/tree_filters.h"
#include "backbones-generator/tree_data.h"
#include "backbones-generator/tree_helper.h"

/*
 *
 * debug en modo cadenas:
 * cadenas = ceil(nres / chain_size)
 * soluciones = fragmentos^cadenas * angulos^(cadenas - 1)
 *
 * debug en modo simple:
 * equivalencia con modo cadenas: fragmentos = 1, chain_size = 1, cadenas = nres
 * soluciones = angulos^(cadenas - 2)// menos 2 xq la semilla aca cuenta
 *
 */

struct TestHelperSingleMode
{
    static unsigned int count(unsigned int angles, unsigned int residues)
    {
        return (residues < 3) ? 1 : std::pow(angles, (residues - 2));
    }

    static unsigned int testCount(unsigned int nres, std::istream& inputFile)
    {
        const float radius = 0.5;
        const float scal = 1;
        //const std::string outputFile("testSimpleMode");

        TreeFilters treeFilters;
        Grillado grilla(100, 100, 100);
        prot_filer::AnglesMapping anglesMapping(nres);
        prot_filer::AnglesData anglesData(nres, &anglesMapping);
        IncrementalBackbone incrementalBackbone(nres, grilla, anglesData, anglesMapping, treeFilters);
        TreeData treeData(nres);
        BareBackbone::treeData = &treeData;
        readData(inputFile, treeData, anglesMapping);
        treeFilters.setr(radius, radius, radius, scal, scal);
        TreeHelper treeHelper(treeData, incrementalBackbone);

        IGeneratorSimple* const generatorPtr = IGeneratorSimple::Factory::new_class("compressed");
        std::auto_ptr<IGeneratorSimple> g(generatorPtr);

        const  unsigned int count = g->generate(treeHelper);

        return count;
    }

};

struct TestHelperChainMode
{
    static unsigned int testCount(unsigned int nres, std::istream& inputFile, const std::string& residuesInput)
    {
        const float radius = 0.5;
        const float scal = 1;
        //const std::string outputFile("testChainMode");

        TreeFilters treeFilters;
        Grillado grilla(100, 100, 100);
        prot_filer::AnglesMapping anglesMapping(nres);
        prot_filer::AnglesData anglesData(nres, &anglesMapping);
        IncrementalBackbone incrementalBackbone(nres, grilla, anglesData, anglesMapping, treeFilters);
        TreeData treeData(nres);
        BareBackbone::treeData = &treeData;
        readData(inputFile, treeData, anglesMapping);
        treeFilters.setr(radius, radius, radius, scal, scal);
        TreeHelper treeHelper(treeData, incrementalBackbone);

        FullCachedAnglesSeqReader* const readerPtr = FactoryReaderChains::new_class("compressed", residuesInput, "");
        IGeneratorChains* const generatorPtr = IGeneratorChains::Factory::new_class("xtc");

        std::auto_ptr<FullCachedAnglesSeqReader> db(readerPtr);
        std::auto_ptr<IGeneratorChains> g(generatorPtr);

        const unsigned int count = g->generate(treeHelper, db.get());


        return count;
    }

    static void generateSimple(const std::string& /*outputFile*/, std::istream& inputFile, unsigned int nres)
    {
        const float radius = 0.5;
        const float scal = 1;

        TreeFilters treeFilters;
        Grillado grilla(100, 100, 100);
        prot_filer::AnglesMapping anglesMapping(nres);
        prot_filer::AnglesData anglesData(nres, &anglesMapping);
        IncrementalBackbone incrementalBackbone(nres, grilla, anglesData, anglesMapping, treeFilters);
        TreeData treeData(nres);
        BareBackbone::treeData = &treeData;
        readData(inputFile, treeData, anglesMapping);
        treeFilters.setr(radius, radius, radius, scal, scal);
        TreeHelper treeHelper(treeData, incrementalBackbone);

        IGeneratorSimple* const generatorPtr = IGeneratorSimple::Factory::new_class("compressed");
        std::auto_ptr<IGeneratorSimple> g(generatorPtr);
        g->generate(treeHelper);
    }

    static unsigned int count(unsigned int fragments,  unsigned int chainSize, unsigned int angles, unsigned int residues)
    {
        const unsigned int amountChains = (residues / chainSize);//ceil?
        return std::pow(fragments, amountChains) * std::pow(angles, (amountChains - 1));
    }
};

TEST(TestTreeGenerator, CountSingleMode_NRES_2)
{
    const unsigned int ANGLES = 2;
    const unsigned int NRES = 2;

    std::stringstream inputFile("0 0\n0 90\n");

    const unsigned int amountGenerated = TestHelperSingleMode::testCount(NRES, inputFile);
    ASSERT_EQ(amountGenerated, TestHelperSingleMode::count(ANGLES, NRES));
}

TEST(TestTreeGenerator, CountSingleMode_NRES_4)
{
    const unsigned int ANGLES = 2;
    const unsigned int NRES = 4;

    std::stringstream inputFile("0 0\n0 90\n");

    const unsigned int amountGenerated = TestHelperSingleMode::testCount(NRES, inputFile);
    ASSERT_EQ(amountGenerated, TestHelperSingleMode::count(ANGLES, NRES));
}

TEST(TestTreeGenerator, CountSingleMode_NRES_10)
{
    const unsigned int ANGLES = 2;
    const unsigned int NRES = 10;

    std::stringstream inputFile("0 0\n0 90\n");

    const unsigned int amountGenerated = TestHelperSingleMode::testCount(NRES, inputFile);
    ASSERT_EQ(amountGenerated, TestHelperSingleMode::count(ANGLES, NRES));
}

TEST(TestTreeGenerator, CountChainMode_NRES_9)
{
    const unsigned int NRES = 9;
    const unsigned int ANGLES = 2;

    const unsigned int NRES_SINGLE = 4;
    const unsigned int ANGLES_SINGLE = 2;

    const unsigned int CHAIN_SIZE = NRES_SINGLE - 1;
    const unsigned int FRAGMENTS = TestHelperSingleMode::count(ANGLES_SINGLE, NRES_SINGLE);

    const std::string outputFileSingleMode("traj");

    std::stringstream inputFileSimple("0 0\n0 90\n");

    TestHelperChainMode::generateSimple(outputFileSingleMode, inputFileSimple, NRES_SINGLE);

    const std::string residuesInput = outputFileSingleMode + ".cps";

    std::stringstream inputFile("0 0\n0 90\n");
    const unsigned int amountGenerated = TestHelperChainMode::testCount(NRES, inputFile, residuesInput);

    ASSERT_EQ(amountGenerated, TestHelperChainMode::count(FRAGMENTS, CHAIN_SIZE, ANGLES, NRES));
}

TEST(TestTreeGenerator, CountChainMode_NRES_9_WHIT_DIFERENTS_ANGLES)
{
    const unsigned int ANGLES = 3;
    const unsigned int NRES = 9;

    const unsigned int ANGLES_SINGLE = 2;
    const unsigned int NRES_SINGLE = 4;

    const unsigned int CHAIN_SIZE = NRES_SINGLE - 1;
    const unsigned int FRAGMENTS = TestHelperSingleMode::count(ANGLES_SINGLE, NRES_SINGLE);

    const std::string outputFileSingleMode("traj");
    std::stringstream inputFileSimple("0 0\n0 90\n");

    TestHelperChainMode::generateSimple(outputFileSingleMode, inputFileSimple, NRES_SINGLE);

    const std::string residuesInput = outputFileSingleMode + ".cps";

    std::stringstream inputFile("0 0\n0 90\n90 0\n");
    const unsigned int amountGenerated = TestHelperChainMode::testCount(NRES, inputFile, residuesInput);

    ASSERT_EQ(amountGenerated, TestHelperChainMode::count(FRAGMENTS, CHAIN_SIZE, ANGLES, NRES));
}
