#include <cassert>
#include <iostream>
#include <memory>
#include <fstream>

#include "getoptpp/getopt_pp.h"
#include "mili/mili.h"

#include "backbones-generator/petu.h"
#include "backbones-generator/generator.h"
#include "backbones-generator/factory_reader_chains.h"
#include "backbones-generator/tree_filters.h"
#include "backbones-generator/tree_data.h"
#include "backbones-generator/tree_helper.h"

int main(int argc, char** argv)
{
    CommandLineOptions o;

    std::cout << mili::getGPL3LicenseText("Backbones-generator", "2.0", "Facundo Muñoz, Daniel Gutson, Hugo Arregui, Marcos Villarreal and Rodrigo Castaño", "2006-2012");

    if (o.parse(argc, argv))
    {
        std::ifstream filer(o.data.c_str());

        TreeFilters treeFilters;
        Grillado grilla(o.m, o.n, o.z);

        prot_filer::AnglesMapping anglesMapping(o.Nres);
        prot_filer::AnglesData anglesData(o.Nres, &anglesMapping);
        
        IncrementalBackbone incrementalBackbone(o.Nres, grilla, anglesData, anglesMapping, treeFilters);

        TreeData treeData(o.Nres, incrementalBackbone);

        BareBackbone::treeData = &treeData;

        readData(filer, treeData, anglesMapping);//read angles

        // Fill r[][][] with the minimun squared distance between atoms
        treeFilters.setr(o.RN, o.RCa, o.RC, o.Scal_1_4, o.Scal_1_5);

        //o.outputFile
        TreeHelper treeHelper(treeData, treeFilters);

        std::cout << "Number of fi-si combinations in file=" << treeHelper.getNAngles() << std::endl;

        if (o.residuesInput.empty())
        {
            IGeneratorSimple* const generatorPtr = IGeneratorSimple::Factory::new_class(o.writeFormat);
            std::auto_ptr<IGeneratorSimple> g(generatorPtr);

            g->generate(treeHelper);
        }
        else
        {
            FullCachedAnglesSeqReader* const readerPtr = FactoryReaderChains::new_class(o.inputFormat, o.residuesInput, o.fragmentsFile);
            IGeneratorChains* const generatorPtr = IGeneratorChains::Factory::new_class(o.writeFormat);

            std::auto_ptr<FullCachedAnglesSeqReader> db(readerPtr);
            std::auto_ptr<IGeneratorChains> g(generatorPtr);

            g->generate(treeHelper, db.get());
        }

        std::cout << "Number of chains generated = " << treeData.cont << std::endl;

        prot_filer::Coord3DReaderFactory::destroy_instance();
        prot_filer::Coord3DSeqReaderFactory::destroy_instance();
        prot_filer::AnglesReaderFactory::destroy_instance();

        return EXIT_SUCCESS;
    }
    else
    {
        o.show_usage();

        return EXIT_FAILURE;
    }
}

void CommandLineOptions::show_usage()
{
    const std::string indent = "                  ";
    std::cerr << "Usage: " << std::endl;
    std::cerr << indent << "[ -i <data_file> ], default = data" << std::endl;
    std::cerr << indent << "[ -o <output file name> ], default = traj" << std::endl;
    std::cerr << indent << "[ {-w, --write_format} <xtc|compressed|fragments> ], default = xtc" << std::endl;
    std::cerr << indent << "[ --chains_input [fragments_file] <input_file> ], default = No using chains" << std::endl;
    std::cerr << indent << "if using chains: [ {-f, --input_format} <compressed|fragments> ], default = compressed " << std::endl;
}



//TODO:
// Hay que decidir si efectivamente estos valores queremos que se puedan configurar desde
// la linea de comandos.
bool CommandLineOptions::parse(int argc, char** argv)
{
    GetOpt::GetOpt_pp ops(argc, argv);
    std::vector<std::string> inputFiles;

    //loadDefaultOptions
    ops
            >> GetOpt::Option('n', "Rn",           RN,          1.5f)
            >> GetOpt::Option('a', "Rca",          RCa,         1.7f)
            >> GetOpt::Option('c', "Rc",           RC,          1.6f)
            >> GetOpt::Option('s', "Scal_1_4",     Scal_1_4,    0.85f)
            >> GetOpt::Option('l', "Scal_1_5",     Scal_1_5,    1.0f)
            >> GetOpt::Option('i', "input_file",   data,        std::string("data"))
            >> GetOpt::Option('o', "output_file",  outputFile,  std::string("traj"))
            >> GetOpt::Option('N', "rows",         n,           static_cast<size_t>(100))
            >> GetOpt::Option('M', "cols",         m,           static_cast<size_t>(100))
            >> GetOpt::Option('Z', "depth",        z,           static_cast<size_t>(100))
            >> GetOpt::Option('w', "write_format", writeFormat, std::string("xtc"));

    bool result;
    result = ops >> GetOpt::Option('r', "Nres", Nres);

    result = result && Nres != 0;
    if (Nres == 0)
        std::cerr << "Error: the amount of residues must be greater than zero" << std::endl;

    const bool chains = ops >> GetOpt::Option("chains_input", inputFiles);

    result = result && (chains || writeFormat != "fragments");
    if (!chains && writeFormat == "fragments")
        std::cerr << "Error: fragments output format cannot be used without chains input" << std::endl;

    if (result && chains)
    {
        ops >> GetOpt::Option('f', "input_format", inputFormat, "compressed");
        if (inputFormat == "compressed")
        {
            if (inputFiles.size() == 0)
            {
                std::cerr << "Error: Compressed input file required" << std::endl;
                result = false;
            }
            else
                residuesInput = inputFiles[0];
        }
        else if (inputFormat == "fragments")
        {
            if (inputFiles.size() < 2)
            {
                std::cerr << "Error: Fragments input format, requires both fragments files" << std::endl;
                result = false;
            }
            else
            {
                fragmentsFile = inputFiles[0];
                residuesInput = inputFiles[1];
            }
        }
    }//end chain

    return result;
}
