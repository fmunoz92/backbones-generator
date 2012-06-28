#include <cassert>
#include <iostream>
#include <memory>
#include <fstream>

#include "getoptpp/getopt_pp.h"

#include "petu.h"

#include "generator.h"
#include "readdata.h"
#include "utils.h"
#include "tree_data.h"

int main(int argc, char** argv)
{
    CommandLineOptions o;

    if (o.parse(argc, argv))
    {
        std::ifstream filer(o.data.c_str());
        TreeData tree_data(o.Nres, o.m, o.n, o.z, filer);

        // Fill r[][][] with the minimun squared distance between atoms
        setr(o.RN, o.RCa, o.RC, o.Scal_1_4, o.Scal_1_5);

        std::cout << "Number of fi-si combinations in file=" << tree_data.cossi.size() << std::endl;

        if (o.residues_input.empty())
        {
            IGeneratorSimple* const generatorPtr = IGeneratorSimple::Factory::new_class(o.write_format);
            std::auto_ptr<IGeneratorSimple> g(generatorPtr);

            g->generate(tree_data);
        }
        else
        {
            FullCachedAnglesSeqReader* readerPtr = read_chains(o.input_format, o.residues_input, o.fragments_file);
            IGeneratorChains* const generatorPtr = IGeneratorChains::Factory::new_class(o.write_format);

            std::auto_ptr<FullCachedAnglesSeqReader> db(readerPtr);
            std::auto_ptr<IGeneratorChains> g(generatorPtr);

            g->generate(tree_data, db.get());
        }

        std::cout << "Number of chains generated=" << tree_data.cont << std::endl;

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
    std::string indent = "                  ";
    std::cerr << "Usage: " << std::endl;
    std::cerr << indent << "[ -i <data_file> ], default = data" << std::endl;
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
    std::vector<std::string> input_files;

    if (ops >> GetOpt::Option('r', "Nres", Nres))
    {
        ops
                >> GetOpt::Option('n', "Rn", RN, 1.5f)
                >> GetOpt::Option('a', "Rca", RCa, 1.7f)
                >> GetOpt::Option('c', "Rc", RC, 1.6f)
                >> GetOpt::Option('s', "Scal_1_4", Scal_1_4, 0.85f)
                >> GetOpt::Option('l', "Scal_1_5", Scal_1_5, 1.0f)
                >> GetOpt::Option('i', "input_file", data)
                >> GetOpt::Option('N', "rows", n, static_cast<size_t>(100))
                >> GetOpt::Option('M', "cols", m, static_cast<size_t>(100))
                >> GetOpt::Option('Z', "depth", z, static_cast<size_t>(100))
                >> GetOpt::Option('w', "write_format", write_format, std::string("xtc"));

        const bool chains = ops >> GetOpt::Option("chains_input", input_files);

        if (!chains && write_format == "fragments")
        {
            std::cerr << "Error: fragments output format cannot be used without chains input" << std::endl;
            return false;
        }
        if (chains)
        {
            ops >> GetOpt::Option('f', "input_format", input_format, "compressed");
            if (input_format == "compressed")
            {
                if (input_files.size() == 0)
                {
                    std::cerr << "Error: Compressed input file required" << std::endl;
                    return false;
                }
                else
                {
                    residues_input = input_files[0];
                }
            }
            else if (input_format == "fragments")
            {
                if (input_files.size() < 2)
                {
                    std::cerr << "Error: Fragments input format, requires both fragments files" << std::endl;
                    return false;
                }
                else
                {
                    fragments_file = input_files[0];
                    residues_input = input_files[1];
                }
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}
