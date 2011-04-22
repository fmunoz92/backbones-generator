#include <cassert>
#include <iostream>
#include <fstream>
#include "getopt_pp.h"

#include <mili/mili.h>
#include "tree_generator.h"
#include "readdata.h"
#include "utils.h"

static void show_usage();

using namespace GetOpt;
using std::string;

struct CommandLineOptions
{
    int Nres;   // Number of amino acids in the chains to build.
    float RN;         // Radius of the Nitrogen atom.
    float RCa;        // Radius of the Carbon atom.
    float RC;         // Radius of the Carbon atom.
    float Scal_1_4;   // Scaling factor for the radii. Used to check for 1-4 clashes.
    float Scal_1_5;   // Scaling factor for the radii. Used to check for 1-5 clashes.
    string data;   // Name of the input file.
    string write_format;
    string residues_input; //Name of the residue chains file.
    string input_format;
    string fragments_file;
    size_t m; // They indicate the size of each dimention of the grid.
    size_t n;
    size_t z;

    bool parse(int arg, char** argv);
};

int main(int argc, char** argv)
{
    CommandLineOptions o;

    if (o.parse(argc, argv))
    {
        std::ifstream filer;
        filer.open(o.data.c_str());
        Grillado* grilla = new Grillado(o.m, o.n, o.z);
        TreeData tree_data(o.Nres, grilla);

        // Fill r[][][] with the minimun squared distance between atoms
        setr(o.RN, o.RCa, o.RC, o.Scal_1_4, o.Scal_1_5);

        readdata(filer, tree_data);
        tree_data.angles_data = new AnglesData(tree_data.nres, *tree_data.angles_mapping);

        cout << "Number of fi-si combinations in file=" << tree_data.cossi.size() << endl;

        if (o.residues_input.empty())
        {
            Generate<SimpleTreeGenerator> g;
            g(o.write_format, tree_data, NULL);
        }
        else
        {
            FullCachedAnglesSeqReader* db = read_chains(o.input_format, o.residues_input, o.fragments_file);
            Generate<ChainsTreeGenerator> g;
            g(o.write_format, tree_data, db);
        }

        cout << "Number of chains generated=" << tree_data.cont << endl;

        Coord3DReaderFactory::destroy_instance();
        Coord3DSeqReaderFactory::destroy_instance();
        AnglesReaderFactory::destroy_instance();

        filer.close();
        return EXIT_SUCCESS;
    }
    else
    {
        show_usage();
        return EXIT_FAILURE;
    }
}

void show_usage()
{
    string indent = "                  ";
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
    GetOpt_pp ops(argc, argv);
    std::vector<string> input_files;

    if (ops >> Option('r', "Nres", Nres))
    {
        ops
            >> Option('n', "Rn", RN, 1.5f)
            >> Option('a', "Rca", RCa, 1.7f)
            >> Option('c', "Rc", RC, 1.6f)
            >> Option('s', "Scal_1_4", Scal_1_4, 0.85f)
            >> Option('l', "Scal_1_5", Scal_1_5, 1.0f)
            >> Option('i', "input_file", data)
            >> Option('N', "rows", n, static_cast<size_t>(100))
            >> Option('M', "cols", m, static_cast<size_t>(100))
            >> Option('Z', "depth", z, static_cast<size_t>(100))
            >> Option('w', "write_format", write_format, string("xtc"));

        bool chains = ops >> Option("chains_input", input_files);

        if (!chains && write_format == "fragments")
        {
            cerr << "Error: fragments output format cannot be used without chains input" << endl;
            return false;
        }
        if (chains)
        {
            ops >> Option('f', "input_format", input_format, "compressed");
            if (input_format == "compressed")
            {
                if (input_files.size() == 0)
                {
                    cerr << "Error: Compressed input file required" << endl;
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
                    cerr << "Error: Fragments input format, requires both fragments files" << endl;
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