#include <cassert>
#include <iostream>
#include <fstream>
#include "getopt_pp.h"

#include <mili/mili.h>
#include "tree_generator.h"
#include "readdata.h"
#include "utils.h"

static void show_usage();

int main(int argc, char** argv)
{
    using namespace GetOpt;

    int     Nres;   // Number of amino acids in the chains to build.

    float   RN,         // Radius of the Nitrogen atom.
            RCa,        // Radius of the Carbon atom.
            RC,         // Radius of the Carbon atom.
            Scal_1_4,   // Scaling factor for the radii. Used to check for 1-4 clashes.
            Scal_1_5;   // Scaling factor for the radii. Used to check for 1-5 clashes.

    std::string data;   // Name of the input file.
    std::string write_format;
    std::string residues_input; //Name of the residue chains file.

    size_t m; // They indicate the size of each dimention of the grid.
    size_t n;
    size_t z;

    GetOpt_pp ops(argc, argv);

    if (ops >> Option('r', "Nres", Nres))
    {
        ops
                >> Option('n', "Rn", RN, 1.5f)
                >> Option('a', "Rca", RCa, 1.7f)
                >> Option('c', "Rc", RC, 1.6f)
                >> Option('s', "Scal_1_4", Scal_1_4, 0.85f)
                >> Option('l', "Scal_1_5", Scal_1_5, 1.0f)
                >> Option('i', "input_file", data, "ramachandran.dat")
                >> Option('N', "rows", n, static_cast<size_t>(100))
                >> Option('M', "cols", m, static_cast<size_t>(100))
                >> Option('Z', "depth", z, static_cast<size_t>(100))
                >> Option('w', "write_format", write_format, string("xtc"))
                >> Option("chains_input", residues_input);
        //TODO:
        // Hay que decidir si efectivamente estos valores queremos que se puedan configurar desde
        // la linea de comandos.
        std::ifstream filer;
        filer.open(data.c_str());
        Grillado* grilla = new Grillado(m, n, z);
        TreeData tree_data(Nres, grilla);

        // Fill r[][][] with the minimun squared distance between atoms
        setr(RN, RCa, RC, Scal_1_4, Scal_1_5);

        readdata(filer, tree_data);
        tree_data.angles_data = new AnglesData(tree_data.nres, *tree_data.angles_mapping);

        cout << "Number of fi-si combinations in file=" << tree_data.cossi.size() << endl;

        TreeGenerator* generator = NULL;

        if (residues_input.empty())
        {
            generator = new SimpleTreeGenerator(tree_data, write_format);
        }
        else
        {
            FullCachedAnglesSeqReader* residue_chain_database = read_chains(residues_input);
            generator = new ChainsTreeGenerator(tree_data, write_format, residue_chain_database);
        }
        generator->generate();

        cout << "Number of chains generated=" << tree_data.cont << endl;

        Coord3DReaderFactory::destroy_instance();
        Coord3DSeqReaderFactory::destroy_instance();
        AnglesReaderFactory::destroy_instance();

        //residue_chain_database is deleted by TreeGenerator
        filer.close();
        delete generator;
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
    std::cerr << "Invalid arguments." << std::endl;
}
