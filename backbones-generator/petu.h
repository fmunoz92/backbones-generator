#ifndef PETU_H
#define PETU_H

#include <string>

struct CommandLineOptions
{
    unsigned int Nres;// Number of amino acids in the chains to build.
    float RN;         // Radius of the Nitrogen atom.
    float RCa;        // Radius of the Carbon atom.
    float RC;         // Radius of the Carbon atom.
    float Scal_1_4;   // Scaling factor for the radii. Used to check for 1-4 clashes.
    float Scal_1_5;   // Scaling factor for the radii. Used to check for 1-5 clashes.
    std::string data;      // Name of the input file.
    std::string outputFile; // Name of the output file.
    std::string writeFormat;
    std::string residuesInput; //Name of the residue chains file.
    std::string inputFormat;
    std::string fragmentsFile;
    size_t m; // They indicate the size of each dimention of the grid.
    size_t n;
    size_t z;

    bool parse(int arg, char** argv);
    void show_usage();
};

#endif
