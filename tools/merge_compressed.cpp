
#include <iostream>
#include "prot-filer/angles.h"
#include "getopt_pp_standalone.h"
#include "prot-filer/read_utils.h"
#include <fstream>
#include <vector>


#define DIM 3
typedef float matrix[DIM][DIM];
using namespace std;
using namespace GetOpt;
void print_help();
void copy_structures(istream& input, ostream& output, const AnglesMapping& angles_mapping);
void read_box(matrix& box, istream& file);
void write_box(const matrix& box, ostream& file);
bool boxes_are_equal(const matrix& box1, const matrix& box2);
int main(int argc , char** argv)
{

    string output_filename;
    vector<string> filenames;
    vector<string> rejected;
    GetOpt_pp ops(argc, argv);
    ops >> Option('i', "input", filenames);
    ops >> Option('o', "output", output_filename, "merged.comp");
    AnglesMapping header(1); // this nres will be overwritten with read_mapping();
    fstream output_file;
    output_file.open(output_filename.c_str(), ios::binary | ios::out);
    if (filenames.size() > 0)
    {

        fstream input_file;
        matrix original_box;
        matrix box;
        input_file.open(filenames[0].c_str(), ios::binary | ios::in);
        header.read_mapping(input_file);
        read_box(original_box, input_file);
        header.write_mapping(output_file);
        write_box(original_box, output_file);
        copy_structures(input_file, output_file, header);
        input_file.close();
        AnglesMapping angles_mapping(1); // this nres will be overwritten with read_mapping();
        for (unsigned int i(1); i < filenames.size(); ++i)
        {
            input_file.open(filenames[i].c_str(), ios::binary | ios::in);
            angles_mapping.read_mapping(input_file);
            read_box(box, input_file);
            if (angles_mapping == header && boxes_are_equal(box, original_box))
            {
                copy_structures(input_file, output_file, angles_mapping);
            }
            else
            {
                rejected.push_back(filenames[i]);
            }
            input_file.close();
        }
        output_file.close();
        if (rejected.size() > 0)
        {
            cerr << "The following files were rejected because their header didn't match the header in " << filenames[0] << "." << endl;
            cerr << "Rejected files: " << endl;
            for (vector<string>::iterator it = rejected.begin(); it != rejected.end(); ++it)
            {
                cerr << *it << endl;
            }
        }
    }
    else
    {
        print_help();
    }

    return 0;
}


void print_help()
{
    cout << "Instrucciones " << endl;
}


void copy_structures(istream& input, ostream& output, const AnglesMapping& angles_mapping)
{
    AnglesData structure(angles_mapping.get_nres(), angles_mapping);
    while (!structure.read_structure(input))
    {
        structure.write_structure(output);
    }
}


bool boxes_are_equal(const matrix& box1, const matrix& box2)
{
    bool result = true;
    for (size_t i(0); result && i < DIM; ++i)
    {
        for (size_t j(0); result && j < DIM; ++j)
        {
            result = result && (box1[i][j] == box2[i][j]);
        }
    }
    return result;
}
void read_box(matrix& _box, istream& file)
{
    //DIM is defined in xdrfile.h (should be 3)
    for (size_t i(0); i < DIM; ++i)
    {
        for (size_t j(0); j < DIM; ++j)
        {
            float temp;
            read_var(file, temp);
            _box[i][j] = temp;
        }
    }
}
void write_box(const matrix& box, ostream& file)
{
    //DIM is defined in xdrfile.h (should be 3)
    for (size_t i(0); i < DIM; ++i)
    {
        for (size_t j(0); j < DIM; ++j)
        {
            write_var(file, box[i][j]);
        }
    }
}
