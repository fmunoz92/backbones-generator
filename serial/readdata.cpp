
#include "readdata.h"

using namespace mili;

void readdata(std::ifstream& filer, std::vector<float> &cosfi, std::vector<float> &sinfi, std::vector<float> &cossi, std::vector<float> &sinsi, AnglesMapping* angles_mapping)
{

    float fi = 0.0f;
    float si = 0.0f;
    unsigned int i = 0;
    while (filer.good())
    {

        if (filer >> fi && filer >> si)
        {
            cosfi.push_back(cos(deg2rad(fi)));
            sinfi.push_back(sin(deg2rad(fi)));
            cossi.push_back(cos(deg2rad(si)));
            sinsi.push_back(sin(deg2rad(si)));
            angles_mapping->set_mapping(fi, si);
        }
        ++i;
    }

    // Nota a futuro: se deberia lanzar una Excepcion si el formato del archivo fuera equivocado.
}

AnglesDatabase* read_chains(const string& input_file, const string& read_format, const string& cache_type)
{
    FormatFiler* reader;
    reader = FilerFactory::get_instance()->create(read_format);
    reader->open_read(input_file);

    AnglesDatabase* residue_chain_database = new AnglesDatabase(reader, cache_type);

    int i = 0;

    while (!residue_chain_database->finished_reading())
    {
        const IncompleteAnglesData some_protein = (*residue_chain_database)[i];
        i += !residue_chain_database->finished_reading();
    }

    // this cannot be done here, since the reader is kept in the residue_chain_database,
    // and therefore it's his responsibility from now on. (and it's used later).
    //reader->close();
    //delete reader;

    return residue_chain_database;
}
