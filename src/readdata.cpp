#include "readdata.h"

using namespace mili;

void readdata(std::istream& filer, TreeData& tree_data)
{
    float fi = 0.0f;
    float si = 0.0f;
    unsigned int i = 0;
    while (filer.good())
    {
        if (filer >> fi && filer >> si)
        {
            tree_data.cosfi.push_back(cos(deg2rad(fi)));
            tree_data.sinfi.push_back(sin(deg2rad(fi)));
            tree_data.cossi.push_back(cos(deg2rad(si)));
            tree_data.sinsi.push_back(sin(deg2rad(si)));
            tree_data.angles_mapping->set_mapping(fi, si);
        }
        ++i;
    }
    // Nota a futuro: se deberia lanzar una Excepcion si el formato del archivo fuera equivocado.
}

FullCachedAnglesSeqReader* read_chains(const string& input_file)
{
    AnglesReader* r = AnglesReaderFactory::get_instance()->create("compressed");
    FullCachedAnglesSeqReader* reader = new FullCachedAnglesSeqReader(r);
    reader->open(input_file);
    return reader;
}
