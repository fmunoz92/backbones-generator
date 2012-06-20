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

FullCachedAnglesSeqReader* read_chains(const string& format, const string& input_file, const string& fragments_file)
{
    if (format == "compressed")
    {
        prot_filer::AnglesReader* r = prot_filer::AnglesReaderFactory::get_instance()->create(format);
        r->open(input_file);
        return new FullCachedAnglesSeqReader(r);
    }
    else if (format == "fragments")
    {
        prot_filer::AnglesReader* r = prot_filer::AnglesReaderFactory::get_instance()->create("compressed");
        prot_filer::FragmentsFromReader* fragments = new prot_filer::FragmentsFromReader(r);
        r->open(fragments_file);

        prot_filer::FragmentsAnglesReader* ar = prot_filer::FragmentsAnglesReaderFactory::get_instance()->create(format);
        ar->open(fragments, input_file);
        FullCachedAnglesSeqReader* reader = new FullCachedAnglesSeqReader(ar);
        return reader;
    }
    else
    {
        throw runtime_error("wrong input format");
    }
}
