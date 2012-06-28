#include "readdata.h"

FullCachedAnglesSeqReader* read_chains(const std::string& format, const std::string& input_file, const std::string& fragments_file)
{
    FullCachedAnglesSeqReader* reader = NULL;

    if (format == "compressed")
    {
        prot_filer::AnglesReader* ar = prot_filer::AnglesReaderFactory::get_instance()->create(format);
        ar->open(input_file);
        reader = new FullCachedAnglesSeqReader(ar);
    }
    else if (format == "fragments")
    {
        prot_filer::AnglesReader* r = prot_filer::AnglesReaderFactory::get_instance()->create("compressed");
        prot_filer::FragmentsFromReader* fragments = new prot_filer::FragmentsFromReader(r);
        r->open(fragments_file);

        prot_filer::FragmentsAnglesReader* ar = prot_filer::FragmentsAnglesReaderFactory::get_instance()->create(format);
        ar->open(fragments, input_file);
        reader = new FullCachedAnglesSeqReader(ar);
    }
    else
    {
        throw std::runtime_error("wrong input format");
    }

    return reader;
}
