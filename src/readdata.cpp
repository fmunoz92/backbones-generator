#include "readdata.h"

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
