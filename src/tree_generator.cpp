#include "tree_generator.h"

XtcWriterHelper::XtcWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

XtcWriterHelper::~XtcWriterHelper()
{
    writer.close();
}

void XtcWriterHelper::write(TreeData& tree_data)
{
    writer.write(tree_data.atm, *tree_data.angles_data);
}

CompressedWriterHelper::CompressedWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

void CompressedWriterHelper::write(TreeData& tree_data)
{
    writer.write(*tree_data.angles_data);
}

FragmentsWriterHelper::FragmentsWriterHelper(FullCachedAnglesSeqReader* reader) :
    output_file(output_file),
    reader(reader)
{
    writer.open(output_file);
}

FragmentsWriterHelper::~FragmentsWriterHelper()
{
    writer.close();
}

void FragmentsWriterHelper::write(TreeData& tree_data)
{
    size_t fragment_nres = reader->get_reader().get_atom_number() / 3;
    writer.write(fragment_nres, tree_data.fragment_ids, *tree_data.angles_data);
}

/* REGISTER_FACTORIZABLE_CLASS TAKE A TYPEDEF IF DERIVED CLASS IS A TEMPLATE-CLASS*/
typedef GeneratorChains<CompressedWriterHelper> ChainsFormatGeneratorCompressedWriter;
typedef GeneratorChains<XtcWriterHelper>        ChainsFormatGeneratorXtcWriter;
typedef GeneratorChains<FragmentsWriterHelper>  ChainsFormatGeneratorFragmentsWriter;

typedef GeneratorSimple<XtcWriterHelper>        SimpleFormatGeneratorXtcWriter;
typedef GeneratorSimple<CompressedWriterHelper> SimpleFormatGeneratorCompressedWriter;


REGISTER_FACTORIZABLE_CLASS(IGeneratorChains, ChainsFormatGeneratorXtcWriter,        std::string, "xtc");
REGISTER_FACTORIZABLE_CLASS(IGeneratorChains, ChainsFormatGeneratorFragmentsWriter,  std::string, "fragments");
REGISTER_FACTORIZABLE_CLASS(IGeneratorChains, ChainsFormatGeneratorCompressedWriter, std::string, "compressed");

REGISTER_FACTORIZABLE_CLASS(IGeneratorSimple, SimpleFormatGeneratorXtcWriter,        std::string, "xtc");
REGISTER_FACTORIZABLE_CLASS(IGeneratorSimple, SimpleFormatGeneratorCompressedWriter, std::string, "compressed");
