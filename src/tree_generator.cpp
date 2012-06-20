#include "tree_generator.h"

XtcWriterHelper::XtcWriterHelper(TreeData& tree_data) :
    tree_data(tree_data)
{
    writer.open(tree_data.output_file);
}

XtcWriterHelper::~XtcWriterHelper()
{
    writer.close();
}

void XtcWriterHelper::write()
{
    writer.write(tree_data.atm, *tree_data.angles_data);
}

CompressedWriterHelper::CompressedWriterHelper(TreeData& tree_data) :
    tree_data(tree_data)
{
    writer.open(tree_data.output_file);
}

CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

void CompressedWriterHelper::write()
{
    writer.write(*tree_data.angles_data);
}

FragmentsWriterHelper::FragmentsWriterHelper(TreeData& tree_data, FullCachedAnglesSeqReader* reader) :
    tree_data(tree_data),
    reader(reader)
{
    writer.open(tree_data.output_file);
}

FragmentsWriterHelper::~FragmentsWriterHelper()
{
    writer.close();
}

void FragmentsWriterHelper::write()
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
