#include "tree_generator.h"

void ChainsFormatGeneratorFragmentsWriter::generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader)
{
    TreeGenerator<ChainsTreeOperator<FragmentsWriterHelper> > generator(tree_data, reader);
    generator.generate();
}

void ChainsFormatGeneratorXtcWriter::generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader)
{
    TreeGenerator<ChainsTreeOperator<XtcWriterHelper> > generator(tree_data, reader);
    generator.generate();
}

void ChainsFormatGeneratorCompressedWriter::generate(TreeData& tree_data, FullCachedAnglesSeqReader* const reader)
{
    TreeGenerator<ChainsTreeOperator<CompressedWriterHelper> > generator(tree_data, reader);
    generator.generate();
}

void SimpleFormatGeneratorCompressedWriter::generate(TreeData& tree_data)
{
    TreeGenerator<SimpleTreeOperator<CompressedWriterHelper> > generator(tree_data, NULL);
    generator.generate();
}

void SimpleFormatGeneratorXtcWriter::generate(TreeData& tree_data)
{
    TreeGenerator<SimpleTreeOperator<XtcWriterHelper> > generator(tree_data, NULL);
    generator.generate();
}

REGISTER_FACTORIZABLE_CLASS(IGeneratorChains, ChainsFormatGeneratorXtcWriter,        std::string, "xtc");
REGISTER_FACTORIZABLE_CLASS(IGeneratorChains, ChainsFormatGeneratorFragmentsWriter,  std::string, "fragments");
REGISTER_FACTORIZABLE_CLASS(IGeneratorChains, ChainsFormatGeneratorCompressedWriter, std::string, "compressed");

REGISTER_FACTORIZABLE_CLASS(IGeneratorSimple, SimpleFormatGeneratorXtcWriter,        std::string, "xtc");
REGISTER_FACTORIZABLE_CLASS(IGeneratorSimple, SimpleFormatGeneratorCompressedWriter, std::string, "compressed");

inline XtcWriterHelper::XtcWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

inline XtcWriterHelper::~XtcWriterHelper()
{
    writer.close();
}

inline void XtcWriterHelper::write(TreeData& tree_data)
{
    writer.write(tree_data.atm, *tree_data.angles_data);
}

inline CompressedWriterHelper::CompressedWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

inline CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

inline void CompressedWriterHelper::write(TreeData& tree_data)
{
    writer.write(*tree_data.angles_data);
}

inline FragmentsWriterHelper::FragmentsWriterHelper(FullCachedAnglesSeqReader* reader) :
    output_file(output_file),
    reader(reader)
{
    writer.open(output_file);
}

inline FragmentsWriterHelper::~FragmentsWriterHelper()
{
    writer.close();
}

inline void FragmentsWriterHelper::write(TreeData& tree_data)
{
    size_t fragment_nres = reader->get_reader().get_atom_number() / 3;
    writer.write(fragment_nres, tree_data.fragment_ids, *tree_data.angles_data);
}
