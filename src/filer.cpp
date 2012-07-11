#include <string>

#include "filer.h"

XtcWriterHelper::XtcWriterHelper(TreeHelper& tree_helper)
    : tree_helper(tree_helper)
{
    const std::string XTC_EXTENSION = ".xtc";
    writer.open(tree_helper.getOutputFile() + XTC_EXTENSION);
}

XtcWriterHelper::~XtcWriterHelper()
{
    writer.close();
}

void XtcWriterHelper::write()
{
    writer.write(tree_helper.getAtm(), tree_helper.getAnglesData());
}

CompressedWriterHelper::CompressedWriterHelper(TreeHelper& tree_helper)
    : tree_helper(tree_helper)
{
    const std::string COMPRESSED_EXTENSION = ".cps";
    writer.open(tree_helper.getOutputFile() + COMPRESSED_EXTENSION);
}

CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

void CompressedWriterHelper::write()
{
    writer.write(tree_helper.getAnglesData());
}

FragmentsWriterHelper::FragmentsWriterHelper(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader)
    : tree_helper(tree_helper),
      reader(reader)
{
    const std::string FRAGMENSTS_EXTENSION = ".fgs";
    writer.open(tree_helper.getOutputFile() + FRAGMENSTS_EXTENSION);
}

FragmentsWriterHelper::~FragmentsWriterHelper()
{
    writer.close();
}

void FragmentsWriterHelper::write()
{
    size_t fragment_nres = reader->get_reader().get_atom_number() / 3;

    writer.write(fragment_nres, tree_helper.getFragmentIds(), tree_helper.getAnglesData());
}
