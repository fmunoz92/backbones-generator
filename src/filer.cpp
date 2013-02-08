#include <string>

#include "backbones-generator/filer.h"

//TODO: usar el nombre de la salida que quiere el usuario

XtcWriterHelper::XtcWriterHelper(TreeHelper& tree_helper)
    : tree_helper(tree_helper)
{
    const std::string XTC_EXTENSION = ".xtc";
    writer.open("traj" + XTC_EXTENSION);
}

XtcWriterHelper::~XtcWriterHelper()
{
    writer.close();
}

void XtcWriterHelper::write()
{
    writer.write(tree_helper.getAtm(), tree_helper.getAtm().getAnglesData());
}

CompressedWriterHelper::CompressedWriterHelper(TreeHelper& tree_helper)
    : tree_helper(tree_helper)
{
    const std::string COMPRESSED_EXTENSION = ".cps";
    writer.open("traj" + COMPRESSED_EXTENSION);
}

CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

void CompressedWriterHelper::write()
{
    writer.write(tree_helper.getAtm().getAnglesData());
}

FragmentsWriterHelper::FragmentsWriterHelper(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader)
    : tree_helper(tree_helper),
      reader(reader)
{
    const std::string FRAGMENSTS_EXTENSION = ".fgs";
    writer.open("traj" + FRAGMENSTS_EXTENSION);
}

FragmentsWriterHelper::~FragmentsWriterHelper()
{
    writer.close();
}

void FragmentsWriterHelper::write()
{
    const size_t fragment_nres = reader->get_reader().get_atom_number() / 3;

    writer.write(fragment_nres, tree_helper.getAtm().getFragmentIds(), tree_helper.getAtm().getAnglesData());
}
