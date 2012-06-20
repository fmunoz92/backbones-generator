#include "filer.h"

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
