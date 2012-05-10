#ifndef INTERNAL_FILER_H
#error Internal header file, DO NOT include this.
#endif

const string output_f = "traj.xtc";

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

inline CompressedWriterHelper::~CompressedWriterHelper()
{
    writer.close();
}

inline CompressedWriterHelper::CompressedWriterHelper() :
    output_file(output_f)
{
    writer.open(output_file);
}

inline void CompressedWriterHelper::write(TreeData& tree_data)
{
    writer.write(*tree_data.angles_data);
}
