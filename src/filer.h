#ifndef FILER_H
#define FILER_H

#include "tree_data.h"

class XtcWriterHelper
{
public:
    XtcWriterHelper(TreeData& tree_data);
    ~XtcWriterHelper();
    void write();
private:
    TreeData& tree_data;
    prot_filer::XtcWriter writer;
};

class CompressedWriterHelper
{
public:
    CompressedWriterHelper(TreeData& tree_data);
    ~CompressedWriterHelper();
    void write();
private:
    TreeData& tree_data;
    prot_filer::CompressedWriter writer;
};

class FragmentsWriterHelper
{
public:
    FragmentsWriterHelper(TreeData& tree_data, FullCachedAnglesSeqReader* reader);//Adapter
    ~FragmentsWriterHelper();
    void write();
private:
    TreeData& tree_data;
    const FullCachedAnglesSeqReader* reader;
    prot_filer::FragmentsWriter writer;
};
#endif
