#ifndef FILER_H
#define FILER_H

#include "poneres.h"//TreeHelper

class XtcWriterHelper
{
public:
    XtcWriterHelper(TreeHelper& tree_helper);
    ~XtcWriterHelper();
    void write();
private:
    TreeHelper& tree_helper;
    prot_filer::XtcWriter writer;
};

class CompressedWriterHelper
{
public:
    CompressedWriterHelper(TreeHelper& tree_helper);
    ~CompressedWriterHelper();
    void write();
private:
    TreeHelper& tree_helper;
    prot_filer::CompressedWriter writer;
};

class FragmentsWriterHelper
{
public:
    FragmentsWriterHelper(TreeHelper& tree_helper, FullCachedAnglesSeqReader* reader);//Adapter
    ~FragmentsWriterHelper();
    void write();
private:
    TreeHelper& tree_helper;
    const FullCachedAnglesSeqReader* reader;
    prot_filer::FragmentsWriter writer;
};
#endif
