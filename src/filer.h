#include "petu.h"
#include "prot-filer/compressed_filer.h"
#include "prot-filer/xtc_filer.h"

class WriterAdapter
{
public:
    virtual bool open(const string& name) = 0;
    virtual void write(BasicProtein& protein, const AnglesData& angles_data) = 0;
    virtual void close() = 0;
    virtual ~WriterAdapter() {};
};

class XtcWriterAdapter : public WriterAdapter, private XtcWriter
{
public:
    virtual bool open(const string& name)
    {
        return XtcWriter::open(name);
    }
    virtual void write(BasicProtein& protein, const AnglesData& angles_data)
    {
        XtcWriter::write(protein, angles_data);
    }
    void close()
    {
        XtcWriter::close();
    }
    virtual ~XtcWriterAdapter()
    {}
};

class CompressedWriterAdapter : public WriterAdapter, private CompressedWriter
{
public:
    virtual bool open(const string& name)
    {
        return CompressedWriter::open(name);
    }
    virtual void write(BasicProtein& protein, const AnglesData& angles_data)
    {
        CompressedWriter::write(angles_data);
    }
    void close()
    {
        CompressedWriter::close();
    }
    virtual ~CompressedWriterAdapter()
    {}
};
