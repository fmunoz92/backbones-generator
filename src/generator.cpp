#include "generator.h"
#include "tree_generator.h"
#include "filer.h"

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
