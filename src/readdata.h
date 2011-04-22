#include "petu.h"
#include <sstream>
#include <fstream>

void readdata(std::istream& filer, TreeData& tree_data);

FullCachedAnglesSeqReader* read_chains(const string& format, const string& input_file, const string& fragments_file);
