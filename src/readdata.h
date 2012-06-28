#ifndef READDATA_H
#define READDATA_H

#include <sstream>
#include <fstream>
#include <string>

#include "tree_data.h"// for global typedef FullCachedAnglesSeqReader

FullCachedAnglesSeqReader* read_chains(const std::string& format, const std::string& input_file, const std::string& fragments_file);

#endif
