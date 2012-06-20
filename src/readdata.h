#ifndef READDATA_H
#define READDATA_H

#include <sstream>
#include <fstream>

#include "tree_data.h"

FullCachedAnglesSeqReader* read_chains(const string& format, const string& input_file, const string& fragments_file);

#endif
