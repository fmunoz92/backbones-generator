#include "petu.h"
#include <sstream>
#include <fstream>

void readdata(std::ifstream& filer, std::vector<float> &cosfi, std::vector<float> &sinfi, std::vector<float> &cossi, std::vector<float> &sinsi, AnglesMapping* angles_info);

AnglesDatabase* read_chains(const string& input_file, const string& read_format, const string& cache_type);
