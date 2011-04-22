#include "petu.h"

FilterResultType poneres(float* pR, const unsigned int resN, TreeData& tree_data, Residuo& residuo, unsigned int si_index, unsigned int fi_index);

FilterResultType addChain(float* pR, unsigned int resN, TreeData& tree_data, vector<Residuo> &residuos, const IncompleteAnglesData& chain, unsigned int chain_index);