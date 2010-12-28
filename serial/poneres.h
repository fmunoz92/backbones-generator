#include "petu.h"

FilterResultType poneres(float* pR, const unsigned int resN, TreeData* arbol_data, Residuo&  residuo, unsigned int si_index, unsigned int fi_index);

FilterResultType addNRes(float* pR, unsigned int resN, TreeData* tree_data, vector<Residuo> &residuos, const IncompleteAnglesData& residue_chain);
