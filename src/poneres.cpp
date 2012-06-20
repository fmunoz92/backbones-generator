#include "poneres.h"

ClashFilter::ClashFilter(const TreeData& tree_data) :
    tree_data(tree_data)
{}

bool ClashFilter::operator()(unsigned int index, const Atoms& patm, int at) const
{
    const bool clash = isclash(patm, at) == FILTER_FAIL;
    return !clash && (index != 2 || (islong(patm, at, tree_data.dmax2) != FILTER_FAIL));
}

