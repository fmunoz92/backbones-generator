#include <vector>
#include <iostream>
#include <sstream>
#include <mili/mili.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "prot-filer/angles.h"
#include "petu.h"
#include "readdata.h"

using namespace std;
using namespace prot_filer;

TEST(TestReadData, read_data)
{
    TreeData tree_data(1, NULL);
    istringstream f("-60 -40\n");
    readdata(f, tree_data);
    tree_data.angles_data = new AnglesData(tree_data.nres, tree_data.angles_mapping);
    ASSERT_EQ(1, tree_data.angles_mapping->get_mapping_size());
    ASSERT_EQ(-60, tree_data.angles_mapping->get_fi_value(0));
    ASSERT_EQ(-40, tree_data.angles_mapping->get_si_value(0));
}
