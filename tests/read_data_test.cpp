#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <mili/mili.h>
#include <gtest/gtest.h>


#include "backbones-generator/tree_data.h"

using namespace std;

TEST(TestReadData, read_data)
{
    stringstream inputFile("-60 -40\n0 90\n");

    TreeFilters treeFilters;
    IncrementalBackbone incrementalBackbone(0, treeFilters);
    TreeData treeData(1, 4, 4, 4, incrementalBackbone);
    treeData.readData(inputFile);

    ASSERT_EQ(2, treeData.anglesMapping.get_mapping_size());
    ASSERT_EQ(-60, treeData.anglesMapping.get_fi_value(0));
    ASSERT_EQ(-40, treeData.anglesMapping.get_si_value(0));

    ASSERT_EQ(0, treeData.anglesMapping.get_fi_value(1));
    ASSERT_EQ(90, treeData.anglesMapping.get_si_value(1));
}
