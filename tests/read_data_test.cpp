#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <mili/mili.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "backbones-generator/tree_data.h"

using namespace std;

TEST(TestReadData, read_data)
{
    istringstream inputFile("-60 -40\n");

    string outputFile("");

    TreeData treeData(1, 4, 4, 4, inputFile, outputFile);

    ASSERT_EQ(1, treeData.anglesMapping.get_mapping_size());
    ASSERT_EQ(-60, treeData.anglesMapping.get_fi_value(0));
    ASSERT_EQ(-40, treeData.anglesMapping.get_si_value(0));
}
