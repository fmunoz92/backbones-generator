#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <mili/mili.h>
#include <gtest/gtest.h>


#include "backbones-generator/tree_data.h"
#include "backbones-generator/tree_filters.h"
#include "backbones-generator/grillado.h"

using namespace std;

TEST(TestReadData, read_data)
{
    stringstream inputFile("-60 -40\n0 90\n");

    TreeFilters treeFilters;
    Grillado grilla(10, 10, 10);
    prot_filer::AnglesMapping anglesMapping(10);
    prot_filer::AnglesData anglesData(10, &anglesMapping); 
    IncrementalBackbone incrementalBackbone(10, grilla, anglesData, anglesMapping, treeFilters);
    TreeData treeData(10, incrementalBackbone);
    BareBackbone::treeData = &treeData;
    readData(inputFile, treeData, anglesMapping);

    ASSERT_EQ(2, anglesMapping.get_mapping_size());
    ASSERT_EQ(-60, anglesMapping.get_fi_value(0));
    ASSERT_EQ(-40, anglesMapping.get_si_value(0));

    ASSERT_EQ(0, anglesMapping.get_fi_value(1));
    ASSERT_EQ(90, anglesMapping.get_si_value(1));
}
