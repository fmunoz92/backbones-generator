#include <vector>
#include <iostream>
#include <sstream>
#include <mili/mili.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "tree_data.h"

using namespace std;

TEST(TestReadData, read_data)
{
	istringstream f("-60 -40\n");
	TreeData tree_data(1, 4, 4, 4, f);

    ASSERT_EQ(1, tree_data.angles_mapping->get_mapping_size());
    ASSERT_EQ(-60, tree_data.angles_mapping->get_fi_value(0));
    ASSERT_EQ(-40, tree_data.angles_mapping->get_si_value(0));
}
