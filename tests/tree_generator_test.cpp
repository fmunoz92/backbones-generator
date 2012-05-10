#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <iostream>
#include <mili/mili.h>
#include "prot-filer/angles.h"
#include "petu.h"
#include "tree_generator.h"
#include "utils.h"
#include "readdata.h"

using namespace std;
using namespace prot_filer;
using mili::insert_into;
using namespace testing;
using testing::A;
using testing::_;
using ::testing::Expectation;

class MockWriteHelper : public XtcWriterHelper
{
public:
    MOCK_METHOD1(write, void(TreeData& tree_data));
    virtual ~MockWriteHelper()
    {};
};

bool eq(const AnglesData& d1, const AnglesData& d2)
{
    AngleIdPair* expected_pair = d1.angles;
    AngleIdPair* data_pair = d2.angles;
    unsigned int i = 0;
    bool match;
    do
    {
        match = (expected_pair[i].fi == data_pair[i].fi) && (expected_pair[i].si == data_pair[i].si);
        ++i;
    }
    while ((i < d1.nres - 1) && match);
    return match;
}

MATCHER_P(CheckData, d, "")
{
    return eq(d, *(arg.angles_data));
}

TEST(Test, simple_generator)
{
    const unsigned int nres = 3;

    AnglesData d1(nres);
    d1.angles[0] = AngleIdPair(0, 0);
    d1.angles[1] = AngleIdPair(0, 0);
    AnglesData d2(nres);
    d2.angles[0] = AngleIdPair(2, 0);
    d2.angles[1] = AngleIdPair(0, 2);

    Grillado* grilla = new Grillado(100, 100, 100);
    TreeData tree_data(nres, grilla);
    istringstream f("-60  -40\n-60  140\n-130 140\n60   30");
    readdata(f, tree_data);
    tree_data.angles_data = new AnglesData(tree_data.nres, tree_data.angles_mapping);
    MockWriteHelper mock_helper;

    Expectation e1 = EXPECT_CALL(mock_helper, write(CheckData(d1))).Times(1);
    EXPECT_CALL(mock_helper, write(CheckData(d2))).Times(1).After(e1);

    SimpleTreeOperator st(tree_data);
    TreeGenerator g(tree_data, mock_helper, st);
    g.generate();
    ASSERT_EQ(2, tree_data.cont);
}
