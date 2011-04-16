#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <iostream>
#include <mili/mili.h>
#include "prot-filer/angles.h"
#include "petu.h"
#include "tree_generator.h"
#include "utils.h"
#include "readdata.h"
#include "angles_data_matcher.h"

using namespace std;
using namespace prot_filer;
using mili::insert_into;
using namespace testing;
using testing::A;
using testing::_;
using ::testing::Expectation;

class MockWriteAdapter : public WriterAdapter
{
public:
    MOCK_METHOD1(open, bool(const string&));
    MOCK_METHOD2(write, void(BasicProtein& protein, const AnglesData& angles_data));
    MOCK_METHOD0(close, void(void));
    virtual ~MockWriteAdapter()
    {};
};

TEST(Test, simple_generator)
{
    const unsigned int nres = 3;

    IncompleteAnglesData d1(nres);
    d1.angles[0] = AngleIdPair(0, 0);
    d1.angles[1] = AngleIdPair(0, 0);
    IncompleteAnglesData d2(nres);
    d2.angles[0] = AngleIdPair(2, 0);
    d2.angles[1] = AngleIdPair(0, 2);

    Grillado* grilla = new Grillado(100, 100, 100);
    TreeData tree_data(nres, grilla);
    istringstream f("-60  -40\n-60  140\n-130 140\n60   30");
    readdata(f, tree_data);
    tree_data.angles_data = new AnglesData(tree_data.nres, *tree_data.angles_mapping);
    MockWriteAdapter* mock_writer = new MockWriteAdapter();

    EXPECT_CALL(*mock_writer, open(A<const string&>())).Times(1);
    Expectation e1 = EXPECT_CALL(*mock_writer, write(_, IncompleteAnglesDataEq(d1))).Times(1);
    EXPECT_CALL(*mock_writer, write(_, IncompleteAnglesDataEq(d2))).Times(1).After(e1);
    EXPECT_CALL(*mock_writer, close()).Times(1);

    SimpleTreeGenerator generator(tree_data, mock_writer);
    generator.generate();
    ASSERT_EQ(2, tree_data.cont);
}

