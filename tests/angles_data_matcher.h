#include <gmock/gmock.h>

using ::testing::Matcher;
using ::testing::MakeMatcher;
using ::testing::MatcherInterface;
using ::testing::MatchResultListener;

class IncompleteAnglesDataMatcher: public MatcherInterface<const AnglesData&>
{
public:
    explicit IncompleteAnglesDataMatcher(const IncompleteAnglesData& data) :
        expected_angles_data(data)
    {};

    virtual bool MatchAndExplain(const AnglesData& angles_data, MatchResultListener* listener) const
    {
        if (expected_angles_data.nres == angles_data.nres)
        {
            AngleIdPair* expected_pair = expected_angles_data.angles;
            AngleIdPair* data_pair = angles_data.angles;
            unsigned int i = 0;
            bool match;
            do
            {
                match = (expected_pair[i].fi == data_pair[i].fi) && (expected_pair[i].si == data_pair[i].si);
                ++i;
            }
            while ((i < angles_data.nres - 1) && match);
            return match;
        }
        else
        {
            return false;
        }
    }

    virtual void DescribeTo(::std::ostream* os) const
    {
        *os << "nres: " << expected_angles_data.nres << " (...)";
    }

    virtual void DescribeNegationTo(::std::ostream* os) const
    {
        *os << "angles data not match";
    }
private:
    const IncompleteAnglesData& expected_angles_data;
};

inline Matcher<const AnglesData&> IncompleteAnglesDataEq(const IncompleteAnglesData& data)
{
    return MakeMatcher(new IncompleteAnglesDataMatcher(data));
}

