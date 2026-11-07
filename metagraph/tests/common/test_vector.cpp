#include "gtest/gtest.h"

#include "common/vector.hpp"


namespace {

template <class Container>
void check_growth_and_shrink() {
    Container values;
    for (int i = 0; i < 10000; ++i)
        values.push_back(i);

    values.resize(1000);
    values.shrink_to_fit();
    for (int i = 1000; i < 20000; ++i)
        values.push_back(i);

    for (int i = 0; i < 20000; ++i)
        ASSERT_EQ(i, values[i]);
}

TEST(Vector, GrowthAndShrink) {
    check_growth_and_shrink<Vector<int>>();
    check_growth_and_shrink<SmallVector<int>>();
}

} // namespace
