#include <vector>

#include "gtest/gtest.h"

#include "cli/query.hpp"
#include "common/vector.hpp"


namespace {

using namespace mtg;
using namespace mtg::cli;

TEST(get_collapsed_coord_ranges, sample) {
    SmallVector<uint64_t> v1;
    v1.push_back(1);
    v1.push_back(2);
    v1.push_back(3);

    SmallVector<uint64_t> v2;
    v1.push_back(1);
    v1.push_back(3);

    std::vector<SmallVector<uint64_t> v = { v1, v2 };
}

}