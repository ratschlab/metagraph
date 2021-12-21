#include <vector>

#include "gtest/gtest.h"

#include "cli/query.hpp"
#include "common/vector.hpp"


namespace {

using namespace mtg;

TEST(collapse_coord_ranges, sample) {
    std::vector<SmallVector<uint64_t>> tuples = { { 1, 2, 3 },
                                                  { 2, 3, 4 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1-2", "0-2-3", "0-3-4" }),
              cli::collapse_coord_ranges(tuples));
}

} // namespace
