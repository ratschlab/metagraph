#include <gtest/gtest.h>

#include "common/configured_pair_id.hpp"

namespace {
using mtg::common::configured_pair_id;

TEST(ConfiguredPairId, PythonUuid5Vectors) {
    // Shared with Python uuid.uuid5; lengths count UTF-8 bytes, not code points.
    EXPECT_EQ("n1-0c8ac32d-6ec3-5dc6-85c0-6993b9c02020",
              configured_pair_id("refseq", "/graphs/a.dbg", "/annotations/a.annodbg"));
    EXPECT_EQ("n1-f3db2c00-cc7c-5044-a87d-a679fd98328e",
              configured_pair_id("refseq", "/graphs/a.dbg", "/annotations/b.annodbg"));
    EXPECT_EQ("n1-b226283c-5004-52d2-a5f4-ccc494acf745",
              configured_pair_id("réf", "/图/a", ""));
    EXPECT_EQ("n1-d600c46d-30c3-5f3d-a520-7d231fcb1141",
              configured_pair_id("refseq", std::nullopt, std::nullopt));
    EXPECT_EQ("n1-5ac59ea8-6d4a-5b01-aa08-e763a91371d2",
              configured_pair_id("g:12", "12:/a:b", "/b:34"));
}

TEST(ConfiguredPairId, DistinguishesFieldsAndAbsence) {
    EXPECT_NE(configured_pair_id("g", "ab", "c"), configured_pair_id("g", "a", "bc"));
    EXPECT_NE(configured_pair_id("g", "a", "b"), configured_pair_id("h", "a", "b"));
    EXPECT_NE(configured_pair_id("g", "a", "b"), configured_pair_id("g", "c", "b"));
    EXPECT_NE(configured_pair_id("g", "", "b"), configured_pair_id("g", std::nullopt, "b"));
}
} // namespace
