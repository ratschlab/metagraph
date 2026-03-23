#include <vector>
#include <unordered_set>

#include "gtest/gtest.h"

#include "cli/query.hpp"
#include "common/vector.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg::cli {
/**
 * Split long contigs into overlapping segments for better work balancing.
 */
void split_contigs_for_rebalancing(size_t k,
                                   size_t kmers_per_seq,
                                   std::vector<std::pair<std::string, std::vector<uint64_t>>> *contigs);
} // namespace mtg::cli

namespace {

using namespace mtg;

TEST(collapse_coord_ranges, empty) {
    std::vector<SmallVector<uint64_t>> tuples = {};

    EXPECT_EQ(std::vector<std::string>(),
              cli::collapse_coord_ranges(tuples));
}

TEST(collapse_coord_ranges, single) {
    std::vector<SmallVector<uint64_t>> tuples = { { 1, 3, 5 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1", "0-3", "0-5" }),
              cli::collapse_coord_ranges(tuples));
}

TEST(collapse_coord_ranges, standard) {
    std::vector<SmallVector<uint64_t>> tuples = { { 1, 2, 3 },
                                                  { 2, 3, 4 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1-2", "0-2-3", "0-3-4" }),
              cli::collapse_coord_ranges(tuples));

    // test longer
    std::vector<SmallVector<uint64_t>> longer = { { 1, 2, 3 },
                                                  { 2, 3, 4 },
                                                  { 3, 4, 5 },
                                                  { 4, 5, 6 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1-4", "0-2-5", "0-3-6" }),
              cli::collapse_coord_ranges(longer));
}

TEST(collapse_coord_ranges, staggered) {
    std::vector<SmallVector<uint64_t>> tuples = { { 1, 2, 3 },
                                                  { 3, 4, 5 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1", "0-2-3", "0-3-4", "1-5" }),
              cli::collapse_coord_ranges(tuples));
}

TEST(collapse_coord_ranges, disjoint) {
    std::vector<SmallVector<uint64_t>> tuples = { { 1, 2, 3 },
                                                  { 5, 6, 7 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1", "0-2", "0-3", "1-5", "1-6", "1-7" }),
              cli::collapse_coord_ranges(tuples));

    std::vector<SmallVector<uint64_t>> longer = { { 1, 2, 3 },
                                                  { 5, 6, 7 },
                                                  { 6, 7, 8 },
                                                  { 7, 8, 9 } };

    EXPECT_EQ(std::vector<std::string>({ "0-1", "0-2", "0-3", "1-5-7", "1-6-8", "1-7-9" }),
              cli::collapse_coord_ranges(longer));
}


void expect_contig_string_set_eq(
        const std::vector<std::pair<std::string, std::vector<uint64_t>>> &contigs,
        std::initializer_list<std::string> expected) {
    auto contig_strings = utils::get_firsts<std::vector<std::string>>(contigs);
    EXPECT_EQ(std::unordered_set<std::string>(expected.begin(), expected.end()),
              std::unordered_set<std::string>(contig_strings.begin(), contig_strings.end()));
    for (const auto &[_, path] : contigs) {
        EXPECT_TRUE(path.empty());
    }
}

TEST(split_contigs_for_rebalancing, no_split_when_path_fits) {
    std::vector<std::pair<std::string, std::vector<uint64_t>>> contigs = {
        { "abcdefgh", {} } // for k=4, path size is 5
    };

    cli::split_contigs_for_rebalancing(4, 5, &contigs);

    ASSERT_EQ(1u, contigs.size());
    expect_contig_string_set_eq(contigs, { "abcdefgh" });
}

TEST(split_contigs_for_rebalancing, splits_and_appends_extra_segments) {
    std::vector<std::pair<std::string, std::vector<uint64_t>>> contigs = {
        { "abcdefghijklmnop", {} },
        { "qrstuvwxyz", {} }
    };

    cli::split_contigs_for_rebalancing(4, 5, &contigs);

    ASSERT_EQ(5u, contigs.size());
    expect_contig_string_set_eq(contigs, { "abcdefgh", "qrstuvwx", "fghijklm", "klmnop", "vwxyz" });
}

TEST(split_contigs_for_rebalancing, split_when_exceeds_by_one_node) {
    std::vector<std::pair<std::string, std::vector<uint64_t>>> contigs = {
        { "abcdefghi", {} } // k=4, path size = 6 (> max_path_size=5 by one)
    };

    cli::split_contigs_for_rebalancing(4, 5, &contigs);

    ASSERT_EQ(2u, contigs.size());
    expect_contig_string_set_eq(contigs, { "abcdefgh", "fghi" });
}

TEST(split_contigs_for_rebalancing, k_equals_one_has_no_overlap) {
    std::vector<std::pair<std::string, std::vector<uint64_t>>> contigs = {
        { "abcdefghijk", {} }
    };

    cli::split_contigs_for_rebalancing(1, 4, &contigs);

    ASSERT_EQ(3u, contigs.size());
    expect_contig_string_set_eq(contigs, { "abcd", "efgh", "ijk" });
}

TEST(split_contigs_for_rebalancing, exact_multiple_of_chunk_steps) {
    std::vector<std::pair<std::string, std::vector<uint64_t>>> contigs = {
        { "abcdefghijklmnopq", {} } // k=3, path size=15, max_path_size=5
    };

    cli::split_contigs_for_rebalancing(3, 5, &contigs);

    ASSERT_EQ(3u, contigs.size());
    expect_contig_string_set_eq(contigs, { "abcdefg", "fghijkl", "klmnopq" });
}

} // namespace
