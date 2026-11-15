#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "graph/search_cut.hpp"

namespace {
using Rows = std::vector<std::pair<size_t, size_t>>;

TEST(SearchCut, PreservesExistingSelectionAndCountsTies) {
    std::mt19937 random(31979);
    for (size_t size = 0; size < 80; ++size) {
        Rows input;
        for (size_t i = 0; i < size; ++i)
            input.emplace_back(i, 1 + random() % 9);
        std::shuffle(input.begin(), input.end(), random);
        for (size_t limit = 0; limit <= size + 1; ++limit) {
            auto expected = input;
            auto rank = [](const auto &a, const auto &b) {
                return std::make_pair(b.second, a.first) < std::make_pair(a.second, b.first);
            };
            if (size > limit) {
                std::nth_element(expected.begin(), expected.begin() + limit, expected.end(), rank);
                expected.resize(limit);
                std::sort(expected.begin(), expected.end(), rank);
            }
            auto actual = input;
            mtg::graph::SearchCut cut;
            mtg::graph::top_n_sorted(actual, limit, &cut);
            ASSERT_EQ(expected, actual);
            ASSERT_EQ(size, cut.matched);
            ASSERT_EQ(std::min(limit, size), cut.returned);
            if (size <= limit || limit == 0) {
                EXPECT_FALSE(cut.nth_score);
                EXPECT_EQ(0u, cut.tied);
                EXPECT_EQ(0u, cut.tied_returned);
            } else {
                const auto score = expected.back().second;
                EXPECT_EQ(score, *cut.nth_score);
                auto tied = [score](auto row) { return row.second == score; };
                EXPECT_EQ(size_t(std::count_if(input.begin(), input.end(), tied)), cut.tied);
                EXPECT_EQ(size_t(std::count_if(actual.begin(), actual.end(), tied)), cut.tied_returned);
            }
        }
    }
}

TEST(SearchCut, HeaderIdentityPairsAndReusedOutput) {
    std::vector<std::pair<std::pair<size_t, size_t>, size_t>> rows{
        {{1, 9}, 5}, {{0, 7}, 5}, {{0, 2}, 5}, {{3, 0}, 4}};
    mtg::graph::SearchCut cut;
    mtg::graph::top_n_sorted(rows, 2, &cut);
    EXPECT_EQ(std::make_pair(size_t(0), size_t(2)), rows.front().first);
    EXPECT_EQ(3u, cut.tied);
    EXPECT_EQ(2u, cut.tied_returned);
    rows.clear();
    mtg::graph::top_n_sorted(rows, 2, &cut);
    EXPECT_EQ(0u, cut.matched);
    EXPECT_FALSE(cut.nth_score);
}
} // namespace
