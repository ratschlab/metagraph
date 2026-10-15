#include "gtest/gtest.h"
#include "coordinate_test_graph.hpp"

namespace mtg::test {
TEST(HeaderCounts, CountsMixedKmersAndTheirOccurrences) {
    auto graph = coordinate_test_graph(3, {{"ACGACGACG", "F", "first"},
                                            {"ACGTTT", "F", "second"},
                                            {"ACGACG", "F", "third"}});
    // ACGACG contains ACG twice, CGA once and GAC once. A record matches
    // each query position once, regardless of how many target occurrences
    // it has. Abundance counts all of those occurrences instead.
    const std::vector<std::pair<std::string, size_t>> counts{
        {"first", 4}, {"second", 2}, {"third", 4}
    };
    const std::vector<std::pair<std::string, size_t>> abundances{
        {"first", 10}, {"second", 2}, {"third", 6}
    };
    for (bool weighted : {false, true}) {
        auto actual = graph->get_top_labels("ACGACG", 100, 0, 0, weighted);
        std::sort(actual.begin(), actual.end());
        EXPECT_EQ(weighted ? abundances : counts, actual);
    }
    auto filtered = graph->get_top_labels("ACGACG", 100, 0.75, 0);
    std::sort(filtered.begin(), filtered.end());
    EXPECT_EQ((std::vector<std::pair<std::string, size_t>>{{"first", 4}, {"third", 4}}),
              filtered);
}

TEST(HeaderCounts, MatchesCoordinatesForThresholdsLimitsAndRepeats) {
    auto graph = coordinate_test_graph(3, {{"ACGACGACG", "file1", "same header"},
                                            {"ACGACG", "file1", "same header"},
                                            {"ACGTTT", "file2", "third"}});
    for (const auto &query : {"ACGACGACG", "ACGTTT", "ACGNNNACG", "TTTACG", "TTTTTT", "AC"}) {
        for (double discovery : {0., 0.5, 1.}) {
            for (double presence : {0., 0.5, 1.}) {
                for (size_t top : {size_t(0), size_t(1), size_t(2), size_t(100)}) {
                    auto coords = graph->get_kmer_coordinates(query, top, discovery, presence);
                    for (bool weighted : {false, true}) {
                        std::vector<std::pair<std::string, size_t>> expected;
                        for (const auto &[label, count, positions, length] : coords) {
                            size_t value = count;
                            if (weighted) {
                                value = 0;
                                for (const auto &tuple : positions)
                                    value += tuple.size();
                            }
                            expected.emplace_back(label, value);
                        }
                        EXPECT_EQ(expected, graph->get_top_labels(query, top, discovery, presence, weighted));
                    }
                }
                auto counts = graph->get_top_labels(query, 100, discovery, presence);
                std::vector<std::string> expected;
                for (const auto &[label, count] : counts)
                    expected.push_back(label);
                EXPECT_EQ(expected, graph->get_labels(query, discovery, presence));
            }
        }
    }
}

TEST(HeaderCounts, NonAnchorRowDiffCoordinates) {
    std::vector<CoordinateRecord> records{{"ACGACGACGTT", "F", "first"},
                                          {"ACGTTACG", "F", "second"}};
    auto plain = coordinate_test_graph(3, records);
    auto diffs = coordinate_test_graph(3, records, true);
    for (const auto &query : {"ACGACGACGTT", "ACGTT", "TTTTACG"}) {
        for (bool weighted : {false, true}) {
            for (size_t top : {size_t(1), size_t(100)})
                EXPECT_EQ(plain->get_top_labels(query, top, 0, 0, weighted),
                          diffs->get_top_labels(query, top, 0, 0, weighted));
        }
    }
}

TEST(HeaderCounts, BRWTRowsAndRepeatedQueryPositions) {
    std::vector<CoordinateRecord> records{{"ACGACGACGTT", "F0", "first"},
                                          {"ACGTTACG", "F1", "second"},
                                          {"ACGACGTTT", "F2", "third"}};
    for (bool row_diff : {false, true}) {
        auto graph = coordinate_test_graph(3, records, row_diff, true);
        for (const auto &query : {"ACGACGACGTT", "ACGTT", "TTTTACG"}) {
            auto coords = graph->get_kmer_coordinates(query, 100, 0, 0);
            for (bool weighted : {false, true}) {
                std::vector<std::pair<std::string, size_t>> expected;
                for (const auto &[label, count, positions, length] : coords) {
                    size_t value = count;
                    if (weighted) {
                        value = 0;
                        for (const auto &tuple : positions)
                            value += tuple.size();
                    }
                    expected.emplace_back(label, value);
                }
                EXPECT_EQ(expected, graph->get_top_labels(query, 100, 0, 0, weighted));
            }
        }
    }
}

TEST(HeaderCounts, PreservesSelectionBeforeAbundanceWeighting) {
    auto graph = coordinate_test_graph(3, {{"AAAA", "F", "short"}, {"AAAAAAAA", "F", "long"}});
    // Both records have one matched query position. The existing tie break
    // selects the first record before replacing that count with its abundance.
    const std::vector<std::pair<std::string, size_t>> expected{{"short", 2}};
    EXPECT_EQ(expected, graph->get_top_labels("AAA", 1, 0, 0, true));
}

TEST(HeaderCounts, CountsQueryPositionsInsteadOfDistinctKmers) {
    auto graph = coordinate_test_graph(3, {{"AAAAAAAA", "F", "long"}, {"AAAA", "F", "short"}});
    const auto counts = graph->get_top_labels("AAAAA", 2, 1., 1.);
    const auto weighted = graph->get_top_labels("AAAAA", 2, 1., 1., true);
    ASSERT_EQ(counts.size(), 2u);
    for (const auto &[label, count] : counts)
        EXPECT_EQ(count, 3u);
    for (const auto &[label, count] : weighted)
        EXPECT_EQ(count, label == "long" ? 18u : 6u);
}
} // namespace mtg::test
