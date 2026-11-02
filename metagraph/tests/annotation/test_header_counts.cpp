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

TEST(HeaderCounts, PositionalResultsMatchCoordinatesAcrossBackends) {
    const std::vector<CoordinateRecord> records{
        {"ACGACGACGTT", "F0", "same header"},
        {"ACGTTACG", "F0", "same header"},
        {"TTTACGACG", "F1", "third"},
        {"ACGACGTTT", "F2", "fourth"}
    };
    for (bool row_diff : {false, true}) {
        for (bool brwt : {false, true}) {
            auto graph = coordinate_test_graph(3, records, row_diff, brwt);
            for (const auto &query : {"ACGACGACGTT", "ACGNNNACG", "NNNACG", "ACGNNN",
                                      "TTTTTT", "NNNNN", "AC"}) {
                for (double discovery : {0., 0.5, 1.}) {
                    for (double presence : {0., 0.5, 1.}) {
                        for (size_t top : {size_t(0), size_t(1), size_t(2), size_t(100)}) {
                            SCOPED_TRACE(::testing::Message() << row_diff << ' ' << brwt
                                         << ' ' << query << ' ' << discovery << ' '
                                         << presence << ' ' << top);
                            auto coords = graph->get_kmer_coordinates(query, top, discovery, presence);
                            decltype(graph->get_kmer_counts(query, top, discovery, presence)) expected_counts;
                            decltype(graph->get_top_label_signatures(query, top, discovery, presence)) expected_masks;
                            for (const auto &[label, count, positions, length] : coords) {
                                std::vector<size_t> counts;
                                sdsl::bit_vector mask(positions.size(), 0);
                                for (const auto &tuple : positions) {
                                    mask[counts.size()] = !tuple.empty();
                                    counts.push_back(tuple.size());
                                }
                                expected_counts.emplace_back(label, count, std::move(counts));
                                expected_masks.emplace_back(label, count, std::move(mask));
                            }
                            EXPECT_EQ(expected_counts, graph->get_kmer_counts(query, top, discovery, presence));
                            EXPECT_EQ(expected_masks, graph->get_top_label_signatures(query, top, discovery, presence));
                        }
                    }
                }
            }
        }
    }
}

TEST(HeaderCounts, PositionalResultsPreserveGapsAndRepeatedKmers) {
    auto graph = coordinate_test_graph(3, {{"ACGACGACG", "F", "first"},
                                          {"ACGTTT", "F", "second"}});
    // Query positions are ACG, CGN, GNN, NNN, NNA, NAC, ACG.
    const std::string query = "ACGNNNACG";
    auto counts = graph->get_kmer_counts(query, 100, 0, 0);
    std::sort(counts.begin(), counts.end());
    ASSERT_EQ(counts.size(), 2u);
    EXPECT_EQ((std::tuple<std::string, size_t, std::vector<size_t>>{
        "first", 2, {3, 0, 0, 0, 0, 0, 3}}), counts[0]);
    EXPECT_EQ((std::tuple<std::string, size_t, std::vector<size_t>>{
        "second", 2, {1, 0, 0, 0, 0, 0, 1}}), counts[1]);
    const auto signatures = graph->get_top_label_signatures(query, 100, 0, 0);
    ASSERT_EQ(signatures.size(), 2u);
    for (const auto &[label, count, mask] : signatures) {
        EXPECT_EQ(count, 2u);
        ASSERT_EQ(mask.size(), 7u);
        for (size_t i = 0; i < mask.size(); ++i)
            EXPECT_EQ(mask[i], i == 0 || i == 6);
    }
}
} // namespace mtg::test
