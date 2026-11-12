#include <algorithm>
#include <limits>

#include "gtest/gtest.h"
#include "coordinate_test_graph.hpp"
#include "dense_search_values.hpp"

namespace mtg::test {
namespace {
// The records of |coords| (every record get_kmer_coordinates returns, num_top_labels unlimited)
// as get_top_labels returns them: the score reported (the coordinates matched if |weighted|,
// the query positions matched otherwise), the top |top| in the contract's order (C5,
// SPEC-multigraph-envelope.md section 3: score descending, then the header's bytes; a
// coordinate annotation holds no abundances)
template <class Coordinates>
std::vector<std::pair<std::string, size_t>>
contract_top(const Coordinates &coords, size_t top, bool weighted) {
    std::vector<std::pair<std::string, size_t>> rows;
    for (const auto &[label, count, positions, length] : coords) {
        size_t value = count;
        if (weighted) {
            value = 0;
            for (const auto &tuple : positions)
                value += tuple.size();
        }
        rows.emplace_back(label, value);
    }
    std::stable_sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
        return a.second != b.second ? a.second > b.second : a.first < b.first;
    });
    rows.resize(std::min(rows.size(), top));
    return rows;
}
} // namespace

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
                // every record above the thresholds
                auto all = dense(graph->get_kmer_coordinates(
                        query, std::numeric_limits<size_t>::max(), discovery, presence));
                for (size_t top : {size_t(0), size_t(1), size_t(2), size_t(100)}) {
                    // the top N by the score reported, with the counts by the coordinates
                    // matched (the contract's order: a record is not selected by its k-mer
                    // matches and then reported with its coordinates)
                    for (bool weighted : {false, true}) {
                        EXPECT_EQ(contract_top(all, top, weighted),
                                  graph->get_top_labels(query, top, discovery, presence, weighted));
                    }
                    // without the counts, the records and the order get_kmer_coordinates returns
                    auto coords = dense(graph->get_kmer_coordinates(query, top, discovery, presence));
                    std::vector<std::pair<std::string, size_t>> expected;
                    for (const auto &[label, count, positions, length] : coords)
                        expected.emplace_back(label, count);
                    EXPECT_EQ(expected, graph->get_top_labels(query, top, discovery, presence));
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
            auto coords = dense(graph->get_kmer_coordinates(query, 100, 0, 0));
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

TEST(HeaderCounts, SelectsByTheReportedScore) {
    auto graph = coordinate_test_graph(3, {{"AAAA", "F", "short"}, {"AAAAAAAA", "F", "long"}});
    // Both records have one matched query position, "long" six coordinates of AAA and "short"
    // two. The top N is taken by the score reported (C5, SPEC-multigraph-envelope.md
    // section 3): with the counts, by the coordinates matched; without them, the tie of the
    // positions matched by the header's bytes
    EXPECT_EQ((std::vector<std::pair<std::string, size_t>>{{"long", 6}}),
              graph->get_top_labels("AAA", 1, 0, 0, true));
    EXPECT_EQ((std::vector<std::pair<std::string, size_t>>{{"long", 6}, {"short", 2}}),
              graph->get_top_labels("AAA", 2, 0, 0, true));
    EXPECT_EQ((std::vector<std::pair<std::string, size_t>>{{"long", 1}}),
              graph->get_top_labels("AAA", 1, 0, 0, false));
}

TEST(HeaderCounts, OneHeaderInTwoColumnsTiedAtTheCut) {
    // the header "same" in two columns: in F0 (column 0) the record ACGT, every query k-mer
    // once; in F1 (column 1) the record ACGTACGT, every query k-mer twice. The query ACGT
    // (k-mers ACG, CGT) matches both records at both positions: the same k-mer score and the
    // same header, so the tie at a top-1 cut is broken by (column, seq_id), and column 0's
    // record is kept; with the counts the score is the coordinates matched (2 against 4) and
    // column 1's record ranks first. Every expectation is written out by hand.
    auto graph = coordinate_test_graph(3, {{"ACGT", "F0", "same"},
                                            {"ACGTACGT", "F1", "same"}});
    using Counts = std::vector<std::tuple<std::string, size_t, std::vector<size_t>>>;
    using Pairs = std::vector<std::pair<std::string, size_t>>;
    graph::SearchCut cut;

    // top 1, without the counts: column 0's record (one coordinate per position)
    EXPECT_EQ((Counts{{"same", 2, {1, 1}}}), dense(graph->get_kmer_counts("ACGT", 1, 0, 0, &cut)));
    EXPECT_EQ(1u, cut.returned);
    EXPECT_FALSE(cut.complete);
    EXPECT_EQ(2u, cut.nth_score);
    EXPECT_EQ(2u, cut.tied);
    EXPECT_EQ(1u, cut.tied_returned);
    EXPECT_EQ((Pairs{{"same", 2}}), graph->get_top_labels("ACGT", 1, 0, 0, false, &cut));
    EXPECT_EQ(2u, cut.tied);
    EXPECT_EQ(1u, cut.tied_returned);
    auto coords = dense(graph->get_kmer_coordinates("ACGT", 1, 0, 0, &cut));
    ASSERT_EQ(1u, coords.size());
    EXPECT_EQ("same", std::get<0>(coords[0]));
    EXPECT_EQ(2u, std::get<1>(coords[0]));
    // ACG at its record's position 0, CGT at 1; the record has two k-mers
    EXPECT_EQ((std::vector<SmallVector<uint64_t>>{{0}, {1}}), std::get<2>(coords[0]));
    EXPECT_EQ(2u, std::get<3>(coords[0]));
    EXPECT_EQ(2u, cut.tied);
    EXPECT_EQ(1u, cut.tied_returned);

    // top 2: both records, column 0's first; nothing cut
    EXPECT_EQ((Counts{{"same", 2, {1, 1}}, {"same", 2, {2, 2}}}),
              dense(graph->get_kmer_counts("ACGT", 2, 0, 0, &cut)));
    EXPECT_EQ(2u, cut.returned);
    EXPECT_TRUE(cut.complete);

    // with the counts: the coordinates matched rank, column 1's record first, no tie
    EXPECT_EQ((Pairs{{"same", 4}}), graph->get_top_labels("ACGT", 1, 0, 0, true, &cut));
    EXPECT_EQ(1u, cut.returned);
    EXPECT_FALSE(cut.complete);
    EXPECT_EQ(4u, cut.nth_score);
    EXPECT_EQ(1u, cut.tied);
    EXPECT_EQ(1u, cut.tied_returned);
    EXPECT_EQ((Pairs{{"same", 4}, {"same", 2}}), graph->get_top_labels("ACGT", 2, 0, 0, true));
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
                            auto coords = dense(graph->get_kmer_coordinates(query, top, discovery, presence));
                            decltype(dense(graph->get_kmer_counts(query, top, discovery, presence))) expected_counts;
                            decltype(dense(graph->get_top_label_signatures(query, top, discovery, presence))) expected_masks;
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
                            EXPECT_EQ(expected_counts, dense(graph->get_kmer_counts(query, top, discovery, presence)));
                            EXPECT_EQ(expected_masks, dense(graph->get_top_label_signatures(query, top, discovery, presence)));
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
    auto counts = dense(graph->get_kmer_counts(query, 100, 0, 0));
    std::sort(counts.begin(), counts.end());
    ASSERT_EQ(counts.size(), 2u);
    EXPECT_EQ((std::tuple<std::string, size_t, std::vector<size_t>>{
        "first", 2, {3, 0, 0, 0, 0, 0, 3}}), counts[0]);
    EXPECT_EQ((std::tuple<std::string, size_t, std::vector<size_t>>{
        "second", 2, {1, 0, 0, 0, 0, 0, 1}}), counts[1]);
    const auto signatures = dense(graph->get_top_label_signatures(query, 100, 0, 0));
    ASSERT_EQ(signatures.size(), 2u);
    for (const auto &[label, count, mask] : signatures) {
        EXPECT_EQ(count, 2u);
        ASSERT_EQ(mask.size(), 7u);
        for (size_t i = 0; i < mask.size(); ++i)
            EXPECT_EQ(mask[i], i == 0 || i == 6);
    }
}
} // namespace mtg::test
