#include "gtest/gtest.h"
#include "tests/annotation/coordinate_test_graph.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/traversal/context_cursor.hpp"
#include <map>
#include <random>

#if defined(_DNA_GRAPH)
namespace mtg::graph::traversal {
namespace {
using test::coordinate_test_graph;
using Stop = ContextCursor::StopReason;
using Closure = ContextCursor::Closure;

TEST(ContextCursor, PreservesJointFlanksAfterSelection) {
    auto graph = coordinate_test_graph(3, {{"CACGAG", "F", "first"},
                                            {"TACGAA", "F", "second"}});
    ContextCursor cursor(*graph, "ACGA");
    EXPECT_EQ(cursor.close(100), (Closure{0, 0, Stop::BRANCH, Stop::BRANCH}));
    auto choices = cursor.choices();
    ASSERT_EQ(choices.size(), 2u);
    EXPECT_EQ(choices[0].left, 'C'); EXPECT_EQ(choices[0].right, 'G');
    EXPECT_EQ(choices[1].left, 'T'); EXPECT_EQ(choices[1].right, 'A');
    EXPECT_TRUE(cursor.extend(ContextCursor::Side::LEFT, 'C'));
    EXPECT_EQ(cursor.close(100), (Closure{0, 1, Stop::RECORD_END, Stop::RECORD_END}));
    EXPECT_EQ(cursor.sequence(), "CACGAG");
    ASSERT_EQ(cursor.occurrences().size(), 1u);
    EXPECT_EQ(cursor.occurrences()[0].record, 0u);
    EXPECT_EQ(cursor.occurrences()[0].offset, 0u);
    EXPECT_FALSE(cursor.extend(ContextCursor::Side::RIGHT, 'A'));
    EXPECT_EQ(cursor.sequence(), "CACGAG");
}

TEST(ContextCursor, CopiesCanExploreAlternativesIndependently) {
    auto graph = coordinate_test_graph(3, {{"CACGAG", "F", "one"}, {"TACGAA", "F", "two"}});
    ContextCursor original(*graph, "ACGA");
    original.close(100);
    auto first = original;
    auto second = original;
    EXPECT_TRUE(first.extend(ContextCursor::Side::LEFT, 'C'));
    EXPECT_TRUE(second.extend(ContextCursor::Side::LEFT, 'T'));
    first.close(100);
    second.close(100);
    EXPECT_EQ(first.sequence(), "CACGAG");
    EXPECT_EQ(second.sequence(), "TACGAA");
    EXPECT_EQ(original.sequence(), "ACGA");
    EXPECT_EQ(original.occurrences().size(), 2u);
    EXPECT_EQ(original.choices().size(), 2u);
}

TEST(ContextCursor, RejectsMatchesSpanningAdjacentRecords) {
    auto graph = coordinate_test_graph(3, {{"AACG", "F", "one"}, {"ACGTT", "F", "two"}});
    ContextCursor spanning(*graph, "AACGTT");
    EXPECT_TRUE(spanning.occurrences().empty());
    ContextCursor common(*graph, "ACG");
    common.close(100);
    EXPECT_EQ(common.sequence(), "ACG");
    EXPECT_EQ(common.occurrences().size(), 2u);
    EXPECT_EQ(common.choices().size(), 2u);
}

TEST(ContextCursor, RetainsAllRepeatOccurrencesAndRecordEnds) {
    auto graph = coordinate_test_graph(3, {{"AAAAAAAA", "F", "repeat"}});
    ContextCursor cursor(*graph, "AAA");
    EXPECT_EQ(cursor.occurrences().size(), 6u);
    cursor.close(100);
    EXPECT_EQ(cursor.sequence(), "AAA");
    EXPECT_TRUE(cursor.extend(ContextCursor::Side::LEFT, 'A'));
    EXPECT_EQ(cursor.occurrences().size(), 5u);
    EXPECT_EQ(cursor.sequence(), "AAAA");
    cursor.close(100);
    EXPECT_EQ(cursor.sequence(), "AAAA");
    EXPECT_THROW((ContextCursor(*graph, "AAA", {}, 5)), std::length_error);
}

TEST(ContextCursor, FiltersColumnsAndBoundsClosure) {
    auto graph = coordinate_test_graph(3, {{"TACGAG", "A", "one"}, {"CACGAA", "B", "two"}});
    ContextCursor cursor(*graph, "acg", {"A"});
    EXPECT_EQ(cursor.close(1), (Closure{1, 1, Stop::RADIUS_LIMIT, Stop::RADIUS_LIMIT}));
    EXPECT_EQ(cursor.sequence(), "TACGA");
    EXPECT_EQ(cursor.close(1), (Closure{0, 1, Stop::RECORD_END, Stop::RADIUS_LIMIT}));
    EXPECT_EQ(cursor.sequence(), "TACGAG");
    EXPECT_THROW((ContextCursor(*graph, "ACN")), std::invalid_argument);
    EXPECT_THROW((ContextCursor(*graph, "AC")), std::invalid_argument);
    EXPECT_THROW((ContextCursor(*graph, "ACG", {"absent"})), std::exception);
}

TEST(ContextCursor, ReconstructsNonAnchorCoordinateRows) {
    std::vector<test::CoordinateRecord> records{{"GACGATT", "F", "one"},
                                                {"TACGACC", "F", "two"}};
    auto plain = coordinate_test_graph(3, records);
    auto diffs = coordinate_test_graph(3, records, true);
    const auto &matrix = dynamic_cast<const annot::RowDiffCoordAnnotator::binary_matrix_type &>(
        diffs->get_annotator().get_matrix());
    EXPECT_LT(matrix.anchor().num_set_bits(), matrix.num_rows());
    for (const auto &seed : {"GAC", "CGA", "ACG", "TACGA"}) {
        ContextCursor a(*plain, seed), b(*diffs, seed);
        EXPECT_EQ(a.close(100), b.close(100));
        EXPECT_EQ(a.sequence(), b.sequence());
        EXPECT_EQ(a.occurrences().size(), b.occurrences().size());
        auto x = a.choices(), y = b.choices();
        ASSERT_EQ(x.size(), y.size());
        for (size_t i = 0; i < x.size(); ++i) {
            EXPECT_EQ(x[i].left, y[i].left);
            EXPECT_EQ(x[i].right, y[i].right);
            EXPECT_EQ(x[i].count, y[i].count);
        }
    }
}

TEST(ContextCursor, BRWTRowsMayHaveUnsortedColumns) {
    std::vector<test::CoordinateRecord> records{{"CACGAG", "F0", "one"},
                                                {"TACGAA", "F1", "two"},
                                                {"GACGAT", "F2", "three"}};
    for (bool row_diff : {false, true}) {
        auto graph = coordinate_test_graph(3, records, row_diff, true);
        ContextCursor cursor(*graph, "ACGA");
        ASSERT_EQ(cursor.occurrences().size(), 3u);
        EXPECT_EQ(cursor.close(100), (Closure{0, 0, Stop::BRANCH, Stop::BRANCH}));
        EXPECT_EQ(cursor.choices().size(), 3u);
        EXPECT_TRUE(cursor.extend(ContextCursor::Side::LEFT, 'C'));
        cursor.close(100);
        EXPECT_EQ(cursor.sequence(), "CACGAG");
    }
}

TEST(ContextCursor, StopsAtUnindexedSourceContext) {
    auto graph = coordinate_test_graph(3, {{"AACGNTTTA", "F", "gapped"}});
    ContextCursor cursor(*graph, "ACG");
    cursor.close(100);
    EXPECT_EQ(cursor.sequence(), "AACG");
    auto choices = cursor.choices();
    ASSERT_EQ(choices.size(), 1u);
    EXPECT_EQ(choices[0].left, '\0');
    EXPECT_EQ(choices[0].right, '?');
}

TEST(ContextCursor, RequiresCoordinatesAndRecordBoundaries) {
    auto graph = std::make_shared<DBGSuccinct>(3);
    graph->add_sequence("ACGTT");
    graph->mask_dummy_kmers(1, false);
    auto annotation = std::make_unique<annot::ColumnCompressed<>>(graph->max_index());
    AnnotatedDBG unindexed(graph, std::move(annotation));
    EXPECT_THROW((ContextCursor(unindexed, "ACG")), std::invalid_argument);
}

TEST(ContextCursor, LongerSeedResolvesTerminalKmerBranch) {
    auto graph = coordinate_test_graph(3, {{"GACGATT", "F", "one"}, {"TACGACC", "F", "two"}});
    ContextCursor short_seed(*graph, "CGA");
    short_seed.close(100);
    EXPECT_EQ(short_seed.sequence(), "ACGA");
    ContextCursor long_seed(*graph, "GACGA");
    long_seed.close(100);
    EXPECT_EQ(long_seed.sequence(), "GACGATT");
}

TEST(ContextCursor, ReportsEmptySetsLimitsAndBoundaryStops) {
    auto graph = coordinate_test_graph(3, {{"TACGAG", "F", "one"}});
    ContextCursor absent(*graph, "GGG");
    EXPECT_EQ(absent.close(0), (Closure{0, 0, Stop::NO_OCCURRENCES, Stop::NO_OCCURRENCES}));
    EXPECT_TRUE(absent.choices().empty());
    EXPECT_FALSE(absent.extend(ContextCursor::Side::LEFT, 'A'));
    EXPECT_THROW((ContextCursor(*graph, "GGG", {}, 0)), std::invalid_argument);

    ContextCursor cursor(*graph, "ACG", {}, 1);
    EXPECT_EQ(cursor.close(0), (Closure{0, 0, Stop::RADIUS_LIMIT, Stop::RADIUS_LIMIT}));
    EXPECT_EQ(cursor.sequence(), "ACG");
    EXPECT_EQ(cursor.close(1), (Closure{1, 1, Stop::RADIUS_LIMIT, Stop::RADIUS_LIMIT}));
    EXPECT_EQ(cursor.close(1), (Closure{0, 1, Stop::RECORD_END, Stop::RADIUS_LIMIT}));
    EXPECT_EQ(cursor.close(1), (Closure{0, 0, Stop::RECORD_END, Stop::RECORD_END}));

    auto gapped = coordinate_test_graph(3, {{"AACGNTTTA", "F", "gap"}});
    EXPECT_EQ(ContextCursor(*gapped, "ACG").close(10),
              (Closure{1, 0, Stop::RECORD_END, Stop::UNINDEXED}));
    EXPECT_EQ(ContextCursor(*gapped, "TTT").close(10),
              (Closure{0, 1, Stop::UNINDEXED, Stop::RECORD_END}));
}

TEST(ContextCursor, DoesNotDropTruncatedOccurrencesDuringClosure) {
    auto graph = coordinate_test_graph(3, {{"ACGA", "F", "end"},
                                          {"TACGAN", "F", "gap"},
                                          {"CACGAG", "F", "long"}});
    ContextCursor cursor(*graph, "ACGA");
    EXPECT_EQ(cursor.close(100), (Closure{0, 0, Stop::RECORD_END, Stop::RECORD_END}));
    EXPECT_EQ(cursor.occurrences().size(), 3u);
    EXPECT_TRUE(cursor.extend(ContextCursor::Side::LEFT, 'T'));
    EXPECT_EQ(cursor.close(100), (Closure{0, 0, Stop::RECORD_END, Stop::UNINDEXED}));
    EXPECT_EQ(cursor.occurrences().size(), 1u);
}

TEST(ContextCursor, MatchesLiteralRecordsThroughRepeatedChoices) {
    // Independent model: look directly at source strings, with no graph or
    // annotation reads. Include repeats, gaps, multiple records and columns.
    std::mt19937 random(670);
    for (size_t trial = 0; trial < 24; ++trial) {
        std::vector<test::CoordinateRecord> records{
            {"CACGAG", "F0", "paired"}, {"TACGAA", "F0", "paired"},
            {"AAAANAAAAAAA", "F1", "repeat"}
        };
        for (size_t r = 0; r < 3; ++r) {
            std::string sequence(12 + random() % 12, 'A');
            for (char &c : sequence)
                c = "ACGTN"[random() % 5];
            records.push_back({sequence, "F1", "random"});
        }
        for (bool row_diff : {false, true}) {
            for (bool brwt : {false, true}) {
                auto graph = coordinate_test_graph(3, records, row_diff, brwt);
                for (const std::string seed : {"ACGA", "AAA", "ACG", "TTT", "GGG"}) {
                    SCOPED_TRACE(::testing::Message() << trial << ' ' << row_diff << ' ' << brwt << ' ' << seed);
                    ContextCursor cursor(*graph, seed);
                    std::string sequence = seed;
                    // (source record index, zero-based start of current sequence)
                    std::vector<std::pair<size_t, size_t>> hits;
                    for (size_t r = 0; r < records.size(); ++r)
                        for (size_t p = 0; p + seed.size() <= records[r].sequence.size(); ++p)
                            if (records[r].sequence.compare(p, seed.size(), seed) == 0)
                                hits.emplace_back(r, p);

                    auto boundary = [&](const auto &hit, size_t side) -> char {
                        const auto &[r, p] = hit;
                        const auto &source = records[r].sequence;
                        if ((!side && !p) || (side && p + sequence.size() == source.size()))
                            return '\0';
                        char c = source[side ? p + sequence.size() : p - 1];
                        return c == 'N' ? '?' : c;
                    };
                    auto extend = [&](size_t side, char base) {
                        auto selected = hits;
                        std::erase_if(selected, [&](const auto &hit) { return boundary(hit, side) != base; });
                        if (selected.empty())
                            return false;
                        hits = std::move(selected);
                        if (side) {
                            sequence.push_back(base);
                        } else {
                            sequence.insert(sequence.begin(), base);
                            for (auto &[r, p] : hits) --p;
                        }
                        return true;
                    };
                    auto close = [&](size_t radius) {
                        std::array<size_t, 2> added{};
                        std::array<Stop, 2> stops{};
                        for (size_t side = 0; side < 2; ++side) {
                            while (true) {
                                if (hits.empty()) { stops[side] = Stop::NO_OCCURRENCES; break; }
                                if (added[side] == radius) { stops[side] = Stop::RADIUS_LIMIT; break; }
                                std::string bases;
                                for (const auto &hit : hits) bases.push_back(boundary(hit, side));
                                if (bases.find('\0') != std::string::npos) { stops[side] = Stop::RECORD_END; break; }
                                if (bases.find('?') != std::string::npos) { stops[side] = Stop::UNINDEXED; break; }
                                if (bases.find_first_not_of(bases.front()) != std::string::npos) {
                                    stops[side] = Stop::BRANCH;
                                    break;
                                }
                                extend(side, bases.front());
                                ++added[side];
                            }
                        }
                        return Closure{added[0], added[1], stops[0], stops[1]};
                    };
                    auto check = [&] {
                        EXPECT_EQ(cursor.sequence(), sequence);
                        using Hit = std::tuple<uint64_t, size_t, uint64_t>;
                        std::vector<Hit> expected, actual;
                        for (const auto &[r, p] : hits)
                            expected.emplace_back(r < 2 ? 0 : 1, r < 2 ? r : r - 2, p);
                        for (const auto &hit : cursor.occurrences())
                            actual.emplace_back(hit.column, hit.record, hit.offset);
                        std::ranges::sort(expected);
                        std::ranges::sort(actual);
                        EXPECT_EQ(actual, expected);
                        std::map<std::pair<char, char>, size_t> choices, observed;
                        for (const auto &hit : hits) ++choices[{boundary(hit, 0), boundary(hit, 1)}];
                        for (const auto &c : cursor.choices()) observed[{c.left, c.right}] = c.count;
                        EXPECT_EQ(observed, choices);
                    };
                    check();
                    for (size_t round = 0; round < 6; ++round) {
                        size_t radius = random() % 8;
                        EXPECT_EQ(cursor.close(radius), close(radius));
                        check();
                        size_t side = random() % 2;
                        char base = "ACGT"[random() % 4];
                        auto branch = cursor;
                        bool extended = extend(side, base);
                        EXPECT_EQ(branch.extend(side ? ContextCursor::Side::RIGHT : ContextCursor::Side::LEFT, base), extended);
                        EXPECT_EQ(cursor.extend(side ? ContextCursor::Side::RIGHT : ContextCursor::Side::LEFT, base), extended);
                        EXPECT_EQ(cursor.sequence(), branch.sequence());
                        check();
                    }
                }
            }
        }
    }
}
} // namespace
} // namespace mtg::graph::traversal

#endif // _DNA_GRAPH
