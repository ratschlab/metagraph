#include "gtest/gtest.h"
#include "tests/annotation/coordinate_test_graph.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/traversal/context_cursor.hpp"

#if defined(_DNA_GRAPH)
namespace mtg::graph::traversal {
namespace {
using test::coordinate_test_graph;

TEST(ContextCursor, PreservesJointFlanksAfterSelection) {
    auto graph = coordinate_test_graph(3, {{"CACGAG", "F", "first"},
                                            {"TACGAA", "F", "second"}});
    ContextCursor cursor(*graph, "ACGA");
    EXPECT_EQ(cursor.close(100), (std::pair<size_t, size_t>{0, 0}));
    auto choices = cursor.choices();
    ASSERT_EQ(choices.size(), 2u);
    EXPECT_EQ(choices[0].left, 'C'); EXPECT_EQ(choices[0].right, 'G');
    EXPECT_EQ(choices[1].left, 'T'); EXPECT_EQ(choices[1].right, 'A');
    EXPECT_TRUE(cursor.extend(ContextCursor::Side::LEFT, 'C'));
    EXPECT_EQ(cursor.close(100), (std::pair<size_t, size_t>{0, 1}));
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
    EXPECT_EQ(cursor.close(1), (std::pair<size_t, size_t>{1, 1}));
    EXPECT_EQ(cursor.sequence(), "TACGA");
    EXPECT_EQ(cursor.close(1), (std::pair<size_t, size_t>{0, 1}));
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
        EXPECT_EQ(cursor.close(100), (std::pair<size_t, size_t>{0, 0}));
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
} // namespace
} // namespace mtg::graph::traversal

#endif // _DNA_GRAPH
