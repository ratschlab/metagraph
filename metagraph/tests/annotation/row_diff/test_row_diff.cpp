#include <gtest/gtest.h>

// work around clang-related bug in GMock
namespace testing {
namespace internal {
using UInt64 = uint64_t;
using Int64 = int64_t;
using Int32 = int32_t;
} // namespace internal
} // namespace testing

#include <gmock/gmock.h>
#include <graph/representation/succinct/dbg_succinct.hpp>

#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "common/utils/file_utils.hpp"

namespace {
using namespace mtg;
using namespace testing;
using ::testing::_;

TEST(RowDiff, Empty) {
    annot::binmat::RowDiff rowdiff;
    EXPECT_EQ(0, rowdiff.diffs().size());
    EXPECT_EQ(0, rowdiff.boundary().size());
    EXPECT_EQ(0, rowdiff.terminal().size());
    EXPECT_EQ(0, rowdiff.num_relations());
}

TEST(RowDiff, Serialize) {
    Vector<uint64_t> diffs = { 1, 2, 3, 4 };
    sdsl::bit_vector boundary = { 0, 1, 0, 0, 1, 1 };
    sdsl::bit_vector terminal = { 0, 0, 0, 1 };
    sdsl::enc_vector<> ediffs(diffs);
    annot::binmat::RowDiff annot(2, 5, nullptr, ediffs, boundary, terminal);

    utils::TempFile tempfile;
    std::ofstream &out = tempfile.ofstream();
    annot.serialize(out);

    annot::binmat::RowDiff loaded;
    EXPECT_TRUE(loaded.load(tempfile.ifstream()));
    EXPECT_EQ(annot.diffs().size(), loaded.diffs().size());
    for (uint32_t i = 0; i < annot.diffs().size(); ++i) {
        EXPECT_EQ(annot.diffs()[i], loaded.diffs()[i]);
    }

    EXPECT_EQ(annot.boundary().size(), loaded.boundary().size());
    for (uint32_t i = 0; i < annot.boundary().size(); ++i) {
        EXPECT_EQ(annot.boundary()[i], loaded.boundary()[i]);
    }

    EXPECT_EQ(annot.terminal().size(), loaded.terminal().size());
    for (uint32_t i = 0; i < annot.terminal().size(); ++i) {
        EXPECT_EQ(annot.terminal()[i], loaded.terminal()[i]);
    }
    EXPECT_EQ(2, loaded.num_columns());
    EXPECT_EQ(5, loaded.num_relations());
    EXPECT_EQ(4, loaded.num_rows());
}

TEST(RowDiff, GetDiff) {
    Vector<uint64_t> diffs = { 1, 2, 3, 4 };
    sdsl::bit_vector boundary = { 1, 0, 1, 0, 0, 1, 0, 1 };
    sdsl::bit_vector terminal = { 0, 0, 0, 1 };
    annot::binmat::RowDiff annot(2, 5, nullptr, sdsl::enc_vector<>(diffs), boundary,
                                 terminal);

    EXPECT_TRUE(annot.get_diff(0).empty());
    ASSERT_THAT(annot.get_diff(1), ElementsAre(1));
    ASSERT_THAT(annot.get_diff(2), ElementsAre(2, 3));
    ASSERT_THAT(annot.get_diff(3), ElementsAre(4));
}

TEST(RowDiff, GetDiff2) {
    Vector<uint64_t> diffs = { 0, 1, 1, 1 };
    sdsl::bit_vector boundary = { 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1 };
    sdsl::bit_vector terminal = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0 };
    annot::binmat::RowDiff annot(2, 5, nullptr, sdsl::enc_vector<>(diffs), boundary, terminal);

    std::vector<uint32_t> empty_annot = { 0, 1, 2, 4, 5, 6, 7, 9 };
    for (const uint32_t v : empty_annot) {
        EXPECT_TRUE(annot.get_diff(v).empty());
    }
    ASSERT_THAT(annot.get_diff(3), ElementsAre(0));
    ASSERT_THAT(annot.get_diff(8), ElementsAre(1));
    ASSERT_THAT(annot.get_diff(10), ElementsAre(1));
    ASSERT_THAT(annot.get_diff(11), ElementsAre(1));
}

/**
 * Tests annotations on the graph in
 * https://docs.google.com/document/d/1e0MFgZRJfmDUSvmDPuC_lvnnWA0VKm5hPdzM8mdrHMM/edit#bookmark=id.ciri4266pkc4
 */
TEST(RowDiff, GetAnnotation) {
    graph::DBGSuccinct graph(4);
    graph.add_sequence("ACTAGCTAGCTAGCTAGCTAGC");
    graph.add_sequence("ACTCTAG");

    Vector<uint64_t> diffs = { 0, 1, 1, 1 };
    sdsl::bit_vector boundary = { 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1 };
    sdsl::bit_vector terminal = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0 };
    annot::binmat::RowDiff annot(2, 5, &graph, sdsl::enc_vector<>(diffs), boundary, terminal);

    EXPECT_EQ("CTAG", graph.get_node_sequence(4));
    ASSERT_THAT(annot.get_row(3), ElementsAre(0, 1));

    EXPECT_EQ("AGCT", graph.get_node_sequence(6));
    ASSERT_THAT(annot.get_row(5), ElementsAre(1));

    EXPECT_EQ("CTCT", graph.get_node_sequence(7));
    ASSERT_THAT(annot.get_row(6), ElementsAre(0));

    EXPECT_EQ("TAGC", graph.get_node_sequence(8));
    ASSERT_THAT(annot.get_row(7), ElementsAre(1));

    EXPECT_EQ("ACTA", graph.get_node_sequence(9));
    ASSERT_THAT(annot.get_row(8), ElementsAre(1));

    EXPECT_EQ("ACTC", graph.get_node_sequence(10));
    ASSERT_THAT(annot.get_row(9), ElementsAre(0));

    EXPECT_EQ("GCTA", graph.get_node_sequence(11));
    ASSERT_THAT(annot.get_row(10), ElementsAre(1));

    EXPECT_EQ("TCTA", graph.get_node_sequence(12));
    ASSERT_THAT(annot.get_row(11), ElementsAre(0));
}

/**
 * Tests annotations on the graph in
 * https://docs.google.com/document/d/1e0MFgZRJfmDUSvmDPuC_lvnnWA0VKm5hPdzM8mdrHMM/edit#bookmark=id.ciri4266pkc4
 * after having removed dummy nodes.
 */
TEST(RowDiff, GetAnnotationMasked) {
    graph::DBGSuccinct graph(4);
    graph.add_sequence("ACTAGCTAGCTAGCTAGCTAGC");
    graph.add_sequence("ACTCTAG");
    graph.mask_dummy_kmers(1, false);

    Vector<uint64_t> diffs = { 0, 1, 1, 1 };
    sdsl::bit_vector boundary = { 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1 };
    sdsl::bit_vector terminal = { 0, 0, 0, 0, 1, 0, 1, 0 };
    annot::binmat::RowDiff annot(2, 9, &graph, sdsl::enc_vector<>(diffs), boundary, terminal);

    EXPECT_EQ("CTAG", graph.get_node_sequence(1));
    ASSERT_THAT(annot.get_row(0), ElementsAre(0, 1));

    EXPECT_EQ("AGCT", graph.get_node_sequence(2));
    ASSERT_THAT(annot.get_row(1), ElementsAre(1));

    EXPECT_EQ("CTCT", graph.get_node_sequence(3));
    ASSERT_THAT(annot.get_row(2), ElementsAre(0));

    EXPECT_EQ("TAGC", graph.get_node_sequence(4));
    ASSERT_THAT(annot.get_row(3), ElementsAre(1));

    EXPECT_EQ("ACTA", graph.get_node_sequence(5));
    ASSERT_THAT(annot.get_row(4), ElementsAre(1));

    EXPECT_EQ("ACTC", graph.get_node_sequence(6));
    ASSERT_THAT(annot.get_row(5), ElementsAre(0));

    EXPECT_EQ("GCTA", graph.get_node_sequence(7));
    ASSERT_THAT(annot.get_row(6), ElementsAre(1));

    EXPECT_EQ("TCTA", graph.get_node_sequence(8));
    ASSERT_THAT(annot.get_row(7), ElementsAre(0));
}

/**
 *  Tests that annotations for the graph in
 *  https://docs.google.com/document/d/1siWApHWBDtiYCsetb6vHPuT7WwBkVFJp5fgti_AJK-s/edit
 *  are correctly retrieved.
 */
TEST(RowDiff, GetAnnotationBifurcation) {
    graph::DBGSuccinct graph(4);
    graph.add_sequence("TACTAGCTAGCTAGCTAGCTAGC");
    graph.add_sequence("ACTCTAGCTAT");

    Vector<uint64_t> diffs = { 1, 0, 1, 0, 0, 1 };
    sdsl::bit_vector boundary
            = { 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };
    sdsl::bit_vector terminal = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0 };
    annot::binmat::RowDiff annot(2, 9, &graph, sdsl::enc_vector<>(diffs), boundary, terminal);

    EXPECT_EQ("CTAG", graph.get_node_sequence(4));
    ASSERT_THAT(annot.get_row(3), ElementsAre(0, 1));

    EXPECT_EQ("CTAT", graph.get_node_sequence(5));
    ASSERT_THAT(annot.get_row(4), ElementsAre(1));

    EXPECT_EQ("TACT", graph.get_node_sequence(6));
    ASSERT_THAT(annot.get_row(5), ElementsAre(0));

    EXPECT_EQ("AGCT", graph.get_node_sequence(7));
    ASSERT_THAT(annot.get_row(6), ElementsAre(0, 1));

    EXPECT_EQ("CTCT", graph.get_node_sequence(8));
    ASSERT_THAT(annot.get_row(7), ElementsAre(1));

    EXPECT_EQ("TAGC", graph.get_node_sequence(9));
    ASSERT_THAT(annot.get_row(8), ElementsAre(0, 1));

    EXPECT_EQ("ACTA", graph.get_node_sequence(12));
    ASSERT_THAT(annot.get_row(11), ElementsAre(0));

    EXPECT_EQ("ACTC", graph.get_node_sequence(13));
    ASSERT_THAT(annot.get_row(12), ElementsAre(1));

    EXPECT_EQ("GCTA", graph.get_node_sequence(14));
    ASSERT_THAT(annot.get_row(13), ElementsAre(0, 1));

    EXPECT_EQ("TCTA", graph.get_node_sequence(15));
    ASSERT_THAT(annot.get_row(14), ElementsAre(1));
}

TEST(RowDiff, GetAnnotationBifurcationMasked) {
    graph::DBGSuccinct graph(4);
    graph.add_sequence("TACTAGCTAGCTAGCTAGCTAGC");
    graph.add_sequence("ACTCTAGCTAT");
    graph.mask_dummy_kmers(1, false);

    Vector<uint64_t> diffs = { 1, 0, 1, 0, 0, 1 };
    sdsl::bit_vector boundary = { 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1 };
    sdsl::bit_vector terminal = { 0, 1, 0, 0, 0, 0, 1, 0, 1, 0 };
    annot::binmat::RowDiff annot(2, 9, &graph, sdsl::enc_vector<>(diffs), boundary, terminal);

    EXPECT_EQ("CTAG", graph.get_node_sequence(1));
    ASSERT_THAT(annot.get_row(0), ElementsAre(0, 1));

    EXPECT_EQ("CTAT", graph.get_node_sequence(2));
    ASSERT_THAT(annot.get_row(1), ElementsAre(1));

    EXPECT_EQ("TACT", graph.get_node_sequence(3));
    ASSERT_THAT(annot.get_row(2), ElementsAre(0));

    EXPECT_EQ("AGCT", graph.get_node_sequence(4));
    ASSERT_THAT(annot.get_row(3), ElementsAre(0, 1));

    EXPECT_EQ("CTCT", graph.get_node_sequence(5));
    ASSERT_THAT(annot.get_row(4), ElementsAre(1));

    EXPECT_EQ("TAGC", graph.get_node_sequence(6));
    ASSERT_THAT(annot.get_row(5), ElementsAre(0, 1));

    EXPECT_EQ("ACTA", graph.get_node_sequence(7));
    ASSERT_THAT(annot.get_row(6), ElementsAre(0));

    EXPECT_EQ("ACTC", graph.get_node_sequence(8));
    ASSERT_THAT(annot.get_row(7), ElementsAre(1));

    EXPECT_EQ("GCTA", graph.get_node_sequence(9));
    ASSERT_THAT(annot.get_row(8), ElementsAre(0, 1));

    EXPECT_EQ("TCTA", graph.get_node_sequence(10));
    ASSERT_THAT(annot.get_row(9), ElementsAre(1));
}

} // namespace
