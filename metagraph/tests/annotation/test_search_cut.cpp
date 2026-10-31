#include "gtest/gtest.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "../test_helpers.hpp"
#include "../graph/all/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "graph/annotated_dbg.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;
using namespace mtg::graph;
using mtg::annot::matrix::CSRMatrix;

// The C5 order and the cut of a search's top N (SPEC-multigraph-envelope.md section 3)

TEST(SearchCut, LowerMedian) {
    // odd n: the middle element
    EXPECT_EQ(7u, lower_median({ 100, 7, 3 }));
    EXPECT_EQ(5u, lower_median({ 5 }));
    // even n: the lower of the two middle elements, index floor((n - 1) / 2)
    EXPECT_EQ(3u, lower_median({ 100, 7, 3, 1 }));
    EXPECT_EQ(1u, lower_median({ 1, 100 }));
    // zeros excluded: positions where the label is absent are not abundances
    EXPECT_EQ(100u, lower_median({ 0, 0, 0, 100 }));
    EXPECT_EQ(1u, lower_median({ 0, 1, 0, 100, 0 }));
    EXPECT_EQ(7u, lower_median({ 0, 100, 100, 3, 0, 0, 7, 7, 7, 7, 7 }));
    // repeated values, one per position
    EXPECT_EQ(2u, lower_median({ 2, 2, 2, 9 }));
    // no value left
    EXPECT_EQ(0u, lower_median({}));
    EXPECT_EQ(0u, lower_median({ 0, 0 }));
}

TEST(SearchCut, RowOrder) {
    // score descending first
    EXPECT_TRUE(search_row_before({ 5, 0, "z" }, { 4, 100, "a" }));
    EXPECT_FALSE(search_row_before({ 4, 100, "a" }, { 5, 0, "z" }));
    // then the median descending
    EXPECT_TRUE(search_row_before({ 5, 100, "z" }, { 5, 1, "a" }));
    EXPECT_FALSE(search_row_before({ 5, 1, "a" }, { 5, 100, "z" }));
    // then the label ascending, by its bytes
    EXPECT_TRUE(search_row_before({ 5, 7, "a" }, { 5, 7, "b" }));
    EXPECT_TRUE(search_row_before({ 5, 7, "Z" }, { 5, 7, "a" }));
    EXPECT_TRUE(search_row_before({ 5, 7, "a" }, { 5, 7, "ab" }));
    EXPECT_TRUE(search_row_before({ 5, 7, "NC_1" }, { 5, 7, "NZ_0" }));
    // a strict weak order: no row before itself
    EXPECT_FALSE(search_row_before({ 5, 7, "a" }, { 5, 7, "a" }));
}

// search_top_n over |labels| with |scores| (and |medians|, if any), the labels decoded counted
struct TopN {
    std::vector<std::pair<size_t, std::string>> rows;
    SearchCut cut;
    std::vector<size_t> decoded;
};

TopN top_n(const std::vector<std::string> &labels,
           const std::vector<size_t> &scores,
           size_t n,
           const std::vector<uint64_t> &medians = {}) {
    TopN out;
    std::function<std::vector<uint64_t>(const std::vector<size_t> &)> medians_of;
    if (medians.size()) {
        medians_of = [&](const std::vector<size_t> &indexes) {
            std::vector<uint64_t> values;
            for (size_t i : indexes) {
                values.push_back(medians.at(i));
            }
            return values;
        };
    }
    out.rows = search_top_n(scores, n, [&](size_t i) {
        out.decoded.push_back(i);
        return labels.at(i);
    }, medians_of, &out.cut);
    std::sort(out.decoded.begin(), out.decoded.end());
    return out;
}

std::vector<std::string> labels_of(const TopN &top) {
    std::vector<std::string> labels;
    for (const auto &[i, label] : top.rows) {
        labels.push_back(label);
    }
    return labels;
}

TEST(SearchCut, Complete) {
    // fewer matches than N: every one returned, in the C5 order, no N-th row
    auto top = top_n({ "c", "a", "b" }, { 3, 5, 3 }, 10);
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c" }), labels_of(top));
    EXPECT_EQ(3u, top.cut.returned);
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.tied);
    EXPECT_EQ(0u, top.cut.tied_returned);
    EXPECT_FALSE(top.cut.abundances);

    // exactly N matches: nothing is cut
    top = top_n({ "c", "a", "b" }, { 3, 5, 3 }, 3);
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c" }), labels_of(top));
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.tied);
}

TEST(SearchCut, Empty) {
    auto top = top_n({}, {}, 10);
    EXPECT_TRUE(top.rows.empty());
    EXPECT_EQ(0u, top.cut.returned);
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.nth_score);
    EXPECT_EQ(0u, top.cut.tied);
    EXPECT_EQ(0u, top.cut.tied_returned);
    EXPECT_TRUE(top.decoded.empty());
}

TEST(SearchCut, TiesStraddleTheCut) {
    // scores 9, then four labels at 5, then 1: N = 3 cuts the tie group at 5, which
    // competes for two places by label
    std::vector<std::string> labels { "e", "low", "d", "top", "b", "c" };
    std::vector<size_t> scores {       5,   1,    5,    9,    5,   5  };
    auto top = top_n(labels, scores, 3);
    EXPECT_EQ((std::vector<std::string>{ "top", "b", "c" }), labels_of(top));
    EXPECT_EQ(3u, top.cut.returned);
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(5u, top.cut.nth_score);
    EXPECT_EQ(4u, top.cut.tied);
    EXPECT_EQ(2u, top.cut.tied_returned);
    EXPECT_FALSE(top.cut.abundances);
    // the rows below the N-th score are never decoded
    EXPECT_EQ((std::vector<size_t>{ 0, 2, 3, 4, 5 }), top.decoded);

    // the returned rows point at their candidates
    for (const auto &[i, label] : top.rows) {
        EXPECT_EQ(labels[i], label);
    }

    // a cut between two scores: the N-th score has one label, returned
    top = top_n(labels, scores, 5);
    EXPECT_EQ((std::vector<std::string>{ "top", "b", "c", "d", "e" }), labels_of(top));
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(5u, top.cut.nth_score);
    EXPECT_EQ(4u, top.cut.tied);
    EXPECT_EQ(4u, top.cut.tied_returned);
    EXPECT_EQ((std::vector<size_t>{ 0, 2, 3, 4, 5 }), top.decoded);

    // N = 1: the single top row, no tie
    top = top_n(labels, scores, 1);
    EXPECT_EQ((std::vector<std::string>{ "top" }), labels_of(top));
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(9u, top.cut.nth_score);
    EXPECT_EQ(1u, top.cut.tied);
    EXPECT_EQ(1u, top.cut.tied_returned);
    EXPECT_EQ((std::vector<size_t>{ 3 }), top.decoded);
}

TEST(SearchCut, MediansOrderTheTies) {
    // a counts annotation: the tie group at the cut is ordered by median, then by label
    std::vector<std::string> labels   { "a", "b", "c", "d", "top" };
    std::vector<size_t> scores        {  12,  12,  12,  12,  20   };
    std::vector<uint64_t> medians     {   1,   7,   7, 100,  1    };
    auto top = top_n(labels, scores, 3, medians);
    EXPECT_EQ((std::vector<std::string>{ "top", "d", "b" }), labels_of(top));
    EXPECT_EQ(3u, top.cut.returned);
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(12u, top.cut.nth_score);
    EXPECT_EQ(4u, top.cut.tied);
    EXPECT_EQ(2u, top.cut.tied_returned);
    EXPECT_TRUE(top.cut.abundances);
    EXPECT_EQ(7u, top.cut.nth_median);
    EXPECT_EQ(2u, top.cut.tied_same_median);
    // the tie member of a lower median than the N-th row's is not decoded
    EXPECT_EQ((std::vector<size_t>{ 1, 2, 3, 4 }), top.decoded);

    // N = 1 of two labels of equal k-mer counts: the higher median wins over the label
    top = top_n({ "SRR0999", "SRR1001" }, { 12, 12 }, 1, { 1, 7 });
    EXPECT_EQ((std::vector<std::string>{ "SRR1001" }), labels_of(top));
    EXPECT_EQ(12u, top.cut.nth_score);
    EXPECT_EQ(2u, top.cut.tied);
    EXPECT_EQ(1u, top.cut.tied_returned);
    EXPECT_EQ(7u, top.cut.nth_median);
    EXPECT_EQ(1u, top.cut.tied_same_median);

    // returned rows of one score are ordered by median also without a cut
    top = top_n(labels, scores, 10, medians);
    EXPECT_EQ((std::vector<std::string>{ "top", "d", "b", "c", "a" }), labels_of(top));
    EXPECT_TRUE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.nth_median);
    EXPECT_EQ(0u, top.cut.tied_same_median);
}

TEST(SearchCut, NoneRequested) {
    // N = 0: nothing returned, the matches cut, no N-th row
    auto top = top_n({ "a", "b" }, { 3, 2 }, 0);
    EXPECT_TRUE(top.rows.empty());
    EXPECT_EQ(0u, top.cut.returned);
    EXPECT_FALSE(top.cut.complete);
    EXPECT_EQ(0u, top.cut.tied);
    EXPECT_TRUE(top.decoded.empty());
}


// An index whose annotation holds k-mer counts: one sequence, every k-mer of it annotated
// with |abundances[j]| in column j
std::unique_ptr<AnnotatedDBG> counts_index(const std::string &sequence, size_t k,
                                           const std::vector<std::string> &labels,
                                           const std::vector<uint64_t> &abundances) {
    auto graph = build_graph_batch<DBGHashOrdered>(k, { sequence });
    Vector<CSRMatrix::RowValues> rows(graph->max_index());
    graph->map_to_nodes(sequence, [&](DeBruijnGraph::node_index node) {
        auto &row = rows[AnnotatedDBG::graph_to_anno_index(node)];
        row.clear();
        for (size_t j = 0; j < labels.size(); ++j) {
            if (abundances[j])
                row.emplace_back(j, abundances[j]);
        }
    });
    annot::LabelEncoder<std::string> encoder;
    for (const auto &label : labels) {
        encoder.insert_and_encode(label);
    }
    return std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<annot::IntRowAnnotator>(
            std::make_unique<CSRMatrix>(std::move(rows), labels.size()), encoder));
}

TEST(SearchCut, CountsTieBrokenByMedianAbundance) {
    const size_t k = 5;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAAC";
    const size_t num_kmers = sequence.size() - k + 1;
    // two labels on the same k-mers, abundance 1 and 100: equal k-mer counts. The label
    // order alone would return "a_low"; the median abundance returns "b_high".
    auto index = counts_index(sequence, k, { "a_low", "b_high" }, { 1, 100 });
    ASSERT_TRUE(AnnotatedDBG::holds_kmer_abundances(index->get_annotator()));
    ASSERT_TRUE(index->has_kmer_abundances());

    auto check_cut = [&](const SearchCut &cut) {
        EXPECT_EQ(1u, cut.returned);
        EXPECT_FALSE(cut.complete);
        EXPECT_EQ(num_kmers, cut.nth_score);
        EXPECT_EQ(2u, cut.tied);
        EXPECT_EQ(1u, cut.tied_returned);
        EXPECT_TRUE(cut.abundances);
        EXPECT_EQ(100u, cut.nth_median);
        EXPECT_EQ(1u, cut.tied_same_median);
    };

    SearchCut cut;
    auto matches = index->get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("b_high", matches[0].first);
    EXPECT_EQ(num_kmers, matches[0].second);
    check_cut(cut);

    auto counts = index->get_kmer_counts(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, counts.size());
    EXPECT_EQ("b_high", std::get<0>(counts[0]));
    EXPECT_EQ(num_kmers, std::get<1>(counts[0]));
    EXPECT_EQ(std::vector<size_t>(num_kmers, 100), std::get<2>(counts[0]));
    check_cut(cut);

    auto signatures = index->get_top_label_signatures(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, signatures.size());
    EXPECT_EQ("b_high", std::get<0>(signatures[0]));
    check_cut(cut);

    // without a cut both are returned, the higher median first
    matches = index->get_top_labels(sequence, 2, 0.0, 0.0, false, &cut);
    ASSERT_EQ(2u, matches.size());
    EXPECT_EQ("b_high", matches[0].first);
    EXPECT_EQ("a_low", matches[1].first);
    EXPECT_TRUE(cut.complete);
    EXPECT_EQ(2u, cut.returned);
    EXPECT_EQ(0u, cut.tied);

    // the abundance sums differ: no tie
    matches = index->get_top_labels(sequence, 1, 0.0, 0.0, true, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("b_high", matches[0].first);
    EXPECT_EQ(100 * num_kmers, matches[0].second);
    EXPECT_EQ(100 * num_kmers, cut.nth_score);
    EXPECT_EQ(1u, cut.tied);
    EXPECT_EQ(1u, cut.tied_same_median);

    // no k-mer of the query in the graph: nothing returned, complete
    matches = index->get_top_labels("TTTTTTTTTT", 1, 0.0, 0.0, false, &cut);
    EXPECT_TRUE(matches.empty());
    EXPECT_EQ(0u, cut.returned);
    EXPECT_TRUE(cut.complete);
    EXPECT_TRUE(cut.abundances);
}

TEST(SearchCut, BinaryTieOrderedByLabel) {
    // a binary annotation: ties by label (byte order), not by column
    const size_t k = 5;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAAC";
    auto index = build_anno_graph<DBGHashOrdered, annot::ColumnCompressed<>>(
            k, { sequence, sequence, sequence }, { "b", "c", "a" });
    ASSERT_FALSE(index->has_kmer_abundances());

    SearchCut cut;
    auto matches = index->get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("a", matches[0].first);
    EXPECT_EQ(1u, cut.returned);
    EXPECT_FALSE(cut.complete);
    EXPECT_EQ(3u, cut.tied);
    EXPECT_EQ(1u, cut.tied_returned);
    EXPECT_FALSE(cut.abundances);

    matches = index->get_top_labels(sequence, 3, 0.0, 0.0, false, &cut);
    ASSERT_EQ(3u, matches.size());
    EXPECT_EQ("a", matches[0].first);
    EXPECT_EQ("b", matches[1].first);
    EXPECT_EQ("c", matches[2].first);
    EXPECT_TRUE(cut.complete);
}

TEST(SearchCut, CoordinateRecordsTiedByHeader) {
    // a coordinate annotation with its record mapping (.seqs): two records of one column,
    // the same sequence, so every k-mer matches both. The records are ranked by header; the
    // coordinates are positions, not abundances, and order no tie.
    const size_t k = 11;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAACGGATCCATGCAAGT";
    const uint64_t n = sequence.size() - k + 1;
    auto coords = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, { sequence, sequence }, { "F", "F" }, DeBruijnGraph::BASIC, true, { 0, n });
    ASSERT_FALSE(AnnotatedDBG::holds_kmer_abundances(coords->get_annotator()));

    std::filesystem::path dir = utils::create_temp_dir("", "test_search_cut");
    const std::string prefix = dir / "coords";
    coords->get_annotator().serialize(prefix);
    auto annotator = std::make_unique<annot::ColumnCoordAnnotator>();
    ASSERT_TRUE(annotator->load(prefix));
    utils::remove_temp_dir(dir);

    // the first record's header sorts after the second's
    auto mapping = std::make_shared<annot::CoordToHeader>(
            std::vector<std::vector<std::string>>{ { "z_record", "a_record" } },
            std::vector<std::vector<uint64_t>>{ { n, n } });
    auto graph = std::dynamic_pointer_cast<DeBruijnGraph>(
            std::const_pointer_cast<SequenceGraph>(coords->get_graph_ptr()));
    ASSERT_TRUE(graph);
    AnnotatedDBG index(graph, std::move(annotator), false, mapping);
    ASSERT_FALSE(index.has_kmer_abundances());

    auto check_cut = [&](const SearchCut &cut) {
        EXPECT_EQ(1u, cut.returned);
        EXPECT_FALSE(cut.complete);
        EXPECT_EQ(2u, cut.tied);
        EXPECT_EQ(1u, cut.tied_returned);
        EXPECT_FALSE(cut.abundances);
    };

    SearchCut cut;
    auto records = index.get_kmer_coordinates(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, records.size());
    EXPECT_EQ("a_record", std::get<0>(records[0]));
    EXPECT_EQ(n, std::get<1>(records[0]));
    EXPECT_EQ(n, cut.nth_score);
    check_cut(cut);

    auto matches = index.get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("a_record", matches[0].first);
    check_cut(cut);

    // the coordinate counts: ranked by the score reported
    matches = index.get_top_labels(sequence, 1, 0.0, 0.0, true, &cut);
    ASSERT_EQ(1u, matches.size());
    EXPECT_EQ("a_record", matches[0].first);
    EXPECT_EQ(n, matches[0].second);
    check_cut(cut);

    auto counts = index.get_kmer_counts(sequence, 1, 0.0, 0.0, &cut);
    ASSERT_EQ(1u, counts.size());
    EXPECT_EQ("a_record", std::get<0>(counts[0]));
    check_cut(cut);

    // both records: in header order
    records = index.get_kmer_coordinates(sequence, 2, 0.0, 0.0, &cut);
    ASSERT_EQ(2u, records.size());
    EXPECT_EQ("a_record", std::get<0>(records[0]));
    EXPECT_EQ("z_record", std::get<0>(records[1]));
    EXPECT_TRUE(cut.complete);
    EXPECT_EQ(2u, cut.returned);
}

} // namespace
