#include <vector>
#include <unordered_set>

#include <sdsl/int_vector.hpp>
#include "gtest/gtest.h"

#include "../graph/all/test_dbg_helpers.hpp"

#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/config/config.hpp"
#include "cli/query.hpp"
#include "common/vector.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg::cli {
void split_contigs_for_rebalancing(size_t k,
                                   size_t kmers_per_seq,
                                   std::vector<std::pair<std::string, std::vector<uint64_t>>> *contigs);
std::string encode_presence_mask(const sdsl::bit_vector &mask);
std::unique_ptr<graph::AnnotatedDBG>
construct_query_graph(const graph::AnnotatedDBG &full_anno_graph,
                      StringGenerator call_sequences,
                      size_t num_threads,
                      const Config &config);
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


sdsl::bit_vector make_bv(std::initializer_list<uint8_t> bits) {
    sdsl::bit_vector bv(bits.size());
    size_t i = 0;
    for (uint8_t b : bits) bv[i++] = b;
    return bv;
}

TEST(encode_presence_mask, all_zeros) {
    EXPECT_EQ("o8", cli::encode_presence_mask(make_bv({0,0,0,0,0,0,0,0})));
}

TEST(encode_presence_mask, all_ones) {
    EXPECT_EQ("x8", cli::encode_presence_mask(make_bv({1,1,1,1,1,1,1,1})));
}

TEST(encode_presence_mask, leading_run) {
    EXPECT_EQ("x3o2x2o1", cli::encode_presence_mask(make_bv({1,1,1,0,0,1,1,0})));
}

TEST(encode_presence_mask, leading_gap) {
    EXPECT_EQ("o3x5", cli::encode_presence_mask(make_bv({0,0,0,1,1,1,1,1})));
}

TEST(encode_presence_mask, single_bits_alternating) {
    EXPECT_EQ("x1o1x1o1x1o1x1o1", cli::encode_presence_mask(make_bv({1,0,1,0,1,0,1,0})));
}

TEST(encode_presence_mask, single_bit_at_each_end) {
    EXPECT_EQ("x1o6x1", cli::encode_presence_mask(make_bv({1,0,0,0,0,0,0,1})));
}

TEST(encode_presence_mask, empty) {
    EXPECT_EQ("", cli::encode_presence_mask(sdsl::bit_vector(0)));
}

// The cut of a sequence's search as /search states it (SPEC-multigraph-envelope.md section 3)
cli::SeqSearchResult with_cut(const graph::SearchCut &cut) {
    return cli::SeqSearchResult(cli::QuerySequence { 0, "q", "ACGT" },
                                cli::SeqSearchResult::LabelCountVec(), {}, cut);
}

TEST(SeqSearchResult, cut_to_json_complete) {
    graph::SearchCut cut;
    cut.returned = 2;
    Json::Value json = with_cut(cut).cut_to_json();
    EXPECT_EQ(2u, json["returned"].asUInt64());
    EXPECT_TRUE(json["complete"].asBool());
    EXPECT_TRUE(json.isMember("nth_score"));
    EXPECT_TRUE(json["nth_score"].isNull());
    EXPECT_EQ(0u, json["tied"].asUInt64());
    EXPECT_EQ(0u, json["tied_returned"].asUInt64());
    // the median members on counts annotations only
    EXPECT_FALSE(json.isMember("nth_median"));
    EXPECT_FALSE(json.isMember("tied_same_median"));

    cut.abundances = true;
    json = with_cut(cut).cut_to_json();
    EXPECT_TRUE(json.isMember("nth_median"));
    EXPECT_TRUE(json["nth_median"].isNull());
    EXPECT_EQ(0u, json["tied_same_median"].asUInt64());

    // an empty result: nothing returned, complete
    json = with_cut(graph::SearchCut()).cut_to_json();
    EXPECT_EQ(0u, json["returned"].asUInt64());
    EXPECT_TRUE(json["complete"].asBool());
    EXPECT_TRUE(json["nth_score"].isNull());
}

TEST(SeqSearchResult, cut_to_json_cut) {
    graph::SearchCut cut;
    cut.returned = 1;
    cut.complete = false;
    cut.nth_score = 12;
    cut.tied = 2;
    cut.tied_returned = 1;
    cut.abundances = true;
    cut.nth_median = 7;
    cut.tied_same_median = 1;
    Json::Value json = with_cut(cut).cut_to_json();
    EXPECT_EQ(1u, json["returned"].asUInt64());
    EXPECT_FALSE(json["complete"].asBool());
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_EQ(2u, json["tied"].asUInt64());
    EXPECT_EQ(1u, json["tied_returned"].asUInt64());
    EXPECT_EQ(7u, json["nth_median"].asUInt64());
    EXPECT_EQ(1u, json["tied_same_median"].asUInt64());
    EXPECT_EQ(7u, json.size());

    cut.abundances = false;
    json = with_cut(cut).cut_to_json();
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_FALSE(json.isMember("nth_median"));
    EXPECT_EQ(5u, json.size());
}


// A batched search's query graph on an annotation with k-mer counts (R10-SEARCH-COST-1): the
// modes that read no counts keep the binary slice (a fraction of the integer one's memory),
// and the ties of their top N are still ordered by the median abundance, read from the full
// index only where a tie needs it

// a Config made from a query command line, as main() makes it (the swap path restored)
std::unique_ptr<cli::Config> query_config(const std::string &mode) {
    const auto swap = utils::get_swap_path();
    std::vector<std::string> args { "metagraph", "query", "-i", "graph.dbg", "-a",
                                    "annotation.annodbg", "--query-mode", mode, "query.fa" };
    std::vector<char*> argv;
    for (std::string &a : args) {
        argv.push_back(a.data());
    }
    argv.push_back(nullptr);
    auto config = std::make_unique<cli::Config>(static_cast<int>(args.size()), argv.data());
    utils::set_swap_path(swap);
    return config;
}

TEST(construct_query_graph, counts_annotation_binary_slice_reads_tie_medians) {
    const size_t k = 5;
    const std::string sequence = "ACGTTGCATTGACCGTAGGCTTAAC";
    const size_t num_kmers = sequence.size() - k + 1;
    // two labels on the same k-mers, abundance 1 and 100: equal k-mer counts, the median
    // abundance returns "b_high"
    auto graph = test::build_graph_batch<graph::DBGHashOrdered>(k, { sequence });
    Vector<annot::matrix::CSRMatrix::RowValues> rows(graph->max_index());
    graph->map_to_nodes(sequence, [&](graph::DeBruijnGraph::node_index node) {
        rows[graph::AnnotatedDBG::graph_to_anno_index(node)] = { { 0, 1 }, { 1, 100 } };
    });
    annot::LabelEncoder<std::string> encoder;
    encoder.insert_and_encode("a_low");
    encoder.insert_and_encode("b_high");
    graph::AnnotatedDBG full(graph, std::make_unique<annot::IntRowAnnotator>(
            std::make_unique<annot::matrix::CSRMatrix>(std::move(rows), 2), encoder));
    ASSERT_TRUE(full.has_kmer_abundances());

    auto check_cut = [&](const graph::SearchCut &cut) {
        EXPECT_EQ(1u, cut.returned);
        EXPECT_FALSE(cut.complete);
        EXPECT_EQ(num_kmers, cut.nth_score);
        EXPECT_EQ(2u, cut.tied);
        EXPECT_EQ(1u, cut.tied_returned);
        EXPECT_TRUE(cut.abundances);
        EXPECT_EQ(100u, cut.nth_median);
        EXPECT_EQ(1u, cut.tied_same_median);
    };

    for (const std::string mode : { "matches", "signature", "labels" }) {
        SCOPED_TRACE(mode);
        auto config = query_config(mode);
        auto query_graph = cli::construct_query_graph(
                full, [&](const auto &callback) { callback(sequence); }, 1, *config);
        // the binary slice: no integer matrix of the batch's rows
        EXPECT_FALSE(dynamic_cast<const annot::matrix::IntMatrix *>(
                &query_graph->get_annotator().get_matrix()));
        EXPECT_TRUE(query_graph->has_kmer_abundances());

        graph::SearchCut cut;
        auto matches = query_graph->get_top_labels(sequence, 1, 0.0, 0.0, false, &cut);
        ASSERT_EQ(1u, matches.size());
        EXPECT_EQ("b_high", matches[0].first);
        EXPECT_EQ(num_kmers, matches[0].second);
        check_cut(cut);

        auto signatures = query_graph->get_top_label_signatures(sequence, 1, 0.0, 0.0, &cut);
        ASSERT_EQ(1u, signatures.size());
        EXPECT_EQ("b_high", std::get<0>(signatures[0]));
        check_cut(cut);

        // without a cut both are returned, the higher median first
        matches = query_graph->get_top_labels(sequence, 2, 0.0, 0.0, false, &cut);
        ASSERT_EQ(2u, matches.size());
        EXPECT_EQ("b_high", matches[0].first);
        EXPECT_EQ("a_low", matches[1].first);
        EXPECT_TRUE(cut.complete);
    }

    // the counts modes read the integer slice, as before
    for (const std::string mode : { "counts", "counts-sum" }) {
        SCOPED_TRACE(mode);
        auto config = query_config(mode);
        auto query_graph = cli::construct_query_graph(
                full, [&](const auto &callback) { callback(sequence); }, 1, *config);
        EXPECT_TRUE(dynamic_cast<const annot::matrix::IntMatrix *>(
                &query_graph->get_annotator().get_matrix()));
        graph::SearchCut cut;
        auto counts = query_graph->get_kmer_counts(sequence, 1, 0.0, 0.0, &cut);
        ASSERT_EQ(1u, counts.size());
        EXPECT_EQ("b_high", std::get<0>(counts[0]));
        check_cut(cut);
    }
}

} // namespace
