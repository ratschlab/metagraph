#include <algorithm>
#include <atomic>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_set>

#include <sdsl/int_vector.hpp>
#include "gtest/gtest.h"

#include "../graph/all/test_dbg_helpers.hpp"

#include "../annotation/test_matrix_helpers.hpp"
#include "../annotation/coordinate_test_graph.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/align.hpp"
#include "cli/config/config.hpp"
#include "cli/query.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"
#include "common/vector.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
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
SeqSearchResult execute_query(QuerySequence&& sequence,
                              QueryMode query_mode,
                              size_t num_top_labels,
                              double discovery_fraction,
                              double presence_fraction,
                              const graph::AnnotatedDBG &anno_graph);
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
    cut.matched = 3;
    cut.returned = 1;
    cut.complete = false;
    cut.nth_score = 12;
    cut.tied = 2;
    cut.tied_returned = 1;
    cut.abundances = true;
    cut.nth_median = 7;
    cut.tied_same_median = 1;
    cut.nth_median_read = true;
    Json::Value json = with_cut(cut).cut_to_json();
    // matched (train 8, master's member): the labels before the top_labels cut
    EXPECT_EQ(3u, json["matched"].asUInt64());
    EXPECT_EQ(1u, json["returned"].asUInt64());
    EXPECT_FALSE(json["complete"].asBool());
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_EQ(2u, json["tied"].asUInt64());
    EXPECT_EQ(1u, json["tied_returned"].asUInt64());
    EXPECT_EQ(7u, json["nth_median"].asUInt64());
    EXPECT_EQ(1u, json["tied_same_median"].asUInt64());
    EXPECT_EQ(8u, json.size());

    // the N-th score of one label: no median read, none stated
    cut.tied = 1;
    cut.nth_median = 0;
    cut.nth_median_read = false;
    json = with_cut(cut).cut_to_json();
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_TRUE(json.isMember("nth_median"));
    EXPECT_TRUE(json["nth_median"].isNull());
    EXPECT_EQ(1u, json["tied_same_median"].asUInt64());
    EXPECT_EQ(8u, json.size());

    cut.abundances = false;
    json = with_cut(cut).cut_to_json();
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_FALSE(json.isMember("nth_median"));
    EXPECT_EQ(6u, json.size());
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
    for (graph::DeBruijnGraph::node_index node : graph->node_range(sequence)) {
        rows[graph::AnnotatedDBG::graph_to_anno_index(node)] = { { 0, 1 }, { 1, 100 } };
    }
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

// Batched against unbatched search on an annotation with k-mer counts whose abundances differ
// per k-mer, with query k-mers repeated (GPT B-gap): every query mode that reads such an
// annotation answers the same rows, in the same order, with the same cut, whether the query
// runs on the full index or on a batch's query graph (whose rows are the batch's distinct
// k-mers: a k-mer repeated in a query must count once per position, also in the medians that
// order the ties). The C5 order is checked against medians computed per query position here.

// the integer matrices of a counts annotation the search reads: the CSR form, and an
// integer CSC over a BRWT (int_brwt), read in batches of rows
enum class CountsForm { CSR, INT_BRWT };

std::unique_ptr<graph::AnnotatedDBG>
counts_annotation(std::shared_ptr<graph::DeBruijnGraph> graph,
                  const std::vector<std::string> &labels,
                  const std::function<uint64_t(size_t label, const std::string &kmer)> &value,
                  CountsForm form) {
    using namespace annot::matrix;
    const uint64_t num_rows = graph->max_index();
    Vector<CSRMatrix::RowValues> rows(num_rows);
    graph->call_kmers([&](graph::DeBruijnGraph::node_index node, const std::string &kmer) {
        auto &row = rows[graph::AnnotatedDBG::graph_to_anno_index(node)];
        for (size_t j = 0; j < labels.size(); ++j) {
            if (uint64_t v = value(j, kmer))
                row.emplace_back(j, v);
        }
    });
    annot::LabelEncoder<std::string> encoder;
    for (const auto &label : labels) {
        encoder.insert_and_encode(label);
    }
    if (form == CountsForm::CSR) {
        return std::make_unique<graph::AnnotatedDBG>(graph, std::make_unique<annot::IntRowAnnotator>(
                std::make_unique<CSRMatrix>(std::move(rows), labels.size()), encoder));
    }
    test::BitVectorPtrArray columns;
    std::vector<sdsl::int_vector<>> values(labels.size());
    std::vector<sdsl::bit_vector> bits(labels.size(), sdsl::bit_vector(num_rows, false));
    std::vector<std::vector<uint64_t>> column_values(labels.size());
    for (uint64_t r = 0; r < num_rows; ++r) {
        for (const auto &[j, v] : rows[r]) {
            bits[j][r] = true;
            column_values[j].push_back(v);
        }
    }
    std::vector<annot::CountsVector> counts;
    for (size_t j = 0; j < labels.size(); ++j) {
        columns.push_back(std::make_unique<bit_vector_stat>(std::move(bits[j])));
        sdsl::int_vector<> v(column_values[j].size());
        for (size_t i = 0; i < v.size(); ++i) {
            v[i] = column_values[j][i];
        }
        counts.emplace_back(v);
    }
    auto brwt = test::build_matrix_from_columns<BRWT>(columns, num_rows);
    return std::make_unique<graph::AnnotatedDBG>(graph, std::make_unique<annot::IntMultiBRWTAnnotator>(
            std::make_unique<CSCMatrix<BRWT, annot::CountsVector>>(std::move(brwt),
                                                                   std::move(counts)),
            encoder));
}

template <class Graph>
void check_batched_parity(graph::DeBruijnGraph::Mode graph_mode, CountsForm form) {
    const size_t k = 11;
    std::mt19937 rng(11);
    auto random_seq = [&](size_t length) {
        std::string s(length, 'A');
        for (char &c : s) {
            c = "ACGT"[rng() % 4];
        }
        return s;
    };
    // A has 20 k-mers, B 30, C 15: the query has A three times, so 60 of A's positions and 30
    // of B's, and a median that counted A's k-mers once (20 against 30) would land on the
    // other side
    const std::string A = random_seq(20 + k - 1);
    const std::string B = random_seq(30 + k - 1);
    const std::string C = random_seq(15 + k - 1);
    auto graph = test::build_graph_batch<Graph>(k, { A, B, C }, graph_mode);
    std::map<std::string, std::pair<char, size_t>> where;
    for (const auto &[name, segment] : { std::make_pair('A', A), std::make_pair('B', B),
                                         std::make_pair('C', C) }) {
        for (size_t i = 0; i + k <= segment.size(); ++i) {
            where[segment.substr(i, k)] = { name, i };
        }
    }
    // the segment and offset of a k-mer of the graph, on either strand in a canonical one
    auto locate = [&](const std::string &kmer) -> std::pair<char, size_t> {
        if (auto it = where.find(kmer); it != where.end())
            return it->second;
        std::string rc = kmer;
        reverse_complement(rc.begin(), rc.end());
        if (auto it = where.find(rc); it != where.end())
            return it->second;
        return { '-', 0 };
    };
    const std::vector<std::string> labels { "u_wide", "v_a_only", "w_like_x", "x_low_a",
                                            "y_high_a", "z_flat", "t_part_b" };
    auto value = [&](size_t label, const std::string &kmer) -> uint64_t {
        const auto [segment, i] = locate(kmer);
        const bool a = segment == 'A', b = segment == 'B', c = segment == 'C';
        switch (label) {
            case 0: return a || b || c ? 7 + i % 2 : 0;           // u: 7 or 8 everywhere
            case 1: return a ? 3 + i % 4 : 0;                     // v: A only, uneven
            case 2: case 3: return a ? 1 + i % 2 : b ? 100 + i : 0;  // w = x: low on A
            case 4: return a ? 100 + i : b ? 1 + i % 3 : 0;      // y: high on A
            case 5: return a || b ? 50 : 0;                       // z: flat
            case 6: return b && i >= 10 && i < 20 ? 1 + i : 0;    // t: part of B
        }
        return 0;
    };
    auto full = counts_annotation(graph, labels, value, form);
    ASSERT_TRUE(full->has_kmer_abundances());

    std::string reversed = A + B + A + A;
    reverse_complement(reversed.begin(), reversed.end());
    const std::vector<std::string> queries {
        A + B + A + A,
        // A's k-mers only, five times
        A + A + A + A + A,
        B + "T" + A + C + A,
        reversed,
        random_seq(60),
    };

    // the C5 order of the labels by the positions of |query|, each position once: the
    // reference the engine's rows are checked against (MATCHES)
    auto expected_matches = [&](const std::string &query, size_t n) {
        std::vector<std::tuple<size_t, uint64_t, std::string>> rows;
        for (size_t j = 0; j < labels.size(); ++j) {
            std::vector<uint64_t> values;
            for (graph::DeBruijnGraph::node_index node : graph->node_range(query)) {
                if (!node)
                    continue;
                std::string kmer = graph->get_node_sequence(node);
                if (uint64_t v = value(j, kmer))
                    values.push_back(v);
            }
            if (values.size())
                rows.emplace_back(values.size(), graph::lower_median(values), labels[j]);
        }
        std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
            if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) > std::get<0>(b);
            if (std::get<1>(a) != std::get<1>(b))
                return std::get<1>(a) > std::get<1>(b);
            return std::get<2>(a) < std::get<2>(b);
        });
        std::vector<std::string> names;
        for (size_t r = 0; r < std::min(n, rows.size()); ++r) {
            names.push_back(std::get<2>(rows[r]));
        }
        return names;
    };

    for (const std::string mode : { "labels", "matches", "counts-sum", "counts", "signature" }) {
        auto config = query_config(mode);
        auto query_graph = cli::construct_query_graph(*full, [&](const auto &callback) {
            for (const auto &query : queries) {
                callback(query);
            }
        }, 1, *config);
        for (size_t q = 0; q < queries.size(); ++q) {
            for (size_t n : { 0, 1, 2, 3, 4, 5, 6, 7, 100 }) {
                SCOPED_TRACE(mode + " query " + std::to_string(q) + " top " + std::to_string(n));
                auto run = [&](const graph::AnnotatedDBG &index) {
                    return cli::execute_query(cli::QuerySequence { q, "q", queries[q] },
                                              config->query_mode, n, 0.0, 0.0, index);
                };
                const cli::SeqSearchResult direct = run(*full);
                const cli::SeqSearchResult batched = run(*query_graph);
                EXPECT_EQ(direct.to_json(true, k), batched.to_json(true, k));
                EXPECT_EQ(direct.cut_to_json(), batched.cut_to_json());
                EXPECT_EQ(direct.to_string("\t", false, true, k),
                          batched.to_string("\t", false, true, k));
                if (std::string(mode) == "matches") {
                    std::vector<std::string> names;
                    for (const auto &[label, count]
                            : std::get<cli::SeqSearchResult::LabelCountVec>(direct.get_result())) {
                        names.push_back(label);
                    }
                    EXPECT_EQ(expected_matches(queries[q], n), names);
                }
            }
        }
    }
}

TEST(construct_query_graph, batched_counts_search_answers_as_unbatched) {
    check_batched_parity<graph::DBGHashOrdered>(graph::DeBruijnGraph::BASIC, CountsForm::CSR);
#if ! _PROTEIN_GRAPH
    // a Protein build has no reverse complement, so no canonical graph
    check_batched_parity<graph::DBGHashOrdered>(graph::DeBruijnGraph::CANONICAL, CountsForm::CSR);
#endif
    check_batched_parity<graph::DBGSuccinct>(graph::DeBruijnGraph::BASIC, CountsForm::CSR);
    check_batched_parity<graph::DBGHashOrdered>(graph::DeBruijnGraph::BASIC, CountsForm::INT_BRWT);
#if ! _PROTEIN_GRAPH
    check_batched_parity<graph::DBGSuccinct>(graph::DeBruijnGraph::CANONICAL, CountsForm::INT_BRWT);
#endif
}

// An int_brwt annotation's matrix that counts the integer values read from it
class CountingIntBRWT : public annot::matrix::CSCMatrix<annot::matrix::BRWT, annot::CountsVector> {
  public:
    using CSCMatrix::CSCMatrix;
    std::vector<RowValues> get_row_values(std::span<const Row> rows,
                                          size_t num_threads = 1) const override {
        ++reads;
        rows_read += rows.size();
        return CSCMatrix::get_row_values(rows, num_threads);
    }
    mutable size_t reads = 0;
    mutable size_t rows_read = 0;
};

// GPT B3: a batch on a counts annotation keeps the binary slice in the modes that read no counts,
// and reads no integer value where no median orders a tie: N = 0, every score its own, and a
// tie of one at the cut (only a shared score reads the medians, once, a repeated k-mer's row
// once). Unbatched alike
TEST(construct_query_graph, no_values_read_without_a_tie) {
    const size_t k = 11;
    std::mt19937 rng(5);
    std::string segment(60 + k - 1, 'A');
    for (char &c : segment) {
        c = "ACGT"[rng() % 4];
    }
    auto graph = test::build_graph_batch<graph::DBGHashOrdered>(k, { segment });
    std::map<std::string, size_t> position;
    for (size_t i = 0; i + k <= segment.size(); ++i) {
        position[segment.substr(i, k)] = i;
    }
    // a: 60 k-mers, b: 40, c: 20, uneven; d: b's k-mers with other values (the tie with b)
    const std::vector<std::string> labels { "a", "b", "c", "d" };
    const size_t covered[] = { 60, 40, 20, 40 };
    auto value = [&](size_t label, const std::string &kmer) -> uint64_t {
        const size_t i = position.at(kmer);
        return i < covered[label] ? 1 + (i * (label + 3)) % 7 : 0;
    };
    const uint64_t num_rows = graph->max_index();
    std::vector<sdsl::bit_vector> bits(labels.size(), sdsl::bit_vector(num_rows, false));
    std::vector<std::vector<uint64_t>> values(labels.size());
    std::vector<std::pair<uint64_t, std::string>> kmers;
    graph->call_kmers([&](graph::DeBruijnGraph::node_index node, const std::string &kmer) {
        kmers.emplace_back(graph::AnnotatedDBG::graph_to_anno_index(node), kmer);
    });
    std::sort(kmers.begin(), kmers.end());
    for (const auto &[row, kmer] : kmers) {
        for (size_t j = 0; j < labels.size(); ++j) {
            if (uint64_t v = value(j, kmer)) {
                bits[j][row] = true;
                values[j].push_back(v);
            }
        }
    }
    test::BitVectorPtrArray columns;
    std::vector<annot::CountsVector> counts;
    for (size_t j = 0; j < labels.size(); ++j) {
        columns.push_back(std::make_unique<bit_vector_stat>(std::move(bits[j])));
        sdsl::int_vector<> v(values[j].size());
        for (size_t i = 0; i < v.size(); ++i) {
            v[i] = values[j][i];
        }
        counts.emplace_back(v);
    }
    auto matrix = std::make_unique<CountingIntBRWT>(
            test::build_matrix_from_columns<annot::matrix::BRWT>(columns, num_rows),
            std::move(counts));
    const CountingIntBRWT &counted = *matrix;
    annot::LabelEncoder<std::string> encoder;
    for (const auto &label : labels) {
        encoder.insert_and_encode(label);
    }
    graph::AnnotatedDBG full(graph, std::make_unique<annot::IntMultiBRWTAnnotator>(
            std::move(matrix), encoder));
    ASSERT_TRUE(full.has_kmer_abundances());

    // the segment, its first 20 k-mers three times more
    const std::string query = segment + "T" + segment.substr(0, 20 + k - 1)
                            + "T" + segment.substr(0, 20 + k - 1) + "T"
                            + segment.substr(0, 20 + k - 1);
    for (const std::string mode : { "matches", "signature" }) {
        auto config = query_config(mode);
        auto query_graph = cli::construct_query_graph(
                full, [&](const auto &callback) { callback(query); }, 1, *config);
        EXPECT_FALSE(dynamic_cast<const annot::matrix::IntMatrix *>(
                &query_graph->get_annotator().get_matrix())) << mode;
        for (const graph::AnnotatedDBG *index : { query_graph.get(), &full }) {
            SCOPED_TRACE(mode + (index == &full ? " unbatched" : " batched"));
            auto run = [&](size_t n) {
                return cli::execute_query(cli::QuerySequence { 0, "q", query },
                                          config->query_mode, n, 0.0, 0.0, *index);
            };
            const size_t reads = counted.reads;
            const size_t rows_read = counted.rows_read;
            // N = 0, and a cut whose N-th score has one label (positions: a 120, b and d 100,
            // c 80): no value read
            for (size_t n : { 0, 1 }) {
                EXPECT_FALSE(run(n).get_cut().nth_median_read) << n;
                EXPECT_EQ(reads, counted.reads) << n;
            }
            // a cut in the tie of b and d: the medians are read once, each of the query's 60
            // distinct k-mers' rows once (the 60 repeated positions not again)
            size_t searches = 1;
#ifndef NDEBUG
            // (a build with assertions checks a signature against get_top_labels, a second
            // search that reads the medians again)
            if (mode == "signature")
                searches = 2;
#endif
            const cli::SeqSearchResult tie = run(2);
            EXPECT_TRUE(tie.get_cut().nth_median_read);
            EXPECT_EQ(2u, tie.get_cut().tied);
            EXPECT_EQ(reads + searches, counted.reads);
            EXPECT_EQ(rows_read + 60 * searches, counted.rows_read);
        }
    }
}

// REVIEW 9726cb420 [P2], lane T6-W2: query_fasta hands each record's result to its callback as
// an rvalue, so that a callback that keeps it (the server's /search, which moves it into the
// pair's results) takes its rows and their per-position values over instead of copying them:
// a copy was held beside the original, outside the /search answer budget's values account.
// Every mode, the sequential path (a record mapping, or query_coords) and the batched one
TEST(query_fasta, a_record_kept_by_the_callback_is_moved_not_copied) {
    std::mt19937 rng(17);
    auto random_sequence = [&](size_t length) {
        std::string s(length, 'A');
        for (char &c : s) {
            c = "ACGT"[rng() % 4];
        }
        return s;
    };
    const size_t k = 7;
    std::vector<test::CoordinateRecord> records;
    for (size_t i = 0; i < 12; ++i) {
        records.push_back({ random_sequence(40 + rng() % 30),
                            "column_with_a_long_name_" + std::to_string(i % 3),
                            "record_with_a_long_name_" + std::to_string(i) });
    }
    utils::TempFile fasta;
    for (size_t q = 0; q < 6; ++q) {
        fasta.ofstream() << ">q" << q << "\n" << records[q].sequence << records[q + 3].sequence
                         << records[q].sequence << "\n";
    }
    fasta.ofstream().close();

    for (bool record_mapping : { true, false }) {
        auto index = test::coordinate_test_graph(k, records, false, false, false, record_mapping);
        for (const std::string mode : { "labels", "matches", "signature", "coords" }) {
            for (unsigned long long batch : { 0ull, 100'000'000ull }) {
                SCOPED_TRACE(::testing::Message() << "record_mapping " << record_mapping
                             << " mode " << mode << " batch " << batch);
                auto config = query_config(mode);
                config->discovery_fraction = 0;
                config->num_top_labels = 1000;
                config->query_batch_size = batch;
                std::mutex mutex;
                std::vector<cli::SeqSearchResult> kept;
                size_t moved = 0;
                // a callback written as the server's: it keeps each result it is handed
                cli::query_fasta(fasta.name(), [&](auto &&result) {
                    auto rows_buffer = [](const cli::SeqSearchResult &r) {
                        return std::visit([](const auto &rows) {
                            return static_cast<const void *>(rows.data());
                        }, r.get_result());
                    };
                    const void *before = rows_buffer(result);
                    std::lock_guard<std::mutex> lock(mutex);
                    kept.emplace_back(std::move(result));
                    // the rows' buffer is the one the query made: taken over, not copied
                    if (before && rows_buffer(kept.back()) == before)
                        ++moved;
                }, *config, *index, nullptr);
                ASSERT_EQ(6u, kept.size());
                size_t with_rows = 0;
                for (const auto &result : kept) {
                    std::visit([&](const auto &rows) { with_rows += !rows.empty(); },
                               result.get_result());
                }
                EXPECT_EQ(6u, with_rows);
                EXPECT_EQ(with_rows, moved);
            }
        }
    }
}

// Lane T7-P1, the /search partial refusal per unit (DESIGN-search-partial-refusal.md section
// 7, tests 4, 5 and 23): query_fasta with a request's budget makes each record's values account
// on its task; with `partial_refusal` a unit past the bound is answered as its refused stub
// (no row, the cut the search computed, marked refused; the alignment and the record's members
// kept), and every other unit as without a budget, byte for byte.

std::string compact_json(const Json::Value &value) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    return Json::writeString(builder, value);
}

// the results of |fasta| on |index| by record, with |budget| (null: none)
std::map<size_t, cli::SeqSearchResult>
search_units(const std::string &fasta, const cli::Config &config,
             const graph::AnnotatedDBG &index, graph::SearchBudget *budget,
             const graph::align::DBGAlignerConfig *aligner = nullptr) {
    std::mutex mutex;
    std::map<size_t, cli::SeqSearchResult> results;
    cli::query_fasta(fasta, [&](cli::SeqSearchResult &&result) {
        std::lock_guard<std::mutex> lock(mutex);
        const size_t id = result.get_sequence().id;
        results.emplace(id, std::move(result));
    }, config, index, aligner, budget, "pair (p, p.dbg, p.annodbg)");
    return results;
}

// each unit either answered as |base| or refused alone on its values: its stub
void expect_answered_or_stub(const std::map<size_t, cli::SeqSearchResult> &base,
                             const std::map<size_t, cli::SeqSearchResult> &got,
                             const graph::SearchBudget &budget, size_t *stubs) {
    ASSERT_EQ(base.size(), got.size());
    *stubs = 0;
    for (const auto &[id, expected] : base) {
        const cli::SeqSearchResult &result = got.at(id);
        if (result.get_cut().refused == graph::SearchCut::Refused::NO) {
            EXPECT_EQ(compact_json(expected.to_json(false, 7)), compact_json(result.to_json(false, 7)))
                    << id;
            EXPECT_EQ(compact_json(expected.cut_to_json()), compact_json(result.cut_to_json()))
                    << id;
            continue;
        }
        ++*stubs;
        EXPECT_EQ(graph::SearchCut::Refused::VALUES, result.get_cut().refused) << id;
        // the member set of the cut the search computed
        EXPECT_EQ(expected.get_cut().abundances, result.get_cut().abundances) << id;
        cli::SeqSearchResult stub = expected;
        stub.refuse(graph::SearchCut::Refused::VALUES, budget.max_bytes() >> 20);
        EXPECT_EQ(compact_json(stub.cut_to_json()), compact_json(result.cut_to_json())) << id;
        // the record's members (its alignment too), no row
        EXPECT_EQ(compact_json(stub.to_json(false, 7)), compact_json(result.to_json(false, 7)))
                << id;
        std::visit([&](const auto &rows) { EXPECT_TRUE(rows.empty()) << id; },
                   result.get_result());
    }
    // one refusal noted per stub: a stub is never refused twice (its handover not charged)
    EXPECT_EQ(*stubs, budget.unit_refusals().size());
    for (const auto &refusal : budget.unit_refusals()) {
        EXPECT_EQ(graph::SearchCut::Refused::VALUES, got.at(refusal.record).get_cut().refused);
        EXPECT_EQ("record " + std::to_string(refusal.record) + " on pair (p, p.dbg, p.annodbg)",
                  refusal.where);
        EXPECT_EQ("values", std::string(refusal.account));
        EXPECT_LE(refusal.held, budget.max_bytes());
        EXPECT_GT(refusal.held + refusal.charge, budget.max_bytes());
    }
}

// Test 4: a refusal at each charge of a unit's search (every prefix of its charges, found from
// the last charge back: the bound one byte below what the unit held before the charge that
// refused it), on every path (the header paths with a record mapping, the column paths
// (filter_and_aggregate, W2's store scratch) without one, a signature at discovery_fraction 1,
// the handover in a plain mode, the batched path), each answered as its stub with the computed
// cut's member set, every other unit as without a budget. Test 5: with `align` the stub keeps
// the alignment
TEST(query_fasta, a_unit_refused_at_any_charge_is_its_stub) {
    std::mt19937 rng(73);
    auto random_sequence = [&](size_t length) {
        std::string s(length, 'A');
        for (char &c : s) {
            c = "ACGT"[rng() % 4];
        }
        return s;
    };
    const size_t k = 7;
    std::vector<test::CoordinateRecord> records;
    for (size_t i = 0; i < 12; ++i) {
        records.push_back({ random_sequence(40 + rng() % 30), "column_" + std::to_string(i % 4),
                            "record_" + std::to_string(i) });
    }
    // record 2, the largest unit: three records and one again
    utils::TempFile fasta;
    for (size_t q = 0; q < 4; ++q) {
        fasta.ofstream() << ">q" << q << "\n" << records[q].sequence;
        if (q == 2)
            fasta.ofstream() << records[5].sequence << records[6].sequence << records[5].sequence;
        fasta.ofstream() << "\n";
    }
    fasta.ofstream().close();

    // the indexes: coordinates with a record mapping (the header paths), without one (the
    // column paths), and k-mer counts (the column paths of an annotation with abundances: the
    // medians' members)
    enum Index { HEADERS, COLUMNS, COUNTS };
    std::vector<std::string> sequences;
    for (const auto &record : records) {
        sequences.push_back(record.sequence);
    }
    std::unique_ptr<graph::AnnotatedDBG> indexes[] = {
        test::coordinate_test_graph(k, records, false, false, false, true),
        test::coordinate_test_graph(k, records, false, false, false, false),
        counts_annotation(test::build_graph_batch<graph::DBGHashOrdered>(k, sequences),
                          { "l0", "l1", "l2", "l3" },
                          [&](size_t label, const std::string &kmer) -> uint64_t {
                              const size_t h = std::hash<std::string>()(kmer);
                              return (h >> label) % 3 ? 1 + (h >> 8) % (label + 2) : 0;
                          }, CountsForm::CSR)
    };
    ASSERT_TRUE(indexes[COUNTS]->has_kmer_abundances());
    struct Case {
        Index index;
        std::string mode;
        double discovery;
        unsigned long long batch;
        bool align;
    };
    std::vector<Case> cases;
    for (Index index : { HEADERS, COLUMNS, COUNTS }) {
        for (const std::string mode : { "matches", "signature", "counts", "coords" }) {
            if (index == COUNTS && mode == "coords")
                continue;
            for (unsigned long long batch : { 0ull, 100'000'000ull }) {
                cases.push_back({ index, mode, 0.0, batch, false });
            }
        }
        cases.push_back({ index, "signature", 1.0, 0, false });
        cases.push_back({ index, index == COUNTS ? "counts" : "coords", 0.0, 0, true });
        cases.push_back({ index, "signature", 0.0, 100'000'000ull, true });
    }
    size_t plain_sites = 0, value_sites = 0, with_medians = 0;
    for (const Case &c : cases) {
        SCOPED_TRACE(::testing::Message() << c.mode << " index " << c.index
                     << " discovery " << c.discovery << " batch " << c.batch << " align "
                     << c.align);
        const auto &index = indexes[c.index];
        auto config = query_config(c.mode);
        config->discovery_fraction = c.discovery;
        config->num_top_labels = 1000;
        config->query_batch_size = c.batch;
        std::unique_ptr<graph::align::DBGAlignerConfig> aligner;
        if (c.align) {
            aligner = std::make_unique<graph::align::DBGAlignerConfig>(
                    cli::initialize_aligner_config(*config, index->get_graph()));
        }
        const auto base = search_units(fasta.name(), *config, *index, nullptr, aligner.get());
        ASSERT_EQ(4u, base.size());
        with_medians += base.at(2).cut_to_json().isMember("nth_median");
        if (c.align) {
            ASSERT_TRUE(base.at(2).get_alignment());
        }
        // with a budget that fits: the same results, the units' values noted
        graph::SearchBudget large(uint64_t(1) << 40, true);
        size_t stubs = 0;
        expect_answered_or_stub(base, search_units(fasta.name(), *config, *index, &large,
                                                   aligner.get()),
                                large, &stubs);
        EXPECT_EQ(0u, stubs);
        const uint64_t most = large.values_max();
        ASSERT_LT(0u, most);
        // the bound just below the largest unit, then below each charge of it, back to its
        // first
        uint64_t bound = most - 1;
        std::optional<size_t> target;
        size_t sites = 0;
        for (; sites < 400; ++sites) {
            graph::SearchBudget budget(bound, true);
            const auto got = search_units(fasta.name(), *config, *index, &budget, aligner.get());
            EXPECT_FALSE(budget.passed());
            expect_answered_or_stub(base, got, budget, &stubs);
            ASSERT_LT(0u, stubs);
            if (!target)
                target = budget.unit_refusals().front().record;
            uint64_t held = 0;
            bool found = false;
            for (const auto &refusal : budget.unit_refusals()) {
                if (refusal.record == *target) {
                    held = refusal.held;
                    found = true;
                }
            }
            ASSERT_TRUE(found) << "the largest unit answered at " << bound;
            if (!held)
                break;
            bound = held - 1;
        }
        ASSERT_LT(sites, 400u);
        ++sites;
        // a plain mode charges its rows once, as they are handed over; the per-position modes
        // their values as they grow, then their rows
        if (c.mode == "matches") {
            EXPECT_EQ(1u, sites);
            plain_sites += sites;
        } else {
            EXPECT_LE(2u, sites);
            value_sites += sites;
        }
        // without the field: the request refused on its values, nothing more refused alone
        graph::SearchBudget whole(most - 1, false);
        EXPECT_THROW(search_units(fasta.name(), *config, *index, &whole, aligner.get()),
                     graph::SearchAnswerTooLarge);
        EXPECT_TRUE(whole.passed());
        EXPECT_EQ("values", std::string(whole.account()));
        EXPECT_TRUE(whole.unit_refusals().empty());
    }
    EXPECT_LT(0u, plain_sites);
    EXPECT_LT(0u, value_sites);
    EXPECT_LT(0u, with_medians);
}

// Without a budget query_fasta is as before (the CLI's): no account, no refusal; a stub's
// handover is never charged (a unit refused on its values is noted once)
TEST(query_fasta, a_unit_without_a_budget_has_no_account) {
    const size_t k = 7;
    std::vector<test::CoordinateRecord> records { { "ACGTTGCATTGACCGTAGGCTTAACGGT", "c", "r" } };
    auto index = test::coordinate_test_graph(k, records);
    utils::TempFile fasta;
    fasta.ofstream() << ">q\n" << records[0].sequence << "\n";
    fasta.ofstream().close();
    auto config = query_config("coords");
    config->discovery_fraction = 0;
    const auto results = search_units(fasta.name(), *config, *index, nullptr);
    ASSERT_EQ(1u, results.size());
    EXPECT_EQ(graph::SearchCut::Refused::NO, results.at(0).get_cut().refused);
    EXPECT_FALSE(results.at(0).cut_to_json().isMember("refused"));
}

// Lane T7-P2, a /search whose client is gone: query_fasta with a stop starts no record once the
// stop answered true (asked before each record's work, and on the batched path before each
// batch and each of its records), lets the records running finish, notes how many it started
// and throws SearchStopped. A stop that never fires, or none, answers as before.

// the results of |fasta| by record (the server's callback, keeping each), with |budget|,
// |stop| (null: none) and |aligner| (null: no alignment); |callbacks| counts the results handed
// over
std::map<size_t, cli::SeqSearchResult>
search_with_stop(const std::string &fasta, const cli::Config &config,
                 const graph::AnnotatedDBG &index, graph::SearchBudget *budget,
                 cli::SearchStop *stop, std::atomic<size_t> *callbacks,
                 const graph::align::DBGAlignerConfig *aligner = nullptr) {
    std::mutex mutex;
    std::map<size_t, cli::SeqSearchResult> results;
    cli::query_fasta(fasta, [&](cli::SeqSearchResult &&result) {
        callbacks->fetch_add(1);
        std::lock_guard<std::mutex> lock(mutex);
        const size_t id = result.get_sequence().id;
        results.emplace(id, std::move(result));
    }, config, index, aligner, budget, "the graph", 0, stop);
    return results;
}

TEST(SearchStop, NoRecordStartsOnceStopped) {
    std::mt19937 rng(29);
    auto random_sequence = [&](size_t length) {
        std::string s(length, 'A');
        for (char &c : s) {
            c = "ACGT"[rng() % 4];
        }
        return s;
    };
    const size_t k = 7;
    std::vector<test::CoordinateRecord> records;
    for (size_t i = 0; i < 12; ++i) {
        records.push_back({ random_sequence(40 + rng() % 30), "column_" + std::to_string(i % 4),
                            "record_" + std::to_string(i) });
    }
    const size_t num_records = 40;
    utils::TempFile fasta;
    for (size_t q = 0; q < num_records; ++q) {
        fasta.ofstream() << ">q" << q << "\n" << records[q % 12].sequence
                         << records[(q + 5) % 12].sequence << "\n";
    }
    fasta.ofstream().close();
    // the records as the server counts them for its log line
    std::ifstream in(fasta.name());
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    ASSERT_EQ(num_records, cli::count_query_records(text, false));
    ASSERT_EQ(2 * num_records, cli::count_query_records(text, true));

    const unsigned int threads = get_num_threads();
    // the sequential path (a record mapping; query_coords) and the batched one (no record
    // mapping, a mode other than coords): one batch, its records queried after its graph
    struct Case {
        bool record_mapping;
        std::string mode;
        unsigned long long batch;
    };
    for (const Case &c : { Case { true, "coords", 0 }, Case { true, "signature", 0 },
                           Case { false, "matches", 100'000'000 },
                           Case { false, "signature", 100'000'000 } }) {
        auto index = test::coordinate_test_graph(k, records, false, false, false,
                                                 c.record_mapping);
        auto config = query_config(c.mode);
        config->discovery_fraction = 0;
        config->num_top_labels = 1000;
        config->query_batch_size = c.batch;
        for (unsigned int t : { 1u, 4u }) {
            SCOPED_TRACE(::testing::Message() << c.mode << " record_mapping " << c.record_mapping
                         << " batch " << c.batch << " threads " << t);
            set_num_threads(t);
            std::atomic<size_t> callbacks { 0 };
            const auto base = search_with_stop(fasta.name(), *config, *index, nullptr, nullptr,
                                               &callbacks);
            ASSERT_EQ(num_records, base.size());
            ASSERT_EQ(num_records, callbacks.load());

            // a stop that never fires (with a budget): every record, as without one
            graph::SearchBudget large(uint64_t(1) << 40, true);
            std::atomic<size_t> asked { 0 };
            cli::SearchStop never([&]() -> const char * {
                asked.fetch_add(1);
                return nullptr;
            });
            callbacks = 0;
            const auto all = search_with_stop(fasta.name(), *config, *index, &large, &never,
                                              &callbacks);
            ASSERT_EQ(num_records, all.size());
            EXPECT_EQ(num_records, callbacks.load());
            for (const auto &[id, result] : base) {
                EXPECT_EQ(compact_json(result.to_json(false, k)),
                          compact_json(all.at(id).to_json(false, k))) << id;
                EXPECT_EQ(compact_json(result.cut_to_json()),
                          compact_json(all.at(id).cut_to_json())) << id;
            }
            EXPECT_FALSE(never.stopped());
            EXPECT_EQ(num_records, never.records_started());
            // asked once per record (and on the batched path once per batch read and started)
            EXPECT_EQ(num_records + (c.batch ? 2 : 0), asked.load());

            // a stop that fires at its |n|-th question, from the first to past the last. With a
            // budget whose bound nothing reaches, so that the stop is asked beside the values
            // accounts the tasks make; nothing about the budget is asserted (this callback
            // charges no answer and nothing passes 1 TiB, so it would hold without the stop,
            // and the budget is the request's own: no count across requests to leak into)
            for (size_t n : { size_t(1), size_t(2), size_t(6), size_t(num_records) }) {
                SCOPED_TRACE(::testing::Message() << "stop at question " << n);
                graph::SearchBudget budget(uint64_t(1) << 40, true);
                std::atomic<size_t> questions { 0 };
                cli::SearchStop stop([&]() -> const char * {
                    return questions.fetch_add(1) + 1 >= n ? "the test's" : nullptr;
                });
                callbacks = 0;
                std::map<size_t, cli::SeqSearchResult> got;
                EXPECT_THROW(got = search_with_stop(fasta.name(), *config, *index, &budget,
                                                    &stop, &callbacks),
                             cli::SearchStopped);
                EXPECT_TRUE(stop.stopped());
                EXPECT_STREQ("the test's", stop.reason());
                // the predicate is not asked again once it answered true
                const size_t asked_then = questions.load();
                EXPECT_TRUE(stop());
                EXPECT_EQ(asked_then, questions.load());
                // no record started after the stop: those started (each handed over) are fewer
                // than all, and at most the questions answered false
                const size_t started = stop.records_started();
                EXPECT_EQ(started, callbacks.load());
                EXPECT_LT(started, num_records);
                EXPECT_LE(started, n - 1);
                if (t == 1) {
                    // one thread: the records start in order, each after one question (the
                    // batched path asks twice before its batch's records: when it reads the
                    // batch and before its query graph)
                    const size_t before = c.batch ? 2 : 0;
                    EXPECT_EQ(n - 1 > before ? n - 1 - before : 0, started);
                }
            }
        }
    }
    set_num_threads(threads);
}

// a stop that would fire only after every record has started is never asked again by
// query_fasta: the record is answered, nothing thrown (the server's next pair would see it)
TEST(SearchStop, AStopAfterTheLastRecordStopsNothing) {
    const size_t k = 7;
    std::vector<test::CoordinateRecord> records { { "ACGTTGCATTGACCGTAGGCTTAACGGT", "c", "r" } };
    auto index = test::coordinate_test_graph(k, records);
    utils::TempFile fasta;
    fasta.ofstream() << ">q\n" << records[0].sequence << "\n";
    fasta.ofstream().close();
    auto config = query_config("coords");
    config->discovery_fraction = 0;
    std::atomic<size_t> callbacks { 0 };
    std::atomic<size_t> questions { 0 };
    cli::SearchStop stop([&]() -> const char * {
        return questions.fetch_add(1) >= 1 ? "the test's" : nullptr;
    });
    const auto results = search_with_stop(fasta.name(), *config, *index, nullptr, &stop,
                                          &callbacks);
    ASSERT_EQ(1u, results.size());
    EXPECT_EQ(1u, questions.load());
    EXPECT_FALSE(stop.stopped());
    EXPECT_EQ(1u, stop.records_started());
    // a later question (the server's next pair) fires; this call is over
    EXPECT_TRUE(stop());
}


// the batched path with alignments on the full graph and one batch or several (Fable review
// m2): each alignment asks the stop first and the batch asks again before its query graph, so a
// stop during a batch's alignments makes neither the alignments left nor its query graph; the
// batches are read one at a time, each asked for. The records started are those whose query on
// their batch's query graph began
TEST(SearchStop, TheBatchedPathStopsInItsAlignmentsAndBetweenBatches) {
    std::mt19937 rng(31);
    auto random_sequence = [&](size_t length) {
        std::string s(length, 'A');
        for (char &c : s) {
            c = "ACGT"[rng() % 4];
        }
        return s;
    };
    const size_t k = 7;
    std::vector<test::CoordinateRecord> records;
    for (size_t i = 0; i < 12; ++i) {
        records.push_back({ random_sequence(40 + rng() % 30), "column_" + std::to_string(i % 4),
                            "record_" + std::to_string(i) });
    }
    const size_t num_records = 24;
    std::vector<size_t> lengths;
    utils::TempFile fasta;
    for (size_t q = 0; q < num_records; ++q) {
        // a substitution, so that the alignment is not the exact match
        std::string sequence = records[q % 12].sequence + records[(q + 5) % 12].sequence;
        sequence[sequence.size() / 2] = sequence[sequence.size() / 2] == 'A' ? 'C' : 'A';
        fasta.ofstream() << ">q" << q << "\n" << sequence << "\n";
        lengths.push_back(sequence.size());
    }
    fasta.ofstream().close();
    // no record mapping and a mode other than coords: the batched path
    auto index = test::coordinate_test_graph(k, records, false, false, false, false);
    auto config = query_config("signature");
    config->discovery_fraction = 0;
    config->num_top_labels = 1000;
    const auto aligner = std::make_unique<graph::align::DBGAlignerConfig>(
            cli::initialize_aligner_config(*config, index->get_graph()));

    const unsigned int threads = get_num_threads();
    for (unsigned long long batch : { 100'000'000ull, 250ull }) {
        config->query_batch_size = batch;
        // the batches as batched_query_fasta reads them: records while the bytes read are
        // within the batch size
        size_t batches = 0;
        for (size_t q = 0; q < num_records; ++batches) {
            for (size_t bytes = 0; q < num_records && bytes <= batch; ++q) {
                bytes += lengths[q];
            }
        }
        if (batch > 1000) {
            ASSERT_EQ(1u, batches);
        } else {
            // two or three records a batch
            ASSERT_LE(num_records / 3, batches);
            ASSERT_GE(num_records / 2, batches);
        }
        for (unsigned int t : { 1u, 4u }) {
            SCOPED_TRACE(::testing::Message() << "batch " << batch << " threads " << t);
            set_num_threads(t);
            std::atomic<size_t> callbacks { 0 };
            const auto base = search_with_stop(fasta.name(), *config, *index, nullptr, nullptr,
                                               &callbacks, aligner.get());
            ASSERT_EQ(num_records, base.size());
            for (const auto &[id, result] : base) {
                ASSERT_TRUE(result.get_alignment()) << id;
            }

            // a stop that never fires: every record and its alignment, as without one; asked
            // when each batch is read, before each alignment, before each query graph and
            // before each record's query
            std::atomic<size_t> asked { 0 };
            cli::SearchStop never([&]() -> const char * {
                asked.fetch_add(1);
                return nullptr;
            });
            callbacks = 0;
            const auto all = search_with_stop(fasta.name(), *config, *index, nullptr, &never,
                                              &callbacks, aligner.get());
            ASSERT_EQ(num_records, all.size());
            for (const auto &[id, result] : base) {
                EXPECT_EQ(compact_json(result.to_json(false, k)),
                          compact_json(all.at(id).to_json(false, k))) << id;
                EXPECT_EQ(compact_json(result.cut_to_json()),
                          compact_json(all.at(id).cut_to_json())) << id;
                ASSERT_TRUE(all.at(id).get_alignment()) << id;
                EXPECT_EQ(result.get_alignment()->cigar, all.at(id).get_alignment()->cigar)
                        << id;
                EXPECT_EQ(result.get_alignment()->score, all.at(id).get_alignment()->score)
                        << id;
            }
            EXPECT_FALSE(never.stopped());
            EXPECT_EQ(num_records, never.records_started());
            EXPECT_EQ(2 * batches + 2 * num_records, asked.load());

            // a stop that fires at its |n|-th question
            std::vector<size_t> ns { 1, 2, 3, 2 * batches + 2 * num_records };
            if (batches == 1) {
                // one batch: question 1 reads it, 2 to num_records + 1 are its alignments,
                // num_records + 2 its query graph, then its records' queries. A stop anywhere
                // up to its first record starts none
                ns = { 2, 1 + num_records / 2, num_records + 1, num_records + 2,
                       num_records + 3 };
            }
            for (size_t n : ns) {
                SCOPED_TRACE(::testing::Message() << "stop at question " << n);
                std::atomic<size_t> questions { 0 };
                cli::SearchStop stop([&]() -> const char * {
                    return questions.fetch_add(1) + 1 >= n ? "the test's" : nullptr;
                });
                callbacks = 0;
                std::map<size_t, cli::SeqSearchResult> got;
                EXPECT_THROW(got = search_with_stop(fasta.name(), *config, *index, nullptr,
                                                    &stop, &callbacks, aligner.get()),
                             cli::SearchStopped);
                EXPECT_TRUE(stop.stopped());
                const size_t started = stop.records_started();
                EXPECT_EQ(started, callbacks.load());
                EXPECT_LT(started, num_records);
                EXPECT_LE(started, n - 1);
                if (batches == 1 || n == 1) {
                    EXPECT_EQ(0u, started);
                }
            }
        }
    }
    set_num_threads(threads);
}

} // namespace
