#include <algorithm>
#include <functional>
#include <map>
#include <random>
#include <tuple>
#include <vector>
#include <unordered_set>

#include <sdsl/int_vector.hpp>
#include "gtest/gtest.h"

#include "../graph/all/test_dbg_helpers.hpp"

#include "../annotation/test_matrix_helpers.hpp"

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/config/config.hpp"
#include "cli/query.hpp"
#include "common/seq_tools/reverse_complement.hpp"
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
    EXPECT_EQ(1u, json["returned"].asUInt64());
    EXPECT_FALSE(json["complete"].asBool());
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_EQ(2u, json["tied"].asUInt64());
    EXPECT_EQ(1u, json["tied_returned"].asUInt64());
    EXPECT_EQ(7u, json["nth_median"].asUInt64());
    EXPECT_EQ(1u, json["tied_same_median"].asUInt64());
    EXPECT_EQ(7u, json.size());

    // the N-th score of one label: no median read, none stated
    cut.tied = 1;
    cut.nth_median = 0;
    cut.nth_median_read = false;
    json = with_cut(cut).cut_to_json();
    EXPECT_EQ(12u, json["nth_score"].asUInt64());
    EXPECT_TRUE(json.isMember("nth_median"));
    EXPECT_TRUE(json["nth_median"].isNull());
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
            graph->map_to_nodes(query, [&](graph::DeBruijnGraph::node_index node) {
                if (!node)
                    return;
                std::string kmer = graph->get_node_sequence(node);
                if (uint64_t v = value(j, kmer))
                    values.push_back(v);
            });
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
    check_batched_parity<graph::DBGHashOrdered>(graph::DeBruijnGraph::CANONICAL, CountsForm::CSR);
    check_batched_parity<graph::DBGSuccinct>(graph::DeBruijnGraph::BASIC, CountsForm::CSR);
    check_batched_parity<graph::DBGHashOrdered>(graph::DeBruijnGraph::BASIC, CountsForm::INT_BRWT);
    check_batched_parity<graph::DBGSuccinct>(graph::DeBruijnGraph::CANONICAL, CountsForm::INT_BRWT);
}

// An int_brwt annotation's matrix that counts the integer values read from it
class CountingIntBRWT : public annot::matrix::CSCMatrix<annot::matrix::BRWT, annot::CountsVector> {
  public:
    using CSCMatrix::CSCMatrix;
    std::vector<RowValues> get_row_values(const std::vector<Row> &rows,
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
            const cli::SeqSearchResult tie = run(2);
            EXPECT_TRUE(tie.get_cut().nth_median_read);
            EXPECT_EQ(2u, tie.get_cut().tied);
            EXPECT_EQ(reads + 1, counted.reads);
            EXPECT_EQ(rows_read + 60, counted.rows_read);
        }
    }
}

} // namespace
