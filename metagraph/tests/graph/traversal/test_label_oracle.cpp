#include "gtest/gtest.h"

#include <algorithm>
#include <random>

#include "tests/test_helpers.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include "graph/traversal/label_oracle.hpp"
#include "graph/traversal/resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;

const size_t kTestK = 11;

// Fixed random sequences: long enough for every k-mer to be distinct with high
// probability, short enough for RowDiff fixtures to build quickly.
std::vector<std::string> make_sequences(size_t n, size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::vector<std::string> seqs;
    for (size_t i = 0; i < n; ++i) {
        std::string s(len, 'A');
        for (char &c : s) {
            c = "ACGT"[gen() % 4];
        }
        seqs.push_back(s);
    }
    return seqs;
}

#if ! _PROTEIN_GRAPH
std::string rc(std::string s) { ::reverse_complement(s); return s; }
#endif

std::vector<DeBruijnGraph::Mode> all_modes() {
    return {
        DeBruijnGraph::BASIC,
#if ! _PROTEIN_GRAPH
        DeBruijnGraph::CANONICAL,
        DeBruijnGraph::PRIMARY,
#endif
    };
}

template <typename Pair>
class LabelOracleTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGHashOrdered, annot::ColumnCompressed<>>,
    std::pair<DBGBitmap, annot::ColumnCompressed<>>,
    std::pair<DBGSSHash, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowFlatAnnotator>,
    std::pair<DBGHashFast, annot::RowFlatAnnotator>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>,
    std::pair<DBGHashFast, annot::RowDiffColumnAnnotator>,
    std::pair<DBGSuccinct, annot::RowDiffDiskAnnotator>
> OracleTypes;
TYPED_TEST_SUITE(LabelOracleTest, OracleTypes);

Regime expected_regime(const DeBruijnGraph &graph, DeBruijnGraph::Mode mode) {
    if (dynamic_cast<const CanonicalDBG*>(&graph))
        return Regime::PRIMARY;
    // DBGSSHash builds PRIMARY as native CANONICAL and is not wrapped
    if (graph.get_mode() == DeBruijnGraph::CANONICAL)
        return Regime::CANONICAL;
    EXPECT_EQ(DeBruijnGraph::BASIC, mode);
    return Regime::BASIC;
}

// regime detection and oriented-node -> annotation-key mapping
TYPED_TEST(LabelOracleTest, RegimeAndKeyMapping) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto seqs = make_sequences(3, 80, 1);
    std::vector<std::string> labels { "A", "B", "C" };

    for (DeBruijnGraph::Mode mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, mode);
        LabelOracle oracle(*anno);
        EXPECT_EQ(expected_regime(anno->get_graph(), mode), oracle.regime());
        EXPECT_EQ(kTestK, oracle.get_k());
        EXPECT_EQ(3u, oracle.num_columns());

        const DeBruijnGraph &graph = oracle.graph();
        for (size_t s = 0; s < seqs.size(); ++s) {
            const std::string &seq = seqs[s];
            auto keys = oracle.keys_of_sequence(seq);
            auto nodes = map_to_nodes_sequentially(graph, seq);
            ASSERT_EQ(seq.size() - kTestK + 1, keys.size());
            ASSERT_EQ(keys.size(), nodes.size());

            LabelQuery query(oracle, { oracle.resolve_label(labels[s]) });
            auto hits = query.fetch(keys);
            for (size_t i = 0; i < keys.size(); ++i) {
                ASSERT_NE(npos, nodes[i]);
                ASSERT_NE(npos, keys[i]);
                std::string_view kmer(seq.data() + i, kTestK);
                EXPECT_EQ(keys[i], oracle.key_of(nodes[i], kmer));
                // the label written on this k-mer is visible through the key
                EXPECT_EQ((LabelQuery::NodeHits{ 0 }), hits[i]) << "k-mer " << kmer;
            }

#if ! _PROTEIN_GRAPH
            // reverse-complement orientation
            std::string seq_rc = rc(seq);
            auto nodes_rc = map_to_nodes_sequentially(graph, seq_rc);
            auto keys_rc = oracle.keys_of_sequence(seq_rc);
            std::reverse(nodes_rc.begin(), nodes_rc.end());
            std::reverse(keys_rc.begin(), keys_rc.end());
            auto hits_rc = query.fetch(keys_rc);
            if (oracle.regime() == Regime::BASIC) {
                // forward strand only: an RC k-mer is in the graph only if it occurs on the
                // forward strand of some sequence, and then it carries that sequence's label
                for (size_t i = 0; i < keys_rc.size(); ++i) {
                    std::string rc_kmer = rc(seq.substr(i, kTestK));
                    bool on_forward_strand = false;
                    for (const auto &other : seqs) {
                        on_forward_strand |= other.find(rc_kmer) != std::string::npos;
                    }
                    EXPECT_EQ(on_forward_strand, keys_rc[i] != npos) << "rc k-mer " << rc_kmer;
                    bool in_this_seq = seq.find(rc_kmer) != std::string::npos;
                    EXPECT_EQ(in_this_seq, !hits_rc[i].empty()) << "rc k-mer " << rc_kmer;
                }
            } else {
                // both orientations map to the same annotation key and see the label
                for (size_t i = 0; i < keys.size(); ++i) {
                    ASSERT_NE(npos, nodes_rc[i]);
                    EXPECT_EQ(keys[i], keys_rc[i]);
                    std::string_view kmer_rc(seq_rc.data() + seq_rc.size() - kTestK - i, kTestK);
                    EXPECT_EQ(keys[i], oracle.key_of(nodes_rc[i], kmer_rc));
                    EXPECT_EQ((LabelQuery::NodeHits{ 0 }), hits_rc[i]);
                }
            }
#endif
        }
    }
}

// the hits agree with AnnotatedDBG::get_labels for every label subset
TYPED_TEST(LabelOracleTest, HitsAgreeWithGetLabels) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    // overlapping sequences so that k-mers carry several labels
    auto base = make_sequences(1, 120, 7)[0];
    std::vector<std::string> seqs { base, base.substr(20, 70), base.substr(50) + "ACGTACGTACGTAAACCC",
                                    make_sequences(1, 60, 9)[0] };
    std::vector<std::string> labels { "L0", "L1", "L2", "L3" };

    for (DeBruijnGraph::Mode mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, mode);
        LabelOracle oracle(*anno);

        std::vector<std::vector<std::string>> label_sets {
            { "L0" }, { "L1", "L3" }, { "L0", "L1", "L2", "L3" }
        };
        for (const auto &names : label_sets) {
            std::vector<LabelRef> refs;
            for (const auto &n : names) {
                refs.push_back(oracle.resolve_label(n));
            }
            LabelQuery query(oracle, refs);
            for (const auto &seq : seqs) {
                auto keys = oracle.keys_of_sequence(seq);
                auto hits = query.fetch(keys);
                ASSERT_EQ(keys.size(), hits.size());
                for (size_t i = 0; i < keys.size(); ++i) {
                    std::vector<std::string> expected;
                    for (const auto &l : anno->get_labels(seq.substr(i, kTestK), 1.0)) {
                        if (std::find(names.begin(), names.end(), l) != names.end())
                            expected.push_back(l);
                    }
                    std::vector<std::string> got;
                    for (LabelId h : hits[i]) {
                        got.push_back(names[h]);
                    }
                    std::sort(expected.begin(), expected.end());
                    std::sort(got.begin(), got.end());
                    EXPECT_EQ(expected, got) << "k-mer " << seq.substr(i, kTestK);
                }
            }
        }
        EXPECT_THROW(oracle.resolve_label("no_such_label"), std::invalid_argument);
        EXPECT_THROW(LabelQuery(oracle, {}), std::invalid_argument);
        EXPECT_THROW(LabelQuery(oracle, { oracle.resolve_label("L0"), oracle.resolve_label("L0") }),
                     std::invalid_argument);
    }
}

// Rows read with a RowDiffCache are the rows read without one: one key at a time, as a
// walk reads them, in both directions, and in batches. The bound is small, so that the
// cache drops generations.
TYPED_TEST(LabelOracleTest, RowDiffCacheGivesTheSameRows) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto seqs = make_sequences(3, 80, 5);
    auto anno = build_anno_graph<Graph, Annotation>(kTestK, seqs, { "A", "B", "C" });
    LabelOracle oracle(*anno);
    annot::matrix::RowDiffCache<annot::matrix::BinaryMatrix::SetBitPositions> cache(2048);
    for (const std::string &seq : seqs) {
        std::vector<Row> rows;
        for (node_index key : oracle.keys_of_sequence(seq)) {
            ASSERT_NE(npos, key);
            rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
        }
        for (Row row : rows) {
            EXPECT_EQ(oracle.get_rows({ row }), oracle.get_rows({ row }, &cache));
        }
        for (auto it = rows.rbegin(); it != rows.rend(); ++it) {
            EXPECT_EQ(oracle.get_rows({ *it }), oracle.get_rows({ *it }, &cache));
        }
        EXPECT_EQ(oracle.get_rows(rows), oracle.get_rows(rows, &cache));
    }
    if (oracle.is_row_diff()) {
        EXPECT_GT(cache.hits, 0u);
    } else {
        EXPECT_EQ(0u, cache.hits);
    }
}

// Unmasked DBGSuccinct (CLI-built graphs keep dummy k-mers): keys never point at dummies
TEST(LabelOracle, UnmaskedSuccinctKeys) {
    auto seqs = make_sequences(2, 60, 3);
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        auto graph = std::make_shared<DBGSuccinct>(kTestK, mode);
        for (const auto &s : seqs) {
            graph->add_sequence(s);
        }
        ASSERT_EQ(nullptr, graph->get_mask());
        auto anno = std::make_unique<AnnotatedDBG>(
            graph, std::make_unique<annot::ColumnCompressed<>>(graph->max_index()));
        anno->annotate_sequence(seqs[0], { "A" });
        anno->annotate_sequence(seqs[1], { "B" });

        LabelOracle oracle(*anno);
        EXPECT_EQ(mode == DeBruijnGraph::BASIC ? Regime::BASIC : Regime::CANONICAL, oracle.regime());
        LabelQuery query(oracle, { oracle.resolve_label("A"), oracle.resolve_label("B") });
        auto keys = oracle.keys_of_sequence(seqs[0]);
        auto hits = query.fetch(keys);
        for (size_t i = 0; i < keys.size(); ++i) {
            ASSERT_NE(npos, keys[i]);
            EXPECT_LE(keys[i], graph->max_index());
            ASSERT_GE(hits[i].size(), 1u);
            EXPECT_EQ(0u, hits[i][0]);
        }
    }
}

// A batch that mixes cached and uncached keys must survive eviction: the cache is cleared
// when it would overflow, so every key of the current call has to be refetched
TEST(LabelOracle, EvictionKeepsTheCurrentBatchAnswerable) {
    auto seqs = make_sequences(1, 60, 21);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kTestK, seqs, { "A" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    auto keys = oracle.keys_of_sequence(seqs[0]);
    ASSERT_GE(keys.size(), 3u);

    for (size_t capacity : { size_t(1), size_t(2), size_t(3) }) {
        LabelQuery query(oracle, { oracle.resolve_label("A") }, capacity);
        auto first = query.fetch({ keys[0], keys[1] });
        ASSERT_EQ(2u, first.size());
        // keys[0] is a hit, keys[2] a miss: the mixed batch forces an eviction
        auto second = query.fetch({ keys[0], keys[2] });
        ASSERT_EQ(2u, second.size());
        EXPECT_EQ(first[0], second[0]) << "capacity " << capacity;
        for (const auto &hits : second) {
            EXPECT_EQ((LabelQuery::NodeHits{ 0 }), hits);
        }
        // and a batch larger than the whole cache is still answered
        auto all = query.fetch(keys);
        ASSERT_EQ(keys.size(), all.size());
        for (const auto &hits : all) {
            EXPECT_EQ((LabelQuery::NodeHits{ 0 }), hits);
        }
    }
}


// Header labels via CoordToHeader (the refseq33m shape)
template <typename Pair>
class LabelOracleCoordTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>,
    std::pair<DBGHashFast, annot::RowDiffColumnAnnotator>
> CoordOracleTypes;
TYPED_TEST_SUITE(LabelOracleCoordTest, CoordOracleTypes);

TYPED_TEST(LabelOracleCoordTest, HeaderLabels) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    // one column "F" holding two sequences, like a FASTA file with two records,
    // plus a second column "G" with one sequence that shares a 30 bp block with s1
    auto rnd = make_sequences(3, 70, 11);
    std::string s1 = rnd[0], s2 = rnd[1], s3 = rnd[2].substr(0, 20) + s1.substr(10, 30) + rnd[2].substr(50);
    std::vector<std::string> seqs { s1, s2, s3 };
    std::vector<std::string> labels { "F", "F", "G" };
    uint64_t n1 = s1.size() - kTestK + 1, n2 = s2.size() - kTestK + 1, n3 = s3.size() - kTestK + 1;

    auto anno = build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, DeBruijnGraph::BASIC,
                                                    true, { 0, n1, 0 });
    ASSERT_TRUE(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&anno->get_annotator().get_matrix()));

    // CoordToHeader: column order is the label encoder order
    Column col_f = anno->get_annotator().get_label_encoder().encode("F");
    Column col_g = anno->get_annotator().get_label_encoder().encode("G");
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[col_f] = { "acc1", "acc2" }; num_kmers[col_f] = { n1, n2 };
    headers[col_g] = { "acc3" };         num_kmers[col_g] = { n3 };
    auto anno_h = with_coord_to_header(
            build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, DeBruijnGraph::BASIC,
                                                true, { 0, n1, 0 }),
            std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno_h);
    ASSERT_TRUE(oracle.has_coordinates());
    ASSERT_TRUE(oracle.coord_to_header());

    auto acc1 = oracle.resolve_label("acc1");
    auto acc2 = oracle.resolve_label("acc2");
    auto acc3 = oracle.resolve_label("acc3");
    auto col_label_f = oracle.resolve_label("F");
    EXPECT_EQ(LabelKind::HEADER, acc1.kind);
    EXPECT_EQ(col_f, acc1.column);
    EXPECT_EQ(0u, acc1.seq_id);
    EXPECT_EQ(1u, acc2.seq_id);
    EXPECT_EQ(LabelKind::COLUMN, col_label_f.kind);
    EXPECT_EQ(std::make_pair(n1, n1 + n2 - 1), oracle.sequence_coords(col_f, 1));
    EXPECT_EQ(std::make_pair(uint64_t(1), uint64_t(3)), oracle.map_coord(col_f, n1 + 3));

    LabelQuery query(oracle, { acc1, acc2, acc3, col_label_f });

    // s2's k-mers: header acc2 and column F, no acc1
    auto hits2 = query.fetch(oracle.keys_of_sequence(s2));
    for (size_t i = 0; i < hits2.size(); ++i) {
        EXPECT_EQ((LabelQuery::NodeHits{ 1, 3 }), hits2[i]) << "s2 k-mer " << i;
    }

    // s1's k-mers: acc1 everywhere; acc3 on the k-mers shared with s3
    auto hits1 = query.fetch(oracle.keys_of_sequence(s1));
    for (size_t i = 0; i < hits1.size(); ++i) {
        bool shared = s3.find(s1.substr(i, kTestK)) != std::string::npos;
        const LabelQuery::NodeHits expected = shared ? LabelQuery::NodeHits{ 0, 2, 3 }
                                                     : LabelQuery::NodeHits{ 0, 3 };
        EXPECT_EQ(expected, hits1[i]) << "s1 k-mer " << i;
    }

    // a name of both a column and a header, or the header of several records, is ambiguous
    std::vector<std::vector<std::string>> clash(2);
    std::vector<std::vector<uint64_t>> clash_kmers(2);
    clash[col_f] = { "G", "acc2" }; clash_kmers[col_f] = { n1, n2 };
    clash[col_g] = { "acc2" };      clash_kmers[col_g] = { n3 };
    auto clash_anno = with_coord_to_header(
            build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, DeBruijnGraph::BASIC,
                                                true, { 0, n1, 0 }),
            std::move(clash), std::move(clash_kmers));
    LabelOracle clash_oracle(*clash_anno);
    EXPECT_THROW(clash_oracle.resolve_label("G"), std::invalid_argument);
    EXPECT_THROW(clash_oracle.resolve_label("acc2"), std::invalid_argument);
    EXPECT_EQ(LabelKind::COLUMN, clash_oracle.resolve_label("F").kind);
    // discovery leaves those names out, so that every discovered label can be resubmitted
    ResolveOptions discover_clash;
    discover_clash.discover = true;
    discover_clash.discover_kind = LabelKind::HEADER;
    EXPECT_TRUE(resolve_support(clash_oracle, s1, discover_clash).labels.empty());
    discover_clash.discover_kind = LabelKind::COLUMN;
    auto columns = resolve_support(clash_oracle, s1, discover_clash).labels;
    ASSERT_EQ(1u, columns.size());
    EXPECT_EQ("F", columns[0].label.name);

    // a column without sequence records has no header labels: header discovery skips it
    std::vector<std::vector<std::string>> partial(2);
    std::vector<std::vector<uint64_t>> partial_kmers(2);
    partial[col_f] = { "acc1", "acc2" }; partial_kmers[col_f] = { n1, n2 };
    auto partial_anno = with_coord_to_header(
            build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, DeBruijnGraph::BASIC,
                                                true, { 0, n1, 0 }),
            std::move(partial), std::move(partial_kmers));
    LabelOracle partial_oracle(*partial_anno);
    EXPECT_TRUE(partial_oracle.has_header_records(col_f));
    EXPECT_FALSE(partial_oracle.has_header_records(col_g));
    ResolveOptions discover;
    discover.discover = true;
    discover.discover_kind = LabelKind::HEADER;
    auto profile = resolve_support(partial_oracle, s3, discover);
    ASSERT_FALSE(profile.labels.empty());
    for (const auto &support : profile.labels) {
        EXPECT_EQ(col_f, support.label.column);
    }

    // tuple rows read with a RowDiffCache are the rows read without one
    annot::matrix::RowDiffCache<annot::matrix::MultiIntMatrix::RowTuples> cache(2048);
    for (const std::string &seq : seqs) {
        std::vector<Row> rows;
        for (node_index key : oracle.keys_of_sequence(seq)) {
            rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
        }
        for (Row row : rows) {
            EXPECT_EQ(oracle.get_row_tuples({ row }), oracle.get_row_tuples({ row }, &cache));
        }
        for (auto it = rows.rbegin(); it != rows.rend(); ++it) {
            EXPECT_EQ(oracle.get_row_tuples({ *it }), oracle.get_row_tuples({ *it }, &cache));
        }
        EXPECT_EQ(oracle.get_row_tuples(rows), oracle.get_row_tuples(rows, &cache));
    }
    EXPECT_EQ(oracle.is_row_diff(), cache.hits > 0);

    // header labels need a CoordToHeader
    LabelOracle no_headers(*anno);
    EXPECT_THROW(no_headers.resolve_label("acc1"), std::invalid_argument);
    EXPECT_THROW(LabelQuery(no_headers, { acc1 }), std::invalid_argument);
}

// The reverse header index answers each header's first record, also after a load()
TEST(CoordToHeaderIndex, FindHeader) {
    for (const std::vector<std::string> &order : { std::vector<std::string>{ "A", "X", "Y", "Z" },
                                                   std::vector<std::string>{ "A", "Z", "Y", "X" } }) {
        auto anno = with_coord_to_header(
                build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                        5, { "ACGTAGATCGAA", "GTCTAACGTTGA", "TATGCCAGATCA", "AGCGTCTTGCGA" },
                        { "F", "F", "F", "F" }, DeBruijnGraph::BASIC, true, { 0, 8, 16, 24 }),
                { order }, { { 8, 8, 8, 8 } });
        LabelOracle oracle(*anno);
        for (size_t s = 0; s < order.size(); ++s) {
            auto ref = oracle.find_header(order[s]);
            ASSERT_TRUE(ref) << order[s];
            EXPECT_EQ(s, ref->seq_id) << order[s];
            EXPECT_EQ(order[s], oracle.header_name(ref->column, ref->seq_id));
            EXPECT_EQ(ref->seq_id, oracle.resolve_label(order[s]).seq_id);
        }
        EXPECT_FALSE(oracle.find_header("W"));
    }

    // a header repeated in a second column resolves to its first column, and is marked
    annot::CoordToHeader two({ { "A", "B" }, { "B", "C" } }, { { 1, 1 }, { 1, 1 } });
    two.index_headers();
    EXPECT_EQ(std::make_pair(Column(0), size_t(1)), *two.find_header("B"));
    EXPECT_EQ(std::make_pair(Column(1), size_t(1)), *two.find_header("C"));
    EXPECT_TRUE(two.is_repeated_header("B"));
    EXPECT_FALSE(two.is_repeated_header("C"));
    EXPECT_FALSE(two.is_repeated_header("W"));

    // load() replaces the headers and the index with them
    const std::string base = test_dump_dir() + "/coord_to_header_index";
    two.serialize(base);
    annot::CoordToHeader loaded({ { "Q" } }, { { 1 } });
    ASSERT_TRUE(loaded.find_header("Q"));
    ASSERT_TRUE(loaded.load(base));
    EXPECT_FALSE(loaded.find_header("Q"));
    EXPECT_EQ(std::make_pair(Column(1), size_t(1)), *loaded.find_header("C"));
}

} // namespace
