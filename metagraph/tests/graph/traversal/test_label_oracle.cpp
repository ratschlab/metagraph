#include "gtest/gtest.h"

#include <algorithm>
#include <random>

#include "tests/test_helpers.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include "graph/traversal/label_oracle.hpp"
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

std::string rc(std::string s) { ::reverse_complement(s); return s; }

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

// T2: regime detection and oriented-node -> annotation-key mapping
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
            auto path_keys = oracle.keys_of_path(nodes, seq);
            ASSERT_EQ(keys, path_keys);

            LabelQuery query(oracle, { oracle.resolve_label(labels[s]) }, false);
            for (size_t i = 0; i < keys.size(); ++i) {
                ASSERT_NE(npos, nodes[i]);
                ASSERT_NE(npos, keys[i]);
                std::string_view kmer(seq.data() + i, kTestK);
                EXPECT_EQ(keys[i], oracle.key_of(nodes[i], kmer));
                // the label written on this k-mer is visible through the key
                const auto &hits = query.fetch(keys[i]);
                ASSERT_EQ(1u, hits.size()) << "k-mer " << kmer;
                EXPECT_EQ(0u, hits[0].label);
            }

#if ! _PROTEIN_GRAPH
            // reverse-complement orientation
            std::string seq_rc = rc(seq);
            auto nodes_rc = map_to_nodes_sequentially(graph, seq_rc);
            auto keys_rc = oracle.keys_of_sequence(seq_rc);
            std::reverse(nodes_rc.begin(), nodes_rc.end());
            std::reverse(keys_rc.begin(), keys_rc.end());
            if (oracle.regime() == Regime::BASIC) {
                // forward strand only: an RC k-mer is in the graph only if it occurs
                // on the forward strand of some sequence (e.g. an inverted repeat),
                // and then it carries that sequence's label, not ours by orientation
                for (size_t i = 0; i < keys_rc.size(); ++i) {
                    std::string rc_kmer = rc(seq.substr(i, kTestK));
                    bool on_forward_strand = false;
                    for (const auto &other : seqs) {
                        on_forward_strand |= other.find(rc_kmer) != std::string::npos;
                    }
                    EXPECT_EQ(on_forward_strand, keys_rc[i] != npos) << "rc k-mer " << rc_kmer;
                    if (keys_rc[i] == npos)
                        continue;
                    const auto &hits = query.fetch(keys_rc[i]);
                    bool in_this_seq = seq.find(rc_kmer) != std::string::npos;
                    EXPECT_EQ(in_this_seq, !hits.empty()) << "rc k-mer " << rc_kmer;
                }
            } else {
                // both orientations map to the same annotation key and see the label
                for (size_t i = 0; i < keys.size(); ++i) {
                    ASSERT_NE(npos, nodes_rc[i]);
                    EXPECT_EQ(keys[i], keys_rc[i]);
                    std::string_view kmer_rc(seq_rc.data() + seq_rc.size() - kTestK - i, kTestK);
                    EXPECT_EQ(keys[i], oracle.key_of(nodes_rc[i], kmer_rc));
                    const auto &hits = query.fetch(keys_rc[i]);
                    ASSERT_EQ(1u, hits.size());
                    EXPECT_EQ(0u, hits[0].label);
                }
            }
#endif
        }
    }
}

// T1: direct-cell and batched-row predicates agree, and both agree with AnnotatedDBG
TYPED_TEST(LabelOracleTest, AccessPathsAgree) {
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
            LabelQuery rows(oracle, refs, false, LabelOracle::Access::ROWS);
            EXPECT_STREQ("rows", rows.access_path());

            std::unique_ptr<LabelQuery> direct;
            if (oracle.supports_direct()) {
                direct = std::make_unique<LabelQuery>(oracle, refs, false, LabelOracle::Access::DIRECT);
                EXPECT_STREQ("direct", direct->access_path());
            } else {
                EXPECT_THROW(LabelQuery(oracle, refs, false, LabelOracle::Access::DIRECT),
                             std::invalid_argument);
            }

            for (const auto &seq : seqs) {
                auto keys = oracle.keys_of_sequence(seq);
                auto hits_rows = rows.fetch(keys);
                ASSERT_EQ(keys.size(), hits_rows.size());
                if (direct) {
                    auto hits_direct = direct->fetch(keys);
                    EXPECT_EQ(hits_rows, hits_direct);
                }
                // ground truth: labels of each k-mer via the public API
                for (size_t i = 0; i < keys.size(); ++i) {
                    std::vector<std::string> expected;
                    for (const auto &l : anno->get_labels(seq.substr(i, kTestK), 1.0)) {
                        if (std::find(names.begin(), names.end(), l) != names.end())
                            expected.push_back(l);
                    }
                    std::vector<std::string> got;
                    for (const auto &h : hits_rows[i]) {
                        got.push_back(names[h.label]);
                    }
                    std::sort(expected.begin(), expected.end());
                    std::sort(got.begin(), got.end());
                    EXPECT_EQ(expected, got) << "k-mer " << seq.substr(i, kTestK);
                }
            }
        }
        // cache hits are counted
        EXPECT_GT(oracle.counters().rows_requested, 0u);
        EXPECT_THROW(oracle.resolve_label("no_such_label"), std::invalid_argument);
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
        EXPECT_NE(nullptr, oracle.dbg_succ());
        EXPECT_NE(nullptr, oracle.node_first_cache());
        LabelQuery query(oracle, { oracle.resolve_label("A"), oracle.resolve_label("B") }, false);
        auto keys = oracle.keys_of_sequence(seqs[0]);
        auto hits = query.fetch(keys);
        for (size_t i = 0; i < keys.size(); ++i) {
            ASSERT_NE(npos, keys[i]);
            EXPECT_LE(keys[i], graph->max_index());
            ASSERT_GE(hits[i].size(), 1u);
            EXPECT_EQ(0u, hits[i][0].label);
        }
    }
}


// A batch that mixes cached and uncached keys must survive eviction: the cache is
// cleared when it would overflow, so every key of the current call has to be refetched,
// not just the misses. Previously the second call below threw "Couldn't find key".
TEST(LabelOracle, EvictionKeepsTheCurrentBatchAnswerable) {
    auto seqs = make_sequences(1, 60, 21);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kTestK, seqs, { "A" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    auto keys = oracle.keys_of_sequence(seqs[0]);
    ASSERT_GE(keys.size(), 3u);

    for (size_t capacity : { size_t(1), size_t(2), size_t(3) }) {
        LabelQuery query(oracle, { oracle.resolve_label("A") }, false,
                         LabelOracle::Access::AUTO, capacity);
        auto first = query.fetch({ keys[0], keys[1] });
        ASSERT_EQ(2u, first.size());
        // keys[0] is a hit, keys[2] a miss: the mixed batch forces an eviction
        auto second = query.fetch({ keys[0], keys[2] });
        ASSERT_EQ(2u, second.size());
        EXPECT_EQ(first[0], second[0]) << "capacity " << capacity;
        for (const auto &hits : second) {
            ASSERT_EQ(1u, hits.size());
            EXPECT_EQ(0u, hits[0].label);
        }
        // and a batch larger than the whole cache is still answered
        auto all = query.fetch(keys);
        ASSERT_EQ(keys.size(), all.size());
        for (const auto &hits : all) {
            ASSERT_EQ(1u, hits.size());
        }
    }
}


// Coordinate-aware fixtures (the refseq33m shape): header labels via CoordToHeader
template <typename Pair>
class LabelOracleCoordTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>,
    std::pair<DBGHashFast, annot::RowDiffColumnAnnotator>
> CoordOracleTypes;
TYPED_TEST_SUITE(LabelOracleCoordTest, CoordOracleTypes);

TYPED_TEST(LabelOracleCoordTest, HeaderLabelsAndCoordinates) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    // one column "F" holding two sequences, like a FASTA file with two records,
    // plus a second column "G" with one sequence that shares a 30 bp block with s1
    auto rnd = make_sequences(3, 70, 11);
    std::string s1 = rnd[0], s2 = rnd[1], s3 = rnd[2].substr(0, 20) + s1.substr(10, 30) + rnd[2].substr(50);
    std::vector<std::string> seqs { s1, s2, s3 };
    std::vector<std::string> labels { "F", "F", "G" };
    uint64_t n1 = s1.size() - kTestK + 1, n2 = s2.size() - kTestK + 1, n3 = s3.size() - kTestK + 1;
    std::vector<uint64_t> coord_starts { 0, n1, 0 };

    auto anno = build_anno_graph<Graph, Annotation>(kTestK, seqs, labels, DeBruijnGraph::BASIC,
                                                    true, coord_starts);
    ASSERT_TRUE(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&anno->get_annotator().get_matrix()));

    // CoordToHeader: column order is the label encoder order (F then G)
    Column col_f = anno->get_annotator().get_label_encoder().encode("F");
    Column col_g = anno->get_annotator().get_label_encoder().encode("G");
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[col_f] = { "acc1", "acc2" }; num_kmers[col_f] = { n1, n2 };
    headers[col_g] = { "acc3" };         num_kmers[col_g] = { n3 };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno, &cth);
    ASSERT_TRUE(oracle.has_coordinates());
    ASSERT_EQ(&cth, oracle.coord_to_header());

    auto acc1 = oracle.resolve_label("acc1");
    auto acc2 = oracle.resolve_label("acc2");
    auto acc3 = oracle.resolve_label("acc3");
    auto col_label_f = oracle.resolve_label("F");
    EXPECT_EQ(LabelKind::HEADER, acc1.kind);
    EXPECT_EQ(col_f, acc1.column);
    EXPECT_EQ(0u, acc1.seq_id);
    EXPECT_EQ(1u, acc2.seq_id);
    EXPECT_EQ(LabelKind::COLUMN, col_label_f.kind);
    EXPECT_EQ(n2, oracle.num_kmers_in_sequence(col_f, 1));

    LabelQuery query(oracle, { acc1, acc2, acc3, col_label_f }, true);
    EXPECT_STREQ("tuples", query.access_path());
    EXPECT_THROW(LabelQuery(oracle, { acc1 }, false, LabelOracle::Access::DIRECT), std::invalid_argument);

    // s2's k-mers: header acc2 with local coords 0..n2-1, column F, no acc1
    auto keys2 = oracle.keys_of_sequence(s2);
    auto hits2 = query.fetch(keys2);
    for (size_t i = 0; i < keys2.size(); ++i) {
        ASSERT_NE(npos, keys2[i]);
        std::vector<LabelId> got;
        for (const auto &h : hits2[i]) got.push_back(h.label);
        ASSERT_EQ((std::vector<LabelId>{ 1, 3 }), got) << "s2 k-mer " << i;
        EXPECT_EQ((SmallVector<Coord>{ i }), hits2[i][0].coords);          // acc2 local coord
        EXPECT_EQ((SmallVector<Coord>{ n1 + i }), hits2[i][1].coords);     // column F global coord
    }

    // s1's k-mers: acc1 everywhere; acc3 on the shared block with s3's local coords
    auto keys1 = oracle.keys_of_sequence(s1);
    auto hits1 = query.fetch(keys1);
    for (size_t i = 0; i < keys1.size(); ++i) {
        std::vector<LabelId> got;
        for (const auto &h : hits1[i]) got.push_back(h.label);
        // shared with s3 iff the k-mer occurs in s3 (the planted block plus any
        // chance match at its boundaries); its acc3 coordinate is the position in s3
        size_t pos_in_s3 = s3.find(s1.substr(i, kTestK));
        bool shared = pos_in_s3 != std::string::npos;
        std::vector<LabelId> expected = shared ? std::vector<LabelId>{ 0, 2, 3 }
                                               : std::vector<LabelId>{ 0, 3 };
        ASSERT_EQ(expected, got) << "s1 k-mer " << i;
        EXPECT_EQ((SmallVector<Coord>{ i }), hits1[i][0].coords);
        if (shared)
            EXPECT_EQ((SmallVector<Coord>{ pos_in_s3 }), hits1[i][1].coords);
    }
    EXPECT_GT(oracle.counters().tuple_rows_fetched, 0u);
    EXPECT_GT(oracle.counters().coords_mapped, 0u);
}

} // namespace
