#include "gtest/gtest.h"

#include <algorithm>
#include <random>

#include "tests/test_helpers.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/graph/traversal/traversal_test_helpers.hpp"

#include "graph/traversal/label_oracle.hpp"
#include "graph/annotation_key_mapper.hpp"
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
            auto keys = map_to_nodes(oracle.graph(), seq);
            auto nodes = map_to_nodes_sequentially(graph, seq);
            ASSERT_EQ(seq.size() - kTestK + 1, keys.size());
            ASSERT_EQ(keys.size(), nodes.size());

            LabelQuery query(oracle, { oracle.resolve_label(labels[s]) });
            auto hits = query.fetch(keys);
            for (size_t i = 0; i < keys.size(); ++i) {
                ASSERT_NE(npos, nodes[i]);
                ASSERT_NE(npos, keys[i]);
                std::string_view kmer(seq.data() + i, kTestK);
                EXPECT_EQ(keys[i], AnnotationKeyMapper(oracle.graph()).key(nodes[i], kmer));
                // the label written on this k-mer is visible through the key
                EXPECT_EQ((LabelQuery::NodeHits{ 0 }), hits[i]) << "k-mer " << kmer;
            }

#if ! _PROTEIN_GRAPH
            // reverse-complement orientation
            std::string seq_rc = rc(seq);
            auto nodes_rc = map_to_nodes_sequentially(graph, seq_rc);
            auto keys_rc = map_to_nodes(oracle.graph(), seq_rc);
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
                    EXPECT_EQ(keys[i], AnnotationKeyMapper(oracle.graph()).key(nodes_rc[i], kmer_rc));
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
                auto keys = map_to_nodes(oracle.graph(), seq);
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
        for (node_index key : map_to_nodes(oracle.graph(), seq)) {
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
        auto keys = map_to_nodes(oracle.graph(), seqs[0]);
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
    auto keys = map_to_nodes(oracle.graph(), seqs[0]);
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
    EXPECT_EQ(std::make_pair(n1, n1 + n2 - 1), oracle.coord_to_header()->sequence_coords(col_f, 1));
    EXPECT_EQ(std::make_pair(uint64_t(1), uint64_t(3)), oracle.map_coord(col_f, n1 + 3));

    LabelQuery query(oracle, { acc1, acc2, acc3, col_label_f });

    // s2's k-mers: header acc2 and column F, no acc1
    auto hits2 = query.fetch(map_to_nodes(oracle.graph(), s2));
    for (size_t i = 0; i < hits2.size(); ++i) {
        EXPECT_EQ((LabelQuery::NodeHits{ 1, 3 }), hits2[i]) << "s2 k-mer " << i;
    }

    // the coordinates of a header label: s2's k-mer i is coordinate n1 + i of acc2's record
    {
        const auto with_coords = query.fetch_hits(map_to_nodes(oracle.graph(), s2));
        ASSERT_EQ(hits2.size(), with_coords.size());
        for (size_t i = 0; i < hits2.size(); ++i) {
            EXPECT_EQ(hits2[i], with_coords[i]->labels);
            const LabelQuery::HeaderCoord expected { 1, n1 + i };
            EXPECT_TRUE(std::find(with_coords[i]->coords.begin(), with_coords[i]->coords.end(),
                                  expected) != with_coords[i]->coords.end()) << "s2 k-mer " << i;
            EXPECT_TRUE(std::is_sorted(with_coords[i]->coords.begin(),
                                       with_coords[i]->coords.end()));
        }
        EXPECT_EQ(std::make_pair(n1, n1 + n2 - 1), query.record_coords(1));
        EXPECT_EQ(std::make_pair(uint64_t(0), uint64_t(0)), query.record_coords(3));
    }

    // s1's k-mers: acc1 everywhere; acc3 on the k-mers shared with s3
    auto hits1 = query.fetch(map_to_nodes(oracle.graph(), s1));
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
    EXPECT_TRUE(clash_oracle.is_ambiguous("G"));
    EXPECT_TRUE(clash_oracle.is_ambiguous("acc2"));
    EXPECT_FALSE(clash_oracle.is_ambiguous("F"));
    EXPECT_FALSE(clash_oracle.is_ambiguous("nope"));
    EXPECT_THROW(clash_oracle.resolve_label("G"), std::invalid_argument);
    EXPECT_THROW(clash_oracle.resolve_label("acc2"), std::invalid_argument);
    EXPECT_THROW(clash_oracle.resolve_label("nope"), std::invalid_argument);
    EXPECT_EQ(LabelKind::COLUMN, clash_oracle.resolve_label("F").kind);
    // discovery leaves those names out, so that every discovered label can be resubmitted
    ResolveOptions discover_clash;
    discover_clash.discover = true;
    discover_clash.discover_kind = LabelKind::HEADER;
    EXPECT_TRUE(resolve_support(clash_oracle, s1, discover_clash).labels.empty());
    // a cut list counts only the labels returned: here none
    discover_clash.discover_max_labels = 1;
    auto cut = resolve_support(clash_oracle, s1, discover_clash);
    EXPECT_TRUE(cut.labels.empty());
    ASSERT_TRUE(cut.labels_truncated.has_value());
    EXPECT_EQ(0u, cut.labels_truncated->kept);
    EXPECT_EQ(2u, cut.labels_truncated->total);
    discover_clash.discover_max_labels = ResolveOptions().discover_max_labels;
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
        for (node_index key : map_to_nodes(oracle.graph(), seq)) {
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

// The header coordinates of fetch_hits do not depend on the cache: a record with a tandem
// repeat (many coordinates per k-mer) read through a small cache, which is cleared by its
// coordinate bound and by its key bound, and in a batch too large to cache
TYPED_TEST(LabelOracleCoordTest, HeaderCoordinatesSurviveCacheEviction) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto rnd = make_sequences(3, 30, 23);
    std::string record = rnd[0];
    for (size_t i = 0; i < 40; ++i) {
        record += rnd[1].substr(0, 20);
    }
    record += rnd[2];
    const uint64_t n = record.size() - kTestK + 1;
    auto anno = with_coord_to_header(
            build_anno_graph<Graph, Annotation>(kTestK, { record }, { "F" },
                                                DeBruijnGraph::BASIC, true, { 0 }),
            { { "acc" } }, { { n } });
    LabelOracle oracle(*anno);
    const LabelRef acc = oracle.resolve_label("acc");
    ASSERT_EQ(LabelKind::HEADER, acc.kind);
    const std::vector<node_index> keys = map_to_nodes(oracle.graph(), record);
    ASSERT_EQ(n, keys.size());

    // the reference: one fetch through a cache that holds everything
    LabelQuery whole(oracle, { acc });
    std::vector<LabelQuery::Hits> expected;
    size_t max_coords = 0;
    size_t num_coords = 0;
    for (const LabelQuery::Hits *hit : whole.fetch_hits(keys)) {
        expected.push_back(*hit);
        max_coords = std::max(max_coords, hit->coords.size());
        num_coords += hit->coords.size();
    }
    // every k-mer has its own coordinate among those of its key
    for (size_t i = 0; i < keys.size(); ++i) {
        ASSERT_NE(npos, keys[i]);
        EXPECT_EQ((LabelQuery::NodeHits{ 0 }), expected[i].labels);
        const LabelQuery::HeaderCoord own { 0, i };
        EXPECT_TRUE(std::binary_search(expected[i].coords.begin(), expected[i].coords.end(), own))
            << "k-mer " << i;
    }
    // a k-mer of the tandem repeat has a coordinate per copy
    ASSERT_GE(max_coords, 30u);
    ASSERT_GT(num_coords, keys.size());

    auto same = [&](const LabelQuery::Hits &hit, size_t i) {
        return hit.labels == expected[i].labels && hit.coords == expected[i].coords;
    };
    // 8 keys and 8 * 16 coordinates: four k-mers of the repeat exceed the coordinate
    // bound, ten k-mers of a flank the key bound
    LabelQuery small(oracle, { acc }, 8);
    for (size_t round = 0; round < 2; ++round) {
        for (size_t begin = 0; begin < keys.size(); begin += 2) {
            const size_t end = std::min(keys.size(), begin + 2);
            const std::vector<node_index> batch(keys.begin() + begin, keys.begin() + end);
            const auto hits = small.fetch_hits(batch);
            ASSERT_EQ(batch.size(), hits.size());
            for (size_t i = begin; i < end; ++i) {
                EXPECT_TRUE(same(*hits[i - begin], i)) << "round " << round << " k-mer " << i;
            }
            // read again, now from the cache
            const auto again = small.fetch_hits(batch);
            for (size_t i = begin; i < end; ++i) {
                EXPECT_TRUE(same(*again[i - begin], i)) << "round " << round << " k-mer " << i;
            }
        }
        // a batch larger than the cache is answered without caching it
        const auto hits = small.fetch_hits(keys);
        ASSERT_EQ(keys.size(), hits.size());
        for (size_t i = 0; i < keys.size(); ++i) {
            EXPECT_TRUE(same(*hits[i], i)) << "round " << round << " k-mer " << i;
        }
        const auto labels = small.fetch(keys);
        for (size_t i = 0; i < keys.size(); ++i) {
            EXPECT_EQ(expected[i].labels, labels[i]);
        }
    }
    // no cache at all
    LabelQuery none(oracle, { acc }, 0);
    const auto hits = none.fetch_hits(keys);
    for (size_t i = 0; i < keys.size(); ++i) {
        EXPECT_TRUE(same(*hits[i], i)) << "k-mer " << i;
    }
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
