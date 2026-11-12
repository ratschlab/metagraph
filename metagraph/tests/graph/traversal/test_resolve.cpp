#include "gtest/gtest.h"

#include <algorithm>
#include <map>
#include <random>
#include <string>

#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/graph/traversal/traversal_test_helpers.hpp"

#include "graph/traversal/resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;

const size_t kK = 11;

std::vector<KmerInterval> iv(std::initializer_list<std::pair<uint64_t, uint64_t>> l) {
    std::vector<KmerInterval> out;
    for (auto [a, b] : l) {
        out.push_back({ a, b });
    }
    return out;
}

template <typename Pair>
class ResolveTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGSSHash, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowFlatAnnotator>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>
> ResolveTypes;
TYPED_TEST_SUITE(ResolveTest, ResolveTypes);

// Query = 200 bp with a 15 bp foreign insert at [80,95); label A is written on
// query[0,120) without the insert, B on query[100,200).
TYPED_TEST(ResolveTest, ThreeStatesAndRuns) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string q = random_seq(200, 42);
    std::string a_graph = q.substr(0, 80) + q.substr(95, 25);
    std::string b_seq = q.substr(100);

    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC
#if ! _PROTEIN_GRAPH
                                      , DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY
#endif
                                    }) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { a_graph, b_seq }, { "A", "B" }, mode);
        LabelOracle oracle(*anno);
        auto profile = resolve_support(oracle, q, { "A", "B" });
        EXPECT_EQ(190u, profile.num_kmers);
        EXPECT_EQ(kK, profile.k);

        // k-mers overlapping the insert [80,95) are absent from the graph: starts 70..94
        ASSERT_EQ(2u, profile.graph_runs.size()) << "mode " << mode;
        EXPECT_EQ((KmerInterval{ 0, 70 }), profile.graph_runs[0]);
        EXPECT_EQ((KmerInterval{ 95, 190 }), profile.graph_runs[1]);

        ASSERT_EQ(2u, profile.labels.size());
        const auto &A = profile.labels[0], &B = profile.labels[1];
        EXPECT_EQ("A", A.label.name);
        // A: present on [0,70) and on [95,110) (k-mers fully inside q[95,120))
        EXPECT_EQ(iv({ {0, 70}, {95, 110} }), A.runs);
        EXPECT_EQ(85u, A.kmers_supported);
        // B: k-mers fully inside q[100,200): starts 100..189
        EXPECT_EQ(iv({ {100, 190} }), B.runs);
        EXPECT_EQ(90u, B.kmers_supported);

        // agreement with the signature API for the labels present
        auto sigs = anno->get_top_label_signatures(q, 10, 0.0, 0.0);
        for (const auto &[label, count, mask] : sigs) {
            const auto &lp = label == "A" ? A : B;
            EXPECT_EQ(count, lp.kmers_supported);
            std::vector<bool> expected(mask.size());
            for (size_t i = 0; i < mask.size(); ++i) {
                expected[i] = mask[i];
            }
            std::vector<bool> got(profile.num_kmers, false);
            for (const auto &run : lp.runs) {
                for (uint64_t i = run.begin; i < run.end; ++i) {
                    got[i] = true;
                }
            }
            EXPECT_EQ(expected, got) << label;
        }
    }
}

// requests that are refused
TYPED_TEST(ResolveTest, Validation) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string q = random_seq(100, 5);
    auto anno = build_anno_graph<Graph, Annotation>(kK, { q, q.substr(0, 50) }, { "X", "W" },
                                                    DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    EXPECT_THROW(resolve_support(oracle, q, {}), std::invalid_argument);
    EXPECT_THROW(resolve_support(oracle, q, { "nope" }), std::invalid_argument);
    EXPECT_THROW(resolve_support(oracle, q.substr(0, kK - 1), { "X" }), std::invalid_argument);
}


template <typename Pair>
class ResolveCoordTest : public ::testing::Test {};
typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>
> ResolveCoordTypes;
TYPED_TEST_SUITE(ResolveCoordTest, ResolveCoordTypes);

// header labels: acc1 = L + R + M + R (repeat R twice) and acc2 = R, both in column F
TYPED_TEST(ResolveCoordTest, HeaderLabels) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string L = random_seq(40, 1), R = random_seq(40, 2), M = random_seq(40, 3);
    std::string acc1 = L + R + M + R;
    std::string acc2 = R;
    uint64_t n1 = acc1.size() - kK + 1, n2 = acc2.size() - kK + 1;
    auto anno = with_coord_to_header(
            build_anno_graph<Graph, Annotation>(kK, { acc1, acc2 }, { "F", "F" },
                                                DeBruijnGraph::BASIC, true, { 0, n1 }),
            { { "acc1", "acc2" } }, { { n1, n2 } });
    LabelOracle oracle(*anno);

    std::string q = R + M + R;
    auto profile = resolve_support(oracle, q, { "acc1", "acc2", "F" });
    ASSERT_EQ(3u, profile.labels.size());
    EXPECT_TRUE(profile.labels[0].label.kind == LabelKind::HEADER);
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), profile.labels[0].runs);     // acc1
    // acc2 = R: k-mers of the two R copies, i.e. [0,30) and [80,110)
    EXPECT_EQ(iv({ {0, 30}, {80, 110} }), profile.labels[1].runs);
    EXPECT_TRUE(profile.labels[2].label.kind == LabelKind::COLUMN);
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), profile.labels[2].runs);     // column F

    // an unknown header
    EXPECT_THROW(resolve_support(oracle, q, { "acc3" }), std::invalid_argument);
}

// The profiles of a few hundred named labels are those of the signature API, and a label
// absent from the query has none
TYPED_TEST(ResolveCoordTest, NamedLabelsGiveTheSignatureProfiles) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const std::string R = random_seq(40, 21), M = random_seq(40, 22);
    const std::string q = random_seq(150, 23) + R + M + R + M + random_seq(150, 24);
    std::mt19937 gen(25);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 300; ++i) {
        // a stretch of the query, some with a second stretch after it, some with a base changed
        const size_t a = gen() % (q.size() - 40);
        std::string s = q.substr(a, 20 + gen() % std::min<size_t>(200, q.size() - a - 20));
        if (i % 3 == 0) {
            const size_t c = gen() % (q.size() - 40);
            s += q.substr(c, 20 + gen() % 20);
        }
        if (i % 5 == 0)
            s[gen() % s.size()] = "ACGT"[gen() % 4];
        seqs.push_back(s);
        labels.push_back("L" + std::to_string(i));
    }
    seqs.push_back(random_seq(60, 26));
    labels.push_back("absent");
    auto anno = build_anno_graph<Graph, Annotation>(kK, seqs, labels, DeBruijnGraph::BASIC, true);
    LabelOracle oracle(*anno);

    const SupportProfile explicit_ = resolve_support(oracle, q, labels);
    ASSERT_EQ(labels.size(), explicit_.labels.size());
    std::map<std::string, const LabelProfile *> by_name;
    for (const LabelProfile &lp : explicit_.labels) {
        by_name[lp.label.name] = &lp;
    }
    const auto signatures = anno->get_top_label_signatures(q, labels.size(), 0.0, 0.0);
    ASSERT_GT(signatures.size(), 200u);
    for (const auto &[name, count, mask] : signatures) {
        const LabelProfile &e = *by_name.at(name);
        EXPECT_EQ(runs_of(mask.size(), [&](uint64_t i) { return mask[i]; }), e.runs) << name;
        EXPECT_EQ(count, e.kmers_supported) << name;
    }
    const LabelProfile &absent = *by_name.at("absent");
    EXPECT_TRUE(absent.runs.empty());
    EXPECT_EQ(0u, absent.kmers_supported);
}

} // namespace
