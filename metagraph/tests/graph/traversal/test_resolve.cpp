#include "gtest/gtest.h"

#include <algorithm>
#include <map>
#include <random>
#include <string>

#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"

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

std::string random_seq(size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::string s(len, 'A');
    for (char &c : s) {
        c = "ACGT"[gen() % 4];
    }
    return s;
}

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
        const bool canonical = mode != DeBruijnGraph::BASIC;

        ResolveOptions opts;
        opts.labels = { "A", "B" };
        auto profile = resolve_support(oracle, q, opts);
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

        // candidates: identical runs grouped, ordered longer first
        ASSERT_EQ(3u, profile.candidates.size());
        EXPECT_EQ((KmerInterval{ 100, 190 }), profile.candidates[0].kmers);
        EXPECT_EQ((std::vector<LabelId>{ 1 }), profile.candidates[0].labels);
        EXPECT_EQ((KmerInterval{ 0, 70 }), profile.candidates[1].kmers);
        EXPECT_EQ((KmerInterval{ 95, 110 }), profile.candidates[2].kmers);

        EXPECT_EQ("x70o25x15o80", encode_runs(A.runs, profile.num_kmers));

        // longest_first
        SelectionPolicy pol;
        pol.policy = SelectionPolicy::LONGEST_FIRST;
        pol.max_seeds = 2;
        pol.merge_overlapping = false;
        auto sel = select_seeds(profile, q, pol, canonical);
        EXPECT_EQ(3u, sel.num_candidates);
        EXPECT_EQ(3u, sel.num_eligible);
        ASSERT_EQ(2u, sel.seeds.size());
        EXPECT_EQ((KmerInterval{ 100, 190 }), sel.seeds[0].kmers);
        EXPECT_EQ(q.substr(100), sel.seeds[0].sequence);
        EXPECT_EQ((std::vector<std::string>{ "B" }), sel.seeds[0].labels);
        EXPECT_EQ((KmerInterval{ 0, 70 }), sel.seeds[1].kmers);
        EXPECT_EQ(q.substr(0, 80), sel.seeds[1].sequence);
        EXPECT_TRUE(sel.seeds[0].overlaps_with.empty());
        EXPECT_EQ(16u, sel.seeds[0].seed_id.size());
        EXPECT_NE(sel.seeds[0].seed_id, sel.seeds[1].seed_id);

        // max_support with overlap: A and B both cover [100,110)
        pol = SelectionPolicy();
        pol.max_seeds = 1;
        sel = select_seeds(profile, q, pol, canonical);
        ASSERT_EQ(1u, sel.seeds.size());
        EXPECT_EQ((KmerInterval{ 100, 110 }), sel.seeds[0].kmers);
        EXPECT_EQ((std::vector<std::string>{ "A", "B" }), sel.seeds[0].labels);
        EXPECT_EQ(2u, sel.seeds[0].population.supporting_total);

        // max_labels_per_seed drops deterministically and reports the drop
        pol.max_labels_per_seed = 1;
        sel = select_seeds(profile, q, pol, canonical);
        ASSERT_EQ(1u, sel.seeds.size());
        EXPECT_EQ(1u, sel.seeds[0].labels.size());
        EXPECT_EQ(1u, sel.seeds[0].population.dropped_count);
        EXPECT_EQ(1u, sel.seeds[0].population.dropped.size());
        EXPECT_FALSE(sel.seeds[0].population.dropped_digest.empty());
        EXPECT_EQ(sel.seeds[0].labels, select_seeds(profile, q, pol, canonical).seeds[0].labels);

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

// discovery with truncation keeps the top labels by k-mers, ties by column id
TYPED_TEST(ResolveTest, DiscoverTruncation) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string q = random_seq(100, 5);
    // X, Y, Z carry the whole query; W carries half
    auto anno = build_anno_graph<Graph, Annotation>(kK, { q, q, q, q.substr(0, 50) },
                                                    { "X", "Y", "Z", "W" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    ResolveOptions opts;
    opts.discover = true;
    opts.discover_max_labels = 2;
    auto profile = resolve_support(oracle, q, opts);
    ASSERT_TRUE(profile.labels_truncated.has_value());
    EXPECT_EQ(2u, profile.labels_truncated->kept);
    EXPECT_EQ(4u, profile.labels_truncated->total);
    EXPECT_EQ(90u, profile.labels_truncated->min_kept_kmers);
    EXPECT_EQ(90u, profile.labels_truncated->max_dropped_kmers);
    EXPECT_EQ(1u, profile.labels_truncated->dropped_full_length);  // Z dropped, W is not full-length
    ASSERT_EQ(2u, profile.labels.size());
    Column cx = anno->get_annotator().get_label_encoder().encode("X");
    Column cy = anno->get_annotator().get_label_encoder().encode("Y");
    EXPECT_EQ(std::min(cx, cy), profile.labels[0].label.column);
    EXPECT_EQ(std::max(cx, cy), profile.labels[1].label.column);

    opts.discover_max_labels = 10;
    profile = resolve_support(oracle, q, opts);
    EXPECT_FALSE(profile.labels_truncated.has_value());
    ASSERT_EQ(4u, profile.labels.size());
    EXPECT_EQ("W", profile.labels[3].label.name);
    EXPECT_EQ(40u, profile.labels[3].kmers_supported);

    // validation
    ResolveOptions bad;
    EXPECT_THROW(resolve_support(oracle, q, bad), std::invalid_argument);
    bad.labels = { "nope" };
    EXPECT_THROW(resolve_support(oracle, q, bad), std::invalid_argument);
    bad.discover = true;
    EXPECT_THROW(resolve_support(oracle, q, bad), std::invalid_argument);
    ResolveOptions headers;
    headers.discover = true;
    headers.discover_kind = LabelKind::HEADER;
    EXPECT_THROW(resolve_support(oracle, q, headers), std::invalid_argument);  // no CoordToHeader
    EXPECT_THROW(resolve_support(oracle, q.substr(0, kK - 1), opts), std::invalid_argument);
}

// max_support at many labels: a synthetic profile whose runs end at jittered positions. The
// breadth-optimal block is the prefix every label covers
TEST(Resolve, MaxSupportScalesWithManyLabels) {
    const size_t num_labels = 1000;
    const uint64_t num_kmers = 5000;
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = num_kmers;
    profile.graph_runs = iv({ {0, num_kmers} });
    std::mt19937 gen(7);
    for (size_t i = 0; i < num_labels; ++i) {
        LabelProfile lp;
        lp.label.column = static_cast<Column>(i);
        lp.label.name = "L" + std::to_string(i);
        // every label covers [0, 1000); each ends somewhere in [1000, 5000)
        uint64_t end = 1000 + gen() % (num_kmers - 1000);
        lp.runs = iv({ {0, end} });
        lp.kmers_supported = end;
        profile.labels.push_back(lp);
    }
    for (LabelId l = 0; l < profile.labels.size(); ++l) {
        profile.candidates.push_back({ profile.labels[l].runs[0], { l } });
    }

    SelectionPolicy policy;
    policy.max_seeds = 1;
    policy.min_block_bp = kK;   // min_kmers == 1
    const std::string query(num_kmers + kK - 1, 'A');
    auto selection = select_seeds(profile, query, policy, false);
    ASSERT_EQ(1u, selection.seeds.size());
    uint64_t smallest_end = num_kmers;
    for (const auto &lp : profile.labels) {
        smallest_end = std::min(smallest_end, lp.runs[0].end);
    }
    EXPECT_EQ((KmerInterval{ 0, smallest_end }), selection.seeds[0].kmers);
    EXPECT_EQ(num_labels, selection.seeds[0].labels.size());

    // a length floor trades labels for length
    policy.min_block_bp = 3000 + kK - 1;
    auto longer = select_seeds(profile, query, policy, false);
    ASSERT_EQ(1u, longer.seeds.size());
    EXPECT_GE(longer.seeds[0].kmers.size(), 3000u);
    EXPECT_LT(longer.seeds[0].labels.size(), num_labels);
    for (const auto &name : longer.seeds[0].labels) {
        const auto &lp = *std::find_if(profile.labels.begin(), profile.labels.end(),
                                       [&](const LabelProfile &x) { return x.label.name == name; });
        EXPECT_GE(lp.runs[0].end, longer.seeds[0].kmers.end);
    }
}

// Later max_support rounds see the runs CLIPPED by what earlier rounds took: with A on
// [0,10) and B on [0,5), the first round takes [0,5) (2 labels), the second [5,10)
TEST(Resolve, MaxSupportClipsRunsBetweenRounds) {
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = 10;
    profile.graph_runs = iv({ {0, 10} });
    for (auto [name, end] : { std::pair<const char*, uint64_t>{ "A", 10 },
                              std::pair<const char*, uint64_t>{ "B", 5 } }) {
        LabelProfile lp;
        lp.label.column = static_cast<Column>(profile.labels.size());
        lp.label.name = name;
        lp.runs = iv({ {0, end} });
        lp.kmers_supported = end;
        profile.labels.push_back(lp);
    }

    SelectionPolicy policy;
    policy.max_seeds = 2;
    policy.min_block_bp = kK;          // min_kmers == 1
    policy.merge_overlapping = false;  // keep the rounds' intervals as chosen
    auto sel = select_seeds(profile, std::string(10 + kK - 1, 'A'), policy, false);
    ASSERT_EQ(2u, sel.seeds.size());
    EXPECT_EQ((KmerInterval{ 0, 5 }), sel.seeds[0].kmers);
    EXPECT_EQ((std::vector<std::string>{ "A", "B" }), sel.seeds[0].labels);
    EXPECT_EQ((KmerInterval{ 5, 10 }), sel.seeds[1].kmers);
    EXPECT_EQ((std::vector<std::string>{ "A" }), sel.seeds[1].labels);
    EXPECT_TRUE(sel.seeds[0].overlaps_with.empty());
    EXPECT_TRUE(sel.seeds[1].overlaps_with.empty());
}

// seed ids: canonical orientation and label order do not change the id
TEST(Resolve, SeedId) {
    std::string s = "ACGTTGCAAGT";
    std::string rc = "ACTTGCAACGT";
    EXPECT_EQ(make_seed_id(s, true, { "b", "a" }), make_seed_id(rc, true, { "a", "b" }));
    EXPECT_NE(make_seed_id(s, false, { "a" }), make_seed_id(rc, false, { "a" }));
    // length-prefixed labels: {"a,b"} != {"a","b"}
    EXPECT_NE(make_seed_id(s, false, { "a,b" }), make_seed_id(s, false, { "a", "b" }));
    // the id of the labeled-traversal prototype for the same seed (no release name)
    EXPECT_EQ(hex64(fnv1a64("1:a", fnv1a64("\t", fnv1a64(s, fnv1a64("\t"))))),
              make_seed_id(s, false, { "a" }));
    EXPECT_EQ("x3o2x5", encode_runs(iv({ {0, 3}, {5, 10} }), 10));
}


template <typename Pair>
class ResolveCoordTest : public ::testing::Test {};
typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>
> ResolveCoordTypes;
TYPED_TEST_SUITE(ResolveCoordTest, ResolveCoordTypes);

// header labels: acc1 = L + R + M + R (repeat R twice) and acc2 = R, both in column F
TYPED_TEST(ResolveCoordTest, HeaderLabelsAndDiscovery) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string L = random_seq(40, 1), R = random_seq(40, 2), M = random_seq(40, 3);
    std::string acc1 = L + R + M + R;
    std::string acc2 = R;
    uint64_t n1 = acc1.size() - kK + 1, n2 = acc2.size() - kK + 1;
    auto anno = build_anno_graph<Graph, Annotation>(kK, { acc1, acc2 }, { "F", "F" },
                                                    DeBruijnGraph::BASIC, true, { 0, n1 });
    annot::CoordToHeader cth({ { "acc1", "acc2" } }, { { n1, n2 } });
    LabelOracle oracle(*anno, &cth);

    std::string q = R + M + R;
    ResolveOptions opts;
    opts.labels = { "acc1", "acc2", "F" };
    auto profile = resolve_support(oracle, q, opts);
    ASSERT_EQ(3u, profile.labels.size());
    EXPECT_TRUE(profile.labels[0].label.kind == LabelKind::HEADER);
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), profile.labels[0].runs);     // acc1
    // acc2 = R: k-mers of the two R copies, i.e. [0,30) and [80,110)
    EXPECT_EQ(iv({ {0, 30}, {80, 110} }), profile.labels[1].runs);
    EXPECT_TRUE(profile.labels[2].label.kind == LabelKind::COLUMN);
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), profile.labels[2].runs);     // column F

    // header discovery gives the same profiles
    ResolveOptions disc;
    disc.discover = true;
    disc.discover_kind = LabelKind::HEADER;
    auto discovered = resolve_support(oracle, q, disc);
    ASSERT_EQ(2u, discovered.labels.size());
    EXPECT_EQ("acc1", discovered.labels[0].label.name);
    EXPECT_EQ(profile.labels[0].runs, discovered.labels[0].runs);
    EXPECT_EQ("acc2", discovered.labels[1].label.name);
    EXPECT_EQ(60u, discovered.labels[1].kmers_supported);
    EXPECT_EQ(profile.labels[1].runs, discovered.labels[1].runs);

    // an unknown header
    opts.labels = { "acc3" };
    EXPECT_THROW(resolve_support(oracle, q, opts), std::invalid_argument);
}

// The explicit profiles of a few hundred labels are the discovery's, and a label absent from
// the query has none
TYPED_TEST(ResolveCoordTest, ExplicitLabelsGiveTheDiscoverysProfiles) {
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

    ResolveOptions disc;
    disc.discover = true;
    const SupportProfile discovered = resolve_support(oracle, q, disc);
    ResolveOptions opts;
    opts.labels = labels;
    const SupportProfile explicit_ = resolve_support(oracle, q, opts);
    ASSERT_EQ(labels.size(), explicit_.labels.size());
    std::map<std::string, const LabelProfile *> by_name;
    for (const LabelProfile &lp : explicit_.labels) {
        by_name[lp.label.name] = &lp;
    }
    ASSERT_GT(discovered.labels.size(), 200u);
    for (const LabelProfile &d : discovered.labels) {
        const LabelProfile &e = *by_name.at(d.label.name);
        EXPECT_EQ(d.runs, e.runs) << d.label.name;
        EXPECT_EQ(d.kmers_supported, e.kmers_supported) << d.label.name;
    }
    const LabelProfile &absent = *by_name.at("absent");
    EXPECT_TRUE(absent.runs.empty());
    EXPECT_EQ(0u, absent.kmers_supported);
    EXPECT_EQ(discovered.candidates.size(), explicit_.candidates.size());
}

} // namespace
