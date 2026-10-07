#include "gtest/gtest.h"

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "tests/test_helpers.hpp"
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
    for (char &c : s) c = "ACGT"[gen() % 4];
    return s;
}

std::vector<KmerInterval> iv(std::initializer_list<std::pair<uint64_t, uint64_t>> l) {
    std::vector<KmerInterval> out;
    for (auto [a, b] : l) out.push_back({ a, b });
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

// T3b: hand-built query. Query = 200 bp; a 15 bp foreign insert at [80,95);
// label A written on query[0,120) (with the insert), B on query[100,200).
TYPED_TEST(ResolveTest, ThreeStatesAndRuns) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string q = random_seq(200, 42);
    std::string a_seq = q.substr(0, 120);
    std::string b_seq = q.substr(100);
    // the insert region [80,95) is only in A's sequence; make the graph lack it by
    // writing A without it: A = q[0,80) + q[95,120)
    std::string a_graph = q.substr(0, 80) + q.substr(95, 25);

    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC
#if ! _PROTEIN_GRAPH
                                      , DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY
#endif
                                    }) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { a_graph, b_seq }, { "A", "B" }, mode);
        LabelOracle oracle(*anno);

        ResolveOptions opts;
        opts.labels = { "A", "B" };
        auto profile = resolve_support(oracle, q, opts);
        EXPECT_EQ(190u, profile.num_kmers);
        EXPECT_EQ(kK, profile.k);

        // k-mers overlapping the insert [80,95) are absent from the graph:
        // starts 70..94 (a k-mer starting at i covers [i, i+11))
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
        // in graph but without B: e.g. k-mer 0; absent from graph: k-mer 80
        EXPECT_FALSE(profile.graph_runs[0].begin > 0);

        // candidates: identical runs grouped, ordered longer first
        ASSERT_EQ(3u, profile.candidates.size());
        EXPECT_EQ((KmerInterval{ 100, 190 }), profile.candidates[0].kmers);
        EXPECT_EQ((std::vector<LabelId>{ 1 }), profile.candidates[0].labels);
        EXPECT_EQ((KmerInterval{ 0, 70 }), profile.candidates[1].kmers);
        EXPECT_EQ((KmerInterval{ 95, 110 }), profile.candidates[2].kmers);

        EXPECT_EQ("x70o25x15o80", encode_runs(A.runs, profile.num_kmers));

        // T4: disconnected blocks are separate candidates and an explicit gapped seed is rejected
        SelectionPolicy pol;
        pol.policy = SelectionPolicy::EXPLICIT;
        pol.explicit_seeds = { { { 60, 100 }, { "A" } } };
        EXPECT_THROW(select_seeds(profile, q, pol, mode != DeBruijnGraph::BASIC), std::invalid_argument);

        // longest_first
        pol = SelectionPolicy();
        pol.policy = SelectionPolicy::LONGEST_FIRST;
        pol.max_seeds = 2;
        pol.merge_overlapping = false;
        auto sel = select_seeds(profile, q, pol, mode != DeBruijnGraph::BASIC);
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
        pol.policy = SelectionPolicy::MAX_SUPPORT;
        pol.max_seeds = 1;
        sel = select_seeds(profile, q, pol, mode != DeBruijnGraph::BASIC);
        ASSERT_EQ(1u, sel.seeds.size());
        EXPECT_EQ((KmerInterval{ 100, 110 }), sel.seeds[0].kmers);
        EXPECT_EQ((std::vector<std::string>{ "A", "B" }), sel.seeds[0].labels);
        EXPECT_EQ(2u, sel.seeds[0].population.supporting_total);

        // explicit with a label that does not cover the interval
        pol = SelectionPolicy();
        pol.policy = SelectionPolicy::EXPLICIT;
        pol.explicit_seeds = { { { 120, 150 }, { "A", "B" } } };
        sel = select_seeds(profile, q, pol, mode != DeBruijnGraph::BASIC);
        ASSERT_EQ(1u, sel.seeds.size());
        EXPECT_EQ((std::vector<std::string>{ "B" }), sel.seeds[0].labels);
        EXPECT_EQ((std::vector<std::string>{ "A" }), sel.seeds[0].labels_not_covering);

        // max_labels_per_seed drops deterministically and reports the drop
        pol = SelectionPolicy();
        pol.policy = SelectionPolicy::MAX_SUPPORT;
        pol.max_seeds = 1;
        pol.max_labels_per_seed = 1;
        pol.label_order = SelectionPolicy::COLUMN_ID;
        sel = select_seeds(profile, q, pol, mode != DeBruijnGraph::BASIC);
        ASSERT_EQ(1u, sel.seeds.size());
        EXPECT_EQ(1u, sel.seeds[0].labels.size());
        EXPECT_EQ(1u, sel.seeds[0].population.dropped_count);
        EXPECT_EQ(1u, sel.seeds[0].population.dropped.size());
        EXPECT_FALSE(sel.seeds[0].population.dropped_digest.empty());

        // oracle agreement with the signature API for the labels present
        auto sigs = anno->get_top_label_signatures(q, 10, 0.0, 0.0);
        for (const auto &[label, count, mask] : sigs) {
            const auto &lp = label == "A" ? A : B;
            EXPECT_EQ(count, lp.kmers_supported);
            std::vector<bool> expected(mask.size());
            for (size_t i = 0; i < mask.size(); ++i) expected[i] = mask[i];
            std::vector<bool> got(profile.num_kmers, false);
            for (const auto &run : lp.runs) {
                for (uint64_t i = run.begin; i < run.end; ++i) got[i] = true;
            }
            EXPECT_EQ(expected, got) << label;
        }
    }
}

// T4c: discovery with truncation keeps the top labels by k-mers, ties by column id
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
    ResolveOptions trace;
    trace.labels = { "X" };
    trace.support = Support::TRACE;
    EXPECT_THROW(resolve_support(oracle, q, trace), std::invalid_argument);  // no coordinates
}

// max_support at the label counts this is meant for: a synthetic profile (no graph
// needed) with many labels whose runs end at jittered positions. The breadth-optimal
// block is the prefix every label covers; picking it must stay fast and must agree
// with a brute-force count (checked by an assert inside the sweep in debug builds).
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
        lp.label.kind = LabelKind::COLUMN;
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
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.max_seeds = 1;
    policy.min_block_bp = kK;   // min_kmers == 1
    auto selection = select_seeds(profile, std::string(num_kmers + kK - 1, 'A'), policy, false);
    ASSERT_EQ(1u, selection.seeds.size());
    // all 1000 labels cover [0, b) for b = the smallest run end
    uint64_t smallest_end = num_kmers;
    for (const auto &lp : profile.labels) {
        smallest_end = std::min(smallest_end, lp.runs[0].end);
    }
    EXPECT_EQ((KmerInterval{ 0, smallest_end }), selection.seeds[0].kmers);
    EXPECT_EQ(num_labels, selection.seeds[0].labels.size());

    // a length floor trades labels for length: the interval must be at least 3000
    // k-mers, so only the labels reaching that far can support it
    policy.min_block_bp = 3000 + kK - 1;
    auto longer = select_seeds(profile, std::string(num_kmers + kK - 1, 'A'), policy, false);
    ASSERT_EQ(1u, longer.seeds.size());
    EXPECT_GE(longer.seeds[0].kmers.size(), 3000u);
    EXPECT_LT(longer.seeds[0].labels.size(), num_labels);
    for (const auto &name : longer.seeds[0].labels) {
        const auto &lp = *std::find_if(profile.labels.begin(), profile.labels.end(),
                                       [&](const LabelProfile &x) { return x.label.name == name; });
        EXPECT_GE(lp.runs[0].end, longer.seeds[0].kmers.end);
    }
}

// Later max_support rounds must see the runs CLIPPED by what earlier rounds took.
// With A on [0,10) and B on [0,5), the first round takes the breadth-optimal [0,5)
// (2 labels); the second must then return [5,10), not a stub at the original run's end.
TEST(Resolve, MaxSupportClipsRunsBetweenRounds) {
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = 10;
    profile.graph_runs = iv({ {0, 10} });
    for (auto [name, end] : { std::pair<const char*, uint64_t>{ "A", 10 },
                              std::pair<const char*, uint64_t>{ "B", 5 } }) {
        LabelProfile lp;
        lp.label.kind = LabelKind::COLUMN;
        lp.label.column = static_cast<Column>(profile.labels.size());
        lp.label.name = name;
        lp.runs = iv({ {0, end} });
        lp.kmers_supported = end;
        profile.labels.push_back(lp);
    }

    SelectionPolicy policy;
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.max_seeds = 2;
    policy.min_block_bp = kK;          // min_kmers == 1
    policy.merge_overlapping = false;  // keep the rounds' intervals as chosen
    auto sel = select_seeds(profile, std::string(10 + kK - 1, 'A'), policy, false);
    ASSERT_EQ(2u, sel.seeds.size());
    EXPECT_EQ((KmerInterval{ 0, 5 }), sel.seeds[0].kmers);
    EXPECT_EQ((std::vector<std::string>{ "A", "B" }), sel.seeds[0].labels);
    EXPECT_EQ((KmerInterval{ 5, 10 }), sel.seeds[1].kmers);
    EXPECT_EQ((std::vector<std::string>{ "A" }), sel.seeds[1].labels);
    // the two seeds tile the supported region without overlapping
    EXPECT_TRUE(sel.seeds[0].overlaps_with.empty());
    EXPECT_TRUE(sel.seeds[1].overlaps_with.empty());
}

// seed ids: canonical orientation and label order do not change the id; release does
TEST(Resolve, SeedId) {
    std::string s = "ACGTTGCAAGT";
    std::string rc = "ACTTGCAACGT";
    EXPECT_EQ(make_seed_id("r1", s, true, { "b", "a" }), make_seed_id("r1", rc, true, { "a", "b" }));
    EXPECT_NE(make_seed_id("r1", s, false, { "a" }), make_seed_id("r1", rc, false, { "a" }));
    EXPECT_NE(make_seed_id("r1", s, true, { "a" }), make_seed_id("r2", s, true, { "a" }));
    // length-prefixed labels: {"a,b"} != {"a","b"}
    EXPECT_NE(make_seed_id("", s, false, { "a,b" }), make_seed_id("", s, false, { "a", "b" }));
    EXPECT_EQ("x3o2x5", encode_runs(iv({ {0, 3}, {5, 10} }), 10));
}


// T3/T29 (resolve part): trace-consistent runs on a coordinate index (refseq33m shape)
template <typename Pair>
class ResolveCoordTest : public ::testing::Test {};
typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>
> ResolveCoordTypes;
TYPED_TEST_SUITE(ResolveCoordTest, ResolveCoordTypes);

TYPED_TEST(ResolveCoordTest, TraceRunsAndHeaderDiscovery) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    // accession acc1 = L + R + M + R (repeat R twice); acc2 = R alone; one column "F"
    std::string L = random_seq(40, 1), R = random_seq(40, 2), M = random_seq(40, 3);
    std::string acc1 = L + R + M + R;
    std::string acc2 = R;
    uint64_t n1 = acc1.size() - kK + 1, n2 = acc2.size() - kK + 1;
    auto anno = build_anno_graph<Graph, Annotation>(kK, { acc1, acc2 }, { "F", "F" },
                                                    DeBruijnGraph::BASIC, true, { 0, n1 });
    std::vector<std::vector<std::string>> headers { { "acc1", "acc2" } };
    std::vector<std::vector<uint64_t>> num_kmers { { n1, n2 } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);

    // query = R + M + R: under kmer support acc1 covers everything; under trace support
    // the chain is consistent (second occurrence of R continues from M), while a query
    // M + R + L breaks the trace after R (R is followed by M or by the end in acc1).
    std::string q = R + M + R;
    ResolveOptions opts;
    opts.labels = { "acc1", "acc2", "F" };
    opts.support = Support::KMER;
    auto kmer_profile = resolve_support(oracle, q, opts);
    ASSERT_EQ(3u, kmer_profile.labels.size());
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), kmer_profile.labels[0].runs);     // acc1
    // acc2 = R: k-mers of the two R copies, i.e. [0,30) and [80,110)
    EXPECT_EQ(iv({ {0, 30}, {80, 110} }), kmer_profile.labels[1].runs);
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), kmer_profile.labels[2].runs);     // column F

    opts.support = Support::TRACE;
    auto trace_profile = resolve_support(oracle, q, opts);
    EXPECT_EQ(iv({ {0, q.size() - kK + 1} }), trace_profile.labels[0].runs);
    EXPECT_TRUE(trace_profile.labels[0].trace_breaks.empty());

    // M + R + L: the R->L junction k-mers [70,80) are absent from the graph, so both
    // support kinds give two runs separated by a graph gap, and no trace break
    std::string q2 = M + R + L;
    auto trace2 = resolve_support(oracle, q2, opts);
    EXPECT_EQ(iv({ {0, 70}, {80, 110} }), trace2.labels[0].runs);
    EXPECT_TRUE(trace2.labels[0].trace_breaks.empty());
    EXPECT_EQ(iv({ {0, 70}, {80, 110} }), trace2.graph_runs);
    opts.support = Support::KMER;
    EXPECT_EQ(iv({ {0, 70}, {80, 110} }), resolve_support(oracle, q2, opts).labels[0].runs);

    // R + M + R + M: every k-mer and junction exists (R->M after the first R, M->R
    // before the second), so kmer support is one run, but the trace through
    // acc1 = L R M R jumps back from coordinate 149 to 70 at the second R->M
    // junction: trace runs [0,110) and [110,150) with a break at 110
    std::string q3 = R + M + R + M;
    opts.support = Support::KMER;
    EXPECT_EQ(iv({ {0, 150} }), resolve_support(oracle, q3, opts).labels[0].runs);
    opts.support = Support::TRACE;
    auto trace3 = resolve_support(oracle, q3, opts);
    EXPECT_EQ(iv({ {0, 110}, {110, 150} }), trace3.labels[0].runs)
        << encode_runs(trace3.labels[0].runs, trace3.num_kmers);
    EXPECT_EQ((std::vector<uint64_t>{ 110 }), trace3.labels[0].trace_breaks);
    EXPECT_EQ(150u, trace3.labels[0].kmers_supported);

    // header discovery
    ResolveOptions disc;
    disc.discover = true;
    disc.discover_kind = LabelKind::HEADER;
    auto discovered = resolve_support(oracle, q, disc);
    ASSERT_EQ(2u, discovered.labels.size());
    EXPECT_EQ("acc1", discovered.labels[0].label.name);
    EXPECT_EQ("acc2", discovered.labels[1].label.name);
    EXPECT_EQ(60u, discovered.labels[1].kmers_supported);

    // trace is rejected on non-basic regimes
#if ! _PROTEIN_GRAPH
    auto canon = build_anno_graph<Graph, Annotation>(kK, { acc1 }, { "F" }, DeBruijnGraph::CANONICAL, true);
    LabelOracle oracle_c(*canon);
    ResolveOptions t;
    t.labels = { "F" };
    t.support = Support::TRACE;
    EXPECT_THROW(resolve_support(oracle_c, q, t), std::invalid_argument);
#endif
}

// The review of 2026-10-06, X-EFFICIENCY-04: explicit labels with trace scanned every k-mer's
// hit list once per label (O(labels x k-mers x hits): 2,000 labels took 5.9 s where the
// discovery returning the same profiles took 0.2 s). One pass now hands each hit to its label's
// accumulator, the calls the scan made in the same order: the explicit profiles of a few
// hundred labels — runs, trace breaks, k-mers supported, under kmer and trace — are the
// discovery's, and a label absent from the query has none
TYPED_TEST(ResolveCoordTest, ExplicitLabelsGiveTheDiscoverysProfiles) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    // the query holds a repeat (R M R M), so that traces jump where k-mer presence does not
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
    // a label nowhere in the query
    seqs.push_back(random_seq(60, 26));
    labels.push_back("absent");
    auto anno = build_anno_graph<Graph, Annotation>(kK, seqs, labels, DeBruijnGraph::BASIC, true);
    LabelOracle oracle(*anno);
    size_t breaks = 0;
    for (Support support : { Support::KMER, Support::TRACE }) {
        ResolveOptions disc;
        disc.discover = true;
        disc.discover_max_labels = 1000;
        disc.support = support;
        const SupportProfile discovered = resolve_support(oracle, q, disc);
        ResolveOptions opts;
        opts.labels = labels;
        opts.support = support;
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
            EXPECT_EQ(d.trace_breaks, e.trace_breaks) << d.label.name;
            EXPECT_EQ(d.kmers_supported, e.kmers_supported) << d.label.name;
            breaks += e.trace_breaks.size();
        }
        const LabelProfile &absent = *by_name.at("absent");
        EXPECT_TRUE(absent.runs.empty());
        EXPECT_EQ(0u, absent.kmers_supported);
        EXPECT_EQ(discovered.candidates.size(), explicit_.candidates.size());
    }
    // the trace jumped somewhere: the comparison covers trace breaks
    EXPECT_GT(breaks, 0u);
}


// ---- milestone 1b: the deadline of a /resolve (ResolveOptions::time_up)

// what a client reads of a profile, compared field by field
void expect_same_profile(const SupportProfile &expected, const SupportProfile &got,
                         const std::string &where) {
    EXPECT_EQ(expected.num_kmers, got.num_kmers) << where;
    EXPECT_EQ(expected.graph_runs, got.graph_runs) << where;
    ASSERT_EQ(expected.labels.size(), got.labels.size()) << where;
    for (size_t i = 0; i < expected.labels.size(); ++i) {
        const LabelProfile &e = expected.labels[i], &g = got.labels[i];
        EXPECT_EQ(e.label.name, g.label.name) << where;
        EXPECT_TRUE(e.label.kind == g.label.kind) << where;
        EXPECT_EQ(e.kmers_supported, g.kmers_supported) << where << " " << e.label.name;
        EXPECT_EQ(e.runs, g.runs) << where << " " << e.label.name;
        EXPECT_EQ(e.trace_breaks, g.trace_breaks) << where << " " << e.label.name;
    }
    ASSERT_EQ(expected.labels_truncated.has_value(), got.labels_truncated.has_value()) << where;
    if (expected.labels_truncated) {
        const LabelTruncation &e = *expected.labels_truncated, &g = *got.labels_truncated;
        EXPECT_EQ(std::make_tuple(e.kept, e.total, e.min_kept_kmers, e.max_dropped_kmers,
                                  e.dropped_full_length),
                  std::make_tuple(g.kept, g.total, g.min_kept_kmers, g.max_dropped_kmers,
                                  g.dropped_full_length)) << where;
    }
    ASSERT_EQ(expected.candidates.size(), got.candidates.size()) << where;
    for (size_t i = 0; i < expected.candidates.size(); ++i) {
        EXPECT_EQ(expected.candidates[i].kmers, got.candidates[i].kmers) << where;
        EXPECT_EQ(expected.candidates[i].labels, got.candidates[i].labels) << where;
    }
}

// Every stop |opts| can come to, one per poll of the work (the deadline passing at the n-th
// read, n = 1, 2, ... until the work completes before it), each checked against what the
// decision B7 promises: exactly the resolve of the query's first stop->resolved_kmers k-mers,
// resolved_kmers being a k-mer in the graph whose labels were not read, never decreasing in n.
// Without a stop the profile is the unbudgeted one. Returns the stops seen: (phase, k-mers
// resolved)
std::vector<std::pair<ResolveStop::Phase, uint64_t>>
check_every_stop(LabelOracle &oracle, const std::string &q, const ResolveOptions &opts,
                 const std::string &what) {
    const SupportProfile full = resolve_support(oracle, q, opts);
    EXPECT_FALSE(full.stop) << what;
    std::vector<std::pair<ResolveStop::Phase, uint64_t>> stops;
    uint64_t last = 0;
    for (size_t n = 1; ; ++n) {
        size_t polls = 0;
        ResolveOptions timed = opts;
        timed.time_up = [&polls, n]() { return ++polls >= n; };
        const SupportProfile got = resolve_support(oracle, q, timed);
        const std::string where = what + ", deadline at read " + std::to_string(n);
        if (!got.stop) {
            expect_same_profile(full, got, where);
            break;
        }
        const uint64_t x = got.stop->resolved_kmers;
        stops.emplace_back(got.stop->phase, x);
        EXPECT_EQ(x, got.num_kmers) << where;
        EXPECT_EQ(full.num_kmers, got.stop->query_kmers) << where;
        EXPECT_GE(x, last) << where;
        last = x;
        EXPECT_TRUE(std::any_of(full.graph_runs.begin(), full.graph_runs.end(),
                                [x](const KmerInterval &r) { return r.begin <= x && x < r.end; }))
            << where << ": the stop is at k-mer " << x << ", not one in the graph";
        if (x == 0) {
            EXPECT_TRUE(got.graph_runs.empty()) << where;
            EXPECT_TRUE(got.candidates.empty()) << where;
            EXPECT_EQ(opts.discover ? 0u : opts.labels.size(), got.labels.size()) << where;
            for (const LabelProfile &lp : got.labels) {
                EXPECT_EQ(0u, lp.kmers_supported) << where;
                EXPECT_TRUE(lp.runs.empty()) << where;
            }
        } else {
            ResolveOptions untimed = opts;
            expect_same_profile(resolve_support(oracle, q.substr(0, x + kK - 1), untimed), got,
                                where);
        }
        if (n > 100000) {
            ADD_FAILURE() << what << ": no end to the stops";
            break;
        }
    }
    return stops;
}

size_t count_phase(const std::vector<std::pair<ResolveStop::Phase, uint64_t>> &stops,
                   ResolveStop::Phase phase) {
    return std::count_if(stops.begin(), stops.end(),
                         [phase](const auto &s) { return s.first == phase; });
}

// Milestone 1b (decision B7, DESIGN-traverse-graphlet.md §21): a resolve its deadline stops is
// exactly the resolve of a query prefix, wherever the stop falls among the row batches — a
// discovery with and without truncation and explicit labels (primed from rows or tuple rows),
// presence and trace, on a query with a stretch absent from the graph and a repeat (rows kept
// for a later occurrence), its rows read three at a time
TYPED_TEST(ResolveCoordTest, DeadlineStopIsTheResolveOfAPrefix) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const std::string R = random_seq(40, 31), M = random_seq(40, 32);
    const std::string left = random_seq(150, 34) + R + M;
    const std::string right = R + M + random_seq(150, 35);
    // the k-mers touching the middle stretch are in no label's sequence
    const std::string q = left + random_seq(30, 33) + right;
    std::mt19937 gen(36);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 120; ++i) {
        const std::string &part = i % 2 ? left : right;
        const size_t a = gen() % (part.size() - 30);
        std::string s = part.substr(a, 20 + gen() % std::min<size_t>(150, part.size() - a - 20));
        if (i % 5 == 0)
            s[gen() % s.size()] = "ACGT"[gen() % 4];
        seqs.push_back(s);
        labels.push_back("L" + std::to_string(i));
    }
    auto anno = build_anno_graph<Graph, Annotation>(kK, seqs, labels, DeBruijnGraph::BASIC, true);
    LabelOracle oracle(*anno);
    size_t stops = 0;
    for (Support support : { Support::KMER, Support::TRACE }) {
        for (size_t max_labels : { size_t(0), size_t(1000), size_t(10) }) {
            ResolveOptions opts;
            opts.support = support;
            if (max_labels) {
                opts.discover = true;
                opts.discover_max_labels = max_labels;
            } else {
                opts.labels = labels;
            }
            opts.batch_rows = 3;
            const std::string what = std::string(support == Support::TRACE ? "trace" : "kmer")
                    + (max_labels ? ", discover " + std::to_string(max_labels) : ", explicit");
            const auto seen = check_every_stop(oracle, q, opts, what);
            // rows read three at a time: a stop between most of them
            EXPECT_GT(count_phase(seen, ResolveStop::ROWS), 30u) << what;
            stops += seen.size();
        }
    }
    EXPECT_GT(stops, 200u);
}

// Milestone 1b: the explicit labels' hits under a deadline. On the direct path (three labels of
// a column annotation: no rows primed, the hits are cell reads) they are fetched
// kResolveCheckKmers k-mers at a time with the deadline read before each piece: a stop there
// resolves the k-mers before the next one in the graph from the piece on. On the row paths each
// priming batch's k-mers are scattered as it is primed, so every stop falls between batches and
// keeps the rows read (review of 2026-10-07, T3-01/V1-01). Without a stop the profile is the one
// fetch's
TYPED_TEST(ResolveTest, DeadlineStopsTheExplicitHitsPass) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    // three and a half pieces of k-mers, a stretch of them absent from the graph across the
    // second piece's start
    const std::string q = random_seq(3 * kResolveCheckKmers + 2000, 41);
    const std::string a = q.substr(0, 4000), b = q.substr(4300, 9000);
    const std::string c = q.substr(2000, 1000) + q.substr(11000, 1000);
    auto anno = build_anno_graph<Graph, Annotation>(kK, { a, b, c }, { "A", "B", "C" },
                                                    DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    ResolveOptions opts;
    opts.labels = { "A", "B", "C" };
    // the row paths' rows in batches of 1,000: about 13 of them, a stop between each two
    opts.batch_rows = 1000;
    const auto seen = check_every_stop(oracle, q, opts, "explicit A, B, C");
    // The direct path's pieces begin at 4096, 8192 and 12288: a stop before each but the
    // first, the first at the next k-mer in the graph (about 4300: those of [3990, 4300) touch
    // the stretch no label has, but for an 11-mer of the graph met by chance). Its first read
    // (before any hit) stops in this phase, at k-mer 0
    ResolveOptions untimed = opts;
    const SupportProfile full = resolve_support(oracle, q, untimed);
    uint64_t next_in_graph = 0;
    for (const KmerInterval &run : full.graph_runs) {
        if (run.end > kResolveCheckKmers) {
            next_in_graph = std::max<uint64_t>(run.begin, kResolveCheckKmers);
            break;
        }
    }
    EXPECT_GT(next_in_graph, 4000u);
    EXPECT_LE(next_in_graph, 4300u);
    std::set<uint64_t> in_hits;
    uint64_t last_rows_stop = 0;
    for (const auto &[phase, x] : seen) {
        if (phase == ResolveStop::SUPPORT && x)
            in_hits.insert(x);
        if (phase == ResolveStop::ROWS)
            last_rows_stop = std::max(last_rows_stop, x);
    }
    // the path the annotation is read on (LabelQuery::access_path): a column annotation and
    // RowFlat read cells (direct), RowDiff decodes rows
    const bool row_path = count_phase(seen, ResolveStop::ROWS) > 0;
    if (std::is_same_v<Annotation, annot::ColumnCompressed<>>)
        EXPECT_FALSE(row_path);
    if (std::is_same_v<Annotation, annot::RowDiffColumnAnnotator>)
        EXPECT_TRUE(row_path);
    if (!row_path) {
        EXPECT_EQ((std::set<uint64_t>{ next_in_graph, 2 * kResolveCheckKmers,
                                       3 * kResolveCheckKmers }), in_hits);
    } else {
        // no stop in a hits pass, and the late ones keep their rows: before the fix every
        // priming stop past 4,096 k-mers became (support, about 4,300) at the hits pass's
        // second read of the deadline already passed
        EXPECT_TRUE(in_hits.empty());
        EXPECT_GT(count_phase(seen, ResolveStop::ROWS), 10u);
        EXPECT_GT(last_rows_stop, 2 * kResolveCheckKmers);
    }
    // a discovery reads the same rows in batches, never in the hits pass
    ResolveOptions disc;
    disc.discover = true;
    const auto discovered = check_every_stop(oracle, q, disc, "discover");
    EXPECT_EQ(0u, count_phase(discovered, ResolveStop::SUPPORT));
    EXPECT_GT(count_phase(discovered, ResolveStop::ROWS), 3u);
}

// Milestone 1b: the loops over the labels after the work read ResolveOptions::finish_check
// every kResolveCheckLabels labels, and what it throws abandons the request
TEST(Resolve, FinishCheckIsReadInTheLoopsOverTheLabels) {
    const std::string q = random_seq(200, 51);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 2 * kResolveCheckLabels + 10; ++i) {
        seqs.push_back(q.substr(i % 150, 40));
        labels.push_back("L" + std::to_string(i));
    }
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, seqs, labels);
    LabelOracle oracle(*anno);
    for (bool discover : { true, false }) {
        ResolveOptions opts;
        if (discover) {
            opts.discover = true;
            opts.discover_max_labels = labels.size();
        } else {
            opts.labels = labels;
        }
        size_t reads = 0;
        opts.finish_check = [&reads]() { reads++; };
        const SupportProfile profile = resolve_support(oracle, q, opts);
        EXPECT_EQ(labels.size(), profile.labels.size());
        EXPECT_GE(reads, 4u) << discover;
        opts.finish_check = []() { throw std::runtime_error("late"); };
        EXPECT_THROW(resolve_support(oracle, q, opts), std::runtime_error) << discover;
    }
}

} // namespace
