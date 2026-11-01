#include "gtest/gtest.h"

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <limits>
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
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/int_matrix/rank_extended/tuple_csc_matrix.hpp"
#include "common/vectors/bit_vector_sd.hpp"


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

// T4c: discovery with truncation keeps the top labels by k-mers, ties by name
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
    // X, Y and Z tie at the boundary, the first two by name kept
    EXPECT_EQ(3u, profile.labels_truncated->tied_at_boundary);
    EXPECT_EQ(2u, profile.labels_truncated->tied_kept);
    ASSERT_EQ(2u, profile.labels.size());
    EXPECT_EQ("X", profile.labels[0].label.name);
    EXPECT_EQ("Y", profile.labels[1].label.name);
    Column cx = anno->get_annotator().get_label_encoder().encode("X");
    Column cy = anno->get_annotator().get_label_encoder().encode("Y");
    EXPECT_EQ(cx, profile.labels[0].label.column);
    EXPECT_EQ(cy, profile.labels[1].label.column);

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
    ResolveOptions none;
    none.discover = true;
    none.discover_max_labels = 0;
    EXPECT_THROW(resolve_support(oracle, q, none), std::invalid_argument);
    ResolveOptions trace;
    trace.labels = { "X" };
    trace.support = Support::TRACE;
    EXPECT_THROW(resolve_support(oracle, q, trace), std::invalid_argument);  // no coordinates
}

// what a discovery kept, by name, in its order
std::vector<std::string> names_of(const SupportProfile &profile) {
    std::vector<std::string> names;
    for (const LabelProfile &lp : profile.labels) {
        names.push_back(lp.label.name);
    }
    return names;
}

// (total, kept, min_kept_kmers, max_dropped_kmers, dropped_full_length, tied_at_boundary,
// tied_kept) of a discovery's cut
using Cut = std::tuple<size_t, size_t, uint64_t, uint64_t, size_t, size_t, size_t>;
Cut cut_of(const SupportProfile &profile) {
    const LabelTruncation &t = profile.labels_truncated.value();
    return { t.total, t.kept, t.min_kept_kmers, t.max_dropped_kmers, t.dropped_full_length,
             t.tied_at_boundary, t.tied_kept };
}

// SPEC-multigraph-envelope.md §7 (E18): a discovery ranks by k-mers supported, ties by label
// name (bytes), never by the index's column numbering. The columns here are numbered in the
// reverse of their names' order within each support (c, b, a carry the whole query; y, x half
// of it each): the cut keeps the names first, states the tie group it went through
// (tied_at_boundary, tied_kept), and without a cut the equal-support labels are in name order
TYPED_TEST(ResolveTest, DiscoveryTiesByNameNotColumn) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const std::string q = random_seq(100, 7);
    auto anno = build_anno_graph<Graph, Annotation>(
            kK, { q.substr(0, 50), q, q.substr(50), q, q }, { "y", "c", "x", "b", "a" },
            DeBruijnGraph::BASIC);
    const auto &encoder = anno->get_annotator().get_label_encoder();
    // the premise: the columns' order is the reverse of the names' within each support
    ASSERT_LT(encoder.encode("c"), encoder.encode("b"));
    ASSERT_LT(encoder.encode("b"), encoder.encode("a"));
    ASSERT_LT(encoder.encode("y"), encoder.encode("x"));
    LabelOracle oracle(*anno);
    ResolveOptions opts;
    opts.discover = true;

    // no cut: more k-mers first, each support's labels by name
    opts.discover_max_labels = 10;
    SupportProfile profile = resolve_support(oracle, q, opts);
    EXPECT_FALSE(profile.labels_truncated.has_value());
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c", "x", "y" }), names_of(profile));
    ASSERT_EQ(5u, profile.labels.size());
    EXPECT_EQ(90u, profile.labels[2].kmers_supported);
    EXPECT_EQ(40u, profile.labels[3].kmers_supported);
    EXPECT_EQ(40u, profile.labels[4].kmers_supported);
    // exactly as many labels as found: no cut either
    opts.discover_max_labels = 5;
    EXPECT_EQ(names_of(profile), names_of(resolve_support(oracle, q, opts)));
    EXPECT_FALSE(resolve_support(oracle, q, opts).labels_truncated.has_value());

    // a cut through the full-length group keeps its first names: the column-first label (c)
    // is the one dropped first
    opts.discover_max_labels = 1;
    profile = resolve_support(oracle, q, opts);
    EXPECT_EQ((std::vector<std::string>{ "a" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 1, 90, 90, 2, 3, 1 }), cut_of(profile));

    opts.discover_max_labels = 2;
    profile = resolve_support(oracle, q, opts);
    EXPECT_EQ((std::vector<std::string>{ "a", "b" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 2, 90, 90, 1, 3, 2 }), cut_of(profile));

    // a cut ending with the group: all of it kept
    opts.discover_max_labels = 3;
    profile = resolve_support(oracle, q, opts);
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 3, 90, 40, 0, 3, 3 }), cut_of(profile));

    // a cut through the half-length group: x before y, the group before it whole and by name
    opts.discover_max_labels = 4;
    profile = resolve_support(oracle, q, opts);
    EXPECT_EQ((std::vector<std::string>{ "a", "b", "c", "x" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 4, 40, 40, 0, 2, 1 }), cut_of(profile));
    // the kept labels' columns are theirs
    for (const LabelProfile &lp : profile.labels) {
        EXPECT_EQ(encoder.encode(lp.label.name), lp.label.column) << lp.label.name;
    }
}

// max_support at the label counts this is meant for: a synthetic profile (no graph
// needed) with many labels whose runs end at jittered positions. The breadth-optimal
// block is the prefix every label covers; picking it must stay fast and must agree
// with a brute-force count (checked by an assert inside the sweep in debug builds).
// What a single-graph /resolve refuses before its admission of an explicit selection: the
// faults the request and k decide, in select_seeds' words, the fault select_seeds would state
// first also when that is a seed not fully in the graph; a request without such a fault passes,
// whatever select_seeds finds of it after the work
TEST(Resolve, ExplicitSeedsCheckedBeforeTheWork) {
    const std::string seq = random_seq(200, 7);
    auto anno_graph = build_anno_graph<DBGHashFast, annot::ColumnCompressed<>>(
            kK, { seq.substr(0, 100) }, { "A" });
    LabelOracle oracle(*anno_graph);
    // the first 90 k-mers are in the graph, then k - 1 + 100 are not
    const std::string query = seq.substr(0, 150);
    const uint64_t num_kmers = query.size() - kK + 1;
    ResolveOptions options;
    options.labels = { "A" };

    auto text = [&](std::vector<ExplicitSeed> seeds,
                    SelectionPolicy::Policy kind = SelectionPolicy::EXPLICIT,
                    const ResolveOptions &with = ResolveOptions()) -> std::string {
        SelectionPolicy policy;
        policy.policy = kind;
        policy.explicit_seeds = std::move(seeds);
        try {
            check_explicit_seeds(oracle, query,
                                 with.discover || with.drop_labels_not_in_graph ? with : options,
                                 policy);
        } catch (const std::invalid_argument &e) {
            return e.what();
        }
        return "";
    };
    const std::string range = "Explicit seed interval out of range";
    const std::string absent = "Explicit seed interval is not fully present in the graph";
    const std::string label = "Explicit seed label was not profiled: 'nope'";
    const KmerInterval good { 0, 10 };
    const KmerInterval missing { 80, 120 };

    // nothing certain: passed, a seed not fully in the graph left to the work
    EXPECT_EQ("", text({ { good, { "A" } } }));
    EXPECT_EQ("", text({ { missing, { "A" } } }));
    EXPECT_EQ("", text({ { good, { "A" } }, { missing, { "A" } } }));
    // one fault of the request's own
    EXPECT_EQ(range, text({ { { 0, num_kmers + 1 }, { "A" } } }));
    EXPECT_EQ(range, text({ { { 5, 5 }, { "A" } } }));
    EXPECT_EQ(range, text({ { { 6, 5 }, { "A" } } }));
    EXPECT_EQ(label, text({ { good, { "A", "nope" } } }));
    EXPECT_EQ(label, text({ { good, { "A" } }, { good, { "nope" } } }));
    EXPECT_EQ(range, text({ { good, { "A" } }, { { 5, 5 }, { "A" } } }));
    // two faults: select_seeds' first, seed by seed, the interval, its presence, its labels
    EXPECT_EQ(range, text({ { { 5, 5 }, { "nope" } } }));
    EXPECT_EQ(absent, text({ { missing, { "nope" } } }));
    EXPECT_EQ(absent, text({ { missing, { "A" } }, { { 5, 5 }, { "A" } } }));
    EXPECT_EQ(absent, text({ { missing, { "A" } }, { good, { "nope" } } }));
    EXPECT_EQ(range, text({ { { 5, 5 }, { "A" } }, { missing, { "nope" } } }));
    EXPECT_EQ(label, text({ { good, { "nope" } }, { missing, { "A" } } }));
    // not an explicit selection: its seeds are read by nothing
    EXPECT_EQ("", text({ { { 5, 5 }, { "A" } } }, SelectionPolicy::MAX_SUPPORT));
    // a discovery's labels are known after the work
    ResolveOptions discovery;
    discovery.discover = true;
    EXPECT_EQ("", text({ { good, { "nope" } } }, SelectionPolicy::EXPLICIT, discovery));
    EXPECT_EQ(range, text({ { { 5, 5 }, { "nope" } } }, SelectionPolicy::EXPLICIT, discovery));
    // so the scan ends at the first seed that names labels, after its interval: select_seeds
    // checks those labels before a later seed (recheck4: a later seed's interval was stated
    // where the base states the first seed's label). A fault before it is still the first
    const std::vector<ExplicitSeed> label_first { { good, { "nolabel" } },
                                                  { { 0, num_kmers + 1 }, {} } };
    EXPECT_EQ("", text(label_first, SelectionPolicy::EXPLICIT, discovery));
    EXPECT_EQ("", text({ { good, { "nope" } }, { { 5, 5 }, { "A" } } },
                       SelectionPolicy::EXPLICIT, discovery));
    EXPECT_EQ("", text({ { missing, { "nope" } }, { { 5, 5 }, { "A" } } },
                       SelectionPolicy::EXPLICIT, discovery));
    EXPECT_EQ(range, text({ { { 5, 5 }, { "A" } }, { good, { "nope" } } },
                          SelectionPolicy::EXPLICIT, discovery));
    EXPECT_EQ(range, text({ { good, {} }, { { 5, 5 }, { "nope" } } },
                          SelectionPolicy::EXPLICIT, discovery));
    EXPECT_EQ(absent, text({ { missing, {} }, { { 5, 5 }, { "nope" } } },
                           SelectionPolicy::EXPLICIT, discovery));
    // the same for a multi-graph pair's labels (a label it does not hold covers nothing)
    ResolveOptions dropping = options;
    dropping.drop_labels_not_in_graph = true;
    EXPECT_EQ("", text(label_first, SelectionPolicy::EXPLICIT, dropping));
    EXPECT_EQ(range, text({ { good, {} }, { { 5, 5 }, { "nope" } } },
                          SelectionPolicy::EXPLICIT, dropping));

    // what passes and what is stated are what select_seeds does on the profile
    auto select = [&](std::vector<ExplicitSeed> seeds,
                      const ResolveOptions &with = ResolveOptions()) -> std::string {
        SelectionPolicy policy;
        policy.policy = SelectionPolicy::EXPLICIT;
        policy.explicit_seeds = std::move(seeds);
        policy.release_id = "r";
        LabelOracle fresh(*anno_graph);
        try {
            select_seeds(resolve_support(fresh, query, with.discover ? with : options), query,
                         policy, false);
        } catch (const std::invalid_argument &e) {
            return e.what();
        }
        return "";
    };
    // the discovery's first fault is the first seed's label
    EXPECT_EQ("Explicit seed label was not profiled: 'nolabel'", select(label_first, discovery));
    const std::vector<std::vector<ExplicitSeed>> requests {
        { { { 5, 5 }, { "nope" } } },
        { { missing, { "nope" } } },
        { { missing, { "A" } }, { { 5, 5 }, { "A" } } },
        { { missing, { "A" } }, { good, { "nope" } } },
        { { good, { "nope" } }, { missing, { "A" } } },
        { { good, { "A" } }, { { 0, num_kmers + 1 }, { "A" } } },
    };
    for (size_t i = 0; i < requests.size(); ++i) {
        EXPECT_EQ(select(requests[i]), text(requests[i])) << i;
        EXPECT_NE("", text(requests[i])) << i;
    }
}

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

// A header discovery ties by the headers (SPEC-multigraph-envelope.md §7), never by (column,
// seq_id): four sequences carry the whole query, numbered in the reverse of their headers'
// order (column F holds hd, hc; column E holds hb, ha), and a fifth (aa, first by name) carries
// part of it, so ranks last whatever its name. Under presence and trace alike
TYPED_TEST(ResolveCoordTest, HeaderDiscoveryTiesByName) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const std::string q = random_seq(100, 8);
    const std::string part = q.substr(0, 50);
    const uint64_t n = q.size() - kK + 1, m = part.size() - kK + 1;
    auto anno = build_anno_graph<Graph, Annotation>(kK, { q, q, q, q, part },
                                                    { "F", "F", "E", "E", "E" },
                                                    DeBruijnGraph::BASIC, true,
                                                    { 0, n, 0, n, 2 * n });
    const auto &encoder = anno->get_annotator().get_label_encoder();
    ASSERT_EQ(0u, encoder.encode("F"));
    ASSERT_EQ(1u, encoder.encode("E"));
    std::vector<std::vector<std::string>> headers { { "hd", "hc" }, { "hb", "ha", "aa" } };
    std::vector<std::vector<uint64_t>> num_kmers { { n, n }, { n, n, m } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);
    for (Support support : { Support::KMER, Support::TRACE }) {
        const std::string what = support == Support::TRACE ? "trace" : "kmer";
        ResolveOptions opts;
        opts.discover = true;
        opts.discover_kind = LabelKind::HEADER;
        opts.support = support;

        opts.discover_max_labels = 10;
        SupportProfile profile = resolve_support(oracle, q, opts);
        EXPECT_FALSE(profile.labels_truncated.has_value()) << what;
        EXPECT_EQ((std::vector<std::string>{ "ha", "hb", "hc", "hd", "aa" }), names_of(profile))
            << what;
        ASSERT_EQ(5u, profile.labels.size()) << what;
        EXPECT_EQ(n, profile.labels[3].kmers_supported) << what;
        EXPECT_EQ(m, profile.labels[4].kmers_supported) << what;
        // each name's sequence: ha is (E, 1), hd (F, 0)
        const std::vector<std::pair<Column, uint64_t>> ids { { 1, 1 }, { 1, 0 }, { 0, 1 },
                                                             { 0, 0 }, { 1, 2 } };
        for (size_t i = 0; i < ids.size(); ++i) {
            EXPECT_EQ(LabelKind::HEADER, profile.labels[i].label.kind) << what;
            EXPECT_EQ(ids[i], std::make_pair(profile.labels[i].label.column,
                                             profile.labels[i].label.seq_id)) << what << " " << i;
        }

        // the cut through the four full-length sequences keeps the first headers
        opts.discover_max_labels = 1;
        profile = resolve_support(oracle, q, opts);
        EXPECT_EQ((std::vector<std::string>{ "ha" }), names_of(profile)) << what;
        EXPECT_EQ((Cut{ 5, 1, n, n, 3, 4, 1 }), cut_of(profile)) << what;

        opts.discover_max_labels = 3;
        profile = resolve_support(oracle, q, opts);
        EXPECT_EQ((std::vector<std::string>{ "ha", "hb", "hc" }), names_of(profile)) << what;
        EXPECT_EQ((Cut{ 5, 3, n, n, 1, 4, 3 }), cut_of(profile)) << what;

        opts.discover_max_labels = 4;
        profile = resolve_support(oracle, q, opts);
        EXPECT_EQ((std::vector<std::string>{ "ha", "hb", "hc", "hd" }), names_of(profile))
            << what;
        EXPECT_EQ((Cut{ 5, 4, n, m, 0, 4, 4 }), cut_of(profile)) << what;
    }
}

// SPEC-multigraph-envelope.md §7 under support "trace" (E18): a discovery ranks, cuts and
// states its cut (min_kept_kmers, max_dropped_kmers, dropped_full_length, tied_at_boundary,
// tied_kept) by the number it reports as kmers_supported -- under trace the traced k-mers
// (those where the label has a coordinate), under presence the k-mers it is present at. A
// coordinate matrix built by hand gives B every k-mer of the query but coordinates at its
// first 50 only (an empty tuple at the others), so that B has more k-mers present than A
// (q[0, 70): 60) and fewer traced; D and E carry the whole query with coordinates, C q[0, 40)
TEST(ResolveTraceDiscovery, RanksByTracedKmers) {
    const std::string q = random_seq(100, 21);
    auto graph = build_graph_batch<DBGSuccinct>(kK, { q }, DeBruijnGraph::BASIC);
    const uint64_t n = q.size() - kK + 1;
    // the annotation rows of the query's k-mers, in the query's order (each once)
    std::vector<uint64_t> rows;
    graph->map_to_nodes(q, [&](node_index i) {
        ASSERT_NE(DeBruijnGraph::npos, i);
        rows.push_back(AnnotatedDBG::graph_to_anno_index(i));
    });
    ASSERT_EQ(n, rows.size());
    ASSERT_EQ(n, std::set<uint64_t>(rows.begin(), rows.end()).size());
    const uint64_t num_rows = graph->max_index();
    // (label, k-mers present: the first ones, of them with a coordinate: the first ones)
    const std::vector<std::tuple<std::string, uint64_t, uint64_t>> labels {
        { "A", 60, 60 }, { "B", n, 50 }, { "C", 30, 30 }, { "D", n, n }, { "E", n, n } };
    std::vector<std::unique_ptr<bit_vector>> columns;
    std::vector<bit_vector_smart> delimiters;
    std::vector<sdsl::int_vector<>> values;
    annot::LabelEncoder<std::string> encoder;
    for (const auto &[name, present, traced] : labels) {
        encoder.insert_and_encode(name);
        // the label's rows in row order, each with its coordinate (consecutive along the
        // query: one trace run) or none
        std::vector<std::pair<uint64_t, std::vector<uint64_t>>> entries;
        for (uint64_t i = 0; i < present; ++i) {
            entries.emplace_back(rows[i], i < traced ? std::vector<uint64_t>{ 1000 + i }
                                                     : std::vector<uint64_t>{});
        }
        std::sort(entries.begin(), entries.end());
        sdsl::bit_vector column(num_rows, 0);
        std::vector<bool> delims { 1 };
        std::vector<uint64_t> vals;
        for (const auto &[row, coords] : entries) {
            column[row] = 1;
            for (uint64_t c : coords) {
                vals.push_back(c);
                delims.push_back(0);
            }
            delims.push_back(1);
        }
        columns.push_back(std::make_unique<bit_vector_sd>(column));
        sdsl::bit_vector d(delims.size());
        for (size_t i = 0; i < delims.size(); ++i) {
            d[i] = delims[i];
        }
        delimiters.emplace_back(d);
        sdsl::int_vector<> v(vals.size(), 0, 64);
        for (size_t i = 0; i < vals.size(); ++i) {
            v[i] = vals[i];
        }
        values.push_back(std::move(v));
    }
    auto matrix = std::make_unique<annot::matrix::TupleCSCMatrix<annot::matrix::ColumnMajor>>(
            annot::matrix::ColumnMajor(std::move(columns)), std::move(delimiters),
            std::move(values));
    AnnotatedDBG anno(graph, std::make_unique<annot::ColumnCoordAnnotator>(std::move(matrix),
                                                                           encoder));
    LabelOracle oracle(anno);
    ASSERT_TRUE(oracle.has_coordinates());

    auto discover = [&](Support support, size_t max_labels) {
        ResolveOptions opts;
        opts.discover = true;
        opts.support = support;
        opts.discover_max_labels = max_labels;
        return resolve_support(oracle, q, opts);
    };
    auto supported = [](const SupportProfile &profile) {
        std::vector<uint64_t> out;
        for (const LabelProfile &lp : profile.labels) {
            out.push_back(lp.kmers_supported);
        }
        return out;
    };

    // without a cut: the order of the k-mers each reports, ties by name
    SupportProfile all = discover(Support::TRACE, 10);
    EXPECT_FALSE(all.labels_truncated.has_value());
    EXPECT_EQ((std::vector<std::string>{ "D", "E", "A", "B", "C" }), names_of(all));
    EXPECT_EQ((std::vector<uint64_t>{ n, n, 60, 50, 30 }), supported(all));
    all = discover(Support::KMER, 10);
    EXPECT_EQ((std::vector<std::string>{ "B", "D", "E", "A", "C" }), names_of(all));
    EXPECT_EQ((std::vector<uint64_t>{ n, n, n, 60, 30 }), supported(all));

    // under trace B ranks by its 50 traced k-mers: below A, and never full length (traced at
    // every one of the query's present k-mers), which by presence it is
    SupportProfile profile = discover(Support::TRACE, 1);
    EXPECT_EQ((std::vector<std::string>{ "D" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 1, n, n, 1, 2, 1 }), cut_of(profile));
    profile = discover(Support::TRACE, 3);
    EXPECT_EQ((std::vector<std::string>{ "D", "E", "A" }), names_of(profile));
    EXPECT_EQ((std::vector<uint64_t>{ n, n, 60 }), supported(profile));
    EXPECT_EQ((Cut{ 5, 3, 60, 50, 0, 1, 1 }), cut_of(profile));
    profile = discover(Support::TRACE, 4);
    EXPECT_EQ((std::vector<std::string>{ "D", "E", "A", "B" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 4, 50, 30, 0, 1, 1 }), cut_of(profile));
    // by presence, unchanged: B is among the full-length labels
    profile = discover(Support::KMER, 1);
    EXPECT_EQ((std::vector<std::string>{ "B" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 1, n, n, 2, 3, 1 }), cut_of(profile));
    profile = discover(Support::KMER, 3);
    EXPECT_EQ((std::vector<std::string>{ "B", "D", "E" }), names_of(profile));
    EXPECT_EQ((Cut{ 5, 3, n, 60, 0, 3, 3 }), cut_of(profile));
}

// Explicit labels with trace hand each hit to its label's accumulator in one pass, the calls
// a per-label scan would make in the same order: scanning every k-mer's hit list once per
// label is O(labels x k-mers x hits) (2,000 labels: 5.9 s where the discovery returning the
// same profiles takes 0.2 s). The explicit profiles of a few hundred labels — runs, trace
// breaks, k-mers supported, under kmer and trace — are the discovery's, and a label absent
// from the query has none
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


// ---- the deadline of a /resolve (ResolveOptions::time_up)

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
                                  e.dropped_full_length, e.tied_at_boundary, e.tied_kept),
                  std::make_tuple(g.kept, g.total, g.min_kept_kmers, g.max_dropped_kmers,
                                  g.dropped_full_length, g.tied_at_boundary, g.tied_kept)) << where;
    }
    ASSERT_EQ(expected.candidates.size(), got.candidates.size()) << where;
    for (size_t i = 0; i < expected.candidates.size(); ++i) {
        EXPECT_EQ(expected.candidates[i].kmers, got.candidates[i].kmers) << where;
        EXPECT_EQ(expected.candidates[i].labels, got.candidates[i].labels) << where;
    }
}

// Every stop |opts| can come to, one per poll of the work (the deadline passing at the n-th
// read, n = 1, 2, ... until the work completes before it; past |dense| polls every tenth part
// of them, since each check re-runs the whole resolve and a prefix's), each checked against
// what the prefix rule promises: exactly the resolve of the query's first
// stop->resolved_kmers k-mers, resolved_kmers being a k-mer in the graph whose labels were not
// read, never decreasing in n. Without a stop the profile is the unbudgeted one. Returns the
// stops seen: (phase, k-mers resolved)
std::vector<std::pair<ResolveStop::Phase, uint64_t>>
check_every_stop(LabelOracle &oracle, const std::string &q, const ResolveOptions &opts,
                 const std::string &what, size_t dense = std::numeric_limits<size_t>::max()) {
    const SupportProfile full = resolve_support(oracle, q, opts);
    EXPECT_FALSE(full.stop) << what;
    std::vector<std::pair<ResolveStop::Phase, uint64_t>> stops;
    uint64_t last = 0;
    for (size_t n = 1; ; n += n < dense ? 1 : std::max<size_t>(1, n / 10)) {
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

// The prefix rule (DESIGN-traverse-graphlet.md §13): a resolve its deadline stops is exactly
// the resolve of a query prefix, wherever the stop falls among the row batches — a discovery
// with and without truncation and explicit labels (primed from rows or tuple rows), presence
// and trace, on a query with a stretch absent from the graph and a repeat (rows kept for a
// later occurrence), its rows read three at a time
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
    // every stop on the column annotation; on the row-diff one (whose reads cost the most)
    // the first 20 and then every tenth part
    const size_t dense = std::is_same_v<Annotation, annot::RowDiffColumnAnnotator>
            ? 20 : std::numeric_limits<size_t>::max();
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
            const auto seen = check_every_stop(oracle, q, opts, what, dense);
            // rows read three at a time: a stop between most of them
            EXPECT_GT(count_phase(seen, ResolveStop::ROWS), 30u) << what;
            stops += seen.size();
        }
    }
    EXPECT_GT(stops, 200u);
}

// The explicit labels' hits under a deadline. On the direct path (three labels of a column
// annotation: no rows primed, the hits are cell reads) they are fetched kResolveCheckKmers
// k-mers at a time with the deadline read before each piece: a stop there resolves the k-mers
// before the next one in the graph from the piece on. On the row paths each priming batch's
// k-mers are scattered as it is primed, so every stop falls between batches and keeps the rows
// read. Without a stop the profile is the one fetch's
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
    // (braces: GCC 13 -Werror=dangling-else on a gtest macro under an if)
    if (std::is_same_v<Annotation, annot::ColumnCompressed<>>) {
        EXPECT_FALSE(row_path);
    }
    if (std::is_same_v<Annotation, annot::RowDiffColumnAnnotator>) {
        EXPECT_TRUE(row_path);
    }
    if (!row_path) {
        EXPECT_EQ((std::set<uint64_t>{ next_in_graph, 2 * kResolveCheckKmers,
                                       3 * kResolveCheckKmers }), in_hits);
    } else {
        // no stop in a hits pass, and the late ones keep their rows: otherwise every priming
        // stop past 4,096 k-mers would become (support, about 4,300) at the hits pass's
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

// The loops over the labels after the work read ResolveOptions::finish_check every
// kResolveCheckLabels labels, and what it throws abandons the request. The reads are counted
// exactly (a lower bound would let one of a discovery's loops lose its check unseen): a loop
// over n labels reads it at the labels kResolveCheckLabels, 2 kResolveCheckLabels, ... below
// n; a discovery has five such loops (the ranking, the names of the labels the order needs —
// here every label, nothing is cut —, the kept labels, the profiles, the candidates'
// grouping), explicit labels two (the profiles, the grouping)
TEST(Resolve, FinishCheckIsReadInTheLoopsOverTheLabels) {
    const std::string q = random_seq(200, 51);
    for (size_t n : { kResolveCheckLabels, kResolveCheckLabels + 1,
                      2 * kResolveCheckLabels + 10 }) {
        std::vector<std::string> seqs, labels;
        for (size_t i = 0; i < n; ++i) {
            seqs.push_back(q.substr(i % 150, 40));
            labels.push_back("L" + std::to_string(i));
        }
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, seqs, labels);
        LabelOracle oracle(*anno);
        const size_t per_loop = (n - 1) / kResolveCheckLabels;
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
            EXPECT_EQ((discover ? 5 : 2) * per_loop, reads) << discover << " " << n;
            if (!per_loop)
                continue;
            opts.finish_check = []() { throw std::runtime_error("late"); };
            EXPECT_THROW(resolve_support(oracle, q, opts), std::runtime_error) << discover;
        }
    }
}

} // namespace
