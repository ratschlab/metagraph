#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "tests/test_helpers.hpp"
#include "tests/source_path.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include "graph/traversal/resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/traverse.hpp"
#include "common/utils/file_utils.hpp"


// The bounds of a /resolve discovery (SPEC-labeled-traversal-core.md §4.6, GPT review 5 #3):
// the label cap (ResolveOptions::discover_max_tracked) cuts the labels tracked by the first-met
// rule and states what the ranking then covers (SupportProfile::labels_capped); the memory
// bound (ResolveOptions::discover_max_memory) charges the pass before it allocates and stops it
// at a k-mer, the answer then exactly the resolve of the query's prefix before it. Neither
// changes a profile it does not bind.

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

// A query and records on it: |columns| columns, each with 1 to |per_column| records, each a
// stretch of the query (its left or right part: the 25 bp between are in no record), some with
// a base changed; the query repeats a stretch (R M ... R M), so that rows are kept for a later
// occurrence. Every record is laid out in its column's coordinates after the one before it,
// and the mapping names it "r<column>_<record>"
struct Records {
    std::string q;
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
    size_t num_records = 0;
};

template <class Graph, class Annotation>
Records make_records(size_t columns, size_t per_column, uint32_t seed) {
    Records r;
    std::mt19937 gen(seed);
    const std::string R = random_seq(40, seed + 1), M = random_seq(40, seed + 2);
    const std::string left = random_seq(120, seed + 3) + R + M;
    const std::string right = R + M + random_seq(120, seed + 4);
    r.q = left + random_seq(25, seed + 5) + right;
    std::vector<std::string> seqs, labels;
    std::vector<uint64_t> starts;
    std::vector<std::vector<std::string>> headers(columns);
    std::vector<std::vector<uint64_t>> num_kmers(columns);
    for (size_t c = 0; c < columns; ++c) {
        const size_t n = 1 + gen() % per_column;
        uint64_t offset = 0;
        for (size_t i = 0; i < n; ++i) {
            const std::string &part = gen() % 2 ? left : right;
            const size_t a = gen() % (part.size() - 30);
            std::string s = part.substr(a, 20 + gen() % std::min<size_t>(150, part.size() - a - 20));
            if (gen() % 4 == 0)
                s[gen() % s.size()] = "ACGT"[gen() % 4];
            seqs.push_back(s);
            labels.push_back("C" + std::to_string(c));
            starts.push_back(offset);
            headers[c].push_back("r" + std::to_string(c) + "_" + std::to_string(i));
            num_kmers[c].push_back(s.size() - kK + 1);
            offset += s.size() - kK + 1;
        }
        r.num_records += n;
    }
    r.anno = build_anno_graph<Graph, Annotation>(kK, seqs, labels, DeBruijnGraph::BASIC, true,
                                                 starts);
    r.cth = std::make_unique<annot::CoordToHeader>(std::move(headers), std::move(num_kmers));
    return r;
}

// what a client reads of a profile, compared field by field (the cap's statement apart)
void expect_same_profile(const SupportProfile &expected, const SupportProfile &got,
                         const std::string &where) {
    EXPECT_EQ(expected.num_kmers, got.num_kmers) << where;
    EXPECT_EQ(expected.graph_runs, got.graph_runs) << where;
    ASSERT_EQ(expected.labels.size(), got.labels.size()) << where;
    for (size_t i = 0; i < expected.labels.size(); ++i) {
        const LabelProfile &e = expected.labels[i], &g = got.labels[i];
        EXPECT_EQ(e.label.name, g.label.name) << where;
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

void expect_same_cap(const std::optional<LabelCap> &e, const std::optional<LabelCap> &g,
                     const std::string &where) {
    ASSERT_EQ(e.has_value(), g.has_value()) << where;
    if (!e)
        return;
    EXPECT_EQ(std::make_tuple(e->max_tracked, e->tracked, e->cut_kmer, e->new_at_cut,
                              e->untracked_max_kmers, e->ranking_exact, e->kept,
                              e->min_kept_kmers),
              std::make_tuple(g->max_tracked, g->tracked, g->cut_kmer, g->new_at_cut,
                              g->untracked_max_kmers, g->ranking_exact, g->kept,
                              g->min_kept_kmers)) << where;
}

// ranking_exact is checkable from the answer alone (the suite's condition, §4.6): true exactly
// when the kept list has max_labels entries and its last supports more than untracked_max_kmers
void expect_ranking_exact_checkable(const SupportProfile &p, size_t max_labels,
                                    const std::string &where) {
    if (!p.labels_capped)
        return;
    const bool from_answer = p.labels.size() == max_labels
            && p.labels.back().kmers_supported > p.labels_capped->untracked_max_kmers;
    EXPECT_EQ(from_answer, p.labels_capped->ranking_exact) << where;
}

std::vector<std::string> names_of(const SupportProfile &profile) {
    std::vector<std::string> names;
    for (const LabelProfile &lp : profile.labels) {
        names.push_back(lp.label.name);
    }
    return names;
}

uint64_t in_graph_from(const SupportProfile &profile, uint64_t from) {
    uint64_t n = 0;
    for (const KmerInterval &run : profile.graph_runs) {
        if (run.end > from)
            n += run.end - std::max(run.begin, from);
    }
    return n;
}

constexpr size_t kAll = std::numeric_limits<size_t>::max();

std::vector<KmerInterval> iv_of(uint64_t begin, uint64_t end) {
    return { KmerInterval { begin, end } };
}

template <typename Pair>
class ResolveBoundsTest : public ::testing::Test {};
typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>
> ResolveBoundsTypes;
TYPED_TEST_SUITE(ResolveBoundsTest, ResolveBoundsTypes);

struct Case {
    LabelKind kind;
    Support support;
    std::string name() const {
        return std::string(kind == LabelKind::HEADER ? "header" : "column")
                + (support == Support::TRACE ? " trace" : " kmer");
    }
};
const std::vector<Case> kCases {
    { LabelKind::COLUMN, Support::KMER }, { LabelKind::COLUMN, Support::TRACE },
    { LabelKind::HEADER, Support::KMER }, { LabelKind::HEADER, Support::TRACE },
};

ResolveOptions discovery(const Case &c, size_t max_labels) {
    ResolveOptions opts;
    opts.discover = true;
    opts.discover_kind = c.kind;
    opts.support = c.support;
    opts.discover_max_labels = max_labels;
    return opts;
}

// A cap at or above the labels met and a memory bound with room: the profile is the unbounded
// one, its modelled memory included, and nothing is stated
TYPED_TEST(ResolveBoundsTest, BoundsThatDoNotBindChangeNothing) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 101);
    LabelOracle oracle(*r.anno, r.cth.get());
    for (const Case &c : kCases) {
        for (size_t max_labels : { kAll, size_t(5) }) {
            const std::string what = c.name() + " max_labels " + std::to_string(max_labels);
            const SupportProfile plain = resolve_support(oracle, r.q, discovery(c, max_labels));
            ASSERT_FALSE(plain.labels_capped) << what;
            ASSERT_FALSE(plain.stop) << what;
            const size_t met = plain.labels_truncated ? plain.labels_truncated->total
                                                      : plain.labels.size();
            ASSERT_GT(met, 20u) << what;
            ResolveOptions bounded = discovery(c, max_labels);
            bounded.discover_max_tracked = met;
            bounded.discover_max_memory = uint64_t(1) << 30;
            const SupportProfile got = resolve_support(oracle, r.q, bounded);
            expect_same_profile(plain, got, what);
            EXPECT_FALSE(got.labels_capped) << what;
            EXPECT_FALSE(got.stop) << what;
            EXPECT_EQ(plain.memory_peak, got.memory_peak) << what;
            // one label fewer binds
            bounded.discover_max_tracked = met - 1;
            const SupportProfile binds = resolve_support(oracle, r.q, bounded);
            EXPECT_TRUE(binds.labels_capped) << what;
            expect_ranking_exact_checkable(binds, max_labels, what);
        }
    }
}

// The first-met rule, for every cap: the labels tracked are exactly those first met at a k-mer
// before the cut k-mer, the first k-mer whose new labels would pass the cap; each tracked label's
// profile is the uncut one; an untracked label supports at most the in-graph k-mers from the
// cut k-mer on; the kept labels above that bound are the uncut discovery's in its order, and
// the kept list is the uncut top max_labels exactly when ranking_exact says so
TYPED_TEST(ResolveBoundsTest, TheCapTracksTheLabelsFirstMetBeforeTheCut) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 202);
    LabelOracle oracle(*r.anno, r.cth.get());
    size_t exact = 0, inexact = 0;
    for (const Case &c : kCases) {
        // every label met, with its uncut profile and the k-mer it is first met at
        const SupportProfile all = resolve_support(oracle, r.q, discovery(c, kAll));
        ASSERT_FALSE(all.labels_truncated) << c.name();
        std::map<std::string, const LabelProfile *> uncut;
        std::map<uint64_t, size_t> first_met;      // k-mer -> labels first met there
        for (const LabelProfile &lp : all.labels) {
            uncut[lp.label.name] = &lp;
            ASSERT_FALSE(lp.runs.empty()) << c.name() << " " << lp.label.name;
            first_met[lp.runs.front().begin]++;
        }
        const size_t met = all.labels.size();
        for (size_t max_labels : { kAll, size_t(1), size_t(7) }) {
            const SupportProfile top = resolve_support(oracle, r.q, discovery(c, max_labels));
            for (size_t cap = 1; cap <= met; ++cap) {
                const std::string what = c.name() + " max_labels " + std::to_string(max_labels)
                        + " cap " + std::to_string(cap);
                ResolveOptions opts = discovery(c, max_labels);
                opts.discover_max_tracked = cap;
                const SupportProfile got = resolve_support(oracle, r.q, opts);
                // the cut k-mer: the first whose new labels would pass the cap
                size_t before = 0;
                uint64_t cut = 0;
                size_t fresh = 0;
                for (const auto &[kmer, n] : first_met) {
                    if (before + n > cap) {
                        cut = kmer;
                        fresh = n;
                        break;
                    }
                    before += n;
                }
                if (before == met) {
                    EXPECT_FALSE(got.labels_capped) << what;
                    expect_same_profile(top, got, what);
                    continue;
                }
                ASSERT_TRUE(got.labels_capped) << what;
                expect_ranking_exact_checkable(got, max_labels, what);
                const LabelCap &capped = *got.labels_capped;
                EXPECT_EQ(cap, capped.max_tracked) << what;
                EXPECT_EQ(before, capped.tracked) << what;
                EXPECT_EQ(cut, capped.cut_kmer) << what;
                EXPECT_EQ(fresh, capped.new_at_cut) << what;
                EXPECT_EQ(before + fresh, capped.met_at_least()) << what;
                const uint64_t bound = in_graph_from(all, cut);
                EXPECT_EQ(bound, capped.untracked_max_kmers) << what;
                // the whole query, every tracked label kept or counted in the truncation
                EXPECT_EQ(all.num_kmers, got.num_kmers) << what;
                EXPECT_FALSE(got.stop) << what;
                const size_t kept = std::min(before, max_labels);
                ASSERT_EQ(kept, got.labels.size()) << what;
                EXPECT_EQ(before > max_labels, got.labels_truncated.has_value()) << what;
                if (got.labels_truncated) {
                    EXPECT_EQ(before, got.labels_truncated->total) << what;
                }
                for (const LabelProfile &lp : got.labels) {
                    const LabelProfile &u = *uncut.at(lp.label.name);
                    EXPECT_LT(u.runs.front().begin, cut) << what << " " << lp.label.name;
                    EXPECT_EQ(u.kmers_supported, lp.kmers_supported) << what;
                    EXPECT_EQ(u.runs, lp.runs) << what << " " << lp.label.name;
                    EXPECT_EQ(u.trace_breaks, lp.trace_breaks) << what << " " << lp.label.name;
                }
                // every untracked label within the bound
                for (const LabelProfile &lp : all.labels) {
                    if (lp.runs.front().begin >= cut) {
                        EXPECT_LE(lp.kmers_supported, bound) << what << " " << lp.label.name;
                    }
                }
                // above the bound the kept labels are the uncut discovery's, in its order
                std::vector<std::string> above, above_uncut;
                for (const LabelProfile &lp : got.labels) {
                    if (lp.kmers_supported > bound)
                        above.push_back(lp.label.name);
                }
                for (const LabelProfile &lp : top.labels) {
                    if (lp.kmers_supported > bound)
                        above_uncut.push_back(lp.label.name);
                }
                EXPECT_EQ(above_uncut, above) << what;
                const bool exact_here = kept == max_labels && kept
                        && got.labels.back().kmers_supported > bound;
                EXPECT_EQ(exact_here, capped.ranking_exact) << what;
                if (capped.ranking_exact) {
                    EXPECT_EQ(names_of(top), names_of(got)) << what;
                    ++exact;
                } else {
                    ++inexact;
                }
            }
        }
    }
    // both statements occur
    EXPECT_GT(exact, 0u);
    EXPECT_GT(inexact, 0u);
}

// The memory bound stops the pass at a k-mer: the answer is then exactly the resolve of the
// query's prefix before it (its labels, their profiles, the truncation, the candidates, and
// under the cap its cut), the stop in phase rows by reason memory; the k-mer never moves back as
// the bound grows, and without a stop the profile is the unbounded one
TYPED_TEST(ResolveBoundsTest, AMemoryStopIsTheResolveOfAPrefix) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 303);
    LabelOracle oracle(*r.anno, r.cth.get());
    size_t stops = 0;
    for (const Case &c : kCases) {
        for (size_t cap : { size_t(0), size_t(30) }) {
            ResolveOptions opts = discovery(c, 10);
            opts.discover_max_tracked = cap;
            opts.batch_rows = 7;
            const SupportProfile plain = resolve_support(oracle, r.q, opts);
            uint64_t last = 0;
            bool ended = false;
            for (uint64_t bound = 1024; !ended && bound < (uint64_t(64) << 20);
                    bound += bound < 65536 ? 512 : bound / 8) {
                const std::string what = c.name() + " cap " + std::to_string(cap) + " bound "
                        + std::to_string(bound);
                ResolveOptions bounded = opts;
                bounded.discover_max_memory = bound;
                const SupportProfile got = resolve_support(oracle, r.q, bounded);
                expect_ranking_exact_checkable(got, 10, what);
                if (!got.stop) {
                    expect_same_profile(plain, got, what);
                    expect_same_cap(plain.labels_capped, got.labels_capped, what);
                    ended = true;
                    break;
                }
                ++stops;
                EXPECT_EQ(ResolveStop::MEMORY, got.stop->reason) << what;
                EXPECT_EQ(ResolveStop::ROWS, got.stop->phase) << what;
                EXPECT_EQ(bound, got.stop->max_memory) << what;
                const uint64_t x = got.stop->resolved_kmers;
                EXPECT_EQ(x, got.num_kmers) << what;
                EXPECT_EQ(plain.num_kmers, got.stop->query_kmers) << what;
                EXPECT_GE(x, last) << what;
                last = x;
                EXPECT_TRUE(std::any_of(plain.graph_runs.begin(), plain.graph_runs.end(),
                        [x](const KmerInterval &run) { return run.begin <= x && x < run.end; }))
                    << what << ": the stop is at k-mer " << x << ", not one in the graph";
                if (x == 0) {
                    EXPECT_TRUE(got.labels.empty()) << what;
                    continue;
                }
                const SupportProfile prefix
                        = resolve_support(oracle, r.q.substr(0, x + kK - 1), opts);
                expect_same_profile(prefix, got, what);
                expect_same_cap(prefix.labels_capped, got.labels_capped, what);
            }
            EXPECT_TRUE(ended) << c.name();
        }
    }
    EXPECT_GT(stops, 100u);
}

// Where the memory bound stops the pass depends on the k-mers before it and its own row only:
// the same k-mer whatever the batches and the rows kept for repeats
TYPED_TEST(ResolveBoundsTest, AMemoryStopDoesNotDependOnTheBatches) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 404);
    LabelOracle oracle(*r.anno, r.cth.get());
    size_t compared = 0;
    for (const Case &c : kCases) {
        for (uint64_t bound = 4096; bound < (uint64_t(1) << 20); bound *= 2) {
            std::optional<uint64_t> x;
            for (size_t rows : { size_t(1), size_t(3), size_t(64), size_t(4096) }) {
                for (uint64_t kept : { uint64_t(0), kResolveKeptBytes }) {
                    ResolveOptions opts = discovery(c, 10);
                    opts.discover_max_memory = bound;
                    opts.batch_rows = rows;
                    opts.kept_bytes = kept;
                    const SupportProfile got = resolve_support(oracle, r.q, opts);
                    const uint64_t at = got.stop ? got.stop->resolved_kmers : got.num_kmers;
                    const std::string what = c.name() + " bound " + std::to_string(bound)
                            + " rows " + std::to_string(rows) + " kept " + std::to_string(kept);
                    if (x) {
                        EXPECT_EQ(*x, at) << what;
                        ++compared;
                    }
                    x = at;
                }
            }
        }
    }
    EXPECT_GT(compared, 50u);
}

// The cut under a deadline: the stop's prefix holds the cut k-mer, and the statement is the
// prefix's (the untracked labels' bound counted to the k-mers resolved)
TYPED_TEST(ResolveBoundsTest, TheCutUnderADeadlineIsThePrefixs) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 505);
    LabelOracle oracle(*r.anno, r.cth.get());
    size_t both = 0;
    for (const Case &c : kCases) {
        ResolveOptions opts = discovery(c, 5);
        opts.discover_max_tracked = 12;
        opts.batch_rows = 3;
        for (size_t n = 1; n < 400; ++n) {
            size_t polls = 0;
            ResolveOptions timed = opts;
            timed.time_up = [&polls, n]() { return ++polls >= n; };
            const SupportProfile got = resolve_support(oracle, r.q, timed);
            expect_ranking_exact_checkable(got, 5, c.name());
            if (!got.stop)
                break;
            EXPECT_EQ(ResolveStop::TIME, got.stop->reason);
            const uint64_t x = got.stop->resolved_kmers;
            if (!x)
                continue;
            const SupportProfile prefix = resolve_support(oracle, r.q.substr(0, x + kK - 1), opts);
            const std::string what = c.name() + " deadline at read " + std::to_string(n);
            expect_same_profile(prefix, got, what);
            expect_same_cap(prefix.labels_capped, got.labels_capped, what);
            both += got.labels_capped.has_value();
        }
    }
    EXPECT_GT(both, 10u);
}

// A cut at the first k-mer in the graph tracks nothing: no label, no truncation, the statement
// says so (the query is too unspecific for the cap)
TEST(ResolveBounds, ACutAtTheFirstKmerTracksNothing) {
    // A, B and D carry the whole query, C and E its parts from k-mers 20 and 50
    const std::string q = random_seq(100, 606);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { q, q, q.substr(20), q, q.substr(50) }, { "A", "B", "C", "D", "E" });
    LabelOracle oracle(*anno);
    for (size_t cap : { size_t(1), size_t(2) }) {
        const std::string what = "cap " + std::to_string(cap);
        ResolveOptions opts = discovery({ LabelKind::COLUMN, Support::KMER }, 10);
        opts.discover_max_tracked = cap;
        const SupportProfile got = resolve_support(oracle, q, opts);
        ASSERT_TRUE(got.labels_capped) << what;
        expect_ranking_exact_checkable(got, 10, what);
        EXPECT_EQ(0u, got.labels_capped->tracked) << what;
        EXPECT_EQ(3u, got.labels_capped->new_at_cut) << what;
        EXPECT_EQ(0u, got.labels_capped->cut_kmer) << what;
        EXPECT_EQ(q.size() - kK + 1, got.labels_capped->untracked_max_kmers) << what;
        EXPECT_FALSE(got.labels_capped->ranking_exact) << what;
        EXPECT_TRUE(got.labels.empty()) << what;
        EXPECT_FALSE(got.labels_truncated) << what;
        EXPECT_TRUE(got.candidates.empty()) << what;
        EXPECT_EQ(iv_of(0, q.size() - kK + 1), got.graph_runs) << what;
    }
    // with room for the three: tracked, and C (first met at k-mer 20) is the cut
    ResolveOptions opts = discovery({ LabelKind::COLUMN, Support::KMER }, 10);
    opts.discover_max_tracked = 3;
    const SupportProfile three = resolve_support(oracle, q, opts);
    ASSERT_TRUE(three.labels_capped);
    EXPECT_EQ(3u, three.labels_capped->tracked);
    EXPECT_EQ(20u, three.labels_capped->cut_kmer);
    EXPECT_EQ(1u, three.labels_capped->new_at_cut);
    EXPECT_EQ((std::vector<std::string>{ "A", "B", "D" }), names_of(three));
    expect_ranking_exact_checkable(three, 10, "cap 3");
    opts.discover_max_tracked = 4;
    expect_ranking_exact_checkable(resolve_support(oracle, q, opts), 10, "cap 4");
    opts.discover_max_tracked = 5;
    EXPECT_FALSE(resolve_support(oracle, q, opts).labels_capped);
}

// A matrix without the budget-aware read (a row-flat annotation; a coordinate column annotation
// for tuple rows) is read by the default read and its rows counted after it: the same profile,
// and a stop that is the resolve of a prefix
TEST(ResolveBounds, ReadsThatAreNotChargedAreCountedAfterTheRead) {
    const std::string q = random_seq(400, 707);
    std::mt19937 gen(708);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 150; ++i) {
        const size_t a = gen() % 350;
        seqs.push_back(q.substr(a, 20 + gen() % std::min<size_t>(200, q.size() - a - 20)));
        labels.push_back("L" + std::to_string(i));
    }
    auto anno = build_anno_graph<DBGSuccinct, annot::RowFlatAnnotator>(kK, seqs, labels);
    LabelOracle oracle(*anno);
    ASSERT_FALSE(oracle.decode_charged());
    ResolveOptions opts = discovery({ LabelKind::COLUMN, Support::KMER }, 20);
    const SupportProfile plain = resolve_support(oracle, q, opts);
    size_t stops = 0;
    for (uint64_t bound = 1024; bound < (uint64_t(16) << 20); bound += bound / 4) {
        ResolveOptions bounded = opts;
        bounded.discover_max_memory = bound;
        const SupportProfile got = resolve_support(oracle, q, bounded);
        const std::string what = "bound " + std::to_string(bound);
        if (!got.stop) {
            expect_same_profile(plain, got, what);
            break;
        }
        ++stops;
        EXPECT_EQ(ResolveStop::MEMORY, got.stop->reason) << what;
        const uint64_t x = got.stop->resolved_kmers;
        if (x)
            expect_same_profile(resolve_support(oracle, q.substr(0, x + kK - 1), opts), got, what);
    }
    EXPECT_GT(stops, 3u);
}


// ---- the route (process_resolve_request): the members and the byte identity

// An index of |r|'s records with their mapping, as a server loads it (the annotation reloaded
// as a coordinate column annotation, the mapping attached)
std::unique_ptr<AnnotatedDBG> with_mapping(const Records &r) {
    std::filesystem::path dir = utils::create_temp_dir("", "test_resolve_bounds");
    const std::string prefix = dir / "coords";
    r.anno->get_annotator().serialize(prefix);
    auto annotator = std::make_unique<annot::ColumnCoordAnnotator>();
    EXPECT_TRUE(annotator->load(prefix));
    utils::remove_temp_dir(dir);
    auto mapping = std::make_shared<annot::CoordToHeader>(*r.cth);
    auto graph = std::dynamic_pointer_cast<DeBruijnGraph>(
            std::const_pointer_cast<SequenceGraph>(r.anno->get_graph_ptr()));
    return std::make_unique<AnnotatedDBG>(graph, std::move(annotator), false, mapping);
}

std::string compact(const Json::Value &v) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    return Json::writeString(builder, v);
}

Json::Value untimed(Json::Value v) {
    v.removeMember("timing");
    return v;
}

// The server's default bounds bind on none of these requests: every answer byte for byte the
// one without bounds (the CLI's), the discovery's work.memory_bytes included
TEST(ResolveBoundsRoute, AnswersThatBindNeitherAreUnchanged) {
    const Records r = make_records<DBGSuccinct, annot::ColumnCompressed<>>(40, 5, 808);
    const auto index = with_mapping(r);
    const cli::ResolveTimeLimits none;
    const cli::ResolveTimeLimits defaults { 0, cli::kResolveFinalizeMs, 1'000'000, 1024 };
    size_t n = 0;
    for (const char *kind : { "column", "header" }) {
        for (const char *support : { "kmer", "trace" }) {
            for (int max_labels : { 1000, 3 }) {
                Json::Value request;
                request["sequence"] = r.q;
                request["discover"]["kind"] = kind;
                request["discover"]["max_labels"] = max_labels;
                request["support"] = support;
                request["select"]["policy"] = "max_support";
                const std::string what = compact(request).substr(r.q.size());
                const Json::Value a = cli::process_resolve_request(request, *index, "r1", 0,
                                                                    nullptr, nullptr, none);
                const Json::Value b = cli::process_resolve_request(request, *index, "r1", 0,
                                                                    nullptr, nullptr, defaults);
                EXPECT_EQ(compact(untimed(a)), compact(untimed(b))) << what;
                EXPECT_FALSE(b.isMember("labels_capped")) << what;
                EXPECT_FALSE(b.isMember("stop")) << what;
                ++n;
            }
        }
    }
    EXPECT_EQ(8u, n);
}

// The cap states its cut (labels_capped, only where it binds), an explicit selection naming a
// label it did not keep is not made (selection null, the reason in the cut's message), and the
// memory bound's stop is stated without a deadline (reason memory, max_memory_mb)
TEST(ResolveBoundsRoute, TheCutAndTheMemoryStopAreStated) {
    const Records r = make_records<DBGSuccinct, annot::ColumnCompressed<>>(40, 5, 909);
    const auto index = with_mapping(r);
    Json::Value request;
    request["sequence"] = r.q;
    request["discover"]["kind"] = "header";
    request["discover"]["max_labels"] = 4;
    const cli::ResolveTimeLimits capped { 0, cli::kResolveFinalizeMs, 20, 0 };
    const Json::Value out = cli::process_resolve_request(request, *index, "r1", 0, nullptr,
                                                         nullptr, capped);
    ASSERT_TRUE(out.isMember("labels_capped")) << compact(out);
    const Json::Value &c = out["labels_capped"];
    EXPECT_EQ((std::vector<std::string>{ "cut_kmer", "max_tracked", "message", "met_at_least",
                                         "new_at_cut", "ranking_exact", "remainder_from_bp",
                                         "tracked", "untracked_max_kmers" }),
              c.getMemberNames());
    EXPECT_EQ(20u, c["max_tracked"].asUInt64());
    EXPECT_LE(c["tracked"].asUInt64(), 20u);
    EXPECT_EQ(c["cut_kmer"].asUInt64(), c["remainder_from_bp"].asUInt64());
    EXPECT_EQ(c["tracked"].asUInt64() + c["new_at_cut"].asUInt64(),
              c["met_at_least"].asUInt64());
    EXPECT_GT(c["tracked"].asUInt64() + c["new_at_cut"].asUInt64(), 20u);
    EXPECT_TRUE(c["ranking_exact"].isBool());
    // checkable from the answer alone: the kept list has max_labels entries and its last
    // supports more than untracked_max_kmers
    const Json::Value &kept = out["labels"];
    EXPECT_EQ(kept.size() == 4u
                  && kept[kept.size() - 1]["kmers_supported"].asUInt64()
                         > c["untracked_max_kmers"].asUInt64(),
              c["ranking_exact"].asBool());
    const std::string message = c["message"].asString();
    EXPECT_EQ(0u, message.find("--resolve-max-discover-labels 20: the "));
    EXPECT_NE(std::string::npos, message.find("Resolve the sequence from base "
                                              + std::to_string(c["cut_kmer"].asUInt64())));
    EXPECT_FALSE(out.isMember("stop"));

    // an explicit seed naming a label the capped discovery did not keep: not made
    Json::Value select = request;
    const uint64_t first = out["graph_runs"][0][0].asUInt64();
    Json::Value seed;
    seed["kmer_interval"].append(first);
    seed["kmer_interval"].append(first + 1);
    seed["labels"].append("no_such_record");
    select["select"]["policy"] = "explicit";
    select["select"]["seeds"].append(seed);
    const Json::Value blocked = cli::process_resolve_request(select, *index, "r1", 0, nullptr,
                                                             nullptr, capped);
    EXPECT_TRUE(blocked["selection"].isNull());
    EXPECT_NE(std::string::npos, blocked["labels_capped"]["message"].asString().find(
            "; the explicit selection was not made (selection null): select.seeds[0] names a "
            "label the capped discovery did not keep"));
    // without the cap it is the 400 it always was
    EXPECT_THROW(cli::process_resolve_request(select, *index, "r1", 0, nullptr, nullptr, {}),
                 cli::InvalidRequest);

    // the memory bound: a 20 kbp query whose pass does not fit 1 MiB stops before its first
    // row, with the stop block and no deadline
    std::string long_query;
    while (long_query.size() < 20000) {
        long_query += r.q;
    }
    Json::Value big = request;
    big["sequence"] = long_query;
    const cli::ResolveTimeLimits small { 0, cli::kResolveFinalizeMs, 0, 1 };
    const Json::Value stopped = cli::process_resolve_request(big, *index, "r1", 0, nullptr,
                                                             nullptr, small);
    ASSERT_TRUE(stopped.isMember("stop")) << compact(stopped).substr(0, 400);
    const Json::Value &s = stopped["stop"];
    EXPECT_EQ("memory", s["reason"].asString());
    EXPECT_EQ("rows", s["phase"].asString());
    EXPECT_EQ(1u, s["max_memory_mb"].asUInt64());
    EXPECT_EQ(s["resolved_kmers"].asUInt64(), stopped["num_kmers"].asUInt64());
    EXPECT_EQ(s["resolved_kmers"].asUInt64(), s["remainder_from_bp"].asUInt64());
    EXPECT_EQ(0u, s["message"].asString().find("--resolve-max-discover-memory-mb 1 ran out while "
                                               "the annotation rows of the query's k-mers were "
                                               "read: this answer is exactly the resolve of the "
                                               "sequence's first "));
    EXPECT_FALSE(stopped.isMember("limits"));
    // with a deadline the stop is stated as before, beside the limits
    big["bounds"]["time_budget_ms"] = 100000;
    const Json::Value timed = cli::process_resolve_request(big, *index, "r1", 0, nullptr,
                                                           nullptr, small);
    EXPECT_EQ("memory", timed["stop"]["reason"].asString());
    EXPECT_TRUE(timed.isMember("limits"));
    EXPECT_EQ(compact(s), compact(timed["stop"]));
}

// The capabilities' `discovery` member of the `resolve` block (§4.6, §10.3)
TEST(ResolveBoundsRoute, CapabilitiesBlock) {
    const Json::Value c = cli::resolve_capabilities_json({ 30000, cli::kResolveFinalizeMs,
                                                           1'000'000, 1024 });
    EXPECT_EQ((std::vector<std::string>{ "discovery", "time_budget" }), c.getMemberNames());
    EXPECT_EQ("{\"max_memory_mb\":1024,\"max_tracked_labels\":1000000,\"rule\":"
              "\"SPEC-labeled-traversal-core.md section 4.6\",\"stop_reasons\":[\"time\","
              "\"memory\"]}", compact(c["discovery"]));
    EXPECT_TRUE(c["discovery"]["max_memory_mb"].isUInt64());
    EXPECT_EQ("{\"max_memory_mb\":0,\"max_tracked_labels\":0,\"rule\":"
              "\"SPEC-labeled-traversal-core.md section 4.6\",\"stop_reasons\":[\"time\","
              "\"memory\"]}", compact(cli::resolve_capabilities_json({})["discovery"]));
}

// The reference samples of the bounds (SPEC-labeled-traversal-core.md §4.6), real runs: 30 a cut by the label cap
// (labels_capped with the members the route writes, ranking_exact checkable from the entry alone), 31 a memory
// stop without a deadline (the members of a time stop's block, reason memory and max_memory_mb, the entry the
// resolve of its first resolved_kmers k-mers)
TEST(ResolveBoundsRoute, Samples30And31) {
    auto sample = [](const std::string &name) {
        const std::string path = source_path("tests/data/multigraph_contract/" + name);
        std::ifstream in(path, std::ios::binary);
        EXPECT_TRUE(in.good()) << path << " not found";
        Json::Value doc;
        std::string errors;
        Json::CharReaderBuilder builder;
        EXPECT_TRUE(Json::parseFromStream(builder, in, &doc, &errors)) << errors;
        return doc;
    };
    const Json::Value s30 = sample("30_resolve_labels_capped.json");
    EXPECT_EQ(200, s30["status"].asInt());
    const uint64_t max_labels = s30["request"]["discover"]["max_labels"].asUInt64();
    ASSERT_EQ(2u, s30["body"]["answers"].size());
    for (const Json::Value &entry : s30["body"]["answers"]) {
        const Json::Value &c = entry["labels_capped"];
        EXPECT_EQ((std::vector<std::string>{ "cut_kmer", "max_tracked", "message", "met_at_least",
                                             "new_at_cut", "ranking_exact", "remainder_from_bp",
                                             "tracked", "untracked_max_kmers" }),
                  c.getMemberNames());
        const Json::Value &kept = entry["labels"];
        EXPECT_EQ(kept.size() == max_labels
                      && kept[kept.size() - 1]["kmers_supported"].asUInt64()
                             > c["untracked_max_kmers"].asUInt64(),
                  c["ranking_exact"].asBool());
        EXPECT_TRUE(c["ranking_exact"].asBool());
        EXPECT_EQ(c["cut_kmer"].asUInt64(), c["remainder_from_bp"].asUInt64());
        EXPECT_EQ(c["tracked"].asUInt64() + c["new_at_cut"].asUInt64(),
                  c["met_at_least"].asUInt64());
        EXPECT_FALSE(entry.isMember("stop"));
    }
    const Json::Value s31 = sample("31_resolve_memory_stop.json");
    EXPECT_EQ(200, s31["status"].asInt());
    ASSERT_EQ(1u, s31["body"]["answers"].size());
    const Json::Value &entry = s31["body"]["answers"][0];
    const Json::Value &stop = entry["stop"];
    EXPECT_EQ((std::vector<std::string>{ "max_memory_mb", "message", "phase", "query_bp", "query_kmers",
                                         "reason", "remainder_from_bp", "resolved_bp",
                                         "resolved_kmers" }),
              stop.getMemberNames());
    EXPECT_EQ("memory", stop["reason"].asString());
    EXPECT_EQ("rows", stop["phase"].asString());
    EXPECT_EQ(stop["resolved_kmers"].asUInt64(), entry["num_kmers"].asUInt64());
    EXPECT_EQ(stop["resolved_kmers"].asUInt64(), stop["remainder_from_bp"].asUInt64());
    EXPECT_FALSE(entry.isMember("limits"));
}

} // namespace
