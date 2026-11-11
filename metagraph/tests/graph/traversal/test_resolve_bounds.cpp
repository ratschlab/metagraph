#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
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
#include "graph/traversal/charged_table.hpp"
#include "graph/traversal/kept_rows.hpp"
#include "graph/traversal/seeded_hash.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/traverse.hpp"
#include "common/utils/file_utils.hpp"

#if defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(__GLIBC__)
#include <malloc.h>
#endif

#if USE_JEMALLOC
extern "C" int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen);
#endif


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
                // the labels present at the cut k-mer: above the cap, a discovery of the sequence
                // from it cuts at its first k-mer and tracks nothing (the message says so)
                if (max_labels == kAll) {
                    const SupportProfile rest = resolve_support(oracle, r.q.substr(cut), opts);
                    ASSERT_TRUE(capped.labels_at_cut <= cap || rest.labels_capped) << what;
                    if (capped.labels_at_cut > cap) {
                        EXPECT_EQ(0u, rest.labels_capped->tracked) << what;
                        EXPECT_EQ(rest.graph_runs.front().begin, rest.labels_capped->cut_kmer)
                            << what;
                    } else {
                        EXPECT_TRUE(!rest.labels_capped || rest.labels_capped->tracked) << what;
                    }
                }
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

// A client following remainder_from_bp either progresses or is told it cannot: every memory stop
// past the query's first k-mer in the graph resolves at least that k-mer, and a stop at it states
// cannot_progress -- and then a resolve of the sequence from there stops there again
TYPED_TEST(ResolveBoundsTest, TheRemainderProgressesOrSaysItCannot) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 313);
    LabelOracle oracle(*r.anno, r.cth.get());
    size_t walks = 0, stuck = 0, steps = 0;
    for (const Case &c : kCases) {
        for (uint64_t bound = 1024; bound < (uint64_t(1) << 20);
                bound += bound < 65536 ? 1024 : bound / 4) {
            ResolveOptions opts = discovery(c, 10);
            opts.discover_max_memory = bound;
            std::string rest = r.q;
            ++walks;
            for (size_t step = 0; step < 1000 && rest.size() >= kK; ++step) {
                const std::string what = c.name() + " bound " + std::to_string(bound) + " step "
                        + std::to_string(step);
                const SupportProfile got = resolve_support(oracle, rest, opts);
                if (!got.stop)
                    break;
                ++steps;
                const uint64_t x = got.stop->resolved_kmers;
                const uint64_t first = got.stop->query_graph_runs.empty()
                        ? got.stop->query_kmers : got.stop->query_graph_runs.front().begin;
                EXPECT_EQ(x == first, got.stop->cannot_progress) << what;
                if (got.stop->cannot_progress) {
                    // the claim: the sequence from x stops at its first k-mer in the graph again
                    const SupportProfile again = resolve_support(oracle, rest.substr(x), opts);
                    ASSERT_TRUE(again.stop) << what;
                    EXPECT_TRUE(again.stop->cannot_progress) << what;
                    EXPECT_EQ(again.stop->query_graph_runs.front().begin,
                              again.stop->resolved_kmers) << what;
                    ++stuck;
                    break;
                }
                EXPECT_GT(x, first) << what;
                rest = rest.substr(x);
            }
        }
    }
    EXPECT_GT(stuck, 0u);
    EXPECT_GT(steps, stuck);
    EXPECT_GT(walks, 10u);
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
    // the remainder's discovery from the cut: resolve it, or (the cut k-mer alone holding more
    // labels than the cap) it cannot track them
    const std::string from = "from base " + std::to_string(c["cut_kmer"].asUInt64());
    EXPECT_NE(message.find("Resolve the sequence " + from) == std::string::npos,
              message.find("A discovery of the sequence " + from + " cannot track them under "
                           "this cap: its first k-mer alone holds ") == std::string::npos)
        << message;
    EXPECT_FALSE(out.isMember("stop"));

    // an explicit seed naming a label the capped discovery did not keep: not made where the
    // ranking is not exact (here fewer labels tracked than max_labels 1000); where it is
    // exact, the uncut answer's 400 (ExplicitSeedsUnderAnExactCut)
    Json::Value select = request;
    select["discover"]["max_labels"] = 1000;
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
    // at the query's first k-mer in the graph: the sequence from there cannot progress either
    // (the graph runs cut at x: none before the query's first k-mer in the graph)
    EXPECT_EQ(0u, stopped["graph_runs"].size());
    EXPECT_NE(std::string::npos, s["message"].asString().find(
            "; it cannot be resolved past base " + std::to_string(s["resolved_kmers"].asUInt64())
            + " under this bound: its first k-mer in the graph does not fit beside the query's "
              "own setup, and a resolve of the sequence from base "
            + std::to_string(s["resolved_kmers"].asUInt64()) + " stops there again"))
        << s["message"].asString();
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
                                                           1'000'000, 1024, 10'000, 10'000,
                                                           1024 });
    EXPECT_EQ((std::vector<std::string>{ "discovery", "max_explicit_labels",
                                         "max_explicit_seeds", "seed_order", "time_budget" }),
              c.getMemberNames());
    EXPECT_TRUE(c["max_explicit_labels"].isUInt64());
    EXPECT_EQ(10000u, c["max_explicit_labels"].asUInt64());
    EXPECT_EQ(0u, cli::resolve_capabilities_json({})["max_explicit_labels"].asUInt64());
    // the cap on a request's explicit seeds (GPT review 7, finding 1; 0: none)
    EXPECT_TRUE(c["max_explicit_seeds"].isUInt64());
    EXPECT_EQ(1024u, c["max_explicit_seeds"].asUInt64());
    EXPECT_EQ(0u, cli::resolve_capabilities_json({})["max_explicit_seeds"].asUInt64());
    // the order of a selection's seeds (§4.3), whatever the bounds
    EXPECT_EQ("rank", c["seed_order"].asString());
    EXPECT_EQ("rank", cli::resolve_capabilities_json({})["seed_order"].asString());
    EXPECT_EQ("{\"max_answer_labels\":10000,\"max_memory_mb\":1024,\"max_tracked_labels\":1000000,"
              "\"rule\":\"SPEC-labeled-traversal-core.md section 4.6\",\"stop_reasons\":[\"time\","
              "\"memory\"]}", compact(c["discovery"]));
    EXPECT_TRUE(c["discovery"]["max_memory_mb"].isUInt64());
    EXPECT_EQ("{\"max_answer_labels\":0,\"max_memory_mb\":0,\"max_tracked_labels\":0,\"rule\":"
              "\"SPEC-labeled-traversal-core.md section 4.6\",\"stop_reasons\":[\"time\","
              "\"memory\"]}", compact(cli::resolve_capabilities_json({})["discovery"]));
}

// A discover.max_labels above the server's maximum (--resolve-max-discover-answer-labels) is
// lowered to it and stated in limits.clamped, with or without a deadline; the answer is that of
// the maximum; at the maximum, or with none, nothing is stated
TEST(ResolveBoundsRoute, TheAnswersLabelsAreClamped) {
    const Records r = make_records<DBGSuccinct, annot::ColumnCompressed<>>(40, 5, 707);
    const auto index = with_mapping(r);
    Json::Value request;
    request["sequence"] = r.q;
    request["discover"]["kind"] = "header";
    request["discover"]["max_labels"] = 50;
    cli::ResolveTimeLimits limits;
    limits.max_discover_answer_labels = 7;
    const Json::Value clamped = cli::process_resolve_request(request, *index, "r1", 0, nullptr,
                                                             nullptr, limits);
    EXPECT_EQ("{\"clamped\":[{\"effective\":7,\"field\":\"discover.max_labels\",\"requested\":50}]}",
              compact(clamped["limits"]));
    EXPECT_EQ(7u, clamped["labels"].size());
    EXPECT_FALSE(clamped.isMember("stop"));
    // the answer of max_labels 7
    Json::Value seven = request;
    seven["discover"]["max_labels"] = 7;
    Json::Value plain = cli::process_resolve_request(seven, *index, "r1", 0, nullptr, nullptr,
                                                     limits);
    EXPECT_FALSE(plain.isMember("limits"));
    Json::Value got = clamped;
    got.removeMember("limits");
    got.removeMember("timing");
    plain.removeMember("timing");
    EXPECT_EQ(compact(plain), compact(got));
    // with a deadline: the clamp first, beside the budget
    Json::Value timed = request;
    timed["bounds"]["time_budget_ms"] = 100000;
    limits.max_time_ms = 60000;
    const Json::Value both = cli::process_resolve_request(timed, *index, "r1", 0, nullptr,
                                                          nullptr, limits);
    EXPECT_EQ("[{\"effective\":7,\"field\":\"discover.max_labels\",\"requested\":50},"
              "{\"effective\":60000,\"field\":\"bounds.time_budget_ms\",\"requested\":100000}]",
              compact(both["limits"]["clamped"]));
    // none without the maximum or under it
    EXPECT_FALSE(cli::process_resolve_request(request, *index, "r1", 0, nullptr, nullptr, {})
                         .isMember("limits"));
    limits = {};
    limits.max_discover_answer_labels = 50;
    EXPECT_FALSE(cli::process_resolve_request(request, *index, "r1", 0, nullptr, nullptr, limits)
                         .isMember("limits"));
}

// N3 of the review: a cut whose ranking is exact kept the uncut discovery's labels, so an
// explicit seed naming a label it did not keep is refused as the uncut answer refuses it (A, B
// and D carry the whole query and are kept; C, first met at k-mer 20, is not kept by either)
TEST(ResolveBoundsRoute, ExplicitSeedsUnderAnExactCut) {
    const std::string q = random_seq(100, 606);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { q, q, q.substr(20), q, q.substr(50) }, { "A", "B", "C", "D", "E" });
    Json::Value request;
    request["sequence"] = q;
    request["discover"]["max_labels"] = 3;
    Json::Value seed;
    seed["kmer_interval"].append(30);
    seed["kmer_interval"].append(40);
    seed["labels"].append("A");
    seed["labels"].append("C");
    request["select"]["policy"] = "explicit";
    request["select"]["seeds"].append(seed);
    cli::ResolveTimeLimits limits;
    limits.max_discover_labels = 3;
    Json::Value discover_only = request;
    discover_only.removeMember("select");
    const Json::Value cut = cli::process_resolve_request(discover_only, *anno, "r1", 0, nullptr,
                                                         nullptr, limits);
    ASSERT_TRUE(cut["labels_capped"]["ranking_exact"].asBool()) << compact(cut);
    std::string capped_error, uncut_error;
    try {
        cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, limits);
    } catch (const cli::InvalidRequest &e) {
        capped_error = e.what();
    }
    try {
        cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, {});
    } catch (const cli::InvalidRequest &e) {
        uncut_error = e.what();
    }
    EXPECT_EQ("Explicit seed label was not profiled: 'C'", uncut_error);
    EXPECT_EQ(uncut_error, capped_error);
    // a seed on the kept labels is selected as without the cap
    request["select"]["seeds"][0]["labels"] = Json::Value(Json::arrayValue);
    request["select"]["seeds"][0]["labels"].append("A");
    request["select"]["seeds"][0]["labels"].append("D");
    Json::Value a = cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, limits);
    Json::Value b = cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, {});
    EXPECT_EQ(compact(b["selection"]), compact(a["selection"]));
}

// NF1 of the verification: an explicit seed naming a label the answer-labels clamp may have left
// out is refused naming the clamp (A, B and D carry the whole query; C, from k-mer 20, ranks
// fourth: kept with max_labels 5, not with the clamp at 3); where the clamp dropped nothing, or
// did not lower max_labels, the 400 is the plain one
TEST(ResolveBoundsRoute, ExplicitSeedsUnderTheAnswerClamp) {
    const std::string q = random_seq(100, 606);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { q, q, q.substr(20), q, q.substr(50) }, { "A", "B", "C", "D", "E" });
    Json::Value request;
    request["sequence"] = q;
    request["discover"]["max_labels"] = 5;
    Json::Value seed;
    seed["kmer_interval"].append(30);
    seed["kmer_interval"].append(40);
    seed["labels"].append("C");
    request["select"]["policy"] = "explicit";
    request["select"]["seeds"].append(seed);
    auto refusal = [&](const Json::Value &body, const cli::ResolveTimeLimits &limits) {
        try {
            cli::process_resolve_request(body, *anno, "r1", 0, nullptr, nullptr, limits);
        } catch (const cli::InvalidRequest &e) {
            return std::string(e.what());
        }
        return std::string();
    };
    // without the clamp C is kept and the seed made
    EXPECT_EQ("", refusal(request, {}));
    cli::ResolveTimeLimits limits;
    limits.max_discover_answer_labels = 3;
    EXPECT_EQ("Explicit seed label was not profiled: 'C'; discover.max_labels was lowered from 5 "
              "to 3 (--resolve-max-discover-answer-labels), and the discovery kept 3 of the 5 "
              "labels it ranked: the label may be one ranked past 3. Ask for at most 3 labels "
              "and choose the seeds among them, or resolve with the seed's label named in labels",
              refusal(request, limits));
    // a label the discovery never met, the clamp having dropped labels: the same words (the
    // answer cannot tell the two apart)
    Json::Value unknown = request;
    unknown["select"]["seeds"][0]["labels"][0] = "Z";
    EXPECT_EQ(0u, refusal(unknown, limits).find(
            "Explicit seed label was not profiled: 'Z'; discover.max_labels was lowered from 5 "));
    // the clamp lowering max_labels but dropping nothing: the plain 400
    limits.max_discover_answer_labels = 6;
    unknown["discover"]["max_labels"] = 7;
    EXPECT_EQ("Explicit seed label was not profiled: 'Z'", refusal(unknown, limits));
    // max_labels at the clamp (not lowered): the plain 400, C ranked past 3 by the request
    limits.max_discover_answer_labels = 3;
    request["discover"]["max_labels"] = 3;
    EXPECT_EQ("Explicit seed label was not profiled: 'C'", refusal(request, limits));
}

// ---- the rows a pass keeps (KeptRows, SPEC §4.6; GPT review 6, finding 8)

// The table of kept rows never shrinks while it lives (erase and clear keep its buckets), so
// what the budget holds for it covers its buckets until the table is freed: a last use releases
// the row's own charge and keeps its entry's (a high-water count, charged again only past the
// most entries held at once), and drop() releases everything and frees the buckets. Before,
// a last use and a drop released the entries' charge while the buckets stayed: 4.5 M entries
// left 448 MiB of buckets uncharged
TEST(ResolveBounds, KeptRowsChargeTheirTableUntilItIsFreed) {
    using annot::matrix::DecodeBudget;
    using Row = annot::matrix::BinaryMatrix::SetBitPositions;
    using Kept = KeptRows<Row>;
    // the table's own bytes: its buckets (a neighbourhood bitmap and the value each, the
    // overflow buckets of the last neighbourhood beside them)
    constexpr uint64_t kBucket = sizeof(typename Kept::Table::value_type) + sizeof(uint64_t);
    auto table_bytes = [&](const Kept &k) {
        return k.bucket_count() ? (k.bucket_count() + 61) * kBucket : 0;
    };
    DecodeBudget budget;
    Kept kept(budget);
    constexpr size_t n = 100'000;
    for (uint64_t key = 0; key < n; ++key) {
        ASSERT_TRUE(budget.charge(8));   // the row's own charge, which the table takes over
        ASSERT_TRUE(kept.keep(key, Row { 7 }, 8, 8));
    }
    EXPECT_EQ(n * (Kept::kEntryBytes + 8), budget.held());
    EXPECT_GE(budget.held(), table_bytes(kept));
    // every row used for the last time: the rows' charges released, the entries' kept while
    // the buckets are
    for (uint64_t key = 0; key < n; ++key) {
        kept.last_use(kept.find(key));
    }
    EXPECT_TRUE(kept.empty());
    EXPECT_TRUE(kept.holds());
    EXPECT_GT(kept.bucket_count(), n);
    EXPECT_EQ(n * Kept::kEntryBytes, budget.held());
    EXPECT_GE(budget.held(), table_bytes(kept));
    // kept again up to the high water: no entry charged twice; one past it: one more
    for (uint64_t key = 0; key <= n; ++key) {
        ASSERT_TRUE(budget.charge(8));
        ASSERT_TRUE(kept.keep(n + key, Row { 7 }, 8, 8));
    }
    EXPECT_EQ((n + 1) * (Kept::kEntryBytes + 8), budget.held());
    EXPECT_GE(budget.held(), table_bytes(kept));
    // dropped: everything released and the table freed
    kept.drop();
    EXPECT_FALSE(kept.holds());
    EXPECT_EQ(0u, kept.bucket_count());
    EXPECT_EQ(0u, budget.held());
    // an entry whose charge does not fit is not kept, and nothing is charged
    DecodeBudget small(Kept::kEntryBytes - 1);
    Kept none(small);
    EXPECT_FALSE(none.keep(1, Row { 1 }, 8, 0));
    EXPECT_FALSE(none.holds());
    EXPECT_EQ(0u, small.held());
}

// The pass's node-keyed tables are seeded (SeededHash; SPEC §4.6, T3 review R2, F1). A table
// takes a key's bucket from the hash's low bits and doubles, down to 10% load, when a key's
// neighbourhood is full: under the identity, a query's node indices chosen so that 63 of them
// share their low bits at each size grow a table far past its keys. Seeded, the same keys
// leave it as natural ones do, its buckets within the per-entry charge (hopscotch_entry_bytes:
// 3 buckets an entry at a growth at 80% load or above, so the table may hold up to 2.5 buckets
// a key). R2 measured 1 M such kept rows at 1,344 MiB against 183 MiB charged
TEST(ResolveBounds, ChosenNodeIndicesDoNotGrowTheTablesPastTheirCharge) {
    using annot::matrix::DecodeBudget;
    using Row = annot::matrix::BinaryMatrix::SetBitPositions;
    using Kept = KeptRows<Row>;
    namespace hop = tsl::detail_hopscotch_hash;
    using KeptBucket = hop::hopscotch_bucket<Kept::Table::value_type, 62, false>;
    using Positions = std::pair<uint64_t, size_t>;
    using Record = std::pair<std::pair<uint64_t, uint64_t>, uint64_t>;
    // each charge covers 3 of the pinned table's buckets at 80% load: a kept row of set bits
    // 180 B (a 48 B bucket), of tuples 210 B (56 B), a setup key 90 B, a record 120 B
    using KeptTuples = KeptRows<annot::matrix::MultiIntMatrix::RowTuples>;
    EXPECT_EQ(180u, Kept::kEntryBytes);
    EXPECT_EQ(210u, KeptTuples::kEntryBytes);
    EXPECT_GE(Kept::kEntryBytes * 4, 3 * sizeof(KeptBucket) * 5);
    EXPECT_GE(KeptTuples::kEntryBytes * 4,
              3 * sizeof(hop::hopscotch_bucket<KeptTuples::Table::value_type, 62, false>) * 5);
    EXPECT_EQ(90u, hopscotch_entry_bytes(sizeof(Positions)));
    EXPECT_GE(90u * 4, 3 * sizeof(hop::hopscotch_bucket<Positions, 62, false>) * 5);
    EXPECT_EQ(kDiscoverIndexEntryBytes, hopscotch_entry_bytes(sizeof(Record)));
    EXPECT_GE(kDiscoverIndexEntryBytes * 4,
              3 * sizeof(hop::hopscotch_bucket<Record, 62, false>) * 5);
    // one seed a process, and not the identity
    EXPECT_EQ(SeededHash()(12345), SeededHash()(12345));
    EXPECT_NE(12345u, SeededHash()(12345));

    // natural node indices, then 63 at each size sharing their low bits (R2's hop.cpp)
    std::mt19937_64 rng(1);
    constexpr size_t filler = 50'000;
    std::vector<uint64_t> keys;
    for (size_t i = 0; i < filler; ++i) {
        keys.push_back(rng() & ((1ull << 36) - 1));
    }
    size_t level = 1;
    while (0.9 * level <= filler) {
        level <<= 1;
    }
    for ( ; filler >= level / 10; level <<= 1) {
        const uint64_t c = rng() & (level - 1);
        for (uint64_t t = 0; t < 63; ++t) {
            keys.push_back(((2 * t + (t & 1)) * level + c) | (1ull << 40));
        }
    }
    const double most_buckets = 2.5 * keys.size();
    // identity: the table outgrows the charge (the train's 80 B a key and more)
    tsl::hopscotch_map<uint64_t, size_t> identity;
    for (uint64_t key : keys) {
        identity.emplace(key, 0);
    }
    EXPECT_GT(identity.bucket_count(), most_buckets);
    // seeded: within it, the setup's table as the kept rows'
    tsl::hopscotch_map<uint64_t, size_t, SeededHash> seeded;
    for (uint64_t key : keys) {
        seeded.emplace(key, 0);
    }
    EXPECT_LE(seeded.bucket_count(), most_buckets);
    DecodeBudget budget;
    Kept kept(budget);
    for (uint64_t key : keys) {
        ASSERT_TRUE(kept.keep(key, Row { 7 }, 8, 0));
    }
    EXPECT_LE(kept.bucket_count(), most_buckets);
    // what is held covers the table at its last growth: its buckets and half as many before
    EXPECT_GE(budget.held(), kept.bucket_count() * sizeof(KeptBucket) * 3 / 2);
    kept.drop();
}

// ---- the tables' account (TableAccount, charged_table.hpp; SPEC §4.6; GPT review 7, finding 7)

// A hash whose seed is known (0): keys can be chosen for it by inverting mix64 (V7-B's crafting)
struct KnownSeedHash {
    size_t operator()(uint64_t key) const { return static_cast<size_t>(mix64(key)); }
};

uint64_t inverse_odd(uint64_t a) {
    uint64_t x = a;
    for (int i = 0; i < 6; ++i) {
        x *= 2 - a * x;
    }
    return x;
}

// the key KnownSeedHash maps to |h|
uint64_t unmix64(uint64_t h) {
    static const uint64_t i1 = inverse_odd(0xFF51AFD7ED558CCDULL);
    static const uint64_t i2 = inverse_odd(0xC4CEB9FE1A85EC53ULL);
    h ^= h >> 33;
    h *= i2;
    h ^= h >> 33;
    h *= i1;
    h ^= h >> 33;
    return h;
}

// Keys chosen for a known seed (V7-B): random keys to just past 10% of the buckets, then groups
// of 63 whose hashes agree modulo the bucket count but differ in its next bit, so that tsl grows
// the table on a full neighbourhood at a tenth of its load. |insert(key)| false ends the stream
template <class Insert, class Buckets>
void chosen_keys(size_t max_buckets, uint64_t seed, const Buckets &buckets, const Insert &insert) {
    std::mt19937_64 rng(seed);
    uint64_t group = 0;
    while (buckets() < max_buckets) {
        const size_t n = std::max<size_t>(buckets(), 1);
        for (size_t i = 0; i < n / 10 + 2; ++i) {
            if (!insert(rng()))
                return;
        }
        const uint64_t low = rng() & (n - 1);
        const size_t before = buckets();
        for (uint64_t i = 0; i < 63 && buckets() == before; ++i) {
            if (!insert(unmix64(low | ((i & 1) * n) | ((++group) << 40))))
                return;
        }
    }
}

// The kept rows' table is charged as allocated beyond its model, whatever the keys: keys chosen
// for the hash's seed double it at a tenth of its load, its buckets passing the model's
// kEntryBytes an entry several times (V7-B: up to 8 times at a growth), and the budget holds
// them all the same — at every step at least the table's buckets (with tsl's tail of 61) and its
// overflow nodes, at every growth the old and the new buckets at once. Under a bound a growth
// that does not fit is refused before it is made: the row is not kept (it stays with the
// caller) and the budget never passes the bound. Before, the model alone was charged
TEST(ResolveBounds, ChosenKeysCannotGrowTheKeptRowsPastTheirAccount) {
    using annot::matrix::DecodeBudget;
    using annot::matrix::alloc_bytes;
    using Row = annot::matrix::BinaryMatrix::SetBitPositions;
    using Kept = KeptRows<Row, uint64_t, KnownSeedHash>;
    namespace hop = tsl::detail_hopscotch_hash;
    constexpr uint64_t kBucket = sizeof(hop::hopscotch_bucket<Kept::Table::value_type, 62, false>);
    auto table_bytes = [&](size_t buckets, size_t overflow) {
        return (buckets ? alloc_bytes((buckets + 61) * kBucket) : 0)
             + overflow * alloc_bytes(sizeof(Kept::Table::value_type) + 2 * sizeof(void*));
    };
    DecodeBudget budget;
    double most_over_model = 0;
    {
        Kept kept(budget);
        chosen_keys(1 << 16, 7, [&]() { return kept.bucket_count(); }, [&](uint64_t key) {
            const size_t before = kept.bucket_count();
            if (kept.contains(key))
                return true;
            Row row { 7 };
            EXPECT_TRUE(kept.keep(key, row, 8, 0));
            const uint64_t held = table_bytes(kept.bucket_count(), kept.overflow_size());
            EXPECT_GE(budget.held(), held) << kept.size();
            EXPECT_EQ(held, kept.account().held()) << kept.size();
            if (kept.bucket_count() != before && before) {
                // the growth: the old buckets and the new ones held at once
                EXPECT_GE(budget.peak(), table_bytes(before, 0)
                                         + table_bytes(kept.bucket_count(), 0)) << kept.size();
                most_over_model = std::max(most_over_model,
                        double(table_bytes(before, 0) + table_bytes(kept.bucket_count(), 0))
                            / (kept.size() * Kept::kEntryBytes));
            }
            return true;
        });
        EXPECT_GE(kept.bucket_count(), size_t(1) << 16);
        kept.drop();
        EXPECT_EQ(0u, budget.held());
    }
    // the model alone was passed several times over
    EXPECT_GT(most_over_model, 3.0);

    // under a bound: a growth that does not fit is refused before it is made
    for (uint64_t bound : { uint64_t(4) << 10, uint64_t(64) << 10, uint64_t(1) << 20 }) {
        DecodeBudget small(bound);
        Kept kept(small);
        size_t refused = 0;
        chosen_keys(1 << 20, 11, [&]() { return kept.bucket_count(); }, [&](uint64_t key) {
            if (kept.contains(key))
                return true;
            Row row { 7 };
            const size_t size = kept.size(), buckets = kept.bucket_count();
            if (!kept.keep(key, row, 8, 0)) {
                // nothing kept, the table as it was, the row still the caller's
                EXPECT_EQ(size, kept.size());
                EXPECT_EQ(buckets, kept.bucket_count());
                EXPECT_EQ(Row { 7 }, row);
                return ++refused < 3;
            }
            EXPECT_LE(small.held(), bound);
            EXPECT_GE(small.held(), table_bytes(kept.bucket_count(), kept.overflow_size()));
            return true;
        });
        EXPECT_GT(refused, 0u) << bound;
        EXPECT_LE(small.peak(), bound) << bound;
    }
    // the smallest bound an entry's model fits in: its table's first buckets (63) do not
    DecodeBudget one(Kept::kEntryBytes);
    Kept first(one);
    Row row { 1 };
    EXPECT_FALSE(first.keep(1, row, 8, 0));
    EXPECT_EQ(0u, first.size());
    EXPECT_EQ(Kept::kEntryBytes, one.held());
}

// A table whose growth cannot be refused where it happens (the record index, inside a k-mer's
// accumulation) is charged as allocated past the budget's maximum; the budget then refuses
// every charge, so the next check stops, and freeing the table releases it all
TEST(ResolveBounds, AnUnrefusedGrowthIsChargedPastTheMaximum) {
    using annot::matrix::DecodeBudget;
    using Table = tsl::hopscotch_map<uint64_t, uint64_t, KnownSeedHash, std::equal_to<uint64_t>,
                                     ChargedAllocator<std::pair<uint64_t, uint64_t>>>;
    DecodeBudget budget(64 << 10);
    {
        TableAccount account(budget, /* refuse */ false);
        Table table(0, KnownSeedHash(), std::equal_to<uint64_t>(),
                    ChargedAllocator<std::pair<uint64_t, uint64_t>>(&account));
        uint64_t key = 0;
        while (budget.held() <= budget.max_bytes()) {
            table.emplace(key++, 0);
        }
        EXPECT_GT(budget.held(), budget.max_bytes());
        EXPECT_EQ(account.held(), budget.held());
        EXPECT_FALSE(budget.charge(1));
        // the model covering the entries: what passed it is released, the rest held
        account.cover_beyond(table.size() * hopscotch_entry_bytes(sizeof(Table::value_type)));
        EXPECT_EQ(account.charged(), budget.held());
        EXPECT_EQ(std::max(account.held(), account.covered()), account.charged());
        Table(0, KnownSeedHash(), std::equal_to<uint64_t>(),
              ChargedAllocator<std::pair<uint64_t, uint64_t>>(&account)).swap(table);
        EXPECT_EQ(0u, account.held());
    }
    EXPECT_EQ(0u, budget.held());
    EXPECT_TRUE(budget.charge(1));
}

// ---- the cap on a request's explicit label lists (--resolve-max-explicit-labels, SPEC §4.1)

std::string refusal_of(const Json::Value &request, uint64_t cap) {
    try {
        cli::parse_resolve_request(request, cap);
    } catch (const cli::InvalidRequest &e) {
        return e.what();
    }
    return "";
}

// The request's `labels` and each explicit seed's, checked before anything else of the request
// (a list of N passes, N + 1 is refused, naming the list, the cap and the flag), `labels` first,
// then the seeds in order; a list of another type is left to the parse; 0 is no cap
TEST(ResolveBoundsRoute, TheExplicitLabelsCap) {
    auto names = [](size_t n, const std::string &prefix) {
        Json::Value list(Json::arrayValue);
        for (size_t i = 0; i < n; ++i) {
            list.append(prefix + std::to_string(i));
        }
        return list;
    };
    auto seed = [&](size_t n) {
        Json::Value s;
        s["kmer_interval"].append(0);
        s["kmer_interval"].append(1);
        s["labels"] = names(n, "s");
        return s;
    };
    Json::Value request;
    request["sequence"] = random_seq(100, 1);
    request["labels"] = names(4, "l");
    request["select"]["policy"] = "explicit";
    request["select"]["seeds"].append(seed(4));
    request["select"]["seeds"].append(seed(4));
    EXPECT_EQ("", refusal_of(request, 4));
    EXPECT_NO_THROW(cli::parse_resolve_request(request, 4));
    request["labels"] = names(5, "l");
    EXPECT_EQ("request.labels: at most 4 labels (--resolve-max-explicit-labels)",
              refusal_of(request, 4));
    EXPECT_EQ("", refusal_of(request, 5));
    EXPECT_EQ("", refusal_of(request, 0));
    request["labels"] = names(4, "l");
    request["select"]["seeds"][1] = seed(5);
    EXPECT_EQ("request.select.seeds[1].labels: at most 4 labels per seed "
              "(--resolve-max-explicit-labels)", refusal_of(request, 4));
    // `labels` first, then the seeds in order
    request["select"]["seeds"][0] = seed(6);
    EXPECT_EQ("request.select.seeds[0].labels: at most 4 labels per seed "
              "(--resolve-max-explicit-labels)", refusal_of(request, 4));
    request["labels"] = names(9, "l");
    EXPECT_EQ("request.labels: at most 4 labels (--resolve-max-explicit-labels)",
              refusal_of(request, 4));
    // before every other fault of the request's fields
    request["bogus"] = 1;
    request.removeMember("sequence");
    EXPECT_EQ("request.labels: at most 4 labels (--resolve-max-explicit-labels)",
              refusal_of(request, 4));
    // a list of another type is the parse's to refuse, in its words
    Json::Value typed;
    typed["sequence"] = random_seq(100, 1);
    typed["labels"] = "l0";
    EXPECT_EQ(refusal_of(typed, 0), refusal_of(typed, 1));
    EXPECT_NE("", refusal_of(typed, 1));
    // the route reads its cap from the limits: refused before the index is read
    const Records r = make_records<DBGSuccinct, annot::ColumnCompressed<>>(4, 2, 808);
    Json::Value over;
    over["sequence"] = r.q;
    over["labels"] = names(3, "C");
    cli::ResolveTimeLimits limits;
    limits.max_explicit_labels = 2;
    try {
        cli::process_resolve_request(over, *r.anno, "r1", 0, nullptr, nullptr, limits);
        ADD_FAILURE() << "not refused";
    } catch (const cli::InvalidRequest &e) {
        EXPECT_EQ("request.labels: at most 2 labels (--resolve-max-explicit-labels)",
                  std::string(e.what()));
    }
}

// The explicit seeds' names are looked up by hash (review 5 #6, the H1 review's S3): a pair of a
// multi-graph server, 100,000 names it does not hold beside one it holds, and an explicit seed
// naming all of them, answered in well under a second each way -- the selection made from the
// profile (select_seeds), and under a stop the check of what the prefix can tell
// (explicit_selection_blocked) -- where a linear search per name took minutes; the names covering
// nothing are listed in the request's order. And the single-graph check before the work
// (check_explicit_seeds) over 100,000 listed names
TEST(ResolveBoundsRoute, ExplicitSeedNamesAreLookedUpByHash) {
    const std::string q = random_seq(200, 909);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, { q }, { "held" });
    const size_t n = 100'000;
    Json::Value request;
    request["sequence"] = q;
    Json::Value seed;
    seed["kmer_interval"].append(0);
    seed["kmer_interval"].append(10);
    for (size_t i = 0; i < n; ++i) {
        const std::string name = "absent" + std::to_string(i);
        request["labels"].append(name);
        seed["labels"].append(name);
    }
    request["labels"].append("held");
    seed["labels"].append("held");
    request["select"]["policy"] = "explicit";
    request["select"]["seeds"].append(seed);
    auto timed = [](const std::function<Json::Value()> &f, double *seconds) {
        const auto start = std::chrono::steady_clock::now();
        Json::Value out = f();
        *seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        return out;
    };
    double seconds = 0;
    const Json::Value out = timed([&]() {
        return cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, {},
                                            nullptr, nullptr, /* drop */ true);
    }, &seconds);
    EXPECT_LT(seconds, 10.0) << "select_seeds over " << n << " names";
    ASSERT_EQ(1u, out["selection"]["seeds"].size());
    const Json::Value &picked = out["selection"]["seeds"][0];
    EXPECT_EQ("[\"held\"]", compact(picked["labels"]));
    ASSERT_EQ(n, picked["labels_not_covering"].size());
    EXPECT_EQ("absent0", picked["labels_not_covering"][0].asString());
    EXPECT_EQ("absent99999", picked["labels_not_covering"][static_cast<Json::ArrayIndex>(n - 1)].asString());
    EXPECT_EQ(n, out["unknown_labels"].size());

    // under a stop (the deadline passed before the first row): the prefix cannot tell the seed
    Json::Value budgeted = request;
    budgeted["bounds"]["time_budget_ms"] = 251;
    auto base = std::chrono::steady_clock::now();
    uint64_t reads = 0;
    std::function<std::chrono::steady_clock::time_point()> clock = [&]() {
        return base + std::chrono::milliseconds(reads++ ? 10 : 0);
    };
    const Json::Value stopped = timed([&]() {
        return cli::process_resolve_request(budgeted, *anno, "r1", 0, nullptr, nullptr,
                                            { 0, cli::kResolveFinalizeMs }, nullptr, clock,
                                            /* drop */ true);
    }, &seconds);
    EXPECT_LT(seconds, 10.0) << "explicit_selection_blocked over " << n << " names";
    ASSERT_TRUE(stopped["stop"].isObject()) << compact(stopped).substr(0, 300);
    EXPECT_TRUE(stopped["selection"].isNull());

    // the check before the work: the request lists every name a seed names
    ResolveOptions options;
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::EXPLICIT;
    ExplicitSeed ex;
    ex.kmers = { 0, 10 };
    for (size_t i = 0; i < n; ++i) {
        options.labels.push_back("l" + std::to_string(i));
    }
    ex.labels = options.labels;
    std::reverse(ex.labels.begin(), ex.labels.end());
    policy.explicit_seeds.push_back(ex);
    LabelOracle oracle(*anno);
    const auto start = std::chrono::steady_clock::now();
    EXPECT_NO_THROW(check_explicit_seeds(oracle, q, options, policy));
    EXPECT_LT(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(),
              10.0) << "check_explicit_seeds over " << n << " names";
    // a name the request does not list is still the one refused, the first in the seed
    policy.explicit_seeds[0].labels.push_back("unlisted");
    policy.explicit_seeds[0].labels.push_back("unlisted2");
    try {
        check_explicit_seeds(oracle, q, options, policy);
        ADD_FAILURE() << "not refused";
    } catch (const std::invalid_argument &e) {
        EXPECT_EQ("Explicit seed label was not profiled: 'unlisted'", std::string(e.what()));
    }
}

// The reference samples of the bounds (SPEC-labeled-traversal-core.md §4.1, §4.6), real runs: 30 a cut by the label cap
// (labels_capped with the members the route writes, ranking_exact checkable from the entry alone), 31 a memory
// stop without a deadline (the members of a time stop's block, reason memory and max_memory_mb, the entry the
// resolve of its first resolved_kmers k-mers)
TEST(ResolveBoundsRoute, Samples30To32) {
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

    // 32: the explicit labels' cap, the route's own text for the sample's request and cap
    const Json::Value s32 = sample("32_resolve_400_explicit_labels.json");
    EXPECT_EQ(400, s32["status"].asInt());
    EXPECT_EQ((std::vector<std::string>{ "error" }), s32["body"].getMemberNames());
    EXPECT_EQ(s32["body"]["error"].asString(), refusal_of(s32["request"], 4));
    EXPECT_EQ("", refusal_of(s32["request"], 5));
}


// ---- the cap on a request's explicit seeds (--resolve-max-explicit-seeds; GPT review 7 #1)

// A request of |n| explicit seeds on k-mer 0, each naming |label| (the request's one label)
Json::Value seeds_request(const std::string &sequence, const std::string &label, size_t n) {
    Json::Value request;
    request["sequence"] = sequence;
    request["labels"].append(label);
    request["select"]["policy"] = "explicit";
    Json::Value &seeds = request["select"]["seeds"];
    seeds = Json::Value(Json::arrayValue);
    for (size_t i = 0; i < n; ++i) {
        Json::Value seed;
        seed["kmer_interval"].append(0);
        seed["kmer_interval"].append(1);
        seed["labels"].append(label);
        seeds.append(seed);
    }
    return request;
}

std::string seeds_refusal_of(const Json::Value &request, uint64_t labels_cap, uint64_t seeds_cap) {
    try {
        cli::parse_resolve_request(request, labels_cap, seeds_cap);
    } catch (const cli::TooManyExplicitSeeds &e) {
        EXPECT_EQ(seeds_cap, e.max_seeds());
        return std::string("seeds: ") + e.what();
    } catch (const cli::InvalidRequest &e) {
        return e.what();
    }
    return "";
}

// More explicit seeds than the cap are refused before anything of the request is copied: after
// the request's `labels`, before each seed's labels and every other fault of the request's
// fields; at the cap the request is parsed; 0 is no cap. The 400's body: error, code, the cap
TEST(ResolveBoundsRoute, TheExplicitSeedsCap) {
    const std::string q = random_seq(200, 1717);
    Json::Value request = seeds_request(q, "held", 4);
    EXPECT_EQ("", seeds_refusal_of(request, 10, 4));
    EXPECT_EQ("", seeds_refusal_of(request, 10, 0));
    EXPECT_EQ("", seeds_refusal_of(request, 0, 0));
    request = seeds_request(q, "held", 5);
    const std::string refused
            = "seeds: request.select.seeds: at most 4 seeds (--resolve-max-explicit-seeds)";
    EXPECT_EQ(refused, seeds_refusal_of(request, 10, 4));
    EXPECT_EQ(refused, seeds_refusal_of(request, 0, 4));
    EXPECT_EQ("", seeds_refusal_of(request, 10, 5));
    EXPECT_EQ("", seeds_refusal_of(request, 10, 0));
    // an InvalidRequest (the CLI's exit status 1; it sets no cap)
    try {
        cli::parse_resolve_request(request, 0, 4);
        ADD_FAILURE() << "not refused";
    } catch (const cli::InvalidRequest &e) {
        EXPECT_EQ("request.select.seeds: at most 4 seeds (--resolve-max-explicit-seeds)",
                  std::string(e.what()));
    }
    // the request's `labels` first, then the count, then each seed's labels
    Json::Value labelled = request;
    labelled["labels"].append("a");
    labelled["labels"].append("b");
    EXPECT_EQ("request.labels: at most 2 labels (--resolve-max-explicit-labels)",
              seeds_refusal_of(labelled, 2, 4));
    labelled = request;
    for (const char *name : { "a", "b", "c" }) {
        labelled["select"]["seeds"][2]["labels"].append(name);
    }
    EXPECT_EQ(refused, seeds_refusal_of(labelled, 2, 4));
    EXPECT_EQ("request.select.seeds[2].labels: at most 2 labels per seed "
              "(--resolve-max-explicit-labels)", seeds_refusal_of(labelled, 2, 5));
    // before every other fault of the request's fields, whatever the policy
    Json::Value faulty = request;
    faulty["bogus"] = 1;
    faulty.removeMember("sequence");
    faulty["select"]["policy"] = "max_support";
    EXPECT_EQ(refused, seeds_refusal_of(faulty, 0, 4));
    // seeds of another type are the parse's to refuse, in its words
    Json::Value typed = request;
    typed["select"]["seeds"] = "s";
    EXPECT_EQ(seeds_refusal_of(typed, 0, 0), seeds_refusal_of(typed, 0, 1));
    EXPECT_NE("", seeds_refusal_of(typed, 0, 1));

    // the 400's body, written error first (as the intake refusals)
    const cli::TooManyExplicitSeeds e(1024);
    const Json::Value body = cli::too_many_explicit_seeds_json(e);
    EXPECT_EQ((std::vector<std::string>{ "code", "error", "max_explicit_seeds" }),
              body.getMemberNames());
    EXPECT_EQ("{\"error\":\"request.select.seeds: at most 1024 seeds "
              "(--resolve-max-explicit-seeds)\",\"code\":\"too_many_explicit_seeds\","
              "\"max_explicit_seeds\":1024}", cli::error_first_json_text(body));
    EXPECT_TRUE(body["max_explicit_seeds"].isUInt64());

    // the route reads its cap from the limits: refused before the index is read; at the cap
    // every seed is selected
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, { q }, { "held" });
    cli::ResolveTimeLimits limits;
    limits.max_explicit_seeds = 4;
    try {
        cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, limits);
        ADD_FAILURE() << "not refused";
    } catch (const cli::TooManyExplicitSeeds &e) {
        EXPECT_EQ(4u, e.max_seeds());
    }
    const Json::Value at_cap = cli::process_resolve_request(seeds_request(q, "held", 4), *anno,
                                                            "r1", 0, nullptr, nullptr, limits);
    EXPECT_EQ(4u, at_cap["selection"]["counts"]["selected"].asUInt64());
    // (a cap below select.max_seeds' default, 10, lowers the default too, stated as any clamp)
    EXPECT_EQ("{\"clamped\":[{\"effective\":4,\"field\":\"select.max_seeds\","
              "\"requested\":10}]}", compact(at_cap["limits"]));
    // 0: none (the CLI)
    EXPECT_EQ(5u, cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, {})
                          ["selection"]["counts"]["selected"].asUInt64());
}

// A select.max_seeds above the cap is lowered to it and stated in limits.clamped (as the
// answer-labels clamp of discover.max_labels), with or without a deadline; the answer is that of
// the cap; at the cap, or with none, nothing is stated
TEST(ResolveBoundsRoute, MaxSeedsIsClampedToTheSeedsCap) {
    const Records r = make_records<DBGSuccinct, annot::ColumnCompressed<>>(40, 5, 717);
    const auto index = with_mapping(r);
    Json::Value request;
    request["sequence"] = r.q;
    request["discover"]["kind"] = "header";
    request["select"]["policy"] = "longest_first";
    request["select"]["merge_overlapping"] = false;
    request["select"]["max_seeds"] = 50;
    // the request picks more than 2 seeds without the cap
    const Json::Value uncapped = cli::process_resolve_request(request, *index, "r1", 0, nullptr,
                                                              nullptr, {});
    ASSERT_GT(uncapped["selection"]["counts"]["selected"].asUInt64(), 2u);
    EXPECT_FALSE(uncapped.isMember("limits"));
    cli::ResolveTimeLimits limits;
    limits.max_explicit_seeds = 2;
    const Json::Value clamped = cli::process_resolve_request(request, *index, "r1", 0, nullptr,
                                                             nullptr, limits);
    EXPECT_EQ("{\"clamped\":[{\"effective\":2,\"field\":\"select.max_seeds\",\"requested\":50}]}",
              compact(clamped["limits"]));
    EXPECT_EQ(2u, clamped["selection"]["counts"]["selected"].asUInt64());
    EXPECT_EQ(2u, clamped["selection"]["policy"]["max_seeds"].asUInt64());
    EXPECT_FALSE(clamped.isMember("stop"));
    // the answer of max_seeds 2
    Json::Value two = request;
    two["select"]["max_seeds"] = 2;
    Json::Value plain = cli::process_resolve_request(two, *index, "r1", 0, nullptr, nullptr,
                                                     limits);
    EXPECT_FALSE(plain.isMember("limits"));
    Json::Value got = clamped;
    got.removeMember("limits");
    EXPECT_EQ(compact(untimed(plain)), compact(untimed(got)));
    // with a deadline and the answer-labels clamp: discover.max_labels, select.max_seeds, then
    // the budget
    Json::Value timed = request;
    timed["bounds"]["time_budget_ms"] = 100000;
    timed["discover"]["max_labels"] = 50;
    limits.max_time_ms = 60000;
    limits.max_discover_answer_labels = 7;
    const Json::Value all = cli::process_resolve_request(timed, *index, "r1", 0, nullptr,
                                                         nullptr, limits);
    EXPECT_EQ("[{\"effective\":7,\"field\":\"discover.max_labels\",\"requested\":50},"
              "{\"effective\":2,\"field\":\"select.max_seeds\",\"requested\":50},"
              "{\"effective\":60000,\"field\":\"bounds.time_budget_ms\",\"requested\":100000}]",
              compact(all["limits"]["clamped"]));
    // none at the cap or without it, nor without a selection
    limits = {};
    limits.max_explicit_seeds = 50;
    EXPECT_FALSE(cli::process_resolve_request(request, *index, "r1", 0, nullptr, nullptr, limits)
                         .isMember("limits"));
    limits.max_explicit_seeds = 2;
    Json::Value unselected = request;
    unselected.removeMember("select");
    EXPECT_FALSE(cli::process_resolve_request(unselected, *index, "r1", 0, nullptr, nullptr,
                                              limits).isMember("limits"));
}

// The selection reads the answer's deadline every kResolveCheckLabels steps (select_seeds'
// |check|): explicit seeds frozen and their overlapping pairs listed, and each max_support round;
// a check that throws ends it there, after at most kResolveCheckLabels steps more; the selection
// does not depend on the check
TEST(ResolveBounds, TheSelectionReadsTheDeadline) {
    const uint64_t num_kmers = 10000;
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = num_kmers;
    profile.graph_runs = { { 0, num_kmers } };
    LabelProfile lp;
    lp.label.kind = LabelKind::COLUMN;
    lp.label.name = "L";
    lp.runs = { { 0, num_kmers } };
    lp.kmers_supported = num_kmers;
    profile.labels.push_back(lp);
    const std::string query(num_kmers + kK - 1, 'A');
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::EXPLICIT;
    const size_t n = 1024;
    for (size_t i = 0; i < n; ++i) {
        policy.explicit_seeds.push_back({ { 0, 1 }, { "L" } });
    }
    const SeedSelection plain = select_seeds(profile, query, policy, false);
    ASSERT_EQ(n, plain.seeds.size());
    EXPECT_EQ(n - 1, plain.seeds[0].overlaps_with.size());
    size_t calls = 0;
    const SeedSelection checked = select_seeds(profile, query, policy, false,
                                               [&]() { ++calls; });
    // n (n - 1) pairs, n seeds frozen and n named: about n^2 / 4096 readings
    EXPECT_GE(calls, n * n / kResolveCheckLabels);
    EXPECT_LE(calls, (n * n + 4 * n) / kResolveCheckLabels + 1);
    ASSERT_EQ(n, checked.seeds.size());
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(plain.seeds[i].seed_id, checked.seeds[i].seed_id);
        EXPECT_EQ(plain.seeds[i].overlaps_with, checked.seeds[i].overlaps_with);
    }
    // a deadline that passes stops it at that reading
    struct Late {};
    calls = 0;
    EXPECT_THROW(select_seeds(profile, query, policy, false, [&]() {
        if (++calls == 3)
            throw Late();
    }), Late);
    EXPECT_EQ(3u, calls);

    // each max_support round sweeps every run (here 5,000 runs of one label), and reads the
    // deadline inside it: a step a run clipped, a candidate start swept (about 10,000), the runs
    // before each of its sorts, a label and its runs for the winner -- about 25,000 steps, six
    // readings a round
    LabelProfile gapped;
    gapped.label.kind = LabelKind::COLUMN;
    gapped.label.name = "G";
    for (uint64_t b = 0; b + 1 < num_kmers; b += 2) {
        gapped.runs.push_back({ b, b + 1 });
    }
    gapped.kmers_supported = gapped.runs.size();
    profile.labels = { gapped };
    SelectionPolicy rounds;
    rounds.policy = SelectionPolicy::MAX_SUPPORT;
    rounds.max_seeds = 5;
    rounds.min_block_bp = kK;
    calls = 0;
    const SeedSelection swept = select_seeds(profile, query, rounds, false, [&]() { ++calls; });
    EXPECT_EQ(5u, swept.seeds.size());
    EXPECT_GE(calls, 5u * 5);
    EXPECT_LE(calls, 5u * 7);
    calls = 0;
    EXPECT_THROW(select_seeds(profile, query, rounds, false, [&]() {
        if (++calls == 1)
            throw Late();
    }), Late);
}

// The max_support selection of the base (6872dd1fd), verbatim but for the deadline: each round
// clips every run by every pick so far, then sweeps the candidate starts (the reference of
// TheMaxSupportClipIsTheBases)
std::vector<std::pair<KmerInterval, std::vector<LabelId>>>
base_max_support(const SupportProfile &profile, size_t max_seeds, uint64_t min_kmers) {
    struct Run { uint64_t begin, end; };
    std::vector<Run> runs;
    for (const auto &lp : profile.labels) {
        for (const auto &run : lp.runs) {
            if (run.size() >= min_kmers)
                runs.push_back({ run.begin, run.end });
        }
    }
    std::sort(runs.begin(), runs.end(), [](const Run &a, const Run &b) { return a.begin < b.begin; });
    const std::vector<Run> all_runs = runs;
    std::vector<KmerInterval> taken;
    std::vector<std::pair<KmerInterval, std::vector<LabelId>>> picked;
    while (picked.size() < max_seeds) {
        runs.clear();
        for (const Run &run : all_runs) {
            uint64_t begin = run.begin;
            for (const auto &t : taken) {
                if (t.end <= begin || t.begin >= run.end)
                    continue;
                if (t.begin > begin && t.begin - begin >= min_kmers)
                    runs.push_back({ begin, t.begin });
                begin = std::max(begin, t.end);
            }
            if (begin < run.end && run.end - begin >= min_kmers)
                runs.push_back({ begin, run.end });
        }
        if (runs.empty())
            break;
        // the best interval by brute force over the candidate starts and every end: the most
        // runs containing it, then longer, then the smaller start
        std::vector<uint64_t> starts;
        for (const Run &run : runs) {
            starts.push_back(run.begin);
            if (run.end >= min_kmers)
                starts.push_back(run.end - min_kmers);
        }
        std::sort(starts.begin(), starts.end());
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        std::optional<KmerInterval> best;
        size_t best_count = 0;
        for (uint64_t a : starts) {
            if (a + min_kmers > profile.num_kmers)
                break;
            size_t count = 0;
            uint64_t b = std::numeric_limits<uint64_t>::max();
            for (const Run &run : runs) {
                if (run.begin <= a && run.end >= a + min_kmers) {
                    ++count;
                    b = std::min(b, run.end);
                }
            }
            if (!count)
                continue;
            const KmerInterval iv { a, b };
            if (!best || count > best_count
                    || (count == best_count && (iv.size() > best->size()
                            || (iv.size() == best->size() && iv.begin < best->begin)))) {
                best = iv;
                best_count = count;
            }
        }
        if (!best)
            break;
        std::vector<LabelId> labels;
        for (LabelId l = 0; l < profile.labels.size(); ++l) {
            for (const auto &run : profile.labels[l].runs) {
                if (run.begin <= best->begin && run.end >= best->end) {
                    labels.push_back(l);
                    break;
                }
            }
        }
        taken.push_back(*best);
        std::sort(taken.begin(), taken.end());
        picked.emplace_back(*best, std::move(labels));
    }
    return picked;
}

// A random profile: |labels| labels, each with runs of 1 to |longest| k-mers separated by gaps
// of 1 to |gap| k-mers, over |num_kmers| k-mers
SupportProfile random_runs_profile(size_t labels, uint64_t num_kmers, uint64_t longest,
                                   uint64_t gap, uint32_t seed) {
    std::mt19937 gen(seed);
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = num_kmers;
    profile.graph_runs = { { 0, num_kmers } };
    for (size_t l = 0; l < labels; ++l) {
        LabelProfile lp;
        lp.label.kind = LabelKind::COLUMN;
        lp.label.name = "L" + std::to_string(l);
        uint64_t at = gen() % (gap + 1);
        while (at < num_kmers) {
            const uint64_t end = std::min(num_kmers, at + 1 + gen() % longest);
            lp.runs.push_back({ at, end });
            lp.kmers_supported += end - at;
            at = end + 1 + gen() % gap;
        }
        profile.labels.push_back(std::move(lp));
    }
    return profile;
}

// The clip by a binary search into the picks (T4-L1) picks what the base's scan of every pick
// picked: the same intervals with the same labels, on random profiles of many labels and runs
TEST(ResolveBounds, TheMaxSupportClipIsTheBases) {
    size_t compared = 0;
    for (uint32_t seed = 0; seed < 30; ++seed) {
        const uint64_t num_kmers = 40 + seed * 3;
        const SupportProfile profile = random_runs_profile(1 + seed % 6, num_kmers,
                                                           1 + seed % 13, 1 + seed % 5, seed);
        const std::string query(num_kmers + kK - 1, 'A');
        for (uint64_t min_bp : { kK, kK + 3, kK + 20 }) {
            SelectionPolicy policy;
            policy.policy = SelectionPolicy::MAX_SUPPORT;
            policy.max_seeds = 1000;
            policy.min_block_bp = min_bp;
            const SeedSelection got = select_seeds(profile, query, policy, false);
            const auto expected = base_max_support(profile, 1000, min_bp - kK + 1);
            ASSERT_EQ(expected.size(), got.seeds.size()) << seed << " " << min_bp;
            for (size_t i = 0; i < expected.size(); ++i) {
                EXPECT_EQ(expected[i].first, got.seeds[i].kmers) << seed << " " << i;
                std::vector<std::string> names;
                for (LabelId l : expected[i].second) {
                    names.push_back(profile.labels[l].label.name);
                }
                std::sort(names.begin(), names.end());
                EXPECT_EQ(names, got.seeds[i].labels) << seed << " " << i;
            }
            compared += expected.size();
        }
    }
    EXPECT_GT(compared, 500u);
}

// A tiny time budget stops a max_support selection of many runs inside its first round (V7-A:
// the rounds are at most --resolve-max-explicit-seeds, but a round sweeps every run): the
// deadline is read about once every kResolveCheckLabels runs swept, so a deadline passing at the
// selection's second reading ends it within a few thousand steps, a small part of one round;
// before, it was read once a round, after every run was clipped by every pick and swept
TEST(ResolveBounds, AManyRunMaxSupportSelectionStopsInTime) {
    // 400 labels of about 2,500 runs each: about a million runs
    const SupportProfile profile = random_runs_profile(400, 10000, 2, 2, 2024);
    size_t num_runs = 0;
    for (const LabelProfile &lp : profile.labels) {
        num_runs += lp.runs.size();
    }
    ASSERT_GT(num_runs, 800'000u);
    const std::string query(profile.num_kmers + kK - 1, 'A');
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.min_block_bp = kK;
    auto timed = [&](size_t max_seeds, const std::function<void()> &check) {
        policy.max_seeds = max_seeds;
        const auto t0 = std::chrono::steady_clock::now();
        try {
            select_seeds(profile, query, policy, false, check);
        } catch (...) {}
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    // no round: the runs gathered and sorted (read before the sort); three rounds, every reading
    // counted: at least one per kResolveCheckLabels runs clipped in each round
    const double setup = timed(0, {});
    size_t calls = 0;
    const double three_rounds = timed(3, [&]() { ++calls; });
    const double round = (three_rounds - setup) / 3;
    EXPECT_GE(calls, 3 * num_runs / kResolveCheckLabels);
    // the deadline passes at the third reading, the first inside the first round: stopped
    // there, a small part of the round (before, the next reading came after the whole round)
    struct Late {};
    calls = 0;
    const double stopped = timed(1024, [&]() {
        if (++calls == 3)
            throw Late();
    });
    EXPECT_EQ(3u, calls);
    EXPECT_LT((stopped - setup) * 8, round) << stopped - setup << " s into the first round of "
                                            << round << " s";
}

// The route (process_resolve_request): under bounds.time_budget_ms the selection of as many
// explicit seeds as the cap (1024, every pair overlapping) reads the deadline about 256 times, and
// a deadline that passes while it runs ends it at the next reading, the 503 deadline of an answer
// that could not be built in time (nothing partial); before, the selection was read once it was
// done, after all of its 1024 x 1023 overlaps were listed
TEST(ResolveBoundsRoute, ATimeStopInTheSelectionIsTheDeadline) {
    const std::string q = random_seq(200, 1818);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, { q }, { "held" });
    const cli::ResolveTimeLimits limits { 0, cli::kResolveFinalizeMs, 0, 0, 0, 0, 1024 };
    const auto base = std::chrono::steady_clock::now();
    uint64_t reads = 0;
    uint64_t late_from = std::numeric_limits<uint64_t>::max();
    std::function<std::chrono::steady_clock::time_point()> clock = [&]() {
        return base + std::chrono::hours(reads++ >= late_from ? 1 : 0);
    };
    auto run = [&](size_t seeds) {
        Json::Value request = seeds_request(q, "held", seeds);
        request["bounds"]["time_budget_ms"] = 1000;
        reads = 0;
        return cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, limits,
                                            nullptr, clock);
    };
    // within the budget: the readings of one seed and of 1024
    const Json::Value one = run(1);
    const uint64_t reads_one = reads;
    ASSERT_TRUE(one["stop"].isNull());
    ASSERT_EQ(1u, one["selection"]["counts"]["selected"].asUInt64());
    const Json::Value all = run(1024);
    const uint64_t reads_all = reads;
    ASSERT_EQ(1024u, all["selection"]["counts"]["selected"].asUInt64());
    EXPECT_GE(reads_all - reads_one, 1024u * 1023u / kResolveCheckLabels);
    // the deadline passes at the third reading inside the selection (the work before it reads
    // the clock as one seed's does, which reads it twice after its selection): the 503 there
    late_from = reads_one;
    try {
        run(1024);
        ADD_FAILURE() << "not stopped";
    } catch (const cli::ResolveDeadline &e) {
        EXPECT_EQ("deadline", cli::resolve_deadline_body(e)["code"].asString());
    }
    // stopped at the reading that saw it, far before the selection's end
    EXPECT_EQ(late_from + 1, reads);
    EXPECT_LT(reads + 200, reads_all);
}

// Sample 38 (SPEC-labeled-traversal-core.md §4.1), a real run: more explicit seeds than the cap
// (2 on the recording server), the route's own body for the sample's request and cap
TEST(ResolveBoundsRoute, Sample38) {
    const std::string path
            = source_path("tests/data/multigraph_contract/38_resolve_400_explicit_seeds.json");
    std::ifstream in(path, std::ios::binary);
    ASSERT_TRUE(in.good()) << path << " not found";
    Json::Value doc;
    std::string errors;
    Json::CharReaderBuilder builder;
    ASSERT_TRUE(Json::parseFromStream(builder, in, &doc, &errors)) << errors;
    EXPECT_EQ(400, doc["status"].asInt());
    const uint64_t cap = doc["body"]["max_explicit_seeds"].asUInt64();
    EXPECT_EQ(2u, cap);
    EXPECT_GT(doc["request"]["select"]["seeds"].size(), cap);
    EXPECT_EQ(compact(cli::too_many_explicit_seeds_json(cli::TooManyExplicitSeeds(cap))),
              compact(doc["body"]));
    EXPECT_EQ("seeds: " + doc["body"]["error"].asString(),
              seeds_refusal_of(doc["request"], 10'000, cap));
    EXPECT_EQ("", seeds_refusal_of(doc["request"], 10'000, cap + 1));
}

// ---- the explicit labels' pass and the answer account (T4-L1, GPT review 7, finding 2)

// the names a discovery of every label of |c|'s kind finds on |r|'s query: the explicit labels
std::vector<std::string> all_labels(LabelOracle &oracle, const std::string &q, const Case &c) {
    ResolveOptions opts = discovery(c, kAll);
    return names_of(resolve_support(oracle, q, opts));
}

ResolveOptions explicit_labels(const Case &c, const std::vector<std::string> &labels) {
    ResolveOptions opts;
    opts.labels = labels;
    opts.support = c.support;
    return opts;
}

// The explicit labels' pass on the row paths is a discovery's (for_each_row, charged to the
// memory bound): under a bound that binds it stops at a k-mer, the answer exactly the resolve of
// the query's prefix before it (reason memory, phase rows), the k-mer never moving back as the
// bound grows; without a stop the profile is the unbounded one. Before, the explicit labels' rows
// and hits were outside any account and never stopped
TYPED_TEST(ResolveBoundsTest, ExplicitLabelsStopAtAPrefix) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 1919);
    LabelOracle oracle(*r.anno, r.cth.get());
    size_t stops = 0, cannot = 0;
    for (const Case &c : kCases) {
        const std::vector<std::string> labels = all_labels(oracle, r.q, c);
        ASSERT_GT(labels.size(), 16u) << c.name() << ": not the row paths";
        for (size_t rows : { size_t(7), kResolveBatchRows }) {
            ResolveOptions opts = explicit_labels(c, labels);
            opts.batch_rows = rows;
            const SupportProfile plain = resolve_support(oracle, r.q, opts);
            ASSERT_FALSE(plain.stop);
            uint64_t last = 0;
            bool ended = false;
            for (uint64_t bound = 1024; !ended && bound < (uint64_t(64) << 20);
                    bound += bound < 65536 ? 1024 : bound / 8) {
                const std::string what = c.name() + " rows " + std::to_string(rows) + " bound "
                        + std::to_string(bound);
                ResolveOptions bounded = opts;
                bounded.discover_max_memory = bound;
                const SupportProfile got = resolve_support(oracle, r.q, bounded);
                if (!got.stop) {
                    expect_same_profile(plain, got, what);
                    ended = true;
                    break;
                }
                ++stops;
                EXPECT_EQ(ResolveStop::MEMORY, got.stop->reason) << what;
                EXPECT_EQ(ResolveStop::ROWS, got.stop->phase) << what;
                EXPECT_EQ(bound, got.stop->max_memory) << what;
                const uint64_t x = got.stop->resolved_kmers;
                EXPECT_EQ(x, got.num_kmers) << what;
                EXPECT_GE(x, last) << what;
                last = x;
                // every explicit label is profiled, the prefix's support
                ASSERT_EQ(labels.size(), got.labels.size()) << what;
                const uint64_t first = got.stop->query_graph_runs.front().begin;
                EXPECT_EQ(x == first, got.stop->cannot_progress) << what;
                cannot += got.stop->cannot_progress;
                if (x == 0)
                    continue;
                const SupportProfile prefix
                        = resolve_support(oracle, r.q.substr(0, x + kK - 1), opts);
                expect_same_profile(prefix, got, what);
            }
            EXPECT_TRUE(ended) << c.name();
        }
    }
    EXPECT_GT(stops, 50u);
    EXPECT_GT(cannot, 0u);
}

// The explicit labels' pass holds one k-mer's hits at a time: under a deadline that never passes
// and without one, and whatever the batches and the rows kept, the profile is the same; the
// smallest bound that answers the whole query answers it as without a bound
TYPED_TEST(ResolveBoundsTest, ExplicitLabelsHoldOneKmersHits) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const Records r = make_records<Graph, Annotation>(40, 5, 2020);
    LabelOracle oracle(*r.anno, r.cth.get());
    for (const Case &c : kCases) {
        const std::vector<std::string> labels = all_labels(oracle, r.q, c);
        const SupportProfile plain = resolve_support(oracle, r.q, explicit_labels(c, labels));
        for (size_t rows : { size_t(1), size_t(3), size_t(64) }) {
            for (uint64_t kept : { uint64_t(0), kResolveKeptBytes }) {
                ResolveOptions opts = explicit_labels(c, labels);
                opts.batch_rows = rows;
                opts.kept_bytes = kept;
                opts.time_up = []() { return false; };
                const std::string what = c.name() + " rows " + std::to_string(rows);
                expect_same_profile(plain, resolve_support(oracle, r.q, opts), what);
            }
        }
        // the smallest bound that answers the whole query, bisected: above it every bound does
        uint64_t lo = 1, hi = uint64_t(64) << 20;
        while (lo + 1 < hi) {
            ResolveOptions opts = explicit_labels(c, labels);
            opts.discover_max_memory = (lo + hi) / 2;
            (resolve_support(oracle, r.q, opts).stop ? lo : hi) = opts.discover_max_memory;
        }
        ResolveOptions opts = explicit_labels(c, labels);
        opts.discover_max_memory = hi;
        const SupportProfile got = resolve_support(oracle, r.q, opts);
        ASSERT_FALSE(got.stop) << c.name();
        expect_same_profile(plain, got, c.name());
    }
}

// The answer's JSON is priced before each part is built (profile_to_json's answer account): an
// account smaller than the price refuses it at that part (AnswerTooLarge), a larger one answers
// it byte for byte as without an account, and the price covers what jsoncpp's tree and its text
// hold (the allocator's live bytes, where they can be read)
std::optional<int64_t> live_heap_bytes() {
#if USE_JEMALLOC
    uint64_t allocated = 0, deallocated = 0;
    size_t sz = sizeof(uint64_t);
    if (!mallctl("thread.allocated", &allocated, &sz, nullptr, 0)
            && !mallctl("thread.deallocated", &deallocated, &sz, nullptr, 0)
            && (allocated || deallocated)) {
        return int64_t(allocated) - int64_t(deallocated);
    }
#endif
#if defined(__APPLE__)
    malloc_statistics_t stats;
    malloc_zone_statistics(nullptr, &stats);
    return int64_t(stats.size_in_use);
#elif defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    return int64_t(mallinfo2().uordblks);
#else
    return std::nullopt;
#endif
}

// A profile of |labels| labels of |runs| runs each (names of |name_bytes| bytes), its candidates
// (a run of its own per label and one interval all share) and a selection of |seeds| seeds
// naming every label, each overlapping the others
std::pair<SupportProfile, SeedSelection> priced_profile(size_t labels, size_t runs,
                                                        size_t name_bytes, size_t seeds,
                                                        bool trace) {
    SupportProfile p;
    p.k = 31;
    p.num_kmers = 3 * runs + 2000;
    p.support = trace ? Support::TRACE : Support::KMER;
    p.graph_runs = { { 0, p.num_kmers } };
    for (size_t l = 0; l < labels; ++l) {
        LabelProfile lp;
        lp.label.kind = LabelKind::HEADER;
        lp.label.name = std::string(name_bytes, 'a' + l % 26) + std::to_string(l);
        for (size_t i = 0; i < runs; ++i) {
            lp.runs.push_back({ 3 * i + l % 2, 3 * i + 2 });
            if (trace)
                lp.trace_breaks.push_back(3 * i + 1);
        }
        lp.runs.push_back({ 3 * runs + 7, 3 * runs + 1000 });
        lp.kmers_supported = 2 * runs + 993;
        p.labels.push_back(std::move(lp));
    }
    std::map<KmerInterval, std::vector<LabelId>> groups;
    for (LabelId l = 0; l < labels; ++l) {
        for (const auto &run : p.labels[l].runs) {
            groups[run].push_back(l);
        }
    }
    for (auto &[iv, ls] : groups) {
        p.candidates.push_back({ iv, ls });
    }
    SeedSelection sel;
    for (size_t i = 0; i < seeds; ++i) {
        FrozenSeed seed;
        seed.seed_id = "0123456789abcdef";
        seed.sequence = std::string(1000, 'A');
        seed.kmers = { 3 * runs + 7, 3 * runs + 977 };
        for (const LabelProfile &lp : p.labels) {
            seed.labels.push_back(lp.label.name);
        }
        seed.population.supporting_total = labels;
        seed.population.included = labels;
        for (size_t j = 0; j < seeds; ++j) {
            if (j != i)
                seed.overlaps_with.push_back(j);
        }
        sel.seeds.push_back(std::move(seed));
    }
    return { std::move(p), std::move(sel) };
}

TEST(ResolveAnswer, ThePriceCoversTheTree) {
    size_t measured = 0;
    for (bool trace : { false, true }) {
        for (const char *format : { "intervals", "rle" }) {
            for (auto [labels, runs, name, seeds] : { std::make_tuple(2000, 3, 12, 2),
                                                       std::make_tuple(20, 2000, 40, 0),
                                                       std::make_tuple(300, 100, 100, 30) }) {
                const auto [p, sel] = priced_profile(labels, runs, name, seeds, trace);
                const std::string what = std::string(format) + (trace ? " trace " : " kmer ")
                        + std::to_string(labels) + "x" + std::to_string(runs);
                AnswerAccount account(uint64_t(1) << 40);
                const std::optional<int64_t> before = live_heap_bytes();
                uint64_t tree_and_text = 0;
                {
                    const Json::Value out = cli::profile_to_json(p, &sel, format, nullptr,
                                                                 &account);
                    const std::string text = compact(out);
                    const std::optional<int64_t> after = live_heap_bytes();
                    if (before && after && *after > *before) {
                        tree_and_text = *after - *before;
                        ++measured;
                    }
                    // the same answer as without an account
                    EXPECT_EQ(compact(cli::profile_to_json(p, &sel, format)), text) << what;
                    // the text, three times its bytes in the price, at least once
                    EXPECT_GE(account.held(), text.size()) << what;
                }
                EXPECT_GE(account.held(), tree_and_text) << what << ": priced " << account.held()
                                                         << ", the heap held " << tree_and_text;
                // an account a byte short of the price refuses it, at the answer's part
                AnswerAccount short_one(account.held() - 1);
                try {
                    cli::profile_to_json(p, &sel, format, nullptr, &short_one);
                    ADD_FAILURE() << what << ": not refused";
                } catch (const AnswerTooLarge &e) {
                    EXPECT_EQ("answer", e.part()) << what;
                }
            }
        }
    }
    if (!measured)
        GTEST_SKIP() << "no allocator counters to compare the prices with";
}

// A max_support selection's scratch is charged as its rounds allocate it (T4 review R1, F1): a
// round holds a piece for each part of a run that the picks clip, so long runs split by many
// picks hold far more than one entry a gathered run, the price before. Here 100 labels with a run
// over the whole query and three labels of short runs every 32 k-mers, every pick inside the long
// runs: at every reading of the selection's check the heap the selection added is within what its
// account holds then (a few KiB aside: the picks' and the claimed intervals' vectors), and the
// selection is the one made without an account
TEST(ResolveAnswer, TheSelectionChargesTheClippedPieces) {
    SupportProfile profile;
    profile.k = kK;
    profile.num_kmers = 4000;
    profile.graph_runs = { { 0, profile.num_kmers } };
    for (size_t l = 0; l < 100; ++l) {
        LabelProfile lp;
        lp.label.kind = LabelKind::HEADER;
        lp.label.name = "long" + std::to_string(l);
        lp.runs.push_back({ 0, profile.num_kmers });
        lp.kmers_supported = profile.num_kmers;
        profile.labels.push_back(std::move(lp));
    }
    for (uint64_t offset : { 7, 17, 27 }) {
        LabelProfile lp;
        lp.label.kind = LabelKind::HEADER;
        lp.label.name = "short" + std::to_string(offset);
        for (uint64_t s = offset; s < profile.num_kmers; s += 32) {
            lp.runs.push_back({ s, s + 1 });
        }
        lp.kmers_supported = lp.runs.size();
        profile.labels.push_back(std::move(lp));
    }
    const std::string query(profile.num_kmers + kK - 1, 'A');
    SelectionPolicy policy;
    policy.policy = SelectionPolicy::MAX_SUPPORT;
    policy.min_block_bp = kK;
    policy.max_seeds = 250;
    const SeedSelection plain = select_seeds(profile, query, policy, false);
    ASSERT_EQ(250u, plain.seeds.size());

    AnswerAccount account(uint64_t(1) << 40);
    const std::optional<int64_t> base = live_heap_bytes();
    int64_t worst = std::numeric_limits<int64_t>::min();
    int64_t heap_at_worst = 0;
    uint64_t held_at_worst = 0;
    size_t readings = 0;
    auto check = [&]() {
        const std::optional<int64_t> now = live_heap_bytes();
        if (!base || !now)
            return;
        ++readings;
        const int64_t gap = (*now - *base) - int64_t(account.held());
        if (gap > worst) {
            worst = gap;
            heap_at_worst = *now - *base;
            held_at_worst = account.held();
        }
    };
    const SeedSelection got = select_seeds(profile, query, policy, false, check, &account);
    ASSERT_EQ(plain.seeds.size(), got.seeds.size());
    for (size_t i = 0; i < got.seeds.size(); ++i) {
        EXPECT_EQ(plain.seeds[i].kmers, got.seeds[i].kmers) << i;
        EXPECT_EQ(plain.seeds[i].labels, got.seeds[i].labels) << i;
    }
    // the rounds' pieces: 100 long runs, each split by up to 250 picks
    EXPECT_GT(account.peak(), 100u * 250 * sizeof(uint64_t) * 4);
    if (!readings)
        GTEST_SKIP() << "no allocator counters to compare the prices with";
    EXPECT_GT(readings, 1000u);
    EXPECT_LE(worst, 64 << 10) << "the heap held " << heap_at_worst << " B while the account held "
                               << held_at_worst << " B";
    std::cerr << "max_support selection: at its worst reading of " << readings << " the heap held "
              << heap_at_worst << " B, the account " << held_at_worst << " B (peak "
              << account.peak() << " B)" << std::endl;
}

// The route: an answer priced past --resolve-max-discover-memory-mb is refused, 400
// resolve_answer_too_large with the bound (members error, code, max_memory_mb; the server writes
// them in that order); the same request with its runs as rle strings and no candidates of single
// k-mers, or under a larger bound, is answered as without a bound. Labels that differ from the
// query every 12 bases have about a run a dozen k-mers: an answer of about 450 B a k-mer
TEST(ResolveBoundsRoute, TheAnswerIsPricedBeforeItIsBuilt) {
    const std::string q = random_seq(8000, 2121);
    std::vector<std::string> seqs { q }, names { "A" };
    std::mt19937 gen(2122);
    for (size_t b = 0; b < 4; ++b) {
        std::string s = q;
        for (size_t i = gen() % 12; i < s.size(); i += 12) {
            s[i] = s[i] == 'A' ? 'C' : 'A';
        }
        seqs.push_back(s);
        names.push_back("B" + std::to_string(b));
    }
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, seqs, names);
    Json::Value request;
    request["sequence"] = q;
    for (const std::string &name : names) {
        request["labels"].append(name);
    }
    const cli::ResolveTimeLimits none;
    const cli::ResolveTimeLimits one_mib { 0, cli::kResolveFinalizeMs, 0, 1 };
    const cli::ResolveTimeLimits large { 0, cli::kResolveFinalizeMs, 0, 1024 };
    const Json::Value plain = cli::process_resolve_request(request, *anno, "r1", 0, nullptr,
                                                           nullptr, none);
    ASSERT_GT(compact(plain).size(), 100'000u);
    try {
        cli::process_resolve_request(request, *anno, "r1", 0, nullptr, nullptr, one_mib);
        ADD_FAILURE() << "not refused";
    } catch (const AnswerTooLarge &e) {
        const Json::Value body = cli::resolve_answer_too_large_json(e);
        EXPECT_EQ((std::vector<std::string>{ "code", "error", "max_memory_mb" }),
                  body.getMemberNames());
        EXPECT_EQ("resolve_answer_too_large", body["code"].asString());
        EXPECT_EQ(1u, body["max_memory_mb"].asUInt64());
        EXPECT_EQ("resolve: the answer would hold more than --resolve-max-discover-memory-mb 1 MiB ("
                  + e.part() + "): nothing partial is sent; ask for run_format \"rle\", a larger "
                  "min_block_kmers, fewer labels or seeds, or a shorter sequence",
                  body["error"].asString());
    }
    EXPECT_EQ(compact(untimed(plain)),
              compact(untimed(cli::process_resolve_request(request, *anno, "r1", 0, nullptr,
                                                           nullptr, large))));
    Json::Value small = request;
    small["run_format"] = "rle";
    small["min_block_kmers"] = 12;
    EXPECT_EQ(compact(untimed(cli::process_resolve_request(small, *anno, "r1", 0, nullptr,
                                                           nullptr, none))),
              compact(untimed(cli::process_resolve_request(small, *anno, "r1", 0, nullptr,
                                                           nullptr, one_mib))));
}

// The route's explicit labels on the row paths under the bound: a memory stop without a deadline,
// stated as a discovery's (reason memory, max_memory_mb, the prefix), the message the memory
// stop's
TEST(ResolveBoundsRoute, AnExplicitLabelsMemoryStopIsStated) {
    const Records r = make_records<DBGSuccinct, annot::ColumnCompressed<>>(40, 5, 2323);
    const auto index = with_mapping(r);
    LabelOracle oracle(*r.anno, r.cth.get());
    Json::Value request;
    request["sequence"] = r.q;
    for (const std::string &name : all_labels(oracle, r.q, { LabelKind::HEADER, Support::KMER })) {
        request["labels"].append(name);
    }
    // a bound of 1 MiB with the query repeated: its setup alone does not fit
    std::string q = r.q;
    while (q.size() < 20000) {
        q += r.q;
    }
    request["sequence"] = q;
    const cli::ResolveTimeLimits one_mib { 0, cli::kResolveFinalizeMs, 0, 1 };
    const Json::Value out = cli::process_resolve_request(request, *index, "r1", 0, nullptr,
                                                         nullptr, one_mib);
    ASSERT_TRUE(out["stop"].isObject()) << compact(out).substr(0, 300);
    EXPECT_EQ("memory", out["stop"]["reason"].asString());
    EXPECT_EQ("rows", out["stop"]["phase"].asString());
    EXPECT_EQ(1u, out["stop"]["max_memory_mb"].asUInt64());
    EXPECT_EQ(out["num_kmers"], out["stop"]["resolved_kmers"]);
    EXPECT_EQ(0u, out["stop"]["message"].asString().find(
            "--resolve-max-discover-memory-mb 1 ran out while the annotation rows of the query's "
            "k-mers were read: this answer is exactly the resolve of the sequence's first "));
    EXPECT_FALSE(out.isMember("limits"));
}

} // namespace
