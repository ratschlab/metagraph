/**
 * The segment search at the graph level (PatternSearch::count with Request::segments) against
 * the oracles of pattern_segment_oracles.hpp, which read the pattern, the windows and the graph
 * on their own (strings and k-mer sets, breadth-first layers):
 *  - random small graphs (k = 4 .. 11; BASIC, CANONICAL, wrapped PRIMARY), DNA, IUPAC and
 *    peptide patterns made of record pieces, spliced, mutated, on either strand, every window
 *    choice, windows and steps at their bounds, max_anchors low and high: the windows (offsets,
 *    states, anchor counts, roots, owned), the reported extents in their order, the segments
 *    grouped with their roots, and the work (edges of both arms, the bwd steps of the left arm,
 *    the states memoised, the horizon cuts) equal the oracle's exactly;
 *  - named junction cases (two exons, a short exon, the reverse strand, a record ending inside
 *    the pattern, two records, a palindrome, a mosaic) with their segments written out, and the
 *    window/step cases (S = 1, S = W, S beyond L - W as "ends", W = L, W > L, the guarantee
 *    W + S - 1 and its tightness, windows not admitted, below the floor, too many windows);
 *  - graphs without the mask answer the masked twin's segments;
 *  - stop sweeps (max_steps, the clock at every reading, the sink's memory, max_segments):
 *    the extents handed over are a prefix of the full run's, the windows before the stop as in
 *    the full run, the one it touched STOPPED, the later ones NOT_REACHED, the memory released.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"
#include "pattern_segment_oracles.hpp"
#include "pattern_test_support.hpp"

#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

// the nucleotide builds the engine serves, as in test_pattern_search.cpp
#if _DNA_GRAPH || _DNA5_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using namespace mtg::test::segments;
using mtg::test::build_graph;
using mtg::test::unbounded_deadline;

typedef DeBruijnGraph::node_index node_index;

constexpr uint64_t kManySteps = 1'000'000'000;
constexpr uint64_t kManyBytes = std::numeric_limits<uint64_t>::max() / 2;


// ---------------------------------------------------------------- running the engine

// one extent as the sink received it
struct Recorded {
    Orientation orientation;
    uint32_t window;
    node_index root;
    uint32_t from;
    uint32_t to;

    bool operator==(const Recorded &other) const {
        return std::tie(orientation, window, root, from, to)
                == std::tie(other.orientation, other.window, other.root, other.from, other.to);
    }
};

std::ostream& operator<<(std::ostream &out, const Recorded &r) {
    return out << orientation_key(r.orientation) << " w" << r.window << " root " << r.root << " ["
               << r.from << ", " << r.to << ")";
}

// records the extents in order and passes everything to a GraphSegmentStore
class TeeSink : public SegmentSink {
  public:
    explicit TeeSink(GraphSegmentStore &store) : store_(store) {}

    bool admit(uint64_t bytes) override { return store_.admit(bytes); }
    void release(uint64_t bytes) override { store_.release(bytes); }
    bool extent(const SegmentExtent &e) override {
        recorded.push_back(Recorded { e.orientation, e.window, e.root.node, e.from, e.to });
        return store_.extent(e);
    }
    const char* stop_reason() const override { return store_.stop_reason(); }

    std::vector<Recorded> recorded;

  private:
    GraphSegmentStore &store_;
};

struct SegmentRunResult {
    Result result;
    std::vector<Recorded> recorded;
    std::map<OracleSegmentKey, uint64_t> segments;
    uint64_t held = 0;
    uint64_t peak = 0;
    uint64_t store_size = 0;
    uint64_t steps_used = 0;
    std::string reason;
};

struct Limits {
    uint64_t max_steps = kManySteps;
    uint64_t max_bytes = kManyBytes;
    uint64_t max_segments = std::numeric_limits<uint64_t>::max();
    bool stop_at_threshold = false;
    // a deadline with an injected clock; unset: none
    std::function<Deadline()> deadline;
};

Request segment_request(const SegmentPlan &plan, Strands strands, uint64_t max_anchors) {
    Request request;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_anchors = max_anchors;
    request.segments = plan;
    request.segments.enabled = true;
    return request;
}

SegmentRunResult run_segments(const DeBruijnGraph &graph, const Pattern &pattern, Request request,
                 const Limits &limits = Limits()) {
    GraphSegmentStore store(static_cast<uint32_t>(pattern.length()), limits.max_bytes,
                            limits.max_segments, limits.stop_at_threshold);
    TeeSink sink(store);
    request.segment_sink = &sink;
    Budget budget(limits.max_steps,
                  limits.deadline ? limits.deadline() : unbounded_deadline());
    PatternSearch search(graph);
    SegmentRunResult run;
    run.result = search.count(pattern, request, budget);
    run.recorded = sink.recorded;
    store.for_each([&](const GraphSegment &s) {
        EXPECT_EQ(0u, run.segments.count({ s.orientation, s.from, s.to }));
        run.segments[{ s.orientation, s.from, s.to }] = s.roots;
    });
    run.held = store.held();
    run.peak = store.peak();
    run.store_size = store.size();
    run.steps_used = budget.steps_used();
    run.reason = store.stop_reason() ? store.stop_reason() : "";
    return run;
}

std::string describe_segments(const std::map<OracleSegmentKey, uint64_t> &segments) {
    std::ostringstream out;
    for (const auto &[key, roots] : segments) {
        const auto &[o, from, to] = key;
        out << orientation_key(o) << "[" << from << "," << to << ")x" << roots << " ";
    }
    return out.str();
}


// ---------------------------------------------------------------- the comparison

/**
 * Runs the engine and the oracle on one case and compares everything the graph level states:
 * windows, extents in order, segments, work, the memory released, the first windows' counts.
 */
SegmentRunResult check_against_oracle(const DeBruijnGraph &graph, PatternKind kind, const std::string &text,
                         Strands strands, const SegmentPlan &plan, uint64_t max_anchors,
                         const std::string &what) {
    const Pattern pattern = Pattern::parse(kind, text);
    const OracleGraph g = oracle_graph(graph);
    const std::vector<OracleQ> searched = oracle_orientations(kind, text, strands);
    const OracleRun e = oracle_segments(g, searched,
                                        OraclePlan { plan.window, plan.step, plan.choice,
                                                     max_anchors, 0 });
    SegmentRunResult run = run_segments(graph, pattern, segment_request(plan, strands, max_anchors));
    const Result &r = run.result;
    EXPECT_FALSE(r.refusal) << what << ": " << r.refusal->message;
    if (r.refusal)
        return run;
    EXPECT_FALSE(r.stop) << what;
    EXPECT_TRUE(r.segments) << what;
    if (!r.segments)
        return run;
    const SegmentCounts &c = *r.segments;
    EXPECT_EQ(e.window, c.window) << what;
    EXPECT_EQ(searched.size(), c.windows.size()) << what;

    uint64_t roots = 0;
    uint64_t owned = 0;
    bool extended = false;
    bool not_admitted = false;
    for (const OracleQ &q : searched) {
        const auto it = c.windows.find(q.orientation);
        EXPECT_TRUE(it != c.windows.end()) << what;
        if (it == c.windows.end())
            continue;
        const std::vector<OracleWindow> &expected = e.windows.at(q.orientation);
        const std::vector<WindowCount> &got = it->second;
        EXPECT_EQ(expected.size(), got.size()) << what;
        for (size_t i = 0; i < std::min(expected.size(), got.size()); ++i) {
            const OracleWindow &x = expected[i];
            const WindowCount &w = got[i];
            const std::string at = what + " " + orientation_key(q.orientation) + " window "
                    + std::to_string(x.offset);
            EXPECT_EQ(x.offset, w.offset) << at;
            const uint64_t n = x.anchors.size();
            if (x.admitted) {
                EXPECT_EQ(WindowState::EXTENDED, w.state) << at;
                EXPECT_EQ(Relation::EXACT, w.anchors.relation) << at;
                EXPECT_EQ(n, w.anchors.value) << at;
                extended = true;
            } else {
                EXPECT_EQ(WindowState::NOT_ADMITTED, w.state) << at;
                // its discovery stopped above max_anchors: AT_LEAST, never more than there are;
                // EXACT only where the running count can lag (even-k wrapped PRIMARY)
                const bool lags = dynamic_cast<const CanonicalDBG*>(&graph)
                        && !(graph.get_k() % 2);
                if (!lags || w.anchors.relation == Relation::AT_LEAST) {
                    EXPECT_EQ(Relation::AT_LEAST, w.anchors.relation) << at;
                    EXPECT_LT(max_anchors, w.anchors.value) << at;
                    EXPECT_GE(n, w.anchors.value) << at;
                } else {
                    EXPECT_EQ(Relation::EXACT, w.anchors.relation) << at;
                    EXPECT_EQ(n, w.anchors.value) << at;
                }
                not_admitted = true;
            }
            EXPECT_EQ(x.roots, w.roots) << at;
            EXPECT_EQ(x.owned, w.owned) << at;
            roots += x.roots;
            owned += x.owned;
        }
        // counts.anchors keeps its meaning: the first window's
        EXPECT_TRUE(r.anchors) << what;
        if (r.anchors && !got.empty()) {
            const Count &first = r.anchors->by_orientation.at(q.orientation);
            EXPECT_EQ(got.front().anchors.relation, first.relation) << what;
            EXPECT_EQ(got.front().anchors.value, first.value) << what;
        }
    }
    EXPECT_EQ(!extended && not_admitted ? Extension::NOT_ADMITTED : Extension::COMPLETED,
              c.search) << what;

    std::vector<Recorded> expected_order;
    uint64_t below = 0;
    for (const OracleExtent &x : e.extents) {
        if (x.reported) {
            expected_order.push_back(Recorded { x.orientation, x.window, x.root, x.from, x.to });
        } else {
            ++below;
        }
    }
    EXPECT_EQ(expected_order, run.recorded) << what;
    EXPECT_EQ(e.segments, run.segments) << what << "\nexpected " << describe_segments(e.segments)
                                        << "\ngot      " << describe_segments(run.segments);
    EXPECT_EQ(expected_order.size(), c.reported) << what;

    const Work &w = r.work;
    EXPECT_EQ(roots, w.segment_anchors_extended) << what;
    EXPECT_EQ(roots, w.extension_anchors) << what;
    EXPECT_EQ(owned, w.segment_anchors_owned) << what;
    EXPECT_EQ(below, w.segment_occurrences_below_window) << what;
    EXPECT_EQ(e.states, w.segment_states_memoised) << what;
    EXPECT_EQ(e.right_edges, w.extension_edges) << what;
    EXPECT_EQ(e.right_branches, w.extension_branches) << what;
    EXPECT_EQ(e.left_edges, w.left_edges) << what;
    EXPECT_EQ(e.left_cuts, w.left_horizon_cuts) << what;
    // a root's left walks never reach the previous extended window (its anchors' right arms
    // would have reached the root): the horizon guards, it never cuts
    EXPECT_EQ(0u, w.left_horizon_cuts) << what;

    const size_t k = graph.get_k();
    if (dynamic_cast<const CanonicalDBG*>(&graph)) {
        EXPECT_EQ(e.left_edges * (k - 2), w.left_backward_steps) << what;
    } else {
        // the incoming edge a node's bwd names costs 2, any other K - 1, and every root at an
        // offset > 0 its own pair, K - 1
        const auto &boss = oracle_base_dbg(graph).get_boss();
        const uint64_t K = boss.get_k();
        const uint64_t other = K >= 2 ? K - 1 : 1;
        uint64_t bwds = e.left_roots * other;
        for (const std::string &u : e.left_expanded) {
            const node_index node = g.kmers.at(u);
            for (const std::string &x : g.incoming(u)) {
                bwds += K >= 2 && g.kmers.at(x) == boss.bwd(node) ? 2 : other;
            }
        }
        EXPECT_EQ(bwds, w.left_backward_steps) << what;
    }
    EXPECT_EQ(run.steps_used, w.steps) << what;
    EXPECT_EQ(w.ranges_visited + w.mask_scans + w.extension_edges + w.left_edges
                    + w.left_backward_steps, w.steps) << what;
    EXPECT_EQ(w.segment_windows, [&]() {
        uint64_t n = 0;
        for (const auto &[o, windows] : e.windows) {
            n += windows.size();
        }
        return n;
    }()) << what;
    // the search's states released; the store holds its own segments only
    EXPECT_EQ(run.store_size * GraphSegmentStore::kGraphSegmentBytes, run.held) << what;
    EXPECT_GE(run.peak, run.held) << what;
    return run;
}


// ---------------------------------------------------------------- random cases

std::string random_dna(std::mt19937 &rng, size_t n) {
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        s.push_back(kOracleBases[rng() % 4]);
    }
    return s;
}

// a base other than |base|
char other_base(char base, std::mt19937 &rng) {
    char c;
    do {
        c = kOracleBases[rng() % 4];
    } while (c == base);
    return c;
}

// an IUPAC code admitting |base|, wider sometimes
char iupac_of(char base, std::mt19937 &rng) {
    static const std::string codes = "RYSWKMBDHVN";
    std::vector<char> options;
    for (char c : codes) {
        if (oracle_iupac(c).find(base) != std::string::npos)
            options.push_back(c);
    }
    return options[rng() % options.size()];
}

struct Case {
    size_t k;
    std::vector<std::string> records;
    DeBruijnGraph::Mode mode;
    PatternKind kind;
    std::string text;
    Strands strands;
    SegmentPlan plan;
    uint64_t max_anchors;
};

std::string describe(const Case &c) {
    std::ostringstream out;
    out << "k=" << c.k << " mode=" << c.mode << " kind=" << to_string(c.kind) << " text="
        << c.text << " strands=" << to_string(c.strands) << " W=" << c.plan.window
        << " S=" << c.plan.step << " windows=" << to_string(c.plan.choice)
        << " max_anchors=" << c.max_anchors << " records:";
    for (const std::string &r : c.records) {
        out << " " << r;
    }
    return out.str();
}

// a piece of |source| of |n| bases, or its reverse complement
std::string piece_of(const std::string &source, size_t n, std::mt19937 &rng) {
    const std::string piece = source.substr(rng() % (source.size() - n + 1), n);
    return rng() % 3 ? piece : oracle_rev_comp(piece);
}

Case random_case(std::mt19937 &rng) {
    Case c;
    static const size_t ks[] = { 4, 5, 6, 7, 9, 11 };
    c.k = ks[rng() % 6];
    const size_t records = 2 + rng() % 4;
    for (size_t i = 0; i < records; ++i) {
        c.records.push_back(random_dna(rng, 3 * c.k + rng() % (6 * c.k)));
    }
    // shared pieces between records: branches and repeats
    for (size_t r = 0; r < 3; ++r) {
        const std::string &from = c.records[rng() % records];
        std::string &to = c.records[rng() % records];
        const size_t n = c.k - 1 + rng() % 6;
        if (from.size() > n)
            to.insert(rng() % to.size(), piece_of(from, n, rng));
    }
    const DeBruijnGraph::Mode modes[]
        = { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY };
    c.mode = modes[rng() % 3];
    const Strands strands[] = { Strands::BOTH, Strands::FORWARD, Strands::REVERSE };
    c.strands = strands[rng() % 3];

    // the pattern: one to three record pieces, spliced, some bases mutated
    std::string bases;
    const size_t pieces = 1 + rng() % 3;
    for (size_t i = 0; i < pieces; ++i) {
        const std::string &source = c.records[rng() % records];
        const size_t n = std::min(source.size(), c.k + rng() % (2 * c.k));
        bases += piece_of(source, n, rng);
    }
    while (bases.size() <= c.k + 1) {
        bases += random_dna(rng, c.k);
    }
    for (size_t m = rng() % 3; m > 0; --m) {
        const size_t i = rng() % bases.size();
        bases[i] = other_base(bases[i], rng);
    }

    const uint32_t roll = rng() % 10;
    if (roll < 3) {
        // a peptide of the bases at a random frame, a few residues ambiguous or mutated
        c.kind = PatternKind::PROTEIN;
        const size_t frame = rng() % 3;
        for (size_t i = frame; i + 3 <= bases.size(); i += 3) {
            char residue = oracle_translate(1, bases.substr(i, 3));
            const uint32_t r = rng() % 12;
            if (r == 0 && residue != '*') {
                residue = 'X';
            } else if (r == 1 && (residue == 'D' || residue == 'N')) {
                residue = 'B';
            } else if (r == 1 && (residue == 'E' || residue == 'Q')) {
                residue = 'Z';
            } else if (r == 1 && (residue == 'I' || residue == 'L')) {
                residue = 'J';
            } else if (r == 2) {
                residue = "ACDEFGHIKLMNPQRSTVWY"[rng() % 20];
            }
            c.text.push_back(residue);
        }
        if (3 * c.text.size() <= c.k)
            return random_case(rng);
    } else {
        c.kind = roll < 6 ? PatternKind::DNA : PatternKind::IUPAC;
        c.text = bases;
        if (c.kind == PatternKind::IUPAC) {
            for (char &b : c.text) {
                if (!(rng() % 6))
                    b = iupac_of(b, rng);
            }
        }
    }

    const size_t L = c.kind == PatternKind::PROTEIN ? 3 * c.text.size() : c.text.size();
    const uint32_t choice = rng() % 20;
    c.plan.choice = choice < 12 ? WindowChoice::STEPS
                                : choice < 17 ? WindowChoice::ENDS : WindowChoice::FIRST;
    // W: k (0 or k), a little more, sometimes beyond L
    const uint32_t w = rng() % 6;
    c.plan.window = w == 0 ? 0
                  : w < 4 ? static_cast<uint32_t>(c.k + rng() % 6)
                  : w == 4 ? static_cast<uint32_t>(c.k)
                           : static_cast<uint32_t>(L + rng() % 3);
    const uint32_t effective = std::min<uint32_t>(c.plan.window ? c.plan.window : c.k, L);
    c.plan.step = 1 + rng() % (effective + 2);
    c.max_anchors = rng() % 4 ? 1000 : rng() % 4;
    return c;
}

TEST(PatternSegments, RandomCasesAgainstTheOracle) {
    std::mt19937 rng(20261009);
    uint64_t segments = 0;
    uint64_t with_owned = 0;
    uint64_t with_left = 0;
    uint64_t not_admitted = 0;
    for (size_t round = 0; round < 600; ++round) {
        const Case c = random_case(rng);
        const std::string what = "round " + std::to_string(round) + ": " + describe(c);
        auto graph = build_graph<DBGSuccinct>(c.k, c.records, c.mode);
        const SegmentRunResult run = check_against_oracle(*graph, c.kind, c.text, c.strands, c.plan,
                                             c.max_anchors, what);
        if (::testing::Test::HasFailure()) {
            FAIL() << what;
        }
        segments += run.segments.size();
        with_owned += run.result.work.segment_anchors_owned > 0;
        with_left += run.result.work.left_edges > 0;
        if (run.result.segments) {
            for (const auto &[o, windows] : run.result.segments->windows) {
                for (const WindowCount &w : windows) {
                    not_admitted += w.state == WindowState::NOT_ADMITTED;
                }
            }
        }
    }
    // the cases exercise what they are meant to
    EXPECT_GT(segments, 600u);
    EXPECT_GT(with_owned, 100u);
    EXPECT_GT(with_left, 100u);
    EXPECT_GT(not_admitted, 30u);
}

TEST(PatternSegments, GraphsWithoutTheMaskAnswerTheMaskedTwins) {
    std::mt19937 rng(77);
    for (size_t round = 0; round < 150; ++round) {
        Case c = random_case(rng);
        // admission compares the upper bound without the mask: room for every window
        c.max_anchors = 1000;
        const std::string what = "round " + std::to_string(round) + ": " + describe(c);
        auto masked = build_graph<DBGSuccinct>(c.k, c.records, c.mode);
        auto unmasked = build_graph<DBGSuccinct>(c.k, c.records, c.mode);
        const_cast<DBGSuccinct&>(oracle_base_dbg(*unmasked)).reset_mask();
        ASSERT_FALSE(oracle_base_dbg(*unmasked).get_mask());
        const Pattern pattern = Pattern::parse(c.kind, c.text);
        const Request request = segment_request(c.plan, c.strands, c.max_anchors);
        const SegmentRunResult a = run_segments(*masked, pattern, request);
        const SegmentRunResult b = run_segments(*unmasked, pattern, request);
        ASSERT_TRUE(a.result.segments && b.result.segments) << what;
        EXPECT_EQ(a.segments, b.segments) << what;
        EXPECT_EQ(a.recorded, b.recorded) << what;
        for (const auto &[o, windows] : a.result.segments->windows) {
            const auto &other = b.result.segments->windows.at(o);
            ASSERT_EQ(windows.size(), other.size()) << what;
            for (size_t i = 0; i < windows.size(); ++i) {
                EXPECT_EQ(windows[i].state, other[i].state) << what;
                EXPECT_EQ(windows[i].roots, other[i].roots) << what;
                EXPECT_EQ(windows[i].owned, other[i].owned) << what;
                // listed with the source dummies dropped: exact
                EXPECT_EQ(Relation::EXACT, other[i].anchors.relation) << what;
                EXPECT_EQ(windows[i].anchors.value, other[i].anchors.value) << what;
            }
        }
        EXPECT_EQ(a.result.work.segment_states_memoised, b.result.work.segment_states_memoised)
            << what;
    }
}


// ---------------------------------------------------------------- named junction cases

SegmentPlan plan_of(uint32_t window, uint32_t step, WindowChoice choice = WindowChoice::STEPS) {
    SegmentPlan plan;
    plan.window = window;
    plan.step = step;
    plan.choice = choice;
    return plan;
}

// a random sequence whose first base is not |avoid_first| and last base not |avoid_last|
std::string flank(std::mt19937 &rng, size_t n, char avoid_first, char avoid_last) {
    std::string s = random_dna(rng, n);
    if (s.front() == avoid_first)
        s.front() = other_base(avoid_first, rng);
    if (s.back() == avoid_last)
        s.back() = other_base(avoid_last, rng);
    return s;
}

std::map<OracleSegmentKey, uint64_t> segments_of(
        std::initializer_list<std::tuple<Orientation, uint32_t, uint32_t, uint64_t>> list) {
    std::map<OracleSegmentKey, uint64_t> result;
    for (const auto &[o, from, to, roots] : list) {
        result[{ o, from, to }] = roots;
    }
    return result;
}

struct TwoExons {
    std::string exon1;
    std::string exon2;
    std::string intron;
};

// two exons of |n| bases and an intron of |intron| bases whose ends cannot extend either exon
TwoExons two_exons(std::mt19937 &rng, size_t n, size_t intron) {
    TwoExons g;
    g.exon1 = random_dna(rng, n);
    g.exon2 = random_dna(rng, n);
    // the intron's first base differs from exon 2's first (exon 1 cannot run on) and its last
    // from exon 1's last (exon 2 cannot run back)
    g.intron = flank(rng, intron, g.exon2.front(), g.exon1.back());
    return g;
}

TEST(PatternSegments, JunctionTwoExons) {
    // J1, J2: the spliced gene's two halves, each a segment holding its windows; the intron's
    // length (multiple of 3 or not) does not matter at the graph level
    std::mt19937 rng(1);
    for (size_t intron : { 120, 121 }) {
        const TwoExons g = two_exons(rng, 60, intron);
        auto graph = build_graph<DBGSuccinct>(11, { g.exon1 + g.intron + g.exon2 },
                                              DeBruijnGraph::BASIC);
        const SegmentRunResult run = check_against_oracle(*graph, PatternKind::DNA, g.exon1 + g.exon2,
                                             Strands::BOTH, plan_of(0, 15), 1000, "J1");
        EXPECT_EQ(segments_of({ { Orientation::FORWARD, 0, 60, 1 },
                                { Orientation::FORWARD, 60, 120, 1 } }), run.segments);
        // windows 0, 15, ..., 105, 109 of 11 bases: those across the junction (50 .. 59 would
        // be) have none; 0 and 60 are the roots, 15 .. 45 and 75 .. 109 owned
        const auto &windows = run.result.segments->windows.at(Orientation::FORWARD);
        ASSERT_EQ(9u, windows.size());
        EXPECT_EQ(2u, run.result.work.segment_anchors_extended);
        EXPECT_EQ(1u, windows[0].roots);
        EXPECT_EQ(1u, windows[4].roots);
        EXPECT_EQ(60u, windows[4].offset);
        // the reverse strand holds nothing
        for (const WindowCount &w : run.result.segments->windows.at(Orientation::REVERSE)) {
            EXPECT_EQ(0u, w.anchors.value);
        }
    }
}

TEST(PatternSegments, JunctionReverseStrand) {
    // J7: the gene on the record's reverse strand: the same intervals of P on "-"
    std::mt19937 rng(2);
    const TwoExons g = two_exons(rng, 60, 120);
    auto graph = build_graph<DBGSuccinct>(
            11, { oracle_rev_comp(g.exon1 + g.intron + g.exon2) }, DeBruijnGraph::BASIC);
    const SegmentRunResult run = check_against_oracle(*graph, PatternKind::DNA, g.exon1 + g.exon2,
                                         Strands::BOTH, plan_of(0, 15), 1000, "J7");
    EXPECT_EQ(segments_of({ { Orientation::REVERSE, 0, 60, 1 },
                            { Orientation::REVERSE, 60, 120, 1 } }), run.segments);
}

TEST(PatternSegments, JunctionShortExon) {
    // J5: an exon of 9 bases (< k = 11) between two long ones: uncovered; the two others found
    std::mt19937 rng(3);
    const std::string e1 = random_dna(rng, 40);
    const std::string e2 = random_dna(rng, 9);
    const std::string e3 = random_dna(rng, 40);
    const std::string i1 = flank(rng, 100, e2.front(), e1.back());
    const std::string i2 = flank(rng, 130, e3.front(), e2.back());
    auto graph = build_graph<DBGSuccinct>(11, { e1 + i1 + e2 + i2 + e3 },
                                          DeBruijnGraph::BASIC);
    const SegmentRunResult run = check_against_oracle(*graph, PatternKind::DNA, e1 + e2 + e3,
                                         Strands::FORWARD, plan_of(0, 7), 1000, "J5");
    EXPECT_EQ(segments_of({ { Orientation::FORWARD, 0, 40, 1 },
                            { Orientation::FORWARD, 49, 89, 1 } }), run.segments);
}

TEST(PatternSegments, JunctionRecordEnds) {
    // J11: the record ends inside the pattern; J13: the two halves in two records
    std::mt19937 rng(4);
    const std::string gene = random_dna(rng, 100);
    auto ends = build_graph<DBGSuccinct>(11, { gene.substr(0, 70) }, DeBruijnGraph::BASIC);
    SegmentRunResult run = check_against_oracle(*ends, PatternKind::DNA, gene, Strands::FORWARD,
                                   plan_of(0, 15), 1000, "J11");
    EXPECT_EQ(segments_of({ { Orientation::FORWARD, 0, 70, 1 } }), run.segments);

    const std::string a = random_dna(rng, 60);
    const std::string b = random_dna(rng, 60);
    auto two = build_graph<DBGSuccinct>(11, { random_dna(rng, 20) + a, b + random_dna(rng, 20) },
                                        DeBruijnGraph::BASIC);
    run = check_against_oracle(*two, PatternKind::DNA, a + b, Strands::FORWARD,
                               plan_of(0, 15), 1000, "J13");
    // (the record pieces around them may extend either half by chance: the oracle decides;
    // both halves are found whole)
    bool first = false;
    bool second = false;
    for (const auto &[key, roots] : run.segments) {
        first |= std::get<1>(key) == 0 && std::get<2>(key) >= 60;
        second |= std::get<1>(key) <= 60 && std::get<2>(key) == 120;
    }
    EXPECT_TRUE(first && second) << describe_segments(run.segments);
}

TEST(PatternSegments, PalindromeOnce) {
    // J15: a palindromic pattern is one orientation
    std::mt19937 rng(5);
    const std::string half = random_dna(rng, 15);
    const std::string p = half + oracle_rev_comp(half);
    auto graph = build_graph<DBGSuccinct>(11, { random_dna(rng, 10) + p + random_dna(rng, 10) },
                                          DeBruijnGraph::BASIC);
    const SegmentRunResult run = check_against_oracle(*graph, PatternKind::DNA, p, Strands::BOTH,
                                         plan_of(11, 5), 1000, "J15");
    ASSERT_EQ(1u, run.result.segments->windows.size());
    EXPECT_TRUE(run.result.segments->windows.count(Orientation::PALINDROMIC));
    EXPECT_EQ(1u, run.segments.count({ Orientation::PALINDROMIC, 0, 30 }));
}

TEST(PatternSegments, MosaicIsOneExtentAtTheGraphLevel) {
    // J17: the pattern's k-mers spread over two records that share 10 bases: no record holds
    // more than its half, the graph's walks switch records in the shared part
    std::mt19937 rng(6);
    const std::string a = random_dna(rng, 30);
    const std::string shared = random_dna(rng, 10);
    const std::string d = random_dna(rng, 30);
    auto graph = build_graph<DBGSuccinct>(
            11, { a + shared + random_dna(rng, 30), random_dna(rng, 30) + shared + d },
            DeBruijnGraph::BASIC);
    const SegmentRunResult run = check_against_oracle(*graph, PatternKind::DNA, a + shared + d,
                                         Strands::FORWARD, plan_of(0, 15), 1000, "J17");
    EXPECT_EQ(1u, run.segments.count({ Orientation::FORWARD, 0, 70 }))
        << describe_segments(run.segments);
}

TEST(PatternSegments, PeptideAcrossAJunction) {
    // J3: a peptide spliced from two exons, the junction at every phase: two segments, their
    // bases as the oracle's (a codon split by the junction may be held by either side by
    // chance)
    std::mt19937 rng(8);
    for (size_t phase = 0; phase < 3; ++phase) {
        const std::string exon1 = random_dna(rng, 45 + phase);
        const std::string exon2 = random_dna(rng, 45 - phase);
        const std::string intron = flank(rng, 99, exon2.front(), exon1.back());
        std::string spliced = exon1 + exon2;
        std::string peptide;
        for (size_t i = 0; i + 3 <= spliced.size(); i += 3) {
            const char residue = oracle_translate(1, spliced.substr(i, 3));
            peptide.push_back(residue == '*' ? 'W' : residue);
        }
        auto graph = build_graph<DBGSuccinct>(11, { exon1 + intron + exon2 },
                                              DeBruijnGraph::CANONICAL);
        for (uint32_t step : { 1u, 5u, 15u }) {
            const SegmentRunResult run = check_against_oracle(*graph, PatternKind::PROTEIN, peptide,
                                                 Strands::BOTH, plan_of(0, step), 1000,
                                                 "J3 phase " + std::to_string(phase));
            EXPECT_FALSE(run.segments.empty());
        }
    }
}


// ---------------------------------------------------------------- windows and steps

TEST(PatternSegments, WindowOffsetsAgainstTheOracle) {
    // every (L, W, S) with k <= W <= L <= 4k and 1 <= S <= W + 2 on k = 7
    const uint32_t k = 7;
    for (uint32_t L = k + 1; L <= 4 * k; ++L) {
        for (uint32_t w = k; w <= L; ++w) {
            for (uint32_t s = 1; s <= w + 2; ++s) {
                EXPECT_EQ(oracle_window_grid(L, w, s, WindowChoice::STEPS),
                          segment_window_offsets(L, w, s, WindowChoice::STEPS))
                    << L << " " << w << " " << s;
            }
            EXPECT_EQ(oracle_window_grid(L, w, 1, WindowChoice::ENDS),
                      segment_window_offsets(L, w, 1, WindowChoice::ENDS));
            EXPECT_EQ(std::vector<uint32_t> { 0 },
                      segment_window_offsets(L, w, 1, WindowChoice::FIRST));
        }
    }
    // one window when 0 = L - W
    EXPECT_EQ(std::vector<uint32_t> { 0 }, segment_window_offsets(20, 20, 3, WindowChoice::ENDS));
    EXPECT_THROW(segment_window_offsets(20, 21, 3, WindowChoice::STEPS), std::invalid_argument);
    EXPECT_THROW(segment_window_offsets(20, 10, 0, WindowChoice::STEPS), std::invalid_argument);
}

// a gene, a record holding it, and the plain DNA search of it at k = 11
struct Gene {
    std::string gene;
    std::shared_ptr<DeBruijnGraph> graph;
};

Gene gene_case(std::mt19937 &rng, size_t n) {
    Gene g;
    g.gene = random_dna(rng, n);
    g.graph = build_graph<DBGSuccinct>(11, { random_dna(rng, 15) + g.gene + random_dna(rng, 15) },
                                       DeBruijnGraph::BASIC);
    return g;
}

TEST(PatternSegments, StepBeyondTheLastWindowIsEnds) {
    // S > L - W: the first and the last window only, the answer "ends"'s
    std::mt19937 rng(9);
    for (size_t round = 0; round < 10; ++round) {
        const Gene g = gene_case(rng, 40 + rng() % 30);
        const std::string text = g.gene.substr(0, 20) + random_dna(rng, 3)
                + g.gene.substr(23);
        const uint32_t L = static_cast<uint32_t>(text.size());
        const SegmentRunResult steps = check_against_oracle(*g.graph, PatternKind::DNA, text, Strands::BOTH,
                                               plan_of(13, L), 1000, "steps");
        const SegmentRunResult ends = check_against_oracle(*g.graph, PatternKind::DNA, text, Strands::BOTH,
                                              plan_of(13, 1, WindowChoice::ENDS), 1000, "ends");
        EXPECT_EQ(ends.segments, steps.segments);
        EXPECT_EQ(ends.recorded, steps.recorded);
        EXPECT_EQ(2u, steps.result.segments->windows.at(Orientation::FORWARD).size());
    }
}

TEST(PatternSegments, CompleteForLengthNearTheStepMaximum) {
    // W' + S - 1 is computed without wrapping: a step near the type's maximum (the windows are
    // "ends"'s) states a length at least the window, saturated at the maximum; the step that
    // reaches the maximum exactly states it
    std::mt19937 rng(13);
    const Gene g = gene_case(rng, 40);
    const Pattern pattern = Pattern::parse(PatternKind::DNA, g.gene);
    const uint32_t L = 40;
    const uint32_t k = 11;
    const uint32_t max = std::numeric_limits<uint32_t>::max();
    const std::vector<std::pair<uint32_t, uint32_t>> cases {
        { 15, 25 }, { L - k, L - 1 }, { max - 10, max }, { max - 9, max },
        { 4294967290u, max }, { max, max } };
    for (const auto &[step, stated] : cases) {
        SCOPED_TRACE("step " + std::to_string(step));
        const SegmentRunResult run = run_segments(
                *g.graph, pattern, segment_request(plan_of(0, step), Strands::BOTH, 1000));
        ASSERT_FALSE(run.result.refusal);
        ASSERT_TRUE(run.result.segments);
        const SegmentCounts &c = *run.result.segments;
        EXPECT_EQ(k, c.window);
        EXPECT_EQ(step, c.step);
        ASSERT_TRUE(c.complete_for_length);
        EXPECT_EQ(stated, *c.complete_for_length);
        EXPECT_LE(c.window, *c.complete_for_length);
        std::vector<uint32_t> offsets;
        for (const WindowCount &w : c.windows.at(Orientation::FORWARD)) {
            offsets.push_back(w.offset);
        }
        EXPECT_EQ(segment_window_offsets(L, k, step, WindowChoice::STEPS), offsets);
        if (step >= L - k) {
            EXPECT_EQ((std::vector<uint32_t> { 0, L - k }), offsets);
        }
    }
}

TEST(PatternSegments, WindowOfTheWholePatternIsFullLength) {
    // W = L: a segment is reported only when it is the whole pattern; W > L: clipped to L
    std::mt19937 rng(10);
    const Gene g = gene_case(rng, 60);
    const std::string whole = g.gene.substr(5, 40);
    for (uint32_t w : { 40u, 41u, 100u }) {
        const SegmentRunResult run = check_against_oracle(*g.graph, PatternKind::DNA, whole, Strands::BOTH,
                                             plan_of(w, 3), 1000, "W " + std::to_string(w));
        EXPECT_EQ(40u, run.result.segments->window);
        EXPECT_EQ(segments_of({ { Orientation::FORWARD, 0, 40, 1 } }), run.segments);
    }
    // a mismatch inside: no full-length instance, nothing reported, segments below the window
    std::string broken = whole;
    broken[20] = other_base(broken[20], rng);
    const SegmentRunResult run = check_against_oracle(*g.graph, PatternKind::DNA, broken, Strands::BOTH,
                                         plan_of(40, 3), 1000, "W = L broken");
    EXPECT_TRUE(run.segments.empty());
    EXPECT_GT(run.result.work.segment_occurrences_below_window, 0u);
}

TEST(PatternSegments, GuaranteeAndItsTightness) {
    // a match of W + S - 1 bases holds a window wherever it lies; one of W + S - 2 does not at
    // the grid's worst alignment
    std::mt19937 rng(12);
    const uint32_t W = 13;
    const uint32_t S = 6;
    const std::string p = random_dna(rng, 70);
    for (uint32_t M : { W + S - 1, W + S - 2 }) {
        size_t missed = 0;
        for (uint32_t a = 0; a + M <= p.size(); ++a) {
            // the record: P[a, a + M) between flanks that cannot extend it
            const std::string left = flank(rng, 20, 'A', a ? p[a - 1] : 'A');
            const std::string right = flank(rng, 20, a + M < p.size() ? p[a + M] : 'A', 'A');
            auto graph = build_graph<DBGSuccinct>(11, { left + p.substr(a, M) + right },
                                                  DeBruijnGraph::BASIC);
            const SegmentRunResult run = check_against_oracle(*graph, PatternKind::DNA, p, Strands::FORWARD,
                                                 plan_of(W, S), 1000,
                                                 "M " + std::to_string(M) + " at "
                                                    + std::to_string(a));
            missed += !run.segments.count({ Orientation::FORWARD, a, a + M });
        }
        if (M == W + S - 1) {
            EXPECT_EQ(0u, missed);
        } else {
            EXPECT_GT(missed, 0u);
        }
    }
}

TEST(PatternSegments, WindowsNotAdmitted) {
    // the middle windows ambiguous (an R that two records instantiate) and max_anchors 1: not
    // admitted; the segments through them found from their neighbours, the left arm passing
    // them; the first window ambiguous: counts.anchors AT_LEAST, the prefix found from the next.
    // Without the mask the cap compares the candidates counted (source dummies unchecked among
    // them): a window it stopped states them apart, AT_LEAST above max_anchors, its anchors only
    // what its discovery proved
    std::mt19937 rng(13);
    const std::string gene = random_dna(rng, 90);
    for (size_t at : { 3, 40 }) {
        std::string variant = gene;
        variant[at] = variant[at] == 'A' ? 'G' : variant[at] == 'G' ? 'A'
                    : variant[at] == 'C' ? 'T' : 'C';
        std::string text = gene;
        text[at] = variant[at] == 'A' || variant[at] == 'G' ? 'R' : 'Y';
        const std::vector<std::string> records { gene,
                                                 variant.substr(at > 10 ? at - 10 : 0, 21) };
        auto graph = build_graph<DBGSuccinct>(11, records, DeBruijnGraph::BASIC);
        const SegmentRunResult run = check_against_oracle(*graph, PatternKind::IUPAC, text, Strands::FORWARD,
                                             plan_of(11, 10), 1, "R at " + std::to_string(at));
        const auto &windows = run.result.segments->windows.at(Orientation::FORWARD);
        size_t not_admitted = 0;
        for (const WindowCount &w : windows) {
            EXPECT_FALSE(w.candidates);
            if (w.state == WindowState::NOT_ADMITTED) {
                ++not_admitted;
                EXPECT_LE(w.offset, at);
                EXPECT_GT(w.offset + 11, at);
            }
        }
        EXPECT_GT(not_admitted, 0u);

        auto unmasked = build_graph<DBGSuccinct>(11, records, DeBruijnGraph::BASIC);
        const_cast<DBGSuccinct&>(oracle_base_dbg(*unmasked)).reset_mask();
        ASSERT_FALSE(oracle_base_dbg(*unmasked).get_mask());
        const SegmentRunResult twin = run_segments(*unmasked, Pattern::parse(PatternKind::IUPAC,
                                                                             text),
                                                   segment_request(plan_of(11, 10),
                                                                   Strands::FORWARD, 1));
        ASSERT_TRUE(twin.result.segments);
        const auto &twin_windows = twin.result.segments->windows.at(Orientation::FORWARD);
        ASSERT_EQ(windows.size(), twin_windows.size());
        size_t capped = 0;
        for (size_t i = 0; i < windows.size(); ++i) {
            const WindowCount &w = twin_windows[i];
            // the masked twin's count is the truth: EXACT, or 2 for a window over the R
            const uint64_t truth = windows[i].anchors.relation == Relation::EXACT
                    ? windows[i].anchors.value : 2;
            if (w.state != WindowState::NOT_ADMITTED) {
                EXPECT_FALSE(w.candidates) << w.offset;
                continue;
            }
            ++capped;
            EXPECT_EQ(Relation::AT_LEAST, w.anchors.relation) << w.offset;
            EXPECT_GE(truth, w.anchors.value) << w.offset;
            ASSERT_TRUE(w.candidates) << w.offset;
            EXPECT_EQ(Unit::ANCHORS, w.candidates->unit);
            EXPECT_EQ(Relation::AT_LEAST, w.candidates->relation) << w.offset;
            EXPECT_LT(1u, w.candidates->value) << w.offset;
        }
        // at least the windows over the R
        EXPECT_LE(not_admitted, capped);
        // the whole gene is still one segment, from the first admitted window on
        EXPECT_EQ(1u, run.segments.count({ Orientation::FORWARD, 0, 90 }))
            << describe_segments(run.segments);
        if (at < 11) {
            EXPECT_EQ(Relation::AT_LEAST, run.result.anchors->total.relation);
            EXPECT_EQ(2u, run.result.anchors->total.value);
        }
    }
}

TEST(PatternSegments, NoWindowAdmitted) {
    std::mt19937 rng(14);
    const Gene g = gene_case(rng, 60);
    // (a window with 0 anchors is extended: the reverse strand's would make the search
    // completed)
    const SegmentRunResult run = check_against_oracle(*g.graph, PatternKind::DNA, g.gene, Strands::FORWARD,
                                         plan_of(0, 15), 0, "max_anchors 0");
    EXPECT_EQ(Extension::NOT_ADMITTED, run.result.segments->search);
    EXPECT_TRUE(run.segments.empty());
    EXPECT_EQ(0u, run.result.work.segment_anchors_extended);
}

TEST(PatternSegments, WindowWithoutAnchorsIsExtended) {
    // a window with 0 anchors was searched completely: EXTENDED (EXACT 0, no root), never
    // NOT_ADMITTED, so that such a window among windows not admitted makes the search
    // COMPLETED, not "none_admitted". The graph holds the gene's forward strand only and
    // max_anchors is 0: every forward window is not admitted, every reverse one has 0 anchors
    std::mt19937 rng(14);
    const Gene g = gene_case(rng, 60);
    const Pattern pattern = Pattern::parse(PatternKind::DNA, g.gene);
    const SegmentRunResult both = run_segments(*g.graph, pattern,
                                               segment_request(plan_of(0, 15), Strands::BOTH, 0));
    ASSERT_TRUE(both.result.segments);
    EXPECT_FALSE(both.result.stop);
    const SegmentCounts &c = *both.result.segments;
    for (const WindowCount &w : c.windows.at(Orientation::FORWARD)) {
        EXPECT_EQ(WindowState::NOT_ADMITTED, w.state) << w.offset;
    }
    for (const WindowCount &w : c.windows.at(Orientation::REVERSE)) {
        EXPECT_EQ(WindowState::EXTENDED, w.state) << w.offset;
        EXPECT_EQ(Relation::EXACT, w.anchors.relation) << w.offset;
        EXPECT_EQ(0u, w.anchors.value) << w.offset;
        EXPECT_EQ(0u, w.roots) << w.offset;
    }
    EXPECT_EQ(Extension::COMPLETED, c.search);
    EXPECT_TRUE(both.segments.empty());
    EXPECT_EQ(0u, both.result.work.segment_anchors_extended);
    // the forward strand alone: none admitted
    const SegmentRunResult forward_only = run_segments(
            *g.graph, pattern, segment_request(plan_of(0, 15), Strands::FORWARD, 0));
    ASSERT_TRUE(forward_only.result.segments);
    EXPECT_EQ(Extension::NOT_ADMITTED, forward_only.result.segments->search);
}

TEST(PatternSegments, PeptideWithAStopItsTableCannotCode) {
    // table 27 has no codon for '*' (its stop codons are sense codons there): a peptide holding
    // one has no instance, yet it is searched as any other. The windows whose anchored bases
    // hold a base of the '*' codon have EXACT 0 anchors (extended); the others are searched,
    // and the parts before and after the '*' are segments; the note no_stop_codon
    std::mt19937 rng(19);
    const int table = 27;
    std::string gene;
    std::string peptide;
    while (peptide.size() < 24) {
        const std::string codon = random_dna(rng, 3);
        const char residue = oracle_translate(table, codon);
        if (residue == '*')
            continue;
        gene += codon;
        peptide.push_back(residue);
    }
    const uint32_t star = 10;
    peptide[star] = '*';
    const uint32_t from = 3 * star;
    const uint32_t to = from + 3;
    const uint32_t L = 3 * static_cast<uint32_t>(peptide.size());
    auto graph = build_graph<DBGSuccinct>(11, { random_dna(rng, 15) + gene
                                                    + random_dna(rng, 15) },
                                          DeBruijnGraph::BASIC);
    const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, peptide,
                                           GeneticCode::get(table));
    ASSERT_FALSE(pattern.has_instances());
    const SegmentRunResult run = run_segments(*graph, pattern,
                                              segment_request(plan_of(0, 1), Strands::FORWARD,
                                                              1000));
    ASSERT_FALSE(run.result.refusal);
    ASSERT_TRUE(run.result.segments);
    EXPECT_FALSE(run.result.stop);
    EXPECT_EQ(Extension::COMPLETED, run.result.segments->search);
    const auto &windows = run.result.segments->windows.at(Orientation::FORWARD);
    EXPECT_EQ(L - 11 + 1, windows.size());
    for (const WindowCount &w : windows) {
        EXPECT_EQ(WindowState::EXTENDED, w.state) << w.offset;
        EXPECT_EQ(Relation::EXACT, w.anchors.relation) << w.offset;
        if (w.offset < to && w.offset + 11 > from) {
            EXPECT_EQ(0u, w.anchors.value) << w.offset;
        } else {
            EXPECT_LE(1u, w.anchors.value) << w.offset;
        }
    }
    EXPECT_EQ(1u, run.segments.count({ Orientation::FORWARD, 0, from }))
        << describe_segments(run.segments);
    EXPECT_EQ(1u, run.segments.count({ Orientation::FORWARD, to, L }))
        << describe_segments(run.segments);
    EXPECT_EQ(2u, run.segments.size()) << describe_segments(run.segments);
    EXPECT_NE(run.result.notes.end(), std::find(run.result.notes.begin(),
                                                run.result.notes.end(),
                                                std::string(kNoteNoStopCodon)));
}

TEST(PatternSegments, WindowsBelowTheFloor) {
    // an N run over a window's anchored bases: BELOW_FLOOR, not discovered; every window below:
    // the refusal, nothing charged
    std::mt19937 rng(15);
    const Gene g = gene_case(rng, 70);
    std::string text = g.gene;
    for (size_t i = 30; i < 45; ++i) {
        text[i] = 'N';
    }
    const Pattern pattern = Pattern::parse(PatternKind::IUPAC, text);
    Request request = segment_request(plan_of(11, 10), Strands::FORWARD, 1000);
    request.min_information_bits = 14;
    SegmentRunResult run = run_segments(*g.graph, pattern, request);
    ASSERT_FALSE(run.result.refusal);
    const auto &windows = run.result.segments->windows.at(Orientation::FORWARD);
    size_t below = 0;
    for (const WindowCount &w : windows) {
        if (w.state == WindowState::BELOW_FLOOR) {
            ++below;
            EXPECT_EQ(Relation::UNKNOWN, w.anchors.relation);
            EXPECT_LT(w.information_bits, 14);
        } else {
            EXPECT_EQ(WindowState::EXTENDED, w.state);
        }
    }
    EXPECT_EQ(2u, below);  // the windows at 30 and 40
    EXPECT_EQ(windows.size() - below, run.result.work.segment_windows);
    // the N run is admitted by the arms: the gene is one segment
    EXPECT_EQ(1u, run.segments.count({ Orientation::FORWARD, 0, 70 }));

    const Pattern low = Pattern::parse(PatternKind::IUPAC, std::string(40, 'N'));
    request.min_information_bits = 1;
    run = run_segments(*g.graph, low, request);
    ASSERT_TRUE(run.result.refusal);
    EXPECT_EQ("information_below_floor", run.result.refusal->code);
    EXPECT_EQ(0u, run.steps_used);
}

TEST(PatternSegments, TooManyWindows) {
    std::mt19937 rng(16);
    const Gene g = gene_case(rng, 60);
    SegmentPlan plan = plan_of(11, 1);
    plan.max_windows = 10;
    SegmentRunResult run = run_segments(*g.graph, Pattern::parse(PatternKind::DNA, g.gene),
                           segment_request(plan, Strands::BOTH, 1000));
    ASSERT_TRUE(run.result.refusal);
    EXPECT_EQ("too_many_windows", run.result.refusal->code);
    // 50 bases of starts over 9 gaps: a step of 6 fits (0, 6, ..., 48, 49 are 10 windows)
    EXPECT_NE(std::string::npos, run.result.refusal->message.find("a step of 6 bases"))
        << run.result.refusal->message;
    EXPECT_EQ(0u, run.steps_used);
    plan.step = 6;
    run = run_segments(*g.graph, Pattern::parse(PatternKind::DNA, g.gene),
                       segment_request(plan, Strands::BOTH, 1000));
    EXPECT_FALSE(run.result.refusal);
    plan.step = 5;
    run = run_segments(*g.graph, Pattern::parse(PatternKind::DNA, g.gene),
                       segment_request(plan, Strands::BOTH, 1000));
    EXPECT_TRUE(run.result.refusal);
}

TEST(PatternSegments, RequestsTheRouteMustRefuseFirst) {
    std::mt19937 rng(17);
    const Gene g = gene_case(rng, 40);
    const Pattern pattern = Pattern::parse(PatternKind::DNA, g.gene);
    PatternSearch search(*g.graph);
    GraphSegmentStore store(40, kManyBytes);
    Budget budget(kManySteps, unbounded_deadline());
    Request request = segment_request(plan_of(0, 15), Strands::BOTH, 1000);
    // no sink
    EXPECT_THROW(search.count(pattern, request, budget), std::invalid_argument);
    request.segment_sink = &store;
    // a window below k, a step of 0
    request.segments.window = 10;
    EXPECT_THROW(search.count(pattern, request, budget), std::invalid_argument);
    request.segments.window = 0;
    request.segments.step = 0;
    EXPECT_THROW(search.count(pattern, request, budget), std::invalid_argument);
    request.segments.step = 15;
    // mismatches: a segment holds no substitution (the route refuses the combination first)
    request.max_mismatches = 1;
    EXPECT_THROW(search.count(pattern, request, budget), std::invalid_argument);
    EXPECT_EQ(0u, budget.steps_used());
    // a pattern of at most k bases is answered as without segments, its mismatches served
    const Result short_pattern = search.count(Pattern::parse(PatternKind::DNA,
                                                             g.gene.substr(0, 9)),
                                              request, budget);
    EXPECT_FALSE(short_pattern.segments);
    ASSERT_TRUE(short_pattern.contexts);
    EXPECT_EQ(2u, short_pattern.contexts->by_mismatches.size());
    request.max_mismatches = 0;
    // enumerate() does not serve them
    request.mode = Mode::PARTIAL;
    EXPECT_THROW(search.enumerate(pattern, request, budget, [](const Context&) {}),
                 std::invalid_argument);
    EXPECT_NO_THROW(search.count(pattern, request, budget));
    EXPECT_EQ(0u, store.held() - store.size() * GraphSegmentStore::kGraphSegmentBytes);
}

TEST(PatternSegments, ShortPatternsAreAnsweredAsWithout) {
    // L <= k: segments do not apply; the answer is the plain one, step for step
    std::mt19937 rng(18);
    const Gene g = gene_case(rng, 40);
    const Pattern pattern = Pattern::parse(PatternKind::DNA, g.gene.substr(3, 9));
    PatternSearch search(*g.graph);
    Request plain;
    plain.min_information_bits = 0;
    Request with = segment_request(plan_of(0, 15), Strands::BOTH, 1000);
    with.scope = plain.scope;
    Budget a(kManySteps, unbounded_deadline());
    Budget b(kManySteps, unbounded_deadline());
    const Result x = search.count(pattern, plain, a);
    const Result y = search.count(pattern, with, b);
    ASSERT_TRUE(x.contexts && y.contexts);
    EXPECT_FALSE(y.segments);
    EXPECT_EQ(x.contexts->total.value, y.contexts->total.value);
    EXPECT_EQ(x.contexts->total.relation, y.contexts->total.relation);
    EXPECT_EQ(a.steps_used(), b.steps_used());
}


// ---------------------------------------------------------------- stops

/**
 * The invariants of a run stopped somewhere against the full run: the extents handed over a
 * prefix of the full run's; the windows in window-major order the same up to the one the stop
 * touched (STOPPED), NOT_REACHED after it (but those below the floor); the memory released.
 */
void check_stopped(const SegmentRunResult &full, const SegmentRunResult &stopped, const std::string &what) {
    const Result &r = stopped.result;
    ASSERT_TRUE(r.segments) << what;
    ASSERT_LE(stopped.recorded.size(), full.recorded.size()) << what;
    for (size_t i = 0; i < stopped.recorded.size(); ++i) {
        ASSERT_EQ(full.recorded[i], stopped.recorded[i]) << what << " extent " << i;
    }
    EXPECT_EQ(stopped.store_size * GraphSegmentStore::kGraphSegmentBytes, stopped.held) << what;
    if (!r.stop) {
        EXPECT_EQ(full.recorded, stopped.recorded) << what;
        return;
    }
    if (r.segments->search == Extension::NOT_STARTED) {
        // a stop at the pattern's first reading: no window reached
        EXPECT_TRUE(stopped.recorded.empty()) << what;
        EXPECT_EQ(0u, r.work.segment_windows) << what;
        for (const auto &[o, windows] : r.segments->windows) {
            for (const WindowCount &w : windows) {
                EXPECT_TRUE(w.state == WindowState::NOT_REACHED
                                || w.state == WindowState::BELOW_FLOOR) << what;
            }
        }
        return;
    }
    EXPECT_EQ(Extension::STOPPED, r.segments->search) << what;
    // window-major order
    const auto &a = full.result.segments->windows;
    const auto &b = r.segments->windows;
    size_t n = 0;
    for (const auto &[o, windows] : a) {
        n = std::max(n, windows.size());
    }
    bool touched = false;
    for (size_t i = 0; i < n; ++i) {
        for (const auto &[o, windows] : a) {
            const WindowCount &x = windows[i];
            const WindowCount &y = b.at(o)[i];
            if (x.state == WindowState::BELOW_FLOOR) {
                EXPECT_EQ(WindowState::BELOW_FLOOR, y.state) << what;
                continue;
            }
            if (touched) {
                EXPECT_EQ(WindowState::NOT_REACHED, y.state) << what;
                EXPECT_EQ(Relation::UNKNOWN, y.anchors.relation) << what;
                continue;
            }
            if (y.state == WindowState::STOPPED) {
                touched = true;
                EXPECT_LE(y.roots, x.roots) << what;
                continue;
            }
            EXPECT_EQ(x.state, y.state) << what << " window " << i;
            EXPECT_EQ(x.roots, y.roots) << what;
            EXPECT_EQ(x.owned, y.owned) << what;
            EXPECT_EQ(x.anchors.value, y.anchors.value) << what;
        }
    }
    EXPECT_TRUE(touched) << what;
}

// the stop cases: a gene with variants on both strands (DNA), its peptide, an IUPAC mosaic
struct StopCase {
    std::shared_ptr<DeBruijnGraph> graph;
    Pattern pattern;
    Request request;
};

std::vector<StopCase> stop_cases() {
    std::mt19937 rng(19);
    std::vector<StopCase> cases;
    const TwoExons g = two_exons(rng, 60, 120);
    std::string variant = g.exon1 + g.exon2;
    variant[30] = other_base(variant[30], rng);
    auto dna = build_graph<DBGSuccinct>(
            11, { g.exon1 + g.intron + g.exon2, oracle_rev_comp(variant) },
            DeBruijnGraph::BASIC);
    cases.push_back({ dna, Pattern::parse(PatternKind::DNA, g.exon1 + g.exon2),
                      segment_request(plan_of(0, 7), Strands::BOTH, 1000) });
    std::string peptide;
    const std::string spliced = g.exon1 + g.exon2;
    for (size_t i = 0; i + 3 <= spliced.size(); i += 3) {
        const char residue = oracle_translate(1, spliced.substr(i, 3));
        peptide.push_back(residue == '*' ? 'X' : residue);
    }
    cases.push_back({ dna, Pattern::parse(PatternKind::PROTEIN, peptide),
                      segment_request(plan_of(0, 5), Strands::BOTH, 1000) });
    // an IUPAC pattern over records that branch: many walks, few states
    std::vector<std::string> records;
    const std::string base = random_dna(rng, 80);
    for (size_t i = 0; i < 4; ++i) {
        std::string r = base;
        for (size_t j = 10 + i; j < 80; j += 13) {
            r[j] = other_base(r[j], rng);
        }
        records.push_back(r);
    }
    std::string text = base;
    for (size_t j = 10; j < 80; j += 13) {
        for (size_t i = 0; i < 4 && j + i < 80; ++i) {
            text[j + i] = 'N';
        }
    }
    auto mosaic = build_graph<DBGSuccinct>(7, records, DeBruijnGraph::PRIMARY);
    cases.push_back({ mosaic, Pattern::parse(PatternKind::IUPAC, text),
                      segment_request(plan_of(9, 4, WindowChoice::STEPS), Strands::BOTH,
                                      1000) });
    return cases;
}

TEST(PatternSegments, StopSweepMaxSteps) {
    size_t case_id = 0;
    for (const StopCase &c : stop_cases()) {
        const SegmentRunResult full = run_segments(*c.graph, c.pattern, c.request);
        ASSERT_FALSE(full.result.stop);
        ASSERT_FALSE(full.recorded.empty());
        const uint64_t total = full.steps_used;
        const uint64_t stride = std::max<uint64_t>(1, total / 1500);
        for (uint64_t n = 0; n <= total; n += stride) {
            const std::string what = "case " + std::to_string(case_id) + " max_steps "
                    + std::to_string(n);
            Limits limits;
            limits.max_steps = n;
            const SegmentRunResult run = run_segments(*c.graph, c.pattern, c.request, limits);
            check_stopped(full, run, what);
            EXPECT_LE(run.steps_used, n) << what;
            if (n < total) {
                ASSERT_TRUE(run.result.stop) << what;
                EXPECT_EQ(StopReason::MAX_STEPS, run.result.stop->reason) << what;
            }
            // deterministic
            const SegmentRunResult again = run_segments(*c.graph, c.pattern, c.request, limits);
            EXPECT_EQ(run.recorded, again.recorded) << what;
            EXPECT_EQ(run.result.work.steps, again.result.work.steps) << what;
        }
        ++case_id;
    }
}

// a clock that reads the start until its |expire|-th reading, then the start plus an hour
struct CountingClock {
    Deadline::Clock::time_point start = Deadline::Clock::now();
    std::shared_ptr<uint64_t> readings = std::make_shared<uint64_t>(0);
    uint64_t expire = std::numeric_limits<uint64_t>::max();

    Deadline deadline() const {
        auto readings = this->readings;
        auto start = this->start;
        const uint64_t expire = this->expire;
        return Deadline(start, 1000, 0, [=]() {
            return ++*readings >= expire ? start + std::chrono::hours(1) : start;
        });
    }
};

TEST(PatternSegments, StopSweepClock) {
    size_t case_id = 0;
    for (const StopCase &c : stop_cases()) {
        CountingClock clock;
        Limits limits;
        limits.deadline = [&]() { return clock.deadline(); };
        const SegmentRunResult full = run_segments(*c.graph, c.pattern, c.request, limits);
        ASSERT_FALSE(full.result.stop);
        const uint64_t readings = *clock.readings;
        ASSERT_GT(readings, 3u);
        for (uint64_t r = 1; r <= readings; ++r) {
            const std::string what = "case " + std::to_string(case_id) + " reading "
                    + std::to_string(r);
            CountingClock expiring;
            expiring.expire = r;
            limits.deadline = [&]() { return expiring.deadline(); };
            const SegmentRunResult run = run_segments(*c.graph, c.pattern, c.request, limits);
            check_stopped(full, run, what);
            // a reading after the search (the low-complexity diagnostic) leaves it complete
            if (run.result.stop) {
                EXPECT_EQ(StopReason::TIME, run.result.stop->reason) << what;
                EXPECT_TRUE(run.result.time_limited) << what;
            }
        }
        ++case_id;
    }
}

TEST(PatternSegments, StopSweepMemoryAndMaxSegments) {
    size_t case_id = 0;
    for (const StopCase &c : stop_cases()) {
        const SegmentRunResult full = run_segments(*c.graph, c.pattern, c.request);
        const uint64_t stride = std::max<uint64_t>(16, full.peak / 600);
        for (uint64_t bytes = 0; bytes <= full.peak; bytes += stride) {
            const std::string what = "case " + std::to_string(case_id) + " bytes "
                    + std::to_string(bytes);
            Limits limits;
            limits.max_bytes = bytes;
            const SegmentRunResult run = run_segments(*c.graph, c.pattern, c.request, limits);
            check_stopped(full, run, what);
            EXPECT_LE(run.peak, bytes) << what;
            if (run.result.stop) {
                EXPECT_EQ(StopReason::EXTERNAL, run.result.stop->reason) << what;
                EXPECT_EQ(StopPhase::EXTENSION, run.result.stop->phase) << what;
                EXPECT_EQ("max_memory", run.reason) << what;
            }
        }
        for (uint64_t m = 0; m <= full.store_size; ++m) {
            const std::string what = "case " + std::to_string(case_id) + " max_segments "
                    + std::to_string(m);
            Limits limits;
            limits.max_segments = m;
            limits.stop_at_threshold = true;
            const SegmentRunResult run = run_segments(*c.graph, c.pattern, c.request, limits);
            check_stopped(full, run, what);
            if (m < full.store_size) {
                ASSERT_TRUE(run.result.stop) << what;
                EXPECT_EQ(StopReason::EXTERNAL, run.result.stop->reason) << what;
                EXPECT_EQ("max_segments", run.reason) << what;
                EXPECT_EQ(m + 1, run.store_size) << what;
            } else {
                EXPECT_FALSE(run.result.stop) << what;
            }
        }
        ++case_id;
    }
}

TEST(PatternSegments, LongLeftArmReadsTheClock) {
    // a root whose left arm walks 2,989 k-mers (the first window below the floor, so that the
    // last window's anchor is a root): the clock is read before every 64th node expanded, far
    // more often than the 4,096-step stride of the charges
    std::mt19937 rng(21);
    const std::string record = random_dna(rng, 3000);
    std::string text = record;
    for (size_t i = 0; i < 11; ++i) {
        text[i] = 'N';
    }
    auto graph = build_graph<DBGSuccinct>(11, { record }, DeBruijnGraph::BASIC);
    const Pattern pattern = Pattern::parse(PatternKind::IUPAC, text);
    Request request = segment_request(plan_of(11, 1, WindowChoice::ENDS), Strands::FORWARD, 1000);
    request.min_information_bits = 10;
    CountingClock clock;
    Limits limits;
    limits.deadline = [&]() { return clock.deadline(); };
    const SegmentRunResult run = run_segments(*graph, pattern, request, limits);
    ASSERT_FALSE(run.result.stop);
    EXPECT_EQ(segments_of({ { Orientation::FORWARD, 0, 3000, 1 } }), run.segments);
    // one incoming edge per k-mer but where a 10-mer repeats by chance
    EXPECT_GE(run.result.work.left_edges, 2989u);
    EXPECT_LE(run.result.work.left_edges, 2989u + 20);
    // every left node costs 1 + 2 steps: by the steps alone about two readings
    EXPECT_GE(*clock.readings, 2989u / Budget::kReleaseClockStride);
}

TEST(PatternSegments, EveryAnchorReadsTheClock) {
    // a window of ~2,000 anchors whose arms end after a step or two (N^7 then three bases on a
    // random record): the clock is read before every anchor, so at least once per anchor
    // listed; by the arms' and the listing's strides alone about one in 64
    std::mt19937 rng(22);
    const std::string record = random_dna(rng, 2000);
    auto graph = build_graph<DBGSuccinct>(7, { record }, DeBruijnGraph::BASIC);
    const Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NNNNNNNACG");
    const Request request = segment_request(plan_of(0, 1, WindowChoice::FIRST), Strands::FORWARD,
                                            5000);
    CountingClock clock;
    Limits limits;
    limits.deadline = [&]() { return clock.deadline(); };
    const SegmentRunResult run = run_segments(*graph, pattern, request, limits);
    ASSERT_FALSE(run.result.stop);
    const WindowCount &w = run.result.segments->windows.at(Orientation::FORWARD).front();
    ASSERT_EQ(WindowState::EXTENDED, w.state);
    EXPECT_GT(w.anchors.value, 1800u);
    EXPECT_GE(*clock.readings, w.roots + w.owned);
}

#endif // _DNA_GRAPH || _DNA5_GRAPH

} // namespace
