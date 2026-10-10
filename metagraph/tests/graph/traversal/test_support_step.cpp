#include "gtest/gtest.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include "graph/traversal/support_step.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"


// The shared support step (TESTS §2.2, and the strand rule). Every expectation comes from an
// oracle that never asks the module:
//  - the WALKER'S RULE, copied from walker.cpp (process_item): a coordinate x of the next
//    k-mer continues a chain when x - 1 is live (right arm) or x + 1 (left arm), one binary
//    search per coordinate — the rule the run merges must reproduce;
//  - a RECORD SCAN over explicit records in labelled columns: a label carries a walk when
//    every k-mer of the walk is a substring of one of its records (label level), and holds it
//    at the column coordinates where one of its records contains the whole walk (record
//    level) — string searches, never consecutive coordinates. In the reverse-complement
//    orientation the same scan of the reverse-complement walk.
// The rows the frames read are spelled from the same records (each record's k-mers numbered
// from its column's running count, as MetaGraph numbers coordinates), and once from a real
// annotated graph read through LabelQuery.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::graph::traversal::support_step;

// ------------------------------------------------------------------ the test's own tables

char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return c;
    }
}

std::string rc(const std::string &s) {
    std::string out(s.rbegin(), s.rend());
    for (char &c : out) {
        c = complement(c);
    }
    return out;
}

std::string random_dna(std::mt19937 &gen, size_t n) {
    std::string s(n, 'A');
    for (char &c : s) {
        c = "ACGT"[gen() % 4];
    }
    return s;
}

// the coordinates first..last; a run that is empty or wrapped (a start moved below 0) is a
// failure, and stops the expansion instead of looping until the coordinates wrap around
void expand_run(Coord first, Coord last, std::vector<Coord> *out) {
    if (first > last || last == std::numeric_limits<Coord>::max()) {
        ADD_FAILURE() << "a run [" << first << ", " << last << "]";
        return;
    }
    for (Coord c = first; c <= last; ++c) {
        out->push_back(c);
    }
}

std::vector<Coord> expand(const std::vector<CoordRun> &runs) {
    std::vector<Coord> out;
    for (const CoordRun &r : runs) {
        expand_run(r.first, r.last, &out);
    }
    return out;
}

std::vector<Coord> expand(const ChainRun *runs, size_t n) {
    std::vector<Coord> out;
    for (size_t i = 0; i < n; ++i) {
        expand_run(runs[i].first, runs[i].last, &out);
    }
    return out;
}

// a random ascending set of coordinates in [0, n), clustered into runs
std::vector<Coord> random_coords(std::mt19937 &gen, Coord n, double density) {
    std::vector<Coord> out;
    std::bernoulli_distribution start(density), stay(0.7);
    bool in = false;
    for (Coord c = 0; c < n; ++c) {
        in = in ? stay(gen) : start(gen);
        if (in)
            out.push_back(c);
    }
    return out;
}

// The walker's step rule (walker.cpp, process_item), the oracle of the run merge
std::vector<Coord> walker_rule(Arm side, const std::vector<Coord> &live,
                               const std::vector<Coord> &coords) {
    std::vector<Coord> out;
    for (Coord x : coords) {
        bool ok = side == Arm::RIGHT
            ? (x > 0 && std::binary_search(live.begin(), live.end(), x - 1))
            : std::binary_search(live.begin(), live.end(), x + 1);
        if (ok)
            out.push_back(x);
    }
    return out;
}

// ------------------------------------------------------------------ the records

/**
 * Explicit records in labelled columns (label l: column l). The rows are spelled from them:
 * every k-window of every record, at the column coordinate of its record's first k-mer plus
 * its position; the record ranges are the records' coordinate ranges. The oracles below scan
 * the same records as strings.
 */
struct Records {
    size_t k = 0;
    std::vector<std::vector<std::string>> of;       // label -> its records
    std::vector<std::vector<Coord>> offset;         // label -> record -> its first coordinate
    std::map<std::string, std::map<LabelId, std::vector<Coord>>> rows;

    Records(size_t k, std::vector<std::vector<std::string>> records)
          : k(k), of(std::move(records)), offset(of.size()) {
        for (LabelId l = 0; l < of.size(); ++l) {
            Coord next = 0;
            for (const std::string &r : of[l]) {
                EXPECT_GE(r.size(), k);
                offset[l].push_back(next);
                for (size_t i = 0; i + k <= r.size(); ++i) {
                    rows[r.substr(i, k)][l].push_back(next + i);
                }
                next += r.size() - k + 1;
            }
        }
    }

    size_t num_labels() const { return of.size(); }

    // the row of |kmer| (no labels when no record holds it)
    RowRuns row(const std::string &kmer, bool coords) const {
        RowRuns out(coords);
        auto it = rows.find(kmer);
        if (it == rows.end())
            return out;
        for (const auto &[label, c] : it->second) {
            out.add(label, c.data(), c.size());
        }
        return out;
    }

    RecordRange record(LabelId label, Coord c) const {
        for (size_t i = 0; i < of[label].size(); ++i) {
            const Coord first = offset[label][i];
            const Coord last = first + of[label][i].size() - k;
            if (first <= c && c <= last)
                return RecordRange { first, last };
        }
        ADD_FAILURE() << "no record of label " << label << " holds coordinate " << c;
        return RecordRange { c, c };
    }

    RecordOf record_of() const {
        return [this](LabelId label, Coord c) { return record(label, c); };
    }

    // every k-mer of the records, with its reverse complement (the walks of a BASIC graph
    // and of its mirror)
    std::vector<std::string> kmers_both_strands() const {
        std::set<std::string> out;
        for (const auto &[kmer, labels] : rows) {
            out.insert(kmer);
            out.insert(rc(kmer));
        }
        return { out.begin(), out.end() };
    }
};

// The labels carrying every k-mer of |walk|: each k-window a substring of one of their records
std::vector<LabelId> carried(const Records &idx, const std::string &walk) {
    std::vector<LabelId> out;
    for (LabelId l = 0; l < idx.num_labels(); ++l) {
        bool all = true;
        for (size_t i = 0; all && i + idx.k <= walk.size(); ++i) {
            const std::string w = walk.substr(i, idx.k);
            all = std::any_of(idx.of[l].begin(), idx.of[l].end(), [&](const std::string &r) {
                return r.find(w) != std::string::npos;
            });
        }
        if (all)
            out.push_back(l);
    }
    return out;
}

// The column coordinates at which a record of |label| holds |walk| whole (its first k-mer's)
std::vector<Coord> occurrences(const Records &idx, LabelId label, const std::string &walk) {
    std::vector<Coord> out;
    for (size_t i = 0; i < idx.of[label].size(); ++i) {
        const std::string &r = idx.of[label][i];
        for (size_t p = r.find(walk); p != std::string::npos; p = r.find(walk, p + 1)) {
            out.push_back(idx.offset[label][i] + p);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

// the k-mer whose row a frame of |orientation| reads where its walk spells |kmer|
std::string row_kmer(Orientation orientation, const std::string &kmer) {
    return orientation == Orientation::SPELLED ? kmer : rc(kmer);
}

/**
 * Walks |walk| in |orientation| at |level| on |arm|: the frame opened at the arm's first
 * k-mer, then stepped base by base, each with the row of the k-mer its orientation names;
 * |check| sees every frame with the part of the walk spelled so far, and every step's units
 * are checked against their price and against the rule of SPEC-DRAFT §20.7 (the labels of both
 * lists, the runs of the labels on both). Returns the last frame.
 */
std::unique_ptr<Frame> walk_frames(const Records &idx, const std::string &walk,
                                   Orientation orientation, Arm arm, Support level,
                                   const std::function<void(const Frame&, const std::string&)> &check) {
    const size_t k = idx.k;
    const bool coords = level == Support::TRACE;
    const size_t n = walk.size() - k + 1;
    auto frame = std::make_unique<Frame>();
    auto next = std::make_unique<Frame>();
    std::string kmer = arm == Arm::RIGHT ? walk.substr(0, k) : walk.substr(n - 1);
    RowRuns row = idx.row(row_kmer(orientation, kmer), coords);
    frame->open(orientation, arm, level, kmer, row_kmer(orientation, kmer), row, idx.record_of());
    EXPECT_TRUE(frame->complete());
    EXPECT_EQ(0u, frame->depth());
    check(*frame, kmer);
    for (size_t s = 1; s < n; ++s) {
        const char base = arm == Arm::RIGHT ? walk[k - 1 + s] : walk[n - 1 - s];
        kmer = arm == Arm::RIGHT ? kmer.substr(1) + base : base + kmer.substr(0, k - 1);
        row = idx.row(row_kmer(orientation, kmer), coords);
        const Price price = frame->price(row);
        uint64_t expected = frame->labels().size() + row.num_labels();
        if (coords) {
            for (size_t i = 0; i < frame->labels().size(); ++i) {
                const auto &b = row.labels();
                auto it = std::lower_bound(b.begin(), b.end(), frame->labels()[i]);
                if (it != b.end() && *it == frame->labels()[i])
                    expected += frame->num_chain_runs(i) + row.num_runs(it - b.begin());
            }
        }
        const uint64_t u = frame->step(base, row_kmer(orientation, kmer), row, next.get());
        EXPECT_EQ(expected, u);
        EXPECT_LE(u, price.units);
        EXPECT_LE(next->bytes(), price.bytes);
        EXPECT_TRUE(next->complete());
        EXPECT_EQ(kmer, next->kmer());
        EXPECT_EQ(s, next->depth());
        std::swap(frame, next);
        check(*frame, arm == Arm::RIGHT ? walk.substr(0, k + s) : walk.substr(n - 1 - s));
    }
    return frame;
}

// the frame against the record scan, for the part of the walk it has spelled
void check_against_scan(const Records &idx, const Frame &f, const std::string &part,
                        const std::string &context) {
    const std::string held = f.orientation() == Orientation::SPELLED ? part : rc(part);
    ASSERT_EQ(carried(idx, held), f.labels()) << context << " " << part;
    bool any = false;
    if (f.level() == Support::TRACE) {
        for (size_t i = 0; i < f.labels().size(); ++i) {
            std::vector<CoordRun> starts;
            EXPECT_EQ(f.num_chain_runs(i), f.starts(i, &starts));
            const std::vector<Coord> expected = occurrences(idx, f.labels()[i], held);
            EXPECT_EQ(expected, expand(starts)) << context << " " << part << " label "
                                                << f.labels()[i];
            EXPECT_EQ(expected.size(), f.num_alive(i)) << context << " " << part << " label "
                                                       << f.labels()[i];
            any |= !expected.empty();
        }
    } else {
        any = !f.labels().empty();
    }
    EXPECT_EQ(any, f.supported()) << context << " " << part;
}

// combine() of the two last frames of |walk| against the scan: every label carried in either
// orientation, ascending, with how it holds the walk as spelled and as its reverse complement
// (record_verified where one of its records holds that strand whole and the level is TRACE)
void check_combine(const Records &idx, const std::string &walk, Support level,
                   const Frame &spelled, const Frame &reverse_complement,
                   const std::string &context) {
    auto held = [&](LabelId label, const std::string &w) {
        const std::vector<LabelId> c = carried(idx, w);
        if (!std::binary_search(c.begin(), c.end(), label))
            return Held::NONE;
        return level == Support::TRACE && !occurrences(idx, label, w).empty()
                ? Held::RECORD_VERIFIED : Held::LABEL_INTERSECTION;
    };
    std::vector<LabelSupport> expected;
    for (LabelId l = 0; l < idx.num_labels(); ++l) {
        const LabelSupport s { l, held(l, walk), held(l, rc(walk)) };
        if (s.spelled != Held::NONE || s.reverse_complement != Held::NONE)
            expected.push_back(s);
    }
    std::vector<LabelSupport> both;
    EXPECT_EQ(spelled.labels().size() + reverse_complement.labels().size(),
              combine(spelled, reverse_complement, &both)) << context;
    ASSERT_EQ(expected.size(), both.size()) << context;
    for (size_t i = 0; i < both.size(); ++i) {
        EXPECT_EQ(expected[i].label, both[i].label) << context;
        EXPECT_STREQ(to_string(expected[i].spelled), to_string(both[i].spelled))
                << context << " label " << both[i].label;
        EXPECT_STREQ(to_string(expected[i].reverse_complement),
                     to_string(both[i].reverse_complement))
                << context << " label " << both[i].label;
    }
}

const char* name(Arm arm) { return arm == Arm::RIGHT ? "right" : "left"; }

// ------------------------------------------------------------------ the primitives

TEST(SupportStep, IntersectLabelsAgainstSetIntersection) {
    std::mt19937 gen(5);
    for (size_t t = 0; t < 2000; ++t) {
        std::set<LabelId> sa, sb;
        const size_t na = gen() % 30, nb = gen() % 30;
        const LabelId universe = 1 + gen() % 60;
        while (sa.size() < std::min<size_t>(na, universe)) {
            sa.insert(gen() % universe);
        }
        while (sb.size() < std::min<size_t>(nb, universe)) {
            sb.insert(gen() % universe);
        }
        std::vector<LabelId> a(sa.begin(), sa.end()), b(sb.begin(), sb.end()), expected, got;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                              std::back_inserter(expected));
        EXPECT_EQ(a.size() + b.size(), intersect_labels(a.data(), a.size(), b.data(), b.size(), &got));
        EXPECT_EQ(expected, got);
    }
}

TEST(SupportStep, CoordinateRunsAreMaximal) {
    std::mt19937 gen(7);
    for (size_t t = 0; t < 500; ++t) {
        const std::vector<Coord> coords = random_coords(gen, 1 + gen() % 200, 0.2);
        std::vector<CoordRun> runs;
        EXPECT_EQ(coords.size(), coordinate_runs(coords.data(), coords.size(), &runs));
        EXPECT_EQ(coords, expand(runs));
        for (size_t i = 1; i < runs.size(); ++i) {
            EXPECT_GT(runs[i].first, runs[i - 1].last + 1);     // maximal: a gap between
        }
    }
    std::vector<CoordRun> runs;
    const std::vector<Coord> duplicate { 1, 2, 2, 3 }, descending { 5, 4 };
    EXPECT_THROW(coordinate_runs(duplicate.data(), duplicate.size(), &runs), std::invalid_argument);
    EXPECT_THROW(coordinate_runs(descending.data(), descending.size(), &runs), std::invalid_argument);
    // a homopolymer's row is one run
    std::vector<Coord> homopolymer(29970);
    std::iota(homopolymer.begin(), homopolymer.end(), 0);
    runs.clear();
    coordinate_runs(homopolymer.data(), homopolymer.size(), &runs);
    ASSERT_EQ(1u, runs.size());
    EXPECT_EQ((CoordRun { 0, 29969 }), runs[0]);
}

TEST(SupportStep, DirectionsOfTheArmsAndOrientations) {
    EXPECT_EQ(Direction::UP, direction(Arm::RIGHT, Orientation::SPELLED));
    EXPECT_EQ(Direction::DOWN, direction(Arm::LEFT, Orientation::SPELLED));
    EXPECT_EQ(Direction::DOWN, direction(Arm::RIGHT, Orientation::REVERSE_COMPLEMENT));
    EXPECT_EQ(Direction::UP, direction(Arm::LEFT, Orientation::REVERSE_COMPLEMENT));
}

// Without record bounds the merge of runs is the walker's rule, on both arms
TEST(SupportStep, ContinueChainsIsTheWalkersRule) {
    std::mt19937 gen(11);
    for (size_t t = 0; t < 3000; ++t) {
        const Coord n = 1 + gen() % 120;
        const std::vector<Coord> live = random_coords(gen, n, 0.15 + 0.1 * (t % 5));
        const std::vector<Coord> coords = random_coords(gen, n, 0.15 + 0.1 * (t % 7));
        std::vector<CoordRun> live_runs, row;
        coordinate_runs(live.data(), live.size(), &live_runs);
        coordinate_runs(coords.data(), coords.size(), &row);
        for (Arm side : { Arm::RIGHT, Arm::LEFT }) {
            const Direction d = direction(side, Orientation::SPELLED);
            std::vector<ChainRun> chains, next;
            open_chains(d, 0, live_runs.data(), live_runs.size(), RecordOf(), &chains);
            EXPECT_EQ(live, expand(chains.data(), chains.size()));
            EXPECT_EQ(chains.size() + row.size(),
                      continue_chains(d, chains.data(), chains.size(), row.data(), row.size(),
                                      &next));
            EXPECT_EQ(walker_rule(side, live, coords), expand(next.data(), next.size()))
                    << name(side) << " trial " << t;
            for (size_t i = 1; i < next.size(); ++i) {
                EXPECT_GT(next[i].first, next[i - 1].last);
            }
        }
    }
}

// With record bounds a chain continues only inside its record: the walker's rule restricted
// to x - 1 (x + 1) and x in one record, the records a random partition of the coordinates
TEST(SupportStep, ContinueChainsStopsAtRecordEnds) {
    std::mt19937 gen(13);
    for (size_t t = 0; t < 3000; ++t) {
        const Coord n = 1 + gen() % 120;
        // the records: [ends[i - 1] + 1, ends[i]]
        std::vector<Coord> ends;
        for (Coord c = 0; c < n; ++c) {
            if (c + 1 == n || gen() % 6 == 0)
                ends.push_back(c);
        }
        auto record_index = [&](Coord c) {
            return std::lower_bound(ends.begin(), ends.end(), c) - ends.begin();
        };
        RecordOf record_of = [&](LabelId, Coord c) {
            const size_t i = record_index(c);
            return RecordRange { i ? ends[i - 1] + 1 : 0, ends[i] };
        };
        const std::vector<Coord> live = random_coords(gen, n, 0.3);
        const std::vector<Coord> coords = random_coords(gen, n, 0.3);
        std::vector<CoordRun> live_runs, row;
        coordinate_runs(live.data(), live.size(), &live_runs);
        coordinate_runs(coords.data(), coords.size(), &row);
        for (Arm side : { Arm::RIGHT, Arm::LEFT }) {
            const Direction d = direction(side, Orientation::SPELLED);
            std::vector<ChainRun> chains, next;
            const uint64_t opened = open_chains(d, 0, live_runs.data(), live_runs.size(),
                                                record_of, &chains);
            // one run per (run, record) piece, one lookup per piece
            EXPECT_EQ(live_runs.size() + chains.size(), opened);
            EXPECT_EQ(live, expand(chains.data(), chains.size()));
            for (const ChainRun &c : chains) {
                EXPECT_EQ(record_index(c.first), record_index(c.last));
                const RecordRange r = record_of(0, c.first);
                EXPECT_EQ(d == Direction::UP ? r.last : r.first, c.record_end);
            }
            continue_chains(d, chains.data(), chains.size(), row.data(), row.size(), &next);
            std::vector<Coord> expected;
            for (Coord x : walker_rule(side, live, coords)) {
                const Coord from = side == Arm::RIGHT ? x - 1 : x + 1;
                if (record_index(from) == record_index(x))
                    expected.push_back(x);
            }
            EXPECT_EQ(expected, expand(next.data(), next.size())) << name(side) << " trial " << t;
        }
    }
}

// ------------------------------------------------------------------ frames against the scan

// Random records with planted repeats (copies of other records' segments on either strand,
// homopolymers): the walks are record substrings on both strands, mosaics through shared
// k-mers and random strings; every frame of every prefix (right arm) or suffix (left arm), in
// both orientations and at both levels, equals the record scan
TEST(SupportStep, RandomRecordsAgainstAScanOfTheRecords) {
    std::mt19937 gen(17);
    size_t frames_checked = 0, supported_traces = 0, carried_not_verified = 0, mixed_frames = 0;
    for (size_t t = 0; t < 150; ++t) {
        const size_t k = 3 + gen() % 5;
        const size_t num_labels = 2 + gen() % 4;
        std::vector<std::vector<std::string>> records(num_labels);
        std::vector<std::string> all;
        for (auto &column : records) {
            const size_t num_records = 1 + gen() % 4;
            for (size_t r = 0; r < num_records; ++r) {
                std::string s;
                while (s.size() < k + gen() % 25) {
                    switch (gen() % 4) {
                        case 0:
                            s += random_dna(gen, 1 + gen() % 6);
                            break;
                        case 1:
                            s += std::string(2 + gen() % (k + 2), "ACGT"[gen() % 4]);
                            break;
                        default:
                            if (all.empty()) {
                                s += random_dna(gen, k);
                            } else {
                                const std::string &src = all[gen() % all.size()];
                                const size_t len = std::min(src.size(), k + gen() % (2 * k));
                                const std::string seg = src.substr(gen() % (src.size() - len + 1), len);
                                s += gen() % 2 ? seg : rc(seg);
                            }
                    }
                }
                column.push_back(s);
                all.push_back(s);
            }
        }
        const Records idx(k, records);

        std::vector<std::string> walks;
        for (size_t w = 0; w < 6; ++w) {
            const std::string &src = all[gen() % all.size()];
            const size_t len = std::min(src.size(), k + gen() % (2 * k + 1));
            const std::string seg = src.substr(gen() % (src.size() - len + 1), len);
            walks.push_back(gen() % 2 ? seg : rc(seg));
        }
        const std::vector<std::string> kmers = idx.kmers_both_strands();
        const std::set<std::string> kmer_set(kmers.begin(), kmers.end());
        for (size_t w = 0; w < 6; ++w) {
            std::string walk = kmers[gen() % kmers.size()];
            const size_t target = k + gen() % (2 * k + 1);
            while (walk.size() < target) {
                std::string options;
                for (char b : std::string("ACGT")) {
                    if (kmer_set.count(walk.substr(walk.size() - k + 1) + b))
                        options.push_back(b);
                }
                if (options.empty())
                    break;
                walk.push_back(options[gen() % options.size()]);
            }
            walks.push_back(walk);
        }
        walks.push_back(random_dna(gen, k + gen() % (2 * k)));

        for (const std::string &walk : walks) {
            for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
                for (Support level : { Support::KMER, Support::TRACE }) {
                    std::unique_ptr<Frame> last[2];
                    for (Orientation o : { Orientation::SPELLED, Orientation::REVERSE_COMPLEMENT }) {
                        const std::string context = std::string(name(arm)) + " "
                                + to_string(o) + " " + to_string(level) + " k=" + std::to_string(k)
                                + " trial " + std::to_string(t) + " walk " + walk;
                        const size_t side = o == Orientation::SPELLED ? 0 : 1;
                        last[side] = walk_frames(idx, walk, o, arm, level,
                                                 [&](const Frame &f, const std::string &part) {
                            ++frames_checked;
                            check_against_scan(idx, f, part, context);
                            if (level == Support::TRACE && f.supported())
                                ++supported_traces;
                            if (level == Support::TRACE && !f.supported() && !f.labels().empty())
                                ++carried_not_verified;
                        });
                    }
                    // the two orientations meet at the walk's end, label by label
                    check_combine(idx, walk, level, *last[0], *last[1],
                                  std::string(name(arm)) + " " + to_string(level) + " k="
                                      + std::to_string(k) + " trial " + std::to_string(t)
                                      + " walk " + walk);
                    // a frame with a verified label beside a label carried only
                    for (const Frame *f : { last[0].get(), last[1].get() }) {
                        if (level != Support::TRACE)
                            break;
                        bool verified = false, carried_only = false;
                        for (size_t i = 0; i < f->labels().size(); ++i) {
                            (f->num_chain_runs(i) ? verified : carried_only) = true;
                        }
                        mixed_frames += verified && carried_only;
                    }
                }
            }
        }
    }
    // the trials reach every case: supported walks, and labels carried but held whole by no
    // record (mosaics), also beside a label of the same frame that is verified
    EXPECT_GT(frames_checked, 10000u);
    EXPECT_GT(supported_traces, 1000u);
    EXPECT_GT(carried_not_verified, 100u);
    EXPECT_GT(mixed_frames, 10u);
}

// The design's k = 3 example (SPEC-DRAFT §20.9): a label holding the records ACG and CGT
// supports the path ACGT at the label level, but no record holds it whole — the chain from
// ACG's coordinate 0 to CGT's 1 crosses from one record into the next, and only the record end
// tells
TEST(SupportStep, TwoRecordsOfOneLabelK3) {
    const Records idx(3, { { "ACG", "CGT" } });
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        auto label = walk_frames(idx, "ACGT", Orientation::SPELLED, arm, Support::KMER,
                                 [](const Frame&, const std::string&) {});
        EXPECT_EQ((std::vector<LabelId>{ 0 }), label->labels());
        EXPECT_TRUE(label->supported());

        auto trace = walk_frames(idx, "ACGT", Orientation::SPELLED, arm, Support::TRACE,
                                 [&](const Frame &f, const std::string &part) {
            check_against_scan(idx, f, part, name(arm));
        });
        EXPECT_EQ((std::vector<LabelId>{ 0 }), trace->labels());     // carried
        EXPECT_EQ(0u, trace->num_chain_runs(0));                     // not verified
        EXPECT_EQ(0u, trace->num_alive(0));
        EXPECT_FALSE(trace->supported());
    }
    // the same chain without the record mapping survives: the coordinates 0 and 1 are
    // consecutive, the records are not one
    Frame open, next;
    open.open(Orientation::SPELLED, Arm::RIGHT, Support::TRACE, "ACG", "ACG",
              idx.row("ACG", true));
    ASSERT_EQ(1u, open.anchor()->runs.size());
    EXPECT_EQ((ChainRun { 0, 0, unbounded(Direction::UP) }), open.anchor()->runs[0]);
    open.step('T', "CGT", idx.row("CGT", true), &next);
    EXPECT_TRUE(next.supported());
    EXPECT_EQ(1u, next.num_alive(0));
    std::vector<CoordRun> starts;
    EXPECT_EQ(1u, next.starts(0, &starts));
    EXPECT_EQ((std::vector<CoordRun>{ { 0, 0 } }), starts);
}

// The strand rule: a label holding a walk's first k-mers on + and its last on - does not
// support it in either orientation, although every k-mer of it is on one of the label's
// strands (the per-k-mer union of a k-mer's row and its reverse complement's, which the rule
// forbids, would carry it)
TEST(SupportStep, FirstKmersOnPlusLastOnMinusSupportNeither) {
    const std::string walk = "AACGG";                // AAC, ACG, CGG
    // AAC and ACG as deposited; CGG only as its reverse complement CCG
    const Records idx(3, { { "AACG", "CCG" } });
    // the per-k-mer union carries label 0 on every k-mer
    for (size_t i = 0; i + 3 <= walk.size(); ++i) {
        const std::string x = walk.substr(i, 3);
        EXPECT_TRUE(idx.rows.count(x) || idx.rows.count(rc(x))) << x;
    }
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        for (Support level : { Support::KMER, Support::TRACE }) {
            auto f = walk_frames(idx, walk, Orientation::SPELLED, arm, level,
                                 [&](const Frame &f, const std::string &part) {
                check_against_scan(idx, f, part, name(arm));
            });
            auto r = walk_frames(idx, walk, Orientation::REVERSE_COMPLEMENT, arm, level,
                                 [&](const Frame &f, const std::string &part) {
                check_against_scan(idx, f, part, name(arm));
            });
            EXPECT_FALSE(f->supported());
            EXPECT_TRUE(f->labels().empty());
            EXPECT_FALSE(r->supported());
            EXPECT_TRUE(r->labels().empty());
            std::vector<LabelSupport> both;
            combine(*f, *r, &both);
            EXPECT_TRUE(both.empty());
        }
    }
    // the control: a label holding the whole walk on its - strand only supports it as the
    // reverse complement, at both levels, and combine states which orientation
    const Records minus(3, { { "AACG", "CCG" }, { "T" + rc(walk) + "T" } });
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        for (Support level : { Support::KMER, Support::TRACE }) {
            auto f = walk_frames(minus, walk, Orientation::SPELLED, arm, level,
                                 [&](const Frame &f, const std::string &part) {
                check_against_scan(minus, f, part, name(arm));
            });
            auto r = walk_frames(minus, walk, Orientation::REVERSE_COMPLEMENT, arm, level,
                                 [&](const Frame &f, const std::string &part) {
                check_against_scan(minus, f, part, name(arm));
            });
            EXPECT_FALSE(f->supported());
            EXPECT_TRUE(r->supported());
            std::vector<LabelSupport> both;
            EXPECT_EQ(f->labels().size() + r->labels().size(), combine(*f, *r, &both));
            const Held held = level == Support::TRACE ? Held::RECORD_VERIFIED
                                                      : Held::LABEL_INTERSECTION;
            EXPECT_EQ((std::vector<LabelSupport>{ { 1, Held::NONE, held } }), both);
            if (level == Support::TRACE) {
                // rc(walk) starts at position 1 of the record
                std::vector<CoordRun> starts;
                r->starts(0, &starts);
                EXPECT_EQ((std::vector<Coord>{ 1 }), expand(starts));
            }
        }
    }
}

// The records of PatternPaths.CrossRecordPathIsNotRecordVerified (k = 5): B's coordinates 3,
// 4, 5 of ACGTACC are consecutive but cross from b0 into b1, so B carries the walk and no
// record of it holds it; A's chain stays inside a0. A walk inside b1 is held by B
TEST(SupportStep, CrossRecordChainDiesAtTheBoundary) {
    const Records idx(5, { { "GGACGTACCAA" }, { "TTTACGTA", "CGTACCGG" }, { "TTGCATGCAT" } });
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        auto f = walk_frames(idx, "ACGTACC", Orientation::SPELLED, arm, Support::TRACE,
                             [&](const Frame &f, const std::string &part) {
            check_against_scan(idx, f, part, name(arm));
        });
        ASSERT_EQ((std::vector<LabelId>{ 0, 1 }), f->labels());
        std::vector<CoordRun> a;
        f->starts(0, &a);
        EXPECT_EQ((std::vector<Coord>{ 2 }), expand(a));            // a0, 1-based 3-9
        EXPECT_EQ(0u, f->num_chain_runs(1));                         // B: carried only

        auto g = walk_frames(idx, "CGTACC", Orientation::SPELLED, arm, Support::TRACE,
                             [&](const Frame &f, const std::string &part) {
            check_against_scan(idx, f, part, name(arm));
        });
        ASSERT_EQ((std::vector<LabelId>{ 0, 1 }), g->labels());
        std::vector<CoordRun> b;
        g->starts(1, &b);
        EXPECT_EQ((std::vector<Coord>{ 4 }), expand(b));            // b1's first k-mer

        // combine states each label's own support: A verified, B beside it carried only
        auto r = walk_frames(idx, "ACGTACC", Orientation::REVERSE_COMPLEMENT, arm,
                             Support::TRACE, [](const Frame&, const std::string&) {});
        std::vector<LabelSupport> both;
        combine(*f, *r, &both);
        EXPECT_EQ((std::vector<LabelSupport>{ { 0, Held::RECORD_VERIFIED, Held::NONE },
                                              { 1, Held::LABEL_INTERSECTION, Held::NONE } }),
                  both);
        check_combine(idx, "ACGTACC", Support::TRACE, *f, *r, name(arm));
    }
}

// Two labels of one frame with different support (k = 4, walk ACGTACC): label 0 carries every
// k-mer in two records (ACGTA, GTACC) and holds the walk in none; label 1 holds it whole. Each
// gets its own Held, in both orientations of the walk (its reverse complement GGTACGT is held
// by no record, its k-mers by none)
TEST(SupportStep, CombineStatesEachLabelsOwnSupport) {
    const Records idx(4, { { "ACGTA", "GTACC" }, { "TTACGTACCTT" } });
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        for (Support level : { Support::KMER, Support::TRACE }) {
            auto f = walk_frames(idx, "ACGTACC", Orientation::SPELLED, arm, level,
                                 [&](const Frame &f, const std::string &part) {
                check_against_scan(idx, f, part, name(arm));
            });
            auto r = walk_frames(idx, "ACGTACC", Orientation::REVERSE_COMPLEMENT, arm, level,
                                 [&](const Frame &f, const std::string &part) {
                check_against_scan(idx, f, part, name(arm));
            });
            std::vector<LabelSupport> both;
            combine(*f, *r, &both);
            const Held verified = level == Support::TRACE ? Held::RECORD_VERIFIED
                                                          : Held::LABEL_INTERSECTION;
            EXPECT_EQ((std::vector<LabelSupport>{ { 0, Held::LABEL_INTERSECTION, Held::NONE },
                                                  { 1, verified, Held::NONE } }), both)
                    << name(arm) << " " << to_string(level);
            check_combine(idx, "ACGTACC", level, *f, *r, name(arm));
        }
    }
    // the mirror: walking the reverse complement, the same labels hold it as reverse_complement
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        auto f = walk_frames(idx, rc("ACGTACC"), Orientation::SPELLED, arm, Support::TRACE,
                             [](const Frame&, const std::string&) {});
        auto r = walk_frames(idx, rc("ACGTACC"), Orientation::REVERSE_COMPLEMENT, arm,
                             Support::TRACE, [](const Frame&, const std::string&) {});
        std::vector<LabelSupport> both;
        combine(*f, *r, &both);
        EXPECT_EQ((std::vector<LabelSupport>{ { 0, Held::NONE, Held::LABEL_INTERSECTION },
                                              { 1, Held::NONE, Held::RECORD_VERIFIED } }), both)
                << name(arm);
    }
}

// A homopolymer: one record of 30,000 As. Its one k-mer's 29,970 coordinates are one run,
// built once with the row; every step then compares one label and one run on each side
TEST(SupportStep, HomopolymerIsOneRunPerStep) {
    const size_t k = 31;
    const Records idx(k, { { std::string(30000, 'A') } });
    const std::string kmer(k, 'A');
    const RowRuns row = idx.row(kmer, true);
    ASSERT_EQ(1u, row.total_runs());
    const size_t steps = 500;
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        Frame frame, next;
        // 1 label, 1 run and its one record lookup
        EXPECT_EQ(3u, frame.open(Orientation::SPELLED, arm, Support::TRACE, kmer, kmer, row,
                                 idx.record_of()));
        uint64_t units = 0;
        for (size_t s = 0; s < steps; ++s) {
            const uint64_t u = frame.step('A', kmer, row, &next);
            EXPECT_EQ(4u, u);            // 1 + 1 labels, 1 + 1 runs: O(1), whatever the length
            EXPECT_EQ(1u, next.total_chain_runs());
            EXPECT_EQ(29970u - (s + 1), next.num_alive(0));
            units += u;
            std::swap(frame, next);
        }
        EXPECT_EQ(4 * steps, units);
        std::vector<CoordRun> starts;
        frame.starts(0, &starts);
        const std::string walk(k + steps, 'A');
        EXPECT_EQ(occurrences(idx, 0, walk), expand(starts));
        EXPECT_EQ((std::vector<CoordRun>{ { 0, 30000 - (k + steps) } }), starts);
        // the other strand: TTT...T is no k-mer of the record
        Frame mirror;
        mirror.open(Orientation::REVERSE_COMPLEMENT, arm, Support::TRACE, kmer, rc(kmer),
                    idx.row(rc(kmer), true), idx.record_of());
        EXPECT_FALSE(mirror.supported());
    }
}

// ------------------------------------------------------------------ the frames below the anchor

// A chain dies at its record's end at depth 5 while its label stays carried (the walk's last
// k-mer lies in another record of the label), in both directions and both orientations: the
// walk ACGTACCAAT (k = 5) in the record GGACGTACCAA up to its last k-mer, which CCAATG holds
// (UP: right arm as spelled, left arm as the reverse complement), and in CGTACCAATGG from its
// second k-mer, its first held by TTACGTAGG (DOWN: left arm as spelled, right arm as the
// reverse complement)
TEST(SupportStep, ChainDiesAtItsRecordEndAtDepth) {
    const std::string walk = "ACGTACCAAT";
    const Records up(5, { { "GGACGTACCAA", "CCAATG" } });
    const Records down(5, { { "CGTACCAATGG", "TTACGTAGG" } });
    struct Case {
        const Records *idx;
        Arm arm;
        Orientation orientation;
        Coord anchor;       // the chain's coordinate when opened
    };
    for (const Case &c : { Case { &up, Arm::RIGHT, Orientation::SPELLED, 2 },
                           Case { &up, Arm::LEFT, Orientation::REVERSE_COMPLEMENT, 2 },
                           Case { &down, Arm::LEFT, Orientation::SPELLED, 4 },
                           Case { &down, Arm::RIGHT, Orientation::REVERSE_COMPLEMENT, 4 } }) {
        const std::string context = std::string(name(c.arm)) + " " + to_string(c.orientation);
        const std::string w = c.orientation == Orientation::SPELLED ? walk : rc(walk);
        auto last = walk_frames(*c.idx, w, c.orientation, c.arm, Support::TRACE,
                                [&](const Frame &f, const std::string &part) {
            check_against_scan(*c.idx, f, part, context);
            ASSERT_EQ((std::vector<LabelId>{ 0 }), f.labels()) << context << " " << part;
            if (f.depth() < 5) {
                EXPECT_EQ(1u, f.num_alive(0)) << context << " " << part;
                EXPECT_EQ(1u, f.num_chain_runs(0)) << context << " " << part;
                EXPECT_TRUE(f.supported()) << context << " " << part;
                std::vector<CoordRun> starts;
                f.starts(0, &starts);
                const Coord start = f.direction() == Direction::UP ? c.anchor
                                                                   : c.anchor - f.depth();
                EXPECT_EQ((std::vector<CoordRun>{ { start, start } }), starts)
                        << context << " " << part;
            } else {
                EXPECT_EQ(0u, f.num_alive(0)) << context << " " << part;
                EXPECT_EQ(0u, f.num_chain_runs(0)) << context << " " << part;
                EXPECT_FALSE(f.supported()) << context << " " << part;
            }
        });
        EXPECT_EQ(5u, last->depth());
        EXPECT_FALSE(last->supported());
    }
}

// A homopolymer run shrinks by one chain per step and to one chain when the walk leaves it:
// the record of 300 As and a C (k = 31) holds the walk A^(31 + s) C once; the frames below
// the anchor hold one bit per chain, so their bytes do not grow with the 270 chains
TEST(SupportStep, HomopolymerRunShrinksToOneChain) {
    const size_t k = 31, s = 100;
    // the record and the walk read from the arm's end: on the left arm C A^(31 + s), read from A^31
    const Records right(k, { { std::string(300, 'A') + "C" } });
    const Records left(k, { { "C" + std::string(300, 'A') } });
    const std::string kmer(k, 'A');
    ASSERT_EQ(1u, right.row(kmer, true).total_runs());
    for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
        const std::string walk = arm == Arm::RIGHT ? std::string(k + s, 'A') + "C"
                                                   : "C" + std::string(k + s, 'A');
        const Records &records = arm == Arm::RIGHT ? right : left;
        auto last = walk_frames(records, walk, Orientation::SPELLED, arm, Support::TRACE,
                                [&](const Frame &f, const std::string &part) {
            check_against_scan(records, f, part, name(arm));
            ASSERT_EQ(1u, f.labels().size());
            EXPECT_EQ(1u, f.num_chain_runs(0));
            EXPECT_EQ(f.depth() <= s ? 270 - f.depth() : 1, f.num_alive(0)) << name(arm) << " " << part;
            if (f.depth())
                EXPECT_EQ(Frame::step_bytes(Support::TRACE, 1, 270, k), f.bytes());
        });
        EXPECT_EQ(s + 1, last->depth());
        std::vector<CoordRun> starts;
        EXPECT_EQ(1u, last->starts(0, &starts));
        EXPECT_EQ(1u, last->num_alive(0));
        // the right arm: A^31 at 0..269, A^30 C at 270, the walk's first k-mer at 270 - (s + 1);
        // the left arm: C A^30 at 0, A^31 at 1..270, the walk's first k-mer (C A^30) at 0
        const Coord start = arm == Arm::RIGHT ? 270 - (s + 1) : 0;
        EXPECT_EQ((std::vector<CoordRun>{ { start, start } }), starts) << name(arm);
    }
}

// The survival bits of every depth equal the run-based merge applied step by step from the
// opened chains: random rows (labels, coordinates), with and without record bounds, in both
// directions; the frames' chain runs are the reference's runs, one by one, their units the
// reference's. The two frames are used in turn, so the anchor's object is overwritten at
// depth 2 while the deeper frames still stand on its chains
TEST(SupportStep, SurvivalBitsEqualTheRunBasedMerge) {
    std::mt19937 gen(23);
    size_t compared = 0, bounded = 0;
    for (size_t t = 0; t < 400; ++t) {
        const size_t num_labels = 1 + gen() % 6;
        const Coord n = 1 + gen() % 150;
        const bool records = gen() % 2;
        // the records: per label, [ends[i - 1] + 1, ends[i]]
        std::vector<std::vector<Coord>> ends(num_labels);
        for (auto &e : ends) {
            for (Coord c = 0; c < n; ++c) {
                if (c + 1 == n || gen() % 8 == 0)
                    e.push_back(c);
            }
        }
        RecordOf record_of;
        if (records) {
            record_of = [&](LabelId label, Coord c) {
                const auto &e = ends[label];
                const size_t i = std::lower_bound(e.begin(), e.end(), c) - e.begin();
                return RecordRange { i ? e[i - 1] + 1 : 0, e[i] };
            };
            ++bounded;
        }
        const size_t steps = 1 + gen() % 12;
        std::vector<RowRuns> rows;
        for (size_t d = 0; d <= steps; ++d) {
            RowRuns row(true);
            for (LabelId l = 0; l < num_labels; ++l) {
                if (gen() % 5 == 0)
                    continue;
                const std::vector<Coord> coords = random_coords(gen, n, 0.2 + 0.1 * (t % 5));
                if (!coords.empty())
                    row.add(l, coords.data(), coords.size());
            }
            rows.push_back(std::move(row));
        }
        for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
            const Direction d = direction(arm, Orientation::SPELLED);
            // the reference: each label's chain runs, continued with continue_chains
            std::vector<LabelId> labels = rows[0].labels();
            std::vector<std::vector<ChainRun>> chains(labels.size());
            for (size_t i = 0; i < labels.size(); ++i) {
                open_chains(d, labels[i], rows[0].runs(i), rows[0].num_runs(i), record_of,
                            &chains[i]);
            }
            Frame f, next;
            f.open(Orientation::SPELLED, arm, Support::TRACE, "AAA", "AAA", rows[0], record_of);
            for (size_t s = 1; s <= steps; ++s) {
                const RowRuns &row = rows[s];
                std::vector<LabelId> next_labels;
                std::vector<std::vector<ChainRun>> next_chains;
                uint64_t expected_units = labels.size() + row.num_labels();
                for (size_t i = 0; i < labels.size(); ++i) {
                    const auto &b = row.labels();
                    auto it = std::lower_bound(b.begin(), b.end(), labels[i]);
                    if (it == b.end() || *it != labels[i])
                        continue;
                    const size_t j = it - b.begin();
                    next_labels.push_back(labels[i]);
                    next_chains.emplace_back();
                    expected_units += continue_chains(d, chains[i].data(), chains[i].size(),
                                                      row.runs(j), row.num_runs(j),
                                                      &next_chains.back());
                }
                EXPECT_EQ(expected_units, f.step('A', "AAA", row, &next)) << name(arm) << " trial " << t;
                std::swap(f, next);
                labels = std::move(next_labels);
                chains = std::move(next_chains);
                ASSERT_EQ(labels, f.labels()) << name(arm) << " trial " << t << " depth " << s;
                const Coord back = d == Direction::UP ? s : 0;
                for (size_t i = 0; i < labels.size(); ++i) {
                    std::vector<CoordRun> starts;
                    EXPECT_EQ(chains[i].size(), f.starts(i, &starts));
                    EXPECT_EQ(chains[i].size(), f.num_chain_runs(i));
                    ASSERT_EQ(chains[i].size(), starts.size()) << name(arm) << " trial " << t;
                    for (size_t r = 0; r < starts.size(); ++r) {
                        EXPECT_EQ((CoordRun { chains[i][r].first - back, chains[i][r].last - back }),
                                  starts[r]) << name(arm) << " trial " << t << " depth " << s;
                    }
                    EXPECT_EQ(expand(chains[i].data(), chains[i].size()).size(), f.num_alive(i));
                    ++compared;
                }
                size_t total = 0;
                for (const auto &c : chains) {
                    total += c.size();
                }
                EXPECT_EQ(total, f.total_chain_runs());
                EXPECT_EQ(total > 0, f.supported());
            }
        }
    }
    EXPECT_GT(compared, 2000u);
    EXPECT_GT(bounded, 100u);
}

// ------------------------------------------------------------------ the API's guards

// A frame takes only the row of the k-mer its orientation names: the spelled k-mer, or its
// reverse complement — never the other one, at the anchor or at a step
TEST(SupportStep, OrientationIsChecked) {
    const Records idx(3, { { "AACGT", "ACGTT" } });
    Frame f, next;
    EXPECT_THROW(f.open(Orientation::SPELLED, Arm::RIGHT, Support::KMER, "AAC", "GTT",
                        idx.row("GTT", false)), std::invalid_argument);
    EXPECT_THROW(f.open(Orientation::REVERSE_COMPLEMENT, Arm::RIGHT, Support::KMER, "AAC", "AAC",
                        idx.row("AAC", false)), std::invalid_argument);
    EXPECT_THROW(f.open(Orientation::SPELLED, Arm::RIGHT, Support::TRACE, "AAC", "AAC",
                        idx.row("AAC", false)), std::invalid_argument);   // no coordinates
    // a row with coordinates whose label has none (read without them) is refused
    RowRuns coordinates(true);
    EXPECT_THROW(coordinates.add(0), std::invalid_argument);
    const std::vector<Coord> one { 7 };
    EXPECT_EQ(2u, coordinates.add(0, one.data(), one.size()));
    f.open(Orientation::SPELLED, Arm::RIGHT, Support::KMER, "AAC", "AAC", idx.row("AAC", false));
    f.step('G', "ACG", idx.row("ACG", false), &next);
    EXPECT_TRUE(next.complete());
    // the walk's next k-mer is ACG: its reverse complement's row is refused, as is a k-mer
    // that does not follow; a refused step leaves no usable frame behind
    EXPECT_THROW(f.step('G', "CGT", idx.row("CGT", false), &next), std::invalid_argument);
    EXPECT_FALSE(next.complete());
    EXPECT_THROW(f.step('G', "ACT", idx.row("ACT", false), &next), std::invalid_argument);
    EXPECT_THROW(f.step('G', "ACG", idx.row("ACG", false), &f), std::invalid_argument);
    f.step('G', "ACG", idx.row("ACG", false), &next);
    EXPECT_TRUE(next.complete());
    EXPECT_EQ("ACG", next.kmer());

    Frame r, r_next;
    r.open(Orientation::REVERSE_COMPLEMENT, Arm::RIGHT, Support::KMER, "AAC", "GTT",
           idx.row("GTT", false));
    EXPECT_THROW(r.step('G', "ACG", idx.row("ACG", false), &r_next), std::invalid_argument);
    r.step('G', "CGT", idx.row("CGT", false), &r_next);     // rc(ACG)
    EXPECT_EQ("ACG", r_next.kmer());     // the frame spells the walk, not its mirror
    // a palindromic k-mer is its own reverse complement: both orientations read its row
    Frame p;
    p.open(Orientation::REVERSE_COMPLEMENT, Arm::LEFT, Support::KMER, "ACGT", "ACGT",
           Records(4, { { "ACGT" } }).row("ACGT", false));
    EXPECT_TRUE(p.supported());

    // combine: a spelled and a reverse-complement frame of one walk, nothing else
    std::vector<LabelSupport> out;
    EXPECT_THROW(combine(r_next, next, &out), std::invalid_argument);
    EXPECT_THROW(combine(next, next, &out), std::invalid_argument);
    EXPECT_THROW(combine(f, r_next, &out), std::invalid_argument);   // depth 0 and 1
    EXPECT_NO_THROW(combine(next, r_next, &out));
}

// The clock is read once every kPollIterations rounds, before the round; a step it stops
// leaves an incomplete frame, which cannot be stepped
TEST(SupportStep, TheClockStopsALongStep) {
    const size_t n = 10000;
    RowRuns row(false);
    for (LabelId l = 0; l < n; ++l) {
        row.add(l);
    }
    Frame f, next, after;
    f.open(Orientation::SPELLED, Arm::RIGHT, Support::KMER, "AAA", "AAA", row);

    Clock patient;
    size_t asked = 0;
    patient.expired = [&]() { ++asked; return false; };
    EXPECT_EQ(2 * n, f.step('A', "AAA", row, &next, &patient));
    EXPECT_TRUE(next.complete());
    EXPECT_EQ(n, next.labels().size());
    EXPECT_EQ(n / kPollIterations, asked);       // 10,000 rounds of the label merge: 2 readings
    EXPECT_EQ(asked, patient.readings);

    Clock expired;
    expired.expired = []() { return true; };
    const uint64_t done = f.step('A', "AAA", row, &next, &expired);
    EXPECT_TRUE(expired.stopped);
    EXPECT_FALSE(next.complete());
    EXPECT_EQ(kPollIterations - 1, next.labels().size());    // the 4,096th round is not done
    EXPECT_EQ(2 * (kPollIterations - 1), done);
    EXPECT_LT(done, f.price(row).units);
    EXPECT_THROW(next.step('A', "AAA", row, &after), std::logic_error);
    // stopped stays stopped
    EXPECT_EQ(0u, f.step('A', "AAA", row, &after, &expired));
    EXPECT_FALSE(after.complete());
}

// ------------------------------------------------------------------ rows of a real index

// The same frames on the rows of an annotated graph with coordinates (row-diff, BASIC), read
// through LabelQuery and turned into runs by RowRuns::assign, with the index's own record
// mapping (LabelOracle::sequence_range): the coordinates and the record ranges are what the
// records say
TEST(SupportStep, RowsOfARealIndex) {
    const size_t k = 5;
    std::vector<std::vector<std::string>> records {
        { "GGACGTACCAA" }, { "TTTACGTA", "CGTACCGG" }, { "TTGCATGCAT" },
    };
    // a fourth column with repeats on both strands of the others
    std::string d = "AAAAAAAAACGTACC" + rc("GGACGTACCAA") + "CGTACG";
    records.push_back({ d, rc(d).substr(3), "CCCCCCCC" });
    const Records idx(k, records);

    std::vector<std::string> seqs, labels;
    std::vector<uint64_t> starts;
    std::vector<std::vector<std::string>> headers(records.size());
    std::vector<std::vector<uint64_t>> num_kmers(records.size());
    for (size_t l = 0; l < records.size(); ++l) {
        for (size_t r = 0; r < records[l].size(); ++r) {
            seqs.push_back(records[l][r]);
            labels.push_back("L" + std::to_string(l));
            starts.push_back(idx.offset[l][r]);
        }
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            k, seqs, labels, DeBruijnGraph::BASIC, true, starts);
    const auto &encoder = anno->get_annotator().get_label_encoder();
    headers.assign(encoder.size(), {});
    num_kmers.assign(encoder.size(), {});
    for (size_t l = 0; l < records.size(); ++l) {
        const size_t c = encoder.encode("L" + std::to_string(l));
        for (size_t r = 0; r < records[l].size(); ++r) {
            headers[c].push_back("r" + std::to_string(l) + "_" + std::to_string(r));
            num_kmers[c].push_back(records[l][r].size() - k + 1);
        }
    }
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);
    std::vector<LabelRef> refs;
    for (size_t l = 0; l < records.size(); ++l) {
        refs.push_back(oracle.resolve_label("L" + std::to_string(l)));
    }
    LabelQuery query(oracle, refs, true);
    auto real_row = [&](const std::string &kmer) {
        RowRuns out(true);
        const node_index key = oracle.keys_of_sequence(kmer)[0];
        if (key != npos)
            out.assign(query.fetch(key), true);
        return out;
    };
    const RecordOf record_of = [&](LabelId l, Coord c) {
        const LabelOracle::SeqRange r = oracle.sequence_range(refs[l].column, c);
        return RecordRange { r.first, r.last };
    };

    // the walks: every substring of 5 to 12 bases of every record, on both strands
    std::set<std::string> walks;
    for (const auto &column : records) {
        for (const std::string &r : column) {
            for (size_t len = k; len <= std::min<size_t>(12, r.size()); ++len) {
                for (size_t p = 0; p + len <= r.size(); ++p) {
                    walks.insert(r.substr(p, len));
                    walks.insert(rc(r.substr(p, len)));
                }
            }
        }
    }
    size_t verified = 0;
    for (const std::string &walk : walks) {
        for (Orientation o : { Orientation::SPELLED, Orientation::REVERSE_COMPLEMENT }) {
            for (Arm arm : { Arm::RIGHT, Arm::LEFT }) {
                const size_t n = walk.size() - k + 1;
                std::string kmer = arm == Arm::RIGHT ? walk.substr(0, k) : walk.substr(n - 1);
                Frame f, next;
                f.open(o, arm, Support::TRACE, kmer, row_kmer(o, kmer), real_row(row_kmer(o, kmer)),
                       record_of);
                for (size_t s = 1; s < n; ++s) {
                    const char base = arm == Arm::RIGHT ? walk[k - 1 + s] : walk[n - 1 - s];
                    kmer = arm == Arm::RIGHT ? kmer.substr(1) + base : base + kmer.substr(0, k - 1);
                    f.step(base, row_kmer(o, kmer), real_row(row_kmer(o, kmer)), &next);
                    std::swap(f, next);
                }
                check_against_scan(idx, f, walk, std::string(name(arm)) + " " + to_string(o));
                verified += f.supported();
            }
        }
    }
    EXPECT_GT(verified, 100u);
}

} // namespace
