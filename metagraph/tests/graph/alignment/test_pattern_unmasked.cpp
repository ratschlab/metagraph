/**
 * The pattern search on graphs WITHOUT the dummy-edge mask (owner decisions #16 and #17 of
 * 2026-10-08; docs/DESIGN-pattern-search.md §4.4): counts that are true upper bounds U (the
 * BOSS entries of the ranges, source dummies included), exact where the engine can tell,
 * exact lists, and the sampled real fraction f behind the route's estimate U x f.
 *
 * Every graph is built twice from the same records, the second copy's mask reset (a TWIN):
 * the same BOSS, the same node ids, one served with and one without the mask. The oracles:
 *  - the graph-walk oracle on the MASKED twin (its valid edges, and on a wrapped PRIMARY graph
 *    the wrapper's reverse complements): the true contexts, with their node ids and order;
 *  - the record-scan oracle: the records' own k-mers, spelled, never the graph;
 *  - the entry-scan oracle of U on the UNMASKED twin: every BOSS entry (W != $) spelled with
 *    its '$' symbols (get_node_sequence), matched against the oriented window where a '$' may
 *    stand only at a position of the window's leading pattern-N run that the engine does not
 *    search (§4.1, "Cost": on $ACGT the run less its last position); a k-mer holding no '$'
 *    is real. It reads neither the mask nor the engine;
 *  - for f: the spelled entries (every node's symbols read by get_node_seq), BOSS's traversal
 *    of the dummy tree (mark_source_dummy_edges, `stats --count-dummy`), and a re-implementation
 *    of the draws (std::mt19937_64 seeded with the number of edges) with the spelled test.
 * The counts must state true relations (U >= the true count, the lower bound <= it), U must
 * equal the entry-scan count exactly where the design says so, an exact count must equal the
 * masked twin's, and every list must equal the masked twin's.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/genetic_code.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

// the nucleotide builds the engine serves, as in test_pattern_search.cpp
#if _DNA_GRAPH || _DNA5_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using mtg::test::build_graph;
using mtg::test::build_graph_batch;

typedef DeBruijnGraph::node_index node_index;
typedef boss::BOSS BOSS;

constexpr uint64_t kManySteps = 1'000'000'000;


// ---------------------------------------------------------------- helpers

std::string rev_comp(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

// the oracles' IUPAC table (as test_pattern_search.cpp's: the bases each code admits)
std::string iupac_bases(char code) {
    switch (code) {
        case 'A': return "A";
        case 'C': return "C";
        case 'G': return "G";
        case 'T': return "T";
        case 'R': return "AG";
        case 'Y': return "CT";
        case 'S': return "CG";
        case 'W': return "AT";
        case 'K': return "GT";
        case 'M': return "AC";
        case 'B': return "CGT";
        case 'D': return "AGT";
        case 'H': return "ACT";
        case 'V': return "ACG";
        case 'N': return "ACGT";
        default: return "";
    }
}

char complement_base(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return '?';
    }
}

typedef std::vector<std::string> Bases;

Bases oracle_pattern(const std::string &text) {
    Bases q;
    for (char c : text) {
        q.push_back(iupac_bases(c));
    }
    return q;
}

Bases oracle_rc(const Bases &q) {
    Bases rc(q.rbegin(), q.rend());
    for (std::string &bases : rc) {
        for (char &b : bases) {
            b = complement_base(b);
        }
        std::sort(bases.begin(), bases.end());
    }
    return rc;
}

std::vector<std::pair<Orientation, Bases>> orientations(const std::string &text,
                                                        Strands strands) {
    const Bases q = oracle_pattern(text);
    if (q == oracle_rc(q))
        return { { Orientation::PALINDROMIC, q } };
    std::vector<std::pair<Orientation, Bases>> result;
    if (strands != Strands::REVERSE)
        result.emplace_back(Orientation::FORWARD, q);
    if (strands != Strands::FORWARD)
        result.emplace_back(Orientation::REVERSE, oracle_rc(q));
    return result;
}

// a k-mer's piece instantiates |w|: every symbol a base its position admits (never $ or N)
bool matches(const Bases &w, std::string_view s) {
    if (s.size() != w.size())
        return false;
    for (size_t i = 0; i < w.size(); ++i) {
        if (w[i].find(s[i]) == std::string::npos)
            return false;
    }
    return true;
}

// the leading pattern-N run of a window that the engine does not search: on a $ACGT graph
// the run, less the window's last position (§4.1, "Cost"); none on $ACGTN
size_t unsearched_lead(const Bases &w) {
#if _DNA5_GRAPH
    (void)w;
    return 0;
#else
    size_t lead = 0;
    while (lead + 1 < w.size() && w[lead] == "ACGT") {
        ++lead;
    }
    return lead;
#endif
}

// an entry's piece is counted into U (§4.1, no mask): as matches(), except that a '$' may
// stand in the unsearched lead (where only a source dummy has one)
bool matches_entry(const Bases &w, std::string_view s) {
    if (s.size() != w.size())
        return false;
    const size_t lead = unsearched_lead(w);
    for (size_t i = 0; i < w.size(); ++i) {
        if (s[i] == '$' ? i >= lead : w[i].find(s[i]) == std::string::npos)
            return false;
    }
    return true;
}

std::vector<uint32_t> scope_offsets(size_t L, size_t k, Scope scope) {
    if (L > k)
        return { 0 };
    if (scope == Scope::SUFFIX)
        return { static_cast<uint32_t>(k - L) };
    std::vector<uint32_t> offsets;
    for (uint32_t p = 0; p + L <= k; ++p) {
        offsets.push_back(p);
    }
    return offsets;
}

const DBGSuccinct& base_dbg(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

std::shared_ptr<DeBruijnGraph> build(size_t k, const std::vector<std::string> &records,
                                     DeBruijnGraph::Mode mode, bool batch) {
    return batch ? build_graph_batch<DBGSuccinct>(k, records, mode)
                 : build_graph<DBGSuccinct>(k, records, mode);
}

/**
 * The same records built twice; the second copy's mask reset. Checks that the two BOSS
 * tables are the same, entry by entry (so that node ids and answer order carry over).
 */
struct Twin {
    std::shared_ptr<DeBruijnGraph> masked;
    std::shared_ptr<DeBruijnGraph> unmasked;
};

Twin build_twin(size_t k, const std::vector<std::string> &records, DeBruijnGraph::Mode mode,
                bool batch = false) {
    Twin twin { build(k, records, mode, batch), build(k, records, mode, batch) };
    DBGSuccinct &unmasked = const_cast<DBGSuccinct&>(base_dbg(*twin.unmasked));
    unmasked.reset_mask();
    const DBGSuccinct &masked = base_dbg(*twin.masked);
    EXPECT_TRUE(masked.get_mask());
    EXPECT_FALSE(unmasked.get_mask());
    EXPECT_EQ(masked.max_index(), unmasked.max_index());
    for (node_index e = 1; e <= std::min(masked.max_index(), unmasked.max_index()); ++e) {
        EXPECT_EQ(masked.get_boss().get_node_str(e) + std::to_string(masked.get_boss().get_W(e)),
                  unmasked.get_boss().get_node_str(e)
                        + std::to_string(unmasked.get_boss().get_W(e)));
    }
    return twin;
}

Request make_request(Scope scope = Scope::ANY_OFFSET, Strands strands = Strands::BOTH) {
    Request request;
    request.scope = scope;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_contexts = 1'000'000'000;
    request.max_anchors = 1'000'000'000;
    request.max_paths = 1'000'000'000;
    return request;
}

struct Ctx {
    node_index node;
    uint32_t offset;
    Orientation orientation;
    node_index base_node;

    bool operator<(const Ctx &o) const {
        return std::tie(node, offset, orientation) < std::tie(o.node, o.offset, o.orientation);
    }
    bool operator==(const Ctx &o) const {
        return std::tie(node, offset, orientation, base_node)
                == std::tie(o.node, o.offset, o.orientation, o.base_node);
    }
};

std::ostream& operator<<(std::ostream &out, const Ctx &c) {
    return out << "(" << c.node << ", " << c.offset << ", " << orientation_key(c.orientation)
               << ")";
}

Result count_of(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                uint64_t max_steps = kManySteps) {
    Budget budget(max_steps, Deadline::unbounded());
    return PatternSearch(graph).count(pattern, request, budget);
}

std::vector<Ctx> enumerate_of(const DeBruijnGraph &graph, const Pattern &pattern,
                              const Request &request, Result *result,
                              uint64_t max_steps = kManySteps) {
    Budget budget(max_steps, Deadline::unbounded());
    std::vector<Ctx> contexts;
    *result = PatternSearch(graph).enumerate(pattern, request, budget, [&](const Context &c) {
        contexts.push_back(Ctx { c.node, c.offset, c.orientation, c.base_node });
    });
    return contexts;
}

bool has_note(const Result &result, const std::string &note) {
    return std::count(result.notes.begin(), result.notes.end(), note) > 0;
}

// |count| states a true relation to |truth| (§3)
::testing::AssertionResult true_relation(const Count &count, uint64_t truth) {
    bool holds = false;
    switch (count.relation) {
        case Relation::EXACT:
            holds = count.value == truth;
            break;
        case Relation::AT_LEAST:
            holds = count.value <= truth;
            break;
        case Relation::BOUNDS:
            holds = count.lower == count.value && count.lower <= truth && truth <= count.upper;
            break;
        case Relation::UNKNOWN:
            holds = count.value == 0;
            break;
    }
    if (holds)
        return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure()
        << to_string(count.relation) << " " << count.value << " [" << count.lower << ", "
        << count.upper << "] for the true " << truth;
}

// the upper end of an EXACT or BOUNDS count
uint64_t upper_of(const Count &count) {
    return count.relation == Relation::BOUNDS ? count.upper : count.value;
}


// ---------------------------------------------------------------- the oracles

// every k-mer of the masked twin with its node id (on a wrapped PRIMARY graph also the
// wrapper's reverse complements): the true k-mers, ids valid on the unmasked twin too
std::vector<std::pair<node_index, std::string>> true_kmers(const DeBruijnGraph &masked) {
    std::vector<std::pair<node_index, std::string>> kmers;
    const DBGSuccinct &dbg_succ = base_dbg(masked);
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&masked);
    for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
        if (!dbg_succ.in_graph(y))
            continue;
        kmers.emplace_back(y, masked.get_node_sequence(y));
        if (canonical) {
            node_index z = canonical->reverse_complement(y);
            if (z != y)
                kmers.emplace_back(z, masked.get_node_sequence(z));
        }
    }
    return kmers;
}

// the stored node of a node of the served graph, by spelling
node_index stored_node(const DeBruijnGraph &masked, node_index node,
                       const std::map<std::string, node_index> &stored) {
    if (!dynamic_cast<const CanonicalDBG*>(&masked))
        return node;
    const std::string kmer = masked.get_node_sequence(node);
    auto it = stored.find(kmer);
    if (it == stored.end())
        it = stored.find(rev_comp(kmer));
    return it == stored.end() ? DeBruijnGraph::npos : it->second;
}

// the graph-walk oracle on the masked twin: the true contexts (or anchors), answer order
std::vector<Ctx> walk_oracle(const DeBruijnGraph &masked, const std::string &text,
                             const Request &request) {
    const size_t k = masked.get_k();
    std::map<std::string, node_index> stored;
    const DBGSuccinct &dbg_succ = base_dbg(masked);
    for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
        if (dbg_succ.in_graph(y))
            stored.emplace(dbg_succ.get_node_sequence(y), y);
    }
    std::vector<Ctx> result;
    for (const auto &[node, kmer] : true_kmers(masked)) {
        for (const auto &[o, q] : orientations(text, request.strands)) {
            const Bases w(q.begin(), q.begin() + std::min(q.size(), k));
            for (uint32_t p : scope_offsets(text.size(), k, request.scope)) {
                if (matches(w, std::string_view(kmer).substr(p, w.size())))
                    result.push_back(Ctx { node, p, o, stored_node(masked, node, stored) });
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

// the record-scan oracle: the contexts spelled from the records' k-mers (both orientations
// unless BASIC), never from the graph
std::set<std::tuple<std::string, uint32_t, Orientation>>
record_oracle(const std::vector<std::string> &records, size_t k, DeBruijnGraph::Mode mode,
              const std::string &text, const Request &request) {
    std::set<std::string> kmers;
    for (const std::string &record : records) {
        for (size_t i = 0; i + k <= record.size(); ++i) {
            const std::string kmer = record.substr(i, k);
            bool indexed = true;
            for (char c : kmer) {
#if _DNA5_GRAPH
                indexed &= std::string("ACGTN").find(c) != std::string::npos;
#else
                indexed &= std::string("ACGT").find(c) != std::string::npos;
#endif
            }
            if (!indexed)
                continue;
            kmers.insert(kmer);
            if (mode != DeBruijnGraph::BASIC)
                kmers.insert(rev_comp(kmer));
        }
    }
    std::set<std::tuple<std::string, uint32_t, Orientation>> result;
    for (const std::string &kmer : kmers) {
        for (const auto &[o, q] : orientations(text, request.strands)) {
            const Bases w(q.begin(), q.begin() + std::min(q.size(), k));
            for (uint32_t p : scope_offsets(text.size(), k, request.scope)) {
                if (matches(w, std::string_view(kmer).substr(p, w.size())))
                    result.emplace(kmer, p, o);
            }
        }
    }
    return result;
}

/**
 * The entry-scan oracle of U at offset |p| of orientation window |w| on the unmasked twin:
 * the entries (W != $) of the stored graph whose spelled k-mer has |w| at |p|, a '$' allowed
 * in the unsearched lead (matches_entry); on a wrapped PRIMARY graph, the union's
 * A + B - P: plus the entries with rc(w) at the mirrored offset, less the real palindromic
 * k-mers holding w at p (found by both, counted once, §4.1). Exactly the engine's U except on
 * an even-k wrapped PRIMARY graph, whose palindrome scans check some candidates: there an
 * upper bound of it.
 */
uint64_t entry_oracle(const DeBruijnGraph &unmasked, const Bases &w, uint32_t p) {
    const DBGSuccinct &stored = base_dbg(unmasked);
    const size_t k = stored.get_k();
    const bool primary = dynamic_cast<const CanonicalDBG*>(&unmasked);
    const Bases wr = oracle_rc(w);
    const uint32_t mirror = static_cast<uint32_t>(k - w.size()) - p;
    uint64_t count = 0;
    for (node_index e = 1; e <= stored.max_index(); ++e) {
        const std::string s = stored.get_node_sequence(e);
        if (s.back() == '$')
            continue;
        count += matches_entry(w, std::string_view(s).substr(p, w.size()));
        if (primary) {
            count += matches_entry(wr, std::string_view(s).substr(mirror, wr.size()));
            if (s.find('$') == std::string::npos && s == rev_comp(s)
                    && matches(w, std::string_view(s).substr(p, w.size())))
                --count;
        }
    }
    return count;
}


// ---------------------------------------------------------------- the checks

/**
 * One (twin, pattern, request): count() on the unmasked twin states true relations, its U
 * per offset and orientation equals the entry-scan oracle (a bound of it on an even-k
 * PRIMARY graph), an empty block is EXACT 0, an exact count equals the masked twin's;
 * enumerate() lists exactly the masked twin's contexts in its order, with EXACT counts equal
 * to the masked twin's; the ALL_OR_COUNT threshold compares U; PARTIAL's prefixes; the same
 * steps in count() and enumerate(), the same ranges as the masked twin. Returns the number of
 * true contexts (anchors for L > k).
 */
size_t check_unmasked(const Twin &twin, const std::vector<std::string> &records,
                      DeBruijnGraph::Mode mode, const std::string &text, PatternKind kind,
                      const Request &request) {
    const DeBruijnGraph &masked = *twin.masked;
    const DeBruijnGraph &unmasked = *twin.unmasked;
    const size_t k = unmasked.get_k();
    const size_t L = text.size();
    const bool primary = mode == DeBruijnGraph::PRIMARY;
    const Pattern pattern = Pattern::parse(kind, text);
    SCOPED_TRACE("pattern " + text + " k " + std::to_string(k) + " mode "
                 + std::to_string(mode) + " scope " + to_string(request.scope) + " strands "
                 + to_string(request.strands));

    EXPECT_FALSE(PatternSearch::support(unmasked).mask_present);
    EXPECT_TRUE(PatternSearch::support(unmasked).supported);

    const std::vector<Ctx> expected = walk_oracle(masked, text, request);
    if (L <= k) {
        std::set<std::tuple<std::string, uint32_t, Orientation>> spelled;
        for (const Ctx &c : expected) {
            spelled.emplace(masked.get_node_sequence(c.node), c.offset, c.orientation);
        }
        EXPECT_EQ(record_oracle(records, k, mode, text, request), spelled);
    }
    std::map<uint32_t, uint64_t> truth_offset;
    std::map<Orientation, uint64_t> truth_orientation;
    for (uint32_t p : scope_offsets(L, k, request.scope)) {
        truth_offset[p] = 0;
    }
    for (const Ctx &c : expected) {
        ++truth_offset[c.offset];
        ++truth_orientation[c.orientation];
    }

    // U per offset and orientation from the entry-scan oracle
    std::map<uint32_t, uint64_t> u_offset;
    std::map<Orientation, uint64_t> u_orientation;
    uint64_t u_total = 0;
    bool offset0_exact = true;
    for (const auto &[o, q] : orientations(text, request.strands)) {
        const Bases w(q.begin(), q.begin() + std::min(q.size(), k));
        offset0_exact &= !unsearched_lead(w);
        for (uint32_t p : scope_offsets(L, k, request.scope)) {
            const uint64_t u = entry_oracle(unmasked, w, p);
            u_offset[p] += u;
            u_orientation[o] += u;
            u_total += u;
        }
    }
    // the bound is exact except where palindrome scans check candidates
    const bool u_exact = !primary || k % 2;
    auto check_upper = [&](const Count &count, uint64_t u, const std::string &what) {
        if (u_exact) {
            EXPECT_EQ(u, upper_of(count)) << what;
        } else {
            EXPECT_GE(u, upper_of(count)) << what;
        }
    };

    // count(): true relations, U, exact zeros, the masked twin's exact counts
    const Result counted = count_of(unmasked, pattern, request);
    const Result counted_masked = count_of(masked, pattern, request);
    EXPECT_FALSE(counted.refusal);
    if (counted.refusal)
        return 0;
    EXPECT_FALSE(counted.stop);
    EXPECT_FALSE(has_note(counted, kNoteThresholdUpperBound));
    EXPECT_EQ(counted_masked.work.ranges_visited, counted.work.ranges_visited);
    const Count &total = L > k ? counted.anchors->total : counted.contexts->total;
    const Count &total_masked = L > k ? counted_masked.anchors->total
                                      : counted_masked.contexts->total;
    EXPECT_TRUE(true_relation(total, expected.size()));
    EXPECT_TRUE(total.relation == Relation::EXACT || total.relation == Relation::BOUNDS);
    check_upper(total, u_total, "total");
    if (!u_total) {
        // an empty block: absence holds
        EXPECT_EQ(Relation::EXACT, total.relation);
        EXPECT_EQ(0u, total.value);
    }
    if (total.relation == Relation::EXACT) {
        EXPECT_EQ(total_masked.value, total.value);
    }
    EXPECT_EQ(Relation::EXACT, total_masked.relation);
    EXPECT_EQ(expected.size(), total_masked.value);
    const auto &by_orientation = L > k ? counted.anchors->by_orientation
                                       : counted.contexts->by_orientation;
    for (const auto &[o, count] : by_orientation) {
        EXPECT_TRUE(true_relation(count, truth_orientation[o])) << orientation_key(o);
        check_upper(count, u_orientation[o], orientation_key(o));
    }
    if (L <= k) {
        for (const auto &[p, count] : counted.contexts->by_offset) {
            EXPECT_TRUE(true_relation(count, truth_offset[p])) << "offset " << p;
            check_upper(count, u_offset[p], "offset " + std::to_string(p));
            // a range at offset 0 (whole nodes spelled, no lead) holds no source dummy
            if (p == 0 && !primary && offset0_exact && request.scope == Scope::ANY_OFFSET) {
                EXPECT_EQ(Relation::EXACT, count.relation) << "offset 0";
            }
        }
    } else if (!primary && offset0_exact) {
        // a long pattern's anchors: the W rule on whole nodes, exact
        EXPECT_EQ(Relation::EXACT, total.relation);
    }

    // enumerate(): the masked twin's list, in its order, the counts made EXACT
    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    if (L > k)
        all.release_anchors = true;
    Result enumerated;
    const std::vector<Ctx> released = enumerate_of(unmasked, pattern, all, &enumerated);
    EXPECT_EQ(counted.work.steps, enumerated.work.steps);
    EXPECT_EQ(expected, released);
    EXPECT_TRUE(enumerated.extraction && enumerated.extraction->complete);
    EXPECT_EQ(expected.size(), enumerated.extraction->returned);
    const Count &listed = L > k ? enumerated.anchors->total : enumerated.contexts->total;
    EXPECT_EQ(Relation::EXACT, listed.relation);
    EXPECT_EQ(expected.size(), listed.value);
    const auto &listed_orientation = L > k ? enumerated.anchors->by_orientation
                                           : enumerated.contexts->by_orientation;
    for (const auto &[o, count] : listed_orientation) {
        EXPECT_EQ(Relation::EXACT, count.relation);
        EXPECT_EQ(truth_orientation[o], count.value) << orientation_key(o);
    }
    if (L <= k) {
        for (const auto &[p, count] : enumerated.contexts->by_offset) {
            EXPECT_EQ(Relation::EXACT, count.relation);
            EXPECT_EQ(truth_offset[p], count.value) << "offset " << p;
        }
        EXPECT_EQ(truth_offset[k - L], enumerated.contexts->suffix.value);
        EXPECT_EQ(Relation::EXACT, enumerated.contexts->suffix.relation);
    }

    // ALL_OR_COUNT compares U with the threshold: above it, nothing released; the note when
    // the lower bound is not above it
    if (upper_of(total) > 0) {
        Request tight = all;
        const uint64_t cap = upper_of(total) - 1;
        (L > k ? tight.max_anchors : tight.max_contexts) = cap;
        Result withheld;
        EXPECT_TRUE(enumerate_of(unmasked, pattern, tight, &withheld).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, withheld.extraction->withheld);
        EXPECT_EQ(total.lower <= cap && total.relation == Relation::BOUNDS,
                  has_note(withheld, kNoteThresholdUpperBound));
        // and at U itself everything is released
        (L > k ? tight.max_anchors : tight.max_contexts) = upper_of(total);
        EXPECT_EQ(expected, enumerate_of(unmasked, pattern, tight, &withheld));
        EXPECT_TRUE(withheld.extraction->complete);
    }

    // PARTIAL: the first contexts in answer order; one release past the end is complete and
    // makes the counts EXACT; a cut one keeps true relations
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    if (L > k)
        partial.release_anchors = true;
    for (uint64_t cap : { uint64_t(0), uint64_t(1), uint64_t(expected.size() / 2),
                          uint64_t(expected.size()), uint64_t(expected.size() + 1),
                          upper_of(total) + 1 }) {
        (L > k ? partial.max_anchors : partial.max_contexts) = cap;
        Result cut;
        const std::vector<Ctx> prefix = enumerate_of(unmasked, pattern, partial, &cut);
        EXPECT_EQ(std::min<uint64_t>(cap, expected.size()), prefix.size()) << cap;
        EXPECT_TRUE(prefix.size() <= expected.size()
                        && std::equal(prefix.begin(), prefix.end(), expected.begin()))
            << cap;
        const Count &t = L > k ? cut.anchors->total : cut.contexts->total;
        EXPECT_TRUE(true_relation(t, expected.size())) << cap;
        if (cap > expected.size()) {
            EXPECT_TRUE(cut.extraction->complete) << cap;
            EXPECT_EQ(Relation::EXACT, t.relation) << cap;
        }
        if (cut.extraction->complete) {
            EXPECT_EQ(Relation::EXACT, t.relation);
            EXPECT_EQ(expected.size(), t.value);
        } else {
            EXPECT_TRUE(cut.extraction->cut) << cap;
            // the lower bound is raised to what was released
            if (t.relation == Relation::BOUNDS) {
                EXPECT_LE(prefix.size(), t.lower) << cap;
            }
        }
    }

    // stopped by max_steps: true relations; PARTIAL lists real contexts only, in answer
    // order, and its counts are at least what it listed; ALL_OR_COUNT withholds
    const uint64_t steps = counted.work.steps;
    for (uint64_t max_steps : { uint64_t(0), uint64_t(1), uint64_t(2), steps / 3, steps / 2,
                                steps ? steps - 1 : 0 }) {
        const Result stopped = count_of(unmasked, pattern, request, max_steps);
        const Count &t = L > k ? stopped.anchors->total : stopped.contexts->total;
        EXPECT_TRUE(true_relation(t, expected.size())) << max_steps;
        const auto &by_o = L > k ? stopped.anchors->by_orientation
                                 : stopped.contexts->by_orientation;
        for (const auto &[o, count] : by_o) {
            EXPECT_TRUE(true_relation(count, truth_orientation[o])) << max_steps;
        }
        if (L <= k) {
            for (const auto &[p, count] : stopped.contexts->by_offset) {
                EXPECT_TRUE(true_relation(count, truth_offset[p])) << max_steps << " " << p;
            }
        }
        Request partial_all = partial;
        (L > k ? partial_all.max_anchors : partial_all.max_contexts) = 1'000'000'000;
        Result cut;
        const std::vector<Ctx> listed_stopped
            = enumerate_of(unmasked, pattern, partial_all, &cut, max_steps);
        EXPECT_TRUE(std::includes(expected.begin(), expected.end(), listed_stopped.begin(),
                                  listed_stopped.end()))
            << max_steps;
        const Count &ct = L > k ? cut.anchors->total : cut.contexts->total;
        EXPECT_TRUE(true_relation(ct, expected.size())) << max_steps;
        if (ct.relation != Relation::UNKNOWN) {
            EXPECT_LE(listed_stopped.size(), ct.value) << max_steps;
        }
        Result withheld;
        EXPECT_TRUE(enumerate_of(unmasked, pattern, all, &withheld, max_steps).empty()
                        || withheld.extraction->complete)
            << max_steps;
    }

    // L > k: the extension gives the masked twin's paths, with the same counts
    if (L > k) {
        Request extend = request;
        extend.extend_paths = true;
        extend.mode = Mode::ALL_OR_COUNT;
        std::vector<std::tuple<node_index, Orientation, std::string, std::vector<node_index>>>
            paths[2];
        Result results[2];
        for (int m = 0; m < 2; ++m) {
            Budget budget(kManySteps, Deadline::unbounded());
            results[m] = PatternSearch(m ? masked : unmasked)
                .enumerate(pattern, extend, budget, [&](const Context &c) {
                    paths[m].emplace_back(c.node, c.orientation, c.sequence, c.path);
                });
        }
        EXPECT_EQ(paths[1], paths[0]);
        EXPECT_EQ(results[1].anchors->paths.relation, results[0].anchors->paths.relation);
        EXPECT_EQ(results[1].anchors->paths.value, results[0].anchors->paths.value);
        EXPECT_EQ(Relation::EXACT, results[0].anchors->total.relation);
        EXPECT_EQ(expected.size(), results[0].anchors->total.value);
        EXPECT_TRUE(results[0].extraction->complete);
        EXPECT_EQ(results[1].anchors->candidates_examined,
                  results[0].anchors->candidates_examined);
        const Result counted_paths = count_of(unmasked, pattern, extend);
        EXPECT_EQ(results[1].anchors->paths.value, counted_paths.anchors->paths.value);
        EXPECT_EQ(Relation::EXACT, counted_paths.anchors->paths.relation);
    }
    return expected.size();
}


// ---------------------------------------------------------------- graphs and primitives

TEST(PatternUnmasked, GraphWithoutMaskIsServed) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                       DeBruijnGraph::PRIMARY }) {
        Twin twin = build_twin(5, { "ACGTACGGTTAC", "TTGCAGG" }, mode);
        GraphSupport s = PatternSearch::support(*twin.unmasked);
        EXPECT_TRUE(s.supported) << s.reason;
        EXPECT_TRUE(s.reason.empty());
        EXPECT_FALSE(s.mask_present);
        EXPECT_NO_THROW(PatternSearch engine(*twin.unmasked));
        // the masked twin: as before
        s = PatternSearch::support(*twin.masked);
        EXPECT_TRUE(s.supported);
        EXPECT_TRUE(s.mask_present);
    }
}

TEST(PatternUnmasked, NodeHasSentinelMatchesTheSpelling) {
    // every edge of graphs of every mode, both builders, with and without an index of suffix
    // ranges of several lengths: node_has_sentinel equals "the spelled node holds $", and on
    // the entries with W != $ BOSS's traversal of the dummy tree agrees
    std::mt19937 rng(41);
    for (int t = 0; t < 24; ++t) {
        const size_t k = 2 + rng() % 9;
        std::vector<std::string> records(1 + rng() % 5);
        for (std::string &record : records) {
            for (size_t i = 0, n = rng() % 20; i < n; ++i) {
                record.push_back("ACGT"[rng() % 4]);
            }
        }
        auto graph = build(k, records, static_cast<DeBruijnGraph::Mode>(rng() % 2), rng() % 2);
        DBGSuccinct &dbg_succ = const_cast<DBGSuccinct&>(base_dbg(*graph));
        BOSS &boss = dbg_succ.get_boss();
        for (size_t suffix : { size_t(0), size_t(1), size_t(2), boss.get_k() }) {
            if (suffix > boss.get_k())
                continue;
            boss.index_suffix_ranges(suffix);
            sdsl::bit_vector dummies(boss.get_W().size(), false);
            boss.mark_source_dummy_edges(&dummies, 1);
            for (node_index e = 1; e <= boss.num_edges(); ++e) {
                const std::vector<BOSS::TAlphabet> node = boss.get_node_seq(e);
                const bool spelled = std::count(node.begin(), node.end(), 0) > 0;
                EXPECT_EQ(spelled, boss.node_has_sentinel(e))
                    << "k " << k << " edge " << e << " suffix " << suffix;
                if (boss.get_W(e) % boss.alph_size) {
                    EXPECT_EQ(bool(dummies[e]), spelled) << "k " << k << " edge " << e;
                }
            }
        }
    }
}

TEST(PatternUnmasked, PrimitivesWithoutTheMask) {
    // the unmasked primitives against a read of W, on random ranges of whole node groups
    std::mt19937 rng(7);
    for (int t = 0; t < 12; ++t) {
        const size_t k = 3 + rng() % 6;
        std::vector<std::string> records(1 + rng() % 4);
        for (std::string &record : records) {
            for (size_t i = 0, n = k + rng() % 15; i < n; ++i) {
                record.push_back("ACGT"[rng() % 4]);
            }
        }
        Twin twin = build_twin(k, records, DeBruijnGraph::BASIC, rng() % 2);
        const DBGSuccinct &g = base_dbg(*twin.unmasked);
        const BOSS &boss = g.get_boss();
        const node_index n = g.max_index();
        for (node_index first = 1; first <= n; ++first) {
            for (node_index last = first; last <= n; ++last) {
                uint64_t non_sink = 0;
                std::vector<uint64_t> with(boss.alph_size, 0);
                node_index next = 0;
                for (node_index e = first; e <= last; ++e) {
                    const auto w = boss.get_W(e) % boss.alph_size;
                    if (w) {
                        ++non_sink;
                        ++with[w];
                        if (!next)
                            next = e;
                    }
                }
                ASSERT_EQ(non_sink, g.count_non_sink_edges_in_range(first, last));
                ASSERT_EQ(next, g.next_non_sink_edge(first, last));
                for (BOSS::TAlphabet c = 1; c < boss.alph_size; ++c) {
                    ASSERT_EQ(with[c], g.count_edges_with_symbol(first, last, c));
                }
            }
        }
        // the masked twin: valid edges are the non-sink ones less the source dummies
        const DBGSuccinct &m = base_dbg(*twin.masked);
        uint64_t real = 0;
        for (node_index e = 1; e <= n; ++e) {
            real += boss.get_W(e) % boss.alph_size && !boss.node_has_sentinel(e);
        }
        EXPECT_EQ(m.count_valid_edges_in_range(1, n), real);
    }
}


// ---------------------------------------------------------------- named cases

TEST(PatternUnmasked, SourceDummiesCountedNeverReleased) {
    // k = 5, one record ACGTTGCA: its source dummies are $$$$A, $$$AC, $$ACG, $ACGT (and the
    // main dummy). ACG forward: the real k-mer ACGTT at offset 0; the dummies $$ACG (offset 2)
    // and $ACGT (offset 1) are entries too: U = 3, the true count 1
    for (bool batch : { false, true }) {
        Twin twin = build_twin(5, { "ACGTTGCA" }, DeBruijnGraph::BASIC, batch);
        const Pattern acg = Pattern::parse(PatternKind::DNA, "ACG");
        Request forward = make_request(Scope::ANY_OFFSET, Strands::FORWARD);

        Result r = count_of(*twin.unmasked, acg, forward);
        EXPECT_EQ(Relation::BOUNDS, r.contexts->total.relation);
        EXPECT_EQ(1u, r.contexts->total.lower);
        EXPECT_EQ(3u, r.contexts->total.upper);
        EXPECT_EQ(Relation::EXACT, r.contexts->by_offset.at(0).relation);
        EXPECT_EQ(1u, r.contexts->by_offset.at(0).value);
        for (uint32_t p : { 1u, 2u }) {
            EXPECT_EQ(Relation::BOUNDS, r.contexts->by_offset.at(p).relation) << p;
            EXPECT_EQ(0u, r.contexts->by_offset.at(p).lower) << p;
            EXPECT_EQ(1u, r.contexts->by_offset.at(p).upper) << p;
        }
        // the masked twin: exact 1, and the same number once listed without the mask
        EXPECT_EQ(Relation::EXACT, count_of(*twin.masked, acg, forward).contexts->total.relation);

        Request all = forward;
        all.mode = Mode::ALL_OR_COUNT;
        Result e;
        std::vector<Ctx> listed = enumerate_of(*twin.unmasked, acg, all, &e);
        ASSERT_EQ(1u, listed.size());
        EXPECT_EQ("ACGTT", twin.unmasked->get_node_sequence(listed[0].node));
        EXPECT_EQ(0u, listed[0].offset);
        EXPECT_EQ(Relation::EXACT, e.contexts->total.relation);
        EXPECT_EQ(1u, e.contexts->total.value);
        EXPECT_EQ(Relation::EXACT, e.contexts->by_offset.at(2).relation);
        EXPECT_EQ(0u, e.contexts->by_offset.at(2).value);
        EXPECT_TRUE(e.extraction->complete);
        EXPECT_FALSE(has_note(e, kNoteThresholdUpperBound));

        // threshold 2: U = 3 is above it, the lower bound 1 is not: withheld, said so
        all.max_contexts = 2;
        EXPECT_TRUE(enumerate_of(*twin.unmasked, acg, all, &e).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, e.extraction->withheld);
        EXPECT_EQ(Relation::BOUNDS, e.contexts->total.relation);
        EXPECT_TRUE(has_note(e, kNoteThresholdUpperBound));
        // threshold 0: both above it: withheld, no note (the count is above it for certain)
        all.max_contexts = 0;
        EXPECT_TRUE(enumerate_of(*twin.unmasked, acg, all, &e).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, e.extraction->withheld);
        EXPECT_FALSE(has_note(e, kNoteThresholdUpperBound));
        // the masked twin releases at 2
        all.max_contexts = 2;
        EXPECT_EQ(1u, enumerate_of(*twin.masked, acg, all, &e).size());

        // stop_at_threshold at 1 compares U: stopped, AT_LEAST, said so; the masked twin
        // (exact 1, not above 1) is not stopped
        Request stop = forward;
        stop.stop_at_threshold = true;
        stop.max_contexts = 1;
        r = count_of(*twin.unmasked, acg, stop);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopReason::MAX_CONTEXTS, r.stop->reason);
        EXPECT_TRUE(true_relation(r.contexts->total, 1));
        EXPECT_TRUE(has_note(r, kNoteThresholdUpperBound));
        EXPECT_FALSE(count_of(*twin.masked, acg, stop).stop);

        // PARTIAL past the end: complete and EXACT; at 1, the one context
        Request partial = forward;
        partial.mode = Mode::PARTIAL;
        partial.max_contexts = 10;
        EXPECT_EQ(1u, enumerate_of(*twin.unmasked, acg, partial, &e).size());
        EXPECT_TRUE(e.extraction->complete);
        EXPECT_EQ(Relation::EXACT, e.contexts->total.relation);
        partial.max_contexts = 1;
        EXPECT_EQ(1u, enumerate_of(*twin.unmasked, acg, partial, &e).size());
        EXPECT_TRUE(true_relation(e.contexts->total, 1));
        EXPECT_EQ(1u, e.contexts->total.value);

        // suffix scope: the W rule's candidates at offset 2, ACGTT's suffix is TTG..: only
        // the dummy $$ACG ends with ACG: U = 1, true 0
        Request suffix = make_request(Scope::SUFFIX, Strands::FORWARD);
        r = count_of(*twin.unmasked, acg, suffix);
        EXPECT_EQ(Relation::BOUNDS, r.contexts->total.relation);
        EXPECT_EQ(0u, r.contexts->total.lower);
        EXPECT_EQ(1u, r.contexts->total.upper);
        suffix.mode = Mode::ALL_OR_COUNT;
        EXPECT_TRUE(enumerate_of(*twin.unmasked, acg, suffix, &e).empty());
        EXPECT_TRUE(e.extraction->complete);
        EXPECT_EQ(Relation::EXACT, e.contexts->total.relation);
        EXPECT_EQ(0u, e.contexts->total.value);

        check_unmasked(twin, { "ACGTTGCA" }, DeBruijnGraph::BASIC, "ACG", PatternKind::DNA,
                       forward);
        check_unmasked(twin, { "ACGTTGCA" }, DeBruijnGraph::BASIC, "ACG", PatternKind::DNA,
                       make_request(Scope::SUFFIX));
    }
}

TEST(PatternUnmasked, LeadingNRunCountsDummiesWithTheirSentinelUnderTheN) {
    // SPEC §7.4 (review of the mask round, wire finding): U is the masked count plus the
    // source dummies holding the window after their $ run, PLUS, for a window with a leading
    // N run (skipped on $ACGT, §7.8), the dummies whose $ run ends inside the skipped run:
    // those hold '$' under an N, so they do not hold the pattern. k = 5, one record CAGTA:
    // the dummies $CAGT, $$CAG, $$$CA, $$$$C. NC forward: no real context, no dummy holding
    // NC on bases, but each dummy has its last '$' under the N at one offset: U = 1 at every
    // offset, BOUNDS [0, 1], never the EXACT 0 that "masked + dummies holding NC" would give.
    // GN reverse: the same window (the trailing run leads the reverse complement). NN: a
    // window of N only, whose last position is searched: U = 2, 3, 4, 5 (masked 1 each,
    // 0..3 dummies holding NN on bases, one more with '$' under the first N)
#if _DNA5_GRAPH
    GTEST_SKIP() << "$ACGTN: no N run is skipped";
#endif
    Twin twin = build_twin(5, { "CAGTA" }, DeBruijnGraph::BASIC);
    const std::map<std::string, std::vector<uint64_t>> upper {
        { "NC", { 1, 1, 1, 1 } }, { "GN", { 1, 1, 1, 1 } }, { "NN", { 2, 3, 4, 5 } } };
    const std::map<std::string, std::vector<uint64_t>> masked_count {
        { "NC", { 0, 0, 0, 0 } }, { "GN", { 0, 0, 0, 0 } }, { "NN", { 1, 1, 1, 1 } } };
    for (const auto &[text, u] : upper) {
        SCOPED_TRACE(text);
        const Pattern p = Pattern::parse(PatternKind::IUPAC, text);
        const Request request = make_request(Scope::ANY_OFFSET,
                                             text == "GN" ? Strands::REVERSE : Strands::FORWARD);
        const Result r = count_of(*twin.unmasked, p, request);
        const Result m = count_of(*twin.masked, p, request);
        ASSERT_TRUE(r.contexts && m.contexts);
        for (uint32_t offset = 0; offset < 4; ++offset) {
            const Count &c = r.contexts->by_offset.at(offset);
            EXPECT_EQ(Relation::BOUNDS, c.relation) << offset;
            EXPECT_EQ(u[offset], c.upper) << offset;
            EXPECT_EQ(Relation::EXACT, m.contexts->by_offset.at(offset).relation) << offset;
            EXPECT_EQ(masked_count.at(text)[offset], m.contexts->by_offset.at(offset).value)
                << offset;
            EXPECT_TRUE(true_relation(c, masked_count.at(text)[offset])) << offset;
        }
    }

    // the decomposition against the spelled entries, on BASIC and CANONICAL graphs (on an
    // even-k wrapped PRIMARY graph the palindrome scans settle some candidates): U = the
    // masked count + the dummies holding the window on bases + the dummies holding a '$' in
    // the skipped lead; the last term is not empty for a leading N run, empty without one
    const std::vector<std::string> records { "CAGTACCA", "GGCATTG", "TCCAGT" };
    size_t under_the_n = 0;
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        Twin t = build_twin(5, records, mode);
        const DBGSuccinct &stored = base_dbg(*t.unmasked);
        for (const char *text : { "NC", "NNC", "NCA", "NNNG", "NN", "CA", "GNC", "CAN" }) {
            SCOPED_TRACE(std::string(text) + " mode " + std::to_string(mode));
            const Pattern p = Pattern::parse(PatternKind::IUPAC, text);
            const Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH);
            const Result r = count_of(*t.unmasked, p, request);
            const Result m = count_of(*t.masked, p, request);
            ASSERT_TRUE(r.contexts && m.contexts);
            for (const auto &[o, q] : orientations(text, Strands::BOTH)) {
                const size_t lead = unsearched_lead(q);
                uint64_t total_masked = 0, on_bases = 0, sentinel_in_lead = 0;
                for (uint32_t offset : scope_offsets(q.size(), 5, Scope::ANY_OFFSET)) {
                    for (node_index e = 1; e <= stored.max_index(); ++e) {
                        const std::string s = stored.get_node_sequence(e);
                        if (s.back() == '$')
                            continue;
                        const std::string_view piece = std::string_view(s).substr(offset,
                                                                                  q.size());
                        const bool dummy = s.find('$') != std::string::npos;
                        if (matches(q, piece)) {
                            ++(dummy ? on_bases : total_masked);
                        } else if (matches_entry(q, piece)) {
                            // a '$' stands in the skipped lead, and only there
                            EXPECT_TRUE(dummy);
                            EXPECT_LT(piece.rfind('$'), lead);
                            ++sentinel_in_lead;
                        }
                    }
                }
                const Count &c = r.contexts->by_orientation.at(o);
                const Count &cm = m.contexts->by_orientation.at(o);
                EXPECT_EQ(Relation::EXACT, cm.relation);
                EXPECT_EQ(total_masked, cm.value) << orientation_key(o);
                EXPECT_EQ(total_masked + on_bases + sentinel_in_lead, upper_of(c))
                    << orientation_key(o);
                if (!lead) {
                    EXPECT_EQ(0u, sentinel_in_lead) << orientation_key(o);
                }
                under_the_n += sentinel_in_lead;
            }
        }
    }
    // not vacuous: the term the SPEC's first wording left out
    EXPECT_LT(0u, under_the_n);
}

TEST(PatternUnmasked, EmptyBlockIsExactZero) {
    // a pattern no entry holds: U = 0, EXACT 0 (absence holds) without any release
    Twin twin = build_twin(6, { "ACGTTGCAAC", "GGGTTTAAAC" }, DeBruijnGraph::BASIC);
    for (const char *text : { "CCCC", "CACACACAC", "TATATATAT" }) {
        const Pattern p = Pattern::parse(PatternKind::DNA, text);
        Result r = count_of(*twin.unmasked, p, make_request());
        const Count &total = r.contexts ? r.contexts->total : r.anchors->total;
        EXPECT_EQ(Relation::EXACT, total.relation) << text;
        EXPECT_EQ(0u, total.value) << text;
    }
}

TEST(PatternUnmasked, LongPatternsAnchorsAndPaths) {
    // L > k: the anchor window's W rule is on whole nodes (no lead): anchors EXACT without
    // the mask; the extension then gives the masked twin's paths. An anchor window starting
    // with N (IUPAC, $ACGT: a lead) is unchecked: anchors BOUNDS, listed (dummies dropped)
    // and extended all the same
    const std::vector<std::string> records { "ACGTTGCAACGT", "ACGTAGGA", "TTACGTTGCAT" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                       DeBruijnGraph::PRIMARY }) {
        Twin twin = build_twin(5, records, mode);
        for (const char *text : { "ACGTTGCA", "NCGTTGCA", "NNGTTGCAA", "ACGTNGCA" }) {
            const Pattern p = Pattern::parse(PatternKind::IUPAC, text);
            SCOPED_TRACE(std::string(text) + " mode " + std::to_string(mode));
            Request request = make_request();
            request.extend_paths = true;
            Result unmasked = count_of(*twin.unmasked, p, request);
            Result masked = count_of(*twin.masked, p, request);
            EXPECT_EQ(masked.anchors->paths.relation, unmasked.anchors->paths.relation);
            EXPECT_EQ(masked.anchors->paths.value, unmasked.anchors->paths.value);
            EXPECT_EQ(masked.anchors->extension == Extension::NO_ANCHORS
                          ? Extension::NO_ANCHORS : Extension::COMPLETED,
                      unmasked.anchors->extension);
            EXPECT_EQ(Relation::EXACT, unmasked.anchors->total.relation);
            EXPECT_EQ(masked.anchors->total.value, unmasked.anchors->total.value);
            EXPECT_EQ(masked.anchors->candidates_examined,
                      unmasked.anchors->candidates_examined);
            // the paths released: the masked twin's
            request.mode = Mode::ALL_OR_COUNT;
            std::vector<std::pair<std::string, std::vector<node_index>>> paths[2];
            for (int m = 0; m < 2; ++m) {
                Budget budget(kManySteps, Deadline::unbounded());
                PatternSearch(m ? *twin.masked : *twin.unmasked)
                    .enumerate(p, request, budget, [&](const Context &c) {
                        paths[m].emplace_back(c.sequence, c.path);
                    });
            }
            EXPECT_EQ(paths[1], paths[0]);
        }
        // anchors whose upper bound is above max_anchors while their lower bound is not:
        // not admitted, said so
        const Pattern lead = Pattern::parse(PatternKind::IUPAC, "NCGTTGCA");
        Request request = make_request();
        request.extend_paths = true;
        Result r = count_of(*twin.unmasked, lead, request);
        const Count total = count_of(*twin.unmasked, lead, make_request()).anchors->total;
        if (total.relation == Relation::BOUNDS && total.upper > total.lower) {
            request.max_anchors = total.upper - 1;
            r = count_of(*twin.unmasked, lead, request);
            EXPECT_EQ(Extension::NOT_ADMITTED, r.anchors->extension);
            EXPECT_EQ(Relation::UNKNOWN, r.anchors->paths.relation);
            EXPECT_EQ(total.lower <= request.max_anchors,
                      has_note(r, kNoteThresholdUpperBound));
        }
    }
}


// ---------------------------------------------------------------- random cases

TEST(PatternUnmasked, RandomPatternsAgainstOracles) {
    // fixed seeds: graphs of every mode and both builders, k 3 to 12, short records (many
    // sources, so many source dummies), patterns from the records' starts (where the dummies
    // match), IUPAC with leading N runs, short and long, every strand choice and scope
    size_t cases = 0;
    size_t with_dummies = 0;
    size_t nonempty = 0;
    for (uint32_t seed = 1; seed <= 60; ++seed) {
        std::mt19937 rng(9000 + seed);
        const size_t k = 3 + rng() % 10;
        const auto mode = static_cast<DeBruijnGraph::Mode>(rng() % 3);
        const bool batch = rng() % 2;
        std::vector<std::string> records(1 + rng() % 6);
        for (std::string &record : records) {
            for (size_t i = 0, n = k + rng() % 12; i < n; ++i) {
                record.push_back("ACGT"[rng() % 4]);
            }
        }
        const Twin twin = build_twin(k, records, mode, batch);
        for (int t = 0; t < 6; ++t) {
            const std::string &record = records[rng() % records.size()];
            const size_t L = 1 + rng() % std::min(record.size(), k + 3);
            // mostly a record's start (the source dummies spell it)
            const size_t start = rng() % 3 ? 0 : rng() % (record.size() - L + 1);
            std::string text = record.substr(start, L);
            PatternKind kind = PatternKind::DNA;
            if (rng() % 2) {
                kind = PatternKind::IUPAC;
                for (char &c : text) {
                    if (rng() % 4 == 0)
                        c = "RYSWKMBDHVN"[rng() % 11];
                }
                // a leading N run, which the engine does not search on $ACGT
                for (size_t i = 0, n = rng() % 3; i < n && i < text.size(); ++i) {
                    text[i] = 'N';
                }
            }
            Scope scope = Scope::ANY_OFFSET;
            if (text.size() <= k && mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0)
                scope = Scope::SUFFIX;
            const auto strands = static_cast<Strands>(rng() % 3);
            SCOPED_TRACE("seed " + std::to_string(seed) + " batch " + std::to_string(batch));
            const Request request = make_request(scope, strands);
            const size_t found = check_unmasked(twin, records, mode, text, kind, request);
            ++cases;
            nonempty += found > 0;
            const Result r = count_of(*twin.unmasked, Pattern::parse(kind, text), request);
            const Count &total = r.contexts ? r.contexts->total : r.anchors->total;
            with_dummies += upper_of(total) > found;
        }
    }
    EXPECT_EQ(360u, cases);
    // not vacuous: many patterns hit, and many counts hold source dummies
    EXPECT_LT(200u, nonempty);
    EXPECT_LT(60u, with_dummies);
    std::cerr << "unmasked cases: " << cases << ", with hits: " << nonempty
              << ", with source dummies in U: " << with_dummies << std::endl;
}

TEST(PatternUnmasked, PeptidesWithoutTheMask) {
    // peptides (the codon automaton) give the masked twin's contexts and paths
    const std::vector<std::string> records { "ATGTAAGGCTGGTGAATGCCC", "ATGAAATAGTGGCC" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                       DeBruijnGraph::PRIMARY }) {
        for (size_t k : { size_t(4), size_t(6), size_t(9) }) {
            Twin twin = build_twin(k, records, mode);
            for (const char *peptide : { "MK", "MP", "GW", "MX", "WG", "MKX" }) {
                for (int table : { 1, 2, 11, 27 }) {
                    const Pattern p = Pattern::parse(PatternKind::PROTEIN, peptide,
                                                     GeneticCode::get(table));
                    Request request = make_request();
                    request.extend_paths = p.length() > k;
                    request.mode = Mode::ALL_OR_COUNT;
                    Result results[2];
                    std::vector<std::tuple<node_index, uint32_t, Orientation, std::string>>
                        listed[2];
                    for (int m = 0; m < 2; ++m) {
                        Budget budget(kManySteps, Deadline::unbounded());
                        results[m] = PatternSearch(m ? *twin.masked : *twin.unmasked)
                            .enumerate(p, request, budget, [&](const Context &c) {
                                listed[m].emplace_back(c.node, c.offset, c.orientation,
                                                       c.sequence);
                            });
                    }
                    SCOPED_TRACE(std::string(peptide) + " table " + std::to_string(table)
                                 + " k " + std::to_string(k));
                    EXPECT_EQ(listed[1], listed[0]);
                    EXPECT_TRUE(results[0].extraction->complete);
                    const Count &u = results[0].contexts ? results[0].contexts->total
                                                         : results[0].anchors->paths;
                    const Count &m = results[1].contexts ? results[1].contexts->total
                                                         : results[1].anchors->paths;
                    EXPECT_EQ(Relation::EXACT, u.relation);
                    EXPECT_EQ(m.value, u.value);
                }
            }
        }
    }
}


// ---------------------------------------------------------------- the real fraction f

// the draws re-implemented here (the generator, the modulo, the redraw of a W = $ entry) with
// the spelled test of a dummy: what sample_real_fraction must give, count for count
RealFraction oracle_draws(const DBGSuccinct &graph, uint64_t samples) {
    const BOSS &boss = graph.get_boss();
    RealFraction f;
    f.edges = boss.num_edges();
    for (node_index e = 1; e <= f.edges; ++e) {
        f.sentinel_edges += !(boss.get_W(e) % boss.alph_size);
    }
    f.seed = f.edges;
    if (f.sentinel_edges == f.edges)
        return f;
    std::mt19937_64 rng(f.edges);
    while (f.samples < samples) {
        const uint64_t e = 1 + rng() % f.edges;
        if (!(boss.get_W(e) % boss.alph_size))
            continue;
        ++f.samples;
        const auto node = boss.get_node_seq(e);
        f.real += !std::count(node.begin(), node.end(), 0);
    }
    return f;
}

// f over every entry, spelled
std::pair<uint64_t, uint64_t> spelled_fraction(const DBGSuccinct &graph) {
    const BOSS &boss = graph.get_boss();
    uint64_t entries = 0;
    uint64_t real = 0;
    for (node_index e = 1; e <= boss.num_edges(); ++e) {
        if (!(boss.get_W(e) % boss.alph_size))
            continue;
        ++entries;
        const auto node = boss.get_node_seq(e);
        real += !std::count(node.begin(), node.end(), 0);
    }
    return { real, entries };
}

TEST(PatternUnmasked, RealFractionOnSmallGraphs) {
    // f on random small graphs: exact_real_fraction equals the spelled count; the sample is
    // the re-implemented draws, count for count, deterministic, seeded with the edges, and
    // its 95% interval holds the exact value (fixed seeds: a deterministic check)
    std::mt19937 rng(23);
    size_t inside = 0;
    size_t graphs = 0;
    for (int t = 0; t < 30; ++t) {
        const size_t k = 3 + rng() % 20;
        std::vector<std::string> records(1 + rng() % 30);
        for (std::string &record : records) {
            for (size_t i = 0, n = k + rng() % 40; i < n; ++i) {
                record.push_back("ACGT"[rng() % 4]);
            }
        }
        auto graph = build(k, records, static_cast<DeBruijnGraph::Mode>(rng() % 2), rng() % 2);
        DBGSuccinct &dbg_succ = const_cast<DBGSuccinct&>(base_dbg(*graph));
        dbg_succ.reset_mask();

        const auto [real, entries] = spelled_fraction(dbg_succ);
        const RealFraction exact = exact_real_fraction(dbg_succ);
        EXPECT_TRUE(exact.exact);
        EXPECT_EQ(entries, exact.samples);
        EXPECT_EQ(real, exact.real);
        EXPECT_EQ(dbg_succ.max_index() - entries, exact.sentinel_edges);
        EXPECT_DOUBLE_EQ(double(real) / entries, exact.value);
        EXPECT_EQ(exact.value, exact.lower);
        EXPECT_EQ(exact.value, exact.upper);

        for (uint64_t samples : { uint64_t(1), uint64_t(57), kRealFractionSamples }) {
            const RealFraction f = sample_real_fraction(dbg_succ, samples);
            const RealFraction again = sample_real_fraction(dbg_succ, samples);
            const RealFraction oracle = oracle_draws(dbg_succ, samples);
            EXPECT_FALSE(f.exact);
            EXPECT_EQ(samples, f.samples);
            EXPECT_EQ(oracle.real, f.real) << samples;
            EXPECT_EQ(oracle.sentinel_edges, f.sentinel_edges);
            EXPECT_EQ(dbg_succ.max_index(), f.edges);
            EXPECT_EQ(dbg_succ.max_index(), f.seed);
            EXPECT_EQ(f.real, again.real);
            EXPECT_EQ(f.value, again.value);
            EXPECT_EQ(f.lower, again.lower);
            EXPECT_EQ(f.upper, again.upper);
            EXPECT_DOUBLE_EQ(double(f.real) / samples, f.value);
            EXPECT_LE(f.lower, f.value);
            EXPECT_GE(f.upper, f.value);
            EXPECT_LE(0, f.lower);
            EXPECT_GE(1, f.upper);
            if (samples == kRealFractionSamples) {
                ++graphs;
                inside += f.lower <= exact.value && exact.value <= f.upper;
            }
        }
        // the default is 10,000
        EXPECT_EQ(kRealFractionSamples, sample_real_fraction(dbg_succ).samples);
        EXPECT_EQ(10'000u, kRealFractionSamples);
    }
    // a 95% interval: of 30 graphs, nearly all (these fixed seeds: all but at most 3)
    std::cerr << "exact f inside the 95% interval: " << inside << " of " << graphs
              << std::endl;
    EXPECT_LE(graphs - 3, inside);
}

TEST(PatternUnmasked, RealFractionWilsonInterval) {
    // the interval of a known sample: 9,000 of 10,000 real (z = 1.959964)
    const double n = 10'000, p = 0.9, z = 1.959963984540054;
    const double centre = (p + z * z / (2 * n)) / (1 + z * z / n);
    const double half = z / (1 + z * z / n) * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n));
    EXPECT_NEAR(0.894, centre - half, 0.001);
    EXPECT_NEAR(0.906, centre + half, 0.001);
    // a graph without a single entry a pattern can count (only the main dummy): value 1,
    // no samples, the interval [0, 1]
    DBGSuccinct empty(5);
    const RealFraction f = sample_real_fraction(empty);
    EXPECT_EQ(0u, f.samples);
    EXPECT_EQ(1, f.value);
    EXPECT_EQ(0, f.lower);
    EXPECT_EQ(1, f.upper);
}

// build/mini_refseq's graph, when present beside the source tree
std::string mini_refseq_graph() {
    std::string here = __FILE__;
    const std::string marker = "/tests/graph/alignment/";
    const size_t at = here.rfind(marker);
    if (at == std::string::npos)
        return "";
    const std::string path = here.substr(0, at) + "/build/mini_refseq/graph_k31.dbg";
    return std::ifstream(path).good() ? path : "";
}

TEST(PatternUnmasked, RealFractionOnMiniRefseq) {
    // build/mini_refseq (k = 31, no .edgemask): 8,335,760 entries, 375 source and 12 sink
    // dummies by `stats --count-dummy`. The sample's 95% interval holds the exact f
    // (deterministic: the seed is the graph's number of edges), and its cost is measured
    const std::string path = mini_refseq_graph();
    if (path.empty())
        GTEST_SKIP() << "build/mini_refseq/graph_k31.dbg not found";
    DBGSuccinct graph(2);
    ASSERT_TRUE(graph.load(path));
    graph.reset_mask();
    const BOSS &boss = graph.get_boss();
    EXPECT_EQ(8'335'760u, boss.num_edges());

    const RealFraction exact = exact_real_fraction(graph);
    sdsl::bit_vector dummies(boss.get_W().size(), false);
    const uint64_t source = boss.mark_source_dummy_edges(&dummies, 1);
    const uint64_t sink = boss.mark_sink_dummy_edges();
    EXPECT_EQ(375u, source);
    EXPECT_EQ(12u, sink);
    std::cerr << "mini_refseq: edges " << boss.num_edges() << ", W = $ " << exact.sentinel_edges
              << ", source dummies " << source << " (main dummy edge marked: "
              << bool(dummies[1]) << "), sink dummies " << sink << ", exact f "
              << std::setprecision(9) << exact.value << " (" << exact.real << " / "
              << exact.samples << "), indexed suffix length "
              << boss.get_indexed_suffix_length() << std::endl;
    // the entries a pattern can count are the edges with W != $: all but the sinks and the
    // main dummy edge, which `stats --count-dummy` counts among the source dummies; its real
    // edges (edges - source - sink) are the real entries
    EXPECT_TRUE(dummies[1]);
    EXPECT_EQ(sink + 1, exact.sentinel_edges);
    EXPECT_EQ(boss.num_edges() - exact.sentinel_edges, exact.samples);
    EXPECT_EQ(boss.num_edges() - source - sink, exact.real);
    EXPECT_EQ(8'335'373u, exact.real);

    const auto start = std::chrono::steady_clock::now();
    const RealFraction f = sample_real_fraction(graph);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::cerr << "mini_refseq sampled f " << f.value << " [" << f.lower << ", " << f.upper
              << "] from " << f.samples << " entries (" << f.real << " real) in "
              << seconds * 1000 << " ms" << std::endl;
    EXPECT_EQ(kRealFractionSamples, f.samples);
    EXPECT_LE(f.lower, exact.value);
    EXPECT_GE(f.upper, exact.value);
    const RealFraction oracle = oracle_draws(graph, kRealFractionSamples);
    EXPECT_EQ(oracle.real, f.real);
    EXPECT_EQ(oracle.sentinel_edges, f.sentinel_edges);
}

// the sampler's cost on a random graph of |records| x |length| bases at k = 31 (few sources:
// f near 1): 10,000 walks of at most 30 bwd steps, whatever the graph's size but for the
// caches it misses
void measure_real_fraction(size_t records, size_t length) {
    std::mt19937 rng(3);
    std::vector<std::string> sequences(records);
    for (std::string &record : sequences) {
        record.resize(length);
        for (char &c : record) {
            c = "ACGT"[rng() % 4];
        }
    }
    auto graph = build_graph_batch<DBGSuccinct>(31, sequences, DeBruijnGraph::BASIC);
    DBGSuccinct &dbg_succ = const_cast<DBGSuccinct&>(base_dbg(*graph));
    dbg_succ.reset_mask();
    const auto start = std::chrono::steady_clock::now();
    const RealFraction f = sample_real_fraction(dbg_succ);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    const RealFraction exact = exact_real_fraction(dbg_succ);
    std::cerr << "synthetic graph: " << dbg_succ.max_index() << " edges, sampled f "
              << std::setprecision(9) << f.value << " [" << f.lower << ", " << f.upper
              << "] in " << seconds * 1000 << " ms; exact f " << exact.value << std::endl;
    EXPECT_LE(f.lower, exact.value);
    EXPECT_GE(f.upper, exact.value);
    // milliseconds (generous, for slow CI machines)
    EXPECT_GT(2.0, seconds);
}

TEST(PatternUnmasked, RealFractionCost) {
    measure_real_fraction(40, 25'000);
}

// by hand (--gtest_also_run_disabled_tests): about 30 million k-mers
TEST(PatternUnmasked, DISABLED_RealFractionCostOnALargeGraph) {
    measure_real_fraction(300, 100'000);
}

#endif // _DNA_GRAPH || _DNA5_GRAPH

} // namespace
