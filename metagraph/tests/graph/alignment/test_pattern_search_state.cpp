/**
 * The search state of the extension (increment 5s-2; DECISIONS P29 of 2026-10-08): Pattern as
 * the Model, the SupportTracker protocol and the PathSink, on tiny graphs, against a brute-force
 * walk oracle that reads neither the engine's DFS nor its automaton:
 *  - the oracle reads a pattern through this file's own tables: the IUPAC codes as the bases
 *    each admits, and the standard genetic code (NCBI table 1) as a 64-letter string in TCAG
 *    order, with X, B, Z, J and '*' derived from it. An oriented pattern is a predicate on
 *    strings: s is a prefix of an instance (prefix_ok). For rc(P) of a peptide, rc(s) must be
 *    a suffix of an instance of P, codon by codon;
 *  - the walk tree: from every k-mer of the served graph (its valid edges, and the wrapper's
 *    reverse complements on a wrapped PRIMARY graph) that is a prefix of an instance, every
 *    extension by a base whose next k-mer is in the graph and keeps the prefix, to L bases. The
 *    tree's nodes are the branches the DFS can enter, its leaves at L the walks;
 *  - a script for the test tracker (which branches it declares DEAD, at which call it answers
 *    STOPPED) is simulated over the tree in the answer order (anchor node, orientation, then
 *    A < C < G < T): the expected walks, supported paths, prunings, branches entered, edges
 *    examined (the out-degree of every node expanded, over the graph's k-mers) and branchings.
 * The test tracker checks the protocol at every call (frame depth and position, the spelled
 * bases against the graph's node sequence, the stored node against a spelling of every stored
 * k-mer and stored_reverse_complement against the stored k-mer's spelling, the Model's state
 * against the oracle's allowed bases, pops matching pushes, frames unwound after a stop). Each test was run against mutants of the engine (see the progress
 * note of 5s-2): a tracker that prunes nothing must leave every answer of increment 4 as it was,
 * and one that prunes must give exactly the oracle's walks less the pruned subtrees.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"
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
using mtg::test::build_graph;
using mtg::test::unbounded_deadline;

typedef DeBruijnGraph::node_index node_index;
typedef SupportTracker::Verdict Verdict;

constexpr uint64_t kManySteps = 1'000'000'000;
constexpr char kBases[] = "ACGT";


// ---------------------------------------------------------------- the oracle's tables

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

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return '?';
    }
}

std::string rev_comp(std::string_view s) {
    std::string rc(s.rbegin(), s.rend());
    for (char &c : rc) {
        c = complement(c);
    }
    return rc;
}

// the standard genetic code (NCBI table 1), codons in TCAG order (TTT, TTC, TTA, TTG, TCT, ...)
constexpr char kStandardCode[]
        = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

char translate(std::string_view codon) {
    static const std::string order = "TCAG";
    size_t index = 0;
    for (char c : codon) {
        index = 4 * index + order.find(c);
    }
    return kStandardCode[index];
}

// a peptide letter as the set of amino acids it admits ('*' the stops)
bool residue_admits(char residue, char amino) {
    switch (residue) {
        case 'X': return amino != '*';
        case 'B': return amino == 'D' || amino == 'N';
        case 'Z': return amino == 'E' || amino == 'Q';
        case 'J': return amino == 'I' || amino == 'L';
        default: return amino == residue;
    }
}

const std::vector<std::string>& all_codons() {
    static const std::vector<std::string> codons = []() {
        std::vector<std::string> result;
        for (char a : std::string(kBases)) {
            for (char b : std::string(kBases)) {
                for (char c : std::string(kBases)) {
                    result.push_back({ a, b, c });
                }
            }
        }
        return result;
    }();
    return codons;
}

/**
 * An oriented pattern as the oracle reads it. DNA and IUPAC: the bases per position of this
 * orientation. PROTEIN: P's residues; for REVERSE (rc(P)) a string s is a prefix of an
 * instance iff rc(s) is a suffix of an instance of P.
 */
struct Oriented {
    Orientation orientation;
    PatternKind kind;
    size_t length = 0;
    std::vector<std::string> sets;
    std::string residues;
    bool reverse = false;

    bool prefix_ok(std::string_view s) const {
        if (s.size() > length)
            return false;
        if (kind != PatternKind::PROTEIN) {
            for (size_t i = 0; i < s.size(); ++i) {
                if (sets[i].find(s[i]) == std::string::npos)
                    return false;
            }
            return true;
        }
        if (!reverse) {
            // codon c covers [3c, 3c + 3): the part spelled must start a codon of residue c
            for (size_t c = 0; 3 * c < s.size(); ++c) {
                const std::string_view part = s.substr(3 * c, 3);
                if (!codon_with(residues[c], part, 0))
                    return false;
            }
            return true;
        }
        // rc(s) covers P's positions [length - |s|, length)
        const std::string t = rev_comp(s);
        const size_t offset = length - s.size();
        for (size_t c = offset / 3; 3 * c < length; ++c) {
            const size_t first = std::max(3 * c, offset);
            const std::string_view part
                = std::string_view(t).substr(first - offset, 3 * c + 3 - first);
            if (!codon_with(residues[c], part, first - 3 * c))
                return false;
        }
        return true;
    }

    // a codon of |residue| holding |part| from its position |begin| on
    static bool codon_with(char residue, std::string_view part, size_t begin) {
        for (const std::string &codon : all_codons()) {
            if (residue_admits(residue, translate(codon))
                    && std::string_view(codon).substr(begin, part.size()) == part)
                return true;
        }
        return false;
    }

    // the bases b (BaseSet bits A C G T) with s + b a prefix of an instance
    BaseSet next_bases(std::string_view s) const {
        BaseSet set = 0;
        for (size_t b = 0; b < 4; ++b) {
            if (prefix_ok(std::string(s) + kBases[b]))
                set |= BaseSet(1) << b;
        }
        return set;
    }
};

// the oriented patterns a request searches (palindromic DNA and IUPAC once; a peptide of the
// standard code never is)
std::vector<Oriented> orientations(PatternKind kind, const std::string &text, Strands strands) {
    std::vector<Oriented> result;
    if (kind == PatternKind::PROTEIN) {
        const size_t L = 3 * text.size();
        if (strands != Strands::REVERSE)
            result.push_back(Oriented { Orientation::FORWARD, kind, L, {}, text, false });
        if (strands != Strands::FORWARD)
            result.push_back(Oriented { Orientation::REVERSE, kind, L, {}, text, true });
        return result;
    }
    std::vector<std::string> q;
    for (char c : text) {
        q.push_back(iupac_bases(c));
    }
    std::vector<std::string> rc(q.rbegin(), q.rend());
    for (std::string &bases : rc) {
        for (char &b : bases) {
            b = complement(b);
        }
        std::sort(bases.begin(), bases.end());
    }
    if (q == rc)
        return { Oriented { Orientation::PALINDROMIC, kind, q.size(), q, {}, false } };
    if (strands != Strands::REVERSE)
        result.push_back(Oriented { Orientation::FORWARD, kind, q.size(), q, {}, false });
    if (strands != Strands::FORWARD)
        result.push_back(Oriented { Orientation::REVERSE, kind, q.size(), rc, {}, false });
    return result;
}

const Oriented& oriented_of(const std::vector<Oriented> &all, Orientation o) {
    for (const Oriented &x : all) {
        if (x.orientation == o)
            return x;
    }
    throw std::logic_error("orientation not searched");
}


// ---------------------------------------------------------------- graphs

const DBGSuccinct& base_dbg(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

// every k-mer of the served graph with its node id: the valid edges, and on a wrapped PRIMARY
// graph the wrapper's reverse complements
std::map<std::string, node_index> graph_kmers(const DeBruijnGraph &graph) {
    std::map<std::string, node_index> kmers;
    const DBGSuccinct &dbg_succ = base_dbg(graph);
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
        if (!dbg_succ.in_graph(y))
            continue;
        kmers.emplace(graph.get_node_sequence(y), y);
        if (canonical) {
            node_index z = canonical->reverse_complement(y);
            if (z != y)
                kmers.emplace(graph.get_node_sequence(z), z);
        }
    }
    return kmers;
}

// the stored node of every node of the served graph, found by spelling (as in
// test_pattern_search.cpp: never the engine's id arithmetic)
class StoredNodes {
  public:
    explicit StoredNodes(const DeBruijnGraph &graph)
          : graph_(graph), wrapped_(dynamic_cast<const CanonicalDBG*>(&graph)) {
        if (!wrapped_)
            return;
        const DBGSuccinct &dbg_succ = dynamic_cast<const DBGSuccinct&>(wrapped_->get_graph());
        for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
            if (dbg_succ.in_graph(y)) {
                stored_.emplace(dbg_succ.get_node_sequence(y), y);
                spelling_.emplace(y, dbg_succ.get_node_sequence(y));
            }
        }
    }

    node_index of(node_index node) const {
        if (!wrapped_)
            return node;
        const std::string kmer = graph_.get_node_sequence(node);
        auto it = stored_.find(kmer);
        if (it == stored_.end())
            it = stored_.find(rev_comp(kmer));
        return it == stored_.end() ? DeBruijnGraph::npos : it->second;
    }

    // the k-mer a stored node spells (on the wrapper: as the underlying graph stores it)
    std::string kmer_of(node_index stored) const {
        if (!wrapped_)
            return graph_.get_node_sequence(stored);
        auto it = spelling_.find(stored);
        return it == spelling_.end() ? std::string() : it->second;
    }

  private:
    const DeBruijnGraph &graph_;
    const CanonicalDBG *wrapped_;
    std::map<std::string, node_index> stored_;
    std::map<node_index, std::string> spelling_;
};


// ---------------------------------------------------------------- the walk tree

struct Branch {
    Orientation orientation;
    std::string sequence;
    std::vector<node_index> nodes;
    // the branches one base longer, in symbol order
    std::vector<size_t> children;
    // the graph's k-mers following the last one, whatever their base (edges examined when the
    // DFS expands this branch)
    uint64_t outdegree = 0;
};

struct WalkTree {
    size_t k = 0;
    size_t L = 0;
    std::vector<Branch> branches;
    // the roots (anchors), in answer order: (anchor node, orientation)
    std::vector<size_t> anchors;
    std::map<Orientation, uint64_t> anchors_of;
};

WalkTree walk_tree(const DeBruijnGraph &graph, const std::vector<Oriented> &searched) {
    WalkTree tree;
    tree.k = graph.get_k();
    tree.L = searched.front().length;
    const auto kmers = graph_kmers(graph);
    std::function<size_t(Branch)> grow = [&](Branch branch) {
        const size_t id = tree.branches.size();
        tree.branches.push_back(branch);
        if (branch.sequence.size() == tree.L)
            return id;
        const Oriented &o = oriented_of(searched, branch.orientation);
        const std::string suffix = branch.sequence.substr(branch.sequence.size() - tree.k + 1);
        std::vector<size_t> children;
        uint64_t outdegree = 0;
        for (char b : std::string(kBases)) {
            auto it = kmers.find(suffix + b);
            if (it == kmers.end())
                continue;
            ++outdegree;
            if (!o.prefix_ok(branch.sequence + b))
                continue;
            Branch child { branch.orientation, branch.sequence + b, branch.nodes, {}, 0 };
            child.nodes.push_back(it->second);
            children.push_back(grow(child));
        }
        tree.branches[id].children = children;
        tree.branches[id].outdegree = outdegree;
        return id;
    };
    std::vector<std::tuple<node_index, Orientation, std::string>> roots;
    for (const auto &[kmer, node] : kmers) {
        for (const Oriented &o : searched) {
            if (o.prefix_ok(kmer))
                roots.emplace_back(node, o.orientation, kmer);
        }
    }
    std::sort(roots.begin(), roots.end());
    for (const auto &[node, o, kmer] : roots) {
        tree.anchors.push_back(grow(Branch { o, kmer, { node }, {}, 0 }));
        ++tree.anchors_of[o];
    }
    return tree;
}


// ---------------------------------------------------------------- scripts and their simulation

typedef std::pair<Orientation, std::string> BranchKey;

/**
 * What the test tracker answers: DEAD at open() or push() for the branches in |dead|, DEAD at
 * complete() for the walks in |dead_complete|, STOPPED at the |stop_at|-th call (1-based, over
 * open, push and complete in the order the DFS makes them; 0: never).
 */
struct Script {
    std::set<BranchKey> dead;
    std::set<BranchKey> dead_complete;
    uint64_t stop_at = 0;
};

// what the engine must answer under a script, simulated over the walk tree
struct Expected {
    std::map<Orientation, uint64_t> walks;
    std::map<Orientation, uint64_t> supported;
    std::map<Orientation, bool> early;
    std::map<Orientation, uint64_t> extended;
    uint64_t pruned = 0;
    uint64_t candidates = 0;
    uint64_t edges = 0;
    uint64_t branchings = 0;
    uint64_t calls = 0;
    // the anchors whose extension began (Work::extension_anchors)
    uint64_t anchors_begun = 0;
    // a STOPPED verdict, or stop_at_threshold above |max_paths| supported paths
    bool stopped = false;
    bool threshold = false;
    // the supported walks completed, in answer order (branch ids)
    std::vector<size_t> paths;

    uint64_t total(const std::map<Orientation, uint64_t> &m) const {
        uint64_t sum = 0;
        for (const auto &[o, n] : m) {
            sum += n;
        }
        return sum;
    }
};

Expected simulate(const WalkTree &tree, const Script &script,
                  uint64_t threshold = std::numeric_limits<uint64_t>::max()) {
    Expected e;
    auto key = [&](const Branch &b) { return BranchKey { b.orientation, b.sequence }; };
    // false: stopped
    auto call = [&]() {
        ++e.calls;
        if (e.calls == script.stop_at) {
            e.stopped = true;
            return false;
        }
        return true;
    };
    std::function<bool(size_t)> visit = [&](size_t id) {
        const Branch &branch = tree.branches[id];
        e.edges += branch.outdegree;
        e.branchings += branch.children.size() > 1;
        for (size_t c : branch.children) {
            const Branch &child = tree.branches[c];
            const bool at_end = child.sequence.size() == tree.L;
            ++e.candidates;
            if (!call())
                return false;
            if (script.dead.count(key(child))) {
                ++e.pruned;
                if (at_end) {
                    ++e.walks[child.orientation];
                } else {
                    e.early[child.orientation] = true;
                }
                continue;
            }
            if (!at_end) {
                if (!visit(c))
                    return false;
                continue;
            }
            ++e.walks[child.orientation];
            if (!call())
                return false;
            if (script.dead_complete.count(key(child))) {
                ++e.pruned;
                continue;
            }
            ++e.supported[child.orientation];
            e.paths.push_back(c);
            if (e.total(e.supported) > threshold) {
                e.threshold = true;
                return false;
            }
        }
        return true;
    };
    for (size_t root : tree.anchors) {
        const Branch &anchor = tree.branches[root];
        ++e.anchors_begun;
        if (!call())
            break;
        if (script.dead.count(key(anchor))) {
            ++e.pruned;
            e.early[anchor.orientation] = true;
            ++e.extended[anchor.orientation];
            continue;
        }
        if (!visit(root))
            break;
        ++e.extended[anchor.orientation];
    }
    return e;
}


// ---------------------------------------------------------------- the test tracker and sink

std::string describe(const SearchState &s) {
    std::ostringstream out;
    out << "(" << orientation_key(s.orientation) << ", depth " << s.depth << ", position "
        << s.position << ", " << s.spelled << ")";
    return out.str();
}

/**
 * Answers by its Script and checks the protocol of SupportTracker at every call; every
 * deviation is kept in |errors| (the tests assert it empty).
 */
class ScriptedTracker : public SupportTracker {
  public:
    ScriptedTracker(const DeBruijnGraph &graph, const Pattern &pattern,
                    const std::vector<Oriented> &searched, Script script)
          : graph_(graph), stored_(graph), k_(graph.get_k()), pattern_(pattern),
            rc_(pattern.reverse_complement()), searched_(searched), script_(std::move(script)) {}

    Verdict open(const SearchState &s) override {
        ++opens;
        if (!stack_.empty())
            error("open() with frames open: " + describe(s));
        if (s.depth || s.position != k_ || s.base || s.spelled.size() != k_)
            error("open() of no anchor: " + describe(s));
        check(s);
        return verdict(s, script_.dead);
    }

    Verdict push(const SearchState &s) override {
        ++pushes;
        if (stack_.empty()) {
            error("push() without frames: " + describe(s));
        } else {
            const Frame &top = stack_.back();
            if (s.depth != stack_.size() || s.position != top.position + 1
                    || s.orientation != stack_.front().orientation
                    || s.spelled.substr(0, s.spelled.size() - 1) != top.spelled
                    || s.base != s.spelled.back())
                error("push() not one step from its frame: " + describe(s));
        }
        check(s);
        return verdict(s, script_.dead);
    }

    void pop() override {
        ++pops;
        if (stack_.empty()) {
            error("pop() without frames");
            return;
        }
        stack_.pop_back();
    }

    Verdict complete(const PathView &walk) override {
        ++completes;
        const std::string sequence(walk.sequence);
        bool same = walk.path.size() == stack_.size() && !stack_.empty()
                        && sequence == stack_.back().spelled
                        && walk.anchor.node == stack_.front().node
                        && walk.anchor.orientation == stack_.front().orientation
                        && walk.anchor.offset == 0
                        && walk.anchor.base_node == stored_.of(walk.anchor.node);
        for (size_t i = 0; same && i < walk.path.size(); ++i) {
            same = walk.path[i] == stack_[i].node;
        }
        if (!same || sequence.size() != pattern_.length())
            error("complete() of no walk of the frames: " + sequence);
        SearchState last {};
        last.orientation = walk.anchor.orientation;
        last.spelled = walk.sequence;
        last.position = static_cast<uint32_t>(sequence.size());
        return verdict(last, script_.dead_complete, false);
    }

    const char* stop_reason() const override { return stopped_ ? "scripted_stop" : nullptr; }

    size_t frames() const { return stack_.size(); }

    // the calls, and the ALIVE verdicts of open() and push() (the frames pushed)
    uint64_t opens = 0, pushes = 0, pops = 0, completes = 0, calls = 0, alive = 0;
    std::vector<std::string> errors;

  private:
    struct Frame {
        node_index node;
        Orientation orientation;
        uint32_t position;
        std::string spelled;
    };

    const DeBruijnGraph &graph_;
    StoredNodes stored_;
    size_t k_;
    const Pattern &pattern_;
    Pattern rc_;
    std::vector<Oriented> searched_;
    Script script_;
    std::vector<Frame> stack_;
    bool stopped_ = false;

    void error(const std::string &message) {
        if (errors.size() < 20)
            errors.push_back(message);
    }

    // the state of the node entered, against the graph and the oracle
    void check(const SearchState &s) {
        if (s.side != Side::RIGHT)
            error("a left step: " + describe(s));
        if (s.position != s.spelled.size() || s.position != k_ + s.depth)
            error("position: " + describe(s));
        if (s.spelled.size() >= k_
                && graph_.get_node_sequence(s.node) != s.spelled.substr(s.spelled.size() - k_))
            error("the node is not the last k-mer spelled: " + describe(s));
        if (s.base_node != stored_.of(s.node))
            error("base node: " + describe(s));
        // the stored k-mer is the node's k-mer as spelled, or its reverse complement (review of
        // 5s-2: what 5s-3's trackers pass support_step as the row's k-mer)
        if (s.spelled.size() >= k_) {
            const std::string kmer(s.spelled.substr(s.spelled.size() - k_));
            if (stored_.kmer_of(s.base_node)
                    != (s.stored_reverse_complement ? rev_comp(kmer) : kmer)) {
                error("stored_reverse_complement: " + describe(s));
            }
        }
        const Oriented &o = oriented_of(searched_, s.orientation);
        if (!o.prefix_ok(s.spelled))
            error("entered a branch off the pattern: " + describe(s));
        // the Model's state: the bases it allows next are the oracle's
        const Pattern &model = s.orientation == Orientation::REVERSE ? rc_ : pattern_;
        if (s.position < model.length()
                && model.bases(s.model, s.position) != o.next_bases(s.spelled))
            error("the Model's state allows other bases: " + describe(s));
    }

    Verdict verdict(const SearchState &s, const std::set<BranchKey> &dead, bool frame = true) {
        if (++calls == script_.stop_at) {
            stopped_ = true;
            return Verdict::STOPPED;
        }
        if (dead.count(BranchKey { s.orientation, std::string(s.spelled) }))
            return Verdict::DEAD;
        if (frame) {
            stack_.push_back(Frame { s.node, s.orientation, s.position, std::string(s.spelled) });
            ++alive;
        }
        return Verdict::ALIVE;
    }
};

// keeps a copy of every path it is given; answers false at its |stop_at|-th (1-based; 0: never)
class ListSink : public PathSink {
  public:
    explicit ListSink(const ScriptedTracker *tracker = nullptr, uint64_t stop_at = 0)
          : tracker_(tracker), stop_at_(stop_at) {}

    bool accept(const PathView &path, const SupportTracker *support) override {
        if (support != tracker_)
            errors.push_back("the sink was given another tracker");
        // the tracker's frames are still open: one per k-mer of the walk
        if (tracker_ && tracker_->frames() != path.path.size())
            errors.push_back("the frames are not open at the sink");
        if (++accepted == stop_at_) {
            stopped_ = true;
            return false;
        }
        paths.push_back(path.copy());
        return true;
    }

    const char* stop_reason() const override { return stopped_ ? "sink_stop" : nullptr; }

    std::vector<Context> paths;
    uint64_t accepted = 0;
    std::vector<std::string> errors;

  private:
    const ScriptedTracker *tracker_;
    uint64_t stop_at_;
    bool stopped_ = false;
};


// ---------------------------------------------------------------- running the engine

Request path_request(Strands strands = Strands::BOTH) {
    Request request;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_contexts = 1'000'000'000;
    request.max_anchors = 1'000'000'000;
    request.max_paths = 1'000'000'000;
    request.extend_paths = true;
    return request;
}

Budget unbounded_budget() {
    return Budget(kManySteps, unbounded_deadline());
}

struct PathCtx {
    node_index anchor;
    Orientation orientation;
    std::string sequence;
    std::vector<node_index> path;

    bool operator==(const PathCtx &o) const {
        return std::tie(anchor, orientation, sequence, path)
                == std::tie(o.anchor, o.orientation, o.sequence, o.path);
    }
};

std::ostream& operator<<(std::ostream &out, const PathCtx &p) {
    return out << "(" << p.anchor << ", " << orientation_key(p.orientation) << ", " << p.sequence
               << ")";
}

PathCtx path_of(const Context &c) {
    return PathCtx { c.node, c.orientation, c.sequence, c.path };
}

std::vector<PathCtx> paths_of(const WalkTree &tree, const std::vector<size_t> &ids) {
    std::vector<PathCtx> result;
    for (size_t id : ids) {
        const Branch &b = tree.branches[id];
        result.push_back(PathCtx { b.nodes.front(), b.orientation, b.sequence, b.nodes });
    }
    return result;
}

Result enumerate_paths(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                       std::vector<PathCtx> *paths, Budget *budget = nullptr) {
    Budget own = unbounded_budget();
    return PatternSearch(graph).enumerate(pattern, request, budget ? *budget : own,
                                          [&](const Context &c) { paths->push_back(path_of(c)); });
}

Result count_paths(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                   Budget *budget = nullptr) {
    Budget own = unbounded_budget();
    return PatternSearch(graph).count(pattern, request, budget ? *budget : own);
}

void expect_count(const Count &c, Relation relation, uint64_t value, const std::string &what) {
    EXPECT_EQ(Unit::PATHS, c.unit) << what;
    EXPECT_EQ(relation, c.relation) << what << ": " << to_string(c.relation);
    EXPECT_EQ(value, c.value) << what;
}

/**
 * Every count of a long pattern's extension under a tracker against the simulation: walks
 * and supported paths per orientation with their relations, the prunings, the work, the stop.
 */
void check_counts(const Result &r, const WalkTree &tree, const Expected &e,
                  const std::vector<Oriented> &searched, uint64_t discovery_steps,
                  const std::string &what) {
    ASSERT_TRUE(r.anchors) << what;
    const AnchorCounts &a = *r.anchors;
    const bool stopped = e.stopped || e.threshold;
    EXPECT_EQ(Relation::EXACT, a.total.relation) << what;
    EXPECT_EQ(tree.anchors.size(), a.total.value) << what;
    EXPECT_EQ(tree.anchors.empty() ? Extension::NO_ANCHORS
                                   : stopped ? Extension::STOPPED : Extension::COMPLETED,
              a.extension) << what;
    if (e.stopped) {
        ASSERT_TRUE(r.stop) << what;
        EXPECT_EQ(StopPhase::EXTENSION, r.stop->phase) << what;
        EXPECT_EQ(StopReason::EXTERNAL, r.stop->reason) << what;
    } else if (e.threshold) {
        ASSERT_TRUE(r.stop) << what;
        EXPECT_EQ(StopReason::MAX_PATHS, r.stop->reason) << what;
    } else {
        EXPECT_FALSE(r.stop) << what;
    }
    EXPECT_FALSE(r.time_limited) << what;

    bool all_exact = true;
    bool all_supported_exact = true;
    ASSERT_EQ(searched.size(), a.paths_by_orientation.size()) << what;
    ASSERT_EQ(searched.size(), a.supported_by_orientation.size()) << what;
    for (const Oriented &o : searched) {
        const Orientation x = o.orientation;
        const std::string per = what + " " + orientation_key(x);
        const uint64_t anchors = tree.anchors_of.count(x) ? tree.anchors_of.at(x) : 0;
        const bool extended = (e.extended.count(x) ? e.extended.at(x) : 0) == anchors;
        const bool early = e.early.count(x) && e.early.at(x);
        const uint64_t walks = e.walks.count(x) ? e.walks.at(x) : 0;
        const uint64_t supported = e.supported.count(x) ? e.supported.at(x) : 0;
        const Relation paths_relation = !anchors || (extended && !early) ? Relation::EXACT
                                                                         : Relation::AT_LEAST;
        const Relation supported_relation = !anchors || extended ? Relation::EXACT
                                                                 : Relation::AT_LEAST;
        expect_count(a.paths_by_orientation.at(x), paths_relation, walks, per + " walks");
        expect_count(a.supported_by_orientation.at(x), supported_relation, supported,
                     per + " supported");
        all_exact &= paths_relation == Relation::EXACT;
        all_supported_exact &= supported_relation == Relation::EXACT;
    }
    expect_count(a.paths, all_exact ? Relation::EXACT : Relation::AT_LEAST, e.total(e.walks),
                 what + " walks");
    expect_count(a.supported, all_supported_exact ? Relation::EXACT : Relation::AT_LEAST,
                 e.total(e.supported), what + " supported");
    EXPECT_EQ(e.total(e.walks), a.walks) << what;
    EXPECT_EQ(e.pruned, a.branches_pruned) << what;
    bool any_early = false;
    for (const auto &[o, early] : e.early) {
        any_early |= early;
    }
    EXPECT_EQ(any_early, a.pruned_before_completion) << what;
    EXPECT_EQ(e.candidates, a.candidates_examined) << what;
    EXPECT_EQ(e.edges, r.work.extension_edges) << what;
    EXPECT_EQ(e.branchings, r.work.extension_branches) << what;
    // every anchor listed is spelled and handed to the tracker, up to the one it stopped at
    EXPECT_EQ(e.anchors_begun, r.work.extension_anchors) << what;
    EXPECT_EQ(discovery_steps + e.edges, r.work.steps) << what;
}

/**
 * One (graph, pattern, strands, script): count() with the tracker, count() with the tracker and
 * a sink, enumerate() with the tracker in ALL_OR_COUNT and in PARTIAL at several caps, and
 * stop_at_threshold, each against the simulation of the script over the walk tree.
 */
void check_script(const DeBruijnGraph &graph, const Pattern &pattern, PatternKind kind,
                  const std::string &text, Strands strands, const Script &script,
                  const std::string &what) {
    const std::vector<Oriented> searched = orientations(kind, text, strands);
    const WalkTree tree = walk_tree(graph, searched);
    const Expected e = simulate(tree, script);
    const std::vector<PathCtx> expected = paths_of(tree, e.paths);

    Request request = path_request(strands);
    Request anchors_only = request;
    anchors_only.extend_paths = false;
    const uint64_t discovery = count_paths(graph, pattern, anchors_only).work.steps;

    // count(): the tracker alone
    {
        ScriptedTracker tracker(graph, pattern, searched, script);
        Request rq = request;
        rq.support = &tracker;
        Result r = count_paths(graph, pattern, rq);
        check_counts(r, tree, e, searched, discovery, what + " count");
        EXPECT_FALSE(r.extraction);
        EXPECT_EQ(0u, tracker.frames()) << what << ": frames left open";
        EXPECT_EQ(tracker.alive, tracker.pops)
            << what << ": every ALIVE open() and push() popped once";
        EXPECT_TRUE(tracker.errors.empty()) << what << ": " << tracker.errors.front();
        EXPECT_EQ(e.calls, tracker.calls) << what;
        if (e.stopped) {
            EXPECT_STREQ("scripted_stop", tracker.stop_reason());
        }
    }
    // count() with a sink: the supported walks as they complete, in answer order
    {
        ScriptedTracker tracker(graph, pattern, searched, script);
        ListSink sink(&tracker);
        Request rq = request;
        rq.support = &tracker;
        rq.sink = &sink;
        Result r = count_paths(graph, pattern, rq);
        check_counts(r, tree, e, searched, discovery, what + " sink");
        std::vector<PathCtx> got;
        for (const Context &c : sink.paths) {
            got.push_back(path_of(c));
            EXPECT_EQ(0u, c.offset);
        }
        EXPECT_EQ(expected, got) << what;
        EXPECT_TRUE(sink.errors.empty()) << what << ": " << sink.errors.front();
        EXPECT_TRUE(tracker.errors.empty()) << what << ": " << tracker.errors.front();
        EXPECT_EQ(0u, tracker.frames()) << what;
    }
    // enumerate(), ALL_OR_COUNT: the supported paths, or none after a stop
    {
        ScriptedTracker tracker(graph, pattern, searched, script);
        Request rq = request;
        rq.support = &tracker;
        std::vector<PathCtx> got;
        Result r = enumerate_paths(graph, pattern, rq, &got);
        check_counts(r, tree, e, searched, discovery, what + " all_or_count");
        ASSERT_TRUE(r.extraction);
        if (e.stopped) {
            EXPECT_EQ(Withheld::EXTERNAL, r.extraction->withheld) << what;
            EXPECT_TRUE(got.empty()) << what;
            EXPECT_FALSE(r.extraction->complete) << what;
        } else {
            EXPECT_TRUE(r.extraction->complete) << what;
            EXPECT_FALSE(r.extraction->withheld) << what;
            EXPECT_EQ(expected, got) << what;
            EXPECT_EQ(expected.size(), r.extraction->returned) << what;
        }
        EXPECT_EQ(0u, tracker.frames()) << what;
        EXPECT_TRUE(tracker.errors.empty()) << what << ": " << tracker.errors.front();
    }
    // ALL_OR_COUNT at the threshold: the supported count decides
    if (!e.stopped && expected.size()) {
        ScriptedTracker tracker(graph, pattern, searched, script);
        Request rq = request;
        rq.support = &tracker;
        rq.max_paths = expected.size() - 1;
        std::vector<PathCtx> got;
        Result r = enumerate_paths(graph, pattern, rq, &got);
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, r.extraction->withheld) << what;
        EXPECT_TRUE(got.empty()) << what;
        rq.max_paths = expected.size();
        ScriptedTracker again(graph, pattern, searched, script);
        rq.support = &again;
        got.clear();
        r = enumerate_paths(graph, pattern, rq, &got);
        EXPECT_TRUE(r.extraction->complete) << what;
        EXPECT_EQ(expected, got) << what;
    }
    // PARTIAL: the first max_paths supported paths, the cut stated
    for (uint64_t cap : { uint64_t(0), uint64_t(1), uint64_t(expected.size() / 2),
                          uint64_t(expected.size()), uint64_t(expected.size() + 1) }) {
        ScriptedTracker tracker(graph, pattern, searched, script);
        Request rq = request;
        rq.support = &tracker;
        rq.mode = Mode::PARTIAL;
        rq.max_paths = cap;
        std::vector<PathCtx> got;
        Result r = enumerate_paths(graph, pattern, rq, &got);
        const std::string at = what + " partial " + std::to_string(cap);
        check_counts(r, tree, e, searched, discovery, at);
        ASSERT_EQ(std::min<uint64_t>(cap, expected.size()), got.size()) << at;
        EXPECT_TRUE(std::equal(got.begin(), got.end(), expected.begin())) << at;
        if (e.stopped) {
            EXPECT_EQ(StopReason::EXTERNAL, r.extraction->cut) << at;
            EXPECT_FALSE(r.extraction->complete) << at;
        } else if (cap >= expected.size()) {
            EXPECT_TRUE(r.extraction->complete) << at;
            EXPECT_FALSE(r.extraction->cut) << at;
        } else {
            EXPECT_EQ(StopReason::MAX_PATHS, r.extraction->cut) << at;
        }
        EXPECT_EQ(0u, tracker.frames()) << at;
    }
    // stop_at_threshold: once more than max_paths SUPPORTED paths are complete
    if (!e.stopped && expected.size()) {
        const uint64_t threshold = expected.size() / 2;
        const Expected t = simulate(tree, script, threshold);
        ASSERT_TRUE(t.threshold) << what;
        ScriptedTracker tracker(graph, pattern, searched, script);
        Request rq = request;
        rq.support = &tracker;
        rq.stop_at_threshold = true;
        rq.max_paths = threshold;
        Result r = count_paths(graph, pattern, rq);
        check_counts(r, tree, t, searched, discovery, what + " stop_at_threshold");
        EXPECT_EQ(0u, tracker.frames()) << what;
    }
}


// ---------------------------------------------------------------- random cases

struct Case {
    std::vector<std::string> records;
    size_t k;
    DeBruijnGraph::Mode mode;
    PatternKind kind;
    std::string text;
    Strands strands;
};

std::string random_dna(std::mt19937 &rng, size_t n) {
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        s.push_back(kBases[rng() % 4]);
    }
    return s;
}

// a few records with repeats planted (forward and reverse-complemented), so that walks branch
// and join, and a pattern cut from one of them, made degenerate at random positions
Case random_case(std::mt19937 &rng) {
    Case c;
    c.k = 3 + rng() % 4;
    const size_t records = 2 + rng() % 4;
    for (size_t i = 0; i < records; ++i) {
        c.records.push_back(random_dna(rng, 12 + rng() % 20));
    }
    for (size_t r = 0; r < 3; ++r) {
        std::string &from = c.records[rng() % records];
        std::string &to = c.records[rng() % records];
        const size_t n = c.k + rng() % 4;
        if (from.size() <= n)
            continue;
        std::string piece = from.substr(rng() % (from.size() - n), n);
        if (rng() % 2)
            piece = rev_comp(piece);
        to.insert(rng() % to.size(), piece);
    }
    const DeBruijnGraph::Mode modes[]
        = { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY };
    c.mode = modes[rng() % 3];
    const Strands strands[] = { Strands::BOTH, Strands::FORWARD, Strands::REVERSE };
    c.strands = strands[rng() % 3];

    const std::string &source = c.records[rng() % records];
    if (rng() % 3 == 0) {
        // a peptide of m residues, 3m > k, translated from the record in the oracle's code
        size_t m = c.k / 3 + 1 + rng() % 3;
        while (3 * m > source.size()) {
            --m;
        }
        if (3 * m <= c.k)
            return random_case(rng);
        const size_t begin = rng() % (source.size() - 3 * m + 1);
        c.kind = PatternKind::PROTEIN;
        for (size_t i = 0; i < m; ++i) {
            char residue = translate(source.substr(begin + 3 * i, 3));
            const uint32_t roll = rng() % 10;
            if (roll == 0 && residue != '*') {
                residue = 'X';
            } else if (roll == 1 && (residue == 'D' || residue == 'N')) {
                residue = 'B';
            } else if (roll == 1 && (residue == 'I' || residue == 'L')) {
                residue = 'J';
            }
            c.text.push_back(residue);
        }
        return c;
    }
    const size_t L = std::min(source.size(), c.k + 1 + rng() % (2 * c.k));
    if (L <= c.k)
        return random_case(rng);
    const size_t begin = rng() % (source.size() - L + 1);
    c.text = source.substr(begin, L);
    c.kind = PatternKind::DNA;
    for (char &base : c.text) {
        const uint32_t roll = rng() % 8;
        if (roll == 0) {
            base = 'N';
            c.kind = PatternKind::IUPAC;
        } else if (roll == 1) {
            base = base == 'A' || base == 'G' ? 'R' : 'Y';
            c.kind = PatternKind::IUPAC;
        }
    }
    return c;
}

std::string describe(const Case &c) {
    std::ostringstream out;
    out << "k " << c.k << " mode " << c.mode << " " << to_string(c.kind) << " " << c.text
        << " strands " << to_string(c.strands) << " records";
    for (const std::string &r : c.records) {
        out << " " << r;
    }
    return out.str();
}

// every branch of the tree marked DEAD at random (anchors included), every walk DEAD at its
// completion at random
Script random_script(std::mt19937 &rng, const WalkTree &tree, uint32_t per_mille) {
    Script script;
    for (const Branch &b : tree.branches) {
        if (rng() % 1000 < per_mille)
            script.dead.emplace(b.orientation, b.sequence);
        if (b.sequence.size() == tree.L && rng() % 1000 < per_mille)
            script.dead_complete.emplace(b.orientation, b.sequence);
    }
    return script;
}


// ---------------------------------------------------------------- the tests

TEST(PatternSearchState, PatternIsAModel) {
    // along random instances' prefixes, the Model's state allows exactly the oracle's bases
    // (and Pattern::allowed's), from a start window of any length; at L it accepts
    std::mt19937 rng(20261008);
    uint64_t checked = 0;
    for (size_t round = 0; round < 400; ++round) {
        const bool protein = round % 2;
        std::string text;
        if (protein) {
            const std::string letters = "ACDEFGHIKLMNPQRSTVWYXBZJ*";
            for (size_t i = 0, m = 2 + rng() % 5; i < m; ++i) {
                text.push_back(letters[rng() % letters.size()]);
            }
        } else {
            const std::string codes = "ACGTRYSWKMBDHVN";
            for (size_t i = 0, L = 3 + rng() % 10; i < L; ++i) {
                text.push_back(codes[rng() % (rng() % 2 ? 4 : codes.size())]);
            }
        }
        const PatternKind kind = protein ? PatternKind::PROTEIN : PatternKind::IUPAC;
        const Pattern p = Pattern::parse(kind, text);
        const Pattern rc = p.reverse_complement();
        for (const Oriented &o : orientations(kind, text, Strands::BOTH)) {
            const Pattern &model = o.orientation == Orientation::REVERSE ? rc : p;
            ASSERT_EQ(o.length, model.length());
            // a random prefix of an instance, as long as the oracle allows one
            const size_t start = 1 + rng() % (o.length - 1);
            std::string s;
            while (s.size() < start) {
                const BaseSet next = o.next_bases(s);
                if (!next)
                    break;
                std::vector<char> choice;
                for (size_t b = 0; b < 4; ++b) {
                    if (next & (1 << b))
                        choice.push_back(kBases[b]);
                }
                s.push_back(choice[rng() % choice.size()]);
            }
            if (s.size() < start)
                continue;
            Model::State state = model.start(s);
            for (uint32_t position = static_cast<uint32_t>(s.size()); position < o.length;
                    ++position) {
                const BaseSet next = o.next_bases(s);
                ASSERT_EQ(next, model.bases(state, position)) << text << " " << s;
                ASSERT_EQ(next, model.allowed(position, s)) << text << " " << s;
                ++checked;
                if (!next)
                    break;
                std::vector<char> choice;
                for (size_t b = 0; b < 4; ++b) {
                    if (next & (1 << b))
                        choice.push_back(kBases[b]);
                }
                const char base = choice[rng() % choice.size()];
                state = model.next(state, position, base);
                s.push_back(base);
            }
            if (s.size() == o.length) {
                EXPECT_TRUE(o.prefix_ok(s));
                EXPECT_TRUE(model.accepting(state, static_cast<uint32_t>(s.size())));
            }
        }
    }
    EXPECT_LT(1000u, checked);

    // a start window whose codon prefix holds a base other than A, C, G, T admits nothing, as
    // allowed() answers (an N of a DNA5 graph's k-mer)
    const Pattern pep = Pattern::parse(PatternKind::PROTEIN, "MKV");
    EXPECT_EQ(0, pep.bases(pep.start("ATGAN"), 5));
    EXPECT_EQ(0, pep.allowed(5, "ATGAN"));
    EXPECT_NE(0, pep.bases(pep.start("NTGAA"), 5));
    EXPECT_EQ(pep.allowed(5, "NTGAA"), pep.bases(pep.start("NTGAA"), 5));
    // DNA and IUPAC states carry nothing
    const Pattern dna = Pattern::parse(PatternKind::IUPAC, "ACGTN");
    EXPECT_EQ(0u, dna.start("ACG"));
    EXPECT_EQ(kAllBases, dna.bases(dna.next(dna.start("ACG"), 3, 'T'), 4));
}

TEST(PatternSearchState, NullTrackerAndSinkAnswerAsIncrement4) {
    // a request without tracker and sink, and one with a tracker that prunes nothing, give
    // increment 4's paths, counts and work; the latter adds the supported counts (all walks)
    std::mt19937 rng(42);
    uint64_t with_paths = 0;
    for (size_t round = 0; round < 120; ++round) {
        const Case c = random_case(rng);
        const std::string what = describe(c);
        auto graph = build_graph<DBGSuccinct>(c.k, c.records, c.mode);
        const Pattern pattern = Pattern::parse(c.kind, c.text);
        const std::vector<Oriented> searched = orientations(c.kind, c.text, c.strands);
        const WalkTree tree = walk_tree(*graph, searched);
        const Expected e = simulate(tree, Script {});
        const std::vector<PathCtx> expected = paths_of(tree, e.paths);
        with_paths += !expected.empty();

        const Request request = path_request(c.strands);
        std::vector<PathCtx> plain_paths;
        Result plain = enumerate_paths(*graph, pattern, request, &plain_paths);
        ASSERT_EQ(expected, plain_paths) << what;
        ASSERT_TRUE(plain.anchors) << what;
        EXPECT_EQ(e.total(e.walks), plain.anchors->walks) << what;
        // no tracker: no support stated
        EXPECT_EQ(Relation::UNKNOWN, plain.anchors->supported.relation) << what;
        EXPECT_TRUE(plain.anchors->supported_by_orientation.empty()) << what;
        EXPECT_EQ(0u, plain.anchors->branches_pruned) << what;

        ScriptedTracker tracker(*graph, pattern, searched, Script {});
        Request tracked = request;
        tracked.support = &tracker;
        std::vector<PathCtx> tracked_paths;
        Result r = enumerate_paths(*graph, pattern, tracked, &tracked_paths);
        EXPECT_EQ(plain_paths, tracked_paths) << what;
        EXPECT_TRUE(tracker.errors.empty()) << what << ": " << tracker.errors.front();
        EXPECT_EQ(0u, tracker.frames()) << what;
        const AnchorCounts &a = *plain.anchors;
        const AnchorCounts &b = *r.anchors;
        EXPECT_EQ(a.paths.relation, b.paths.relation) << what;
        EXPECT_EQ(a.paths.value, b.paths.value) << what;
        EXPECT_EQ(a.paths.relation, b.supported.relation) << what;
        EXPECT_EQ(a.paths.value, b.supported.value) << what;
        EXPECT_EQ(a.candidates_examined, b.candidates_examined) << what;
        EXPECT_EQ(a.extension, b.extension) << what;
        EXPECT_EQ(0u, b.branches_pruned) << what;
        EXPECT_FALSE(b.pruned_before_completion) << what;
        EXPECT_EQ(plain.work.steps, r.work.steps) << what;
        EXPECT_EQ(plain.work.extension_edges, r.work.extension_edges) << what;
        EXPECT_EQ(plain.work.extension_branches, r.work.extension_branches) << what;
        EXPECT_EQ(plain.extraction->complete, r.extraction->complete) << what;
        for (const auto &[o, count] : a.paths_by_orientation) {
            EXPECT_EQ(count.value, b.paths_by_orientation.at(o).value) << what;
            EXPECT_EQ(count.value, b.supported_by_orientation.at(o).value) << what;
        }
        // the tracker saw every branch of the tree: one open per anchor, one push per branch
        // entered, one complete per walk, every frame popped
        EXPECT_EQ(tree.anchors.size(), tracker.opens) << what;
        EXPECT_EQ(e.candidates, tracker.pushes) << what;
        EXPECT_EQ(e.total(e.walks), tracker.completes) << what;
        EXPECT_EQ(tracker.opens + tracker.pushes, tracker.pops) << what;

        // a sink without a tracker receives every walk, in answer order (count() only)
        ListSink sink;
        Request sunk = request;
        sunk.sink = &sink;
        Result s = count_paths(*graph, pattern, sunk);
        std::vector<PathCtx> sunk_paths;
        for (const Context &x : sink.paths) {
            sunk_paths.push_back(path_of(x));
        }
        EXPECT_EQ(expected, sunk_paths) << what;
        EXPECT_EQ(plain.anchors->paths.value, s.anchors->paths.value) << what;
        EXPECT_EQ(plain.work.steps, s.work.steps) << what;
        EXPECT_TRUE(sink.errors.empty()) << what;
    }
    std::cerr << "random cases with paths: " << with_paths << " of 120" << std::endl;
    EXPECT_LT(40u, with_paths);
}

TEST(PatternSearchState, DeadBranchesAreTheOracleLessTheirSubtrees) {
    // a tracker declaring chosen branches DEAD (anchors, inner branches, last k-mers, and
    // complete walks at complete()) yields exactly the oracle's walks less those subtrees,
    // with every counter and relation of the simulation
    std::mt19937 rng(7);
    uint64_t pruned_cases = 0;
    uint64_t early = 0;
    uint64_t at_end = 0;
    uint64_t supported = 0;
    uint64_t kinds[3] = { 0, 0, 0 };
    uint64_t largest = 0;
    for (size_t round = 0; round < 80; ++round) {
        const Case c = random_case(rng);
        auto graph = build_graph<DBGSuccinct>(c.k, c.records, c.mode);
        const Pattern pattern = Pattern::parse(c.kind, c.text);
        const WalkTree tree = walk_tree(*graph, orientations(c.kind, c.text, c.strands));
        largest = std::max<uint64_t>(largest, tree.branches.size());
        for (uint32_t per_mille : { 100u, 300u, 600u }) {
            const Script script = random_script(rng, tree, per_mille);
            const Expected e = simulate(tree, script);
            pruned_cases += e.pruned > 0;
            early += !e.early.empty();
            at_end += e.pruned > 0 && e.early.empty();
            supported += !e.paths.empty() && e.pruned > 0;
            kinds[static_cast<int>(c.kind)] += e.pruned > 0;
            check_script(*graph, pattern, c.kind, c.text, c.strands, script,
                         describe(c) + " dead " + std::to_string(per_mille));
            if (HasFatalFailure())
                return;
        }
    }
    std::cerr << "scripts with a pruning: " << pruned_cases << " (before L: " << early
              << ", only at L: " << at_end << ", with supported paths left: " << supported
              << "; dna " << kinds[0] << ", iupac " << kinds[1] << ", protein " << kinds[2]
              << "), largest tree " << largest << " branches" << std::endl;
    EXPECT_LT(60u, pruned_cases);
    EXPECT_LT(20u, early);
    EXPECT_LT(5u, at_end);
    EXPECT_LT(20u, supported);
    EXPECT_LT(5u, kinds[2]);
}

TEST(PatternSearchState, StoppedAtEveryCall) {
    // a tracker answering STOPPED at its n-th call, for every n: stop {extension, external},
    // the supported paths completed before it (PARTIAL releases them, cut external; ALL_OR_COUNT
    // withholds), every frame unwound
    std::mt19937 rng(11);
    uint64_t stops = 0;
    for (size_t round = 0; round < 30; ++round) {
        const Case c = random_case(rng);
        auto graph = build_graph<DBGSuccinct>(c.k, c.records, c.mode);
        const Pattern pattern = Pattern::parse(c.kind, c.text);
        const WalkTree tree = walk_tree(*graph, orientations(c.kind, c.text, c.strands));
        Script script = random_script(rng, tree, 150);
        const uint64_t calls = simulate(tree, script).calls;
        for (uint64_t n = 1; n <= calls && n <= 60; ++n) {
            script.stop_at = n;
            check_script(*graph, pattern, c.kind, c.text, c.strands, script,
                         describe(c) + " stop at " + std::to_string(n));
            ++stops;
            if (HasFatalFailure())
                return;
        }
    }
    std::cerr << "stops: " << stops << std::endl;
    EXPECT_LT(100u, stops);
}

TEST(PatternSearchState, RecombinantBubble) {
    // TESTS §4's bubble (k = 4): A ACGTTGCA, B TCGTTGCC; WCGTTGCM has four walks, two of them
    // mosaics. A tracker that kills each mosaic at its last k-mer leaves paths exact 4 (they
    // are complete walks), supported exact 2, two branches pruned, nothing pruned early
    const std::vector<std::string> records { "ACGTTGCA", "TCGTTGCC" };
    auto graph = build_graph<DBGSuccinct>(4, records, DeBruijnGraph::BASIC);
    const Pattern pattern = Pattern::parse(PatternKind::IUPAC, "WCGTTGCM");
    Script script;
    script.dead = { { Orientation::FORWARD, "ACGTTGCC" }, { Orientation::FORWARD, "TCGTTGCA" } };
    const std::vector<Oriented> searched
        = orientations(PatternKind::IUPAC, "WCGTTGCM", Strands::FORWARD);
    ScriptedTracker tracker(*graph, pattern, searched, script);
    Request request = path_request(Strands::FORWARD);
    request.support = &tracker;
    std::vector<PathCtx> got;
    Result r = enumerate_paths(*graph, pattern, request, &got);
    ASSERT_TRUE(r.anchors);
    expect_count(r.anchors->paths, Relation::EXACT, 4, "walks");
    expect_count(r.anchors->supported, Relation::EXACT, 2, "supported");
    EXPECT_EQ(4u, r.anchors->walks);
    EXPECT_EQ(2u, r.anchors->branches_pruned);
    EXPECT_FALSE(r.anchors->pruned_before_completion);
    ASSERT_EQ(2u, got.size());
    EXPECT_EQ("ACGTTGCA", got[0].sequence);
    EXPECT_EQ("TCGTTGCC", got[1].sequence);
    EXPECT_TRUE(r.extraction->complete);
    EXPECT_TRUE(tracker.errors.empty()) << tracker.errors.front();

    // MosaicDiesBeforeL: one base more (A ACGTTGCAA, B TCGTTGCCA, WCGTTGCMA): the mosaics die
    // at k-mer 4 of 6 -- paths at_least 2 (the walks below were not followed), supported exact 2
    const std::vector<std::string> longer { "ACGTTGCAA", "TCGTTGCCA" };
    auto graph2 = build_graph<DBGSuccinct>(4, longer, DeBruijnGraph::BASIC);
    const Pattern p2 = Pattern::parse(PatternKind::IUPAC, "WCGTTGCMA");
    Script s2;
    s2.dead = { { Orientation::FORWARD, "ACGTTGCC" }, { Orientation::FORWARD, "TCGTTGCA" } };
    ScriptedTracker t2(*graph2, p2, orientations(PatternKind::IUPAC, "WCGTTGCMA",
                                                  Strands::FORWARD), s2);
    Request rq2 = path_request(Strands::FORWARD);
    rq2.support = &t2;
    std::vector<PathCtx> got2;
    Result r2 = enumerate_paths(*graph2, p2, rq2, &got2);
    expect_count(r2.anchors->paths, Relation::AT_LEAST, 2, "walks");
    expect_count(r2.anchors->supported, Relation::EXACT, 2, "supported");
    EXPECT_EQ(2u, r2.anchors->branches_pruned);
    EXPECT_TRUE(r2.anchors->pruned_before_completion);
    EXPECT_EQ(Extension::COMPLETED, r2.anchors->extension);
    ASSERT_EQ(2u, got2.size());
    EXPECT_TRUE(r2.extraction->complete);
    EXPECT_FALSE(r2.stop);
}

TEST(PatternSearchState, TrackerTimeIsATimeStop) {
    // a tracker that finds the work time passed through the request's Budget: stop {extension,
    // time}, time_limited, the budget stopped for the next pattern; an external stop leaves the
    // budget to the next pattern
    const std::vector<std::string> records { "ACGTTGCAAGTC", "TCGTTGCCAGTT" };
    auto graph = build_graph<DBGSuccinct>(4, records, DeBruijnGraph::BASIC);
    const Pattern pattern = Pattern::parse(PatternKind::IUPAC, "WCGTTGCMA");
    const Pattern next = Pattern::parse(PatternKind::DNA, "GTTGCAAG");

    auto now = Deadline::Clock::now();
    Deadline::Clock::time_point fake = now;
    Budget budget(kManySteps, Deadline(now, 1000, 10, [&]() { return fake; }));

    class TimeTracker : public SupportTracker {
      public:
        TimeTracker(Budget &budget, Deadline::Clock::time_point &fake,
                    Deadline::Clock::time_point late)
              : budget_(budget), fake_(fake), late_(late) {}
        Verdict open(const SearchState &) override { ++depth; return Verdict::ALIVE; }
        Verdict push(const SearchState &s) override {
            if (s.depth == 3) {
                fake_ = late_;
                if (!budget_.check_time())
                    return Verdict::STOPPED;
            }
            ++depth;
            return Verdict::ALIVE;
        }
        void pop() override { --depth; }
        Verdict complete(const PathView &) override { return Verdict::ALIVE; }
        const char* stop_reason() const override { return "time"; }
        int depth = 0;
      private:
        Budget &budget_;
        Deadline::Clock::time_point &fake_;
        Deadline::Clock::time_point late_;
    } tracker(budget, fake, now + std::chrono::seconds(5));

    Request request = path_request(Strands::FORWARD);
    request.support = &tracker;
    Result r = PatternSearch(*graph).count(pattern, request, budget);
    ASSERT_TRUE(r.stop);
    EXPECT_EQ(StopPhase::EXTENSION, r.stop->phase);
    EXPECT_EQ(StopReason::TIME, r.stop->reason);
    EXPECT_TRUE(r.time_limited);
    EXPECT_EQ(Extension::STOPPED, r.anchors->extension);
    EXPECT_EQ(0, tracker.depth);
    Request plain = path_request(Strands::FORWARD);
    Result after = PatternSearch(*graph).count(next, plain, budget);
    ASSERT_TRUE(after.stop);
    EXPECT_EQ(StopPhase::DISCOVERY, after.stop->phase);
    EXPECT_EQ(StopReason::TIME, after.stop->reason);

    // the same stop without the budget's reading: external, and the next pattern runs
    Budget open_budget = unbounded_budget();
    Script script;
    script.stop_at = 4;
    ScriptedTracker scripted(*graph, pattern,
                             orientations(PatternKind::IUPAC, "WCGTTGCMA", Strands::FORWARD),
                             script);
    request.support = &scripted;
    Result external = PatternSearch(*graph).count(pattern, request, open_budget);
    ASSERT_TRUE(external.stop);
    EXPECT_EQ(StopReason::EXTERNAL, external.stop->reason);
    EXPECT_FALSE(external.time_limited);
    EXPECT_STREQ("external", to_string(external.stop->reason));
    EXPECT_FALSE(open_budget.stopped());
    Result runs = PatternSearch(*graph).count(next, plain, open_budget);
    EXPECT_FALSE(runs.stop);
    EXPECT_EQ(Extension::COMPLETED, runs.anchors->extension);
}

TEST(PatternSearchState, SinkStopsAndRules) {
    const std::vector<std::string> records { "ACGTTGCAAGTC", "TCGTTGCCAGTT", "ACGTTGCCAGTC" };
    auto graph = build_graph<DBGSuccinct>(4, records, DeBruijnGraph::BASIC);
    const std::string text = "WCGTTGCMAGT";
    const Pattern pattern = Pattern::parse(PatternKind::IUPAC, text);
    const std::vector<Oriented> searched = orientations(PatternKind::IUPAC, text, Strands::BOTH);
    const WalkTree tree = walk_tree(*graph, searched);
    const Expected e = simulate(tree, Script {});
    ASSERT_LE(3u, e.paths.size());

    // a sink answering false at its second path: stop {extension, external}, the walks counted
    // up to it
    ScriptedTracker tracker(*graph, pattern, searched, Script {});
    ListSink sink(&tracker, 2);
    Request request = path_request();
    request.support = &tracker;
    request.sink = &sink;
    Result r = count_paths(*graph, pattern, request);
    ASSERT_TRUE(r.stop);
    EXPECT_EQ(StopPhase::EXTENSION, r.stop->phase);
    EXPECT_EQ(StopReason::EXTERNAL, r.stop->reason);
    EXPECT_STREQ("sink_stop", sink.stop_reason());
    EXPECT_EQ(nullptr, tracker.stop_reason());
    EXPECT_EQ(1u, sink.paths.size());
    EXPECT_EQ(2u, sink.accepted);
    EXPECT_EQ(Relation::AT_LEAST, r.anchors->supported.relation);
    EXPECT_EQ(2u, r.anchors->supported.value);
    EXPECT_EQ(0u, tracker.frames());
    EXPECT_TRUE(tracker.errors.empty()) << tracker.errors.front();

    // stop_at_threshold reads the supported paths after the sink's call
    ScriptedTracker t2(*graph, pattern, searched, Script {});
    ListSink s2(&t2);
    request.support = &t2;
    request.sink = &s2;
    request.stop_at_threshold = true;
    request.max_paths = 1;
    Result r2 = count_paths(*graph, pattern, request);
    ASSERT_TRUE(r2.stop);
    EXPECT_EQ(StopReason::MAX_PATHS, r2.stop->reason);
    EXPECT_EQ(2u, s2.paths.size());

    // enumerate() refuses a sink for an extending pattern, not for a short one
    ListSink s3;
    Request with_sink = path_request();
    with_sink.sink = &s3;
    Budget budget = unbounded_budget();
    EXPECT_THROW(PatternSearch(*graph).enumerate(pattern, with_sink, budget,
                                                 [](const Context &) {}),
                 std::invalid_argument);
    const Pattern short_pattern = Pattern::parse(PatternKind::DNA, "CGTT");
    std::vector<PathCtx> none;
    EXPECT_NO_THROW(enumerate_paths(*graph, short_pattern, with_sink, &none));
    EXPECT_TRUE(s3.paths.empty());
    // a long pattern without the extension: the sink is never asked
    Request anchors_only = with_sink;
    anchors_only.extend_paths = false;
    Result a = count_paths(*graph, pattern, anchors_only);
    EXPECT_EQ(Extension::NOT_REQUESTED, a.anchors->extension);
    EXPECT_EQ(0u, s3.accepted);
}

TEST(PatternSearchState, WrappedPrimaryBaseNodes) {
    // on a wrapped PRIMARY graph the states name the wrapper's nodes, virtual ones included,
    // their stored k-mers and whether a stored k-mer is the reverse complement of the k-mer
    // spelled (the tracker checks all three by spelling; stored_reverse_complement always
    // false fails here)
    const std::vector<std::string> records { "ACGTTGCAAGTCAGG", "CCTGACTTGCAACGT" };
    auto graph = build_graph<DBGSuccinct>(5, records, DeBruijnGraph::PRIMARY);
    const std::string text = "ACGTTGCAAG";
    const Pattern pattern = Pattern::parse(PatternKind::DNA, text);
    const std::vector<Oriented> searched = orientations(PatternKind::DNA, text, Strands::BOTH);
    const WalkTree tree = walk_tree(*graph, searched);
    const Expected e = simulate(tree, Script {});
    ASSERT_FALSE(e.paths.empty());
    ScriptedTracker tracker(*graph, pattern, searched, Script {});
    Request request = path_request();
    request.support = &tracker;
    std::vector<PathCtx> got;
    Result r = enumerate_paths(*graph, pattern, request, &got);
    EXPECT_EQ(paths_of(tree, e.paths), got);
    EXPECT_TRUE(tracker.errors.empty()) << tracker.errors.front();
    EXPECT_EQ(tracker.opens + tracker.pushes, tracker.pops);
}

#endif // _DNA_GRAPH || _DNA5_GRAPH

} // namespace
