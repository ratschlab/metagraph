#include "walker.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <string_view>
#include <tuple>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/unix_tools.hpp"


namespace mtg {
namespace graph {
namespace traversal {

const char* to_string(EventType type) {
    switch (type) {
        case EventType::LABEL_END: return "label_end";
        case EventType::SWITCH: return "switch";
        case EventType::BLOCKED: return "blocked";
        case EventType::HAIRPIN: return "hairpin";
        case EventType::TIP: return "tip";
        case EventType::BUBBLE: return "bubble";
        case EventType::RECONVERGE: return "reconverge";
        case EventType::REVISIT: return "revisit";
    }
    return "unknown";
}

const char* to_string(LabelMode mode) {
    return mode == LabelMode::ANNOTATE ? "annotate" : "constrain";
}


namespace {

/*
 * Implementation notes (spec §6, increment I3):
 *
 * - Exploration is level-synchronous: all live path heads of an arm at extension
 *   depth d take one step in one level, the two arms alternate by depth (§6.8).
 *   Labels of every successor of a level are fetched in one LabelQuery batch; along
 *   structurally unbranched runs up to |batch_kmers| nodes are prefetched into the
 *   label cache and a structural lookahead cache only, which never changes which
 *   steps are taken. Lookahead work is not counted in the contract counters until
 *   the walk consumes it (successor_enumerations, keys_mapped, rows_requested are
 *   independent of batch_kmers; rows_fetched / cache_hits are not).
 * - Tip windows and bubble windows (Strategy::tip_window_bp / bubble_window_bp, the
 *   TIP / BUBBLE events and the GrowthBin::tips / bubbles counters) are NOT
 *   implemented in this increment; the fields are kept and stay zero (spec §12).
 * - EndReason has no values for `switched`, `superseded`, `minority`,
 *   `below_min_labels`, `split_limit`, `switch_sources` and `hairpin`: a source
 *   whose lineage only continues under another name ends its run with LABEL_LOST
 *   plus a SWITCH event; quorum and split-limit stops use BRANCH with the text in
 *   Event::text; a label whose only continuation is a skipped hairpin ends with
 *   DEAD_END, text "hairpin"; a label cut from the switch sources by
 *   max_switch_sources ends with LABEL_LOST, text "switch_sources".
 * - Every successor the walker refuses to a label whose lineage would have continued
 *   on it (quorum, split limit, branch limit, loss budget) is stated on the step's
 *   BranchEvent as a Refusal (§7.2): the tuned-run checker takes only that, or a
 *   structural block it can verify itself, as evidence for an omission — a missing
 *   child says nothing about why it is missing. Recording them is linear per
 *   re-minimisation round (ArmResult::refusal_scans), and a loss-budget refusal does
 *   not depend on whether the source goes on elsewhere (review round 4).
 * - Hairpins (§6.5): with skip_hairpins the self-RC step is inadmissible and gets a
 *   HAIRPIN event; without it the step is followed, flagged with a HAIRPIN event
 *   (text "followed") on the parent segment and never counts toward ambiguity
 *   (the split it causes is a divergence).
 * - Runs closed by a reconvergence merge (the lineage continues in the kept
 *   entry's run) have `ended == false` and `to_bp` set: they are not label ends.
 *   The merged segment records per parent which labels entered through it
 *   (Segment::labels_via_parent); an entry taken from a non-first parent carries
 *   the merge position as its route boundary so that continuations (spelled through
 *   first parents) only list labels covering the spelled tail. Under trace support
 *   the live coordinates of a label on both parents are united.
 * - Edge identity: 2-bit packed (k+1)-mer when it fits 64 bits, otherwise a 128-bit
 *   FNV-1a pair; two distinct (k+1)-mers colliding on the pair are treated as the
 *   same edge (conservative: blocks), the full string is not kept.
 * - Annotate mode (§6.9) shares the level loop, the structural admissibility
 *   (check_structure: hairpin, seed re-entry, per-path edge reuse), the caps, merging
 *   and the segment DAG with constrain mode, and nothing of the label state: items
 *   carry the labels PRESENT at their node (LabelRecorder, full rows) and
 *   process_item_annotate() follows every admissible successor. Leaves end with a
 *   path_reason and no label ends; segments carry the present sets as runs.
 * - Completeness boundary (§6.10): a cap trips between two heads of one level, so
 *   ArmState::boundary records the extension depth of the first head NOT expanded (or
 *   the depth of the first pruned head); complete_to_bp is its minimum with the
 *   requested radius. Heads that already reached the radius when the time budget or a
 *   seed-level cap trips are ended with max_extension_bp, not censored, so that an arm
 *   whose every walk is present is reported complete.
 */

constexpr char kSentinel = boss::BOSS::kSentinel;
// the structural lookahead cache is cleared when it grows beyond this many nodes
constexpr size_t kMaxLookahead = 1'000'000;

// deterministic tie-break key of a switch source (§6.3)
using SourceKey = std::tuple<double, uint32_t, Column, uint64_t>;

// One label on one path position (spec §6.3). |pred|, |switched|, |from| and
// |switch_cost| are filled by the recurrence and consumed when a step is committed.
struct Entry {
    LabelId label = 0;
    double loss = 0;
    uint32_t branches = 0;
    uint32_t run = UINT32_MAX;
    // the path spelled through first parents carries this lineage from here on
    // (set to the merge position when the entry came through another parent)
    uint64_t route_bp = 0;
    SmallVector<Coord> coords;      // trace support only: live coordinates at this node
    LabelId pred = 0;
    bool switched = false;
    LabelId from = 0;
    double switch_cost = 0;
};
// sorted by label id
using State = std::vector<Entry>;

// A live path head
struct Item {
    size_t segment = 0;
    node_index node = npos;
    std::string kmer;               // spelling of |node| (uppercase, natural orientation)
    uint64_t ext_bp = 0;
    State state;
    uint32_t splits = 0;
    size_t path_id = 0;
    bool revisiting = false;        // inside a stretch of nodes seen at another distance
    // annotate mode: the labels present at |node| (capped) and their true count
    std::vector<LabelId> present;
    size_t present_total = 0;
};

struct Succ {
    node_index node;
    char ch;
};

// A permitted label at a successor (A(v) of the spec)
struct Target {
    LabelId label;
    SmallVector<Coord> coords;
};

// A successor under evaluation (scratch, reused across steps)
struct Cand {
    const Succ *succ = nullptr;
    std::string kmer;
    std::vector<Target> targets;
    State state;
    bool hairpin = false;           // self-reverse-complementary step (canonical regimes)
    bool skipped = false;           // hairpin not followed (skip_hairpins)
    bool blocked = false;
    EndReason block_reason = EndReason::DEAD_END;
    uint64_t key_lo = 0;
    uint64_t key_hi = 0;
    bool rc = false;
    bool followed = false;
    bool quorum_fail = false;
    const char *quorum_text = "";
    size_t initial_labels = 0;      // |σ_v| before any source was excluded
    bool truncated = false;         // derive() cut the switch sources (max_switch_sources)
    SourceKey cut {};               // key of the last eligible source when truncated
    // the predecessors of the entries whose source the current re-minimisation round
    // excluded (refused "branch"; duplicates removed when the event is emitted)
    std::vector<LabelId> excluded_preds;
    LabelRecorder::NodeLabels present;   // annotate mode: the labels at the successor
    bool admissible() const { return !skipped && !blocked; }
    void reset() {
        succ = nullptr;
        kmer.clear();
        targets.clear();
        state.clear();
        hairpin = skipped = blocked = false;
        block_reason = EndReason::DEAD_END;
        key_lo = key_hi = 0;
        rc = followed = quorum_fail = false;
        quorum_text = "";
        initial_labels = 0;
        truncated = false;
        excluded_preds.clear();
        present.labels.clear();
        present.total = 0;
    }
};

struct EdgeKey {
    uint64_t lo;
    uint64_t hi;
    bool operator==(const EdgeKey &o) const { return lo == o.lo && hi == o.hi; }
};
struct EdgeKeyHash {
    size_t operator()(const EdgeKey &k) const {
        return k.lo ^ (k.hi * 0x9E3779B97F4A7C15ULL) ^ (k.hi >> 29);
    }
};
// one use of an edge: the segment that took it and the orientation bit
struct EdgeUse {
    uint64_t packed = 0;
    EdgeUse() = default;
    EdgeUse(size_t segment, bool rc) : packed((static_cast<uint64_t>(segment) << 1) | rc) {}
    size_t segment() const { return packed >> 1; }
    bool rc() const { return packed & 1; }
};
// the same use keyed the other way round, for probing "did THIS segment take the edge"
struct EdgeSegKey {
    uint64_t lo;
    uint64_t hi;
    size_t segment;
    bool operator==(const EdgeSegKey &o) const {
        return lo == o.lo && hi == o.hi && segment == o.segment;
    }
};
struct EdgeSegKeyHash {
    size_t operator()(const EdgeSegKey &k) const {
        return EdgeKeyHash()(EdgeKey{ k.lo, k.hi }) ^ (k.segment * 0xC2B2AE3D27D4EB4FULL);
    }
};

// structural lookahead: the successors of a node and, when it has exactly one,
// the annotation key of that successor
struct Lookahead {
    std::vector<Succ> succs;
    node_index succ_key = npos;
};

struct LeafInfo {
    bool is_leaf = false;
    std::array<uint32_t, kNumEndReasons> end_reasons {};
    std::optional<EndReason> path_reason;
    std::vector<LabelEnd> end_labels;
    std::optional<Continuation> continuation;
};

struct ArmState {
    Arm arm = Arm::RIGHT;
    ArmResult result;
    std::vector<Item> frontier;
    std::vector<Item> next;
    std::vector<std::string> walk_seq;      // per segment, in walking order
    std::vector<LeafInfo> leaves;           // per segment
    // The per-path edge-reuse record, kept both ways: every use of an edge (scanned
    // when the edge has few uses) and, per (edge, segment), the orientation bits
    // 1 << rc of that segment's uses (probed per ancestor when it has many — in keep
    // mode every live path through a shared region records the same edge, and
    // scanning those uses once per path is quadratic in the number of paths).
    tsl::hopscotch_map<EdgeKey, SmallVector<EdgeUse>, EdgeKeyHash> used_edges;
    tsl::hopscotch_map<EdgeSegKey, uint8_t, EdgeSegKeyHash> used_by_segment;
    tsl::hopscotch_map<node_index, std::pair<size_t, uint64_t>> first_arrival;
    tsl::hopscotch_map<node_index, Lookahead> lookahead;
    size_t finished_leaves = 0;
    size_t next_path_id = 1;
    bool stopped = false;
    // extension depth of the first head not expanded (or first pruned): the level it
    // belongs to is partial and does not count toward complete_to_bp
    uint64_t boundary = std::numeric_limits<uint64_t>::max();
    // ancestor marks of |marked_segment| (ancestors never change after creation), and
    // the marked segments as a list: the proper ancestors, in no particular order
    std::vector<uint32_t> visit_mark;
    uint32_t visit_epoch = 0;
    size_t marked_segment = SIZE_MAX;
    std::vector<size_t> visit_stack;
    std::vector<size_t> ancestors;
};

// Per-label scratch of process_item, sized |label_dict| once and cleared over the
// touched labels only, so that a step costs O(|σ| + |A(v)|) and not O(|label_dict|).
// |marked| is set and cleared within one pass (the sources excluded in one
// re-minimisation round, the predecessors on one successor); |budget_checked| marks the
// sources whose loss-budget refusals were recorded with their label end.
struct Scratch {
    std::vector<uint8_t> trace_broken, excluded, seen, has_cont, stays, only_hairpin,
                         in_quorum_fail, superseded, blocked, is_touched, marked,
                         budget_checked;
    std::vector<uint32_t> cont_count;
    std::vector<const char*> qtext;
    std::vector<LabelId> touched;

    void init(size_t n) {
        for (auto *v : { &trace_broken, &excluded, &seen, &has_cont, &stays, &only_hairpin,
                         &in_quorum_fail, &superseded, &blocked, &is_touched, &marked,
                         &budget_checked }) {
            v->assign(n, 0);
        }
        cont_count.assign(n, 0);
        qtext.assign(n, "");
        touched.clear();
    }
    void touch(LabelId l) {
        if (!is_touched[l]) {
            is_touched[l] = 1;
            touched.push_back(l);
        }
    }
    void reset() {
        for (LabelId l : touched) {
            trace_broken[l] = excluded[l] = seen[l] = has_cont[l] = stays[l] = only_hairpin[l]
                = in_quorum_fail[l] = superseded[l] = blocked[l] = is_touched[l] = marked[l]
                = budget_checked[l] = 0;
            cont_count[l] = 0;
            qtext[l] = "";
        }
        touched.clear();
    }
};
struct ScratchGuard {
    Scratch &scratch;
    ~ScratchGuard() { scratch.reset(); }
};

uint8_t block_rank(EndReason reason) {
    switch (reason) {
        case EndReason::REACHED_SEED: return 3;
        case EndReason::EDGE_REUSE_RC: return 2;
        case EndReason::EDGE_REUSE: return 1;
        default: return 0;
    }
}

EndReason block_reason_of_rank(uint8_t rank) {
    switch (rank) {
        case 3: return EndReason::REACHED_SEED;
        case 2: return EndReason::EDGE_REUSE_RC;
        default: return EndReason::EDGE_REUSE;
    }
}

std::vector<LabelId> labels_of(const State &state) {
    std::vector<LabelId> out;
    out.reserve(state.size());
    for (const Entry &e : state) {
        out.push_back(e.label);
    }
    return out;
}

const Entry* find_entry(const State &state, LabelId label) {
    auto it = std::lower_bound(state.begin(), state.end(), label,
                               [](const Entry &e, LabelId l) { return e.label < l; });
    return it != state.end() && it->label == label ? &*it : nullptr;
}

Entry* find_entry(State &state, LabelId label) {
    return const_cast<Entry*>(find_entry(static_cast<const State&>(state), label));
}

double min_loss(const State &state) {
    double m = kInfiniteLoss;
    for (const Entry &e : state) {
        m = std::min(m, e.loss);
    }
    return m;
}

// sorted set union of live coordinates
void unite_coords(SmallVector<Coord> *into, const SmallVector<Coord> &other) {
    if (other.empty())
        return;
    SmallVector<Coord> out;
    std::set_union(into->begin(), into->end(), other.begin(), other.end(),
                   std::back_inserter(out));
    *into = std::move(out);
}


class Walker {
  public:
    Walker(LabelOracle &oracle,
           const Seed &seed,
           const Strategy &strategy,
           const LabelChangeCost &cost,
           const std::string &release_id)
          : oracle_(oracle), graph_(oracle.graph()), seed_(seed), strategy_(strategy),
            cost_(cost), release_id_(release_id), k_(oracle.get_k()),
            regime_(oracle.regime()), canonical_(oracle.canonical()),
            nfc_(oracle.node_first_cache()),
            trace_(strategy.support == Support::TRACE),
            annotate_(strategy.label_mode == LabelMode::ANNOTATE) {}

    SeedResult run();

  private:
    // ---- setup
    void validate_seed();
    // Permitted set derived from the seed itself (Seed::labels empty): the labels
    // supporting every seed k-mer, with the per-k-mer hits of that same single pass.
    // Returns true when every derived label supports every seed k-mer by construction
    // and |hits| was therefore left EMPTY (the non-trace case: materialising it would
    // be L*M Hit objects with no information in them). Under `support: trace` it fills
    // |hits| with the coordinates found in this pass and returns false.
    bool derive_seed_labels(const std::vector<node_index> &keys,
                            std::vector<LabelRef> *refs,
                            std::vector<LabelQuery::NodeHits> *hits);
    void init_edge_coding();
    void init_arm(ArmState &arm);

    // ---- graph access
    std::vector<Succ> enumerate(ArmState &arm, node_index node, const std::string &kmer);
    std::vector<Succ> successors(ArmState &arm, node_index node, const std::string &kmer,
                                 bool *cached_key, node_index *key);
    node_index key_of_succ(Arm arm, const std::string &kmer, const Succ &s);
    void succ_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const;
    void step_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const;
    bool is_hairpin(const std::string &kmer_u, const std::string &kmer_v,
                    const std::string &step);
    std::pair<EdgeKey, bool> edge_key(const std::string &step);

    // ---- levels
    void run_level(ArmState &arm, uint64_t depth);
    void sort_items(std::vector<Item> &items) const;
    std::optional<EndReason> process_item(ArmState &arm, Item &item,
                                          const std::vector<Succ> &succs,
                                          const LabelQuery::NodeHits *hits,
                                          size_t remaining_in_level);
    // annotate mode: every structurally admissible successor is followed
    std::optional<EndReason> process_item_annotate(ArmState &arm, Item &item,
                                                   const std::vector<Succ> &succs,
                                                   const LabelRecorder::NodeLabels *present,
                                                   size_t remaining_in_level);
    // structural admissibility of the step item -> c.succ (both modes): hairpin, seed
    // re-entry, per-path edge reuse; requires c.kmer
    void check_structure(ArmState &arm, const Item &item, Cand &c);
    // the caps, decided before anything is committed; nf = successors to follow
    std::optional<EndReason> cap_check(const ArmState &arm, size_t nf,
                                       size_t remaining_in_level) const;
    void account_steps(ArmState &arm, uint64_t at, size_t nf);
    void record_edge(ArmState &arm, size_t segment, const Cand &c);
    void prefetch(ArmState &arm, const std::vector<Item> &items,
                  const std::vector<std::vector<Succ>> &succs);
    void arrival(ArmState &arm, Item &item);
    void merge_level(ArmState &arm, uint64_t depth);
    void beam(ArmState &arm, uint64_t depth);
    void trip(ArmState &arm, EndReason reason, std::vector<Item> &items, size_t index);
    void stop_arm(ArmState &arm, EndReason reason, std::vector<Item> *items, size_t from);
    // end the frontier of |arm| with |reason|; heads that already reached the radius
    // are complete (max_extension_bp), not censored
    void stop_frontier(ArmState &arm, EndReason reason);
    bool time_exceeded() const;

    // ---- recorded labels (annotate mode) and bounded label lists (both modes)
    std::vector<LabelId> labels_at(const Item &item) const {
        return annotate_ ? item.present : labels_of(item.state);
    }
    // the first max_labels_per_node ids of a sorted list, counting the cut
    std::vector<LabelId> bounded(ArmState &arm, const std::vector<LabelId> &labels, size_t total);
    // append the labels present at the node entered by the step at |at| to the
    // segment's runs
    void record_present(ArmState &arm, size_t segment, uint64_t at,
                        const LabelRecorder::NodeLabels &nl);
    void summarize_annotate();

    // ---- label state
    void derive(const State &sigma, const std::vector<Target> &targets,
                const std::vector<uint8_t> &excluded, State *out,
                bool *truncated, SourceKey *cut);
    SourceKey source_key(const Entry &e) const {
        const LabelRef &ref = result_.label_dict[e.label];
        return { e.loss, e.branches, ref.column, ref.seq_id };
    }
    // |taken|: runs already continued by an earlier child of the same split; a
    // run taken twice is cloned so that every run belongs to exactly one path
    void commit_entries(ArmState &arm, const Item &item, size_t target_segment,
                        State &state, uint64_t at, std::vector<uint32_t> *taken);

    // ---- ends
    void end_run(ArmState &arm, const Entry &e, uint64_t at, EndReason reason);
    void end_label(ArmState &arm, const Item &item, const Entry &e, EndReason reason,
                   size_t structural, const char *text, double needed = 0);
    void finish_path(ArmState &arm, Item &item, std::optional<EndReason> path_reason);
    void censor_item(ArmState &arm, Item &item, EndReason reason);
    Continuation make_continuation(ArmState &arm, const Item &item);

    // ---- bookkeeping
    GrowthBin& bin(ArmState &arm, uint64_t bp);
    size_t new_segment(ArmState &arm, std::vector<size_t> parents, uint64_t from_bp,
                       std::vector<LabelId> labels_start, size_t labels_start_total);
    uint32_t new_run(ArmState &arm, LabelId label, uint64_t from_bp, bool by_switch,
                     LabelId from, double cost, uint32_t prev);
    void mark_ancestors(ArmState &arm, size_t seg);
    bool is_ancestor_or_self(const ArmState &arm, size_t anc, size_t seg) const;
    void record_live(ArmState &arm, uint64_t bp, const std::vector<Item> &items);
    // distinct labels on the heads a[from_a..] and b below extension depth |below_bp|;
    // |*exact| is cleared when a head carries a list cut by max_labels_per_node
    // (annotate mode), so the count is a lower bound
    size_t distinct_labels(const std::vector<Item> &a, const std::vector<Item> &b,
                           size_t from_a, bool *exact,
                           uint64_t below_bp = std::numeric_limits<uint64_t>::max());
    void finalize(ArmState &arm);
    void summarize();

    LabelOracle &oracle_;
    const DeBruijnGraph &graph_;
    const Seed &seed_;
    const Strategy &strategy_;
    LabelChangeCost cost_;          // remapped to dictionary ids by validate_seed()
    const std::string &release_id_;
    const size_t k_;
    const Regime regime_;
    const CanonicalDBG *canonical_;
    const NodeFirstCache *nfc_;
    const bool trace_;
    const bool annotate_;

    Timer timer_;
    SeedResult result_;
    std::string seed_upper_;
    std::vector<node_index> nodes_;
    tsl::hopscotch_set<node_index> seed_nodes_;
    std::unique_ptr<LabelQuery> query_;         // constrain mode
    std::unique_ptr<LabelRecorder> recorder_;   // annotate mode
    // trace support: live coordinates of each kept seed label at the arm boundaries
    std::vector<SmallVector<Coord>> boundary_coords_[2];

    std::array<uint8_t, 256> code_ {};
    size_t bits_ = 2;
    bool packable_ = true;

    std::array<ArmState, 2> arms_;
    uint64_t steps_total_ = 0;
    uint64_t pair_evaluations_ = 0;
    bool seed_stopped_ = false;

    // scratch reused across steps
    Scratch scratch_;
    std::vector<Cand> cands_;
    std::vector<const Entry*> sw_;
    std::string step_;
    std::string rc_scratch_;
    std::vector<uint32_t> label_stamp_;
    uint32_t label_epoch_ = 0;
};


/********************************** setup ***********************************/

void Walker::validate_seed() {
    const std::string &seq = seed_.sequence;
    if (seq.size() < k_)
        throw std::invalid_argument("Seed shorter than k");

    seed_upper_ = seq;
#if ! _DNA_CASE_SENSITIVE_GRAPH
    for (char &c : seed_upper_) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
#endif
    const std::string &alphabet = graph_.alphabet();
    for (size_t i = 0; i < seed_upper_.size(); ++i) {
        char c = seed_upper_[i];
        if (c == kSentinel || alphabet.find(c) == std::string::npos) {
            throw std::invalid_argument("Invalid character '" + std::string(1, seq[i])
                                        + "' at position " + std::to_string(i) + " of the seed");
        }
    }

    nodes_ = map_to_nodes_sequentially(graph_, seed_upper_);
    assert(nodes_.size() == seq.size() - k_ + 1);
    {
        std::vector<KmerInterval> runs;
        bool missing = false;
        for (uint64_t i = 0; i < nodes_.size(); ++i) {
            if (nodes_[i] == npos) {
                missing = true;
                continue;
            }
            if (runs.empty() || runs.back().end != i) {
                runs.push_back({ i, i + 1 });
            } else {
                runs.back().end = i + 1;
            }
        }
        if (missing) {
            throw std::invalid_argument("Seed is not fully present in the graph; graph runs: "
                                        + encode_runs(runs, nodes_.size()));
        }
    }
    result_.seed_id = seed_.seed_id;
    result_.length_bp = seq.size();
    result_.num_kmers = nodes_.size();

    // ---- seed nodes (either orientation in canonical regimes)
    for (node_index n : nodes_) {
        seed_nodes_.insert(n);
    }
    if (regime_ != Regime::BASIC) {
        std::string rc = seed_upper_;
        ::reverse_complement(rc);
        for (node_index n : map_to_nodes_sequentially(graph_, rc)) {
            if (n != npos)
                seed_nodes_.insert(n);
        }
    }

    // ---- annotate mode: no permitted set, nothing to validate the seed against. The
    // dictionary is filled by the recorder as labels are met; the seed id is the one of
    // the sequence under no labels.
    if (annotate_) {
        if (!seed_.labels.empty()) {
            throw std::invalid_argument("labels.mode \"annotate\" records the labels present and "
                                        "filters by none: omit seeds[].labels ("
                                        + std::to_string(seed_.labels.size()) + " given)");
        }
        result_.num_seed_labels = 0;
        result_.labels_supporting_total = 0;
        result_.validated_seed_id = make_seed_id(release_id_, seed_upper_,
                                                 regime_ != Regime::BASIC, {});
        result_.seed_id_mismatch = !seed_.seed_id.empty()
                                    && seed_.seed_id != result_.validated_seed_id;
        return;
    }

    // ---- labels
    if (trace_) {
        if (!oracle_.has_coordinates())
            throw std::invalid_argument("Trace support requires an annotation with k-mer coordinates");
        if (regime_ != Regime::BASIC) {
            throw std::invalid_argument("Trace support is only defined for BASIC (forward-strand) "
                                        "graphs: coordinates carry no strand in canonical indexes");
        }
        if (strategy_.merge_reconverge) {
            // A merge unites the coordinate sets of both parents but keeps only one
            // parent's (loss, branches, run, route). A later step may then continue on
            // the discarded parent's coordinates while carrying the retained parent's
            // evidence, which can report direct support, a loss and a run that no single
            // occurrence justifies. Doing this correctly needs per-coordinate
            // provenance, so refuse the combination rather than weaken the evidence.
            throw std::invalid_argument("Trace support cannot be combined with reconvergence "
                                        "merging: set branching.on_reconverge to \"keep\" "
                                        "(merging keeps only one parent's evidence per label, "
                                        "which a coordinate trace cannot justify)");
        }
    }

    // support of every seed k-mer. One pass over the seed's annotation rows: either
    // of the labels given, or — when none are given — of every label, whose
    // intersection IS the derived permitted set (§6.1, `permit: per_hit`).
    std::vector<node_index> keys = oracle_.keys_of_sequence(seed_upper_);
    std::vector<LabelRef> seed_refs;
    std::vector<LabelQuery::NodeHits> hits;
    // A set derived from the seed supports every seed k-mer by construction, so unless
    // `trace` still has to find the coordinates consecutive there is nothing for the
    // per-label loop below to look up — and no reason to build an L x M hit matrix.
    bool all_supported = false;
    if (seed_.labels.empty()) {
        all_supported = derive_seed_labels(keys, &seed_refs, &hits);
    } else {
        for (const auto &name : seed_.labels) {
            for (const auto &other : seed_refs) {
                if (other.name == name)
                    throw std::invalid_argument("Duplicate seed label '" + name + "'");
            }
            seed_refs.push_back(oracle_.resolve_label(name));
        }
        LabelQuery validation(oracle_, seed_refs, trace_);
        hits = validation.fetch(keys);
    }

    // request index (seed labels, then extra) -> dictionary id, UINT32_MAX if dropped
    std::vector<LabelId> dict_of_request(seed_refs.size() + strategy_.extra.size(), UINT32_MAX);
    std::vector<std::string> kept_names;
    // Per-label support, gathered in ONE sweep over the k-mers that scatters every hit
    // into its label's run list — O(hits), the size of the data, the shape resolve.cpp
    // uses — instead of looking each label up at each k-mer (which was O(L*M*log L)
    // and, before that, O(L^2*M)). Under trace support the live coordinate set of
    // every label advances in that same sweep. Nothing L x M is materialised: a
    // derived non-trace set left |hits| empty on purpose (all_supported).
    const size_t num_refs = seed_refs.size();
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> runs(num_refs);
    std::vector<uint64_t> supported_kmers(num_refs, 0);
    std::vector<std::vector<Coord>> live(trace_ ? num_refs : 0);
    std::vector<uint8_t> chain_ok(num_refs, 1);
    std::vector<Coord> next_live;
    for (size_t i = 0; !all_supported && i < keys.size(); ++i) {
        for (const LabelQuery::Hit &h : hits[i]) {
            const LabelId l = h.label;
            supported_kmers[l]++;
            if (runs[l].empty() || runs[l].back().second != i) {
                runs[l].emplace_back(i, i + 1);
            } else {
                runs[l].back().second = i + 1;
            }
            if (!trace_)
                continue;
            std::vector<Coord> &lv = live[l];
            if (i == 0) {
                lv.assign(h.coords.begin(), h.coords.end());
            } else {
                next_live.clear();
                for (Coord c : h.coords) {
                    if (c > 0 && std::binary_search(lv.begin(), lv.end(), c - 1))
                        next_live.push_back(c);
                }
                lv.swap(next_live);
            }
            if (lv.empty())
                chain_ok[l] = 0;
        }
    }
    for (LabelId l = 0; l < num_refs; ++l) {
        const bool all = all_supported || supported_kmers[l] == keys.size();
        if (all && chain_ok[l]) {
            dict_of_request[l] = result_.label_dict.size();
            result_.label_dict.push_back(seed_refs[l]);
            kept_names.push_back(seed_refs[l].name);
            if (trace_) {
                SmallVector<Coord> right(live[l].begin(), live[l].end());
                SmallVector<Coord> left;
                for (Coord c : live[l]) {
                    left.push_back(c - (keys.size() - 1));
                }
                boundary_coords_[static_cast<size_t>(Arm::RIGHT)].push_back(right);
                boundary_coords_[static_cast<size_t>(Arm::LEFT)].push_back(left);
            }
        } else {
            DroppedLabel dropped;
            dropped.name = seed_refs[l].name;
            dropped.reason = "seed_unsupported";
            dropped.runs = std::move(runs[l]);
            result_.dropped_labels.push_back(std::move(dropped));
        }
    }
    if (result_.label_dict.empty()) {
        // under `support: trace` a derived label can still be dropped here: the set is
        // derived from k-mer presence and must then also be coordinate-consecutive
        if (result_.labels_from_seed)
            throw SeedDerivationError("No label derived from the seed supports every "
                                      "k-mer of the seed");
        throw std::invalid_argument("No seed label supports every k-mer of the seed");
    }
    result_.num_seed_labels = result_.label_dict.size();
    if (!result_.labels_from_seed) {
        // An explicit list is its own support total. Emitting 0 here would make the one
        // field a client uses to detect truncation
        // (labels_supporting_total > |labels|) read as "nothing supports this seed".
        result_.labels_supporting_total = result_.num_seed_labels;
    }
    // the seed id is defined over the case-mapped sequence (§4.3 / §6.1)
    result_.validated_seed_id = make_seed_id(release_id_, seed_upper_,
                                             regime_ != Regime::BASIC, kept_names);
    result_.seed_id_mismatch = !seed_.seed_id.empty()
                                && seed_.seed_id != result_.validated_seed_id;

    // ---- extra labels (switch targets)
    for (size_t i = 0; i < strategy_.extra.size(); ++i) {
        LabelRef ref = oracle_.resolve_label(strategy_.extra[i]);
        for (const auto &other : result_.label_dict) {
            if (other.same_target(ref)) {
                throw std::invalid_argument("Extra label '" + strategy_.extra[i]
                                            + "' duplicates a seed label");
            }
        }
        dict_of_request[seed_refs.size() + i] = result_.label_dict.size();
        result_.label_dict.push_back(ref);
    }

    // ---- the cost table is authored over request indices: remap it to dictionary
    // ids, dropping entries that name a dropped seed label
    if (cost_.model() == LabelChangeCost::TABLE) {
        std::map<std::pair<LabelId, LabelId>, double> entries;
        for (const auto &[pair, value] : cost_.table()) {
            if (pair.first >= dict_of_request.size() || pair.second >= dict_of_request.size())
                continue;
            LabelId a = dict_of_request[pair.first], b = dict_of_request[pair.second];
            if (a == UINT32_MAX || b == UINT32_MAX)
                continue;
            entries[{ a, b }] = value;
        }
        cost_ = LabelChangeCost::table(std::move(entries), cost_.default_cost());
    }
    for (LabelId id = result_.num_seed_labels; id < result_.label_dict.size(); ++id) {
        bool reachable = false;
        if (cost_.finite()) {
            for (LabelId s = 0; s < result_.num_seed_labels; ++s) {
                reachable |= cost_.cost(s, id) <= strategy_.loss_budget;
            }
        }
        if (!reachable) {
            throw std::invalid_argument("Extra label '" + result_.label_dict[id].name
                                        + "' is unreachable: no seed label can switch to it "
                                        "within the loss budget");
        }
    }
}

bool Walker::derive_seed_labels(const std::vector<node_index> &keys,
                                std::vector<LabelRef> *refs,
                                std::vector<LabelQuery::NodeHits> *hits) {
    // (column, seq_id); seq_id is 0 and meaningless for COLUMN labels
    using Key = std::pair<Column, uint64_t>;

    const LabelKind kind = strategy_.seed_label_kind.value_or(
            oracle_.coord_to_header() ? LabelKind::HEADER : LabelKind::COLUMN);
    if (kind == LabelKind::HEADER && !oracle_.coord_to_header()) {
        throw std::invalid_argument("Deriving sequence header labels from the seed requires a "
                                    "CoordToHeader index; set seed_label_kind to \"column\"");
    }
    // HEADER needs coordinates to tell the indexed sequences of a column apart, trace
    // needs them as evidence; a plain column set needs neither
    const bool with_coords = kind == LabelKind::HEADER || trace_;
    if (with_coords && !oracle_.has_coordinates()) {
        throw std::invalid_argument("Deriving the permitted set from the seed needs an annotation "
                                    "with k-mer coordinates for this label kind");
    }
    if (!strategy_.max_seed_labels)
        throw std::invalid_argument("max_seed_labels must be positive");

    // The derivation is the most expensive phase of a derived seed — one FULL annotation
    // row per seed k-mer, where an explicit list reads only its own columns — and it runs
    // entirely before the walk loop, i.e. before the only other place that looks at the
    // clock. So it watches the budget itself. A non-positive budget means "no extension"
    // (time_exceeded() treats it as already spent), not "no derivation", so it is not an
    // immediate deadline here; the server's own cap is what bounds that case.
    const double budget_ms = strategy_.time_budget_ms;
    auto out_of_time = [&]() {
        return budget_ms > 0 && timer_.elapsed() * 1000.0 >= budget_ms;
    };
    // The candidate set of the FIRST k-mer consumed is a whole annotation row: no
    // intersection has narrowed it yet. A seed of exactly k bases never gets an
    // intersection at all, so its "derived" set would BE that row — on a wide index
    // millions of entries that a cap applied afterwards cannot un-materialise. Refuse
    // such a seed up front instead of building it.
    const size_t max_candidates = std::max<size_t>(
            size_t(1) << 16,
            64 * std::min<size_t>(strategy_.max_seed_labels, size_t(1) << 40));

    // the running intersection, ascending by (column, seq_id), and its distinct columns
    std::vector<Key> live, next;
    std::vector<Column> live_columns;
    // trace support only: the live coordinates of the candidates at every seed k-mer,
    // kept from this one pass so that the per-label validation needs no second read
    std::vector<std::vector<std::pair<Key, SmallVector<Coord>>>> coords_at;
    if (trace_)
        coords_at.resize(keys.size());
    std::vector<size_t> done;         // the k-mers already consumed, in that order
    done.reserve(keys.size());
    size_t compacted_at = 0;          // |live| when |coords_at| was last compacted
    std::vector<std::pair<Key, SmallVector<Coord>>> cur;
    std::vector<std::pair<Key, Coord>> flat;
    auto by_key = [](const std::pair<Key, SmallVector<Coord>> &a,
                     const std::pair<Key, SmallVector<Coord>> &b) { return a.first < b.first; };
    // after the first k-mer only the columns still carrying a candidate matter, so the
    // per-k-mer work shrinks with the intersection instead of scanning the whole row
    auto skip = [&](bool first, Column c) {
        return !first && !std::binary_search(live_columns.begin(), live_columns.end(), c);
    };

    // The derivation fetches in sub-batches of at most this many k-mers, whatever
    // |batch_kmers| the caller asked for: batching the row reconstruction is what makes a
    // long seed fast, but a sub-batch is reconstructed BEFORE the intersection or the
    // clock can stop it, so the work that can be wasted has to stay bounded.
    // The FIRST sub-batch is always kMaxChunk wide, whatever batch_kmers says: it is the
    // window in which the cheapest row is chosen and the candidate guard below is
    // applied, and whether a seed is accepted must not depend on a fetch-size knob.
    constexpr size_t kMaxChunk = 64;
    const size_t chunk = std::clamp<size_t>(strategy_.batch_kmers, 1, kMaxChunk);
    for (size_t begin = 0, width = kMaxChunk; begin < keys.size(); begin += width, width = chunk) {
        const size_t end = std::min(keys.size(), begin + width);
        std::vector<Row> rows;
        rows.reserve(end - begin);
        for (size_t i = begin; i < end; ++i) {
            assert(keys[i] != npos);   // the seed is fully present in the graph
            rows.push_back(AnnotatedDBG::graph_to_anno_index(keys[i]));
        }
        std::vector<Row> distinct = rows;
        std::sort(distinct.begin(), distinct.end());
        distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
        // THE single read of the seed's rows: the derived set is the intersection of
        // the same rows the per-label seed validation consumes below
        std::vector<annot::matrix::BinaryMatrix::SetBitPositions> plain;
        std::vector<annot::matrix::MultiIntMatrix::RowTuples> tuples;
        if (with_coords) {
            tuples = oracle_.get_row_tuples(distinct);
        } else {
            plain = oracle_.get_rows(distinct);
        }
        auto row_of = [&](size_t i) {
            return static_cast<size_t>(
                    std::lower_bound(distinct.begin(), distinct.end(), rows[i - begin])
                    - distinct.begin());
        };

        // In which order the sub-batch is consumed. The k-mer that SEEDS the intersection
        // sets the peak work and the peak memory of the whole derivation, and k-mer 0 is
        // merely where the caller cut the seed: a seed starting in a conserved or
        // repetitive k-mer costs orders of magnitude more than the same seed shifted by
        // two bases. So the first sub-batch starts from its cheapest row instead (its
        // sizes are already known — the batch is fetched before anything is intersected).
        // An intersection is commutative and |live| stays ascending by key, so neither
        // the derived set nor its order depends on this choice.
        std::vector<size_t> order;
        order.reserve(end - begin);
        for (size_t i = begin; i < end; ++i) {
            order.push_back(i);
        }
        if (done.empty()) {
            auto cost_of = [&](size_t i) {
                const size_t r = row_of(i);
                if (!with_coords)
                    return plain[r].size();
                size_t n = 0;    // columns plus coordinates: every one costs a map_coord
                for (const auto &entry : tuples[r]) {
                    n += 1 + entry.second.size();
                }
                return n;
            };
            std::vector<size_t> cost(order.size());
            for (size_t j = 0; j < order.size(); ++j) {
                cost[j] = cost_of(order[j]);
            }
            const size_t best = static_cast<size_t>(
                    std::min_element(cost.begin(), cost.end()) - cost.begin());
            std::iter_swap(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(best));
            if (cost[best] > max_candidates) {
                throw SeedDerivationError(
                        "The permitted set cannot be derived from this seed: the narrowest of "
                        "its first " + std::to_string(end - begin) + " k-mers alone has "
                        + std::to_string(cost[best]) + " annotation entries (columns plus "
                        "k-mer coordinates, an upper bound on its labels), over the limit of "
                        + std::to_string(max_candidates) + ". Start the seed in a less "
                          "repetitive k-mer, or name the labels explicitly.");
            }
        }

        for (size_t i : order) {
            const bool first = done.empty();
            // counted when the row is CONSUMED, not when the sub-batch is fetched, so
            // the reported count is the number of k-mers the derivation actually read
            oracle_.counters().rows_requested++;
            const size_t r = row_of(i);
            cur.clear();
            if (!with_coords) {
                for (Column c : plain[r]) {
                    if (!skip(first, c))
                        cur.emplace_back(Key{ c, 0 }, SmallVector<Coord>());
                }
                std::sort(cur.begin(), cur.end(), by_key);
            } else if (kind == LabelKind::COLUMN) {
                for (const auto &[c, coords] : tuples[r]) {
                    if (skip(first, c))
                        continue;
                    SmallVector<Coord> cc(coords.begin(), coords.end());
                    std::sort(cc.begin(), cc.end());
                    cur.emplace_back(Key{ c, 0 }, std::move(cc));
                }
                std::sort(cur.begin(), cur.end(), by_key);
            } else {
                flat.clear();
                for (const auto &[c, coords] : tuples[r]) {
                    if (skip(first, c))
                        continue;
                    for (Coord coord : coords) {
                        auto [seq_id, local] = oracle_.map_coord(c, coord);
                        flat.emplace_back(Key{ c, seq_id }, local);
                    }
                }
                std::sort(flat.begin(), flat.end());
                for (const auto &[key, local] : flat) {
                    if (cur.empty() || cur.back().first != key)
                        cur.emplace_back(key, SmallVector<Coord>());
                    // presence alone decides the derived set; the coordinates are only
                    // kept when `trace` has to find them consecutive
                    if (trace_ && (cur.back().second.empty() || cur.back().second.back() != local))
                        cur.back().second.push_back(local);
                }
            }

            next.clear();
            if (trace_)
                coords_at[i].clear();
            if (first) {
                for (auto &e : cur) {
                    next.push_back(e.first);
                    if (trace_)
                        coords_at[i].emplace_back(e.first, std::move(e.second));
                }
            } else {
                for (size_t a = 0, b = 0; a < live.size() && b < cur.size(); ) {
                    if (live[a] < cur[b].first) {
                        ++a;
                    } else if (cur[b].first < live[a]) {
                        ++b;
                    } else {
                        next.push_back(live[a]);
                        if (trace_)
                            coords_at[i].emplace_back(live[a], std::move(cur[b].second));
                        ++a;
                        ++b;
                    }
                }
            }
            live.swap(next);
            done.push_back(i);
            // bail out as soon as nothing can support the whole seed. The sub-batch
            // this k-mer came from was already reconstructed, so the waste is bounded by
            // kMaxChunk rows, not by the rest of the seed; |rows_requested| above counts
            // only the rows actually consumed, so it does not claim otherwise.
            if (live.empty()) {
                throw SeedDerivationError("No label supports every k-mer of the seed (the "
                                          "permitted set derived from the seed is empty at "
                                          "k-mer " + std::to_string(i) + " of "
                                          + std::to_string(keys.size()) + ")");
            }
            live_columns.clear();
            for (const Key &key : live) {
                if (live_columns.empty() || live_columns.back() != key.first)
                    live_columns.push_back(key.first);
            }
            if (first) {
                compacted_at = live.size();
            } else if (trace_ && live.size() * 2 <= compacted_at) {
                // |coords_at| is only ever read for the keys that survive to the end, so
                // entries for keys the intersection has already dropped are dead weight —
                // and the first k-mer's entry is a whole row of heap-allocated coordinate
                // vectors. Compacting on every shrink would be quadratic, so compact once
                // the live set has halved: O(log) passes, and the retention stays within a
                // factor of two of the keys that can still matter.
                for (size_t j : done) {
                    auto &v = coords_at[j];
                    size_t w = 0;
                    for (size_t a = 0, b = 0; a < v.size() && b < live.size(); ) {
                        if (v[a].first < live[b]) {
                            ++a;
                        } else if (live[b] < v[a].first) {
                            ++b;
                        } else {
                            if (w != a)
                                v[w] = std::move(v[a]);
                            ++w;
                            ++a;
                            ++b;
                        }
                    }
                    v.resize(w);
                    v.shrink_to_fit();
                }
                compacted_at = live.size();
            }
            if (out_of_time()) {
                throw SeedDerivationError(
                        "The time budget (bounds.time_budget_ms) ran out while deriving the "
                        "permitted set from the seed, after " + std::to_string(done.size())
                        + " of " + std::to_string(keys.size()) + " k-mers; name the labels "
                          "explicitly or shorten the seed");
            }
        }
    }

    result_.labels_from_seed = true;
    result_.labels_supporting_total = live.size();
    auto name_of = [&](const Key &key) -> const std::string& {
        return kind == LabelKind::COLUMN ? oracle_.column_name(key.first)
                                         : oracle_.header_name(key.first, key.second);
    };
    if (kind == LabelKind::HEADER) {
        // A FASTA header is unique within a column, not within the index, and the set is
        // deduplicated by (column, seq_id): two columns holding the same accession both
        // survive. Such a list is not resubmittable — the explicit path rejects duplicate
        // names and a name resolves to the FIRST column holding it — and every name-keyed
        // output (label_dict, label_summary, end_labels, the seed id) would be ambiguous.
        // Refuse it rather than echo a list the caller cannot use.
        std::vector<const std::string*> names;
        names.reserve(live.size());
        for (const Key &key : live) {
            names.push_back(&name_of(key));
        }
        std::sort(names.begin(), names.end(),
                  [](const std::string *a, const std::string *b) { return *a < *b; });
        for (size_t i = 1; i < names.size(); ++i) {
            if (*names[i - 1] == *names[i]) {
                throw SeedDerivationError(
                        "The labels derived from the seed are ambiguous: the sequence header '"
                        + *names[i] + "' occurs in more than one annotation column, so the "
                        "derived list cannot be resubmitted as explicit labels. Set "
                        "seed_label_kind to \"column\" or name the labels explicitly.");
            }
        }
        // ... and a header that an explicit list would resolve to something else is just
        // as unusable. Every derived name must round-trip through the resolver the
        // explicit path uses, resolve_label(): it tries column names FIRST, so a header
        // spelled like a column resolves to that column (whatever find_header() would say:
        // review round 3, finding 3), and a header held by several columns resolves to
        // the first of them, which need not be the derived one.
        for (const Key &key : live) {
            const std::string &name = name_of(key);
            const LabelRef back = oracle_.resolve_label(name);   // the header exists: no throw
            if (back.kind != LabelKind::HEADER || back.column != key.first
                    || back.seq_id != key.second) {
                const std::string resolves_to = back.kind == LabelKind::COLUMN
                    ? "is also the name of an annotation column"
                    : "also occurs in another annotation column ("
                          + oracle_.column_name(back.column) + ")";
                throw SeedDerivationError(
                        "The labels derived from the seed are not resubmittable: the sequence "
                        "header '" + name + "' " + resolves_to + ", which an explicit label list "
                        "would resolve it to. Set seed_label_kind to \"column\" or name the "
                        "labels explicitly.");
            }
        }
    }
    // The cap bounds the traversal state, not the discovery: everything above it was
    // found anyway, so it is counted and digested rather than silently forgotten.
    // Under `exhaustive` it is not cut at all: the preset promises that no walk is
    // dropped, and every walk of a dropped carrier would be. Refuse instead, naming
    // the two levers the caller has.
    if (live.size() > strategy_.max_seed_labels && strategy_.exhaustive) {
        throw SeedDerivationError(
                std::to_string(live.size()) + " labels carry the seed and max_seed_labels is "
                + std::to_string(strategy_.max_seed_labels) + ": under `exhaustive` the "
                  "derived set is not truncated (every walk of a dropped carrier would be "
                  "missing from the trie); raise labels.max_seed_labels or name the labels");
    }
    if (live.size() > strategy_.max_seed_labels) {
        uint64_t digest = kFnvOffsetBasis;
        for (size_t i = strategy_.max_seed_labels; i < live.size(); ++i) {
            digest = fnv1a64(name_of(live[i]) + "\n", digest);
        }
        result_.labels_dropped = live.size() - strategy_.max_seed_labels;
        result_.labels_dropped_digest = hex64(digest);
        live.resize(strategy_.max_seed_labels);
    }
    refs->reserve(live.size());
    for (const Key &key : live) {
        LabelRef ref;
        ref.kind = kind;
        ref.column = key.first;
        ref.seq_id = kind == LabelKind::HEADER ? key.second : 0;
        ref.name = name_of(key);
        refs->push_back(std::move(ref));
    }
    hits->clear();
    if (!trace_) {
        // Every derived label supports every k-mer by construction and there are no
        // coordinates to carry over, so the per-label validation needs no hit matrix at
        // all: materialising one would be |live| x |keys| Hit objects (32 B each) saying
        // nothing but that. The caller is told so instead.
        return true;
    }
    // the hits of this pass in the shape the per-label validation expects: only the
    // coordinates (which `support: trace` still has to find consecutive) are carried over
    hits->assign(keys.size(), LabelQuery::NodeHits());
    for (size_t i = 0; i < keys.size(); ++i) {
        LabelQuery::NodeHits &node = (*hits)[i];
        node.reserve(live.size());
        size_t j = 0;
        for (LabelId l = 0; l < live.size(); ++l) {
            while (j < coords_at[i].size() && coords_at[i][j].first < live[l]) {
                ++j;
            }
            if (j < coords_at[i].size() && coords_at[i][j].first == live[l])
                node.push_back(LabelQuery::Hit{ l, std::move(coords_at[i][j].second) });
        }
    }
    return false;
}

void Walker::init_edge_coding() {
    code_.fill(255);
    size_t n = 0;
    for (char c : graph_.alphabet()) {
        if (c != kSentinel)
            code_[static_cast<unsigned char>(c)] = n++;
    }
    bits_ = 1;
    while ((size_t(1) << bits_) < n) {
        ++bits_;
    }
    packable_ = (k_ + 1) * bits_ <= 64;
}

void Walker::init_arm(ArmState &arm) {
    arm.result.arm = arm.arm;
    arm.result.requested = strategy_.direction == Strategy::BOTH
        || (strategy_.direction == Strategy::LEFT) == (arm.arm == Arm::LEFT);
    if (!arm.result.requested)
        return;

    Item root;
    root.node = arm.arm == Arm::RIGHT ? nodes_.back() : nodes_.front();
    root.kmer = arm.arm == Arm::RIGHT ? seed_upper_.substr(seed_upper_.size() - k_)
                                      : seed_upper_.substr(0, k_);
    for (LabelId l = 0; l < result_.num_seed_labels; ++l) {
        Entry e;
        e.label = l;
        e.pred = l;
        if (trace_)
            e.coords = boundary_coords_[static_cast<size_t>(arm.arm)][l];
        root.state.push_back(e);
    }
    if (annotate_) {
        // the boundary k-mer's own labels: the root's entry node, from which
        // continuous presence (label_summary.direct_bp) is measured
        auto nl = recorder_->fetch({ oracle_.key_of(root.node, root.kmer) });
        root.present = bounded(arm, nl[0].labels, nl[0].total);
        root.present_total = nl[0].total;
    }
    root.segment = new_segment(arm, {}, 0, labels_at(root),
                               annotate_ ? root.present_total : root.state.size());
    for (Entry &e : root.state) {
        e.run = new_run(arm, e.label, 0, false, 0, 0, UINT32_MAX);
    }
    arm.first_arrival.emplace(root.node, std::make_pair(root.segment, uint64_t(0)));
    arm.frontier.push_back(std::move(root));
}


/******************************* graph access *******************************/

std::vector<Succ> Walker::enumerate(ArmState &arm, node_index node, const std::string &kmer) {
    std::vector<Succ> out;
    auto cb = [&](node_index n, char c) {
        if (c != kSentinel && n != npos)
            out.push_back({ n, c });
    };
    try {
        if (arm.arm == Arm::RIGHT) {
            if (canonical_) {
                canonical_->call_outgoing_kmers(node, kmer, cb);
            } else {
                graph_.call_outgoing_kmers(node, cb);
            }
        } else {
            if (canonical_) {
                canonical_->call_incoming_kmers(node, kmer, cb);
            } else if (nfc_) {
                nfc_->call_incoming_kmers(node, cb);
            } else {
                graph_.call_incoming_kmers(node, cb);
            }
        }
    } catch (const std::runtime_error &e) {
        // CanonicalDBG throws an empty runtime_error on an inconsistent primary graph
        throw std::runtime_error("Inconsistent primary graph at node " + std::to_string(node)
                                 + " (" + kmer + ") on arm " + to_string(arm.arm)
                                 + (e.what()[0] ? std::string(": ") + e.what() : std::string()));
    }
    std::sort(out.begin(), out.end(), [](const Succ &a, const Succ &b) { return a.ch < b.ch; });
    return out;
}

// Successors of |node| from the lookahead cache or the graph; |*cached_key| tells
// whether |*key| holds the annotation key of the single successor. Enumerations and
// keys count when the walk consumes them, not when the lookahead computed them.
std::vector<Succ> Walker::successors(ArmState &arm, node_index node, const std::string &kmer,
                                     bool *cached_key, node_index *key) {
    *cached_key = false;
    arm.result.successor_enumerations++;
    auto it = arm.lookahead.find(node);
    if (it != arm.lookahead.end()) {
        Lookahead la = std::move(it.value());
        arm.lookahead.erase(it);
        if (la.succs.size() == 1) {
            *cached_key = true;
            *key = la.succ_key;
            oracle_.counters().keys_mapped++;
        }
        return std::move(la.succs);
    }
    return enumerate(arm, node, kmer);
}

node_index Walker::key_of_succ(Arm arm, const std::string &kmer, const Succ &s) {
    // the spelling is only needed by native canonical graphs
    if (regime_ != Regime::CANONICAL)
        return oracle_.key_of(s.node, std::string_view());
    succ_kmer(arm, kmer, s.ch, &rc_scratch_);
    return oracle_.key_of(s.node, rc_scratch_);
}

// the k-mer of the successor reached by |c| (|out| must not alias |kmer|)
void Walker::succ_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const {
    out->clear();
    if (arm == Arm::RIGHT) {
        out->append(kmer, 1, k_ - 1);
        out->push_back(c);
    } else {
        out->push_back(c);
        out->append(kmer, 0, k_ - 1);
    }
}

// the (k+1)-mer of a step in natural orientation
void Walker::step_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const {
    out->clear();
    if (arm == Arm::RIGHT) {
        out->append(kmer);
        out->push_back(c);
    } else {
        out->push_back(c);
        out->append(kmer);
    }
}

bool Walker::is_hairpin(const std::string &kmer_u, const std::string &kmer_v,
                        const std::string &step) {
    std::string &rc = rc_scratch_;
    rc = step;
    ::reverse_complement(rc);
    if (rc == step)
        return true;
    if (k_ % 2 == 0) {
        rc = kmer_v;
        ::reverse_complement(rc);
        if (rc == kmer_v)
            return true;
        rc = kmer_u;
        ::reverse_complement(rc);
        if (rc == kmer_u)
            return true;
    }
    return false;
}

std::pair<EdgeKey, bool> Walker::edge_key(const std::string &step) {
    std::string_view s = step;
    bool flip = false;
    if (regime_ != Regime::BASIC) {
        rc_scratch_ = step;
        ::reverse_complement(rc_scratch_);
        if (rc_scratch_ < step) {
            s = rc_scratch_;
            flip = true;
        }
    }
    EdgeKey key { 0, 0 };
    if (packable_) {
        for (char c : s) {
            key.lo = (key.lo << bits_) | code_[static_cast<unsigned char>(c)];
        }
    } else {
        key.lo = fnv1a64(s);
        key.hi = fnv1a64(s, 0x84222325cbf29ce4ULL);
    }
    return { key, flip };
}


/******************************** bookkeeping *******************************/

GrowthBin& Walker::bin(ArmState &arm, uint64_t bp) {
    uint64_t width = std::max<uint64_t>(1, strategy_.profile_bin_bp);
    size_t index = bp / width;
    auto &growth = arm.result.growth;
    while (growth.size() <= index) {
        GrowthBin b;
        b.from_bp = growth.size() * width;
        growth.push_back(b);
    }
    return growth[index];
}

size_t Walker::new_segment(ArmState &arm, std::vector<size_t> parents, uint64_t from_bp,
                           std::vector<LabelId> labels_start, size_t labels_start_total) {
    Segment seg;
    seg.id = arm.result.segments.size();
    seg.parents = std::move(parents);
    seg.from_bp = from_bp;
    seg.labels_start = std::move(labels_start);
    seg.labels_start_total = labels_start_total;
    arm.result.segments.push_back(std::move(seg));
    arm.walk_seq.emplace_back();
    arm.leaves.emplace_back();
    return arm.result.segments.size() - 1;
}

uint32_t Walker::new_run(ArmState &arm, LabelId label, uint64_t from_bp, bool by_switch,
                         LabelId from, double cost, uint32_t prev) {
    LabelRun run;
    run.label = label;
    run.from_bp = from_bp;
    run.to_bp = from_bp;
    run.entered_by_switch = by_switch;
    run.from_label = from;
    run.switch_cost = cost;
    run.prev_run = prev;
    arm.result.runs.push_back(run);
    return arm.result.runs.size() - 1;
}

// Mark every ancestor of |seg| (once per segment: the parents of a segment never
// change after its creation, so the marks stay valid until another segment is marked)
void Walker::mark_ancestors(ArmState &arm, size_t seg) {
    if (arm.marked_segment == seg)
        return;
    const auto &segments = arm.result.segments;
    if (arm.visit_mark.size() < segments.size())
        arm.visit_mark.resize(segments.size(), 0);
    ++arm.visit_epoch;
    arm.visit_stack.assign(1, seg);
    arm.ancestors.clear();
    while (!arm.visit_stack.empty()) {
        size_t s = arm.visit_stack.back();
        arm.visit_stack.pop_back();
        for (size_t p : segments[s].parents) {
            if (arm.visit_mark[p] != arm.visit_epoch) {
                arm.visit_mark[p] = arm.visit_epoch;
                arm.visit_stack.push_back(p);
                arm.ancestors.push_back(p);
            }
        }
    }
    arm.marked_segment = seg;
}

// requires mark_ancestors(arm, seg)
bool Walker::is_ancestor_or_self(const ArmState &arm, size_t anc, size_t seg) const {
    assert(arm.marked_segment == seg);
    if (anc == seg)
        return true;
    // segments are numbered in creation order
    return anc < seg && arm.visit_mark[anc] == arm.visit_epoch;
}

size_t Walker::distinct_labels(const std::vector<Item> &a, const std::vector<Item> &b,
                               size_t from_a, bool *exact, uint64_t below_bp) {
    ++label_epoch_;
    size_t n = 0;
    *exact = true;
    auto mark = [&](LabelId l) {
        // the dictionary grows during the walk in annotate mode
        if (label_stamp_.size() <= l)
            label_stamp_.resize(l + 1, 0);
        if (label_stamp_[l] != label_epoch_) {
            label_stamp_[l] = label_epoch_;
            ++n;
        }
    };
    auto count = [&](const Item &item) {
        if (item.ext_bp >= below_bp)
            return;
        if (annotate_) {
            // a cut list hides labels: the count over it is a lower bound
            if (item.present_total > item.present.size())
                *exact = false;
            for (LabelId l : item.present) {
                mark(l);
            }
            return;
        }
        for (const Entry &e : item.state) {
            mark(e.label);
        }
    };
    for (size_t i = from_a; i < a.size(); ++i) {
        count(a[i]);
    }
    for (const Item &item : b) {
        count(item);
    }
    return n;
}

void Walker::record_live(ArmState &arm, uint64_t bp, const std::vector<Item> &items) {
    size_t pairs = 0;
    for (const Item &item : items) {
        pairs += annotate_ ? item.present.size() : item.state.size();
    }
    bool exact = true;
    size_t labels = distinct_labels(items, {}, 0, &exact);
    GrowthBin &b = bin(arm, bp);
    b.max_live_paths = std::max(b.max_live_paths, items.size());
    b.max_live_labels = std::max(b.max_live_labels, labels);
    b.max_live_pairs = std::max(b.max_live_pairs, pairs);
    b.live_labels_exact &= exact;
}


/*********************************** ends ***********************************/

void Walker::end_run(ArmState &arm, const Entry &e, uint64_t at, EndReason reason) {
    LabelRun &run = arm.result.runs[e.run];
    run.to_bp = at;
    run.ended = true;
    run.end_reason = reason;
    bin(arm, at).label_ends[static_cast<size_t>(reason)]++;
}

void Walker::end_label(ArmState &arm, const Item &item, const Entry &e, EndReason reason,
                       size_t structural, const char *text, double needed) {
    end_run(arm, e, item.ext_bp, reason);
    Event ev;
    ev.at_bp = item.ext_bp;
    ev.type = EventType::LABEL_END;
    ev.label = e.label;
    ev.reason = reason;
    ev.structural_successors = structural;
    ev.needed_budget = needed;
    ev.text = text;
    arm.result.segments[item.segment].events.push_back(std::move(ev));
    if (reason == EndReason::LOSS_BUDGET)
        arm.result.needed_budgets.push_back(needed);
}

void Walker::finish_path(ArmState &arm, Item &item, std::optional<EndReason> path_reason) {
    Segment &seg = arm.result.segments[item.segment];
    seg.labels_end = labels_at(item);
    LeafInfo &leaf = arm.leaves[item.segment];
    leaf.is_leaf = true;
    leaf.path_reason = path_reason;
    for (const Entry &e : item.state) {
        const LabelRun &run = arm.result.runs[e.run];
        assert(run.ended);
        leaf.end_reasons[static_cast<size_t>(run.end_reason)]++;
        leaf.end_labels.push_back({ e.label, e.loss, e.branches, e.run, e.route_bp });
    }
    // a continuation is for a walk that could go on: the radius or a cap, never a
    // structural end (annotate mode ends every path with a path_reason)
    if (path_reason && (is_resource_stop(*path_reason) || *path_reason == EndReason::MAX_EXTENSION))
        leaf.continuation = make_continuation(arm, item);
    arm.finished_leaves++;
}

void Walker::censor_item(ArmState &arm, Item &item, EndReason reason) {
    for (const Entry &e : item.state) {
        end_label(arm, item, e, reason, 0, "");
    }
    finish_path(arm, item, reason);
}

// The last max(k, bp since the last switch or route join) bases of seed + flank, up
// to continuation_bp, with the labels covering that whole tail (§7.1). Only the tail
// is spelled: O(continuation_bp) per leaf, not O(path length).
Continuation Walker::make_continuation(ArmState &arm, const Item &item) {
    Continuation c;
    uint64_t boundary = 0;
    for (const Entry &e : item.state) {
        const LabelRun &run = arm.result.runs[e.run];
        if (run.entered_by_switch)
            boundary = std::max(boundary, run.from_bp);
        boundary = std::max(boundary, e.route_bp);
    }
    uint64_t tail = item.ext_bp - boundary;
    uint64_t n = std::min<uint64_t>(strategy_.continuation_bp, std::max<uint64_t>(k_, tail));
    n = std::min<uint64_t>(n, seed_upper_.size() + item.ext_bp);
    const uint64_t from_walk = std::min<uint64_t>(n, item.ext_bp);

    // the last |from_walk| walked bases, collected leaf -> root (reverse walking order)
    std::string rev;
    rev.reserve(from_walk);
    for (size_t s = item.segment; rev.size() < from_walk; ) {
        const std::string &w = arm.walk_seq[s];
        for (auto it = w.rbegin(); it != w.rend() && rev.size() < from_walk; ++it) {
            rev.push_back(*it);
        }
        if (arm.result.segments[s].parents.empty())
            break;
        s = arm.result.segments[s].parents[0];
    }
    assert(rev.size() == from_walk);
    const uint64_t from_seed = n - from_walk;
    if (arm.arm == Arm::RIGHT) {
        c.sequence = seed_upper_.substr(seed_upper_.size() - from_seed);
        c.sequence.append(rev.rbegin(), rev.rend());
    } else {
        // natural orientation of the left flank is the reverse of the walking order
        c.sequence = rev + seed_upper_.substr(0, from_seed);
    }
    // labels whose current run and route cover the whole tail inside the flank
    uint64_t covered_from = n >= item.ext_bp ? 0 : item.ext_bp - n;
    c.loss_used = kInfiniteLoss;
    c.branches_used = UINT32_MAX;
    if (annotate_) {
        // the labels recorded on EVERY node of the tail (an intersection over the
        // runs covering it; a lower bound where a run's list was cut, which the arm's
        // nodes_labels_truncated reports). Such a label validates as a seed label of
        // the continuation, so the tail stays valid /traverse input.
        std::vector<LabelId> alive = item.present;
        for (size_t s = item.segment; !alive.empty(); ) {
            const Segment &seg = arm.result.segments[s];
            for (auto it = seg.label_sets.rbegin(); it != seg.label_sets.rend() && !alive.empty(); ++it) {
                if (it->to_bp <= covered_from)
                    break;
                std::vector<LabelId> still;
                std::set_intersection(alive.begin(), alive.end(), it->labels.begin(),
                                      it->labels.end(), std::back_inserter(still));
                alive.swap(still);
            }
            if (seg.from_bp <= covered_from || seg.parents.empty())
                break;
            s = seg.parents[0];
        }
        c.labels = std::move(alive);
        c.loss_used = 0;
        c.branches_used = 0;
        return c;
    }
    for (const Entry &e : item.state) {
        const LabelRun &run = arm.result.runs[e.run];
        if (run.from_bp <= covered_from && e.route_bp <= covered_from) {
            c.labels.push_back(e.label);
            c.loss_used = std::min(c.loss_used, e.loss);
            c.branches_used = std::min(c.branches_used, e.branches);
        }
    }
    if (c.labels.empty()) {
        // only reachable with n == k: the tail is the current k-mer, which every
        // live label supports
        for (const Entry &e : item.state) {
            c.labels.push_back(e.label);
            c.loss_used = std::min(c.loss_used, e.loss);
            c.branches_used = std::min(c.branches_used, e.branches);
        }
    }
    if (c.labels.empty()) {
        c.loss_used = 0;
        c.branches_used = 0;
    }
    return c;
}


/******************************* label state ********************************/

void Walker::derive(const State &sigma, const std::vector<Target> &targets,
                    const std::vector<uint8_t> &excluded, State *out,
                    bool *truncated, SourceKey *cut) {
    out->clear();
    *truncated = false;
    const double budget = strategy_.loss_budget;
    const bool loss_only = strategy_.switch_on_loss_only;

    auto in_targets = [&](LabelId l) {
        auto it = std::lower_bound(targets.begin(), targets.end(), l,
                                   [](const Target &t, LabelId x) { return t.label < x; });
        return it != targets.end() && it->label == l;
    };

    // sources eligible for a switch
    std::vector<const Entry*> &sw = sw_;
    sw.clear();
    const Entry *b1 = nullptr, *b2 = nullptr;
    if (cost_.finite()) {
        for (const Entry &e : sigma) {
            if (excluded[e.label])
                continue;
            if (loss_only && in_targets(e.label))
                continue;
            sw.push_back(&e);
        }
        if (cost_.model() == LabelChangeCost::CONSTANT) {
            for (const Entry *e : sw) {
                if (!b1 || source_key(*e) < source_key(*b1)) {
                    b2 = b1;
                    b1 = e;
                } else if (!b2 || source_key(*e) < source_key(*b2)) {
                    b2 = e;
                }
            }
        } else {
            std::sort(sw.begin(), sw.end(), [&](const Entry *a, const Entry *b) {
                return source_key(*a) < source_key(*b);
            });
            if (sw.size() > strategy_.max_switch_sources) {
                sw.resize(strategy_.max_switch_sources);
                *truncated = true;
                // no source at all: a key below every real one
                *cut = sw.empty() ? SourceKey{ -1.0, 0, 0, 0 } : source_key(*sw.back());
            }
        }
    }

    for (const Target &t : targets) {
        const Entry *stay = find_entry(sigma, t.label);
        if (stay && excluded[t.label])
            stay = nullptr;
        double stay_loss = stay ? stay->loss : kInfiniteLoss;

        const Entry *best = nullptr;
        double best_loss = kInfiniteLoss;
        if (cost_.model() == LabelChangeCost::CONSTANT) {
            const Entry *s = (b1 && b1->label != t.label) ? b1 : b2;
            if (s) {
                best = s;
                best_loss = s->loss + cost_.cost(s->label, t.label);
            }
            pair_evaluations_ += 1;
        } else if (cost_.model() == LabelChangeCost::TABLE) {
            for (const Entry *s : sw) {
                if (s->label == t.label)
                    continue;
                double v = s->loss + cost_.cost(s->label, t.label);
                if (v < best_loss
                        || (v == best_loss && best && source_key(*s) < source_key(*best))) {
                    best = s;
                    best_loss = v;
                }
            }
            pair_evaluations_ += sw.size();
        }

        Entry e;
        e.label = t.label;
        e.coords = t.coords;
        if (stay && stay_loss <= best_loss) {
            e.loss = stay_loss;
            e.branches = stay->branches;
            e.run = stay->run;
            e.route_bp = stay->route_bp;
            e.pred = t.label;
            e.switched = false;
        } else if (best && best_loss <= budget) {
            e.loss = best_loss;
            e.branches = best->branches;
            e.run = UINT32_MAX;
            e.route_bp = best->route_bp;
            e.pred = best->label;
            e.switched = true;
            e.from = best->label;
            e.switch_cost = cost_.cost(best->label, t.label);
        } else {
            continue;
        }
        out->push_back(std::move(e));
    }
}

void Walker::commit_entries(ArmState &arm, const Item &item, size_t target_segment,
                            State &state, uint64_t at, std::vector<uint32_t> *taken) {
    for (Entry &e : state) {
        if (e.switched) {
            const Entry *src = find_entry(item.state, e.from);
            assert(src);
            e.run = new_run(arm, e.label, at, true, e.from, e.switch_cost, src->run);
            Event ev;
            ev.at_bp = at;
            ev.type = EventType::SWITCH;
            ev.label = e.from;
            ev.to = e.label;
            ev.cost = e.switch_cost;
            arm.result.segments[target_segment].events.push_back(std::move(ev));
        } else if (taken) {
            if (std::find(taken->begin(), taken->end(), e.run) != taken->end()) {
                LabelRun src = arm.result.runs[e.run];
                e.run = new_run(arm, src.label, src.from_bp, src.entered_by_switch,
                                src.from_label, src.switch_cost, src.prev_run);
            } else {
                taken->push_back(e.run);
            }
        }
        e.switched = false;
        e.pred = e.label;
    }
}


/********************************** levels **********************************/

bool Walker::time_exceeded() const {
    return strategy_.time_budget_ms <= 0
        || timer_.elapsed() * 1000.0 >= strategy_.time_budget_ms;
}

void Walker::sort_items(std::vector<Item> &items) const {
    switch (strategy_.order) {
        case Strategy::BREADTH_FIRST:
            std::sort(items.begin(), items.end(),
                      [](const Item &a, const Item &b) { return a.path_id < b.path_id; });
            break;
        case Strategy::LOWEST_LOSS_FIRST:
            std::sort(items.begin(), items.end(), [](const Item &a, const Item &b) {
                return std::make_pair(min_loss(a.state), a.path_id)
                     < std::make_pair(min_loss(b.state), b.path_id);
            });
            break;
        case Strategy::MOST_SUPPORTED_FIRST: {
            // support = the labels alive on the head (constrain) or recorded at its node
            // (annotate), so that a label-free beam follows the majority continuation at a
            // fork instead of whichever head was created first (spec §6.11). The TRUE
            // count, not the recorded list: that list is cut at max_labels_per_node, and
            // ranking by it would tie a 65-label branch with a 10,000-label one
            auto support = [this](const Item &item) {
                return annotate_ ? item.present_total : item.state.size();
            };
            std::sort(items.begin(), items.end(), [&](const Item &a, const Item &b) {
                return std::make_tuple(support(b), min_loss(a.state), a.path_id)
                     < std::make_tuple(support(a), min_loss(b.state), b.path_id);
            });
            break;
        }
    }
}

// Structural lookahead along unbranched runs (§8.3): the successors of up to
// |batch_kmers| nodes go to the lookahead cache, their annotation keys are mapped
// in one sequential call per chunk and their rows warm the label cache. Nothing
// here changes which steps are taken; the counters are charged on consumption.
void Walker::prefetch(ArmState &arm, const std::vector<Item> &items,
                      const std::vector<std::vector<Succ>> &succs) {
    if (!strategy_.batch_kmers)
        return;
    std::vector<node_index> warm_keys;
    std::vector<node_index> chain_curs;      // nodes whose successors were enumerated
    std::vector<Lookahead> chain_entries;    // their lookahead entries
    std::vector<node_index> chain_nodes;     // single successors, in walking order
    std::string window;                      // spelling of |chain_nodes|, natural orientation
    std::string kmer, next_kmer;
    for (size_t i = 0; i < items.size(); ++i) {
        if (succs[i].size() != 1)
            continue;
        node_index cur = succs[i][0].node;
        if (arm.lookahead.count(cur))
            continue;
        succ_kmer(arm.arm, items[i].kmer, succs[i][0].ch, &kmer);
        chain_curs.clear();
        chain_entries.clear();
        chain_nodes.clear();
        window.clear();
        for (size_t n = 0; n < strategy_.batch_kmers; ++n) {
            if (n > 0 && arm.lookahead.count(cur))
                break;
            Lookahead la;
            la.succs = enumerate(arm, cur, kmer);
            chain_curs.push_back(cur);
            const bool single = la.succs.size() == 1;
            const char c = single ? la.succs[0].ch : 0;
            const node_index nxt = single ? la.succs[0].node : npos;
            chain_entries.push_back(std::move(la));
            if (!single)
                break;
            succ_kmer(arm.arm, kmer, c, &next_kmer);
            kmer.swap(next_kmer);
            if (chain_nodes.empty()) {
                window = kmer;
            } else if (arm.arm == Arm::RIGHT) {
                window.push_back(c);
            } else {
                window.insert(window.begin(), c);
            }
            chain_nodes.push_back(nxt);
            cur = nxt;
        }
        if (!chain_nodes.empty()) {
            std::vector<node_index> keys;
            if (arm.arm == Arm::RIGHT) {
                keys = oracle_.keys_of_path(chain_nodes, window);
            } else {
                // the window is in natural orientation: the walked nodes appear reversed
                std::vector<node_index> natural(chain_nodes.rbegin(), chain_nodes.rend());
                keys = oracle_.keys_of_path(natural, window);
                std::reverse(keys.begin(), keys.end());
            }
            // charged when the walk consumes them (successors())
            oracle_.counters().keys_mapped -= keys.size();
            for (size_t j = 0; j < keys.size(); ++j) {
                chain_entries[j].succ_key = keys[j];
                warm_keys.push_back(keys[j]);
            }
        }
        for (size_t j = 0; j < chain_curs.size(); ++j) {
            arm.lookahead.emplace(chain_curs[j], std::move(chain_entries[j]));
        }
    }
    if (!warm_keys.empty()) {
        if (annotate_) {
            recorder_->warm(warm_keys);
        } else {
            query_->warm(warm_keys);
        }
    }
    if (arm.lookahead.size() > kMaxLookahead) {
        for (const auto &kv : arm.lookahead) {
            if (kv.second.succs.size() == 1)
                oracle_.counters().keys_mapped++;
        }
        arm.lookahead.clear();
    }
}

void Walker::arrival(ArmState &arm, Item &item) {
    auto it = arm.first_arrival.find(item.node);
    if (it == arm.first_arrival.end()) {
        arm.first_arrival.emplace(item.node, std::make_pair(item.segment, item.ext_bp));
        item.revisiting = false;
        return;
    }
    if (it->second.second != item.ext_bp) {
        // reported once per stretch of revisited nodes
        if (item.revisiting)
            return;
        item.revisiting = true;
        Event ev;
        ev.at_bp = item.ext_bp;
        ev.type = EventType::REVISIT;
        ev.labels = { static_cast<LabelId>(it->second.first) };
        ev.length_bp = item.ext_bp > it->second.second ? item.ext_bp - it->second.second
                                                        : it->second.second - item.ext_bp;
        arm.result.segments[item.segment].events.push_back(std::move(ev));
    }
    // same distance: resolved in merge_level
}

void Walker::merge_level(ArmState &arm, uint64_t depth) {
    if (arm.next.size() < 2)
        return;
    // group by node, keeping the order of first occurrence
    tsl::hopscotch_map<node_index, size_t> group_of;
    std::vector<std::vector<size_t>> groups;
    for (size_t i = 0; i < arm.next.size(); ++i) {
        auto it = group_of.find(arm.next[i].node);
        if (it == group_of.end()) {
            group_of.emplace(arm.next[i].node, groups.size());
            groups.push_back({ i });
        } else {
            groups[it->second].push_back(i);
        }
    }
    std::vector<Item> merged;
    merged.reserve(arm.next.size());
    for (const auto &g : groups) {
        if (g.size() == 1) {
            merged.push_back(std::move(arm.next[g[0]]));
            continue;
        }
        if (!strategy_.merge_reconverge) {
            for (size_t j = 0; j < g.size(); ++j) {
                Item &item = arm.next[g[j]];
                // reported once per stretch of shared nodes, like arrival()
                if (j > 0 && !item.revisiting) {
                    item.revisiting = true;
                    Event ev;
                    ev.at_bp = depth;
                    ev.type = EventType::REVISIT;
                    ev.labels = { static_cast<LabelId>(arm.next[g[0]].segment) };
                    ev.text = "same_distance";
                    arm.result.segments[item.segment].events.push_back(std::move(ev));
                }
                merged.push_back(std::move(item));
            }
            continue;
        }
        // union of states keeping the min (loss, branches) entry per label; an entry
        // taken from a later parent is routed from here on (route_bp); under trace
        // support the live coordinates of both parents are united
        Item primary = std::move(arm.next[g[0]]);
        State st = primary.state;
        std::vector<size_t> parents { primary.segment };
        size_t path_id = primary.path_id;
        uint32_t splits = primary.splits;
        for (size_t j = 1; j < g.size(); ++j) {
            Item &other = arm.next[g[j]];
            parents.push_back(other.segment);
            path_id = std::min(path_id, other.path_id);
            splits = std::max(splits, other.splits);
            // A lineage routed in from a later parent stops being evidence for the bases
            // the surviving path spells, so the route is recorded on the run too (the
            // earliest such depth wins, since runs are shared across sibling paths).
            auto stamp_route = [&](uint32_t run_id, uint64_t at) {
                if (run_id == UINT32_MAX)
                    return;
                uint64_t &stamped = arm.result.runs[run_id].route_bp;
                stamped = stamped ? std::min(stamped, at) : at;
            };
            // both states are sorted by label: one pass over the two, into a fresh
            // vector (inserting each incoming entry in place would cost their product)
            State out;
            out.reserve(st.size() + other.state.size());
            auto a = st.begin();
            for (Entry &e : other.state) {
                while (a != st.end() && a->label < e.label) {
                    out.push_back(std::move(*a++));
                }
                if (a == st.end() || a->label != e.label) {
                    Entry copy = e;
                    copy.route_bp = depth;
                    stamp_route(copy.run, depth);
                    out.push_back(std::move(copy));
                    continue;
                }
                Entry &cur = *a++;
                bool better = std::make_pair(e.loss, e.branches)
                            < std::make_pair(cur.loss, cur.branches);
                Entry &closed = better ? cur : e;
                // closed by the merge: not a label end, the lineage continues
                LabelRun &run = arm.result.runs[closed.run];
                run.to_bp = depth;
                run.ended = false;
                if (trace_)
                    unite_coords(&e.coords, cur.coords);
                if (better) {
                    Entry kept = e;
                    kept.route_bp = depth;
                    stamp_route(kept.run, depth);
                    out.push_back(std::move(kept));
                } else {
                    if (trace_)
                        cur.coords = e.coords;
                    out.push_back(std::move(cur));
                }
            }
            while (a != st.end()) {
                out.push_back(std::move(*a++));
            }
            st.swap(out);
        }
        // in annotate mode the join is purely structural (same node, same depth) and
        // the merged segment's entry labels are the node's own
        size_t m = new_segment(arm, parents, depth,
                               annotate_ ? primary.present : labels_of(st),
                               annotate_ ? primary.present_total : st.size());
        Segment &mseg = arm.result.segments[m];
        mseg.labels_via_parent.assign(parents.size(), {});
        for (size_t j = 0; j < g.size(); ++j) {
            const Item &item = j == 0 ? primary : arm.next[g[j]];
            // runs are per path, so the kept entry's run identifies its parent
            for (const Entry &e : item.state) {
                const Entry *kept = find_entry(st, e.label);
                if (kept && kept->run == e.run)
                    mseg.labels_via_parent[j].push_back(e.label);
            }
            Segment &seg = arm.result.segments[item.segment];
            seg.children.push_back(m);
            seg.labels_end = labels_at(item);
        }
        Event ev;
        ev.at_bp = depth;
        ev.type = EventType::RECONVERGE;
        for (size_t p : parents) {
            ev.labels.push_back(static_cast<LabelId>(p));
        }
        mseg.events.push_back(std::move(ev));
        bin(arm, depth).reconvergences++;
        arm.first_arrival[primary.node] = { m, depth };

        primary.segment = m;
        primary.state = std::move(st);
        primary.path_id = path_id;
        primary.splits = splits;
        merged.push_back(std::move(primary));
    }
    arm.next.swap(merged);
}

void Walker::beam(ArmState &arm, uint64_t depth) {
    if (strategy_.on_overflow != Strategy::BEAM || arm.next.size() <= strategy_.max_live_paths)
        return;
    // heads at the radius are complete and end as such at the next level: pruning them
    // would report a beam cut on a level that is in fact entirely present
    if (depth >= strategy_.max_extension_bp)
        return;
    std::vector<size_t> order(arm.next.size());
    for (size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    // most supported first, whatever Strategy::order (which only sets the expansion
    // order within a level): the labels alive on the head (constrain) or recorded at
    // its node (annotate, the true count — the recorded list is cut at
    // max_labels_per_node and must not decide the beam) — a label-free beam then
    // follows the majority continuation at a fork rather than whichever head was
    // created first (spec §6.11)
    auto support = [this](const Item &item) {
        return annotate_ ? item.present_total : item.state.size();
    };
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const Item &x = arm.next[a], &y = arm.next[b];
        return std::make_tuple(support(y), min_loss(x.state), x.path_id)
             < std::make_tuple(support(x), min_loss(y.state), y.path_id);
    });
    std::vector<Item> kept;
    bool exact = true;
    const size_t live_labels = distinct_labels(arm.next, {}, 0, &exact);
    for (size_t i = 0; i < order.size(); ++i) {
        Item &item = arm.next[order[i]];
        if (i < strategy_.max_live_paths) {
            kept.push_back(std::move(item));
        } else {
            if (!arm.result.cap_trigger) {
                arm.result.cap_trigger = CapTrigger{ EndReason::BEAM_PRUNED, depth, item.segment,
                                                     arm.next.size(), live_labels, exact };
            }
            // the pruned heads are walks of length |depth| that will not be expanded:
            // walks of that length are all present, longer ones are not
            arm.boundary = std::min(arm.boundary, item.ext_bp);
            censor_item(arm, item, EndReason::BEAM_PRUNED);
        }
    }
    if (arm.result.status == ArmResult::COMPLETE)
        arm.result.status = ArmResult::PRUNED;
    std::sort(kept.begin(), kept.end(),
              [](const Item &a, const Item &b) { return a.path_id < b.path_id; });
    arm.next.swap(kept);
}

// End every live path of |arm| with |reason|: the items of the current level from
// |from| on (if given) and the items already produced for the next level.
void Walker::stop_arm(ArmState &arm, EndReason reason, std::vector<Item> *items, size_t from) {
    const size_t heads = arm.next.size() + (items ? items->size() - from : 0);
    if (!heads)
        return;
    // A head that has already reached the radius is complete, not out of budget (the
    // same rule stop_frontier applies to the frontier): it ends as max_extension_bp
    // below and is no remaining head, so neither frontier_remaining nor the cap
    // trigger counts it or its labels.
    const uint64_t radius = strategy_.max_extension_bp;
    bool exact = true;
    size_t labels = items ? distinct_labels(*items, arm.next, from, &exact, radius)
                          : distinct_labels(arm.next, {}, 0, &exact, radius);
    size_t live = 0;
    for (const Item &item : arm.next) {
        live += item.ext_bp < radius;
    }
    if (items) {
        for (size_t j = from; j < items->size(); ++j) {
            live += (*items)[j].ext_bp < radius;
        }
    }
    const Item &first = items && from < items->size() ? (*items)[from] : arm.next.front();
    if (!arm.result.cap_trigger) {
        arm.result.cap_trigger = CapTrigger{ reason, first.ext_bp, first.segment,
                                             live, labels, exact };
    }
    // Levels are synchronous, so every head of |items| before |from| was expanded and
    // none after it: walks of length first.ext_bp are all present, longer ones are
    // not. (When only |arm.next| is left, first.ext_bp is the completed level's
    // successor depth, which is then complete too.)
    arm.boundary = std::min(arm.boundary, first.ext_bp);
    arm.result.frontier_live_paths = live;
    arm.result.frontier_live_labels = labels;
    arm.result.frontier_live_labels_exact = exact;
    arm.result.status = ArmResult::TRUNCATED;
    auto censor = [&](Item &item) {
        censor_item(arm, item, item.ext_bp >= radius ? EndReason::MAX_EXTENSION : reason);
    };
    if (items) {
        for (size_t j = from; j < items->size(); ++j) {
            censor((*items)[j]);
        }
    }
    for (Item &item : arm.next) {
        censor(item);
    }
    arm.next.clear();
    arm.frontier.clear();
    arm.stopped = true;
}

void Walker::stop_frontier(ArmState &arm, EndReason reason) {
    if (arm.frontier.empty())
        return;
    // heads that already reached the radius are complete, not out of budget: ending
    // them as such keeps "status complete <=> complete_to_bp == max_extension_bp"
    std::vector<Item> rest;
    for (Item &item : arm.frontier) {
        if (item.ext_bp >= strategy_.max_extension_bp) {
            censor_item(arm, item, EndReason::MAX_EXTENSION);
        } else {
            rest.push_back(std::move(item));
        }
    }
    arm.frontier.swap(rest);
    if (!arm.frontier.empty())
        stop_arm(arm, reason, &arm.frontier, 0);
}

void Walker::trip(ArmState &arm, EndReason reason, std::vector<Item> &items, size_t index) {
    stop_arm(arm, reason, &items, index);
    if (reason == EndReason::MAX_STEPS) {
        seed_stopped_ = true;
        for (ArmState &other : arms_) {
            if (&other != &arm)
                stop_frontier(other, reason);
        }
    }
}

void Walker::run_level(ArmState &arm, uint64_t depth) {
    std::vector<Item> items;
    items.swap(arm.frontier);
    sort_items(items);
    record_live(arm, depth, items);

    // structural successors of every item, then one batched label fetch
    std::vector<std::vector<Succ>> succs(items.size());
    std::vector<node_index> keys;
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].ext_bp >= strategy_.max_extension_bp)
            continue;
        bool cached = false;
        node_index key = npos;
        succs[i] = successors(arm, items[i].node, items[i].kmer, &cached, &key);
        if (cached && succs[i].size() == 1) {
            keys.push_back(key);
            continue;
        }
        for (const Succ &s : succs[i]) {
            keys.push_back(key_of_succ(arm.arm, items[i].kmer, s));
        }
    }
    std::vector<LabelQuery::NodeHits> hits;
    std::vector<LabelRecorder::NodeLabels> present;
    if (annotate_) {
        present = recorder_->fetch(keys);
    } else {
        hits = query_->fetch(keys);
    }
    prefetch(arm, items, succs);

    size_t offset = 0;
    for (size_t i = 0; i < items.size(); ++i) {
        Item &item = items[i];
        if (item.ext_bp >= strategy_.max_extension_bp) {
            censor_item(arm, item, EndReason::MAX_EXTENSION);
            continue;
        }
        const size_t remaining = items.size() - i - 1;
        auto cap = annotate_
            ? process_item_annotate(arm, item, succs[i], present.data() + offset, remaining)
            : process_item(arm, item, succs[i], hits.data() + offset, remaining);
        offset += succs[i].size();
        if (cap) {
            trip(arm, *cap, items, i);
            return;
        }
    }
    merge_level(arm, depth + 1);
    beam(arm, depth + 1);
    arm.frontier.swap(arm.next);
}

std::optional<EndReason> Walker::process_item(ArmState &arm, Item &item,
                                              const std::vector<Succ> &succs,
                                              const LabelQuery::NodeHits *hits,
                                              size_t remaining_in_level) {
    const uint64_t at = item.ext_bp;
    const Arm side = arm.arm;
    Scratch &sc = scratch_;
    ScratchGuard guard { sc };      // clears the touched labels on every exit

    if (succs.empty()) {
        for (const Entry &e : item.state) {
            end_label(arm, item, e, EndReason::DEAD_END, 0, "");
        }
        finish_path(arm, item, std::nullopt);
        return std::nullopt;
    }

    for (const Entry &e : item.state) {
        sc.touch(e.label);
    }
    cands_.resize(succs.size());
    for (Cand &c : cands_) {
        c.reset();
    }

    for (size_t i = 0; i < succs.size(); ++i) {
        Cand &c = cands_[i];
        c.succ = &succs[i];
        succ_kmer(side, item.kmer, succs[i].ch, &c.kmer);
        for (const auto &h : hits[i]) {
            sc.touch(h.label);
            Target t { h.label, {} };
            if (trace_) {
                const Entry *e = find_entry(item.state, h.label);
                if (e) {
                    for (Coord x : h.coords) {
                        bool ok = side == Arm::RIGHT
                            ? (x > 0 && std::binary_search(e->coords.begin(), e->coords.end(), x - 1))
                            : std::binary_search(e->coords.begin(), e->coords.end(), x + 1);
                        if (ok)
                            t.coords.push_back(x);
                    }
                    if (t.coords.empty()) {
                        sc.trace_broken[h.label] = 1;
                        continue;
                    }
                } else {
                    t.coords.assign(h.coords.begin(), h.coords.end());
                }
            }
            c.targets.push_back(std::move(t));
        }
        // LabelQuery sorts hits by label id; the recurrence relies on it
        if (!std::is_sorted(c.targets.begin(), c.targets.end(),
                            [](const Target &a, const Target &b) { return a.label < b.label; })) {
            std::sort(c.targets.begin(), c.targets.end(),
                      [](const Target &a, const Target &b) { return a.label < b.label; });
        }
        check_structure(arm, item, c);
    }

    // ---- label recurrence and per-lineage branching (fixpoint over excluded sources)
    for (Cand &c : cands_) {
        if (!c.admissible())
            derive(item.state, c.targets, sc.excluded, &c.state, &c.truncated, &c.cut);
    }
    std::vector<LabelId> ambiguous_over, ambiguous_taken;
    // The explicit per-successor refusals of this step (BranchEvent::refused), one per
    // (successor, cause) with the sources whose lineage it cut; the labels are made
    // distinct and ascending when the event is emitted. Empty on an ordinary step, so
    // nothing is allocated there.
    std::vector<BranchEvent::Refusal> refused;
    // the labels of refusal (ch, cause), created when first needed; a step has a few
    // of them at most (successors x causes). The reference is used before the next call.
    auto refusals = [&](char ch, const char *cause) -> std::vector<LabelId>& {
        for (BranchEvent::Refusal &r : refused) {
            if (r.ch == ch && std::string_view(r.cause) == cause)
                return r.labels;
        }
        refused.push_back({ ch, cause, {} });
        return refused.back().labels;
    };
    // successors with an entry of a source excluded in this round: (the smallest such
    // source, successor index)
    std::vector<std::pair<LabelId, size_t>> excluded_on;
    // Every changing round excludes at least one source for good, so the fixpoint
    // settles within |σ| + 1 rounds (the last one changes nothing). The bound is
    // explicit so that nothing here can spin; it is never what ends the loop.
    const size_t sigma_size = item.state.size();
    size_t rounds = 0;
    bool changed = true;
    while (changed && rounds <= sigma_size) {
        ++rounds;
        for (LabelId l : sc.touched) {
            sc.cont_count[l] = 0;
        }
        for (Cand &c : cands_) {
            if (!c.admissible())
                continue;
            derive(item.state, c.targets, sc.excluded, &c.state, &c.truncated, &c.cut);
            if (rounds == 1)
                c.initial_labels = c.state.size();
            if (c.hairpin)
                continue;   // a followed hairpin never counts toward ambiguity (§6.5)
            for (LabelId l : sc.touched) {
                sc.seen[l] = 0;
            }
            for (const Entry &e : c.state) {
                if (!sc.seen[e.pred]) {
                    sc.seen[e.pred] = 1;
                    sc.cont_count[e.pred]++;
                }
            }
        }
        changed = false;
        const size_t excluded_from = ambiguous_over.size();
        for (const Entry &src : item.state) {
            if (sc.excluded[src.label] || sc.cont_count[src.label] < 2)
                continue;
            if (src.branches + 1 > strategy_.max_label_branches) {
                sc.excluded[src.label] = 1;
                sc.marked[src.label] = 1;
                ambiguous_over.push_back(src.label);
                changed = true;
            }
        }
        if (!changed)
            break;
        // The next derivation removes the excluded sources' entries from every
        // successor: record which successors refused them while those entries are
        // still visible. Each successor's state is scanned once per round for all the
        // sources this round excluded — a scan per excluded source was Θ(|σ|²) at a
        // node where every source is ambiguous (review round 4, finding 2).
        excluded_on.clear();
        for (size_t i = 0; i < cands_.size(); ++i) {
            Cand &c = cands_[i];
            c.excluded_preds.clear();
            if (!c.admissible())
                continue;
            arm.result.refusal_scans += c.state.size();
            for (const Entry &e : c.state) {
                if (sc.marked[e.pred])
                    c.excluded_preds.push_back(e.pred);
            }
            if (!c.excluded_preds.empty()) {
                excluded_on.emplace_back(*std::min_element(c.excluded_preds.begin(),
                                                           c.excluded_preds.end()), i);
            }
        }
        // a new (char, "branch") group is created where the per-source scan created it
        // (by the smallest excluded source on the successor, then successor order), so
        // that the refusals keep their order
        std::sort(excluded_on.begin(), excluded_on.end());
        for (const auto &[first, i] : excluded_on) {
            const Cand &c = cands_[i];
            std::vector<LabelId> &labels = refusals(c.succ->ch, "branch");
            labels.insert(labels.end(), c.excluded_preds.begin(), c.excluded_preds.end());
        }
        for (size_t j = excluded_from; j < ambiguous_over.size(); ++j) {
            sc.marked[ambiguous_over[j]] = 0;
        }
    }
    assert(!changed);
    // the rounds after the first are re-minimisations
    arm.result.reminimisation_rounds += rounds - 1;
    arm.result.max_reminimisation_rounds = std::max(arm.result.max_reminimisation_rounds,
                                                    rounds - 1);
    for (const Entry &src : item.state) {
        if (sc.excluded[src.label] || sc.cont_count[src.label] < 2)
            continue;
        ambiguous_taken.push_back(src.label);
        for (Cand &c : cands_) {
            if (!c.admissible())
                continue;
            for (Entry &e : c.state) {
                if (e.pred == src.label)
                    e.branches = src.branches + 1;
            }
        }
    }

    // ---- quorum, min_live_labels
    std::vector<Cand*> followed;
    for (Cand &c : cands_) {
        if (!c.admissible() || c.state.empty())
            continue;
        size_t n = c.state.size();
        if (n < strategy_.min_successor_labels
                || static_cast<double>(n) < strategy_.min_successor_fraction
                                                * static_cast<double>(sigma_size)) {
            c.quorum_fail = true;
            c.quorum_text = "minority";
        } else if (n < strategy_.min_live_labels) {
            c.quorum_fail = true;
            c.quorum_text = "below_min_labels";
        } else {
            c.followed = true;
            followed.push_back(&c);
        }
    }
    if (followed.size() >= 2 && item.splits + 1 > strategy_.max_splits_per_path) {
        for (Cand *c : followed) {
            c->followed = false;
            c->quorum_fail = true;
            c->quorum_text = "split_limit";
        }
        followed.clear();
    }
    const size_t nf = followed.size();

    // ---- caps, decided before anything is committed
    if (auto cap = cap_check(arm, nf, remaining_in_level))
        return cap;

    // ---- commit: events for inadmissible label-carrying successors and followed hairpins
    for (const Cand &c : cands_) {
        if (c.state.empty())
            continue;
        if (c.admissible() && !(c.hairpin && c.followed))
            continue;
        Event ev;
        ev.at_bp = at;
        ev.ch = c.succ->ch;
        ev.labels = labels_of(c.state);
        ev.labels_total = ev.labels.size();   // the live state is never cut
        if (c.skipped) {
            ev.type = EventType::HAIRPIN;
        } else if (c.blocked) {
            ev.type = EventType::BLOCKED;
            ev.reason = c.block_reason;
            bin(arm, at).blocked_repeat++;
        } else {
            ev.type = EventType::HAIRPIN;
            ev.text = "followed";
        }
        arm.result.segments[item.segment].events.push_back(std::move(ev));
    }

    // What a switch from |src| into the admissible successor |c| would cost (§6.3), for
    // the loss-budget ends and refusals (finite costs only): |present| a target within
    // the budget (kept through a cheaper predecessor, as |src| has no entry there),
    // |cut| a finite switch derive() never priced (max_switch_sources), |over| a target
    // nobody entered that |src| reaches only above the budget, the cheapest at |needed|.
    // A constant cost prices every target but |src|'s own alike, so the targets nobody
    // entered are counted, not scanned (an entry's label is always a target): a test per
    // (source, successor), where scanning the targets for every ending source made a
    // step O(|σ| · Σ|A(v)|).
    struct BudgetLook {
        bool present = false, cut = false, over = false;
        double needed = kInfiniteLoss;
    };
    auto budget_look = [&](const Entry &src, const Cand &c) {
        BudgetLook b;
        const LabelId l = src.label;
        const bool eligible = !c.truncated || !(c.cut < source_key(src));
        if (cost_.model() == LabelChangeCost::CONSTANT) {
            arm.result.refusal_scans++;
            assert(c.state.size() <= c.targets.size());
            auto it = std::lower_bound(c.targets.begin(), c.targets.end(), l,
                                       [](const Target &t, LabelId x) { return t.label < x; });
            const bool own = it != c.targets.end() && it->label == l;
            const double v = src.loss + cost_.default_cost();
            if (c.targets.size() == (own ? 1u : 0u) || v == kInfiniteLoss)
                return b;
            // the targets other than |l| nobody entered
            const size_t absent = c.targets.size() - c.state.size()
                                    - (own && !find_entry(c.state, l) ? 1 : 0);
            if (!eligible) {
                b.cut = true;
            } else if (v <= strategy_.loss_budget) {
                b.present = true;
            } else if (absent) {
                b.over = true;
                b.needed = v;
            }
            return b;
        }
        arm.result.refusal_scans += c.targets.size();
        for (const Target &t : c.targets) {
            if (t.label == l)
                continue;
            double v = src.loss + cost_.cost(l, t.label);
            if (v == kInfiniteLoss)
                continue;
            if (!eligible) {
                b.cut = true;
            } else if (v <= strategy_.loss_budget) {
                b.present = true;
            } else if (!find_entry(c.state, t.label)) {
                b.needed = std::min(b.needed, v);
                b.over = true;
            }
        }
        return b;
    };

    // ---- per-source outcomes
    for (const Cand &c : cands_) {
        // a quorum or split-limit stop refuses the successor to every source on it
        std::vector<LabelId> *quorum_refused = c.quorum_fail && !c.state.empty()
            ? &refusals(c.succ->ch, c.quorum_text) : nullptr;
        for (const Entry &e : c.state) {
            if (c.followed) {
                sc.has_cont[e.pred] = 1;
                if (!e.switched)
                    sc.stays[e.label] = 1;
                if (e.switched && find_entry(item.state, e.label))
                    sc.superseded[e.label] = 1;
            } else if (c.quorum_fail) {
                sc.in_quorum_fail[e.pred] = 1;
                sc.qtext[e.pred] = c.quorum_text;
                quorum_refused->push_back(e.pred);
            } else if (c.blocked) {
                sc.blocked[e.pred] = std::max(sc.blocked[e.pred], block_rank(c.block_reason));
            } else if (c.skipped) {
                sc.only_hairpin[e.pred] = 1;
            }
        }
    }
    std::vector<LabelId> dropped;
    for (const Entry &src : item.state) {
        LabelId l = src.label;
        if (sc.has_cont[l]) {
            if (!sc.stays[l]) {
                // the lineage continues only under other names (switch events on commit)
                end_run(arm, src, at, EndReason::LABEL_LOST);
            }
            continue;
        }
        EndReason reason;
        const char *text = "";
        double needed = 0;
        if (sc.excluded[l]) {
            reason = EndReason::BRANCH;
            dropped.push_back(l);
        } else if (sc.in_quorum_fail[l]) {
            reason = EndReason::BRANCH;
            text = sc.qtext[l];
            dropped.push_back(l);
        } else if (sc.blocked[l]) {
            reason = block_reason_of_rank(sc.blocked[l]);
        } else if (sc.only_hairpin[l]) {
            reason = EndReason::DEAD_END;
            text = "hairpin";
        } else if (sc.trace_broken[l]) {
            reason = EndReason::RECORD_END;
        } else if (sc.superseded[l]) {
            reason = EndReason::LABEL_LOST;
            text = "superseded";
        } else {
            // a switch into a target that nobody else carries, only above the budget?
            // Only the sources derive() actually considered count (max_switch_sources).
            double best_absent = kInfiniteLoss;
            bool present_target = false, cut_finite = false;
            std::string over_budget_on;     // successors refused to l by the budget alone
            if (cost_.finite()) {
                for (const Cand &c : cands_) {
                    if (!c.admissible())
                        continue;
                    const BudgetLook b = budget_look(src, c);
                    present_target |= b.present;
                    cut_finite |= b.cut;
                    if (b.over) {
                        best_absent = std::min(best_absent, b.needed);
                        over_budget_on.push_back(c.succ->ch);
                    }
                }
                sc.budget_checked[l] = 1;
            }
            if (best_absent != kInfiniteLoss) {
                reason = EndReason::LOSS_BUDGET;
                needed = best_absent;
                dropped.push_back(l);
                for (char ch : over_budget_on) {
                    refusals(ch, "loss_budget").push_back(l);
                }
            } else {
                reason = EndReason::LABEL_LOST;
                if (present_target) {
                    text = "superseded";
                } else if (cut_finite) {
                    text = "switch_sources";
                }
            }
        }
        end_label(arm, item, src, reason, succs.size(), text, needed);
    }

    // ---- loss-budget refusals of the sources whose end was not decided by the budget
    // above: a successor on which a lineage could go on only by a switch above the
    // budget is refused to it by the budget whether the lineage continues on another
    // successor or ends for another reason — only the label END depends on that
    // (review round 4, finding 3). A source excluded by the branch limit is left out:
    // the limit took it out of every successor's source set, so no switch of it was
    // priced, and its refusals are the "branch" ones.
    if (cost_.finite()) {
        for (const Cand &c : cands_) {
            // every target entered: the budget refused nothing here
            if (!c.admissible() || c.targets.size() == c.state.size())
                continue;
            for (const Entry &e : c.state) {
                sc.marked[e.pred] = 1;
            }
            for (const Entry &src : item.state) {
                const LabelId l = src.label;
                // ... nor to a lineage that continues on |c| (or that a quorum refused it)
                if (sc.excluded[l] || sc.budget_checked[l] || sc.marked[l])
                    continue;
                if (budget_look(src, c).over)
                    refusals(c.succ->ch, "loss_budget").push_back(l);
            }
            for (const Entry &e : c.state) {
                sc.marked[e.pred] = 0;
            }
        }
    }

    // ---- branch event: an ambiguity taken or any refusal (a quorum or split-limit
    // stop, an excluded source and a successor refused by the budget each recorded one)
    if (!ambiguous_taken.empty() || !refused.empty()) {
        arm.result.branch_events_total++;
        if (arm.result.branch_events.size() < strategy_.max_branch_events) {
            BranchEvent be;
            be.at_bp = at;
            be.segment = item.segment;
            for (const Cand &c : cands_) {
                if (!c.admissible() || !c.initial_labels)
                    continue;
                be.chars.push_back(c.succ->ch);
                be.labels_per_successor.push_back(c.initial_labels);
            }
            be.ambiguous = ambiguous_taken;
            be.ambiguous.insert(be.ambiguous.end(), ambiguous_over.begin(), ambiguous_over.end());
            std::sort(be.ambiguous.begin(), be.ambiguous.end());
            be.dropped = dropped;
            be.labels_affected = be.ambiguous.size() + be.dropped.size();
            for (BranchEvent::Refusal &r : refused) {
                std::sort(r.labels.begin(), r.labels.end());
                r.labels.erase(std::unique(r.labels.begin(), r.labels.end()), r.labels.end());
            }
            be.refused = std::move(refused);
            arm.result.branch_events.push_back(std::move(be));
        }
    }

    if (!nf) {
        finish_path(arm, item, std::nullopt);
        return std::nullopt;
    }

    // ---- steps
    account_steps(arm, at, nf);

    if (nf == 1) {
        Cand &c = *followed[0];
        commit_entries(arm, item, item.segment, c.state, at, nullptr);
        arm.walk_seq[item.segment].push_back(c.succ->ch);
        arm.result.segments[item.segment].length_bp++;
        record_edge(arm, item.segment, c);
        item.node = c.succ->node;
        std::swap(item.kmer, c.kmer);       // keep the buffers for the next step
        item.ext_bp = at + 1;
        std::swap(item.state, c.state);
        arrival(arm, item);
        arm.next.push_back(std::move(item));
        return std::nullopt;
    }

    // ---- split
    Split split { at, item.segment, {}, !ambiguous_taken.empty(), item.state.size(), {} };
    {
        GrowthBin &b = bin(arm, at);
        b.splits++;
        if (split.ambiguous) {
            b.ambiguous_branches++;
        } else {
            b.divergences++;
        }
    }
    arm.result.segments[item.segment].labels_end = labels_of(item.state);
    std::vector<uint32_t> taken;
    for (size_t j = 0; j < nf; ++j) {
        Cand &c = *followed[j];
        std::vector<LabelId> child_labels = labels_of(c.state);
        size_t child = new_segment(arm, { item.segment }, at, child_labels, child_labels.size());
        arm.result.segments[item.segment].children.push_back(child);
        split.children.push_back(child);
        // the trie view: the labels continuing on this branch (bounded list, true count)
        SplitBranch br;
        br.segment = child;
        br.ch = c.succ->ch;
        br.labels_distinct = child_labels.size();
        br.labels = bounded(arm, child_labels, child_labels.size());
        split.branches.push_back(std::move(br));
        commit_entries(arm, item, child, c.state, at, &taken);
        arm.walk_seq[child].push_back(c.succ->ch);
        arm.result.segments[child].length_bp = 1;
        record_edge(arm, child, c);
        Item ni;
        ni.segment = child;
        ni.node = c.succ->node;
        ni.kmer = std::move(c.kmer);
        ni.ext_bp = at + 1;
        ni.state = std::move(c.state);
        ni.splits = item.splits + 1;
        ni.path_id = arm.next_path_id++;
        arrival(arm, ni);
        arm.next.push_back(std::move(ni));
    }
    arm.result.splits.push_back(std::move(split));
    return std::nullopt;
}

void Walker::check_structure(ArmState &arm, const Item &item, Cand &c) {
    step_kmer(arm.arm, item.kmer, c.succ->ch, &step_);
    if (regime_ != Regime::BASIC && is_hairpin(item.kmer, c.kmer, step_)) {
        c.hairpin = true;
        c.skipped = strategy_.skip_hairpins;
    }
    if (c.skipped)
        return;
    if (seed_nodes_.count(c.succ->node)) {
        c.blocked = true;
        c.block_reason = EndReason::REACHED_SEED;
        return;
    }
    auto [key, rc] = edge_key(step_);
    c.key_lo = key.lo;
    c.key_hi = key.hi;
    c.rc = rc;
    auto it = arm.used_edges.find(key);
    if (it == arm.used_edges.end())
        return;
    mark_ancestors(arm, item.segment);
    // Was the edge taken by this segment or an ancestor, and in which orientation?
    // Either scan the edge's uses or probe (edge, segment) for the path's own
    // segments: the same uses, so the same answer, at the cost of whichever is fewer.
    const SmallVector<EdgeUse> &uses = it->second;
    bool same = false, opposite = false;
    if (uses.size() <= arm.ancestors.size() + 1) {
        arm.result.edge_reuse_probes += uses.size();
        for (const EdgeUse &u : uses) {
            if (!is_ancestor_or_self(arm, u.segment(), item.segment))
                continue;
            (u.rc() == c.rc ? same : opposite) = true;
        }
    } else {
        arm.result.edge_reuse_probes += arm.ancestors.size() + 1;
        auto probe = [&](size_t segment) {
            auto jt = arm.used_by_segment.find(EdgeSegKey{ c.key_lo, c.key_hi, segment });
            if (jt == arm.used_by_segment.end())
                return;
            same |= (jt->second >> c.rc) & 1;
            opposite |= (jt->second >> !c.rc) & 1;
        };
        probe(item.segment);
        for (size_t a : arm.ancestors) {
            probe(a);
        }
    }
    if (opposite) {
        c.blocked = true;
        c.block_reason = EndReason::EDGE_REUSE_RC;
    } else if (same) {
        c.blocked = true;
        c.block_reason = EndReason::EDGE_REUSE;
    }
}

std::optional<EndReason> Walker::cap_check(const ArmState &arm, size_t nf,
                                           size_t remaining_in_level) const {
    if (!nf)
        return std::nullopt;
    if (steps_total_ + nf > strategy_.max_steps)
        return EndReason::MAX_STEPS;
    if (arm.result.output_bp + nf > strategy_.max_output_bp)
        return EndReason::MAX_OUTPUT;
    size_t live_after = remaining_in_level + arm.next.size() + nf;
    if (nf >= 2 && arm.finished_leaves + live_after > strategy_.max_paths)
        return EndReason::MAX_PATHS;
    if (strategy_.on_overflow == Strategy::STOP && live_after > strategy_.max_live_paths)
        return EndReason::MAX_LIVE_PATHS;
    return std::nullopt;
}

void Walker::account_steps(ArmState &arm, uint64_t at, size_t nf) {
    bin(arm, at).steps += nf;
    arm.result.steps += nf;
    arm.result.output_bp += nf;
    steps_total_ += nf;
}

void Walker::record_edge(ArmState &arm, size_t segment, const Cand &c) {
    arm.used_edges[EdgeKey{ c.key_lo, c.key_hi }].push_back(EdgeUse(segment, c.rc));
    arm.used_by_segment[EdgeSegKey{ c.key_lo, c.key_hi, segment }]
        |= static_cast<uint8_t>(1u << c.rc);
}

std::vector<LabelId> Walker::bounded(ArmState &arm, const std::vector<LabelId> &labels,
                                     size_t total) {
    arm.result.max_labels_at_node = std::max(arm.result.max_labels_at_node, total);
    if (total > strategy_.max_labels_per_node || labels.size() > strategy_.max_labels_per_node)
        arm.result.nodes_labels_truncated++;
    if (labels.size() <= strategy_.max_labels_per_node)
        return labels;
    return std::vector<LabelId>(labels.begin(), labels.begin() + strategy_.max_labels_per_node);
}

void Walker::record_present(ArmState &arm, size_t segment, uint64_t at,
                            const LabelRecorder::NodeLabels &nl) {
    // the recorder already cut the list to the cap; this counts the cut
    std::vector<LabelId> labels = bounded(arm, nl.labels, nl.total);
    auto &sets = arm.result.segments[segment].label_sets;
    if (!sets.empty() && sets.back().to_bp == at && sets.back().labels_total == nl.total
            && sets.back().labels == labels) {
        sets.back().to_bp = at + 1;
        return;
    }
    LabelSetRun run;
    run.from_bp = at;
    run.to_bp = at + 1;
    run.labels = std::move(labels);
    run.labels_total = nl.total;
    sets.push_back(std::move(run));
}

std::optional<EndReason> Walker::process_item_annotate(ArmState &arm, Item &item,
                                                       const std::vector<Succ> &succs,
                                                       const LabelRecorder::NodeLabels *present,
                                                       size_t remaining_in_level) {
    const uint64_t at = item.ext_bp;
    const Arm side = arm.arm;
    if (succs.empty()) {
        finish_path(arm, item, EndReason::DEAD_END);
        return std::nullopt;
    }

    cands_.resize(succs.size());
    std::vector<Cand*> followed;
    uint8_t blocked = 0;
    for (size_t i = 0; i < succs.size(); ++i) {
        Cand &c = cands_[i];
        c.reset();
        c.succ = &succs[i];
        c.present = present[i];
        succ_kmer(side, item.kmer, succs[i].ch, &c.kmer);
        check_structure(arm, item, c);
        c.followed = c.admissible();
        if (c.followed) {
            followed.push_back(&c);
        } else if (c.blocked) {
            blocked = std::max(blocked, block_rank(c.block_reason));
        }
    }
    const size_t nf = followed.size();
    if (auto cap = cap_check(arm, nf, remaining_in_level))
        return cap;

    // ---- commit: every successor not followed, and every followed hairpin, is reported
    // with the labels present there — the only way a reader can tell what a blocked
    // continuation would have carried. The recorder already cut the list to the cap;
    // for a successor that is NOT followed this event is the only place its list
    // appears, so the cut is counted here (a followed hairpin's node is recorded and
    // counted by record_present below, so its event only repeats the list).
    for (const Cand &c : cands_) {
        if (c.followed && !c.hairpin)
            continue;
        Event ev;
        ev.at_bp = at;
        ev.ch = c.succ->ch;
        ev.labels = c.followed ? c.present.labels
                               : bounded(arm, c.present.labels, c.present.total);
        ev.labels_total = c.present.total;
        if (c.skipped) {
            ev.type = EventType::HAIRPIN;
        } else if (c.blocked) {
            ev.type = EventType::BLOCKED;
            ev.reason = c.block_reason;
            bin(arm, at).blocked_repeat++;
        } else {
            ev.type = EventType::HAIRPIN;
            ev.text = "followed";
        }
        arm.result.segments[item.segment].events.push_back(std::move(ev));
    }

    if (!nf) {
        // a structural end: the path's reason (no label lineage ends in this mode);
        // a successor that was only a skipped hairpin is a dead end with its event
        finish_path(arm, item, blocked ? block_reason_of_rank(blocked) : EndReason::DEAD_END);
        return std::nullopt;
    }

    account_steps(arm, at, nf);

    if (nf == 1) {
        Cand &c = *followed[0];
        arm.walk_seq[item.segment].push_back(c.succ->ch);
        arm.result.segments[item.segment].length_bp++;
        record_present(arm, item.segment, at, c.present);
        record_edge(arm, item.segment, c);
        item.node = c.succ->node;
        std::swap(item.kmer, c.kmer);
        item.ext_bp = at + 1;
        item.present = std::move(c.present.labels);
        item.present_total = c.present.total;
        arrival(arm, item);
        arm.next.push_back(std::move(item));
        return std::nullopt;
    }

    // ---- split: a divergence of walks; "ambiguous" is a lineage notion, and there
    // is no lineage here
    Split split { at, item.segment, {}, false, item.present_total, {} };
    {
        GrowthBin &b = bin(arm, at);
        b.splits++;
        b.divergences++;
    }
    arm.result.segments[item.segment].labels_end = item.present;
    for (size_t j = 0; j < nf; ++j) {
        Cand &c = *followed[j];
        size_t child = new_segment(arm, { item.segment }, at, c.present.labels, c.present.total);
        arm.result.segments[item.segment].children.push_back(child);
        split.children.push_back(child);
        SplitBranch br;
        br.segment = child;
        br.ch = c.succ->ch;
        br.labels_distinct = c.present.total;
        br.labels = c.present.labels;    // already cut to the cap by the recorder
        split.branches.push_back(std::move(br));
        arm.walk_seq[child].push_back(c.succ->ch);
        arm.result.segments[child].length_bp = 1;
        record_present(arm, child, at, c.present);
        record_edge(arm, child, c);
        Item ni;
        ni.segment = child;
        ni.node = c.succ->node;
        ni.kmer = std::move(c.kmer);
        ni.ext_bp = at + 1;
        ni.splits = item.splits + 1;
        ni.path_id = arm.next_path_id++;
        ni.present = std::move(c.present.labels);
        ni.present_total = c.present.total;
        arrival(arm, ni);
        arm.next.push_back(std::move(ni));
    }
    arm.result.splits.push_back(std::move(split));
    return std::nullopt;
}


/********************************* outputs **********************************/

void Walker::finalize(ArmState &arm) {
    ArmResult &res = arm.result;
    // every walk up to the first unexpanded head's depth is present; a complete arm
    // (no boundary) is complete to the radius
    res.complete_to_bp = std::min(arm.boundary, strategy_.max_extension_bp);
    assert((res.status == ArmResult::COMPLETE) == (res.complete_to_bp == strategy_.max_extension_bp)
           || !res.requested);
    // merge_level unites the edge histories of the routes it joins, so under merging
    // the walks present are those admissible under the united history (§6.10)
    res.completeness_scope = strategy_.merge_reconverge ? "united_history" : "per_path";
    for (size_t s = 0; s < res.segments.size(); ++s) {
        Segment &seg = res.segments[s];
        if (strategy_.sequences) {
            seg.sequence = arm.walk_seq[s];
            if (arm.arm == Arm::LEFT)
                std::reverse(seg.sequence.begin(), seg.sequence.end());
        }
        std::stable_sort(seg.events.begin(), seg.events.end(),
                         [](const Event &a, const Event &b) { return a.at_bp < b.at_bp; });
    }
    for (size_t s = 0; s < res.segments.size(); ++s) {
        const LeafInfo &leaf = arm.leaves[s];
        if (!leaf.is_leaf)
            continue;
        PathResult path;
        path.id = res.paths.size();
        size_t cur = s;
        while (true) {
            path.segments.push_back(cur);
            if (res.segments[cur].parents.empty())
                break;
            cur = res.segments[cur].parents[0];
        }
        std::reverse(path.segments.begin(), path.segments.end());
        path.length_bp = res.segments[s].from_bp + res.segments[s].length_bp;
        path.end_reasons = leaf.end_reasons;
        path.path_reason = leaf.path_reason;
        path.end_labels = leaf.end_labels;
        path.continuation = leaf.continuation;
        res.paths.push_back(std::move(path));
    }
}

void Walker::summarize() {
    result_.label_summary.assign(result_.label_dict.size(), {});
    for (ArmState &arm : arms_) {
        size_t a = static_cast<size_t>(arm.arm);
        const auto &runs = arm.result.runs;
        // lineage root of each run: prev_run always points at an earlier run, so one
        // forward pass resolves every switch chain
        std::vector<LabelId> root(runs.size());
        for (size_t r = 0; r < runs.size(); ++r) {
            const LabelRun &run = runs[r];
            if (run.entered_by_switch && run.prev_run != UINT32_MAX) {
                assert(run.prev_run < r);
                root[r] = root[run.prev_run];
            } else {
                root[r] = run.label;
            }
        }
        for (size_t r = 0; r < runs.size(); ++r) {
            const LabelRun &run = runs[r];
            LabelArmSummary &own = result_.label_summary[run.label][a];
            own.runs.push_back(r);
            // direct_bp is measured along the label's OWN route, which is why a merge
            // does not clamp it: a reconvergence joins paths at the same node, so the
            // bases from the merge onward are identical on every route into it, and
            // every k-mer of from_bp..to_bp carries this label on some route. That is
            // label-consistent route support, NOT a contiguous occurrence in the
            // label's sequence (see LabelArmSummary). What a merge invalidates is
            // attributing the label to ANOTHER path's spelled prefix, reported per path
            // by LabelEnd::route_bp and per run by LabelRun::route_bp, not here.
            if (!run.entered_by_switch && run.from_bp == 0)
                own.direct_bp = std::max(own.direct_bp, run.to_bp);
            LabelArmSummary &lineage = result_.label_summary[root[r]][a];
            lineage.reach_bp = std::max(lineage.reach_bp, run.to_bp);
        }
        for (const Segment &seg : arm.result.segments) {
            for (const Event &ev : seg.events) {
                if (ev.type == EventType::SWITCH)
                    result_.label_summary[ev.to][a].reentries++;
            }
        }
    }
}

void Walker::summarize_annotate() {
    result_.label_summary.assign(result_.label_dict.size(), {});
    for (ArmState &arm : arms_) {
        const size_t a = static_cast<size_t>(arm.arm);
        const auto &segs = arm.result.segments;
        if (segs.empty())
            continue;
        // reach_bp: the furthest position at which the label was recorded on any path
        for (const Segment &seg : segs) {
            for (const LabelSetRun &run : seg.label_sets) {
                for (LabelId l : run.labels) {
                    LabelArmSummary &s = result_.label_summary[l][a];
                    s.reach_bp = std::max(s.reach_bp, run.to_bp);
                }
            }
        }
        // direct_bp: continuous presence from the seed boundary along SOME route.
        // Segments are created parents-first (new_segment appends), so one pass in id
        // order that enters each segment with the UNION of its parents' surviving sets
        // visits every segment once; walking the routes instead was exponential on a
        // merged DAG (N nested diamonds: O(N) segments, 2^N routes). The union is exact
        // for an existential claim: a label survives a segment iff it is present on its
        // runs, whichever route brought it in. Exact when no per-node list was cut
        // (nodes_labels_truncated == 0), a lower bound otherwise.
        std::vector<std::vector<LabelId>> alive_end(segs.size());
        std::vector<LabelId> alive, still;
        for (size_t s = 0; s < segs.size(); ++s) {
            const Segment &seg = segs[s];
            alive.clear();
            if (seg.parents.empty()) {
                alive = seg.labels_start;
            } else {
                for (size_t p : seg.parents) {
                    assert(p < s);
                    still.clear();
                    std::set_union(alive.begin(), alive.end(), alive_end[p].begin(),
                                   alive_end[p].end(), std::back_inserter(still));
                    alive.swap(still);
                }
            }
            for (const LabelSetRun &run : seg.label_sets) {
                if (alive.empty())
                    break;
                still.clear();
                std::set_intersection(alive.begin(), alive.end(), run.labels.begin(),
                                      run.labels.end(), std::back_inserter(still));
                for (LabelId l : alive) {
                    if (!std::binary_search(still.begin(), still.end(), l)) {
                        LabelArmSummary &x = result_.label_summary[l][a];
                        x.direct_bp = std::max(x.direct_bp, run.from_bp);
                    }
                }
                alive.swap(still);
            }
            const uint64_t end = seg.from_bp + seg.length_bp;
            for (LabelId l : alive) {
                LabelArmSummary &x = result_.label_summary[l][a];
                x.direct_bp = std::max(x.direct_bp, end);
            }
            alive_end[s] = alive;
        }
    }
}

SeedResult Walker::run() {
    validate_seed();
    init_edge_coding();
    if (annotate_) {
        const LabelKind kind = strategy_.seed_label_kind.value_or(
                oracle_.coord_to_header() ? LabelKind::HEADER : LabelKind::COLUMN);
        recorder_ = std::make_unique<LabelRecorder>(oracle_, kind, strategy_.max_labels_per_node);
    } else {
        query_ = std::make_unique<LabelQuery>(oracle_, result_.label_dict, trace_);
    }
    scratch_.init(result_.label_dict.size());
    label_stamp_.assign(result_.label_dict.size(), 0);

    arms_[static_cast<size_t>(Arm::LEFT)].arm = Arm::LEFT;
    arms_[static_cast<size_t>(Arm::RIGHT)].arm = Arm::RIGHT;
    for (ArmState &arm : arms_) {
        init_arm(arm);
    }

    uint64_t depth = 0;
    while (!seed_stopped_) {
        bool any = false;
        for (const ArmState &arm : arms_) {
            any |= !arm.frontier.empty();
        }
        if (!any)
            break;
        if (depth > 0 && time_exceeded()) {
            for (ArmState &arm : arms_) {
                stop_frontier(arm, EndReason::TIME_BUDGET);
            }
            break;
        }
        for (ArmState &arm : arms_) {
            if (!arm.frontier.empty() && !arm.stopped)
                run_level(arm, depth);
            if (seed_stopped_)
                break;
        }
        ++depth;
    }

    if (annotate_) {
        // the labels met along either arm, in the order first seen
        result_.label_dict = recorder_->labels();
    }
    for (ArmState &arm : arms_) {
        finalize(arm);
        arm.result.pair_evaluations = 0;
    }
    arms_[static_cast<size_t>(Arm::RIGHT)].result.pair_evaluations = pair_evaluations_;
    if (annotate_) {
        summarize_annotate();
    } else {
        summarize();
    }

    result_.arms[static_cast<size_t>(Arm::LEFT)] = std::move(arms_[0].result);
    result_.arms[static_cast<size_t>(Arm::RIGHT)] = std::move(arms_[1].result);
    result_.annotation_counters = oracle_.counters();
    result_.access_path = annotate_ ? recorder_->access_path() : query_->access_path();
    result_.elapsed_seconds = timer_.elapsed();
    return std::move(result_);
}

} // namespace


void validate_strategy(const Strategy &st, const LabelChangeCost &cost) {
    if (st.loss_budget < 0)
        throw std::invalid_argument("Negative loss budget");
    if (st.min_successor_fraction < 0 || st.min_successor_fraction > 1)
        throw std::invalid_argument("min_successor_fraction must be in [0, 1]");
    if (!st.max_live_paths)
        throw std::invalid_argument("max_live_paths must be positive");
    if (!st.max_labels_per_node)
        throw std::invalid_argument("max_labels_per_node must be positive");
    // A conflicting knob is refused, never overridden: a caller who asked for the
    // exhaustive trie and got a pruned one would have no way to tell.
    auto conflict = [](const char *knob, const char *required, const std::string &why) {
        throw std::invalid_argument(std::string(knob) + " " + why + "; set it to " + required
                                    + " or omit it");
    };
    if (st.exhaustive) {
        const std::string why = "conflicts with `exhaustive` (every admissible walk must be "
                                "explored and kept)";
        if (st.max_label_branches != Strategy::kUnlimited)
            conflict("branching.max_label_branches", "\"unlimited\"", why);
        if (st.merge_reconverge) {
            conflict("branching.on_reconverge", "\"keep\"",
                     "conflicts with `exhaustive` (an exhaustive result is a trie of walks, "
                     "not a DAG)");
        }
        if (st.max_splits_per_path != Strategy::kUnlimited)
            conflict("branching.max_splits_per_path", "\"unlimited\"", why);
        if (st.min_successor_labels > 1)
            conflict("branching.min_successor_labels", "1", why + " and a quorum prunes");
        if (st.min_successor_fraction > 0)
            conflict("branching.min_successor_fraction", "0", why + " and a quorum prunes");
        if (st.min_live_labels > 1)
            conflict("bounds.min_live_labels", "1", why + " and a quorum prunes");
        if (st.on_overflow == Strategy::BEAM) {
            conflict("frontier.on_overflow", "\"stop\"",
                     "conflicts with `exhaustive` (a beam prunes; a cap has to stop the arm and "
                     "report complete_to_bp instead)");
        }
        // derive() sorts the switch sources of a pairwise cost and keeps the cheapest
        // max_switch_sources: a target label reachable only from a cut source is not
        // entered, and if it was the only label on that successor the walk into it is
        // pruned with no path-level reason. CONSTANT and FORBID never cut.
        if (cost.model() == LabelChangeCost::TABLE
                && st.max_switch_sources != Strategy::kUnlimited) {
            conflict("labels.max_switch_sources", "\"unlimited\"",
                     "conflicts with `exhaustive` under a \"table\" change cost (a cut source "
                     "list can leave a target label unentered and so prune a label-consistent "
                     "walk)");
        }
    }
    if (st.label_mode == LabelMode::ANNOTATE) {
        const std::string why = "does not apply in labels.mode \"annotate\" (labels are recorded, "
                                "not filtered: there is no permitted set)";
        if (!st.extra.empty())
            conflict("labels.extra", "[]", why);
        if (cost.finite())
            conflict("labels.change_cost", "{\"model\": \"forbid\"}", why);
        if (st.loss_budget > 0)
            conflict("labels.loss_budget", "0", why);
        if (st.min_successor_labels > 1)
            conflict("branching.min_successor_labels", "1", why);
        if (st.min_successor_fraction > 0)
            conflict("branching.min_successor_fraction", "0", why);
        if (st.min_live_labels > 1)
            conflict("bounds.min_live_labels", "1", why);
        if (st.support == Support::TRACE) {
            conflict("support", "\"kmer\"",
                     "cannot be \"trace\" in labels.mode \"annotate\": a coordinate trace is "
                     "per-label evidence and needs a permitted set, which this mode has none of");
        }
    }
}

std::string walk_rule_statement(const Strategy &st, const LabelOracle &oracle) {
    const size_t k = oracle.get_k();
    const bool stranded = oracle.regime() != Regime::BASIC;
    // the edge coding of Walker::init_edge_coding(): a (k+1)-mer is packed exactly
    // when it fits 64 bits, otherwise identified by a 128-bit FNV-1a pair
    size_t symbols = 0;
    for (char c : oracle.graph().alphabet()) {
        symbols += c != kSentinel;
    }
    size_t bits = 1;
    while ((size_t(1) << bits) < symbols) {
        ++bits;
    }
    const bool packable = (k + 1) * bits <= 64;

    std::string s = "every walk of at most complete_to_bp bases from the seed boundary is present "
                    "when it uses no (k+1)-mer edge twice on its own path (";
    s += stranded
        ? "edges are compared as canonical (k+1)-mers, each use keeping its orientation"
        : "edges are compared as (k+1)-mers";
    if (!packable) {
        s += ", identified by a 128-bit FNV-1a hash of the (k+1)-mer, so a hash collision "
             "blocks a step as a reuse";
    }
    s += "), enters no seed node";
    if (stranded)
        s += " (of either strand)";
    s += ", and ";
    if (!stranded) {
        // check_structure() tests for hairpins in stranded regimes only: a basic graph
        // holds one strand, so a self-reverse-complementary step is an ordinary step there
        s += "takes every step regardless of strand (a single-strand graph has no hairpins); ";
    } else {
        std::string hairpin = "a self-reverse-complementary (k+1)-mer";
        if (k % 2 == 0)
            hairpin += " or a step into or out of a self-reverse-complementary k-mer node (even k)";
        s += st.skip_hairpins
            ? "takes no hairpin step (" + hairpin + "); "
            : "may take hairpin steps (" + hairpin + "), which are flagged; ";
    }
    if (st.merge_reconverge) {
        // merge_level unites the edge histories of the routes it joins (conservative), so
        // the walks present are those admissible under the UNITED history, fewer than the
        // per-path set the rule above describes; the certificate must say so
        s += "after a reconvergence the merged walk is also barred from every edge used by "
             "any route merged into it (edge histories are united under on_reconverge: "
             "merge, so fewer walks are present than per path; on_reconverge: keep gives "
             "the per-path set); ";
    }
    s += st.label_mode == LabelMode::ANNOTATE
        ? "every such walk is followed and the labels present at its nodes are recorded"
        : "it is followed while some permitted label supports every node of it under the "
          "strategy's label rules";
    return s;
}

SeedResult traverse_seed(LabelOracle &oracle,
                         const Seed &seed,
                         const Strategy &strategy,
                         const LabelChangeCost &cost,
                         const std::string &release_id) {
    validate_strategy(strategy, cost);
    Walker walker(oracle, seed, strategy, cost, release_id);
    return walker.run();
}

std::string spell_path(const ArmResult &arm, const PathResult &path) {
    std::string out;
    if (arm.arm == Arm::RIGHT) {
        for (size_t s : path.segments) {
            out += arm.segments[s].sequence;
        }
    } else {
        for (auto it = path.segments.rbegin(); it != path.segments.rend(); ++it) {
            out += arm.segments[*it].sequence;
        }
    }
    return out;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
