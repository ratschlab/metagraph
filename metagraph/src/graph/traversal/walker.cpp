#include "walker.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <queue>
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
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
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

const char* to_string(SeedDerivationError::Cause cause) {
    switch (cause) {
        case SeedDerivationError::NO_CARRIER: return "no_carrier";
        case SeedDerivationError::NO_TRACE_CARRIER: return "no_trace_carrier";
        case SeedDerivationError::TOO_WIDE: return "too_wide";
        case SeedDerivationError::TIME_BUDGET: return "time_budget";
        case SeedDerivationError::AMBIGUOUS_HEADER: return "ambiguous_header";
        case SeedDerivationError::OVER_SEED_LABEL_CAP: return "over_seed_label_cap";
        case SeedDerivationError::UNREPRESENTABLE_LABEL_NAME: return "unrepresentable_label_name";
    }
    return "unknown";
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
 *   implemented in this increment; the fields are kept and stay zero (spec §12), and
 *   validate_strategy() refuses a non-zero window instead of ignoring it.
 * - EndReason has no values for `switched`, `superseded`, `minority`,
 *   `below_min_labels`, `split_limit`, `switch_sources` and `hairpin`: a source
 *   whose lineage only continues under another name ends its run with LABEL_LOST
 *   plus a SWITCH event; quorum and split-limit stops use BRANCH with the text in
 *   Event::text; a label whose only continuation is a skipped hairpin ends with
 *   DEAD_END, text "hairpin"; a label cut from the switch sources by
 *   max_switch_sources ends with LABEL_LOST, text "switch_sources". A cut source that
 *   does not end leaves no such trace, so every cut that may have changed a step is
 *   also counted (ArmResult::switch_sources_cut).
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
// Under a §14 budget the level's annotation keys are fetched in calls of at most this
// many (fewer when rows are wide or the work budget is near, fetch_chunk()): a fetch call is
// the unit of charging, so its rows are charged as work when it returns them and what it
// holds is observed before the comparison after it. Without a budget the level is one
// LOGICAL call, as it always was (the cache, and with it the direct_reads counter, depends on
// the batching). A call the deadline may fall into is also decoded in time-sized chunks with
// the deadline checked between them (pass 5, LabelOracle::pacer; a call far from it is one
// piece, as splitting costs what its rows share): a call's counting and cache decisions stay
// whole, so chunking changes nothing it returns, and the uninterruptible unit is one chunk
// (at least one row), not the call.
constexpr size_t kFetchChunk = 8192;

// A budget or the deadline ran out while a head was being PLANNED or a level fetched
// (DESIGN-traverse-graphlet.md §14). Thrown only before anything of the head or level is
// committed, so the caller censors from that head on, exactly as at a cap. Never escapes
// the walker.
struct BudgetTrip {
    ResourceStop::Resource resource;
    double demand;                  // work units, elapsed ms, or bytes
    // annotation_decode: a budget-aware annotation read did not fit (stage 3), with what the
    // walk held then (bytes; -1: the account, as for a refused head) and whether a test hook
    // (WalkerHooks::deny_decode) refused it rather than the budget
    const char *phase = "traversal";
    double used = -1;
    bool injected = false;
    // what did not fit (ResourceStop's cause fields), for the statements
    ResourceStop detail = ResourceStop();
};

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
    // The memory reserved for this head (§14, modelled bytes): |reserve| holds the head
    // itself and what ending it costs (its label ends, leaf, continuation, path chain and
    // their delivery), |merge_reserve| what a reconvergence merge at the end of its level
    // costs (released or consumed by merge_level). Whatever ends the head — its own step,
    // a cap, the radius, a beam or a budget — fits in what it holds, so a prefix that
    // was admitted can always be finished and delivered.
    uint64_t reserve = 0;
    uint64_t merge_reserve = 0;
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
    // derive() cut a source with a finite pair cost into a target of this successor:
    // the cut may have changed a loss or an entry (ArmResult::switch_sources_cut)
    bool cut_reaches = false;
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
        cut_reaches = false;
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
    // depth of the last branch event emitted: levels are synchronous, so it never
    // decreases, which is what makes the stored prefix a depth boundary
    uint64_t last_branch_event_bp = 0;
    // ancestor marks of |marked_segment| (ancestors never change after creation), and
    // the marked segments as a list: the proper ancestors, in no particular order
    std::vector<uint32_t> visit_mark;
    uint32_t visit_epoch = 0;
    size_t marked_segment = SIZE_MAX;
    std::vector<size_t> visit_stack;
    std::vector<size_t> ancestors;
    // per segment, the length of its first-parent chain (root: 1): what a path ending
    // there costs to spell out in JSON paths[].segments, charged per expansion (§14)
    std::vector<uint32_t> chain_len;
    // growth bins already charged to the memory account (bins are created lazily)
    size_t bins_charged = 0;
    // work units not visible as a counter (annotation keys and entries, derivation
    // scans, trace coordinates); ArmResult::work_units adds the weighted counters
    uint64_t work_extra = 0;
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

// What one source (an entry of the head's state) does at a committed step: its lineage
// continues on a followed successor, ends silently (it continues only under other names,
// the switch events are written on commit), or ends with a label_end event.
struct SourcePlan {
    enum Kind : uint8_t { CONTINUES, SILENT_END, LABEL_END };
    Kind kind = CONTINUES;
    EndReason reason = EndReason::DEAD_END;
    const char *text = "";
    double needed = 0;
};

// The plan of one head (DESIGN-traverse-graphlet.md §14, "atomic commit per head"): every
// decision of the step, computed without touching the result, so that the head can be
// refused (censored with resource_limit) after planning and before anything of it is
// written. Reused across heads (member of the walker) to keep its buffers.
struct HeadPlan {
    std::vector<Cand*> followed;
    std::vector<LabelId> ambiguous_over, ambiguous_taken, dropped;
    // the explicit per-successor refusals of the step (BranchEvent::refused), one per
    // (successor, cause); the labels are made distinct and ascending when the event is
    // emitted
    std::vector<BranchEvent::Refusal> refused;
    std::vector<SourcePlan> sources;     // parallel to the head's state
    // The re-minimisation rounds the step needed (the rounds after the first). Kept here
    // and counted when the step is committed (or a cap stops it, as before stage 2), not
    // while planning: max_reminimisation_rounds > 0 is what states greedy_losses, and a
    // head that is not admitted decided nothing greedily (review of stage 2, finding 1).
    size_t reminimisations = 0;
    // what committing the plan costs (§14, modelled bytes): the objects it adds for good,
    // and per followed successor the new head's reservations
    uint64_t committed = 0;
    std::vector<uint64_t> child_reserve, child_merge_reserve;
    // the growth bins the children's level writes into, charged with the first head of a
    // level that creates children (the arm's bins charged up to this count on commit), so
    // that no level starts beyond the budget by its bins
    size_t bins_needed = 0;
    // what the commit will create, checked against it in debug builds
    size_t new_segments = 0, new_runs = 0, new_events = 0;

    void clear() {
        followed.clear();
        ambiguous_over.clear();
        ambiguous_taken.clear();
        dropped.clear();
        refused.clear();
        sources.clear();
        reminimisations = 0;
        committed = 0;
        child_reserve.clear();
        child_merge_reserve.clear();
        bins_needed = 0;
        new_segments = new_runs = new_events = 0;
    }
    uint64_t reserved() const {
        uint64_t total = 0;
        for (size_t j = 0; j < child_reserve.size(); ++j) {
            total += child_reserve[j] + child_merge_reserve[j];
        }
        return total;
    }
};

/**
 * The memory model of §14: bytes per object, retained plus delivered, each a fixed upper
 * bound so that admission is deterministic (never a measurement). A vector grows by
 * doubling, so an element is charged twice its size; an event once more for the stable
 * sort in finalisation (bounded by any one segment's events).
 */
struct CostModel {
    uint64_t segment = 0, seg_label = 0, base = 0, step = 0, run = 0, event = 0,
             event_label = 0, needed = 0, leaf = 0, leaf_label = 0, chain_entry = 0,
             cont_base = 0, split = 0, split_branch = 0, bevent = 0, bevent_entry = 0,
             refusal = 0, presence_run = 0, bin = 0, item = 0, entry = 0, coord = 0,
             present = 0, label = 0, label_name = 0, merge_parent = 0;
};

// the labels of refusal (ch, cause) of |plan|, created when first needed; a step has a
// few of them at most (successors x causes). The reference is used before the next call.
std::vector<LabelId>& refusals_of(HeadPlan &plan, char ch, const char *cause) {
    for (BranchEvent::Refusal &r : plan.refused) {
        if (r.ch == ch && std::string_view(r.cause) == cause)
            return r.labels;
    }
    plan.refused.push_back({ ch, cause, {} });
    return plan.refused.back().labels;
}

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
           const std::string &release_id,
           const WalkerHooks *hooks,
           const AttemptControl *control)
          : oracle_(oracle), graph_(oracle.graph()), seed_(seed), strategy_(strategy),
            cost_(cost), release_id_(release_id), hooks_(hooks),
            control_(control), k_(oracle.get_k()),
            regime_(oracle.regime()), canonical_(oracle.canonical()),
            nfc_(oracle.node_first_cache()),
            trace_(strategy.support == Support::TRACE),
            annotate_(strategy.label_mode == LabelMode::ANNOTATE),
            path_cache_scope_(oracle, strategy.max_memory_bytes > 0) {}

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
    // an index-supplied name as a failure message echoes it; |where| identifies it
    std::string echoed(const std::string &name, const std::string &where) const;
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
    // sets cap_demand_ when a cap trips
    std::optional<EndReason> cap_check(const ArmState &arm, size_t nf,
                                       size_t remaining_in_level);
    // ADMIT (§14): whether the planned head may be committed; false censors it with
    // resource_limit (sets cap_demand_)
    bool admit(const ArmState &arm, const Item &item, size_t followed);
    // the per-source decisions of a step that passed the caps (PLAN-B): which lineages
    // continue, end silently or end with a reason, and the refusals; mutates nothing
    // but scratch and the spent-work counters
    void plan_outcomes(ArmState &arm, const Item &item);
    // COMMIT of a constrained step whose plan was admitted; cannot fail
    void commit_item(ArmState &arm, Item &item, const std::vector<Succ> &succs);
    // the plan's re-minimisation rounds into the arm's counters (commit, or a cap)
    void count_reminimisations(ArmState &arm) {
        arm.result.reminimisation_rounds += plan_.reminimisations;
        arm.result.max_reminimisation_rounds = std::max(arm.result.max_reminimisation_rounds,
                                                        plan_.reminimisations);
    }
    // COST of the planned step (HeadPlan::committed, the children's reservations)
    void plan_cost(const ArmState &arm, const Item &item);
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
    // after a resource stop: keep only the dictionary labels the result records
    void compact_dictionary();
    void summarize_annotate();

    // ---- label state
    // |*truncated|, |*cut|: whether max_switch_sources cut the source list and the key
    // of the last source kept; |*cut_reaches|: whether a cut source had a finite pair
    // cost into one of |targets| (other than its own label). The pairs priced are
    // charged to |arm|, the arm whose head is being processed: the counters (and the
    // greedy_losses limitation read off them) are per arm.
    void derive(ArmState &arm, const State &sigma, const std::vector<Target> &targets,
                const std::vector<uint8_t> &excluded, State *out,
                bool *truncated, SourceKey *cut, bool *cut_reaches);
    // a finite cost(from -> t) for some target t != from
    bool reaches_a_target(LabelId from, const std::vector<Target> &targets) const;
    SourceKey source_key(const Entry &e) const {
        const LabelRef &ref = result_.label_dict[e.label];
        return { e.loss, e.branches, ref.column, ref.seq_id };
    }
    // |split|: the state is one child's of a split, whose runs already continued by an
    // earlier child (taken_run()) are cloned so that every run belongs to exactly one path
    void commit_entries(ArmState &arm, const Item &item, size_t target_segment,
                        State &state, uint64_t at, bool split);
    // The runs continued by the children of one split so far, as stamps indexed by run id:
    // O(1) per entry. A scan of the runs taken so far made a split O(|σ|²), once in its
    // plan and once in its commit — on a 2,876-label locus most of the walk's time, with
    // no budget set (review of stage 2, finding 4). begin_split() starts a split;
    // taken_run() marks |run| and tells whether an earlier child of the split took it.
    void begin_split(const ArmState &arm) {
        if (run_stamp_.size() < arm.result.runs.size())
            run_stamp_.resize(arm.result.runs.size(), 0);
        if (++run_epoch_ == 0) {
            std::fill(run_stamp_.begin(), run_stamp_.end(), 0);
            run_epoch_ = 1;
        }
    }
    bool taken_run(uint32_t run) {
        // a split continues runs that existed when it began (switch-ins get new runs)
        assert(run < run_stamp_.size());
        if (run_stamp_[run] == run_epoch_)
            return true;
        run_stamp_[run] = run_epoch_;
        return false;
    }

    // ---- ends
    // |segment|: where the run ends (LabelRun::segment), the head's own segment
    void end_run(ArmState &arm, const Entry &e, uint64_t at, EndReason reason, size_t segment);
    void end_label(ArmState &arm, const Item &item, const Entry &e, EndReason reason,
                   size_t structural, const char *text, double needed = 0);
    void finish_path(ArmState &arm, Item &item, std::optional<EndReason> path_reason);
    void censor_item(ArmState &arm, Item &item, EndReason reason);
    Continuation make_continuation(ArmState &arm, const Item &item);

    // ---- bookkeeping
    GrowthBin& bin(ArmState &arm, uint64_t bp);
    // every event of a segment is written here (counted for the plan's debug check)
    void push_event(ArmState &arm, size_t segment, Event &&ev) {
        arm.result.segments[segment].events.push_back(std::move(ev));
        events_written_++;
    }
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

    // ---- budgets (DESIGN-traverse-graphlet.md §14)
    void init_budgets();
    // under a memory budget: the seed's row-diff path cache from depth 1 on (PathCacheScope)
    void enable_path_cache();
    uint64_t accounted() const { return base_ + committed_ + reserved_; }
    void note_peak() { peak_ = std::max(peak_, accounted()); }
    // a head's own memory, its stop and its merge reservations (see Item::reserve)
    uint64_t item_bytes(const State &state, size_t present) const;
    uint64_t stop_bytes(size_t state, size_t present, uint64_t ext, uint32_t chain) const;
    uint64_t merge_bytes(size_t labels) const;
    uint64_t stop_bytes(const ArmState &arm, const Item &item) const {
        return stop_bytes(item.state.size(), item.present.size(), item.ext_bp,
                          arm.chain_len[item.segment]);
    }
    // the first stop arrival() would report for a head entering |node| at |ext| (a probe)
    bool would_revisit(const ArmState &arm, node_index node, uint64_t ext, bool revisiting) const;
    // the plan's reservations for a new head (constrain: |state|, annotate: |present|)
    void plan_child(const State &state, size_t present, uint64_t ext, uint32_t chain);
    // book a committed plan: release |held| (the head's reservations), add the rest
    void settle(ArmState &arm, uint64_t held);
    // end a head (cap, radius, beam, budget): commit its stop, release its reservations
    void release_head(const ArmState &arm, Item &item);
    // the growth bins the level at |depth| can touch (its own depth and the next)
    size_t bins_needed(uint64_t depth) const;
    // the plan's share of the bins of its children's level (HeadPlan::bins_needed)
    void plan_bins(const ArmState &arm, uint64_t child_depth);
    // charge the bins a level can touch before it runs (only the roots' level is not
    // charged by an admission)
    void charge_bins(ArmState &arm, uint64_t depth);
    void charge_dictionary();
    // charged work units of an arm (ArmResult::work_units) / of the seed: its seed phase
    // and both arms
    uint64_t work_of(const ArmState &arm) const;
    uint64_t work_used() const { return seed_work_ + work_of(arms_[0]) + work_of(arms_[1]); }
    // The seed phase's work (§14, "rows decoded" count for the locus): the rows read to
    // validate the seed or derive its permitted set, charged as the walk charges its
    // reads, so that a work budget bounds the seed phase too (review of stage 2, finding
    // 6). Over the budget it fails the seed: nothing has been walked yet.
    void charge_seed(uint64_t units);
    std::vector<LabelQuery::NodeHits> fetch_seed_hits(LabelQuery &query,
                                                      const std::vector<node_index> &keys);
    // a budget does not hold the seed itself: throws SeedBudgetError
    [[noreturn]] void fail_seed(ResourceStop::Resource resource, double used, double demand,
                                const std::string &what, const char *phase = "traversal",
                                bool injected = false, const ResourceStop *detail = nullptr);
    // The memory budget does not hold the seed's depth-0 state (the dictionary and both arms'
    // roots with what ending and delivering them costs): fails the seed with |need| bytes.
    // |unread|: the arm whose root row was not read — the state before it already reached the
    // budget, or the |labels| new labels it names (|names| bytes) do not fit beside it — so
    // that |need| is only the least the state needs, and is stated so
    [[noreturn]] void fail_depth0(uint64_t need, const Arm *unread = nullptr,
                                  uint64_t labels = 0, uint64_t names = 0);
    // the work budget's comparison (after every charge: what was charged since the last
    // comparison is the largest_charge_ candidate) and the deadline's (|force|: now,
    // otherwise once the interval has passed, §14), and with the deadline a stop from outside
    // the walk (external_stop). Throws BudgetTrip.
    void checkpoint(bool force);
    // A stop requested from outside the walk (AttemptControl::poll): NONE without a control.
    // A client that is gone throws AttemptAborted here, abandoning the walk.
    ExternalStop external_stop();
    // the seed phase's poll: a cancelled attempt, or one past its bound, fails the seed
    // (SeedBudgetError) — nothing has been walked yet, so no partial result exists
    void seed_external_stop();
    // the attempt's elapsed milliseconds, rounded up: what an external stop states
    double attempt_ms() const;
    // the seed's usage for the caller (AttemptControl::meter), however the walk ends
    void write_meter() const noexcept;
    // what the trip costs the result: the resource stop (the first of the seed) and the
    // cap's demand in the knob's unit; returns the end reason of the censored heads
    EndReason note_stop(const ArmState &arm, const Item *head, ResourceStop::Resource resource,
                        double demand, bool injected = false, const char *phase = "traversal",
                        double used = -1, const ResourceStop *detail = nullptr);
    // the trip of a level's budget-aware read that |refusal| refused (stage 3), with its
    // cause: the row with its dependency rows, or (annotate) the labels it would name first;
    // what the level held then is observed first (memory_bound_soft)
    BudgetTrip read_trip(const FetchRefusal &refusal);
    // Throws the trip of a level whose own lists (keys, successors, its fetch's vectors), held
    // beside the account, already leave none of the budget for its budget-aware read: their
    // cause, not a row's (review of stage 3, F6)
    void level_lists_trip();
    // what a recorded (annotate) dictionary label named |name| costs the account: label_bytes()'s
    // model, read from the name in place
    uint64_t recorded_label_bytes(std::string_view name) const;
    // how a seed-phase or root read that |refusal| refused did not fit (none when a test hook
    // injected the refusal): what its row needed, more than what was left |beside| ...
    static std::string refused_read_text(const FetchRefusal &refusal, bool injected,
                                         const std::string &beside);
    // a level's annotation rows (in calls under a budget, observing what each holds),
    // charged as work to |arm| when a call returns them
    std::vector<LabelQuery::NodeHits> fetch_hits(ArmState &arm, const std::vector<node_index> &keys);
    std::vector<LabelRecorder::NodeLabels> fetch_present(ArmState &arm,
                                                         const std::vector<node_index> &keys);
    // the keys of the next fetch call; |grown| is the call's growth from one key per level
    // under a work budget
    size_t fetch_chunk(size_t grown) const;
    // ---- the chunked deadlines (pass 5, spec §6.8): which deadline a read is paced against —
    // the seed's time budget as the walk reads it (WALK: depth > 0, time_exceeded), as the
    // derivation reads it (DERIVATION: a positive budget only), or none (NONE: depth 0, a
    // validation) — always with the attempt's walk-until when there is one
    enum class Deadline { NONE, WALK, DERIVATION };
    // what a read checks between its chunks; null when reads are not paced (target 0), so that
    // such a read is one piece, as before
    ReadPacing *pacing(Deadline deadline);
    // the stop a paced read saw before a chunk, as a level's trip: the time budget, a cancel,
    // the attempt's walk-until (BudgetTrip), or a gone client (AttemptAborted)
    BudgetTrip paced_trip();
    ReadPacing pacing_;
    // what stopped the last interrupted read: TIME (the seed's time budget) or an external stop
    bool paced_by_time_ = false;
    // the memory bound's soft part (memory_bound_soft, ResourceAccount::soft_overshoot):
    // what is held beyond the admitted account — decoded rows (|scratch|), a cache beyond its
    // allotment, dictionary labels named but not charged yet — observed wherever it is held
    // and before any check after it can throw (GPT review of stage 2, finding 3)
    void observe_soft(uint64_t scratch, uint64_t dictionary = 0);
    // the same for the seed phase, which runs before the account exists
    void observe_seed_scratch(uint64_t scratch);
    // ---- the budget-aware annotation reads (stage 3 of DESIGN-traverse-graphlet.md §14.1;
    // LabelOracle::decode_charged()): what a read may hold — the budget minus the admitted
    // account and what the level's fetch holds beyond it (|level_soft_|); a DecodeBudget of
    // |max| bytes whose charges the test hook can refuse (WalkerHooks::deny_decode)
    uint64_t allowance() const;
    annot::matrix::DecodeBudget decode_budget(DecodeCharge::Where where, Arm arm, uint64_t max);
    // The derivation's window read budget-aware (stage 3): its distinct |rows| in runs within
    // |budget| (a run that does not fit retried in halves), every row admitted against its
    // standalone demand beside the rows read before it, so that where the read stops depends
    // on the rows alone. Every row read: returns rows.size(), |*out| and |*costs| filled and
    // |budget| holding them. Otherwise the position of the first row that does not fit, with
    // |*refusal| saying why; nothing is kept, |budget| is as on entry.
    // |pacing|: its runs at most a paced chunk long, its stop checked before each; a stop
    // returns the position reached with |*refusal| INTERRUPTED and ReadPacing::units the work
    // of the rows read (nothing kept)
    template <class RowT>
    size_t read_window(const std::vector<Row> &rows, annot::matrix::DecodeBudget &budget,
                       std::vector<RowT> *out, std::vector<annot::matrix::RowCost> *costs,
                       FetchRefusal *refusal, ReadPacing *pacing = nullptr);
    // the seed's own bytes, before the account exists (observe_seed_scratch)
    uint64_t seed_bytes() const {
        return seed_upper_.size() * 2 + nodes_.size() * (sizeof(node_index) + 64);
    }
    // the memory a head's admission needs with what the level's fetch holds beside it,
    // observed (memory_bound_soft) whatever the admission decides
    void observe_admission(uint64_t need);

    // what a dictionary label costs the account: its LabelRef and summaries, its name
    // retained twice, and its delivery in the requested detail
    uint64_t label_bytes(LabelId id, const LabelRef &label) const;
    uint64_t dictionary_bytes(size_t *named) const;
    // what the dictionary built so far holds (its labels and names in every copy, the
    // query's or recorder's cache), its delivery apart: held before the depth-0 admission
    uint64_t dictionary_held() const;

    /**
     * The row-diff path cache (LabelOracle::path_cache, the efficiency pass) during this
     * seed. Without a memory budget it is the request's, within the request's bound, kept
     * from seed to seed (only physical work depends on it). Under a memory budget it is the
     * seed's: empty and off until the depth-0 state is admitted (the seed phase and the
     * annotate roots decode as before the efficiency pass, so that what their refusals state
     * is as before too), then within what the label cache leaves of its allotment
     * (enable_path_cache: the allotment is in the account, and the label cache stays as it
     * was), emptied when the seed ends. The scope restores the request's bound when the seed
     * ends.
     */
    struct PathCacheScope {
        PathCacheScope(LabelOracle &oracle, bool memory_budget)
              : oracle(oracle), memory_budget(memory_budget) {
            if (memory_budget) {
                oracle.path_cache().set_bound(0);
                oracle.path_cache().clear();
            }
        }
        ~PathCacheScope() {
            if (memory_budget)
                oracle.path_cache().clear();
            oracle.path_cache().set_bound(oracle.path_cache_max());
        }
        LabelOracle &oracle;
        const bool memory_budget;
    };

    LabelOracle &oracle_;
    const DeBruijnGraph &graph_;
    const Seed &seed_;
    const Strategy &strategy_;
    LabelChangeCost cost_;          // remapped to dictionary ids by validate_seed()
    const std::string &release_id_;
    const WalkerHooks *hooks_;
    // a stop from outside the walk (the server's attempt), null: none
    const AttemptControl *control_;
    // the seed's usage was written to control_->meter (once: at the end of a walk, before its
    // arms move into the result, or on the way out of a failed one)
    mutable bool metered_ = false;
    const size_t k_;
    const Regime regime_;
    const CanonicalDBG *canonical_;
    const NodeFirstCache *nfc_;
    const bool trace_;
    const bool annotate_;
    PathCacheScope path_cache_scope_;

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
    bool seed_stopped_ = false;
    // what the last cap to trip compared against its limit (CapTrigger::demand)
    double cap_demand_ = 0;

    // admissions so far (Admission::ordinal), both arms
    uint64_t admissions_ = 0;

    // ---- budgets (§14). The memory account: base (seed, dictionary, fixed output and,
    // under a budget, the caches' allotments), committed objects, and the reservations
    // of the live heads; its total never exceeds the budget after an admission
    CostModel m_;
    bool budgeted_ = false;          // a §14 budget is set
    uint64_t mem_limit_ = 0;
    uint64_t base_ = 0, committed_ = 0, reserved_ = 0, peak_ = 0, overshoot_ = 0;
    uint64_t fixed_base_ = 0;        // base_ without the dictionary (no label charged)
    uint64_t seed_work_ = 0;         // the seed phase's charged work units
    uint64_t cache_allotment_ = 0;
    // the caches' allotments in the account (the label cache's and the lookahead's), which a
    // memory stop states with its demand (ResourceStop::allotted)
    uint64_t allotted_ = 0;
    size_t dict_charged_ = 0;        // dictionary labels charged (annotate grows it)
    size_t max_lookahead_ = kMaxLookahead;
    uint64_t next_check_ = kWorkCheckInterval;
    // the widest annotation row fetched so far (entries, and coordinates under trace):
    // sizes the fetch calls
    uint64_t widest_row_ = 0;
    // work used at the last comparison with the work budget, and the most charged between
    // two comparisons (ResourceAccount::largest_charge)
    uint64_t compared_at_ = 0;
    uint64_t largest_charge_ = 0;
    uint64_t depth_ = 0;             // the level being processed
    size_t events_written_ = 0;      // events pushed, for the plan's debug check
    // the annotation is read by the budget-aware decode path: a request budget is set and
    // the index has the path (stage 3)
    bool decode_charged_ = false;
    // What the level holds until its heads are processed, beyond the admitted account
    // (bytes): its key and successor lists and what its fetch returned — with the
    // budget-aware reads the held bytes they were charged at (the level's vectors, the rows'
    // hits or label lists), otherwise the estimate stage 2 observed. Subtracted from what a
    // read may hold, and observed at every admission (memory_bound_soft); 0 between levels.
    uint64_t level_soft_ = 0;
    uint64_t decode_charges_ = 0;    // DecodeCharge::ordinal
    // the test hook (WalkerHooks::deny_decode) refused the last charge of a budget-aware read:
    // a refusal whose last charge it denied is injected, one by an admission (a demand or the
    // names that do not fit what is left) never is
    bool decode_denied_ = false;

    // scratch reused across steps
    Scratch scratch_;
    HeadPlan plan_;
    std::vector<Cand> cands_;
    std::vector<const Entry*> sw_;
    std::string step_;
    std::string rc_scratch_;
    std::vector<uint32_t> label_stamp_;
    uint32_t label_epoch_ = 0;
    std::vector<uint32_t> run_stamp_;    // begin_split() / taken_run()
    uint32_t run_epoch_ = 0;
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
        hits = fetch_seed_hits(validation, keys);
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
                left.reserve(live[l].size());
                for (Coord c : live[l]) {
                    left.push_back(c - (keys.size() - 1));
                }
                // moved, not copied: each copy held the label's whole coordinate set once more
                boundary_coords_[static_cast<size_t>(Arm::RIGHT)].push_back(std::move(right));
                boundary_coords_[static_cast<size_t>(Arm::LEFT)].push_back(std::move(left));
            }
        } else {
            DroppedLabel dropped;
            dropped.name = seed_refs[l].name;
            dropped.reason = "seed_unsupported";
            dropped.runs = std::move(runs[l]);
            result_.dropped_labels.push_back(std::move(dropped));
        }
    }
    if (trace_ && mem_limit_) {
        // The validation's peak under trace support: the hits with their coordinates, every
        // label's live coordinate set and both arms' boundary coordinates, all held at once
        // here and none of them charged (the account does not exist yet). Observed before
        // anything can fail the seed, so that a failure and a walk both state them (review
        // of the stage-2 recheck, P2: copies of 9.6 MB under 1 MiB were reported as 3 MiB)
        auto coords_bytes = [](const auto &sets) {
            uint64_t bytes = 0;
            for (const auto &set : sets) {
                bytes += sizeof(set) + set.capacity() * sizeof(Coord);
            }
            return bytes;
        };
        uint64_t held = coords_bytes(live) + next_live.capacity() * sizeof(Coord)
                      + coords_bytes(boundary_coords_[0]) + coords_bytes(boundary_coords_[1]);
        for (const LabelQuery::NodeHits &node : hits) {
            held += sizeof(node) + node.size() * sizeof(LabelQuery::Hit);
            for (const LabelQuery::Hit &h : node) {
                held += h.coords.capacity() * sizeof(Coord);
            }
        }
        observe_seed_scratch(held);
    }
    if (result_.label_dict.empty()) {
        // under `support: trace` a derived label can still be dropped here: the set is
        // derived from k-mer presence and must then also be coordinate-consecutive. The
        // presence carriers exist (the intersection was not empty), so the cause is the
        // trace, and labels the cap cut before this check might have passed it
        if (result_.labels_from_seed) {
            throw SeedDerivationError(SeedDerivationError::NO_TRACE_CARRIER,
                                      "No label derived from the seed supports every "
                                      "k-mer of the seed as one coordinate-consecutive "
                                      "occurrence (support: trace)",
                                      0, static_cast<double>(result_.labels_supporting_total),
                                      "", result_.labels_dropped);
        }
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
    // every extra label must be reachable: entered by some chain of switches from a seed
    // label whose summed cost stays within the loss budget (switch_reach); the walk enforces
    // the cumulative loss switch by switch. One that no chain reaches is refused, all of them
    // named (the first eight, and how many more), as before for one
    const std::vector<double> reach = switch_reach(cost_, result_.label_dict.size(),
                                                   result_.num_seed_labels,
                                                   strategy_.loss_budget);
    std::vector<LabelId> unreachable;
    for (LabelId id = result_.num_seed_labels; id < result_.label_dict.size(); ++id) {
        if (reach[id] == kInfiniteLoss)
            unreachable.push_back(id);
    }
    if (!unreachable.empty()) {
        constexpr size_t kNamed = 8;
        std::string names;
        for (size_t i = 0; i < std::min(kNamed, unreachable.size()); ++i) {
            names += (i ? ", '" : "'") + result_.label_dict[unreachable[i]].name + "'";
        }
        if (unreachable.size() > kNamed)
            names += " and " + std::to_string(unreachable.size() - kNamed) + " more";
        throw std::invalid_argument(
                std::string(unreachable.size() == 1 ? "Extra label " : "Extra labels ") + names
                + (unreachable.size() == 1 ? " is" : " are") + " unreachable: no chain of "
                  "switches from a seed label enters "
                + (unreachable.size() == 1 ? "it" : "them") + " within the loss budget");
    }
}

std::string Walker::echoed(const std::string &name, const std::string &where) const {
    // Under a memory budget a name the INDEX supplies (a header, a column name) is echoed
    // by a failure as a bounded prefix with its length and where it is: a failed seed's
    // result is not admitted, and nothing in the request bounds such a name, so echoing it
    // whole let one failed result exceed the budget (review of the stage-2 fixes, F1: a header
    // of 180,000 control characters, written twice, was 2.16 MB under 1 MiB). Without a
    // budget, or when short, it is echoed whole, as before.
    constexpr size_t kEchoBytes = 256;
    if (!strategy_.max_memory_bytes || name.size() <= kEchoBytes)
        return name;
    size_t cut = kEchoBytes;
    while (cut > 0 && (static_cast<unsigned char>(name[cut]) & 0xC0) == 0x80) {
        --cut;     // never inside a UTF-8 sequence
    }
    return name.substr(0, cut) + "... (" + std::to_string(name.size()) + " bytes; " + where + ")";
}

// the work of one row of the derivation's window as charge_seed charges it: 8 per key, 1 per
// entry and, read with coordinates, 1 per coordinate (its dependency rows apart)
static uint64_t window_row_units(const annot::matrix::BinaryMatrix::SetBitPositions &row) {
    return 8 + row.size();
}

static uint64_t window_row_units(const annot::matrix::MultiIntMatrix::RowTuples &row) {
    uint64_t units = 8;
    for (const auto &entry : row) {
        units += 1 + entry.second.size();
    }
    return units;
}

template <class RowT>
size_t Walker::read_window(const std::vector<Row> &rows, annot::matrix::DecodeBudget &budget,
                           std::vector<RowT> *out, std::vector<annot::matrix::RowCost> *costs,
                           FetchRefusal *refusal, ReadPacing *pacing) {
    using annot::matrix::buffer_bytes;
    using annot::matrix::DecodeStatus;
    using annot::matrix::RowCost;
    const uint64_t at_entry = budget.held();
    const size_t n = rows.size();
    // what the rows admitted so far hold, with the window's vectors
    uint64_t committed = 0;
    auto left = [&]() { return budget.max_bytes() - at_entry - committed; };
    auto refuse = [&](size_t pos, FetchRefusal::Cause cause, uint64_t demand) {
        *refusal = FetchRefusal();
        refusal->cause = cause;
        refusal->position = pos;
        refusal->left = left();
        refusal->held = committed;
        refusal->demand = demand;
        if (cause == FetchRefusal::DECODE && budget.refused_need() > at_entry + committed)
            refusal->need = budget.refused_need() - at_entry - committed;
        std::vector<RowT>().swap(*out);
        std::vector<RowCost>().swap(*costs);
        budget.restore(at_entry);
        return pos;
    };
    out->clear();
    costs->clear();
    const uint64_t vectors = buffer_bytes(n, sizeof(RowT)) + buffer_bytes(n, sizeof(RowCost));
    if (!budget.charge(vectors))
        return refuse(0, FetchRefusal::DECODE, 0);
    committed += vectors;
    out->reserve(n);
    costs->reserve(n);
    DecodePacer &pacer = oracle_.pacer();
    size_t previous = 0;
    double previous_ms = 0;
    size_t run = n;
    for (size_t pos = 0; pos < n; ) {
        size_t len = std::min(run, n - pos);
        if (pacing) {
            if (pacing->stop()) {
                // the rows read so far are decoded work, the caller's to charge
                pacing->interrupted = true;
                pacing->units = 0;
                for (size_t j = 0; j < out->size(); ++j) {
                    pacing->units += window_row_units((*out)[j]) + 8 * (*costs)[j].dependency_rows
                                   + (*costs)[j].dependency_entries;
                }
                return refuse(pos, FetchRefusal::INTERRUPTED, 0);
            }
            len = std::min(len, pacer.next(n - pos, pacing->ms_left(), previous, previous_ms));
        }
        const uint64_t sub_bytes = buffer_bytes(len, sizeof(Row));
        if (!budget.charge(sub_bytes)) {
            if (len == 1)
                return refuse(pos, FetchRefusal::DECODE, 0);
            run = std::max<size_t>(1, len / 2);
            continue;
        }
        const std::vector<Row> sub(rows.begin() + pos, rows.begin() + pos + len);
        std::vector<RowT> got;
        std::vector<RowCost> got_costs;
        std::vector<uint64_t> held;
        DecodeStatus status;
        Timer timer;
        if constexpr(std::is_same_v<RowT, annot::matrix::BinaryMatrix::SetBitPositions>) {
            status = oracle_.get_rows(sub, budget, &got, &got_costs, &held);
        } else {
            status = oracle_.get_row_tuples(sub, budget, &got, &got_costs, &held);
        }
        previous_ms = timer.elapsed() * 1000;
        pacer.record(len, previous_ms);
        previous = len;
        if (status != DecodeStatus::OK) {
            assert(status == DecodeStatus::REFUSED);
            if (len == 1)
                return refuse(pos, FetchRefusal::DECODE, 0);
            budget.release(sub_bytes);
            run = std::max<size_t>(1, len / 2);
            continue;
        }
        for (size_t j = 0; j < len; ++j) {
            if (got_costs[j].demand > left())
                return refuse(pos + j, FetchRefusal::DEMAND, got_costs[j].demand);
            out->push_back(std::move(got[j]));
            costs->push_back(got_costs[j]);
            committed += held[j];
        }
        // the call's three vectors and the sub-batch are freed with this scope
        budget.release(annot::matrix::IRowDiff::output_bytes<RowT>(len) + sub_bytes);
        pos += len;
    }
    assert(budget.held() == at_entry + committed);
    return n;
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
    // With the budget-aware reads (stage 3) a sub-batch is a read that may be refused
    // whole, so its width is fixed whatever batch_kmers says: whether the seed phase fits
    // must not depend on a fetch-size knob (the derived set never does)
    const size_t chunk = decode_charged_ ? kMaxChunk
                                         : std::clamp<size_t>(strategy_.batch_kmers, 1, kMaxChunk);
    // the derivation's own state, as the soft observation below counts it
    auto state_bytes = [&]() {
        uint64_t bytes = (live.capacity() + next.capacity()) * sizeof(Key)
                       + cur.capacity() * sizeof(cur[0]) + flat.capacity() * sizeof(flat[0]);
        for (const auto &at : coords_at) {
            bytes += sizeof(at) + at.size() * sizeof(at[0]);
            for (const auto &entry : at) {
                bytes += entry.second.size() * sizeof(Coord);
            }
        }
        return bytes;
    };
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
        // the work of each distinct row's row-diff dependency rows (budget-aware reads only)
        std::vector<uint64_t> dependency;
        // The window's read is paced under a deadline (pass 5): the seed's time budget as the
        // derivation reads it after every k-mer, and the attempt's walk-until. A stop between
        // its chunks charges what they decoded (work done, stated without a comparison, as a
        // window found too_wide) and ends the derivation as the next k-mer's check would have:
        // time_budget after the k-mers consumed so far, or the attempt's stop failing the seed
        ReadPacing *pace = pacing(Deadline::DERIVATION);
        auto interrupted = [&](uint64_t units) {
            seed_work_ += units;
            if (paced_by_time_) {
                throw SeedDerivationError(SeedDerivationError::TIME_BUDGET,
                        "The time budget (bounds.time_budget_ms) ran out while deriving the "
                        "permitted set from the seed, after " + std::to_string(done.size())
                        + " of " + std::to_string(keys.size()) + " k-mers; name the labels "
                          "explicitly or shorten the seed",
                        budget_ms, timer_.elapsed() * 1000.0);
            }
            seed_external_stop();
            throw std::logic_error("a paced read of the derivation stopped without a stop");
        };
        if (decode_charged_) {
            // The sub-batch — the window the derivation holds at once, to choose its cheapest
            // row — is read within what the memory budget leaves beside the seed and the
            // derivation's state, every dependency row and tuple charged before it is held,
            // and row by row admitted beside the rows before it: a row that does not fit
            // fails the seed (no result exists yet), and the statement says what it needed
            // and what the window's earlier rows held (review of stage 3, F3)
            const uint64_t beside = seed_bytes() + state_bytes();
            const uint64_t max = !mem_limit_ ? std::numeric_limits<uint64_t>::max()
                               : beside < mem_limit_ ? mem_limit_ - beside : 0;
            annot::matrix::DecodeBudget budget = decode_budget(DecodeCharge::SEED, Arm::RIGHT, max);
            std::vector<annot::matrix::RowCost> costs;
            FetchRefusal r;
            const size_t read = with_coords ? read_window(distinct, budget, &tuples, &costs, &r, pace)
                                            : read_window(distinct, budget, &plain, &costs, &r, pace);
            if (pace && pace->interrupted)
                interrupted(pace->units);
            if (read < distinct.size()) {
                size_t kmer = begin;
                while (rows[kmer - begin] != distinct[read]) {
                    ++kmer;
                }
                const bool injected = decode_denied_ && r.cause == FetchRefusal::DECODE;
                ResourceStop detail;
                detail.cause = ResourceStop::READ_ROW;
                detail.where = ResourceStop::DERIVATION;
                detail.left = r.left;
                detail.held = beside + r.held;
                detail.row_demand = r.cause == FetchRefusal::DEMAND ? r.demand : 0;
                detail.lower_bound = r.cause == FetchRefusal::DECODE;
                detail.index = kmer;
                const uint64_t need = beside + r.held
                    + (detail.lower_bound ? std::max<uint64_t>(r.need, r.left + 1) : r.demand);
                fail_seed(ResourceStop::MEMORY, static_cast<double>(beside + r.held),
                          static_cast<double>(need),
                          (injected ? std::string("an injected refusal of an annotation "
                                                        "read (a test hook, not the budget)")
                                          : "the memory budget (bounds.max_memory_mb = "
                                            + std::to_string(mem_limit_ >> 20) + ")")
                          + " did not admit reading the annotation rows of seed k-mers "
                          + std::to_string(begin) + " .. " + std::to_string(end - 1)
                          + " to derive the permitted set (a window the derivation holds at once, "
                            "to choose its cheapest row), at the row of seed k-mer "
                          + std::to_string(kmer)
                          + refused_read_text(r, injected, "the seed phase had left beside the seed "
                                "and the derivation's state (" + std::to_string(beside)
                                + " bytes) and the window's " + std::to_string(read)
                                + " row(s) read before it (" + std::to_string(r.held) + " bytes)")
                          + "; no traversal was made",
                          "annotation_decode", injected, &detail);
            }
            dependency.resize(distinct.size());
            for (size_t r = 0; r < distinct.size(); ++r) {
                dependency[r] = 8 * costs[r].dependency_rows + costs[r].dependency_entries;
            }
        } else {
            DecodePacer &pacer = oracle_.pacer();
            size_t previous = 0;
            double previous_ms = 0;
            for (size_t at = 0; at < distinct.size(); ) {
                size_t piece = distinct.size() - at;
                if (pace) {
                    if (pace->stop()) {
                        uint64_t units = 0;
                        for (const auto &row : plain) {
                            units += window_row_units(row);
                        }
                        for (const auto &row : tuples) {
                            units += window_row_units(row);
                        }
                        interrupted(units);
                    }
                    piece = pacer.next(piece, pace->ms_left(), previous, previous_ms);
                }
                Timer timer;
                if (!at && piece == distinct.size()) {
                    // one piece: the read as it always was
                    if (with_coords) {
                        tuples = oracle_.get_row_tuples(distinct);
                    } else {
                        plain = oracle_.get_rows(distinct);
                    }
                } else {
                    const std::vector<Row> sub(distinct.begin() + at,
                                               distinct.begin() + at + piece);
                    if (with_coords) {
                        for (auto &row : oracle_.get_row_tuples(sub)) {
                            tuples.push_back(std::move(row));
                        }
                    } else {
                        for (auto &row : oracle_.get_rows(sub)) {
                            plain.push_back(std::move(row));
                        }
                    }
                }
                previous_ms = timer.elapsed() * 1000;
                pacer.record(piece, previous_ms);
                previous = piece;
                at += piece;
            }
        }
        // The sub-batch's decoded rows and the running intersection are held before any
        // account exists, the largest scratch of a derived seed: observed as the soft excess,
        // which a seed failed here states too (finding 8) — after the read, and again once
        // the sub-batch is consumed (the intersection grows while the rows are still held)
        auto observe = [&]() {
            if (!strategy_.max_memory_bytes)
                return;
            uint64_t scratch = state_bytes();
            for (const auto &row : plain) {
                scratch += sizeof(row) + row.size() * sizeof(row[0]);
            }
            for (const auto &row : tuples) {
                scratch += sizeof(row) + row.size() * sizeof(row[0]);
                for (const auto &entry : row) {
                    scratch += entry.second.size() * sizeof(entry.second[0]);
                }
            }
            observe_seed_scratch(scratch);
        };
        observe();
        auto row_of = [&](size_t i) {
            return static_cast<size_t>(
                    std::lower_bound(distinct.begin(), distinct.end(), rows[i - begin])
                    - distinct.begin());
        };

        // the entries of a k-mer's row: its columns, plus its coordinates where they are
        // read (each costs a map_coord); also the row's work units beyond its key's 8
        auto cost_of = [&](size_t i) {
            const size_t r = row_of(i);
            if (!with_coords)
                return plain[r].size();
            size_t n = 0;
            for (const auto &entry : tuples[r]) {
                n += 1 + entry.second.size();
            }
            return n;
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
        // The window's rows, charged as ONE charge as soon as they are read, before a
        // comparison can fail the seed: the window was decoded whole, and charging each k-mer's
        // row when it was consumed let the first row's comparison fail the seed with the
        // window's later rows decoded and never charged (review of the stage-2 recheck, P2: two
        // rows of 400,019 units reported as 200,009). Each k-mer's row is charged as the walk
        // charges a fetched row (8 per key, 1 per entry, and its row-diff dependency rows with
        // the budget-aware reads), a row two k-mers share once per k-mer, as before.
        uint64_t window_units = 0;
        for (size_t i = begin; i < end; ++i) {
            window_units += 8 + cost_of(i) + (dependency.empty() ? 0 : dependency[row_of(i)]);
        }
        if (done.empty()) {
            std::vector<size_t> cost(order.size());
            for (size_t j = 0; j < order.size(); ++j) {
                cost[j] = cost_of(order[j]);
            }
            const size_t best = static_cast<size_t>(
                    std::min_element(cost.begin(), cost.end()) - cost.begin());
            std::iter_swap(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(best));
            if (cost[best] > max_candidates) {
                // the window was read whole before it was found too wide: its rows are the
                // seed's work (the usage a ledger reconciles), added without a comparison, so
                // that too_wide stays the cause stated (review of the stage-3 fixes, P3: a
                // window of 400,019 units reported 0)
                seed_work_ += window_units;
                throw SeedDerivationError(SeedDerivationError::TOO_WIDE,
                        "The permitted set cannot be derived from this seed: the narrowest of "
                        "its first " + std::to_string(end - begin) + " k-mers alone has "
                        + std::to_string(cost[best]) + " annotation entries (columns plus "
                        "k-mer coordinates, an upper bound on its labels), over the limit of "
                        + std::to_string(max_candidates) + ". Start the seed in a less "
                          "repetitive k-mer, or name the labels explicitly.",
                        static_cast<double>(max_candidates), static_cast<double>(cost[best]));
            }
        }
        charge_seed(window_units);

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
                    oracle_.map_coords(c, coords.data(), coords.size(),
                                       [&, column = c](Coord, uint64_t seq_id, Coord local) {
                        flat.emplace_back(Key{ column, seq_id }, local);
                    });
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
                throw SeedDerivationError(SeedDerivationError::NO_CARRIER,
                                          "No label supports every k-mer of the seed (the "
                                          "permitted set derived from the seed is empty at "
                                          "k-mer " + std::to_string(i) + " of "
                                          + std::to_string(keys.size()) + ")",
                                          static_cast<double>(keys.size()),
                                          static_cast<double>(done.size()));
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
            // per k-mer, like the deadline: the window was charged when it was read
            seed_external_stop();
            if (out_of_time()) {
                throw SeedDerivationError(SeedDerivationError::TIME_BUDGET,
                        "The time budget (bounds.time_budget_ms) ran out while deriving the "
                        "permitted set from the seed, after " + std::to_string(done.size())
                        + " of " + std::to_string(keys.size()) + " k-mers; name the labels "
                          "explicitly or shorten the seed",
                        budget_ms, timer_.elapsed() * 1000.0);
            }
        }
        observe();
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
        // where a header is, for a failure that cannot echo it whole (echoed())
        auto where = [](const Key &key) {
            return "column " + std::to_string(key.first) + ", sequence " + std::to_string(key.second);
        };
        std::vector<const Key*> names;
        names.reserve(live.size());
        for (const Key &key : live) {
            names.push_back(&key);
        }
        // stable: of the keys sharing a name, the one an echo locates is deterministic
        std::stable_sort(names.begin(), names.end(), [&](const Key *a, const Key *b) {
            return name_of(*a) < name_of(*b);
        });
        for (size_t i = 1; i < names.size(); ++i) {
            if (name_of(*names[i - 1]) == name_of(*names[i])) {
                const std::string name = echoed(name_of(*names[i]), where(*names[i]));
                throw SeedDerivationError(SeedDerivationError::AMBIGUOUS_HEADER,
                        "The labels derived from the seed are ambiguous: the sequence header '"
                        + name + "' occurs in more than one annotation column, so the "
                        "derived list cannot be resubmitted as explicit labels. Set "
                        "seed_label_kind to \"column\" or name the labels explicitly.",
                        0, 0, name);
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
                          + echoed(oracle_.column_name(back.column),
                                   "column " + std::to_string(back.column)) + ")";
                const std::string shown = echoed(name, where(key));
                throw SeedDerivationError(SeedDerivationError::AMBIGUOUS_HEADER,
                        "The labels derived from the seed are not resubmittable: the sequence "
                        "header '" + shown + "' " + resolves_to + ", which an explicit label list "
                        "would resolve it to. Set seed_label_kind to \"column\" or name the "
                        "labels explicitly.",
                        0, 0, shown);
            }
        }
    }
    // The cap bounds the traversal state, not the discovery: everything above it was
    // found anyway, so it is counted and digested rather than silently forgotten.
    // Under `exhaustive` it is not cut at all: the preset promises that no walk is
    // dropped, and every walk of a dropped carrier would be. Refuse instead, naming
    // the two levers the caller has.
    if (live.size() > strategy_.max_seed_labels && strategy_.exhaustive) {
        throw SeedDerivationError(SeedDerivationError::OVER_SEED_LABEL_CAP,
                std::to_string(live.size()) + " labels carry the seed and max_seed_labels is "
                + std::to_string(strategy_.max_seed_labels) + ": under `exhaustive` the "
                  "derived set is not truncated (every walk of a dropped carrier would be "
                  "missing from the trie); raise labels.max_seed_labels or name the labels",
                static_cast<double>(strategy_.max_seed_labels), static_cast<double>(live.size()));
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

void Walker::charge_seed(uint64_t units) {
    seed_work_ += units;
    // the seed phase has no head to stop at: a stop from outside the walk fails the seed at
    // its next charge (a fetch call's rows, a derivation's window)
    seed_external_stop();
    // Compared once every kWorkCheckInterval units (W). A seed phase that has run past the
    // budget by less than W at a comparison is let go on, so that, once it ends, the first
    // head's check stops the walk with a valid result complete to 0 bp (the result is
    // returned and the overrun stated); one that has run past it by W or more fails the seed,
    // where no result exists yet. Comparing only against the budget failed seeds whose overrun
    // was far below W, which the spec lets finish (review of stage 3, F9). A comparison that
    // lets an overrun go on is not one that passed: what is charged after it stays in the
    // stretch since the last comparison that passed, so that a stop still exceeds the budget
    // by at most the stretch it states (largest_charge).
    if (seed_work_ < next_check_)
        return;
    next_check_ = seed_work_ + kWorkCheckInterval;
    if (!strategy_.max_work_units)
        return;
    if (seed_work_ > strategy_.max_work_units
            && seed_work_ - strategy_.max_work_units < kWorkCheckInterval)
        return;
    // a comparison: what was charged since the last one that passed is stated
    largest_charge_ = std::max(largest_charge_, seed_work_ - compared_at_);
    compared_at_ = seed_work_;
    if (seed_work_ > strategy_.max_work_units) {
        fail_seed(ResourceStop::WORK, static_cast<double>(seed_work_),
                  static_cast<double>(seed_work_),
                  "the work budget (bounds.max_work_units = "
                  + std::to_string(strategy_.max_work_units) + ") ran out while the seed was "
                  + (seed_.labels.empty() ? "read to derive its permitted set"
                                          : "validated against its labels")
                  + ", after " + std::to_string(seed_work_) + " work units: no traversal "
                    "was made");
    }
}

std::vector<LabelQuery::NodeHits> Walker::fetch_seed_hits(LabelQuery &query,
                                                          const std::vector<node_index> &keys) {
    std::vector<LabelQuery::NodeHits> hits;
    uint64_t scratch = 0;
    // the work of the row of key |i| (8 units, 1 per hit and coordinate)
    auto row_work = [&](size_t i) {
        uint64_t units = (keys[i] != npos ? 8 : 0) + hits[i].size();
        for (const LabelQuery::Hit &h : hits[i]) {
            units += h.coords.size();
        }
        return units;
    };
    // Every row a call returned is charged as ONE charge, before the comparison that can fail
    // the seed: the call decoded them all, and charging them row by row let the first row's
    // comparison fail the seed with the call's later rows decoded and never charged (review of
    // the stage-2 recheck, P2: two rows of 400,019 units reported as 200,010). A call is one
    // indivisible charge, as a level's fetch call is (largest_charge states it).
    auto charge = [&](size_t from) {
        // What the call holds is observed BEFORE the charge: a charge can fail the seed, and
        // the failure must state what the fetched rows and the query's cache held (review of
        // the stage-2 fixes, F2: observed 0 while the first row alone held 3.2 MB of
        // coordinates under 1 MiB) — the rule of finding 3, observe then check
        uint64_t units = 0;
        for (size_t i = from; i < hits.size(); ++i) {
            scratch += sizeof(hits[i]) + hits[i].size() * sizeof(LabelQuery::Hit);
            for (const LabelQuery::Hit &h : hits[i]) {
                scratch += h.coords.size() * sizeof(Coord);
            }
            units += row_work(i);
        }
        observe_seed_scratch(scratch + query.cache_bytes());
        charge_seed(units);
    };
    // How many keys the next call reads. A call is decoded whole and charged as one, so under
    // a work budget it is the work a failed seed phase can run past its last comparison by:
    // grown from one key, doubling, it holds about one check interval at the widest row read
    // so far (with its hits, coordinates and, budget-aware, dependency rows), so that a failed
    // seed phase ran past the budget by less than two intervals plus one such call. Sized from
    // the query's labels alone, a call of rows wider than that was one charge of 2.5 million
    // units under a budget of 1 once its rows were charged together. Without a work budget the
    // calls are as before.
    // Under the attempt's walk-until the validation's reads are paced (pass 5): a stop between
    // chunks charges the rows the chunks decoded and fails the seed (no result exists yet). The
    // seed's own time budget never stops it (§6.8: the deadline is not checked at depth 0, and
    // a seed validated past its budget still delivers its result complete to 0 bp)
    ReadPacing *pace = pacing(Deadline::NONE);
    auto interrupted = [&]() {
        if (!pace || !pace->interrupted)
            return;
        charge_seed(pace->units);   // its poll fails the seed: the stop is the attempt's
        seed_external_stop();
        throw std::logic_error("a paced read of the seed phase stopped without a stop");
    };
    uint64_t widest = query.labels().size();
    size_t grown = 1;
    auto next_chunk = [&]() -> size_t {
        if (!strategy_.max_work_units)
            return std::max<size_t>(1, kWorkCheckInterval / (8 + query.labels().size()));
        const size_t chunk = std::clamp<uint64_t>(kWorkCheckInterval / (8 + widest), 1, grown);
        grown = std::min(2 * grown, kFetchChunk);
        return chunk;
    };
    if (decode_charged_) {
        // The budget-aware reads (stage 3), in chunks under any budget: each within what the
        // memory budget leaves beside the seed and the hits read so far, every key admitted
        // against its demand; a chunk that does not fit fails the seed. The validation's
        // query caches nothing (the walk's own query is another), and each row's work includes
        // its row-diff dependency rows.
        query.set_max_cache_bytes(0);
        std::vector<KeyCost> costs;
        hits.reserve(keys.size());
        costs.reserve(keys.size());
        uint64_t held = annot::matrix::buffer_bytes(keys.size(), sizeof(LabelQuery::NodeHits))
                      + annot::matrix::buffer_bytes(keys.size(), sizeof(KeyCost));
        for (size_t begin = 0, end = 0; begin < keys.size(); begin = end) {
            end = std::min(keys.size(), begin + next_chunk());
            const uint64_t beside = seed_bytes() + held;
            const uint64_t max = !mem_limit_ ? std::numeric_limits<uint64_t>::max()
                               : beside < mem_limit_ ? mem_limit_ - beside : 0;
            annot::matrix::DecodeBudget budget = decode_budget(DecodeCharge::SEED, Arm::RIGHT, max);
            size_t refused_at = 0;
            if (pace)
                pace->interrupted = false;
            if (!query.fetch(keys.data() + begin, end - begin, budget, &hits, &costs, &refused_at,
                             pace)) {
                interrupted();
                const FetchRefusal &r = query.refusal();
                const bool injected = decode_denied_ && r.cause == FetchRefusal::DECODE;
                ResourceStop detail;
                detail.cause = ResourceStop::READ_ROW;
                detail.where = ResourceStop::SEED;
                detail.left = r.left;
                detail.held = beside + r.held;
                detail.row_demand = r.cause == FetchRefusal::DEMAND ? r.demand : 0;
                detail.lower_bound = r.cause == FetchRefusal::DECODE;
                detail.index = begin + refused_at;
                const uint64_t need = beside + r.held
                    + (detail.lower_bound ? std::max<uint64_t>(r.need, r.left + 1) : r.demand);
                fail_seed(ResourceStop::MEMORY, static_cast<double>(beside + r.held),
                          static_cast<double>(need),
                          (injected ? std::string("an injected refusal of an annotation "
                                                        "read (a test hook, not the budget)")
                                          : "the memory budget (bounds.max_memory_mb = "
                                            + std::to_string(mem_limit_ >> 20) + ")")
                          + " did not admit reading the annotation row of seed k-mer "
                          + std::to_string(begin + refused_at) + " to validate the seed against "
                            "its labels"
                          + refused_read_text(r, injected, "the seed phase had left beside the seed "
                                "and the hits read before it (" + std::to_string(beside + r.held)
                                + " bytes)")
                          + "; no traversal was made",
                          "annotation_decode", injected, &detail);
            }
            held += budget.held();
            observe_seed_scratch(held);
            // the call's rows, each with its row-diff dependency rows, as one charge (above)
            uint64_t units = 0;
            for (size_t i = begin; i < end; ++i) {
                const uint64_t row = row_work(i) + costs[i].dependency_units;
                widest = std::max<uint64_t>(widest, row - (keys[i] != npos ? 8 : 0));
                units += row;
            }
            charge_seed(units);
        }
        return hits;
    }
    if (!strategy_.max_work_units) {
        // one logical call, as always: the fetch's counters (direct_reads) depend on the
        // batching (paced, its decoding is chunked; its counting and cache stay whole)
        hits = query.fetch(keys, pace);
        interrupted();
        charge(0);
        return hits;
    }
    // Under a work budget in calls (next_chunk), with the comparison between them, so that a
    // long seed cannot spend the budget many times over in one call
    hits.reserve(keys.size());
    for (size_t begin = 0, end = 0; begin < keys.size(); begin = end) {
        end = std::min(keys.size(), begin + next_chunk());
        std::vector<LabelQuery::NodeHits> got
            = query.fetch(std::vector<node_index>(keys.begin() + begin, keys.begin() + end), pace);
        interrupted();
        for (auto &h : got) {
            hits.push_back(std::move(h));
        }
        for (size_t i = begin; i < end; ++i) {
            widest = std::max<uint64_t>(widest, row_work(i) - (keys[i] != npos ? 8 : 0));
        }
        charge(begin);
    }
    return hits;
}

void Walker::fail_seed(ResourceStop::Resource resource, double used, double demand,
                       const std::string &what, const char *phase, bool injected,
                       const ResourceStop *detail) {
    ResourceStop q = detail ? *detail : ResourceStop();
    q.resource = resource;
    q.phase = phase;
    q.injected = injected;
    // no arm's head was refused: the seed itself; RIGHT is the arm a one-armed request
    // walks by default, never read for a failed seed
    q.arm = Arm::RIGHT;
    q.at_bp = 0;
    q.limit = resource == ResourceStop::MEMORY ? static_cast<double>(strategy_.max_memory_bytes)
            : resource == ResourceStop::WORK ? static_cast<double>(strategy_.max_work_units)
            : std::ceil(control_ ? control_->bound_ms : 0);   // the attempt's bound
    q.used = used;
    q.demand = demand;
    if (resource == ResourceStop::MEMORY)
        q.allotted = allotted_;
    ResourceAccount account;
    account.memory_limit = strategy_.max_memory_bytes;
    account.memory_peak = peak_;
    account.memory_final = accounted();
    account.soft_overshoot = overshoot_;
    account.work_limit = strategy_.max_work_units;
    account.work_seed = seed_work_;
    account.work_used = work_used();
    account.largest_charge = largest_charge_;
    account.decode_charged = decode_charged_;
    account.row_diff_uncounted = budgeted_ && !decode_charged_ && oracle_.row_diff();
    const size_t labels = annotate_ ? (recorder_ ? recorder_->labels().size() : 0)
                                    : result_.label_dict.size();
    throw SeedBudgetError(what, q, account, !annotate_ && seed_.labels.empty(), labels);
}

void Walker::fail_depth0(uint64_t need, const Arm *unread, uint64_t labels, uint64_t names) {
    const size_t named = annotate_ ? recorder_->labels().size() : result_.label_dict.size();
    std::string what = "the memory budget (bounds.max_memory_mb = "
        + std::to_string(mem_limit_ >> 20) + ") does not hold the seed's depth-0 state: the seed "
          "and " + std::to_string(named) + " label(s), with both arms' roots and what ending and "
          "delivering them costs in the requested detail, need "
        + (unread ? "at least " : "")
        + std::to_string((need + (uint64_t(1) << 20) - 1) >> 20) + " MiB";
    // The need includes the caches' allotments of this budget, which grow with the budget: the
    // knob value that holds the state is stated, not the need (raised to it, the seed failed
    // again with a larger need; review of the stage-3 fixes, P2)
    const uint64_t knob = memory_budget_holding(need, allotted_) >> 20;
    if (allotted_) {
        what += " with the caches' allotments of this budget (5/16 of it, at most 128 MiB), which "
                "grow with it: " + (unread ? "no budget below " + std::to_string(knob)
                                             + " MiB holds them"
                                           : "the smallest budget that holds them is "
                                             + std::to_string(knob) + " MiB");
    }
    if (unread) {
        what += labels
            ? std::string(" (the ") + to_string(*unread) + " arm's root row names " + std::to_string(labels)
              + " further label(s), " + std::to_string(names) + " bytes with their delivery, "
                "which do not fit beside it and were not named"
            : std::string(" (the ") + to_string(*unread) + " arm's root row was not read: the depth-0 state "
              "before it already reached the budget";
        what += strategy_.direction == Strategy::BOTH && *unread == Arm::LEFT
            ? "; the RIGHT arm's root was not read)" : ")";
        // a lower bound is no promise: what was not built may need more (review of stage 3,
        // answer 2)
        what += "; raising the budget to that may still fail: unread roots, labels or later "
                "state may need more";
    }
    ResourceStop detail;
    detail.lower_bound = unread;
    fail_seed(ResourceStop::MEMORY, static_cast<double>(fixed_base_), static_cast<double>(need),
              what + ": no traversal was made", "traversal", false, &detail);
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
        // moved: the root's entry is the only reader of its arm's boundary coordinates, and
        // keeping a second copy for the whole walk held them twice, uncharged
        if (trace_)
            e.coords = std::move(boundary_coords_[static_cast<size_t>(arm.arm)][l]);
        root.state.push_back(std::move(e));
    }
    if (annotate_) {
        // the boundary k-mer's own labels: the root's entry node, from which
        // continuous presence (label_summary.direct_bp) is measured. Its row is charged
        // here and compared with the budget before the first level reads anything
        // (run_level): both roots belong to the depth-0 result, so a stop between them
        // would deliver none, and the stretch they make with the seed phase is stated as
        // a charge of its own (largest_charge). What the row holds beyond the account is
        // observed now: the root's labels are admitted with the depth-0 state, its row and
        // the cache are not
        const node_index key = oracle_.key_of(root.node, root.kmer);
        std::vector<LabelRecorder::NodeLabels> nl;
        if (decode_charged_) {
            // The budget-aware read (stage 3), within what the account leaves. The depth-0
            // state built so far (the other arm's root, its labels and reservation) may already
            // reach the budget: then the state does not fit whatever this row holds, and the
            // seed fails as a depth-0 state, with its levers, without reading the row (review
            // of stage 3, F7: it was reported as this row's decoding)
            if (mem_limit_ && accounted() >= mem_limit_)
                fail_depth0(accounted() + 1, &arm.arm);
            std::vector<KeyCost> costs;
            nl.reserve(1);
            costs.reserve(1);
            annot::matrix::DecodeBudget budget = decode_budget(DecodeCharge::ROOT, arm.arm,
                                                               allowance());
            size_t refused_at = 0;
            auto name_bytes = [this](std::string_view name) { return recorded_label_bytes(name); };
            if (!recorder_->fetch(&key, 1, budget, &nl, &costs, &refused_at, name_bytes)) {
                const FetchRefusal &r = recorder_->refusal();
                if (r.cause == FetchRefusal::NAMES) {
                    // The row fits; the labels it names, with their delivery in the requested
                    // detail, do not: the depth-0 state does not fit (they are not named). Its
                    // need is at least what admitting the read needed (the row read alone, the
                    // names and their provisional naming), the state after it not known.
                    const uint64_t names = r.names_bytes - r.labels * LabelRecorder::kNamingBytes;
                    fail_depth0(accounted() + r.demand + r.names_bytes, &arm.arm, r.labels, names);
                }
                // the row itself, with its row-diff dependency rows: no result complete to 0 bp
                // exists without it
                const bool injected = decode_denied_ && r.cause == FetchRefusal::DECODE;
                ResourceStop detail;
                detail.cause = ResourceStop::READ_ROW;
                detail.where = ResourceStop::ROOT;
                detail.left = r.left;
                detail.held = accounted();
                detail.row_demand = r.cause == FetchRefusal::DEMAND ? r.demand : 0;
                detail.lower_bound = r.cause == FetchRefusal::DECODE;
                // What the depth-0 state built before it holds beyond the fixed state (the other
                // arm's root with its labels and reservation), and its labels: they compete with
                // the read, so the levers that shrink them help when the row would fit without
                // them (traverse.cpp)
                detail.labels = recorder_->labels().size();
                detail.label_bytes = accounted() > fixed_base_ ? accounted() - fixed_base_ : 0;
                const uint64_t need = accounted()
                    + (detail.lower_bound ? std::max<uint64_t>(r.need, r.left + 1) : r.demand);
                fail_seed(ResourceStop::MEMORY, static_cast<double>(accounted()),
                          static_cast<double>(need),
                          (injected ? std::string("an injected refusal of an annotation "
                                                        "read (a test hook, not the budget)")
                                          : "the memory budget (bounds.max_memory_mb = "
                                            + std::to_string(mem_limit_ >> 20) + ")")
                          + " did not admit reading the annotation row of the "
                          + to_string(arm.arm) + " arm's root (labels.mode annotate)"
                          + refused_read_text(r, injected, "the request had left beside the depth-0 "
                                "state built so far (" + std::to_string(accounted())
                                + " bytes: the seed, the fixed records and the caches' allotments"
                                + (detail.labels ? ", and the other arm's root with its "
                                                   + std::to_string(detail.labels) + " label(s)"
                                                 : std::string())
                                + ")")
                          + "; no traversal was made",
                          "annotation_decode", injected, &detail);
            }
            arm.work_extra += costs[0].dependency_units;
        } else {
            nl = recorder_->fetch({ key });
        }
        widest_row_ = std::max<uint64_t>(widest_row_, nl[0].total);
        arm.work_extra += (key != npos ? 8 : 0) + nl[0].total;
        charge_dictionary();
        observe_soft(sizeof(nl[0]) + nl[0].labels.size() * sizeof(LabelId)
                     + (decode_charged_ ? 0 : recorder_->last_call_bytes()));
        root.present = bounded(arm, nl[0].labels, nl[0].total);
        root.present_total = nl[0].total;
    }
    root.segment = new_segment(arm, {}, 0, labels_at(root),
                               annotate_ ? root.present_total : root.state.size());
    for (Entry &e : root.state) {
        e.run = new_run(arm, e.label, 0, false, 0, 0, UINT32_MAX);
    }
    arm.first_arrival.emplace(root.node, std::make_pair(root.segment, uint64_t(0)));
    // the root segment, its runs and its first arrival, and the root head's
    // reservation; admitted with the rest of the depth-0 state in run() (a budget too
    // small for them fails the seed)
    committed_ += m_.segment + labels_at(root).size() * m_.seg_label
                + root.state.size() * m_.run + m_.step;
    root.reserve = item_bytes(root.state, root.present.size()) + stop_bytes(arm, root);
    reserved_ += root.reserve;
    note_peak();
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
    arm.chain_len.push_back(seg.parents.empty() ? 1 : arm.chain_len[seg.parents[0]] + 1);
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

void Walker::end_run(ArmState &arm, const Entry &e, uint64_t at, EndReason reason,
                     size_t segment) {
    LabelRun &run = arm.result.runs[e.run];
    run.to_bp = at;
    run.ended = true;
    run.end_reason = reason;
    run.segment = static_cast<uint32_t>(segment);
    run.branches = e.branches;
    run.loss = e.loss;
    bin(arm, at).label_ends[static_cast<size_t>(reason)]++;
}

void Walker::end_label(ArmState &arm, const Item &item, const Entry &e, EndReason reason,
                       size_t structural, const char *text, double needed) {
    end_run(arm, e, item.ext_bp, reason, item.segment);
    Event ev;
    ev.at_bp = item.ext_bp;
    ev.type = EventType::LABEL_END;
    ev.label = e.label;
    ev.reason = reason;
    ev.structural_successors = structural;
    ev.needed_budget = needed;
    ev.text = text;
    push_event(arm, item.segment, std::move(ev));
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
    release_head(arm, item);
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
        // the labels recorded on EVERY node whose k-mer lies in the tail (an
        // intersection over the runs covering them; a lower bound where a run's list
        // was cut, which the arm's nodes_labels_truncated reports). The tail's n - k + 1
        // k-mers end at its last n - k + 1 bases: the nodes at outward
        // [ext_bp - (n - k + 1), ext_bp). Its first k - 1 bases are no node's newest base
        // inside the tail -- the nodes there start before it, so their labels must not
        // narrow the list (§7.1: the labels covering that whole tail). A tail of at most k
        // bases (continuation_bp 0 included) lies in the head node's k-mer. Seed nodes
        // carry no recorded labels: a tail reaching more than k - 1 bases into the seed
        // is checked on its flank nodes only. Such a label validates as a seed label of
        // the continuation, so the tail stays valid /traverse input.
        const uint64_t kmers = n > k_ ? n - k_ + 1 : 1;
        const uint64_t node_from = item.ext_bp - std::min<uint64_t>(item.ext_bp, kmers);
        std::vector<LabelId> alive = item.present;
        for (size_t s = item.segment; !alive.empty(); ) {
            const Segment &seg = arm.result.segments[s];
            for (auto it = seg.label_sets.rbegin(); it != seg.label_sets.rend() && !alive.empty(); ++it) {
                if (it->to_bp <= node_from)
                    break;
                std::vector<LabelId> still;
                std::set_intersection(alive.begin(), alive.end(), it->labels.begin(),
                                      it->labels.end(), std::back_inserter(still));
                alive.swap(still);
            }
            if (seg.from_bp <= node_from || seg.parents.empty())
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

bool Walker::reaches_a_target(LabelId from, const std::vector<Target> &targets) const {
    assert(cost_.model() == LabelChangeCost::TABLE);
    if (cost_.default_cost() != kInfiniteLoss) {
        // every pair is finite unless an entry forbids it, so this stops at the first
        // or second target unless the table forbids pair after pair
        for (const Target &t : targets) {
            if (t.label != from && cost_.cost(from, t.label) != kInfiniteLoss)
                return true;
        }
        return false;
    }
    // only the entries of |from| are finite: bounded by the table the caller wrote, not
    // by |targets|, so a cut list costs no more to check than the table is long
    const auto &table = cost_.table();
    for (auto it = table.lower_bound({ from, 0 });
            it != table.end() && it->first.first == from; ++it) {
        const LabelId to = it->first.second;
        if (to == from || it->second == kInfiniteLoss)
            continue;
        auto t = std::lower_bound(targets.begin(), targets.end(), to,
                                  [](const Target &x, LabelId l) { return x.label < l; });
        if (t != targets.end() && t->label == to)
            return true;
    }
    return false;
}

void Walker::derive(ArmState &arm, const State &sigma, const std::vector<Target> &targets,
                    const std::vector<uint8_t> &excluded, State *out,
                    bool *truncated, SourceKey *cut, bool *cut_reaches) {
    out->clear();
    *truncated = false;
    *cut_reaches = false;
    const double budget = strategy_.loss_budget;
    const bool loss_only = strategy_.switch_on_loss_only;
    // the derivation's own scan of the sources and targets (its pair evaluations are
    // counted below), charged as work (§14) and compared with the budget before it runs
    arm.work_extra += sigma.size() + targets.size();
    checkpoint(false);

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
                // Whether the cut can matter here: a cut source with a finite switch into
                // a target of this successor might have been the cheapest way in (§6.3).
                // Such a source need not END (it may go on along another successor), so
                // its label end cannot be what reports the cut; this flag is, through
                // ArmResult::switch_sources_cut.
                for (size_t i = strategy_.max_switch_sources; i < sw.size() && !*cut_reaches; ++i) {
                    *cut_reaches = reaches_a_target(sw[i]->label, targets);
                }
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
            arm.result.pair_evaluations += 1;
            checkpoint(false);
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
            arm.result.pair_evaluations += sw.size();
            // one target prices every kept source: the budget is checked within a
            // derivation, not only between heads
            checkpoint(false);
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
                            State &state, uint64_t at, bool split) {
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
            push_event(arm, target_segment, std::move(ev));
        } else if (split && taken_run(e.run)) {
            LabelRun src = arm.result.runs[e.run];
            e.run = new_run(arm, src.label, src.from_bp, src.entered_by_switch,
                            src.from_label, src.switch_cost, src.prev_run);
            // the clone is the same lineage on the other child: where a merge routed it
            // in from a later parent, its displayed evidence starts at that merge too
            // (§7.1: route_bp 0 would claim the whole displayed flank). Every later merge
            // is deeper (merge_level(depth + 1)), so the inherited stamp stays the
            // earliest. |src| is a copy: new_run may reallocate the runs.
            arm.result.runs[e.run].route_bp = src.route_bp;
        }
        e.switched = false;
        e.pred = e.label;
    }
}


/********************************* budgets **********************************/

void Walker::init_budgets() {
    const DeliveryCosts &d = strategy_.delivery;
    // a segment with its walk buffer and leaf record, its visit mark and chain length
    m_.segment = 2 * (sizeof(Segment) + sizeof(std::string) + sizeof(LeafInfo))
               + 2 * 2 * sizeof(uint32_t) + d.segment;
    m_.seg_label = 2 * sizeof(LabelId) + d.segment_label;
    m_.base = 2 + d.base;
    // a step's edge use (used_edges and used_by_segment) and the first arrival at its
    // node, charged whether or not they are new: an upper bound with no probing
    m_.step = 304;
    // a run with its summary root and index (summarize)
    m_.run = 2 * sizeof(LabelRun) + 3 * sizeof(uint32_t) + d.run;
    m_.event = 3 * sizeof(Event) + d.event;
    m_.event_label = 2 * sizeof(LabelId) + d.event_label;
    m_.needed = 2 * sizeof(double);
    m_.leaf = 2 * sizeof(PathResult) + d.leaf;
    // an end label, and the label's ids in labels_end and the continuation
    m_.leaf_label = 2 * sizeof(LabelEnd) + 4 * sizeof(LabelId) + d.leaf_label;
    m_.chain_entry = d.chain_entry;
    m_.cont_base = 1 + d.continuation_base;
    m_.split = 2 * sizeof(Split) + d.split;
    m_.split_branch = 2 * sizeof(SplitBranch) + 2 * 2 * sizeof(size_t) + d.split_branch;
    m_.bevent = 2 * sizeof(BranchEvent) + d.branch_event;
    m_.bevent_entry = 2 * (sizeof(char) + sizeof(size_t)) + d.branch_event_entry;
    m_.refusal = 2 * sizeof(BranchEvent::Refusal) + d.refusal;
    // a merged segment's parent: its id in parents and its labels_via_parent list
    m_.merge_parent = 2 * (sizeof(size_t) + sizeof(std::vector<LabelId>)) + d.merge_parent;
    m_.presence_run = 2 * sizeof(LabelSetRun) + d.presence_run;
    m_.bin = 2 * sizeof(GrowthBin) + d.bin;
    // a live head: the Item (twice: frontier and next), its k-mer when it does not fit
    // the small-string buffer, its label state with trace coordinates or its present list
    m_.item = 2 * sizeof(Item) + (k_ > 22 ? k_ + 1 : 0);
    m_.entry = 2 * sizeof(Entry);
    m_.coord = 2 * sizeof(Coord);
    m_.present = 2 * sizeof(LabelId);
    // a dictionary label: its LabelRef, the per-label scratch and stamps, its summary on
    // both arms and the query's maps; its name twice per byte (the LabelRef and the
    // query's or recorder's copy), its delivery priced by name (label_bytes)
    m_.label = 2 * sizeof(LabelRef) + 2 * 2 * sizeof(LabelArmSummary) + 96 + d.label;
    m_.label_name = 2;

    mem_limit_ = strategy_.max_memory_bytes;
    budgeted_ = strategy_.max_memory_bytes || strategy_.max_work_units;
    // the seed (sequence, nodes, seed node set of both strands), the fixed records of
    // the output and the statement of a stop (resource_stop, Q, K, walk_domain)
    base_ = d.fixed + 4096 + seed_upper_.size() * (2 + d.base)
          + nodes_.size() * (sizeof(node_index) + 64);
    // what the output repeats from the request and the index besides the dictionary: the
    // seed_id and the dropped seed labels with their names and presence runs, each name
    // priced as delivered (escaped) in the requested detail (finding 1)
    base_ += seed_.seed_id.size() + (d.name ? d.name(seed_.seed_id, DeliveryCosts::Name::SEED_ID) : 0);
    for (const DroppedLabel &dl : result_.dropped_labels) {
        base_ += 2 * sizeof(DroppedLabel) + dl.name.size() + dl.reason.size()
               + dl.runs.size() * (2 * sizeof(dl.runs[0]) + d.dropped_run) + d.dropped
               + (d.name ? d.name(dl.name, DeliveryCosts::Name::DROPPED_LABEL) : 0);
    }
    if (mem_limit_) {
        // Fixed allotments of the budget for the caches, which then evict within them:
        // charging what a cache happens to hold would make admission depend on
        // annotation.batch_kmers (the lookahead warms the caches) and on eviction timing.
        cache_allotment_ = std::min<uint64_t>(mem_limit_ / 4, uint64_t(64) << 20);
        const uint64_t lookahead = std::min<uint64_t>(mem_limit_ / 16, uint64_t(64) << 20);
        max_lookahead_ = std::max<size_t>(1, lookahead / 128);
        allotted_ = cache_allotment_ + lookahead;
        assert(allotted_ == memory_allotments(mem_limit_));
        base_ += allotted_;
        if (query_)
            query_->set_max_cache_bytes(cache_allotment_);
        if (recorder_)
            recorder_->set_max_cache_bytes(cache_allotment_);
    }
    fixed_base_ = base_;
    charge_dictionary();
    note_peak();
}

void Walker::enable_path_cache() {
    if (!mem_limit_ || !oracle_.path_cache_max())
        return;
    // The row-diff path cache shares the label cache's allotment: it holds what the label
    // cache leaves of it (re-read at every insert, and trimmed before the label cache grows:
    // LabelOracle::make_room), so that both stay within the allotment the account holds, the
    // label cache evicts as it always did and nothing admitted or observed changes (the
    // lookahead's reads, physical work, are admitted as without the cache: prefetch).
    // Only from here, after both roots were read: a refused annotate root states whether its
    // read alone was refused (a lower bound) or its standalone demand did not fit, and with
    // the cache on, the left root's path could hold the right root's row, whose read then
    // completes where it alone would have been refused — a DEMAND statement where the seed
    // without the cache states a DECODE one (review of the efficiency pass). The levels
    // state only what holds either way (read_trip), so from depth 1 on the cache changes
    // nothing they state
    oracle_.path_cache().set_bound(std::min(oracle_.path_cache_max(), cache_allotment_),
                                   [this]() {
        const uint64_t label = query_ ? query_->cache_bytes()
                             : recorder_ ? recorder_->cache_bytes() : 0;
        return label < cache_allotment_ ? cache_allotment_ - label : 0;
    });
}

uint64_t Walker::label_bytes(LabelId id, const LabelRef &label) const {
    // a seed label's name is written twice by detail full (label_dict and seed.labels),
    // an extra or recorded label's once: the serialiser layer knows, and prices the escaping
    const DeliveryCosts &d = strategy_.delivery;
    const auto use = !annotate_ && id < result_.num_seed_labels ? DeliveryCosts::Name::SEED_LABEL
                                                                : DeliveryCosts::Name::LABEL;
    return m_.label + label.name.size() * m_.label_name + (d.name ? d.name(label.name, use) : 0);
}

void Walker::charge_dictionary() {
    // annotate mode names labels as the walk meets them: charged after the fetch that
    // named them (the overshoot this allows is stated as memory_bound_soft); the seed's
    // dictionary (constrain) and the roots' labels (annotate) are admitted at depth 0
    const std::vector<LabelRef> &dict = annotate_ ? recorder_->labels() : result_.label_dict;
    for (; dict_charged_ < dict.size(); ++dict_charged_) {
        base_ += label_bytes(dict_charged_, dict[dict_charged_]);
    }
    note_peak();
}

uint64_t Walker::item_bytes(const State &state, size_t present) const {
    uint64_t bytes = m_.item + state.size() * m_.entry + present * m_.present;
    if (trace_) {
        for (const Entry &e : state) {
            bytes += e.coords.size() * m_.coord;
        }
    }
    return bytes;
}

// What ending a head at |ext| costs, whatever ends it: a label end and an end label per
// lineage, its labels at the leaf, the leaf with its path chain spelled out (JSON) and a
// continuation of at most min(continuation_bp, |seed| + ext) bases.
uint64_t Walker::stop_bytes(size_t state, size_t present, uint64_t ext, uint32_t chain) const {
    const size_t labels = annotate_ ? present : state;
    const uint64_t continuation = strategy_.continuation_bp
        ? std::min<uint64_t>(strategy_.continuation_bp, seed_upper_.size() + ext) : 0;
    return labels * (m_.seg_label + m_.leaf_label) + state * m_.event + m_.leaf
         + chain * m_.chain_entry + continuation * m_.cont_base;
}

// A head's share of a merge at the end of its level: the merged segment (counted once
// per head, so g heads hold g segments for the one created), its labels_start, the
// partition and the parents' end sets (each at most the head's labels), its parent entry
// (id and partition list), and one label on the reconverge event — or a same_distance
// revisit event under keep.
uint64_t Walker::merge_bytes(size_t labels) const {
    return m_.segment + 3 * labels * m_.seg_label + m_.merge_parent + 32 + m_.event
         + m_.event_label;
}

bool Walker::would_revisit(const ArmState &arm, node_index node, uint64_t ext,
                           bool revisiting) const {
    // arrival() without its side effects
    auto it = arm.first_arrival.find(node);
    return it != arm.first_arrival.end() && it->second.second != ext && !revisiting;
}

void Walker::plan_child(const State &state, size_t present, uint64_t ext, uint32_t chain) {
    plan_.child_reserve.push_back(item_bytes(state, present)
                                  + stop_bytes(state.size(), present, ext, chain));
    plan_.child_merge_reserve.push_back(merge_bytes(annotate_ ? present : state.size()));
}

void Walker::settle(ArmState &arm, uint64_t held) {
    committed_ += plan_.committed;
    reserved_ += plan_.reserved();
    assert(reserved_ >= held);
    reserved_ -= held;
    arm.bins_charged = std::max(arm.bins_charged, plan_.bins_needed);
    note_peak();
}

// The bins of the level at |child_depth| (its own depth and the next, see charge_bins),
// as far as no earlier head charged them: the first head of a level that creates a child
// pays for them in its admission. Charged at the next level's start instead, they could
// put the account over the budget without any admission refusing them, and the stop that
// follows would deliver a result beyond the budget (review of stage 2, finding 3).
void Walker::plan_bins(const ArmState &arm, uint64_t child_depth) {
    const size_t needed = bins_needed(child_depth);
    if (needed > arm.bins_charged) {
        plan_.committed += (needed - arm.bins_charged) * m_.bin;
        plan_.bins_needed = needed;
    }
}

void Walker::release_head(const ArmState &arm, Item &item) {
    const uint64_t stop = stop_bytes(arm, item);
    // the head was admitted with this much reserved for ending it (a merged head with
    // the sum of its parents', which bounds it)
    assert(stop <= item.reserve);
    committed_ += stop;
    assert(reserved_ >= item.reserve + item.merge_reserve);
    reserved_ -= item.reserve + item.merge_reserve;
    item.reserve = item.merge_reserve = 0;
}

size_t Walker::bins_needed(uint64_t depth) const {
    // a level writes into the bins of its own depth and of the next (label ends of its
    // children when they are censored, reconvergences at depth + 1)
    const uint64_t width = std::max<uint64_t>(1, strategy_.profile_bin_bp);
    return (depth + 1) / width + 1;
}

void Walker::charge_bins(ArmState &arm, uint64_t depth) {
    // the roots' level (run()); a later level's bins were charged by the admission of
    // the heads that created it (plan_bins), so this charges nothing there
    const size_t needed = bins_needed(depth);
    if (arm.bins_charged < needed) {
        committed_ += (needed - arm.bins_charged) * m_.bin;
        arm.bins_charged = needed;
        note_peak();
    }
}

uint64_t Walker::work_of(const ArmState &arm) const {
    const ArmResult &r = arm.result;
    return 4 * r.successor_enumerations + r.pair_evaluations + r.edge_reuse_probes
         + r.refusal_scans + r.steps + arm.work_extra;
}

void Walker::checkpoint(bool force) {
    const uint64_t used = work_used();
    // The work budget is compared on every call, and every charge of the walk reaches a
    // checkpoint before the next one (a fetch call's rows, a derivation's scan, each target
    // it prices, each edge-reuse probe, each successor enumeration): a comparison costs
    // nothing, and comparing only once every interval let a level of wide annotation rows
    // run far past the budget (GPT review of stage 2, finding 2: 125,044 units used under a
    // budget of 1). The previous comparison passed, so a stop overruns the budget by what
    // was charged since, at most: the largest such stretch is recorded and stated with its
    // number (review of the stage-2 fixes, F7), rather than a bound that one kind of
    // charge could break.
    if (strategy_.max_work_units) {
        largest_charge_ = std::max(largest_charge_, used - compared_at_);
        compared_at_ = used;
        if (used > strategy_.max_work_units)
            throw BudgetTrip { ResourceStop::WORK, static_cast<double>(used) };
    }
    if (!force && used < next_check_)
        return;
    next_check_ = used + kWorkCheckInterval;
    // The deadline needs the clock, so it is read before every head and at least every
    // interval, so that one wide level cannot overrun it by more than that. Never at depth
    // 0: a zero budget means "no extension", and the boundary check in run() is what ends
    // such a walk after its first level.
    if (depth_ > 0 && time_exceeded())
        throw BudgetTrip { ResourceStop::TIME, timer_.elapsed() * 1000.0 };
    // A stop from outside the walk is read where the deadline is (before every head and at
    // least every interval), at any depth: unlike a zero time budget it says nothing about
    // extension, and a cancelled attempt must stop as soon as it can. It stops the walk like
    // a budget, so the prefix walked so far is delivered (resource_limit)
    switch (external_stop()) {
        case ExternalStop::NONE:
            return;
        case ExternalStop::CANCELLED:
            throw BudgetTrip { ResourceStop::CANCELLED, attempt_ms() };
        case ExternalStop::ATTEMPT_DEADLINE:
        case ExternalStop::CLIENT_GONE:    // thrown as AttemptAborted by external_stop()
            throw BudgetTrip { ResourceStop::ATTEMPT_DEADLINE, attempt_ms() };
    }
}

ReadPacing* Walker::pacing(Deadline deadline) {
    if (!(oracle_.pacer().target_ms > 0))
        return nullptr;
    pacing_ = ReadPacing();
    paced_by_time_ = false;
    pacing_.ms_left = [this, deadline]() {
        double left = std::numeric_limits<double>::infinity();
        const double budget = strategy_.time_budget_ms;
        if (deadline == Deadline::WALK) {
            left = budget > 0 ? budget - timer_.elapsed() * 1000.0 : 0;
        } else if (deadline == Deadline::DERIVATION && budget > 0) {
            left = budget - timer_.elapsed() * 1000.0;
        }
        if (control_ && control_->ms_left)
            left = std::min(left, control_->ms_left());
        return left;
    };
    // the same tests the walk applies after the read (checkpoint, the derivation's
    // out_of_time, the attempt's poll), so that a read stopped here censors the walk exactly
    // where the read run whole would have
    pacing_.stop = [this, deadline]() {
        const double budget = strategy_.time_budget_ms;
        if ((deadline == Deadline::WALK && time_exceeded())
                || (deadline == Deadline::DERIVATION && budget > 0
                        && timer_.elapsed() * 1000.0 >= budget)) {
            paced_by_time_ = true;
            return true;
        }
        return control_ && control_->poll_now && control_->poll_now() != ExternalStop::NONE;
    };
    return &pacing_;
}

BudgetTrip Walker::paced_trip() {
    if (paced_by_time_)
        return BudgetTrip { ResourceStop::TIME, timer_.elapsed() * 1000.0 };
    // the stop flag is set by now: the poll returns it (a gone client throws)
    switch (external_stop()) {
        case ExternalStop::CANCELLED:
            return BudgetTrip { ResourceStop::CANCELLED, attempt_ms() };
        case ExternalStop::NONE:
        case ExternalStop::ATTEMPT_DEADLINE:
        case ExternalStop::CLIENT_GONE:
            break;
    }
    return BudgetTrip { ResourceStop::ATTEMPT_DEADLINE, attempt_ms() };
}

ExternalStop Walker::external_stop() {
    if (!control_ || !control_->poll)
        return ExternalStop::NONE;
    const ExternalStop stop = control_->poll();
    if (stop == ExternalStop::CLIENT_GONE) {
        // nothing of this walk can be delivered: abandoned without finalising anything (the
        // message is logged: it names no request-supplied text, a seed_id can be megabytes)
        throw AttemptAborted("the client closed its connection: the walk was abandoned at "
                             "depth " + std::to_string(depth_) + " bp");
    }
    return stop;
}

void Walker::seed_external_stop() {
    const ExternalStop stop = external_stop();
    if (stop == ExternalStop::NONE)
        return;
    const double ms = attempt_ms();
    const bool cancelled = stop == ExternalStop::CANCELLED;
    fail_seed(cancelled ? ResourceStop::CANCELLED : ResourceStop::ATTEMPT_DEADLINE, ms, ms,
              std::string(cancelled ? "the attempt was cancelled (POST /traverse/cancel)"
                                    : "the attempt reached the time at which the server stops "
                                      "walking it (its duration bound less the larger of half "
                                      "the allowance and the delivery reserve; see "
                                      "usage.bound.walk_until_ms)")
              + " while the seed was "
              + (seed_.labels.empty() ? "read to derive its permitted set"
                                      : "validated against its labels")
              + ", " + std::to_string(static_cast<uint64_t>(ms)) + " ms after the request was "
                "received: no traversal was made");
}

double Walker::attempt_ms() const {
    // whole milliseconds, rounded up: the statements' integers (K, Q and the cap trigger), and
    // never less than what elapsed
    return std::ceil(control_ && control_->elapsed_ms ? control_->elapsed_ms()
                                                      : timer_.elapsed() * 1000.0);
}

void Walker::write_meter() const noexcept {
    if (!control_ || !control_->meter || metered_)
        return;
    metered_ = true;
    AttemptMeter &m = *control_->meter;
    m.work_units = work_used();
    m.work_seed = seed_work_;
    // What the budget admitted: an account past it was refused (the seed failed at depth 0, or
    // the walk stopped) or is held beyond it as the soft excess, so under a memory budget the
    // admitted peak is at most the budget (review of the stage-4 backend, F1: a refused depth-0
    // demand of 21.8 MB was reported as admitted under 1 MiB)
    const uint64_t peak = std::max(peak_, accounted());
    m.memory_peak = mem_limit_ ? std::min(peak, mem_limit_) : peak;
    m.memory_final = accounted();
    m.soft_excess = overshoot_;
    m.walked = true;
}

EndReason Walker::note_stop(const ArmState &arm, const Item *head,
                            ResourceStop::Resource resource, double demand, bool injected,
                            const char *phase, double used, const ResourceStop *detail) {
    if (!result_.resource_stop) {
        ResourceStop q = detail ? *detail : ResourceStop();
        q.resource = resource;
        q.injected = injected;
        q.phase = phase;
        q.arm = arm.arm;
        q.at_bp = head ? head->ext_bp : depth_;
        q.demand = demand;
        switch (resource) {
            case ResourceStop::MEMORY:
                q.limit = static_cast<double>(mem_limit_);
                // a refused read states what the walk held beyond the account too
                q.used = used >= 0 ? used : static_cast<double>(accounted());
                q.allotted = allotted_;
                break;
            case ResourceStop::WORK:
                q.limit = static_cast<double>(strategy_.max_work_units);
                q.used = demand;
                break;
            case ResourceStop::TIME:
                q.limit = strategy_.time_budget_ms;
                q.used = demand;
                break;
            case ResourceStop::CANCELLED:
            case ResourceStop::ATTEMPT_DEADLINE:
                // the attempt's amounts: its bound and its elapsed time, no budget's
                q.limit = std::ceil(control_ ? control_->bound_ms : 0);
                q.used = demand;
                break;
        }
        result_.resource_stop = q;
    }
    // the cap trigger's demand, in the unit of the knob a reader would raise
    // (bounds.max_memory_mb counts whole MiB): the budget that admits it, whose caches'
    // allotments are larger than this budget's (memory_budget_holding; review of the
    // stage-3 fixes, P2: the need at this budget, raised to, failed again)
    switch (resource) {
        case ResourceStop::MEMORY:
            cap_demand_ = static_cast<double>(
                    memory_budget_holding(static_cast<uint64_t>(std::ceil(demand)), allotted_) >> 20);
            return EndReason::RESOURCE_LIMIT;
        case ResourceStop::WORK:
            cap_demand_ = demand;
            return EndReason::RESOURCE_LIMIT;
        case ResourceStop::TIME:
            cap_demand_ = demand;
            return EndReason::TIME_BUDGET;
        case ResourceStop::CANCELLED:
        case ResourceStop::ATTEMPT_DEADLINE:
            // not a time budget that ran out (time_budget would tell a reader to raise one):
            // censored like any stop the request's own budgets did not cause, resource_limit,
            // with the attempt's elapsed milliseconds as the demand
            cap_demand_ = demand;
            return EndReason::RESOURCE_LIMIT;
    }
    return EndReason::RESOURCE_LIMIT;
}

// The bytes the dictionary labels named since |*named| will cost when charge_dictionary()
// charges them (annotate mode names labels inside a fetch): held, and not charged yet.
uint64_t Walker::dictionary_bytes(size_t *named) const {
    uint64_t bytes = 0;
    if (!annotate_)
        return 0;
    const std::vector<LabelRef> &dict = recorder_->labels();
    for (; *named < dict.size(); ++*named) {
        bytes += label_bytes(*named, dict[*named]);
    }
    return bytes;
}

uint64_t Walker::dictionary_held() const {
    const std::vector<LabelRef> &dict = annotate_ ? recorder_->labels() : result_.label_dict;
    uint64_t bytes = annotate_ ? recorder_->cache_bytes() : query_->cache_bytes();
    for (const LabelRef &label : dict) {
        bytes += m_.label - strategy_.delivery.label + label.name.size() * m_.label_name;
    }
    return bytes;
}

void Walker::observe_soft(uint64_t scratch, uint64_t dictionary) {
    if (!mem_limit_)
        return;
    // Only what is held beyond the admitted account is soft: the scratch rows, a cache
    // beyond its allotment and the dictionary labels a fetch named that are not charged
    // yet. The admitted account itself never exceeds the budget (the depth-0 state and
    // every level's bins are admitted), so memory_bound_soft never reports the modelled
    // state's own excess as the decoder's (review of stage 2, finding 3).
    uint64_t cache = 0;
    if (query_) {
        cache = query_->cache_bytes();
    } else if (recorder_) {
        cache = recorder_->cache_bytes();
    }
    const uint64_t soft = scratch + dictionary
                        + (cache > cache_allotment_ ? cache - cache_allotment_ : 0);
    const uint64_t held = accounted() + soft;
    if (held > mem_limit_)
        overshoot_ = std::max(overshoot_, std::min(soft, held - mem_limit_));
}

void Walker::observe_seed_scratch(uint64_t scratch) {
    if (!mem_limit_)
        return;
    // The seed phase runs before the account exists: what it holds is the seed (charged
    // with the depth-0 state, see init_budgets) and the rows it decodes (soft)
    const uint64_t seed = seed_upper_.size() * 2 + nodes_.size() * (sizeof(node_index) + 64);
    const uint64_t held = seed + scratch;
    if (held > mem_limit_)
        overshoot_ = std::max(overshoot_, std::min(scratch, held - mem_limit_));
}

uint64_t Walker::allowance() const {
    if (!mem_limit_)
        return std::numeric_limits<uint64_t>::max();
    const uint64_t held = accounted() + level_soft_;
    return held < mem_limit_ ? mem_limit_ - held : 0;
}

BudgetTrip Walker::read_trip(const FetchRefusal &r) {
    const uint64_t held = accounted() + level_soft_ + r.held;
    // the level's lists and the rows read before the refused one, beside the account
    observe_soft(level_soft_ + r.held);
    BudgetTrip t { ResourceStop::MEMORY, 0, "annotation_decode", static_cast<double>(held),
                   decode_denied_ && r.cause == FetchRefusal::DECODE };
    ResourceStop &d = t.detail;
    d.where = ResourceStop::LEVEL;
    d.left = r.left;
    d.held = held;
    if (r.cause == FetchRefusal::NAMES) {
        // the row fits; the dictionary labels it names first do not: walker state, not decoding
        t.phase = "traversal";
        d.cause = ResourceStop::LABEL_NAMES;
        d.row_demand = r.demand;
        d.labels = r.labels;
        d.label_bytes = r.names_bytes;
        t.demand = static_cast<double>(held + r.demand + r.names_bytes);
    } else {
        d.cause = ResourceStop::READ_ROW;
        // Whether the refused row's exact demand is known depends on whether the lookahead
        // had already read it (DEMAND) or its read was refused (DECODE), and that depends on
        // annotation.batch_kmers. The stop must not (spec §6.8), so a level stop states only
        // what holds either way: the row needed more than the bytes left.
        d.row_demand = 0;
        d.lower_bound = true;
        t.demand = static_cast<double>(held + r.left + 1);
    }
    return t;
}

void Walker::level_lists_trip() {
    if (!mem_limit_ || allowance())
        return;
    observe_soft(level_soft_);
    const uint64_t held = accounted() + level_soft_;
    BudgetTrip t { ResourceStop::MEMORY, static_cast<double>(held + 1), "traversal",
                   static_cast<double>(held), false };
    t.detail.cause = ResourceStop::LEVEL_LISTS;
    t.detail.where = ResourceStop::LEVEL;
    t.detail.held = held;
    t.detail.lower_bound = true;
    throw t;
}

std::string Walker::refused_read_text(const FetchRefusal &r, bool injected,
                                      const std::string &beside) {
    if (injected)
        return "";
    // a demand is an upper bound, so a row refused by it is not said to need more: it is said
    // to be admitted against more
    return std::string(": ")
        + (r.cause == FetchRefusal::DECODE
            ? "the row, read alone with its row-diff dependency rows, needs more"
            : "the row's standalone demand (what reading it alone with its row-diff dependency "
              "rows holds at most, which every read is admitted against) is "
              + std::to_string(r.demand) + " bytes, more")
        + " than the " + std::to_string(r.left) + " bytes " + beside;
}

uint64_t Walker::recorded_label_bytes(std::string_view name) const {
    const DeliveryCosts &d = strategy_.delivery;
    return m_.label + name.size() * m_.label_name
            + (d.name ? d.name(name, DeliveryCosts::Name::LABEL) : 0);
}

annot::matrix::DecodeBudget Walker::decode_budget(DecodeCharge::Where where, Arm arm,
                                                  uint64_t max) {
    annot::matrix::DecodeBudget budget(max);
    decode_denied_ = false;
    if (hooks_ && hooks_->deny_decode) {
        budget.deny = [this, where, arm](uint64_t) {
            const DecodeCharge charge { where, arm, depth_, decode_charges_++ };
            const bool denied = hooks_->deny_decode(charge);
            decode_denied_ = denied;
            return denied;
        };
    }
    return budget;
}

void Walker::observe_admission(uint64_t need) {
    if (!mem_limit_)
        return;
    // The level's fetched rows are held until its heads are processed, beside the admitted
    // account, which grows with every admission: observed here, whatever the admission
    // decides, so that memory_bound_soft states the excess they make with it
    const uint64_t held = need + level_soft_;
    if (held > mem_limit_)
        overshoot_ = std::max(overshoot_, std::min(level_soft_, held - mem_limit_));
}

size_t Walker::fetch_chunk(size_t grown) const {
    // A call's rows are decoded together and charged together, before the comparison that
    // follows the call, so the call is the work a stop can run past the budget by. Sized so
    // that, at the widest row fetched so far, a call holds about one check interval (the
    // deadline is read between calls), and under a work budget no more than the budget has
    // left, so that near the budget a call reads one key; and there grown from one key per
    // level (|grown|), so that rows far wider than any fetched before are met by a small
    // call rather than by a whole chunk of them (review of the stage-2 fixes, F3)
    uint64_t chunk = std::clamp<uint64_t>(kWorkCheckInterval / (8 + widest_row_), 1, kFetchChunk);
    if (strategy_.max_work_units) {
        const uint64_t used = work_used();
        const uint64_t left = strategy_.max_work_units > used ? strategy_.max_work_units - used : 0;
        chunk = std::min<uint64_t>({ chunk, std::max<uint64_t>(1, left / (8 + widest_row_)),
                                     static_cast<uint64_t>(grown) });
    }
    return static_cast<size_t>(chunk);
}

// The work of a fetched row (DESIGN-traverse-graphlet.md §14, "rows decoded" and
// "coordinates mapped"): 8 per key, 1 per entry and 1 per coordinate. A row is charged
// when a fetch returns it — decoded then or by the lookahead ahead of it — whether or not a
// head consumes it: a level cut mid-way decoded its later rows all the same, and charging
// only consumed rows hid that decoding from the budget and from `used` (review of the
// stage-2 fixes, F3: 880,000 of 900,000 decoded entries uncharged)
static uint64_t row_units(node_index key, const LabelQuery::NodeHits &h) {
    uint64_t units = (key != npos ? 8 : 0) + h.size();
    for (const LabelQuery::Hit &hit : h) {
        units += hit.coords.size();
    }
    return units;
}

std::vector<LabelQuery::NodeHits> Walker::fetch_hits(ArmState &arm,
                                                     const std::vector<node_index> &keys) {
    // A call the deadline may fall into is decoded in time-sized chunks with the deadline
    // checked between them (pacing, pass 5): a call's counters, cache and result are those of the
    // whole call, so a level that the deadline does not stop is unchanged; one it stops is
    // censored at its first head, where the whole call would have been (the rows its chunks
    // decoded are charged: decoded work, though no row was returned)
    ReadPacing *pace = pacing(depth_ > 0 ? Deadline::WALK : Deadline::NONE);
    if (!budgeted_) {
        // one logical call, as always: the fetch's counters (direct_reads) depend on the
        // batching (paced, its decoding is chunked; its counting and cache stay whole)
        std::vector<LabelQuery::NodeHits> hits = query_->fetch(keys, pace);
        if (pace && pace->interrupted) {
            arm.work_extra += pace->units;
            throw paced_trip();
        }
        for (size_t i = 0; i < hits.size(); ++i) {
            arm.work_extra += row_units(keys[i], hits[i]);
        }
        return hits;
    }
    if (decode_charged_) {
        // The budget-aware reads (stage 3): each call admits every key it returns against
        // the key's demand, within what the request has left — the budget minus the account
        // and what the level's fetch holds already, so that where a level stops does not
        // depend on how it is cut into calls nor on what the lookahead cached. The level's
        // vectors are held beside the account like its rows; each row's work includes its
        // row-diff dependency rows (KeyCost), whichever read decoded them.
        // A level with no key to read (the radius, where heads only end, or dead ends) reads
        // nothing, so it is not admitted as a read: it finishes within the reservations its
        // heads already hold. Admitting its empty lists stopped a completed radius-0 walk
        // whose depth-0 state filled the budget exactly (review of stage 3, F3).
        if (keys.empty())
            return {};
        std::vector<LabelQuery::NodeHits> hits;
        std::vector<KeyCost> costs;
        hits.reserve(keys.size());
        costs.reserve(keys.size());
        level_soft_ += annot::matrix::buffer_bytes(keys.size(), sizeof(LabelQuery::NodeHits))
                     + annot::matrix::buffer_bytes(keys.size(), sizeof(KeyCost));
        level_lists_trip();
        size_t grown = 1;
        for (size_t begin = 0; begin < keys.size(); ) {
            const size_t end = std::min(keys.size(), begin + fetch_chunk(grown));
            annot::matrix::DecodeBudget budget = decode_budget(DecodeCharge::LEVEL, arm.arm,
                                                               allowance());
            size_t refused_at = 0;
            if (pace)
                pace->interrupted = false;
            if (!query_->fetch(keys.data() + begin, end - begin, budget, &hits, &costs,
                               &refused_at, pace)) {
                if (pace && pace->interrupted) {
                    arm.work_extra += pace->units;
                    throw paced_trip();
                }
                // nothing of the level is committed yet: censored from its first head
                throw read_trip(query_->refusal());
            }
            level_soft_ += budget.held();
            for (size_t i = begin; i < end; ++i) {
                const uint64_t units = row_units(keys[i], hits[i]) + costs[i].dependency_units;
                widest_row_ = std::max<uint64_t>(widest_row_, units - (keys[i] != npos ? 8 : 0));
                arm.work_extra += units;
            }
            observe_soft(level_soft_);
            checkpoint(true);
            grown = std::min(2 * grown, kFetchChunk);
            begin = end;
        }
        return hits;
    }
    // Under a §14 budget in calls (fetch_chunk()), each charged and compared when it
    // returns; the memory a call leaves held beyond the account is observed BEFORE the
    // comparison after it can throw, or a stop inside the fetch would hide it (GPT review
    // of stage 2, finding 3)
    std::vector<LabelQuery::NodeHits> hits;
    hits.reserve(keys.size());
    // the level's key and successor lists, held beside its hits (run_level)
    const uint64_t lists = level_soft_;
    uint64_t scratch = 0;
    size_t grown = 1;
    for (size_t begin = 0; begin < keys.size(); ) {
        const size_t end = std::min(keys.size(), begin + fetch_chunk(grown));
        std::vector<LabelQuery::NodeHits> got
            = query_->fetch(std::vector<node_index>(keys.begin() + begin, keys.begin() + end),
                            pace);
        if (pace && pace->interrupted) {
            arm.work_extra += pace->units;
            observe_soft(lists + scratch + query_->last_call_bytes());
            throw paced_trip();
        }
        for (auto &h : got) {
            const size_t i = hits.size();
            scratch += sizeof(h) + h.size() * sizeof(LabelQuery::Hit);
            const uint64_t units = row_units(keys[i], h);
            for (const auto &hit : h) {
                scratch += hit.coords.size() * sizeof(Coord);
            }
            widest_row_ = std::max<uint64_t>(widest_row_, units - (keys[i] != npos ? 8 : 0));
            arm.work_extra += units;
            hits.push_back(std::move(h));
        }
        // the call's raw rows were held beside the level's hits while it built them
        level_soft_ = lists + scratch;
        observe_soft(level_soft_ + query_->last_call_bytes());
        checkpoint(true);
        grown = std::min(2 * grown, kFetchChunk);
        begin = end;
    }
    return hits;
}

std::vector<LabelRecorder::NodeLabels> Walker::fetch_present(ArmState &arm,
                                                             const std::vector<node_index> &keys) {
    // as fetch_hits; a recorded row costs its whole width (the true count), not the capped
    // list, and the dictionary also grows with what each call names
    auto units = [](node_index key, const LabelRecorder::NodeLabels &nl) {
        return (key != npos ? 8 : 0) + static_cast<uint64_t>(nl.total);
    };
    ReadPacing *pace = pacing(depth_ > 0 ? Deadline::WALK : Deadline::NONE);
    if (!budgeted_) {
        std::vector<LabelRecorder::NodeLabels> present = recorder_->fetch(keys, pace);
        if (pace && pace->interrupted) {
            arm.work_extra += pace->units;
            throw paced_trip();
        }
        for (size_t i = 0; i < present.size(); ++i) {
            arm.work_extra += units(keys[i], present[i]);
        }
        return present;
    }
    if (decode_charged_) {
        // as fetch_hits; the names a call gives are charged inside it (what the account will
        // charge for them), then with the dictionary, so that they are held within the budget;
        // a level with no key to read is not admitted as a read (see fetch_hits)
        if (keys.empty())
            return {};
        std::vector<LabelRecorder::NodeLabels> present;
        std::vector<KeyCost> costs;
        present.reserve(keys.size());
        costs.reserve(keys.size());
        level_soft_ += annot::matrix::buffer_bytes(keys.size(), sizeof(LabelRecorder::NodeLabels))
                     + annot::matrix::buffer_bytes(keys.size(), sizeof(KeyCost));
        level_lists_trip();
        auto name_bytes = [this](std::string_view name) { return recorded_label_bytes(name); };
        size_t grown = 1;
        for (size_t begin = 0; begin < keys.size(); ) {
            const size_t end = std::min(keys.size(), begin + fetch_chunk(grown));
            annot::matrix::DecodeBudget budget = decode_budget(DecodeCharge::LEVEL, arm.arm,
                                                               allowance());
            size_t refused_at = 0;
            if (pace)
                pace->interrupted = false;
            if (!recorder_->fetch(keys.data() + begin, end - begin, budget, &present, &costs,
                                  &refused_at, name_bytes, pace)) {
                if (pace && pace->interrupted) {
                    arm.work_extra += pace->units;
                    throw paced_trip();
                }
                throw read_trip(recorder_->refusal());
            }
            // The lists and the naming charges stay with the level (the names go to the
            // account with the dictionary): the pending labels are freed with the call, but
            // keeping their charge until the level ends makes a later call of the level admit
            // its keys against what one call would have had left
            level_soft_ += budget.held() - recorder_->last_names_bytes();
            charge_dictionary();
            for (size_t i = begin; i < end; ++i) {
                const uint64_t u = units(keys[i], present[i]) + costs[i].dependency_units;
                widest_row_ = std::max<uint64_t>(widest_row_, u - (keys[i] != npos ? 8 : 0));
                arm.work_extra += u;
            }
            observe_soft(level_soft_);
            checkpoint(true);
            grown = std::min(2 * grown, kFetchChunk);
            begin = end;
        }
        return present;
    }
    std::vector<LabelRecorder::NodeLabels> present;
    present.reserve(keys.size());
    const uint64_t lists = level_soft_;
    uint64_t scratch = 0, dictionary = 0;
    size_t named = dict_charged_;
    size_t grown = 1;
    for (size_t begin = 0; begin < keys.size(); ) {
        const size_t end = std::min(keys.size(), begin + fetch_chunk(grown));
        std::vector<LabelRecorder::NodeLabels> got
            = recorder_->fetch(std::vector<node_index>(keys.begin() + begin, keys.begin() + end),
                               pace);
        if (pace && pace->interrupted) {
            arm.work_extra += pace->units;
            observe_soft(lists + scratch + recorder_->last_call_bytes(), dictionary);
            throw paced_trip();
        }
        for (auto &nl : got) {
            const size_t i = present.size();
            scratch += sizeof(nl) + nl.labels.size() * sizeof(LabelId);
            widest_row_ = std::max<uint64_t>(widest_row_, nl.total);
            arm.work_extra += units(keys[i], nl);
            present.push_back(std::move(nl));
        }
        dictionary += dictionary_bytes(&named);
        // the call's raw rows were held beside the level's lists while it built them
        level_soft_ = lists + scratch;
        observe_soft(level_soft_ + recorder_->last_call_bytes(), dictionary);
        checkpoint(true);
        grown = std::min(2 * grown, kFetchChunk);
        begin = end;
    }
    return present;
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
    // Past the seed's deadline no later level runs (run() reads the same clock before the
    // next depth, a head's checkpoint before the other arm's next head), so nothing the
    // lookahead would read can be consumed: its graph steps and reads are skipped. Before
    // the efficiency pass the lookahead ran in full there (2.2 s on a 250 ms budget on
    // refseq33m); only timing and the physical counters can tell the difference
    if (time_exceeded())
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
        // The chain stops at the radius: its n-th node (n from 0, the item's successor at
        // ext_bp + 1 + n) is enumerated and its single successor's row warmed only when that
        // node is below the radius — a head at the radius is ended, not expanded, so its
        // successors are never asked for. Before the efficiency pass the chain ran
        // batch_kmers nodes whatever the radius left (1,047 rows read for 32 consumed on
        // refseq33m); what changes is physical work (timing) and the counters of lookahead
        // work: annotation.direct_reads (a direct-access annotation's lookahead reads single
        // cells), and keys_mapped where the unbudgeted lookahead outgrew kMaxLookahead and
        // its clearing counted the keys
        const uint64_t radius = strategy_.max_extension_bp;
        const uint64_t room = items[i].ext_bp + 1 < radius ? radius - items[i].ext_bp - 1 : 0;
        const size_t chain_max = std::min<uint64_t>(strategy_.batch_kmers, room);
        if (!chain_max)
            continue;
        succ_kmer(arm.arm, items[i].kmer, succs[i][0].ch, &kmer);
        chain_curs.clear();
        chain_entries.clear();
        chain_nodes.clear();
        window.clear();
        for (size_t n = 0; n < chain_max; ++n) {
            if (n > 0 && arm.lookahead.count(cur))
                break;
            // the graph steps of a chain are not free either (about 0.65 ms a first-touch node
            // on refseq33m, 1,000 per chain at batch_kmers 1000): a chain stops at the seed's
            // deadline, as its read does (the entries it made are kept; nothing depends on them)
            if (n > 0 && time_exceeded())
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
        if (budgeted_ && arm.lookahead.size() + chain_curs.size() > max_lookahead_) {
            // Under a budget the lookahead stays within its allotment (charged with the
            // account): it is cleared BEFORE a chain would push it past, and a chain longer
            // than the allotment is not kept at all. A key a cleared or unkept entry carried
            // is mapped again, and counted, when the walk reaches its node — so it is not
            // counted here: keys_mapped then counts each key the walk consumes once, whatever
            // the lookahead did (without a budget a clearing counts it as well, as it always
            // did, which makes the counter depend on annotation.batch_kmers)
            arm.lookahead.clear();
            if (chain_curs.size() > max_lookahead_)
                continue;
        }
        for (size_t j = 0; j < chain_curs.size(); ++j) {
            arm.lookahead.emplace(chain_curs[j], std::move(chain_entries[j]));
        }
    }
    if (!warm_keys.empty()) {
        // Under a deadline the lookahead's read is paced too; a deadline that stops it ends
        // the warming silently, and the next head's checkpoint, which reads the same deadline,
        // stops the walk there, where the whole read would have: nothing depends on the cache.
        // Also at depth 0, where no head checks the deadline: run() reads it before depth 1,
        // so the roots' lookahead read past it was never consumed (the efficiency pass)
        ReadPacing *pace = pacing(Deadline::WALK);
        if (decode_charged_) {
            // within what is left beside the level's rows; a read that does not fit ends the
            // warming silently (budget-aware results never depend on what is cached)
            annot::matrix::DecodeBudget budget = decode_budget(DecodeCharge::WARM, arm.arm,
                                                               allowance());
            // Under a memory budget the lookahead reads until a run does not fit, and with the
            // path cache a run holds less: it reads as if the cache were not there (its rows
            // the cache spared charged; RowDiffCache::admit_as_uncached), so that the cache does
            // not make it read rows no level asks for (review of the efficiency pass)
            struct AdmitAsUncached {
                AdmitAsUncached(LabelOracle &oracle, bool on) : oracle(oracle), on(on) {
                    if (on)
                        oracle.path_cache().set_admit_as_uncached(true);
                }
                ~AdmitAsUncached() {
                    if (on)
                        oracle.path_cache().set_admit_as_uncached(false);
                }
                LabelOracle &oracle;
                const bool on;
            } admit { oracle_, mem_limit_ > 0 };
            if (annotate_) {
                recorder_->warm(warm_keys, budget, pace);
            } else {
                query_->warm(warm_keys, budget, pace);
            }
        } else if (annotate_) {
            recorder_->warm(warm_keys, pace);
        } else {
            query_->warm(warm_keys, pace);
        }
    }
    if (arm.lookahead.size() > max_lookahead_) {
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
        push_event(arm, item.segment, std::move(ev));
    }
    // same distance: resolved in merge_level
}

void Walker::merge_level(ArmState &arm, uint64_t depth) {
    // Every new head holds a merge reservation (§14): consumed below by a merge it takes
    // part in, released otherwise. Nothing here is admitted: the level's heads already
    // were, with this cost included, so a merge cannot fail.
    auto release_merge = [&](Item &item) {
        assert(reserved_ >= item.merge_reserve);
        reserved_ -= item.merge_reserve;
        item.merge_reserve = 0;
    };
    if (arm.next.size() < 2) {
        for (Item &item : arm.next) {
            release_merge(item);
        }
        return;
    }
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
            release_merge(arm.next[g[0]]);
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
                    push_event(arm, item.segment, std::move(ev));
                    committed_ += m_.event + m_.event_label;
                }
                release_merge(item);
                merged.push_back(std::move(item));
            }
            continue;
        }
        // what the merge costs, out of the reservations of the heads it joins (the
        // merged segment, its entry set and partition, the parents' end sets, the ids
        // and the reconverge event); they bound it, as the merged state is at most the
        // union of theirs
        uint64_t merge_cost = m_.segment + m_.event
                            + g.size() * (m_.event_label + m_.merge_parent + 32);
        uint64_t held_merge = 0, held = 0;
        for (size_t j : g) {
            const Item &item = arm.next[j];
            merge_cost += (annotate_ ? item.present.size() : item.state.size()) * m_.seg_label;
            held_merge += item.merge_reserve;
            held += item.reserve;
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
                if (kept && kept->run == e.run) {
                    mseg.labels_via_parent[j].push_back(e.label);
                } else if (e.run != UINT32_MAX) {
                    // closed by this merge on this parent (to_bp and ended set above)
                    LabelRun &closed = arm.result.runs[e.run];
                    closed.segment = static_cast<uint32_t>(item.segment);
                    closed.branches = e.branches;
                    closed.loss = e.loss;
                }
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
        push_event(arm, m, std::move(ev));
        bin(arm, depth).reconvergences++;
        arm.first_arrival[primary.node] = { m, depth };

        primary.segment = m;
        // the merged entry set and its partition (constrain), or the node's own labels
        merge_cost += (annotate_ ? primary.present.size() : 2 * st.size()) * m_.seg_label;
        assert(merge_cost <= held_merge);
        committed_ += merge_cost;
        assert(reserved_ >= held_merge);
        reserved_ -= held_merge;
        // the merged head holds its parents' reservations: they bound its own and its end
        primary.reserve = held;
        primary.merge_reserve = 0;
        primary.state = std::move(st);
        primary.path_id = path_id;
        primary.splits = splits;
        merged.push_back(std::move(primary));
    }
    arm.next.swap(merged);
    note_peak();
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
                                                     arm.next.size(), live_labels, exact,
                                                     static_cast<double>(arm.next.size()) };
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
                                             live, labels, exact, cap_demand_ };
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
    // max_steps, the request's budgets and its deadline bound the whole seed (§6.8 scope
    // table, DESIGN-traverse-graphlet.md §14 "locus" scope): the other arm stops too
    if (reason == EndReason::MAX_STEPS || reason == EndReason::RESOURCE_LIMIT
            || reason == EndReason::TIME_BUDGET) {
        seed_stopped_ = true;
        for (ArmState &other : arms_) {
            if (&other != &arm)
                stop_frontier(other, reason);
        }
    }
}

void Walker::run_level(ArmState &arm, uint64_t depth) {
    depth_ = depth;
    // what the level's fetch holds beside the account, until the level is done
    level_soft_ = 0;
    struct LevelDone {
        uint64_t &soft;
        ~LevelDone() { soft = 0; }
    } level_done { level_soft_ };
    std::vector<Item> items;
    items.swap(arm.frontier);
    sort_items(items);
    record_live(arm, depth, items);

    // structural successors of every item, then one batched label fetch (in chunks with
    // a check between them under a §14 budget). A budget or the deadline that runs out
    // here stops the level before its first head: nothing of it is committed yet.
    std::vector<std::vector<Succ>> succs(items.size());
    std::vector<node_index> keys;
    std::vector<LabelQuery::NodeHits> hits;
    std::vector<LabelRecorder::NodeLabels> present;
    // A level's heads share its depth. At the radius the level only ends them (their
    // reservations hold that), so nothing can be refused there: stopping it would
    // report an arm whose every walk is present as truncated.
    const bool expands = depth < strategy_.max_extension_bp;
    try {
        // what the walk charged since the last head (the roots' rows at depth 0, a seed
        // phase let finish past the budget) is checked before this level reads anything
        if (expands)
            checkpoint(true);
        charge_bins(arm, depth);
        for (size_t i = 0; i < items.size(); ++i) {
            if (items[i].ext_bp >= strategy_.max_extension_bp)
                continue;
            bool cached = false;
            node_index key = npos;
            succs[i] = successors(arm, items[i].node, items[i].kmer, &cached, &key);
            // an enumeration is charged when consumed
            checkpoint(false);
            if (cached && succs[i].size() == 1) {
                keys.push_back(key);
                continue;
            }
            for (const Succ &s : succs[i]) {
                keys.push_back(key_of_succ(arm.arm, items[i].kmer, s));
            }
        }
        if (mem_limit_) {
            // the level's key and successor lists are held until its heads are processed,
            // beside the account like its fetched rows (and observed with them)
            level_soft_ += annot::matrix::buffer_bytes(keys.capacity(), sizeof(node_index))
                         + annot::matrix::buffer_bytes(succs.size(), sizeof(succs[0]));
            for (const auto &s : succs) {
                level_soft_ += annot::matrix::buffer_bytes(s.capacity(), sizeof(Succ));
            }
            // observed now: the fetch may stop the level before it observes anything
            observe_soft(level_soft_);
        }
        if (annotate_) {
            present = fetch_present(arm, keys);
            // the labels the fetch named (observed as soft until here, see fetch_present)
            charge_dictionary();
        } else {
            hits = fetch_hits(arm, keys);
        }
        if (expands) {
            if (mem_limit_ && accounted() > mem_limit_)
                throw BudgetTrip { ResourceStop::MEMORY, static_cast<double>(accounted()) };
            checkpoint(true);
        }
    } catch (const BudgetTrip &t) {
        trip(arm, note_stop(arm, &items[0], t.resource, t.demand, t.injected, t.phase, t.used,
                            &t.detail),
             items, 0);
        return;
    }
    prefetch(arm, items, succs);
    if (mem_limit_) {
        // the lookahead warmed the caches, possibly beyond their allotments (the stage-2
        // reads), while the level's rows are still held: observed before a head's check can
        // throw, with the raw rows the lookahead's read held while it built them
        observe_soft(level_soft_ + (decode_charged_ ? 0
                                    : annotate_ ? recorder_->last_call_bytes()
                                                : query_->last_call_bytes()));
    }

    size_t offset = 0;
    for (size_t i = 0; i < items.size(); ++i) {
        Item &item = items[i];
        if (item.ext_bp >= strategy_.max_extension_bp) {
            censor_item(arm, item, EndReason::MAX_EXTENSION);
            continue;
        }
        const size_t remaining = items.size() - i - 1;
        std::optional<EndReason> cap;
        try {
            cap = annotate_
                ? process_item_annotate(arm, item, succs[i], present.data() + offset, remaining)
                : process_item(arm, item, succs[i], hits.data() + offset, remaining);
        } catch (const BudgetTrip &t) {
            // thrown while the head was planned: censored like a head beyond a cap
            trip(arm, note_stop(arm, &item, t.resource, t.demand, t.injected, t.phase, t.used,
                                &t.detail),
                 items, i);
            return;
        }
        offset += succs[i].size();
        if (cap) {
            trip(arm, *cap, items, i);
            return;
        }
    }
    merge_level(arm, depth + 1);
    beam(arm, depth + 1);
    arm.frontier.swap(arm.next);
    if (hooks_ && hooks_->level)
        hooks_->level(arm.arm, depth, accounted());
    // what the walk will deliver grows with its account: the caller keeps time back for it
    if (control_ && control_->progress)
        control_->progress(accounted());
}

std::optional<EndReason> Walker::process_item(ArmState &arm, Item &item,
                                              const std::vector<Succ> &succs,
                                              const LabelQuery::NodeHits *hits,
                                              size_t remaining_in_level) {
    const Arm side = arm.arm;
    Scratch &sc = scratch_;
    ScratchGuard guard { sc };      // clears the touched labels on every exit
    // PLAN (successors, derived states, label ends, switches, splits) -> ADMIT -> COMMIT
    // (DESIGN-traverse-graphlet.md §14): nothing of the head reaches the result before it
    // is admitted, so that a head refused by a budget leaves no half-written step behind
    HeadPlan &plan = plan_;
    plan.clear();
    // the budgets and the deadline, before anything of the head is planned (its
    // successors' rows were charged when the level's fetch returned them)
    checkpoint(true);

    if (succs.empty()) {
        // a dead end: every lineage ends here, and the path with them — what the head's
        // reservation holds for its end
        plan.committed = stop_bytes(arm, item);
        plan.new_events = item.state.size();
        if (!admit(arm, item, 0))
            return EndReason::RESOURCE_LIMIT;
        const uint64_t held = item.reserve + item.merge_reserve;
        for (const Entry &e : item.state) {
            end_label(arm, item, e, EndReason::DEAD_END, 0, "");
        }
        finish_path(arm, item, std::nullopt);
        settle(arm, held);
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
                    // the coordinates were charged with their row when the fetch returned
                    // it (one indivisible charge with the row, compared after the fetch
                    // call: review of the stage-2 fixes, F5 — charged here as one sum after the
                    // whole row, they ran past the stated bound unchecked)
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
        if (!c.admissible()) {
            derive(arm, item.state, c.targets, sc.excluded, &c.state, &c.truncated, &c.cut,
                   &c.cut_reaches);
        }
    }
    std::vector<LabelId> &ambiguous_over = plan.ambiguous_over;
    std::vector<LabelId> &ambiguous_taken = plan.ambiguous_taken;
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
            derive(arm, item.state, c.targets, sc.excluded, &c.state, &c.truncated, &c.cut,
                   &c.cut_reaches);
            checkpoint(false);
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
            checkpoint(false);
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
            std::vector<LabelId> &labels = refusals_of(plan, c.succ->ch, "branch");
            labels.insert(labels.end(), c.excluded_preds.begin(), c.excluded_preds.end());
        }
        for (size_t j = excluded_from; j < ambiguous_over.size(); ++j) {
            sc.marked[ambiguous_over[j]] = 0;
        }
    }
    assert(!changed);
    // the rounds after the first are re-minimisations; counted by the commit (or the cap
    // below), never by a head the budget refuses
    plan.reminimisations = rounds - 1;
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
    std::vector<Cand*> &followed = plan.followed;
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

    // ---- caps, decided before anything is committed. A capped head still counts its
    // re-minimisation rounds, as it always did: a result without a budget is unchanged
    // by stage 2 (only a head the budget refuses leaves them uncounted)
    if (auto cap = cap_check(arm, nf, remaining_in_level)) {
        count_reminimisations(arm);
        return cap;
    }

    // ---- PLAN-B: what every lineage does at this step, and the refusals
    plan_outcomes(arm, item);
    plan_cost(arm, item);

    // ---- ADMIT
    if (!admit(arm, item, nf))
        return EndReason::RESOURCE_LIMIT;

    // ---- COMMIT (cannot fail)
    const uint64_t held = item.reserve + item.merge_reserve;
#ifndef NDEBUG
    const size_t segments_before = arm.result.segments.size();
    const size_t runs_before = arm.result.runs.size();
    const size_t events_before = events_written_;
#endif
    commit_item(arm, item, succs);
    // the plan foresaw exactly what the commit wrote: the account below is its cost
    assert(arm.result.segments.size() == segments_before + plan.new_segments);
    assert(arm.result.runs.size() == runs_before + plan.new_runs);
    assert(events_written_ == events_before + plan.new_events);
    settle(arm, held);
    return std::nullopt;
}

// COST (§14) of an admitted-to-be constrained step: the objects the commit writes for
// good and the new heads' reservations, from the plan alone (it reads, never writes).
void Walker::plan_cost(const ArmState &arm, const Item &item) {
    HeadPlan &plan = plan_;
    const uint64_t at = item.ext_bp;
    const size_t nf = plan.followed.size();
    uint64_t &cost = plan.committed;
    // blocked / hairpin events
    for (const Cand &c : cands_) {
        if (c.state.empty() || (c.admissible() && !(c.hairpin && c.followed)))
            continue;
        cost += m_.event + c.state.size() * m_.event_label;
        plan.new_events++;
    }
    // the label ends; a head that ends here pays them out of its stop reservation
    for (const SourcePlan &sp : plan.sources) {
        if (sp.kind != SourcePlan::LABEL_END)
            continue;
        cost += (nf ? m_.event : 0) + (sp.reason == EndReason::LOSS_BUDGET ? m_.needed : 0);
        plan.new_events++;
    }
    // the branch event, when it is kept
    if ((!plan.ambiguous_taken.empty() || !plan.refused.empty())
            && arm.result.branch_events.size() < strategy_.max_branch_events) {
        size_t entries = plan.ambiguous_taken.size() + plan.ambiguous_over.size()
                       + plan.dropped.size();
        for (const Cand &c : cands_) {
            entries += c.admissible() && c.initial_labels;
        }
        cost += m_.bevent + entries * m_.bevent_entry;
        for (const BranchEvent::Refusal &r : plan.refused) {
            cost += m_.refusal + r.labels.size() * m_.bevent_entry;
        }
    }
    if (!nf) {
        cost += stop_bytes(arm, item);
        return;
    }
    cost += nf * (m_.step + m_.base);
    plan_bins(arm, at + 1);
    if (nf == 1) {
        const Cand &c = *plan.followed[0];
        for (const Entry &e : c.state) {
            if (e.switched) {
                cost += m_.run + m_.event;
                plan.new_runs++;
                plan.new_events++;
            }
        }
        if (would_revisit(arm, c.succ->node, at + 1, item.revisiting)) {
            cost += m_.event + m_.event_label;
            plan.new_events++;
        }
        plan_child(c.state, 0, at + 1, arm.chain_len[item.segment]);
        return;
    }
    // a split: the split record, the parent's end set, and per child its segment, its
    // branch (labels cut to the cap), its switch-ins and the runs cloned for it
    cost += m_.split + item.state.size() * m_.seg_label;
    begin_split(arm);
    for (const Cand *c : plan.followed) {
        const size_t n = c->state.size();
        cost += m_.split_branch + std::min(n, strategy_.max_labels_per_node) * m_.seg_label
              + m_.segment + n * m_.seg_label;
        plan.new_segments++;
        for (const Entry &e : c->state) {
            if (e.switched) {
                cost += m_.run + m_.event;
                plan.new_runs++;
                plan.new_events++;
            } else if (taken_run(e.run)) {
                cost += m_.run;     // commit_entries clones a run continued twice
                plan.new_runs++;
            }
        }
        if (would_revisit(arm, c->succ->node, at + 1, false)) {
            cost += m_.event + m_.event_label;
            plan.new_events++;
        }
        plan_child(c->state, 0, at + 1, arm.chain_len[item.segment] + 1);
    }
}

// PLAN-B (§14): the per-source decisions of a step that passed the caps — which lineage
// continues, which ends silently (it goes on only under other names) and which ends with
// what reason, and every refusal of the step — in the order the commit replays them, so
// that refusal groups, `dropped` and the label ends keep today's order. Mutates scratch
// and the spent-work counter refusal_scans only.
void Walker::plan_outcomes(ArmState &arm, const Item &item) {
    Scratch &sc = scratch_;
    HeadPlan &plan = plan_;
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
            checkpoint(false);
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
        checkpoint(false);
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
            ? &refusals_of(plan, c.succ->ch, c.quorum_text) : nullptr;
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
    std::vector<LabelId> &dropped = plan.dropped;
    plan.sources.assign(item.state.size(), SourcePlan());
    for (size_t i = 0; i < item.state.size(); ++i) {
        const Entry &src = item.state[i];
        SourcePlan &sp = plan.sources[i];
        LabelId l = src.label;
        if (sc.has_cont[l]) {
            // the lineage continues only under other names (switch events on commit,
            // on the children when this step splits): a run end without an event
            sp.kind = sc.stays[l] ? SourcePlan::CONTINUES : SourcePlan::SILENT_END;
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
                    refusals_of(plan, ch, "loss_budget").push_back(l);
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
        sp.kind = SourcePlan::LABEL_END;
        sp.reason = reason;
        sp.text = text;
        sp.needed = needed;
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
                    refusals_of(plan, c.succ->ch, "loss_budget").push_back(l);
            }
            for (const Entry &e : c.state) {
                sc.marked[e.pred] = 0;
            }
        }
    }
}

// COMMIT (§14) of an admitted constrained step: writes the decisions of PLAN-A/B into the
// result in the order the walker always wrote them (blocked / hairpin events, the label
// ends in state order, the branch event, then the step or the split). Cannot fail.
void Walker::commit_item(ArmState &arm, Item &item, const std::vector<Succ> &succs) {
    const uint64_t at = item.ext_bp;
    HeadPlan &plan = plan_;
    const size_t nf = plan.followed.size();

    // ---- a cut switch-source list that may have changed what this step commits: a loss,
    // an entry of a followed successor, or the labels on a blocked / hairpin event. Once
    // per successor, for the derivation that stands (the last round's); stated as the
    // arm's switch_sources limitation (§7.0)
    for (const Cand &c : cands_) {
        arm.result.switch_sources_cut += c.truncated && c.cut_reaches;
    }
    // ---- the re-minimisations this committed step was decided by (greedy_losses)
    count_reminimisations(arm);

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
        push_event(arm, item.segment, std::move(ev));
    }

    // ---- the label ends decided by the plan, in state order
    for (size_t i = 0; i < item.state.size(); ++i) {
        const Entry &src = item.state[i];
        const SourcePlan &sp = plan.sources[i];
        if (sp.kind == SourcePlan::SILENT_END) {
            end_run(arm, src, at, EndReason::LABEL_LOST, item.segment);
        } else if (sp.kind == SourcePlan::LABEL_END) {
            end_label(arm, item, src, sp.reason, succs.size(), sp.text, sp.needed);
        }
    }

    // ---- branch event: an ambiguity taken or any refusal (a quorum or split-limit
    // stop, an excluded source and a successor refused by the budget each recorded one)
    if (!plan.ambiguous_taken.empty() || !plan.refused.empty()) {
        assert(at >= arm.last_branch_event_bp);
        arm.last_branch_event_bp = at;
        arm.result.branch_events_total++;
        if (arm.result.branch_events.size() >= strategy_.max_branch_events) {
            // the first event not kept: all events below its depth are kept (§7.2)
            arm.result.branch_events_complete_to_bp
                = std::min(arm.result.branch_events_complete_to_bp, at);
        } else {
            BranchEvent be;
            be.at_bp = at;
            be.segment = item.segment;
            for (const Cand &c : cands_) {
                if (!c.admissible() || !c.initial_labels)
                    continue;
                be.chars.push_back(c.succ->ch);
                be.labels_per_successor.push_back(c.initial_labels);
            }
            be.ambiguous = plan.ambiguous_taken;
            be.ambiguous.insert(be.ambiguous.end(), plan.ambiguous_over.begin(),
                                plan.ambiguous_over.end());
            std::sort(be.ambiguous.begin(), be.ambiguous.end());
            be.dropped = plan.dropped;
            be.labels_affected = be.ambiguous.size() + be.dropped.size();
            for (BranchEvent::Refusal &r : plan.refused) {
                std::sort(r.labels.begin(), r.labels.end());
                r.labels.erase(std::unique(r.labels.begin(), r.labels.end()), r.labels.end());
            }
            be.refused = std::move(plan.refused);
            arm.result.branch_events.push_back(std::move(be));
        }
    }

    if (!nf) {
        finish_path(arm, item, std::nullopt);
        return;
    }

    // ---- steps
    account_steps(arm, at, nf);

    if (nf == 1) {
        Cand &c = *plan.followed[0];
        commit_entries(arm, item, item.segment, c.state, at, false);
        arm.walk_seq[item.segment].push_back(c.succ->ch);
        arm.result.segments[item.segment].length_bp++;
        record_edge(arm, item.segment, c);
        item.node = c.succ->node;
        std::swap(item.kmer, c.kmer);       // keep the buffers for the next step
        item.ext_bp = at + 1;
        std::swap(item.state, c.state);
        item.reserve = plan.child_reserve[0];
        item.merge_reserve = plan.child_merge_reserve[0];
        arrival(arm, item);
        arm.next.push_back(std::move(item));
        return;
    }

    // ---- split
    Split split { at, item.segment, {}, !plan.ambiguous_taken.empty(), item.state.size(), {} };
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
    begin_split(arm);
    for (size_t j = 0; j < nf; ++j) {
        Cand &c = *plan.followed[j];
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
        commit_entries(arm, item, child, c.state, at, true);
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
        ni.reserve = plan.child_reserve[j];
        ni.merge_reserve = plan.child_merge_reserve[j];
        arrival(arm, ni);
        arm.next.push_back(std::move(ni));
    }
    arm.result.splits.push_back(std::move(split));
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
        checkpoint(false);
        for (const EdgeUse &u : uses) {
            if (!is_ancestor_or_self(arm, u.segment(), item.segment))
                continue;
            (u.rc() == c.rc ? same : opposite) = true;
        }
    } else {
        arm.result.edge_reuse_probes += arm.ancestors.size() + 1;
        checkpoint(false);
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
                                           size_t remaining_in_level) {
    if (!nf)
        return std::nullopt;
    auto over = [&](uint64_t demand, EndReason reason) {
        cap_demand_ = static_cast<double>(demand);
        return reason;
    };
    if (steps_total_ + nf > strategy_.max_steps)
        return over(steps_total_ + nf, EndReason::MAX_STEPS);
    if (arm.result.output_bp + nf > strategy_.max_output_bp)
        return over(arm.result.output_bp + nf, EndReason::MAX_OUTPUT);
    size_t live_after = remaining_in_level + arm.next.size() + nf;
    if (nf >= 2 && arm.finished_leaves + live_after > strategy_.max_paths)
        return over(arm.finished_leaves + live_after, EndReason::MAX_PATHS);
    if (strategy_.on_overflow == Strategy::STOP && live_after > strategy_.max_live_paths)
        return over(live_after, EndReason::MAX_LIVE_PATHS);
    return std::nullopt;
}

bool Walker::admit(const ArmState &arm, const Item &item, size_t followed) {
    const uint64_t held = item.reserve + item.merge_reserve;
    const uint64_t cost = plan_.committed + plan_.reserved();
    const uint64_t total = accounted();
    const Admission a { arm.arm, admissions_++, item.ext_bp, item.segment, followed,
                        total, held, cost, mem_limit_ };
    // the account after the commit: the head's reservations become its objects and its
    // children's reservations
    const uint64_t need = total - held + cost;
    observe_admission(need);
    const bool over = mem_limit_ && need > mem_limit_;
    if (over || (hooks_ && hooks_->deny && hooks_->deny(a))) {
        note_stop(arm, &item, ResourceStop::MEMORY, static_cast<double>(need), !over);
        return false;
    }
    return true;
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
    HeadPlan &plan = plan_;
    plan.clear();
    checkpoint(true);
    if (succs.empty()) {
        plan.committed = stop_bytes(arm, item);
        if (!admit(arm, item, 0))
            return EndReason::RESOURCE_LIMIT;
        const uint64_t held = item.reserve + item.merge_reserve;
        finish_path(arm, item, EndReason::DEAD_END);
        settle(arm, held);
        return std::nullopt;
    }

    cands_.resize(succs.size());
    std::vector<Cand*> &followed = plan.followed;
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

    // ---- COST: every decision of this mode is structural and made above (PLAN); what
    // the commit below writes is priced from it, so the head is admitted before the
    // first write (§14)
    static const State kNoState;
    uint64_t &cost = plan.committed;
    for (const Cand &c : cands_) {
        if (c.followed && !c.hairpin)
            continue;
        cost += m_.event + c.present.labels.size() * m_.event_label;
        plan.new_events++;
    }
    if (nf)
        plan_bins(arm, at + 1);
    if (!nf) {
        cost += stop_bytes(arm, item);
    } else if (nf == 1) {
        const Cand &c = *followed[0];
        const auto &sets = arm.result.segments[item.segment].label_sets;
        // record_present() extends the segment's last run when the node's set equals it
        const bool extends = !sets.empty() && sets.back().to_bp == at
            && sets.back().labels_total == c.present.total && sets.back().labels == c.present.labels;
        cost += m_.step + m_.base
              + (extends ? 0 : m_.presence_run + c.present.labels.size() * m_.seg_label);
        if (would_revisit(arm, c.succ->node, at + 1, item.revisiting)) {
            cost += m_.event + m_.event_label;
            plan.new_events++;
        }
        plan_child(kNoState, c.present.labels.size(), at + 1, arm.chain_len[item.segment]);
    } else {
        cost += m_.split + item.present.size() * m_.seg_label;
        for (const Cand *c : followed) {
            const size_t n = c->present.labels.size();
            cost += m_.step + m_.base + m_.split_branch + m_.segment + m_.presence_run
                  + 3 * n * m_.seg_label;
            plan.new_segments++;
            if (would_revisit(arm, c->succ->node, at + 1, false)) {
                cost += m_.event + m_.event_label;
                plan.new_events++;
            }
            plan_child(kNoState, n, at + 1, arm.chain_len[item.segment] + 1);
        }
    }
    if (!admit(arm, item, nf))
        return EndReason::RESOURCE_LIMIT;
    const uint64_t held = item.reserve + item.merge_reserve;

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
        push_event(arm, item.segment, std::move(ev));
    }

    if (!nf) {
        // a structural end: the path's reason (no label lineage ends in this mode);
        // a successor that was only a skipped hairpin is a dead end with its event
        finish_path(arm, item, blocked ? block_reason_of_rank(blocked) : EndReason::DEAD_END);
        settle(arm, held);
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
        item.reserve = plan.child_reserve[0];
        item.merge_reserve = plan.child_merge_reserve[0];
        arrival(arm, item);
        arm.next.push_back(std::move(item));
        settle(arm, held);
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
        ni.reserve = plan.child_reserve[j];
        ni.merge_reserve = plan.child_merge_reserve[j];
        arrival(arm, ni);
        arm.next.push_back(std::move(ni));
    }
    arm.result.splits.push_back(std::move(split));
    settle(arm, held);
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
    // The walk's own buffers are MOVED into the result: nothing reads them after this,
    // and a copy would double the retained bases and leaf data at the very end of the
    // request, when the reservation for them is already spent.
    for (size_t s = 0; s < res.segments.size(); ++s) {
        Segment &seg = res.segments[s];
        if (strategy_.sequences) {
            seg.sequence = std::move(arm.walk_seq[s]);
            if (arm.arm == Arm::LEFT)
                std::reverse(seg.sequence.begin(), seg.sequence.end());
        }
        std::stable_sort(seg.events.begin(), seg.events.end(),
                         [](const Event &a, const Event &b) { return a.at_bp < b.at_bp; });
    }
    // One PathResult per leaf, in segment order. The chain is NOT materialised (it is
    // the leaf's first-parent chain, derived when serialised): per leaf this is O(1)
    // plus its own end labels, so a comb-shaped trie finalises in linear time.
    for (size_t s = 0; s < res.segments.size(); ++s) {
        LeafInfo &leaf = arm.leaves[s];
        if (!leaf.is_leaf)
            continue;
        PathResult path;
        path.id = res.paths.size();
        path.leaf = s;
        path.length_bp = res.segments[s].from_bp + res.segments[s].length_bp;
        path.end_reasons = leaf.end_reasons;
        path.path_reason = leaf.path_reason;
        path.end_labels = std::move(leaf.end_labels);
        path.continuation = std::move(leaf.continuation);
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

// Annotate mode names a label when a level's fetch first returns it, before any head of
// the level is admitted. A walk that a budget, a refused admission or the deadline stopped
// after the fetch can then hold labels only the heads it never committed would have
// recorded: in label_dict, its L records and label_summary, but in no segment, split or
// event — a label "met" that appears nowhere in the result (review of stage 2, finding 2).
// Such a result keeps only the labels it records, in their first-seen order, so that the
// kept ids keep their relative order and every sorted list stays sorted. (The labels a
// stop orphans are those first named at or after the first head it censored, i.e. the
// last ids, so the kept ids do not change at all.) A walk that is not stopped records
// every label it named — each fetched successor is followed or stated by an event — so
// this is the identity there. Applied to stops by a request budget or a refused admission
// only: a cap, and a time stop without a budget, keep their dictionary as it was before
// stage 2 (they can hold such a label too; changing them would change results that set
// no budget, which stage 2 leaves byte-identical).
void Walker::compact_dictionary() {
    const size_t n = result_.label_dict.size();
    constexpr LabelId kUnused = std::numeric_limits<LabelId>::max();
    std::vector<LabelId> to(n, kUnused);
    // every list of label ids an annotate result holds (it has no runs, label ends or
    // lineage partitions; REVISIT and RECONVERGE events list segments, not labels)
    auto each_list = [&](auto f) {
        for (ArmState &arm : arms_) {
            ArmResult &r = arm.result;
            assert(r.runs.empty());
            for (Segment &seg : r.segments) {
                f(seg.labels_start);
                f(seg.labels_end);
                for (LabelSetRun &run : seg.label_sets) {
                    f(run.labels);
                }
                for (Event &ev : seg.events) {
                    if (ev.type == EventType::BLOCKED || ev.type == EventType::HAIRPIN)
                        f(ev.labels);
                }
            }
            for (Split &split : r.splits) {
                for (SplitBranch &branch : split.branches) {
                    f(branch.labels);
                }
            }
            for (PathResult &p : r.paths) {
                assert(p.end_labels.empty());
                if (p.continuation)
                    f(p.continuation->labels);
            }
        }
    };
    each_list([&](const std::vector<LabelId> &labels) {
        for (LabelId l : labels) {
            to[l] = 0;
        }
    });
    LabelId kept = 0;
    for (LabelId l = 0; l < n; ++l) {
        if (to[l] != kUnused)
            to[l] = kept++;
    }
    if (kept == n)
        return;
    each_list([&](std::vector<LabelId> &labels) {
        for (LabelId &l : labels) {
            l = to[l];
        }
    });
    std::vector<LabelRef> dict;
    dict.reserve(kept);
    for (LabelId l = 0; l < n; ++l) {
        if (to[l] != kUnused)
            dict.push_back(std::move(result_.label_dict[l]));
    }
    result_.label_dict = std::move(dict);
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
    // A continuation is offered as the seed of the next request (§7.1), and a seed
    // shorter than k is refused: 1 .. k - 1 would hand out continuations that cannot be
    // resubmitted. 0 is "no continuation sequence" and stays allowed.
    if (strategy_.continuation_bp > 0 && strategy_.continuation_bp < k_) {
        throw std::invalid_argument(
                "output.continuation_bp " + std::to_string(strategy_.continuation_bp)
                + " is shorter than k = " + std::to_string(k_) + ": a continuation must be "
                  "valid traverse input; use 0 (no continuation sequence) or at least k");
    }
    // what the walk consumed reaches the caller however it ends: a result, a failed seed or
    // an exception (a ledger reconciles the attempt against it)
    struct Metered {
        const Walker &walker;
        ~Metered() { walker.write_meter(); }
    } metered { *this };
    // the seed phase observes what it holds against the memory budget (observe_seed_scratch)
    mem_limit_ = strategy_.max_memory_bytes;
    budgeted_ = strategy_.max_memory_bytes || strategy_.max_work_units;
    // under a request budget the annotation is read by the budget-aware decode path when the
    // index has it (stage 3 of DESIGN-traverse-graphlet.md §14.1), from the seed phase on
    decode_charged_ = budgeted_ && oracle_.decode_charged();
    try {
        validate_seed();
    } catch (SeedDerivationError &e) {
        // a seed failed in its derivation states memory_bound_soft like any other result
        // under a memory budget, with what the derivation was seen to hold (finding 8)
        e.set_soft_overshoot(overshoot_);
        throw;
    }
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
    init_budgets();

    arms_[static_cast<size_t>(Arm::LEFT)].arm = Arm::LEFT;
    arms_[static_cast<size_t>(Arm::RIGHT)].arm = Arm::RIGHT;
    for (ArmState &arm : arms_) {
        init_arm(arm);
        // the bins a stop before the first level can write into
        if (arm.result.requested)
            charge_bins(arm, 0);
    }
    // the roots hold their coordinates now (charged with the depth-0 state); an arm that is
    // not requested never read its boundary coordinates, which are freed here rather than
    // held, uncharged, for the whole walk
    for (auto &coords : boundary_coords_) {
        std::vector<SmallVector<Coord>>().swap(coords);
    }
    if (annotate_)
        charge_dictionary();
    // What the depth-0 state holds whatever its admission decides — the dictionary's labels
    // with their names in every copy, the query's or recorder's cache — is observed before
    // the admission can fail the seed: a seed the budget does not hold still built it, and
    // its failure states that excess (review of the stage-2 fixes, F6: twelve names of 1 MB
    // under 1 MiB reported memory_bound_soft 0). Within the budget it is no excess.
    observe_seed_scratch(dictionary_held());
    // ---- ADMIT the depth-0 state (§14) like any head: the dictionary and both roots,
    // each reserved with what ending and delivering it costs. A result complete to 0 bp
    // is the shallowest there is, so a budget that does not hold this holds no valid
    // result at all, and the seed fails before anything per label is delivered. Without
    // this the depth-0 result was delivered whatever it cost — ~4 KB per label and arm in
    // detail full, 44 MiB for 2,876 derived labels under a 1 MiB budget, reported as the
    // soft overshoot (review of stage 2, finding 3).
    if (mem_limit_ && accounted() > mem_limit_)
        fail_depth0(accounted());
    enable_path_cache();

    uint64_t depth = 0;
    while (!seed_stopped_) {
        bool any = false;
        for (const ArmState &arm : arms_) {
            any |= !arm.frontier.empty();
        }
        if (!any)
            break;
        if (depth > 0 && time_exceeded()) {
            depth_ = depth;
            // Only a head below the radius is censored by the deadline: stop_frontier ends
            // the others as max_extension_bp, so when every remaining head has reached the
            // radius the walk is complete and nothing was stopped. Stating a time stop then
            // (resource_stop, Q) told the reader to raise a budget that cut nothing (GPT
            // review of stage 2, N1).
            for (const ArmState &arm : arms_) {
                auto head = std::find_if(arm.frontier.begin(), arm.frontier.end(), [&](const Item &item) {
                    return item.ext_bp < strategy_.max_extension_bp;
                });
                if (head != arm.frontier.end()) {
                    note_stop(arm, &*head, ResourceStop::TIME, timer_.elapsed() * 1000.0);
                    break;
                }
            }
            cap_demand_ = timer_.elapsed() * 1000.0;
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
        arm.result.work_units = work_of(arm);
    }
    // a stop by a request budget or a refused admission; a time stop without a budget (and
    // a cap) keeps its dictionary as before stage 2, so that such a result is unchanged
    if (annotate_ && result_.resource_stop
            && (budgeted_ || result_.resource_stop->resource != ResourceStop::TIME))
        compact_dictionary();
    // every head ended within what it held: nothing is reserved any more
    assert(reserved_ == 0);
    ResourceAccount &account = result_.account;
    account.memory_limit = mem_limit_;
    account.memory_peak = peak_;
    account.memory_final = accounted();
    account.soft_overshoot = overshoot_;
    account.work_limit = strategy_.max_work_units;
    account.work_seed = seed_work_;
    account.work_used = work_used();
    account.largest_charge = largest_charge_;
    account.decode_charged = decode_charged_;
    account.row_diff_uncounted = budgeted_ && !decode_charged_ && oracle_.row_diff();
    // before the arms' results move into the seed's
    write_meter();
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


uint64_t memory_allotments(uint64_t budget) {
    return std::min<uint64_t>(budget / 4, uint64_t(64) << 20)
         + std::min<uint64_t>(budget / 16, uint64_t(64) << 20);
}

uint64_t memory_budget_holding(uint64_t need, uint64_t allotted) {
    constexpr uint64_t kMiB = uint64_t(1) << 20;
    // what the account holds whatever the budget: the need without this budget's allotments
    const uint64_t own = need > allotted ? need - allotted : 0;
    if (!allotted || own > std::numeric_limits<uint64_t>::max() / 2)
        return std::max<uint64_t>(1, (own + kMiB - 1) / kMiB) * kMiB;
    // B - memory_allotments(B) never decreases with B (its slope is 11/16, 15/16 or 1), so
    // the budgets that hold |own| beside their allotments are those from the smallest on;
    // one of own + 128 MiB holds it (the allotments never exceed 128 MiB)
    auto holds = [&](uint64_t mib) { return own + memory_allotments(mib * kMiB) <= mib * kMiB; };
    uint64_t lo = 1, hi = (own + (uint64_t(128) << 20) + kMiB - 1) / kMiB;
    while (lo < hi) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (holds(mid)) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return lo * kMiB;
}

// see walker.hpp; validate_seed's rule for labels.extra
std::vector<double> switch_reach(const LabelChangeCost &cost, size_t n, size_t sources,
                                 double budget) {
    std::vector<double> dist(n, kInfiniteLoss);
    for (size_t s = 0; s < std::min(sources, n); ++s) {
        dist[s] = 0;
    }
    if (cost.model() == LabelChangeCost::FORBID)
        return dist;
    if (cost.model() == LabelChangeCost::CONSTANT) {
        const double c = cost.default_cost();
        if (sources && c <= budget) {
            for (size_t l = sources; l < n; ++l) {
                dist[l] = std::min(dist[l], c);
            }
        }
        return dist;
    }
    const auto &table = cost.table();
    const double fallback = cost.default_cost();
    const bool by_default = fallback != kInfiniteLoss && fallback <= budget;
    std::vector<LabelId> pending;
    if (by_default) {
        pending.reserve(n);
        for (size_t l = 0; l < n; ++l) {
            pending.push_back(static_cast<LabelId>(l));
        }
    }
    using Item = std::pair<double, LabelId>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
    for (size_t s = 0; s < std::min(sources, n); ++s) {
        queue.emplace(0.0, static_cast<LabelId>(s));
    }
    std::vector<uint8_t> done(n, 0);
    auto relax = [&](LabelId to, double loss) {
        if (loss <= budget && loss < dist[to]) {
            dist[to] = loss;
            queue.emplace(loss, to);
        }
    };
    while (!queue.empty()) {
        const auto [d, u] = queue.top();
        queue.pop();
        if (done[u] || d > dist[u])
            continue;
        done[u] = 1;
        auto first = table.lower_bound({ u, 0 });
        for (auto it = first; it != table.end() && it->first.first == u; ++it) {
            if (it->first.second < n && it->first.second != u && it->second != kInfiniteLoss)
                relax(it->first.second, d + it->second);
        }
        if (!by_default)
            continue;
        size_t kept = 0;
        for (LabelId v : pending) {
            if (v == u || table.count({ u, v })) {
                // no default from |u| to |v|: a later pop gives it one
                pending[kept++] = v;
            } else {
                relax(v, d + fallback);
            }
        }
        pending.resize(kept);
    }
    for (double &x : dist) {
        if (x > budget)
            x = kInfiniteLoss;
    }
    return dist;
}


void validate_strategy(const Strategy &st, const LabelChangeCost &cost) {
    if (st.loss_budget < 0)
        throw std::invalid_argument("Negative loss budget");
    if (st.min_successor_fraction < 0 || st.min_successor_fraction > 1)
        throw std::invalid_argument("min_successor_fraction must be in [0, 1]");
    if (!st.max_live_paths)
        throw std::invalid_argument("max_live_paths must be positive");
    if (!st.max_labels_per_node)
        throw std::invalid_argument("max_labels_per_node must be positive");
    // Tip and bubble windows (§6.5) are not implemented: accepting a window and walking
    // as if it were 0 would answer a different question than the one asked
    if (st.tip_window_bp) {
        throw std::invalid_argument("branching.tip_window_bp: not implemented (tip windows are not "
                                    "supported yet); set it to 0 or omit it");
    }
    if (st.bubble_window_bp) {
        throw std::invalid_argument("branching.bubble_window_bp: not implemented (bubble windows "
                                    "are not supported yet); set it to 0 or omit it");
    }
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
                         const std::string &release_id,
                         const WalkerHooks *hooks,
                         const AttemptControl *control) {
    validate_strategy(strategy, cost);
    Walker walker(oracle, seed, strategy, cost, release_id, hooks, control);
    return walker.run();
}

std::string spell_path(const ArmResult &arm, const PathResult &path) {
    // walked leaf -> root through first parents: natural orientation is root -> leaf on
    // the right arm (so the pieces are reversed at the end) and leaf -> root on the left
    // (whose segment sequences are stored natural, i.e. reversed walking order)
    std::vector<const std::string*> pieces;
    walk_path_leaf_first(arm, path, [&](size_t s) {
        pieces.push_back(&arm.segments[s].sequence);
        return true;
    });
    if (arm.arm == Arm::RIGHT)
        std::reverse(pieces.begin(), pieces.end());
    std::string out;
    for (const std::string *p : pieces) {
        out += *p;
    }
    return out;
}

std::vector<size_t> path_segments(const ArmResult &arm, const PathResult &path) {
    std::vector<size_t> chain;
    walk_path_leaf_first(arm, path, [&](size_t s) {
        chain.push_back(s);
        return true;
    });
    std::reverse(chain.begin(), chain.end());
    return chain;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
