#ifndef __TRAVERSAL_WALKER_HPP__
#define __TRAVERSAL_WALKER_HPP__

#include <array>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "traversal_types.hpp"
#include "label_oracle.hpp"


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The permitted set could not be DERIVED from a seed (Seed::labels empty): nothing
 * carries every k-mer of it, the derivation ran out of the time budget, the candidate
 * set of the cheapest seed k-mer is too wide to materialise, or the derived names are
 * ambiguous. Unlike a malformed seed or an unknown explicit label, this is a property
 * of the seed and the data, not of the request: the caller has no lever, so a
 * multi-seed request reports it per seed and still traverses the other seeds
 * (spec §6.1). Explicit labels never raise it.
 *
 * The cause is carried as a code with the bound the derivation ran into, so that the
 * per-seed failure can name the request field that would get past it (spec §7.0, a
 * `derivation` limitation) without parsing the message, which is prose and may change.
 */
class SeedDerivationError : public std::invalid_argument {
  public:
    enum Cause {
        // no label carries every k-mer of the seed: |limit| the seed's k-mers,
        // |observed| the k-mers read when the intersection became empty
        NO_CARRIER,
        // `support: trace`: labels carry every k-mer, none as one coordinate-consecutive
        // occurrence: |observed| how many carry every k-mer (|labels_cut| of them were
        // cut by max_seed_labels before the trace check)
        NO_TRACE_CARRIER,
        // the narrowest of the first 64 k-mers has more annotation entries than the
        // derivation materialises: |limit| that bound, |observed| its entries
        TOO_WIDE,
        // bounds.time_budget_ms ran out: |limit| the budget, |observed| the elapsed ms
        TIME_BUDGET,
        // a derived header name is not resubmittable as an explicit label: |subject|
        AMBIGUOUS_HEADER,
        // `exhaustive` refuses to cut the derived set: |limit| max_seed_labels,
        // |observed| the labels carrying the seed
        OVER_SEED_LABEL_CAP,
    };

    SeedDerivationError(Cause cause, const std::string &what, double limit = 0,
                        double observed = 0, std::string subject = "", size_t labels_cut = 0)
          : std::invalid_argument(what), cause_(cause), limit_(limit), observed_(observed),
            subject_(std::move(subject)), labels_cut_(labels_cut) {}

    Cause cause() const { return cause_; }
    double limit() const { return limit_; }
    double observed() const { return observed_; }
    const std::string& subject() const { return subject_; }
    // labels carrying the seed that max_seed_labels cut before the failure (only the
    // trace check runs after the cap): non-zero means raising that knob may succeed
    size_t labels_cut() const { return labels_cut_; }

  private:
    Cause cause_;
    double limit_;
    double observed_;
    std::string subject_;
    size_t labels_cut_;
};
const char* to_string(SeedDerivationError::Cause cause);

/**
 * Cost of switching the supporting label along a path (spec §9). Labels are
 * identified by their index in the *request*: seed labels in the given order,
 * then strategy.extra. A seed label dropped by validation (seed_unsupported)
 * keeps its request index, so a table authored before traverse_seed() runs stays
 * valid; entries naming a dropped label are ignored. SeedResult::label_dict is
 * the compacted dictionary (dropped labels removed) used by all outputs.
 */
class LabelChangeCost {
  public:
    enum Model { FORBID, CONSTANT, TABLE };

    static LabelChangeCost forbid() { return LabelChangeCost(); }
    static LabelChangeCost constant(double value) {
        LabelChangeCost c;
        c.model_ = CONSTANT;
        c.constant_ = value;
        return c;
    }
    // |default_cost| may be kInfiniteLoss
    static LabelChangeCost table(std::map<std::pair<LabelId, LabelId>, double> entries,
                                 double default_cost) {
        LabelChangeCost c;
        c.model_ = TABLE;
        c.table_ = std::move(entries);
        c.constant_ = default_cost;
        return c;
    }

    Model model() const { return model_; }
    bool finite() const { return model_ != FORBID; }
    const std::map<std::pair<LabelId, LabelId>, double>& table() const { return table_; }
    double default_cost() const { return constant_; }
    double cost(LabelId from, LabelId to) const {
        if (from == to)
            return 0;
        switch (model_) {
            case FORBID: return kInfiniteLoss;
            case CONSTANT: return constant_;
            case TABLE: {
                auto it = table_.find({ from, to });
                return it == table_.end() ? constant_ : it->second;
            }
        }
        return kInfiniteLoss;
    }

  private:
    Model model_ = FORBID;
    double constant_ = kInfiniteLoss;
    std::map<std::pair<LabelId, LabelId>, double> table_;
};


/**
 * What labels mean to the walk (spec §6.9).
 *
 * CONSTRAIN: admissibility requires a surviving permitted label (the label state of
 * §6.3, branch limits, quorum, switching). ANNOTATE: admissibility is purely
 * structural — every non-'$' successor is followed, subject only to the per-path
 * edge-reuse rule, hairpin skipping and seed re-entry — and the labels present at
 * every node are RECORDED (Segment::label_sets) instead of filtering anything. There
 * is no permitted set, no loss budget, no switching, no quorum and no branch limit in
 * that mode; it exists to produce an oracle that shares no label logic with the
 * constrained walk, so that one can verify the other.
 */
enum class LabelMode { CONSTRAIN, ANNOTATE };
const char* to_string(LabelMode mode);

struct Strategy {
    enum Direction { BOTH, LEFT, RIGHT };
    enum Order { BREADTH_FIRST, LOWEST_LOSS_FIRST, MOST_SUPPORTED_FIRST };
    enum Overflow { STOP, BEAM };

    // "no limit" for max_label_branches and max_splits_per_path
    static constexpr size_t kUnlimited = std::numeric_limits<size_t>::max();

    Direction direction = BOTH;
    Support support = Support::KMER;

    // The exhaustive preset (spec §6.9): unlimited branch allowance, on_reconverge
    // keep (a trie, not a DAG), no quorum, no beam, no split limit. It sets nothing
    // by itself — validate_strategy() REJECTS a strategy that asks for `exhaustive`
    // and for a conflicting knob, because asking for the exhaustive trie and getting
    // a pruned one is the failure this flag exists to prevent. The JSON layer fills
    // the preset values in for omitted knobs.
    bool exhaustive = false;

    // labels
    LabelMode label_mode = LabelMode::CONSTRAIN;
    // Cap on every per-node label list in the output: the recorded sets of annotate
    // mode (Segment::label_sets, labels_start / labels_end) and the per-branch lists
    // of the trie view in both modes. The true count is always reported beside a cut
    // list (LabelSetRun::labels_total, SplitBranch::labels_distinct) and every cut is
    // counted in ArmResult::nodes_labels_truncated — a silently incomplete label list
    // would defeat the oracle's purpose.
    size_t max_labels_per_node = 64;
    std::vector<std::string> extra;
    double loss_budget = 0;
    bool switch_on_loss_only = true;
    // Pairwise (TABLE) costs only: the cheapest that many sources are considered per
    // step (or kUnlimited). A cut source list can leave a target label unentered, so
    // under `exhaustive` with a TABLE cost it must be kUnlimited (validate_strategy);
    // elsewhere every cut that may have mattered is counted
    // (ArmResult::switch_sources_cut).
    size_t max_switch_sources = 64;
    // Cap on a set DERIVED from the seed (Seed::labels empty). It bounds the
    // traversal state, not the derivation: the labels above the cap are discovered
    // either way and reported as SeedResult::labels_dropped / _digest. Under
    // `exhaustive` a derived set over the cap is REFUSED (SeedDerivationError)
    // instead of cut: the preset promises that no walk is dropped.
    size_t max_seed_labels = 1000;
    // What a label the walker names by itself denotes — a label DERIVED from the seed
    // (Seed::labels empty) or one RECORDED in annotate mode. Unset (the default) means
    // HEADER when the index has a CoordToHeader (labels are the indexed sequences /
    // accessions) and COLUMN otherwise. Ignored when Seed::labels is given.
    std::optional<LabelKind> seed_label_kind;

    // branching
    size_t max_label_branches = 0;     // or kUnlimited
    size_t min_successor_labels = 1;
    double min_successor_fraction = 0;
    // Tip and bubble windows (spec §6.5) are NOT implemented: validate_strategy()
    // refuses a non-zero value rather than accept a knob that would change nothing.
    uint64_t tip_window_bp = 0;
    uint64_t bubble_window_bp = 0;
    bool merge_reconverge = true;
    bool skip_hairpins = true;
    size_t max_splits_per_path = 64;   // or kUnlimited

    // bounds
    uint64_t max_extension_bp = 5000;
    size_t min_live_labels = 1;
    uint64_t max_steps = 1'000'000;
    size_t max_live_paths = 1000;
    size_t max_paths = 10'000;
    uint64_t max_output_bp = 2'000'000;
    double time_budget_ms = 30'000;

    // frontier
    Order order = BREADTH_FIRST;
    Overflow on_overflow = STOP;

    // output
    bool sequences = true;
    uint64_t profile_bin_bp = 100;
    // per arm, the first that many branch events in level order are kept (or
    // kUnlimited); ArmResult::branch_events_complete_to_bp states where a cut starts
    size_t max_branch_events = 100;
    // 0: no continuation sequence (its labels and loss are still reported); otherwise
    // at least k, or the continuation would be shorter than a k-mer and not valid
    // traverse input (traverse_seed() refuses 1 .. k - 1)
    uint64_t continuation_bp = 1000;

    // annotation
    size_t batch_kmers = 64;
};


struct Seed {
    std::string seed_id;
    std::string sequence;
    /**
     * The permitted labels of this seed.
     *
     * MAY BE EMPTY, which means "derive the permitted set from the seed itself" — the
     * design note's `permit: per_hit` default: the permitted set is the resolved hit's
     * own annotation label(s), not all graph labels. The derived set is the labels that
     * support EVERY k-mer of the seed, in ascending `(column, seq_id)` order, of the
     * kind named by Strategy::seed_label_kind, capped at Strategy::max_seed_labels.
     * It costs no extra annotation reads: it is the intersection of the very rows the
     * seed validation reads anyway (spec §6.1). When the set is empty the seed is
     * rejected, exactly as when every given label is dropped.
     *
     * The derived set is the seed's CARRIER set, not one distinguished hit: under the
     * default forbid cost that is `path_common` over exactly the labels that carry the
     * whole seed. To extend under one specific record, name that one label.
     *
     * Deriving replaces the given list everywhere: `label_dict`, the label ids of all
     * outputs and the request indices of a LabelChangeCost TABLE are "seed labels in
     * order, then Strategy::extra" as before, with the derived labels in the place of
     * the given ones (a table authored by name therefore needs an explicit list).
     * An `extra` label that the derivation also produced is rejected as a duplicate.
     *
     * A derived list is resubmittable VERBATIM as an explicit one. Two guards keep that
     * promise, both reported as SeedDerivationError: a HEADER set whose names are not
     * distinct (a FASTA header is unique per column, not per index, so two columns can
     * hold the same accession while the explicit path rejects duplicate names and
     * resolves a name to the first column holding it), and a seed whose cheapest k-mer
     * has more candidates than the derivation will materialise — a seed of exactly k
     * bases has no intersection at all, so its "derived" set would be a whole
     * annotation row.
     */
    std::vector<std::string> labels;
};


enum class EventType : uint8_t {
    LABEL_END,     // label (lineage) leaves the path: reason
    SWITCH,        // lineage continues under another label: label=from, to, cost
    BLOCKED,       // a label-carrying successor was inadmissible: reason, ch, labels
    HAIRPIN,       // skipped self-reverse-complementary step: ch, labels
    TIP,           // alternative not counted as a branch: ch, length_bp
    BUBBLE,        // alternatives reconverging within the window: alleles in |text|
    RECONVERGE,    // another path joined here: labels = segments
    REVISIT,       // a node already visited at another distance: labels = segments
};
const char* to_string(EventType type);

struct Event {
    uint64_t at_bp = 0;
    EventType type = EventType::LABEL_END;
    LabelId label = 0;
    LabelId to = 0;
    double cost = 0;
    double needed_budget = 0;        // LOSS_BUDGET ends
    EndReason reason = EndReason::DEAD_END;
    char ch = 0;
    uint64_t length_bp = 0;
    size_t structural_successors = 0;
    std::vector<LabelId> labels;     // affected labels (or segments for joins)
    // BLOCKED / HAIRPIN: the distinct labels at the successor; > |labels| when the list
    // was cut by Strategy::max_labels_per_node (annotate mode). The list of a not
    // followed successor appears nowhere else, so a cut here is counted in
    // ArmResult::nodes_labels_truncated like any recorded node's.
    size_t labels_total = 0;
    bool truncated() const { return labels_total > labels.size(); }
    std::string text;
};

// Annotate mode: the labels present on every node of a stretch of a segment. Runs are
// maximal stretches with an identical (capped list, total) and are half-open over
// outward base indices, so [from_bp, to_bp) covers the nodes entered by steps
// from_bp + 1 .. to_bp.
struct LabelSetRun {
    uint64_t from_bp = 0;
    uint64_t to_bp = 0;
    std::vector<LabelId> labels;     // ascending, at most Strategy::max_labels_per_node
    size_t labels_total = 0;         // distinct labels at those nodes; > |labels| when cut
    bool truncated() const { return labels_total > labels.size(); }
};

struct Segment {
    size_t id = 0;
    std::vector<size_t> parents;
    // joins only (|parents| > 1): per parent, the labels whose kept entry entered
    // through it, so that per-label routes through the DAG stay reconstructible
    // (constrain mode; empty in annotate mode, which tracks no lineages)
    std::vector<std::vector<LabelId>> labels_via_parent;
    std::vector<size_t> children;
    uint64_t from_bp = 0;
    uint64_t length_bp = 0;
    std::string sequence;            // bases added, natural orientation
    // The labels at the segment's ENTRY node (the seed boundary for the root, the
    // first node of the segment otherwise) and at its last node. In constrain mode
    // these are the live permitted labels; in annotate mode the labels present,
    // capped at max_labels_per_node. |labels_start_total| is the true count at the
    // entry node (> |labels_start| when cut): the root's boundary node is in no run,
    // so this is the only place its count is reported; the last node's count is on
    // the last run.
    std::vector<LabelId> labels_start;
    size_t labels_start_total = 0;
    std::vector<LabelId> labels_end;
    std::vector<Event> events;       // sorted by at_bp
    // annotate mode only: the labels present along the segment, as runs
    std::vector<LabelSetRun> label_sets;
};

struct LabelRun {
    LabelId label = 0;
    uint64_t from_bp = 0;
    uint64_t to_bp = 0;              // half-open; set when the run ends
    // Non-zero when this run's lineage was merged in from another parent at a
    // reconvergence (§6.5) at that extension. A merge joins paths AT A NODE, so the
    // bases from route_bp onward are shared by every route into it, while the earlier
    // [from_bp, route_bp) prefix was travelled on a different route than the one the
    // reported path spells. So of a displayed path, only [route_bp, to_bp) is evidence
    // that this label carries those bases; the prefix is not.
    uint64_t route_bp = 0;
    bool entered_by_switch = false;
    LabelId from_label = 0;
    double switch_cost = 0;
    uint32_t prev_run = UINT32_MAX;
    bool ended = false;
    EndReason end_reason = EndReason::DEAD_END;
};

struct LabelEnd {
    LabelId label;
    double loss;
    uint32_t branches;
    uint32_t run;
    // 0 when this label reached the leaf along the path's own spelled sequence.
    // Non-zero when the label was merged in from another parent at a reconvergence
    // (§6.5) at that extension: the label supports the leaf and the shared bases from
    // route_bp onward, but NOT this path's spelled bases BEFORE route_bp, which belong
    // to the route it did not travel. Only labels with route_bp == 0 are evidence that
    // the label carries the whole spelled flank.
    uint64_t route_bp = 0;
};

struct Continuation {
    std::string sequence;
    std::vector<LabelId> labels;
    double loss_used = 0;
    uint32_t branches_used = 0;
};

struct PathResult {
    size_t id = 0;
    std::vector<size_t> segments;    // root -> leaf (first parent at joins)
    uint64_t length_bp = 0;
    std::array<uint32_t, kNumEndReasons> end_reasons {};
    std::optional<EndReason> path_reason;   // MAX_EXTENSION / resource / beam
    std::vector<LabelEnd> end_labels;
    std::optional<Continuation> continuation;
};

// The trie view of a split (design note §5.1.1): per branch, the labels at the
// branch's first node. A child's count need NOT be a share of the parent's: one label
// may follow several branches (an ambiguous branch in constrain mode; always possible
// in annotate mode), so the children's counts can sum to more than labels_before.
struct SplitBranch {
    size_t segment = 0;
    char ch = 0;                     // the branch's first base, walking direction
    size_t labels_distinct = 0;      // labels at the branch's first node (true count)
    std::vector<LabelId> labels;     // the first max_labels_per_node of them, ascending
};

struct Split {
    uint64_t at_bp;                  // = the shared prefix length from the seed boundary
    size_t segment;
    std::vector<size_t> children;
    bool ambiguous;                  // some label followed several children
    size_t labels_before = 0;        // distinct labels at the split node
    std::vector<SplitBranch> branches;   // parallel to children
};

struct GrowthBin {
    uint64_t from_bp = 0;
    size_t max_live_paths = 0;
    size_t max_live_labels = 0;      // distinct labels alive in the bin
    size_t max_live_pairs = 0;       // (path, label) pairs
    // false when a head counted in this bin carried a list cut by max_labels_per_node
    // (annotate mode): the two counts above are then lower bounds
    bool live_labels_exact = true;
    uint64_t steps = 0;
    size_t divergences = 0;
    size_t ambiguous_branches = 0;
    size_t splits = 0;
    size_t reconvergences = 0;
    size_t bubbles = 0;
    size_t tips = 0;
    size_t blocked_repeat = 0;
    std::array<uint32_t, kNumEndReasons> label_ends {};
};

struct BranchEvent {
    uint64_t at_bp;
    size_t segment;
    std::vector<char> chars;
    std::vector<size_t> labels_per_successor;
    std::vector<LabelId> ambiguous;
    std::vector<LabelId> dropped;
    size_t labels_affected;
    // A successor the walker decided NOT to follow for labels whose lineage would have
    // continued on it: the successor's base (walking direction), why, and those labels
    // (sources alive at the node, ascending). |cause| is "minority", "below_min_labels"
    // or "split_limit" (the quorum and split-limit texts of the BRANCH ends), "branch"
    // (an ambiguous source over its allowance, whose entries on every successor are
    // removed) or "loss_budget" (the only switch into the successor costs more than the
    // budget — stated whether or not the source goes on along another successor, and
    // not for a source the branch limit excluded, which no switch was priced for). Each
    // cause is a decision of one strategy knob, and a checker verifies it against that
    // knob. This is the walker's explicit per-successor refusal evidence: |dropped|
    // and |ambiguous| alone cannot tell a successor that was refused from one that was
    // followed and then deleted from the output (review round 3, finding 1), and a
    // missing child says nothing about why it is missing.
    struct Refusal {
        char ch;
        const char *cause;
        std::vector<LabelId> labels;
    };
    std::vector<Refusal> refused;
};

struct CapTrigger {
    EndReason reason;
    uint64_t at_bp;
    size_t segment;
    size_t live_paths;
    size_t live_labels;
    bool live_labels_exact;          // see GrowthBin::live_labels_exact
    // What the cap compared at the trip, which exceeded its limit: the seed's steps with
    // the step about to be taken (max_steps), the arm's bases likewise (max_output_bp),
    // leaves plus live heads (max_paths), live heads (max_live_paths; the beam: the
    // heads of the level before pruning), elapsed milliseconds (time_budget_ms). So a
    // response can state how far beyond the knob the run needed to go, not only that it
    // stopped (live_paths counts the heads the cap cut, which can be below the limit).
    double demand = 0;
};

struct ArmResult {
    enum Status { COMPLETE, TRUNCATED, PRUNED };

    Arm arm = Arm::RIGHT;
    bool requested = true;
    Status status = COMPLETE;
    size_t frontier_live_paths = 0;  // remaining when stopped
    size_t frontier_live_labels = 0;
    // false when a remaining head carried a cut list (annotate mode): the count above
    // is then a lower bound
    bool frontier_live_labels_exact = true;
    std::optional<CapTrigger> cap_trigger;
    /**
     * The completeness boundary (spec §6.10). Exploration is level-synchronous, so a
     * cap trips between two heads of one level: the heads expanded before it have all
     * their children, the heads after it have none. The level is then PARTIAL and does
     * not count. complete_to_bp is the extension depth up to which EVERY admissible walk
     * is present: every walk of at most complete_to_bp bases from the seed boundary that
     * obeys the per-path edge-reuse rule is in |segments|, and every path that ended
     * before complete_to_bp ended for the reported semantic reason. Walks longer than
     * that may be present (the partial level's children) but nothing is claimed about
     * them. It equals Strategy::max_extension_bp exactly when status == COMPLETE.
     */
    uint64_t complete_to_bp = 0;
    // What "every admissible walk" quantifies over (§6.10): "per_path" under
    // on_reconverge keep (the walk rule's own per-path edge history), "united_history"
    // under merge, where a merge unites the edge histories of the routes it joins and
    // the set of walks present is the smaller one admissible under that union. Stated
    // as a field, not only in walk_rule's prose, so that a checker can pick the
    // termination scope it verifies against (review round 3).
    const char *completeness_scope = "per_path";
    // Per-node label lists cut by Strategy::max_labels_per_node: how many recorded
    // positions (annotate mode: the root's boundary node, every node entered and every
    // not followed successor listed on a BLOCKED / HAIRPIN event; both modes:
    // trie-view branches) lost labels, and the largest true count seen. Non-zero means
    // the recorded sets are incomplete.
    size_t nodes_labels_truncated = 0;
    size_t max_labels_at_node = 0;

    std::vector<Segment> segments;
    std::vector<Split> splits;
    std::vector<PathResult> paths;
    std::vector<LabelRun> runs;
    std::vector<GrowthBin> growth;
    // The first Strategy::max_branch_events branch events in processing order, of
    // branch_events_total. Levels are synchronous, so an arm produces its events in
    // non-decreasing at_bp, and the cut is a depth boundary like complete_to_bp:
    // branch_events_complete_to_bp is the at_bp of the FIRST event not stored
    // (UINT64_MAX when none was dropped). Every branch event at a depth below it —
    // every refusal and every ambiguity the walker decided there — is in
    // |branch_events|; at or beyond it some may be missing. Refusals are the evidence
    // for a successor not taken (§7.2), so this bounds where an omission is
    // guaranteed to carry its reason.
    std::vector<BranchEvent> branch_events;
    size_t branch_events_total = 0;
    uint64_t branch_events_complete_to_bp = std::numeric_limits<uint64_t>::max();
    std::vector<double> needed_budgets;

    uint64_t steps = 0;
    // successor enumerations consumed by the walk (independent of batch_kmers:
    // lookahead enumerations count when, and only when, the walk reaches the node)
    uint64_t successor_enumerations = 0;
    uint64_t output_bp = 0;
    uint64_t pair_evaluations = 0;
    // Edge-use tests made by the per-path edge-reuse check: for every successor whose
    // edge was taken before, the recorded uses of that edge or the path's own segments
    // (itself and its ancestors), whichever are fewer. So a check costs at most the
    // path's depth in segments, not the number of live paths that took the same edge.
    uint64_t edge_reuse_probes = 0;
    // Re-minimisations at ambiguous nodes (§6.3): once a source over its branch
    // allowance is excluded, every successor's state is derived again, and that round
    // can expose another source over the limit. The rounds AFTER the first are counted
    // (0 at an ordinary node), as a total over the arm and as the largest at one node:
    // a node needs at most |σ| of them, and a large maximum marks a pathological locus.
    uint64_t reminimisation_rounds = 0;
    size_t max_reminimisation_rounds = 0;
    // The work of recording refusals (BranchEvent::refused): per round that excludes
    // sources, the successor state entries scanned once for those sources (scanning
    // every successor once per excluded source was Θ(|σ|²) at a node where every
    // source is ambiguous: review round 4), plus the loss-budget tests of a finite
    // change cost — one per (source, successor) under a constant cost, one per target
    // under a table. Linear in |σ| + Σ|σ_v| per round under forbid and constant costs.
    uint64_t refusal_scans = 0;
    // Successor derivations (one per successor of a committed step) in which a TABLE
    // cost's source list was cut by max_switch_sources AND a cut source had a finite
    // pair cost into a target of that successor: the cases in which the cut may have
    // raised a target's loss or left it unentered (§6.3). An over-approximation — the
    // kept sources may still have been the cheapest — but never an under-count, and it
    // covers what a `switch_sources` label end cannot: a cut source that goes on along
    // another successor (or stays) does not end, so no end records the cut.
    uint64_t switch_sources_cut = 0;
};

struct LabelArmSummary {
    // Longest seed-entered run from the boundary, measured along the label's OWN route.
    // It means: there is a graph route of this length from the seed on which EVERY k-mer
    // carries this label. It does NOT mean the label's sequence contains those bases
    // contiguously — a label can hold several contigs or repeat copies, so a route may
    // hop between them (only `support: trace`, or validation against the source record,
    // establishes a contiguous occurrence). It is also not a statement about any
    // particular displayed path: for that, check LabelEnd::route_bp on that path.
    uint64_t direct_bp = 0;
    uint64_t reach_bp = 0;           // furthest extent of the label's lineage
    size_t reentries = 0;
    std::vector<uint32_t> runs;      // indices into ArmResult::runs
};

struct DroppedLabel {
    std::string name;
    std::string reason;
    std::vector<std::pair<uint64_t, uint64_t>> runs;   // k-mer runs on the seed
};

struct SeedResult {
    std::string seed_id;
    std::string validated_seed_id;
    bool seed_id_mismatch = false;
    uint64_t length_bp = 0;
    uint64_t num_kmers = 0;
    // constrain mode: seed labels first, then extra. Annotate mode: every label
    // recorded along either arm, in the order first seen (num_seed_labels is 0).
    std::vector<LabelRef> label_dict;
    size_t num_seed_labels = 0;
    // The seed labels were derived from the seed (Seed::labels was empty). Then
    // |labels_dropped| is how many of them Strategy::max_seed_labels cut, with an
    // FNV-1a-64 hex digest of their names so a truncated set is identifiable (and the
    // complement runnable explicitly). Labels cut by the cap are NOT in
    // |dropped_labels|: they support the seed, they were only not taken.
    // |labels_supporting_total| is how many labels support every seed k-mer — for a
    // derived set the size of the intersection BEFORE the cap, for an explicit list
    // simply |num_seed_labels|, so that `labels_supporting_total > num_seed_labels`
    // means "truncated" in both cases and never reads as "nothing supports the seed".
    bool labels_from_seed = false;
    size_t labels_supporting_total = 0;
    size_t labels_dropped = 0;
    std::string labels_dropped_digest;
    std::vector<DroppedLabel> dropped_labels;
    std::array<ArmResult, 2> arms;             // index by Arm
    std::vector<std::array<LabelArmSummary, 2>> label_summary;   // per label_dict entry
    LabelOracle::Counters annotation_counters;
    double elapsed_seconds = 0;
    const char *access_path = "";
};

/**
 * Reject a strategy whose knobs contradict each other (spec §6.9): `exhaustive` with a
 * branch limit, a split limit, reconvergence merging, a quorum, a beam, or a bounded
 * max_switch_sources under a TABLE cost (a cut source list can leave a target label
 * unentered and so prune a label-consistent walk); annotate mode with label machinery
 * (extra labels, a loss budget, a switch cost, a branch limit, a quorum, trace
 * support). Throws std::invalid_argument naming the knob and the value the preset or
 * mode requires. traverse_seed() calls it; the JSON layer calls it after filling the
 * preset's values in for omitted knobs.
 */
void validate_strategy(const Strategy &strategy, const LabelChangeCost &cost);

/**
 * What "every walk is present" quantifies over (spec §6.10): the per-path edge-reuse
 * rule, hairpin handling and seed re-entry, as the walker enforces them on THIS
 * index — the regime decides whether seed nodes of both strands count, k whether a
 * self-reverse-complementary k-mer node is a hairpin too (even k), and k with the
 * alphabet whether edges are compared exactly or by a 128-bit hash. Stated in the
 * output so that a reader knows which walks complete_to_bp covers — the rule is
 * what makes "all walks" finite in a graph with cycles.
 */
std::string walk_rule_statement(const Strategy &strategy, const LabelOracle &oracle);

/**
 * Validate the seed and extend it in both directions under |strategy|.
 * Throws std::invalid_argument for invalid seeds / labels / strategies,
 * SeedDerivationError (a subclass of it) when an empty |seed.labels| yields no usable
 * permitted set, and std::runtime_error (naming the node) when the primary graph of a
 * PRIMARY index is inconsistent (CanonicalDBG found both strands of a k-mer).
 * |cost| is defined over request indices (seed labels in the given order, then
 * strategy.extra, see LabelChangeCost); SeedResult::label_dict lists the kept seed
 * labels in the given order, then strategy.extra.
 * An empty |seed.labels| derives the permitted set from the seed (see Seed::labels);
 * the derived labels then take the place of the given ones in both orders.
 */
SeedResult traverse_seed(LabelOracle &oracle,
                         const Seed &seed,
                         const Strategy &strategy,
                         const LabelChangeCost &cost,
                         const std::string &release_id = "");

// Reconstruct the flank of |path| in natural orientation from the segments.
std::string spell_path(const ArmResult &arm, const PathResult &path);

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_WALKER_HPP__
