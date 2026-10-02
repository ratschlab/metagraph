#ifndef __TRAVERSAL_WALKER_HPP__
#define __TRAVERSAL_WALKER_HPP__

#include <array>
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
 */
class SeedDerivationError : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

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


struct Strategy {
    enum Direction { BOTH, LEFT, RIGHT };
    enum Order { BREADTH_FIRST, LOWEST_LOSS_FIRST, MOST_SUPPORTED_FIRST };
    enum Overflow { STOP, BEAM };

    Direction direction = BOTH;
    Support support = Support::KMER;

    // labels
    std::vector<std::string> extra;
    double loss_budget = 0;
    bool switch_on_loss_only = true;
    size_t max_switch_sources = 64;
    // Cap on a set DERIVED from the seed (Seed::labels empty). It bounds the
    // traversal state, not the derivation: the labels above the cap are discovered
    // either way and reported as SeedResult::labels_dropped / _digest.
    size_t max_seed_labels = 1000;
    // What a derived label denotes. Unset (the default) means HEADER when the index
    // has a CoordToHeader (labels are the indexed sequences / accessions) and COLUMN
    // otherwise. Ignored when Seed::labels is given.
    std::optional<LabelKind> seed_label_kind;

    // branching
    size_t max_label_branches = 0;
    size_t min_successor_labels = 1;
    double min_successor_fraction = 0;
    uint64_t tip_window_bp = 0;
    uint64_t bubble_window_bp = 0;
    bool merge_reconverge = true;
    bool skip_hairpins = true;
    size_t max_splits_per_path = 64;

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
    size_t max_branch_events = 100;
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
    std::string text;
};

struct Segment {
    size_t id = 0;
    std::vector<size_t> parents;
    // joins only (|parents| > 1): per parent, the labels whose kept entry entered
    // through it, so that per-label routes through the DAG stay reconstructible
    std::vector<std::vector<LabelId>> labels_via_parent;
    std::vector<size_t> children;
    uint64_t from_bp = 0;
    uint64_t length_bp = 0;
    std::string sequence;            // bases added, natural orientation
    std::vector<LabelId> labels_start;
    std::vector<LabelId> labels_end;
    std::vector<Event> events;       // sorted by at_bp
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

struct Split {
    uint64_t at_bp;
    size_t segment;
    std::vector<size_t> children;
    bool ambiguous;                  // some label followed several children
};

struct GrowthBin {
    uint64_t from_bp = 0;
    size_t max_live_paths = 0;
    size_t max_live_labels = 0;      // distinct labels alive in the bin
    size_t max_live_pairs = 0;       // (path, label) pairs
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
};

struct CapTrigger {
    EndReason reason;
    uint64_t at_bp;
    size_t segment;
    size_t live_paths;
    size_t live_labels;
};

struct ArmResult {
    enum Status { COMPLETE, TRUNCATED, PRUNED };

    Arm arm = Arm::RIGHT;
    bool requested = true;
    Status status = COMPLETE;
    size_t frontier_live_paths = 0;  // remaining when stopped
    size_t frontier_live_labels = 0;
    std::optional<CapTrigger> cap_trigger;

    std::vector<Segment> segments;
    std::vector<Split> splits;
    std::vector<PathResult> paths;
    std::vector<LabelRun> runs;
    std::vector<GrowthBin> growth;
    std::vector<BranchEvent> branch_events;
    size_t branch_events_total = 0;
    std::vector<double> needed_budgets;

    uint64_t steps = 0;
    // successor enumerations consumed by the walk (independent of batch_kmers:
    // lookahead enumerations count when, and only when, the walk reaches the node)
    uint64_t successor_enumerations = 0;
    uint64_t output_bp = 0;
    uint64_t pair_evaluations = 0;
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
    std::vector<LabelRef> label_dict;          // seed labels first, then extra
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
