#ifndef __TRAVERSAL_WALKER_HPP__
#define __TRAVERSAL_WALKER_HPP__

#include <array>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
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
        // bounds.time_budget_ms ran out before the first seed k-mer was consumed: |limit| the
        // budget, |observed| the elapsed ms. Once one was, the seed is delivered with the set
        // derived from the k-mers read instead (SeedResult::derivation_partial, decision D3) —
        // unless that set would then fail the seed, which then fails with this, as before D3
        // (derivation_out_of_time)
        TIME_BUDGET,
        // a derived header name is not resubmittable as an explicit label: |subject|
        AMBIGUOUS_HEADER,
        // `exhaustive` refuses to cut the derived set: |limit| max_seed_labels,
        // |observed| the labels carrying the seed
        OVER_SEED_LABEL_CAP,
        // the seed's label dictionary (or its dropped labels) would hold a name that is
        // not valid UTF-8, which no output carries verbatim: never thrown by the walker;
        // the request layer (cli/traverse.cpp) checks the result after the walk (an
        // annotate dictionary is known only then) and reports the seed failed with this
        // cause, in both modes and whether or not the labels were derived (spec §6.1
        // step 4, §7.0)
        UNREPRESENTABLE_LABEL_NAME,
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

    // Under a memory budget, the soft excess the seed phase was seen to hold before it
    // failed (ResourceAccount::soft_overshoot, bytes): the derivation decodes whole rows
    // that no admission charges, so a failed seed states memory_bound_soft like every
    // other result under a budget (GPT review of stage 2, finding 8). Set by the walker.
    uint64_t soft_overshoot() const { return soft_overshoot_; }
    void set_soft_overshoot(uint64_t bytes) { soft_overshoot_ = bytes; }

  private:
    Cause cause_;
    double limit_;
    double observed_;
    std::string subject_;
    size_t labels_cut_;
    uint64_t soft_overshoot_ = 0;
};
const char* to_string(SeedDerivationError::Cause cause);

// The failure of a derived seed whose time budget (|budget_ms|) ran out after |kmers_read| of
// its |num_kmers| k-mers, |elapsed_ms| into the seed: every such seed's before decision D3,
// and still the one of a seed that read none, or whose partial set cannot be delivered — a
// check written for the whole seed's carriers refuses their superset (`exhaustive` over
// max_seed_labels, an ambiguous derived header, an extra label it duplicates, a label name no
// output carries), or a budget does not hold its depth-0 state. Their statements would be false
// or overstated about a set the whole seed may narrow; this one is true (review of W1, finding 1)
SeedDerivationError derivation_out_of_time(uint64_t kmers_read, uint64_t num_kmers,
                                           double budget_ms, double elapsed_ms);

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

/**
 * What the response will cost to deliver, per object, in the requested output detail
 * (bytes; upper bounds of what the serialisers hold at once: the text or the JSON tree
 * with its text, in every copy alive together). The memory budget charges every committed
 * head with the delivery of what it creates (DESIGN-traverse-graphlet.md §14, "delivery is
 * accounted per expansion"), so that a committed prefix stays deliverable within the
 * budget. All zero (the default) charges nothing for the output: a C++ caller that keeps
 * the SeedResult pays only for it. The CLI fills it in from output.detail and
 * output.sequences (cli::delivery_costs, where each bound is derived).
 */
struct DeliveryCosts {
    uint64_t fixed = 0;              // per seed: the envelope, fixed records, limitations
    uint64_t label = 0;              // per dictionary label (its record or objects), name apart
    uint64_t segment = 0;            // with its first parent's id
    uint64_t segment_label = 0;      // per label listed by a segment or a split branch
    // per parent of a merged segment: its id (beyond the first) and its labels_via_parent
    // list (the labels in it are segment_labels)
    uint64_t merge_parent = 0;
    uint64_t base = 0;               // per base of a segment
    uint64_t continuation_base = 0;  // per base of a continuation
    uint64_t run = 0;
    uint64_t event = 0;
    uint64_t event_label = 0;
    uint64_t leaf = 0;               // per path (T / C records, the path object)
    uint64_t leaf_label = 0;         // per end label and continuation label
    uint64_t chain_entry = 0;        // per segment of a path's chain (JSON paths[].segments)
    uint64_t split = 0;
    uint64_t split_branch = 0;
    uint64_t branch_event = 0;
    uint64_t branch_event_entry = 0; // per successor, ambiguous / dropped label, refusal label
    uint64_t refusal = 0;            // per refusal of a branch event, its labels apart
    uint64_t presence_run = 0;       // annotate mode: per recorded LabelSetRun
    uint64_t bin = 0;
    uint64_t dropped = 0;            // per dropped seed label (its record or object), name apart
    uint64_t dropped_run = 0;        // per k-mer presence run of a dropped seed label
    // Record coordinates (Strategy::coordinates; DESIGN-traverse-graphlet.md §18): a run's
    // entry of the coordinates block, a seed label's entry, and one occurrence of either list.
    // Charged only when coordinates are recorded, so that nothing else's account changes
    uint64_t coordinate_run = 0;
    uint64_t coordinate_seed = 0;
    uint64_t occurrence = 0;
    // The part of |fixed| that is the coordinate output's (the block's skeleton with a cut
    // list's limitation and K record, or the null form with its reason; 0 without coordinates):
    // charged with |fixed| as before, and counted again only in the account's coordinate share
    // (ResourceAccount::coordinates), so that the server's delivery reserve can price the whole
    // coordinate output at its own ratio and leave it out of the measured one (plan revision 3)
    uint64_t coordinate_fixed = 0;
    // A seed-level limitation beyond the ones |fixed| holds, charged where a result states one
    // (the derivation limitation of a permitted set derived from part of the seed: Walker,
    // SeedResult::derivation_partial)
    uint64_t extra_limitation = 0;
    /**
     * What delivering one name costs, in every copy the serialisers hold at once: label
     * names, dropped labels' names and the request's seed_id are the one input whose
     * delivered size is no fixed multiple of its length — JSON writes a control character
     * as six bytes, MGT a '%' as three, and a graphlet's MGT text is escaped once more
     * inside its JSON string — so the layer that does the escaping prices each one
     * exactly (GPT review of stage 2, finding 1: a fixed four bytes per name byte let a
     * name of control characters deliver more than the whole budget). Empty: nothing.
     */
    enum class Name { SEED_LABEL, LABEL, DROPPED_LABEL, SEED_ID };
    std::function<uint64_t(std::string_view name, Name use)> name;
};

// The interval, in charged work units, at which the walker reads the clock for its
// deadline at the latest (DESIGN-traverse-graphlet.md §14: "checked at least every W work
// units"), and within which a seed phase that ran past the work budget is let finish (so
// that the result complete to 0 bp is delivered); stated by the server as
// work_check_interval. The work budget itself is compared after every charge of the walk
// (ResourceAccount::largest_charge).
constexpr uint64_t kWorkCheckInterval = 65536;

// The graph steps of the structural lookahead (§8.3) after which it reads the seed's deadline
// and the attempt's stop again: its chains run between two checkpoints, up to
// min(batch_kmers, the radius left) steps per head of a level, and charge no units while they
// are enumerated, so the work interval never polls them (the review of 2026-10-06, W3).
// Stated by the server's deadline_check rule
constexpr size_t kLookaheadPollSteps = 16;

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
    // The request's budgets (DESIGN-traverse-graphlet.md §14), 0 = none, both PER SEED (the
    // design's locus scope): each seed of a request is walked under budgets of its own, so
    // that its result does not depend on the other seeds of the request (§6.8). A request
    // of n seeds can therefore hold up to n times the memory budget — each seed's output is
    // kept until the response is written — and the server bounds n (review of stage 2,
    // finding 7: a request-level ledger would make a seed's stop depend on its position).
    // Memory: the modelled bytes one seed's walk and output retain at once — the walker's
    // state, the label caches (a fixed allotment under a budget) and the output in the
    // requested detail (|delivery|) — charged per head at admission, the depth-0 state
    // included (SeedBudgetError), never measured, so that a stop is reproducible. It is the
    // one bound on the size of the output, and so on the time to serialise it: detail full
    // spells every leaf's chain, quadratic on a comb-shaped trie. Work: charged work units,
    // a weighted sum of the walk's counters (successor enumerations, the annotation rows its
    // fetches return with their entries and coordinates, pair evaluations, refusal scans,
    // edge-reuse probes, derivation scans, steps) and of the seed phase's reads, compared
    // after every charge of the walk (the seed phase: once every kWorkCheckInterval units),
    // so a stop exceeds it by what was charged since the previous comparison — one fetch
    // call's rows, a node's label-state scan, the seed phase with the roots' rows — and the
    // largest such charge is stated (ResourceAccount::largest_charge). Work
    // is the walk's alone, the same in every output detail (finding 5: charging delivery as
    // work would make a work stop depend on the detail). Either budget stops the whole
    // seed, with resource_limit.
    uint64_t max_memory_bytes = 0;
    uint64_t max_work_units = 0;
    DeliveryCosts delivery;

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
    // Record coordinates (DESIGN-traverse-graphlet.md §18, owner decisions C1-C12; opt-in, C1):
    // where each run's sequence and each seed label's seed lies in its label's record (header
    // labels) or column (column labels), recorded only under support trace on an index with
    // k-mer coordinates — everywhere else coordinates_reason says why none are (§18.1). Off (the
    // default), nothing is recorded or charged: every result and every budgeted stop is the same
    // as without the feature. A list keeps the first |max_coordinate_occurrences| occurrences by
    // start (kUnlimited: all; at least 1) and states its true count (RunCoordinates::total).
    bool coordinates = false;
    size_t max_coordinate_occurrences = 16;

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

// NOTE: sizeof(LabelRun) is part of the memory model (Walker::init_budgets charges every run
// twice its size), as are those of the walker's Entry and Item: a field added here changes
// every budgeted account, and with it where every budgeted request stops, whether or not it
// asks for the field. Per-run output that only some requests ask for lives in a side table
// beside |runs| (ArmResult::run_coordinates), charged only when recorded.
struct LabelRun {
    LabelId label = 0;
    uint64_t from_bp = 0;
    uint64_t to_bp = 0;              // half-open; set when the run ends
    // Non-zero when this run's lineage was merged in from another parent at a
    // reconvergence (§6.5): the EARLIEST such merge on its route (runs are shared across
    // sibling paths, so the first stamp wins). A merge joins paths AT A NODE, so the bases
    // before it were travelled on a different route than the one the reported path
    // spells. It is not where displayed support begins: after two such merges d1 < d2 the
    // run keeps d1 while the displayed path spells another route up to d2 (spec §7.1), so
    // displayed support is derived from the merge partitions (evidence_from, design §5.1);
    // only route_bp == 0 says the whole run is displayed support. One exception, also
    // stated in §7.1: at a merge of three or more parents an entry that beat an earlier
    // parent is stamped before a later parent beats it, so the run it closes there (to_bp ==
    // route_bp, `m`) carries that merge, which is not on its route; such a run's displayed
    // support is still [evidence_from, to_bp). Changing the stamp would change R in
    // unbudgeted output, so it stays, stated.
    uint64_t route_bp = 0;
    bool entered_by_switch = false;
    LabelId from_label = 0;
    double switch_cost = 0;
    uint32_t prev_run = UINT32_MAX;
    bool ended = false;
    EndReason end_reason = EndReason::DEAD_END;
    // The segment on which the run ended (its LABEL_END event, if any, is there at
    // to_bp) or was closed by a reconvergence merge (a parent of the merged segment).
    // Neither the events nor the paths determine it: clones made at a split are
    // identical rows, a switch source that goes on only under other names ends
    // silently, and a merge of three or more parents closes several runs at one depth.
    uint32_t segment = UINT32_MAX;
    // The lineage's branch count and loss when the run ended or was closed. Branches
    // are counted before quorum filtering, so they can depend on a successor that left
    // no trace in the output; the leaf records them only for labels alive at a leaf.
    uint32_t branches = 0;
    double loss = 0;
};

/**
 * Record coordinates of one run (Strategy::coordinates; DESIGN-traverse-graphlet.md §18.2):
 * ArmResult::run_coordinates[i] belongs to ArmResult::runs[i]. An occurrence is a coordinate
 * chain of the run's label live at the run's last node — one copy of the run's sequence in the
 * label's record (header label) or column (column label), followed k-mer by k-mer: chains only
 * continue or die along a run, and a split hands each chain to at most one child, since a
 * coordinate is one k-mer. |ends| holds such a chain's k-mer coordinate at that last node; the
 * occurrence's interval follows from the run's length L = to_bp - from_bp and k: [c + k - L,
 * c + k) on the right arm, [c, c + L) on the left (the run's own bases; L = 0 is the empty
 * interval at the seed boundary).
 */
struct RunCoordinates {
    // the first max_coordinate_occurrences chains' coordinates, ascending (so the intervals'
    // starts ascend on both arms)
    std::vector<Coord> ends;
    // the chains live at the last node: the run's true occurrence count (> |ends| when cut)
    uint64_t total = 0;
    // Chains that were live on the run's lineage and continued on no followed path of it
    // before its last node (C5, decision C-N1): a record end, a successor not followed
    // (blocked, a quorum, the branch limit), or a continuation taken over by a switch. Chains
    // partitioned among the children of a split are not ended (each follows one child); a
    // clone made at a split inherits the count of the prefix it shares.
    uint64_t chains_ended = 0;
    // The run was entered by a switch into a label whose own lineage was still live at the
    // switch (revision 1 of the coordinates plan, decision C-N2 option a): its entry kept only
    // the label's chains that continue from before the switch, so chains of the label that
    // start at the switch node are missing — |total| is a lower bound. Stated per run and by
    // the block (complete: false); the walk itself is unchanged.
    bool lower_bound = false;
};

// Record coordinates of one seed label: its occurrences of the whole seed, each the chain's
// coordinate at the first seed k-mer (the interval [c, c + |seed|)); SeedResult::seed_coordinates
struct SeedCoordinates {
    std::vector<Coord> starts;       // the first max_coordinate_occurrences, ascending
    uint64_t total = 0;              // the label's chains over the whole seed
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

/**
 * One leaf of the segment DAG. The path is its leaf's first-parent chain, which is NOT
 * stored: on a comb-shaped trie (one terminating branch per split) the chains sum to
 * Θ(segments²), so materialising one per leaf made finalisation quadratic in the result
 * (DESIGN-traverse-graphlet.md §14, "delivery is accounted per expansion"). The chain is
 * derived on demand from parent pointers — by the serialisers while they write it, and
 * by path_segments() for C++ callers.
 */
struct PathResult {
    size_t id = 0;
    size_t leaf = 0;                 // the leaf segment (root -> leaf: path_segments())
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
    // Record coordinates (SeedResult::coordinates_recorded): one entry per run, parallel to
    // |runs|; empty otherwise. A side table, not a field of LabelRun, whose size the memory
    // model charges for every run of every request (see LabelRun)
    std::vector<RunCoordinates> run_coordinates;
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
    // The charged work units of this arm (Strategy::max_work_units): the weighted sum of
    // the work its heads consumed (successors enumerated 4, annotation keys requested 8,
    // entries returned 1, pair evaluations, refusal scans, edge-reuse probes, derivation
    // scans and steps 1 each). Independent of batch_kmers, like the counters it weighs.
    uint64_t work_units = 0;
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

/**
 * Why the walk of a seed stopped at a budget (DESIGN-traverse-graphlet.md §14; JSON
 * `resource_stop`, MGT `Q`): which resource, in which phase, and the amounts in the
 * resource's raw unit — bytes (memory), work units (work) or milliseconds (time). The
 * first stop of the seed; the stop applies to the whole locus (both arms).
 */
struct ResourceStop {
    // CANCELLED, ATTEMPT_DEADLINE: no budget of the request ran out — the walk was stopped from
    // outside (AttemptControl): the attempt was cancelled, or it reached the duration bound the
    // server enforces for it (DESIGN-traverse-graphlet.md §14 v5.1). Their amounts are the
    // attempt's: |limit| its bound and |used| its elapsed time (whole milliseconds)
    enum Resource { MEMORY, WORK, TIME, CANCELLED, ATTEMPT_DEADLINE };
    Resource resource = MEMORY;
    // traversal | finalisation | serialisation | annotation_decode. Finalisation and
    // serialisation are reserved at admission, so they never stop; annotation_decode: a
    // budget-aware annotation read (LabelOracle::decode_charged()) did not fit the memory
    // the request had left (stage 3 of DESIGN-traverse-graphlet.md §14.1) — a level's read
    // (the level is censored from its first head) or, failing the seed, the seed phase's
    // or an annotate root's
    const char *phase = "traversal";
    Arm arm = Arm::RIGHT;            // the arm whose head was not admitted
    uint64_t at_bp = 0;              // that head's depth
    double limit = 0;                // 0: none (a hook refused the head)
    double used = 0;                 // accounted when the head was refused
    double demand = 0;               // what admitting it needed (used, for work and time)
    // a memory refusal injected by WalkerHooks::deny or deny_decode (test hooks) while the
    // budget, if any, would have admitted the head or the read: the statement must not blame
    // a budget for it
    bool injected = false;

    // ---- what did not fit, for the statements of a memory stop (the message, the actions
    // and the walk_domain they shape; nothing of it is a field of its own on the wire). A stop
    // must name its true cause: a row that does not fit is not a dictionary or a level's lists
    // that do not (review of stage 3, F2, F6, F7)
    enum Cause {
        HEAD,           // a head's admission (traversal), or the seed's depth-0 state
        READ_ROW,       // annotation_decode: an annotation row with its row-diff dependency rows
        LABEL_NAMES,    // traversal: the dictionary labels a level's rows would name first
        LEVEL_LISTS,    // traversal: the level's own lists left nothing to read its rows with
    };
    Cause cause = HEAD;
    // where the refused read was (READ_ROW, LABEL_NAMES): the seed phase (the derivation's
    // window or the validation), an annotate root, or a level's fetch
    enum Where { SEED, DERIVATION, ROOT, LEVEL };
    Where where = LEVEL;
    uint64_t left = 0;               // bytes the read had left at the refused row
    // READ_ROW: the refused row's standalone demand (the bound every returned row is admitted
    // against); 0 when its read alone was refused, which then needed more than |left|
    uint64_t row_demand = 0;
    uint64_t held = 0;               // what was held beside the read (bytes: the account, the
                                     // level's lists and rows read before it, or the seed phase's)
    // LABEL_NAMES: the new labels, and what naming them costs; READ_ROW at an annotate ROOT:
    // the labels already in the dictionary (the other arm's root's), and what the depth-0 state
    // built before the read holds beyond the fixed state (that root with its labels)
    uint64_t labels = 0;
    uint64_t label_bytes = 0;
    uint64_t index = 0;              // SEED, DERIVATION: the seed k-mer of the refused row
    // |demand| is the least the seed was seen to need (part of what it needs was not read
    // because the budget was already exceeded): its statements say "at least"
    bool lower_bound = false;
    // MEMORY: the caches' allotments of the budget that |used| and |demand| include (bytes;
    // 0 in the seed phase, which runs before they are charged). They are fractions of the
    // budget, so a demand measured at this budget is not the budget that holds it: the
    // statements give memory_budget_holding(demand, allotted) instead
    uint64_t allotted = 0;
};

// The caches' allotments of a memory budget of |budget| bytes, charged in the account with the
// depth-0 state (Walker::init_budgets): a quarter of it, at most 64 MiB, for the label cache and
// a sixteenth, at most 64 MiB, for the lookahead
uint64_t memory_allotments(uint64_t budget);

/**
 * The smallest bounds.max_memory_mb, in bytes (a whole number of MiB), whose admitted account
 * holds |need| bytes measured at a budget whose caches' allotments were |allotted| bytes: the
 * allotments grow with the budget (memory_allotments), so a budget of |need| holds less than
 * |need| beside them and the same state fails again (review of the stage-3 fixes, P2: a depth-0
 * failure stated "need 7 MiB" at 1 MiB, then "need 9 MiB" at 7 and "need 10 MiB" at 8 and 9;
 * 10 held it). With |allotted| 0 (nothing allotted in the account yet: the seed phase) it is
 * |need| rounded up to whole MiB, as before.
 */
uint64_t memory_budget_holding(uint64_t need, uint64_t allotted);

// The request's accounts at the end of the seed (memory in modelled bytes)
struct ResourceAccount {
    uint64_t memory_limit = 0;
    uint64_t memory_peak = 0;        // the largest accounted total
    uint64_t memory_final = 0;       // accounted when the walk ended (live reservations 0)
    // The excess over the budget (0 without one) of what was held but could not be
    // charged before it was allocated — a level's decoded annotation rows (stage 3 charges
    // them inside the decoder), a cache beyond its allotment and an annotate dictionary's
    // growth — on top of the admitted account, which never exceeds the budget itself: the
    // memory_bound_soft limitation's observed
    uint64_t soft_overshoot = 0;
    uint64_t work_limit = 0;
    // the seed phase (its validation, or the derivation of its permitted set: 8 per seed
    // k-mer read, 1 per annotation entry and coordinate) and both arms' work_units
    uint64_t work_seed = 0;
    uint64_t work_used = 0;
    uint64_t work_check_interval = kWorkCheckInterval;
    // The most work charged between two comparisons with the work budget (units; 0: no
    // comparison was made, i.e. no work budget). Every comparison follows one that passed,
    // so a work stop exceeds the budget by at most its last such stretch, and so by at most
    // this: the bound a stop states with its number (review of the stage-2 fixes, F7), whatever
    // the charge was — a fetch call's annotation rows decoded whole with their coordinates,
    // a node's label-state scan, or the seed phase with the roots' rows, which are charged
    // before the first comparison so that the result complete to 0 bp can be delivered
    uint64_t largest_charge = 0;
    // How the annotation was read under the request's budgets (stage 3 of
    // DESIGN-traverse-graphlet.md §14.1): |decode_charged| — by the budget-aware decode path
    // (every dependency row and coordinate tuple charged before it is held, a read that does
    // not fit refused whole, each row's dependency rows charged as work); |row_diff_uncounted|
    // — a row-diff annotation without that path, whose dependency rows are neither charged
    // nor counted. Both false without a budget.
    bool decode_charged = false;
    bool row_diff_uncounted = false;
    // the part of the modelled memory that is the coordinate output (Strategy::coordinates:
    // the fixed part, DeliveryCosts::coordinate_fixed, and the seed's and the runs' entries and
    // occurrences as charged at their creation); 0 without them. Their text per account byte
    // differs from the rest of the output's, which the server's delivery reserve tells apart
    uint64_t coordinates = 0;
};

// A seed's deadline record (R8; timing only): its longest uninterruptible piece, and, when a
// stop ended its walk, what stopped it (time_budget, attempt, cancelled, memory, work) and how
// long after its deadline (the seed's time budget, or the attempt's walk-until) it stopped
struct DeadlineRecord {
    UninterruptiblePiece longest;
    std::string stopped_by;
    std::optional<double> after_deadline_ms;
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
    DeadlineRecord deadline;
    const char *access_path = "";
    std::optional<ResourceStop> resource_stop;
    ResourceAccount account;

    // ---- record coordinates (Strategy::coordinates; DESIGN-traverse-graphlet.md §18)
    // Recorded: |seed_coordinates| and both requested arms' run_coordinates hold them.
    // Otherwise, when the strategy asked for them, |coordinates_reason| says why none are
    // (coordinates_reason() before the walk; kCoordinatesPartialDerivation after a derivation
    // the time budget stopped). Both empty when they were not asked for.
    bool coordinates_recorded = false;
    const char *coordinates_reason = nullptr;
    std::vector<SeedCoordinates> seed_coordinates;   // per seed label (id < num_seed_labels)
    // the graph's k: what turns a chain's k-mer coordinate into the interval of a run's bases
    size_t k = 0;

    // The permitted set was derived from the first |kmers_read| of the seed's num_kmers k-mers
    // only: the seed's time budget ran out during the derivation (the owner's decision D3,
    // 2026-10-04). Such a set is the intersection over fewer rows, a superset of the labels
    // carrying the whole seed, so it may hold labels the whole seed would exclude: the result
    // states a `derivation` limitation and its label evidence is qualified (overstated). The
    // walk is the seed itself — stopped by the time budget at depth 0 (complete_to_bp 0) —
    // since the budget is spent. A derivation stopped before its first k-mer still fails the
    // seed (SeedDerivationError::TIME_BUDGET), as before, and so does one whose partial set
    // would fail the seed in any other way (derivation_out_of_time): a result carrying this is
    // always a walked one.
    struct PartialDerivation {
        uint64_t kmers_read = 0;
        double elapsed_ms = 0;
    };
    std::optional<PartialDerivation> derivation_partial;
};

// Why a seed reports no coordinates (SeedResult::coordinates_reason and the JSON's
// coordinates_reason), §18.1: the index has no k-mer coordinates; the support is k-mer
// presence (also in annotate mode, which refuses trace), where a stretch can be stitched from
// several occurrences; the seed was not traversed (failed, refused, not started); or its
// permitted set was derived from part of the seed, so no label's occurrences of the whole seed
// are known (D3, under support trace)
constexpr const char *kCoordinatesNoIndex = "index has no coordinates";
constexpr const char *kCoordinatesSupportKmer = "support kmer";
constexpr const char *kCoordinatesNoTraversal = "no traversal";
constexpr const char *kCoordinatesPartialDerivation = "partial derivation";

// The reason a seed of |strategy| reports no coordinates whatever happens to it, checked in
// the order of §18.1 (the index's, then the support's); nullptr: they are not requested, or
// a walked seed records them
const char* coordinates_reason(const Strategy &strategy, bool index_has_coordinates);

/**
 * A request budget (DESIGN-traverse-graphlet.md §14) does not hold the seed itself: the work
 * budget ran out while the seed was validated or its permitted set derived, or the memory
 * budget does not hold the depth-0 state — the label dictionary and both arms' roots, each
 * root reserved with what ending and delivering its labels costs in the requested detail
 * (an unadmitted depth-0 state was delivered whole, ~4 KB per label and arm in detail full:
 * review of stage 2, finding 3). No valid traversal exists within the budget, not even one
 * complete to 0 bp, so unlike a stop between two heads there is no partial result: the
 * caller reports the seed as failed (outcome.walks: failed) with its resource stop and
 * traverses the other seeds. Raised before anything per label is delivered.
 */
class SeedBudgetError : public std::runtime_error {
  public:
    SeedBudgetError(const std::string &what, const ResourceStop &stop,
                    const ResourceAccount &account, bool labels_from_seed, size_t labels)
          : std::runtime_error(what), stop_(stop), account_(account),
            labels_from_seed_(labels_from_seed), labels_(labels) {}

    // resource MEMORY or WORK, at_bp 0. |used|: memory — what the account held before the
    // depth-0 state (the seed and the fixed allotments); work — the seed phase's units.
    // |demand|: memory — what the depth-0 state needs in all; work — as used
    const ResourceStop& stop() const { return stop_; }
    const ResourceAccount& account() const { return account_; }
    // whether the permitted set was derived from the seed, and its size when the budget
    // refused it (0 when the work budget ran out before the set was known)
    bool labels_from_seed() const { return labels_from_seed_; }
    size_t labels() const { return labels_; }

  private:
    ResourceStop stop_;
    ResourceAccount account_;
    bool labels_from_seed_;
    size_t labels_;
};

/**
 * Reject a strategy whose knobs contradict each other (spec §6.9): `exhaustive` with a
 * branch limit, a split limit, reconvergence merging, a quorum, a beam, or a bounded
 * max_switch_sources under a TABLE cost (a cut source list can leave a target label
 * unentered and so prune a label-consistent walk); annotate mode with label machinery
 * (extra labels, a loss budget, a switch cost, a branch limit, a quorum, trace
 * support), or a cap of 0 on the recorded coordinates' lists (coordinates asked for under
 * k-mer support or in annotate mode are not refused: such a result says why it reports none,
 * coordinates_reason()). Throws std::invalid_argument naming the knob and the value the preset
 * or mode requires. traverse_seed() calls it; the JSON layer calls it after filling the
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
 * What the walker is about to commit for one head (DESIGN-traverse-graphlet.md §14,
 * "atomic commit per head"): every head is planned without touching the result, then
 * ADMITTED against the request's budgets, then committed (which cannot fail). A head
 * that is not admitted is censored like a head beyond a cap, with resource_limit.
 * |ordinal| counts the admissions of the seed so far, over both arms, in processing
 * order, so that a test can deny exactly the n-th one.
 */
struct Admission {
    Arm arm = Arm::RIGHT;
    uint64_t ordinal = 0;
    uint64_t at_bp = 0;              // the head's extension depth
    size_t segment = 0;              // the head's segment
    size_t followed = 0;             // successors the plan follows (0: the head ends)
    // the memory ledger (modelled bytes): the accounted total before the head, the
    // reservation it held (released by committing it), what committing it adds (its
    // objects and its children's reservations), and the budget (0: none)
    uint64_t total = 0;
    uint64_t reserve = 0;
    uint64_t cost = 0;
    uint64_t limit = 0;
};

/**
 * Test and instrumentation hooks of traverse_seed(). |deny| is asked at every admission
 * after the budgets admitted the head: returning true censors the head with
 * resource_limit as if an allocation had been refused — the injection point for the
 * allocation-denial fixtures of §14 (a denial around a switch, a split and a merge must
 * leave runs, events, splits and paths consistent).
 */
/**
 * One charge of a budget-aware annotation read (stage 3 of DESIGN-traverse-graphlet.md §14.1;
 * annot::matrix::DecodeBudget): where the read is — the seed phase (SEED: the derivation's
 * rows or the validation's), an annotate root (ROOT), a level's fetch (LEVEL) or the
 * lookahead (WARM) — the arm and depth, and |ordinal|, which counts the decode charges of the
 * seed so far, so that a test can deny exactly the n-th one.
 */
struct DecodeCharge {
    enum Where { SEED, ROOT, LEVEL, WARM };
    Where where = LEVEL;
    Arm arm = Arm::RIGHT;
    uint64_t at_bp = 0;
    uint64_t ordinal = 0;
};

/**
 * A stop requested from outside the walk (DESIGN-traverse-graphlet.md §14 v5.1: a lease is a
 * valid release point only because the backend stops the attempt itself): the attempt was
 * cancelled, it reached the duration bound the server enforces for it, or its client closed
 * the connection. The core library knows no HTTP (spec §10.1): the caller decides, the walker
 * only asks.
 */
enum class ExternalStop : uint8_t { NONE = 0, CANCELLED, ATTEMPT_DEADLINE, CLIENT_GONE };

// The client of the request is gone (ExternalStop::CLIENT_GONE): nothing can be delivered, so
// the walk is abandoned where it is — no result is finalised — and the caller writes nothing.
// Not an std::invalid_argument: it is no fault of the request.
class AttemptAborted : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// What the walk of one seed consumed, written when traverse_seed() returns or throws, however
// it ends (a result, a failed seed, an invalid request, an abandoned walk): what a ledger
// reconciles a reservation against (DESIGN §14, "usage on every response")
struct AttemptMeter {
    uint64_t work_units = 0;         // ResourceAccount::work_used
    uint64_t work_seed = 0;          // ResourceAccount::work_seed
    uint64_t memory_peak = 0;        // the modelled account's peak (bytes)
    uint64_t memory_final = 0;       // what the result holds when the walk ended
    uint64_t soft_excess = 0;        // ResourceAccount::soft_overshoot (0 without a budget)
    bool walked = false;             // the walker ran (false: the strategy was refused first)
    // the part of the account that is recorded coordinates (ResourceAccount::coordinates)
    uint64_t memory_coordinates = 0;
};

/**
 * The caller's control of a walk from outside (the server's attempt): |poll| is asked where
 * the walk already looks at its budgets and its deadline — before every head, at least every
 * kWorkCheckInterval charged units, at every charge of the seed phase and per k-mer of a
 * derivation — so it must be cheap. CANCELLED and ATTEMPT_DEADLINE stop the walk like a budget
 * (a valid partial result, censored with resource_limit; in the seed phase the seed fails);
 * CLIENT_GONE throws AttemptAborted. |elapsed_ms| is the attempt's clock and |bound_ms| its
 * bound: what such a stop states. |meter|, when given, receives the seed's usage.
 */
struct AttemptControl {
    std::function<ExternalStop()> poll;
    std::function<double()> elapsed_ms;
    double bound_ms = 0;
    AttemptMeter *meter = nullptr;
    // The chunked deadlines (pass 5, spec §6.8): |ms_left| the time left before the walk must
    // stop (the attempt's walk-until; infinity: none), which sizes the chunks of an annotation
    // read, and |poll_now| the poll asked before each chunk — a poll that also reads the clock
    // and the client (not only the stop flag), so that a stop that fell inside a long read is
    // seen within one chunk. Either may be null.
    std::function<double()> ms_left;
    std::function<ExternalStop()> poll_now;
    // the walk's modelled account (bytes, the memory model's, computed with or without a
    // budget) at every level's end, and the part of it that is recorded coordinates (0 without
    // them): what the caller estimates the seed's output from (the server's delivery reserve),
    // whose text per account byte differs between the two parts
    std::function<void(uint64_t account, uint64_t coordinates)> progress;
};

struct WalkerHooks {
    std::function<bool(const Admission&)> deny;
    // asked at every charge of a budget-aware annotation read (only on an index whose reads
    // are budget-aware, and only under a request budget): returning true refuses that charge
    // as if it had not fit — reported like a refusal of |deny|, as an injected memory stop, in
    // phase annotation_decode (a lookahead read gives up silently, as on a real refusal)
    std::function<bool(const DecodeCharge&)> deny_decode;
    // after every level of an arm (its merges and beam done): the arm, the level's depth
    // and the accounted memory total — what a budget must admit to complete that level
    std::function<void(Arm, uint64_t, uint64_t)> level;
    // a derived seed's state is recounted in full at every observation and compared with the
    // running total the derivation keeps (std::logic_error on a difference): what a debug
    // build asserts, so that a Release build's tests check that the total is the sum it
    // replaced (the review of 2026-10-06, W2)
    bool recount_derivation_state = false;
};

/**
 * Validate the seed and extend it in both directions under |strategy|.
 * Throws std::invalid_argument for invalid seeds / labels / strategies,
 * SeedDerivationError (a subclass of it) when an empty |seed.labels| yields no usable
 * permitted set, SeedBudgetError when a request budget does not hold the seed itself
 * (its validation, or its depth-0 state), and std::runtime_error (naming the node) when
 * the primary graph of a PRIMARY index is inconsistent (CanonicalDBG found both strands
 * of a k-mer).
 * |cost| is defined over request indices (seed labels in the given order, then
 * strategy.extra, see LabelChangeCost); SeedResult::label_dict lists the kept seed
 * labels in the given order, then strategy.extra.
 * An empty |seed.labels| derives the permitted set from the seed (see Seed::labels);
 * the derived labels then take the place of the given ones in both orders.
 * |control|: a stop from outside the walk (AttemptControl); null: none, as before.
 */
SeedResult traverse_seed(LabelOracle &oracle,
                         const Seed &seed,
                         const Strategy &strategy,
                         const LabelChangeCost &cost,
                         const std::string &release_id = "",
                         const WalkerHooks *hooks = nullptr,
                         const AttemptControl *control = nullptr);

/**
 * The cheapest loss at which a chain of switches from any of the labels [0, |sources|)
 * enters each label of [0, |n|) under |cost|, for every label it enters within |budget|
 * (+inf for the others): a multi-source shortest path over the cost model, summed the way
 * the walk sums a lineage's loss (left to right along the chain, so a label the walk can
 * enter within the budget is found within it). An extra label is a switch target; the walk
 * enforces cumulative losses switch by switch, so it is accepted when SOME chain reaches it
 * within the budget, not only one switch from a seed label (the stage-2 recheck's design
 * answer: A -> B = 1, B -> C = 1 under a budget of 2 rejected C, which the walk enters).
 *
 * A constant cost reaches every label in one switch (a chain costs more), so only
 * cost <= budget matters. A table is relaxed by Dijkstra: a popped label relaxes its
 * explicit entries, and the default (when finite) is relaxed lazily — the labels no popped
 * label has relaxed by default yet are kept in |pending|, and each pop gives every one of
 * them without an explicit entry from the popped label the default from it, then drops it:
 * pops come in ascending loss, so the first default a label receives is its cheapest. A
 * label is skipped there only for an explicit entry, so the work is O((n + entries) log n),
 * not O(n^2) for a pool of thousands of extra labels.
 */
std::vector<double> switch_reach(const LabelChangeCost &cost, size_t n, size_t sources,
                                 double budget);

// Reconstruct the flank of |path| in natural orientation from the segments.
std::string spell_path(const ArmResult &arm, const PathResult &path);

// The segments of |path|, root -> leaf, through first parents at joins: the chain that
// PathResult no longer stores. O(path depth in segments) per call.
std::vector<size_t> path_segments(const ArmResult &arm, const PathResult &path);

// Visit the segments of |path| leaf -> root (first parents) without allocating; |f|
// returns false to stop early (a continuation needs only the last n bases).
template <class F>
void walk_path_leaf_first(const ArmResult &arm, const PathResult &path, F f) {
    for (size_t s = path.leaf; ; s = arm.segments[s].parents[0]) {
        if (!f(s) || arm.segments[s].parents.empty())
            return;
    }
}

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_WALKER_HPP__
