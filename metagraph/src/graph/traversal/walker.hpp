#ifndef __TRAVERSAL_WALKER_HPP__
#define __TRAVERSAL_WALKER_HPP__

#include <array>
#include <limits>
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
 * Label-constrained extension of a seed.
 *
 * A seed is a sequence whose k-mers are all in the graph, and a set of permitted labels that
 * carry every one of them: named explicitly, or derived as the labels on every seed k-mer.
 * From the seed's last k-mer (right arm) and its first (left arm) the walker extends level
 * by level, one base per level, keeping per path the labels that carry every k-mer walked
 * so far. A label leaves a path when no admissible successor carries it (LABEL_LOST, or
 * DEAD_END without successors), when its continuation is blocked by the walk rule
 * (EDGE_REUSE, EDGE_REUSE_RC, REACHED_SEED), or when it is on several successors beyond
 * its branch allowance (BRANCH). A path splits where its labels continue on several
 * successors. The result is a tree of segments per arm.
 *
 * The walk rule, which makes the set of walks finite on a graph with cycles: a path never
 * takes an edge ((k+1)-mer) twice, in either orientation in canonical graphs, never steps
 * back into the seed, and never takes a self-reverse-complementary step (a hairpin).
 *
 * Completeness: levels are synchronous, so a cap trips between two heads of one level.
 * ArmResult::complete_to_bp is the extension up to which every admissible walk is present.
 */
struct Strategy {
    enum Direction { BOTH, LEFT, RIGHT };

    // "no limit" for max_label_branches
    static constexpr size_t kUnlimited = std::numeric_limits<size_t>::max();

    Direction direction = BOTH;

    // Cap on a permitted set DERIVED from the seed: the labels above it (in ascending
    // (column, seq_id) order) are reported as SeedResult::labels_dropped and its digest.
    size_t max_seed_labels = 1000;
    // What a derived label denotes. Unset: HEADER when the index has a CoordToHeader (the
    // labels are the indexed sequences), COLUMN otherwise.
    std::optional<LabelKind> seed_label_kind;

    // How many times a label may continue on more than one successor (0: it ends with
    // BRANCH at its first fork; kUnlimited: every label-consistent walk is followed)
    size_t max_label_branches = 0;

    uint64_t max_extension_bp = 5000;
    double time_budget_ms = 30'000;

    // output
    bool sequences = true;

    // Safety caps. Not request fields: they keep a request within bounded memory.
    size_t max_splits_per_path = 64;
    uint64_t max_steps = 1'000'000;
    size_t max_live_paths = 1000;
    size_t max_paths = 10'000;
    uint64_t max_output_bp = 2'000'000;
};

struct Seed {
    std::string seed_id;
    std::string sequence;
    // The permitted labels. Empty: derive them from the seed (the labels on every k-mer).
    std::vector<std::string> labels;
};

// A seed with a k-mer the graph does not have
class SeedNotInGraph : public std::invalid_argument {
  public:
    SeedNotInGraph(const std::string &what, uint64_t kmers, uint64_t present)
          : std::invalid_argument(what), kmers_(kmers), present_(present) {}
    uint64_t kmers() const { return kmers_; }
    uint64_t present() const { return present_; }

  private:
    uint64_t kmers_;
    uint64_t present_;
};

/**
 * The permitted set could not be derived from a seed (Seed::labels empty). A property of
 * the seed and the data, not of the request: a multi-seed request reports it per seed.
 */
class SeedDerivationError : public std::invalid_argument {
  public:
    enum Cause {
        // no label carries every k-mer of the seed: |limit| the seed's k-mers, |observed|
        // the k-mers read when the intersection became empty
        NO_CARRIER,
        // the narrowest of the first 64 k-mers has more annotation entries than the
        // derivation materialises: |limit| that bound, |observed| its entries
        TOO_WIDE,
        // the time budget ran out before every k-mer was read: |limit| the budget,
        // |observed| the k-mers read
        TIME_BUDGET,
        // a derived header name is not resubmittable as an explicit label: |subject|
        AMBIGUOUS_HEADER,
    };

    SeedDerivationError(Cause cause, const std::string &what, double limit = 0,
                        double observed = 0, std::string subject = "")
          : std::invalid_argument(what), cause_(cause), limit_(limit), observed_(observed),
            subject_(std::move(subject)) {}

    Cause cause() const { return cause_; }
    double limit() const { return limit_; }
    double observed() const { return observed_; }
    const std::string& subject() const { return subject_; }

  private:
    Cause cause_;
    double limit_;
    double observed_;
    std::string subject_;
};
const char* to_string(SeedDerivationError::Cause cause);

// A label leaving a path at |at_bp| (bases from the seed boundary)
struct LabelEnd {
    uint64_t at_bp = 0;
    LabelId label = 0;
    EndReason reason = EndReason::DEAD_END;
    // why a BRANCH or DEAD_END: "split_limit" (the path's split allowance ran out) or
    // "hairpin" (the label's only continuation was a skipped hairpin); otherwise empty
    const char *text = "";
};

struct Segment {
    size_t id = 0;
    std::vector<size_t> parents;     // none for the root, one otherwise
    std::vector<size_t> children;
    uint64_t from_bp = 0;
    uint64_t length_bp = 0;
    std::string sequence;            // bases added, natural orientation
    // the labels live at the segment's entry node and at its last node
    std::vector<LabelId> labels_start;
    std::vector<LabelId> labels_end;
    std::vector<LabelEnd> label_ends;   // sorted by at_bp
};

struct Split {
    uint64_t at_bp = 0;              // = the shared prefix length from the seed boundary
    size_t segment = 0;
    std::vector<size_t> children;
    bool ambiguous = false;          // some label followed several children
};

struct PathEnd {
    LabelId label = 0;
    EndReason reason = EndReason::DEAD_END;
    uint32_t branches = 0;           // the label's forks with several continuations taken
};

// One leaf of the segment tree; the path is its leaf's chain of parents
struct PathResult {
    size_t id = 0;
    size_t leaf = 0;                 // the leaf segment
    uint64_t length_bp = 0;
    // why the path itself ended when no label decided it: the extension length or a cap
    std::optional<EndReason> path_reason;
    std::vector<PathEnd> end_labels;
};

struct CapTrigger {
    EndReason reason;
    uint64_t at_bp;                  // depth of the first head not expanded
    // what the cap compared at the trip, which exceeded its limit (steps, bases, paths,
    // live paths, or elapsed milliseconds)
    double demand = 0;
};

struct ArmResult {
    enum Status { COMPLETE, TRUNCATED };

    Arm arm = Arm::RIGHT;
    bool requested = true;
    Status status = COMPLETE;
    size_t frontier_live_paths = 0;  // heads below the extension length when stopped
    std::optional<CapTrigger> cap_trigger;
    // Every admissible walk of at most complete_to_bp bases from the seed boundary is in
    // |segments|, and every path that ended before it ended for the reported reason.
    // Equals Strategy::max_extension_bp exactly when status == COMPLETE.
    uint64_t complete_to_bp = 0;

    std::vector<Segment> segments;
    std::vector<Split> splits;
    std::vector<PathResult> paths;

    uint64_t steps = 0;
    uint64_t output_bp = 0;
};

// A seed label that does not carry every seed k-mer
struct DroppedLabel {
    std::string name;
    std::vector<std::pair<uint64_t, uint64_t>> runs;   // the k-mer runs it carries
};

struct SeedResult {
    std::string seed_id;
    std::string validated_seed_id;
    bool seed_id_mismatch = false;
    uint64_t length_bp = 0;
    uint64_t num_kmers = 0;
    // the permitted labels; a label id is an index here
    std::vector<LabelRef> label_dict;
    // The labels were derived from the seed. Then |labels_dropped| is how many of them
    // max_seed_labels cut, with an FNV-1a-64 hex digest of their names.
    // |labels_supporting_total| is how many labels carry every seed k-mer.
    bool labels_from_seed = false;
    size_t labels_supporting_total = 0;
    size_t labels_dropped = 0;
    std::string labels_dropped_digest;
    std::vector<DroppedLabel> dropped_labels;
    std::array<ArmResult, 2> arms;   // indexed by Arm
    LabelOracle::Counters annotation_counters;
    double elapsed_seconds = 0;
};

/**
 * Validate the seed and extend it under |strategy|.
 * Throws SeedNotInGraph for a seed with a k-mer the graph does not have,
 * SeedDerivationError when an empty |seed.labels| yields no permitted set, and
 * std::invalid_argument for invalid seeds, labels or strategies.
 */
SeedResult traverse_seed(const LabelOracle &oracle,
                         const Seed &seed,
                         const Strategy &strategy);

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_WALKER_HPP__
