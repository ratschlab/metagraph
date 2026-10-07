#ifndef __PATTERN_SEARCH_HPP__
#define __PATTERN_SEARCH_HPP__

/**
 * Pattern search: count, and extract without reading any annotation, every graph context of
 * a short motif or an IUPAC pattern (peptides later) on the succinct graph.
 *
 * The design is docs/DESIGN-pattern-search.md (v6 + §22); section numbers below refer to it.
 * This header is the contract between the engine (pattern_search.cpp) and the route
 * (src/cli/pattern.cpp) for increments 0-2: modes count, all_or_count and partial, the last two
 * only with output.labels "none" (the label-free path, §4.3); no annotation read, no placement,
 * no predicate, no extension beyond k, single graph. The JSON the route writes from these
 * types is the contract of §7 (src/cli/pattern.cpp writes it); every JSON name comes from the
 * to_string() family at the end of this file.
 *
 * The owner's guarantee rule holds for every type here: nothing is weakened silently, every
 * count carries its unit and its relation, and a count is never promoted by assumption.
 */

#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct;

namespace pattern {

// The version of the /pattern request and answer JSON (capabilities:
// pattern.pattern_contract_version). Raised when a field changes meaning; a field or value
// that only becomes accepted (a later increment's) does not raise it.
constexpr int kPatternContractVersion = 1;


// ---------------------------------------------------------------- counts

/**
 * What a count counts (§3). Units are never summed into each other: a graph-context count
 * says nothing about placed occurrences, an anchor count nothing about completed paths.
 * PLACED_OCCURRENCES and LABELS need annotation (a later increment): this increment states
 * them UNKNOWN.
 */
enum class Unit { GRAPH_CONTEXTS, ANCHORS, PATHS, PLACED_OCCURRENCES, LABELS };

/**
 * How a count's value relates to the true number (§3, §4.1).
 *  EXACT     the discovery behind it completed: no step, time or threshold stop touched it.
 *  AT_LEAST  a stop interrupted discovery; the value is the sum of the lower bounds of what
 *            was explored. An undiscovered IUPAC branch, offset or strand has no upper bound,
 *            so this is the relation of every interrupted discovery.
 *  BOUNDS    every range of every branch, offset and orientation was discovered and only the
 *            scans of invalid (masked) edges were interrupted: lower <= true <= upper.
 *  UNKNOWN   the phase never ran: not started, after a stop, or not in this increment.
 */
enum class Relation { EXACT, AT_LEAST, BOUNDS, UNKNOWN };

/**
 * One count: {value, relation, unit}, plus lower and upper for BOUNDS (JSON: every count of
 * §7.2). |value| is always the conservative number: the count when EXACT, the
 * lower bound when AT_LEAST or BOUNDS; it is 0 and written as JSON null when UNKNOWN.
 * BOUNDS stays BOUNDS even when lower == upper: EXACT is reserved for completed discovery.
 */
struct Count {
    Unit unit = Unit::GRAPH_CONTEXTS;
    Relation relation = Relation::UNKNOWN;
    uint64_t value = 0;
    // BOUNDS only; lower == value
    uint64_t lower = 0;
    uint64_t upper = 0;

    static Count exact(Unit unit, uint64_t value) {
        return { unit, Relation::EXACT, value, 0, 0 };
    }
    static Count at_least(Unit unit, uint64_t value) {
        return { unit, Relation::AT_LEAST, value, 0, 0 };
    }
    static Count bounds(Unit unit, uint64_t lower, uint64_t upper) {
        assert(lower <= upper);
        return { unit, Relation::BOUNDS, lower, lower, upper };
    }
    static Count unknown(Unit unit) {
        return { unit, Relation::UNKNOWN, 0, 0, 0 };
    }

    /**
     * The sum of two disjoint portions of one search, with the weakest relation (§3):
     *  UNKNOWN  + UNKNOWN            = UNKNOWN
     *  UNKNOWN  + anything else      = AT_LEAST (the explored portion's lower bound: an
     *                                  unexplored portion has no upper bound)
     *  AT_LEAST + anything           = AT_LEAST (sum of the lower bounds)
     *  BOUNDS   + BOUNDS or EXACT    = BOUNDS   (sums of lower and of upper)
     *  EXACT    + EXACT              = EXACT
     * Both operands must have the same unit.
     */
    Count& operator+=(const Count &other) {
        assert(unit == other.unit);
        auto low = [](const Count &c) -> uint64_t {
            return c.relation == Relation::UNKNOWN ? 0 : c.value;
        };
        auto high = [](const Count &c) -> uint64_t {
            return c.relation == Relation::BOUNDS ? c.upper : c.value;
        };
        if (relation == Relation::UNKNOWN && other.relation == Relation::UNKNOWN)
            return *this;
        if (relation == Relation::UNKNOWN || other.relation == Relation::UNKNOWN
                || relation == Relation::AT_LEAST || other.relation == Relation::AT_LEAST) {
            *this = at_least(unit, low(*this) + low(other));
        } else if (relation == Relation::BOUNDS || other.relation == Relation::BOUNDS) {
            *this = bounds(unit, low(*this) + low(other), high(*this) + high(other));
        } else {
            value += other.value;
        }
        return *this;
    }
};


// ---------------------------------------------------------------- patterns

// The kinds of this increment; `protein` (the codon automaton, §6) is a later increment
enum class PatternKind { DNA, IUPAC };

/**
 * The bases allowed at one pattern position: bit 0 A, bit 1 C, bit 2 G, bit 3 T (the IUPAC
 * convention: R = A|G). A set over {A, C, G, T} on purpose: a pattern's N is {A, C, G, T}
 * and never matches a record's N symbol on a DNA5 build (§3), and a set over these four
 * bases cannot express that symbol at all. The engine maps each base to its BOSS code with
 * the graph's encode(); the unconstrained flanks of `any_offset` (§4.1) are not pattern
 * positions and do admit N.
 */
using BaseSet = uint8_t;
constexpr BaseSet kBaseA = 1;
constexpr BaseSet kBaseC = 2;
constexpr BaseSet kBaseG = 4;
constexpr BaseSet kBaseT = 8;
constexpr BaseSet kAllBases = kBaseA | kBaseC | kBaseG | kBaseT;

/**
 * A refusal of one pattern: written as the `error` of that pattern's slot, the other patterns
 * of the request are still answered (§7.2). |code| is the JSON error code of the slot.
 */
class PatternError : public std::invalid_argument {
  public:
    PatternError(std::string code, const std::string &message)
          : std::invalid_argument(message), code_(std::move(code)) {}
    const std::string& code() const { return code_; }

  private:
    std::string code_;
};

/**
 * One oriented pattern: L positions, each a BaseSet (§3).
 */
class Pattern {
  public:
    /**
     * Parses |text| (case-insensitive) as |kind|: DNA over A, C, G, T; IUPAC over the 15
     * codes A C G T R Y S W K M B D H V N. Throws PatternError with code "bad_alphabet" on
     * an empty text or on any other character (U, '-', '.', whitespace included), naming
     * the first offending 0-based position. No length cap: a pattern longer than k costs
     * only its anchor window (§4.1).
     */
    static Pattern parse(PatternKind kind, std::string_view text);

    PatternKind kind() const { return kind_; }
    // the pattern in upper case; for a reverse complement, the IUPAC complement reversed
    const std::string& text() const { return text_; }
    size_t length() const { return positions_.size(); }
    const std::vector<BaseSet>& positions() const { return positions_; }

    /**
     * The engine's one primitive (§3): the bases allowed at |position| after the bases
     * |spelled| at positions [0, position) of this oriented pattern. For DNA and IUPAC it is
     * positions()[position] whatever was spelled; a peptide's codon automaton (§6, a later
     * increment) depends on |spelled|, which is why the range DFS asks through this call.
     */
    BaseSet allowed(size_t position, std::string_view spelled) const;

    // sum of log2(4 / |set_i|) over the positions [begin, end) (§3)
    double information_bits(size_t begin, size_t end) const;
    double information_bits() const { return information_bits(0, length()); }

    // every position is a single base, whatever the kind (the floor is waived for such a
    // pattern in `suffix` scope, §5.3: one range, a few ranks)
    bool is_exact() const;

    /**
     * rc(P): positions reversed, each set complemented (A<->T, C<->G, R<->Y, K<->M, B<->V,
     * D<->H; S, W, N to themselves, COMPL_TAB of reverse_complement.hpp); the kind is kept.
     * Searched completely on its own and never converted to P's ids (§4.1, "Orientation").
     */
    Pattern reverse_complement() const;

    // P == rc(P) position by position (ACGT, RY, NN): searched once, its contexts counted
    // once in every total, strand "=" per context and key "both" in by_strand (§3, "Strand")
    bool is_palindromic() const;

  private:
    Pattern(PatternKind kind, std::string text, std::vector<BaseSet> positions)
          : kind_(kind), text_(std::move(text)), positions_(std::move(positions)) {}

    PatternKind kind_;
    std::string text_;
    std::vector<BaseSet> positions_;
};


// ---------------------------------------------------------------- request and budget

/**
 * What a request answers (§5.2). In this increment the two retrieval modes exist only with
 * output.labels "none" (the route refuses every other projection): their retrieval is the
 * label-free extraction of enumerate(), and no annotation row is ever read.
 *  COUNT         discover and count; nothing retained beyond the DFS frontier; no results.
 *  ALL_OR_COUNT  (the JSON default) the count, then every context — released only when
 *                discovery completed with an EXACT total <= max_contexts — or none, withheld
 *                with the reason.
 *  PARTIAL       the count, then the first max_contexts contexts in answer order among those
 *                discovered, every cut stated (an agent asks for this explicitly).
 */
enum class Mode { COUNT, ALL_OR_COUNT, PARTIAL };

/**
 * What a request counts (§3, "Scope").
 *  SUFFIX      (L <= k) the k-mers whose last L symbols instantiate the pattern: every
 *              occurrence starting at position >= k - L of its retained island, once per
 *              distinct k-mer. Answers absence_scope "suffix_only": a label may still carry
 *              the pattern inside the first k - L bases of an island. BASIC and native
 *              CANONICAL graphs only: on a wrapped PRIMARY graph a virtual suffix is a stored
 *              prefix (§4.1), refused per pattern with "scope_unsupported".
 *  ANY_OFFSET  (L <= k, the default) every k-mer containing the pattern at any offset
 *              p in [0, k - L]: every occurrence inside a retained k-mer. Branching work even
 *              for an exact pattern (the flank ranges), charged to max_steps.
 *  LONG        (L > k) implied by the length whatever the request named, never requested:
 *              this increment counts the anchors (the k-mers instantiating positions
 *              [0, k)), extends nothing and extracts nothing.
 */
enum class Scope { SUFFIX, ANY_OFFSET, LONG };

/**
 * Which oriented patterns a request searches. FORWARD is P as given, REVERSE is rc(P), BOTH
 * both; a palindromic P is searched once whatever this says. On a BASIC graph FORWARD is the
 * + strand and REVERSE the - strand; elsewhere they name orientations of the pattern only.
 */
enum class Strands { BOTH, FORWARD, REVERSE };

/**
 * The per-pattern semantics of a request. The route fills it from the request JSON with the
 * server's defaults and caps applied (§7.1); the steps and the deadline are
 * request-wide and live in Budget.
 */
struct Request {
    // read by enumerate() only; count() answers as COUNT whatever this says
    Mode mode = Mode::ALL_OR_COUNT;
    // SUFFIX or ANY_OFFSET (LONG is never requested; a pattern with L > k is LONG)
    Scope scope = Scope::ANY_OFFSET;
    Strands strands = Strands::BOTH;
    // when set, a pattern's discovery stops as soon as its running lower bound, over all its
    // offsets and orientations, exceeds its threshold (§5.2): contexts against max_contexts
    // (L <= k), anchors against max_anchors (L > k); the count is then AT_LEAST and the stop
    // names the threshold. Ends only that pattern; the next one starts afresh
    bool stop_at_threshold = false;
    /**
     * Per pattern, all offsets and orientations (§5.3): the stop_at_threshold threshold for
     * L <= k, ALL_OR_COUNT's retrieval threshold (released only when the EXACT total is
     * <= it) and PARTIAL's cap on the contexts returned.
     */
    uint64_t max_contexts = 10'000;
    // the stop_at_threshold threshold for L > k (it admits an extension in a later increment)
    uint64_t max_anchors = 1'000;
    /**
     * L > k, enumerate() only: release the anchors (the k-mers instantiating the anchor
     * window, at offset 0) exactly as contexts are released for L <= k, with max_anchors in
     * the place of max_contexts (ALL_OR_COUNT's threshold, PARTIAL's cap). An anchor is not
     * a context of the pattern (its path may not complete, §4.2), so the contract of this
     * increment does not publish them and the route leaves this false: the results of a
     * long pattern are then withheld (PATHS_LATER_INCREMENT). Not a JSON field.
     */
    bool release_anchors = false;
    /**
     * The information floor (§5.3): it gates discovery for every pattern except an exact one
     * (Pattern::is_exact) in SUFFIX scope — IUPAC patterns, ANY_OFFSET, and LONG on the anchor
     * window's bits; a pattern below it is refused with "information_below_floor" and costs
     * nothing. Server policy (--pattern-min-information-bits), not a request field.
     */
    double min_information_bits = 24;
};

/**
 * The request's one deadline (§5.3), started when the request is parsed. Work stops at
 * start + time_budget - finalize_reserve (work_expired), so that the counts kept up to date
 * during the search can be written in the reserve; the answer must be written by
 * start + time_budget (respond_expired), else the route answers 503 "deadline" and sends
 * nothing partial as if whole.
 */
class Deadline {
  public:
    using Clock = std::chrono::steady_clock;

    // |clock| is injectable so that tests stop at a chosen instant
    Deadline(Clock::time_point start, double time_budget_ms, double finalize_reserve_ms,
             std::function<Clock::time_point()> clock = &Clock::now);

    // no deadline (tests that do not test time)
    static Deadline unbounded();

    bool work_expired() const;
    bool respond_expired() const;
    double elapsed_ms() const;
    double time_budget_ms() const { return time_budget_ms_; }
    double finalize_reserve_ms() const { return finalize_reserve_ms_; }

  private:
    Clock::time_point start_;
    double time_budget_ms_;
    double finalize_reserve_ms_;
    bool unbounded_ = false;
    std::function<Clock::time_point()> clock_;
};

// Why a pattern's work stopped, or why PARTIAL's list was cut (JSON stop.reason, cut.reason):
// the knob an agent can turn
enum class StopReason { MAX_STEPS, TIME, MAX_CONTEXTS, MAX_ANCHORS };

/**
 * Where it stopped (JSON stop.phase):
 *  DISCOVERY   the range DFS over node ranges, the W rule, the flank ranges;
 *  MASK_SCAN   the scans deferred until every range is discovered: a W-rule range's invalid
 *              edges or candidates (§4.1), and on an even-k wrapped PRIMARY graph the check
 *              of which contexts are palindromic k-mers — the only stop that, after complete
 *              discovery, can leave BOUNDS; every other discovery stop leaves AT_LEAST;
 *  EXTRACTION  enumerate()'s release of contexts, stopped by the deadline after discovery
 *              completed or stopped (its counts keep their relation: extraction counts
 *              nothing).
 */
enum class StopPhase { DISCOVERY, MASK_SCAN, EXTRACTION };

struct Stop {
    StopPhase phase;
    StopReason reason;
};

/**
 * The request's work budget: max_steps and the deadline, shared by all patterns of the
 * request in request order (§5.3: one share per shard; this increment has one shard), so
 * that the same request on the same index stops at the same step (§5.5). A step is one
 * range evaluation (a tighten_range, successful or not, or a W-rule rank set) or one edge
 * examined by a scan; extraction charges no steps (its work is linear in the contexts it
 * releases, which max_contexts caps, plus the invalid candidates it skips) but reads the
 * clock.
 * Stops are sticky: once a charge is refused, every later charge is refused with the same
 * reason, and the patterns after it answer UNKNOWN counts with that stop.
 */
class Budget {
  public:
    // the clock is read at least every this many steps or examined edges (§5.3), and at
    // every check_time()
    static constexpr uint64_t kClockStride = 4096;

    Budget(uint64_t max_steps, Deadline deadline);

    /**
     * Admits |n| more steps: true when steps_used() + n <= max_steps and the deadline's
     * work time has not passed at the last clock reading (taken when the charge crosses a
     * multiple of kClockStride). False records the reason (unless a stop was recorded
     * before) — MAX_STEPS when the steps would exceed max_steps (checked first: it does not
     * depend on the machine), else TIME — and does not count the steps; the caller must not
     * do the work.
     */
    bool charge(uint64_t n = 1);
    // reads the clock now (phase and pattern boundaries, extraction): false when the work
    // time passed, recording TIME unless a stop was recorded before
    bool check_time();

    std::optional<StopReason> stopped() const { return stopped_; }
    uint64_t steps_used() const { return steps_used_; }
    uint64_t max_steps() const { return max_steps_; }
    const Deadline& deadline() const { return deadline_; }

  private:
    uint64_t max_steps_;
    uint64_t steps_used_ = 0;
    Deadline deadline_;
    std::optional<StopReason> stopped_;
};


// ---------------------------------------------------------------- graphs

/**
 * The three graph modes the engine serves (§4.1, "Three graph modes"), told apart by
 * get_mode() and the wrapper. JSON names as the traversal's regimes.
 *  BASIC      DBGSuccinct, mode BASIC: k-mers as deposited; strand stated; every scope.
 *  CANONICAL  DBGSuccinct, mode CANONICAL (both orientations stored): every scope; no
 *             strand (note "strand_unknown_canonical").
 *  PRIMARY    CanonicalDBG wrapping a PRIMARY DBGSuccinct (as the server loads it). Contexts
 *             are the wrapper graph's: the wrapper k-mers containing Q — stored k-mers
 *             containing Q (base search of Q), and virtual rc(y) of stored y containing rc(Q)
 *             at offset k - L - p (base search of rc(Q), mapped to wrapper id y + offset).
 *             The two are united by (wrapper node id, offset) before anything counts them: a
 *             palindromic stored k-mer (even k only) is found by both and counted once. As on
 *             a native CANONICAL graph, FORWARD and REVERSE then count alike. ANY_OFFSET and
 *             LONG only; no strand.
 */
enum class GraphMode { BASIC, CANONICAL, PRIMARY };

/**
 * What the engine can do on a graph (the graph half of the capabilities' pattern block; the
 * route adds placement, support and annotation from the annotation). Answered without
 * building a PatternSearch.
 */
struct GraphSupport {
    // whether count() and enumerate() can answer on this graph at all
    bool supported = false;
    // when not (also the 400 code of the route): "representation_unsupported" (not a
    // DBGSuccinct, nor a CanonicalDBG over a PRIMARY one), "primary_unwrapped" (a PRIMARY
    // DBGSuccinct not wrapped in CanonicalDBG), "alphabet_unsupported" (the BOSS alphabet is
    // not "$ACGT" or "$ACGTN"), "mask_required" (no valid-edge mask: without it every edge,
    // dummies included, would count as a k-mer, §4)
    std::string reason;
    GraphMode mode = GraphMode::BASIC;
    bool mask_present = false;
    size_t k = 0;
    // the BOSS alphabet, sentinel first: "$ACGT" (DNA4) or "$ACGTN" (DNA5)
    std::string alphabet;
    // true on BASIC only: FORWARD/REVERSE are the + and - strands
    bool strand_stated = false;
    // the requestable scopes: SUFFIX and ANY_OFFSET, or ANY_OFFSET only (PRIMARY)
    std::vector<Scope> scopes;
};


// ---------------------------------------------------------------- results

/**
 * Which oriented pattern a count or a context belongs to: FORWARD = P, REVERSE = rc(P),
 * PALINDROMIC = P searched once because P == rc(P). JSON (§3, "Strand"): on a BASIC graph
 * strand "+", "-", "=" per context and keys "+", "-", "both" in by_strand; on the other modes,
 * where no strand is known, "forward", "reverse", "palindromic" (orientation, by_orientation).
 */
enum class Orientation { FORWARD, REVERSE, PALINDROMIC };

/**
 * counts.contexts of a pattern with L <= k (unit GRAPH_CONTEXTS). A graph context is
 * (orientation, k-mer, offset) (§3): a k-mer containing the pattern twice gives two
 * contexts; a k-mer present in many records is one context.
 */
struct ContextCounts {
    Count total = Count::unknown(Unit::GRAPH_CONTEXTS);
    // the contexts at offset k - L (equal to total in SUFFIX scope)
    Count suffix = Count::unknown(Unit::GRAPH_CONTEXTS);
    // one entry per offset of the scope (every p in [0, k - L] for ANY_OFFSET, only k - L
    // for SUFFIX), zero counts included, so that a missing key never stands for zero
    std::map<uint32_t, Count> by_offset;
    // one entry per orientation searched
    std::map<Orientation, Count> by_orientation;
};

/**
 * counts.anchors and counts.paths of a pattern with L > k. An anchor is a k-mer
 * instantiating positions [0, k) of the oriented pattern (§4.2); on a wrapped PRIMARY graph
 * the anchors of Q[0, k) and of rc(Q[0, k)), mapped, united by wrapper node id (§4.1).
 */
struct AnchorCounts {
    Count total = Count::unknown(Unit::ANCHORS);
    std::map<Orientation, Count> by_orientation;
    // completed paths (§4.2) are a later increment: UNKNOWN, except EXACT 0 when total is
    // EXACT 0 (no path starts without an anchor; a derivation, not a promotion)
    Count paths = Count::unknown(Unit::PATHS);
};

// JSON work: what the search spent, charged against Budget
struct Work {
    // range evaluations (tighten_range calls and W-rule rank sets), one step each
    uint64_t ranges_visited = 0;
    // ranges whose scan began (§4.1: zero in the common case; MASK_SCAN of StopPhase)
    uint64_t mask_scans = 0;
    // every step this pattern charged: ranges_visited plus the edges its scans examined
    uint64_t steps = 0;
};

// A per-pattern refusal decided by the engine (JSON: the slot's error {code, message});
// bad_alphabet comes earlier, from Pattern::parse
struct Refusal {
    // "information_below_floor" | "scope_unsupported"
    std::string code;
    std::string message;
};

// Why enumerate() publishes nothing (JSON withheld.reason, §5.2)
enum class Withheld {
    // ALL_OR_COUNT: discovery completed, EXACT total > max_contexts
    COUNT_ABOVE_THRESHOLD,
    // ALL_OR_COUNT: stop_at_threshold stopped discovery (stop reason MAX_CONTEXTS/MAX_ANCHORS)
    THRESHOLD_CROSSED,
    // ALL_OR_COUNT: discovery or a mask scan stopped at max_steps
    DISCOVERY_BUDGET,
    // ALL_OR_COUNT: discovery or extraction stopped by the deadline (§5.3)
    DEADLINE,
    // either mode, L > k with anchors not EXACT 0: results of a long pattern are its
    // completed paths (§4.2), a later increment; anchors are never released as results
    PATHS_LATER_INCREMENT,
};

/**
 * What enumerate() released: the label-free retrieval (§4.3 "The label-free path first",
 * §5.2 "Projection none"), single shard. Absent after count().
 */
struct Extraction {
    // contexts passed to the callback
    uint64_t returned = 0;
    /**
     * Every context of the pattern in its scope and strands was released: discovery completed
     * with an EXACT count and returned equals it (an EXACT 0 included, and for L > k anchors
     * EXACT 0, hence paths EXACT 0). JSON retrieval_complete: the one flag that licenses an
     * absence claim over graph contexts (§5.1); it claims nothing about labels, which this
     * increment never reads.
     */
    bool complete = false;
    // set: nothing is published ("results": []) and the callback received nothing (an
    // ALL_OR_COUNT release is buffered, so a deadline inside it withholds all of it)
    std::optional<Withheld> withheld;
    // PARTIAL, neither complete nor withheld: why the list is shorter than the pattern's
    // contexts — the reason of the stop that touched the pattern when there is one, else
    // MAX_CONTEXTS (the first max_contexts in answer order were returned; MAX_ANCHORS for
    // the anchors of Request::release_anchors)
    std::optional<StopReason> cut;
};

/**
 * One graph context of a pattern with L <= k, as enumerate() releases it (the route writes
 * one JSON result from it: kmer = graph.get_node_sequence(node) on the graph given to
 * PatternSearch, instance = kmer.substr(offset, L), which instantiates P for FORWARD and
 * PALINDROMIC and rc(P) for REVERSE, row = AnnotatedDBG::graph_to_anno_index(base_node)).
 * A context is (orientation, node, offset) (§3): an IUPAC pattern and its reverse
 * complement can both match one instance (NA and TN both match TA), giving two contexts
 * at one (node, offset) that differ in orientation; an exact DNA pattern cannot.
 */
struct Context {
    Orientation orientation;
    // the oriented pattern's 0-based offset inside the k-mer of |node|, in [0, k - L]
    uint32_t offset;
    // in the graph given to PatternSearch: the BOSS edge index on a DBGSuccinct, the wrapper
    // id on a wrapped PRIMARY graph (a stored k-mer's id, or that plus the wrapper's offset
    // for a virtual reverse complement); never converted between orientations
    DeBruijnGraph::node_index node;
    // the stored k-mer that carries the annotation row: |node| itself, except for a virtual
    // node on a wrapped PRIMARY graph, where it is CanonicalDBG::get_base_node(node)
    DeBruijnGraph::node_index base_node;
};

// JSON notes of a pattern (§7.2), the ones this increment can state:
//  low_complexity_pattern    an exact pattern that sdust flags with the seeder's parameters
//                            (T = 20, W = 64, is_low_complexity): why its counts are large
//  strand_unknown_canonical  graph mode CANONICAL or PRIMARY: orientations, not strands
//  paths_later_increment     L > k: anchors counted, paths neither extended nor extracted
constexpr const char kNoteLowComplexity[] = "low_complexity_pattern";
constexpr const char kNoteStrandUnknown[] = "strand_unknown_canonical";
constexpr const char kNotePathsLater[] = "paths_later_increment";

/**
 * The answer for one pattern. With |refusal| set nothing was searched and only the pattern
 * description (scope, palindromic, information bits) is meaningful.
 */
struct Result {
    std::optional<Refusal> refusal;

    // SUFFIX | ANY_OFFSET as requested, LONG when L > k
    Scope scope = Scope::ANY_OFFSET;
    GraphMode graph_mode = GraphMode::BASIC;
    bool palindromic = false;
    // the orientations searched, in search order (FORWARD before REVERSE)
    std::vector<Orientation> searched;
    double information_bits = 0;
    // L > k: the bits of the anchor window [0, k), stated separately (§3)
    std::optional<double> anchor_information_bits;

    // exactly one of the two is set on an answered pattern: contexts when L <= k, anchors
    // when L > k
    std::optional<ContextCounts> contexts;
    std::optional<AnchorCounts> anchors;

    Work work;
    // the first stop that touched this pattern (discovery, mask scan or extraction), or none
    // when its work completed; a pattern reached after the budget stopped carries that stop
    // (phase DISCOVERY) with UNKNOWN counts
    std::optional<Stop> stop;
    // a TIME stop touched it, in any phase: the answer depends on the machine (JSON
    // determinism "time_limited", else "full"; §5.5)
    bool time_limited = false;
    // set by enumerate() on an answered pattern; never by count()
    std::optional<Extraction> extraction;
    // kNote* values, in the order the constants are declared
    std::vector<std::string> notes;
    double elapsed_ms = 0;
};


// ---------------------------------------------------------------- the engine

/**
 * Counts and extracts the graph contexts of a pattern on the succinct graph by BOSS range
 * narrowing (§4.1), never reading annotation:
 *  - positions 0 .. m-2 (m = min(L, k)) on node ranges, a DFS from the whole edge range
 *    trying only Pattern::allowed() at each depth (the seeder's suffix_to_prefix loop with
 *    its symbol set made a parameter);
 *  - position m-1 on W: per allowed symbol c, the valid edges of the leaf range with
 *    W in {c, c + alph_size} — plain and marked ranks, minus the invalid edges among them,
 *    found by a scan through the mask that is charged one step per edge;
 *  - for ANY_OFFSET, offsets p with p + L <= k - 1 by continuing the DFS from the pattern's
 *    leaves with every non-sentinel symbol (N included on DNA5) for k - 1 - p - L more
 *    steps and counting every valid edge leaving those nodes.
 * No dummy edge is relied on and nothing is scored.
 *
 * Stateless between calls; one instance may serve concurrent requests, each with its own
 * Budget.
 */
class PatternSearch {
  public:
    // What the engine can do on |graph| (see GraphSupport); never throws
    static GraphSupport support(const DeBruijnGraph &graph);

    // Throws std::invalid_argument naming support(graph).reason when it is not supported.
    // |graph| must outlive this object.
    explicit PatternSearch(const DeBruijnGraph &graph);

    const GraphSupport& graph_support() const { return support_; }

    /**
     * Counts |pattern| under |request| (as mode COUNT), charging |budget|. Range descriptors
     * are discarded as they are counted, except those whose scan is deferred until
     * discovery completes (at most one per step charged, none in the common case); the
     * rest of the memory is the DFS frontier, O(k * alphabet). Refusals (information floor, SUFFIX on a wrapped PRIMARY graph) come
     * back as Result::refusal without charging anything. On a budget already stopped (or
     * whose work time has passed at the pattern's check_time) every count is UNKNOWN and
     * stop is {DISCOVERY, the budget's reason}. Never throws for a parsed pattern.
     */
    Result count(const Pattern &pattern, const Request &request, Budget &budget) const;

    /**
     * The label-free path (§4.3): discovery exactly as count() — the same steps charged, the
     * same counts, work, stop and refusals — then the release of request.mode (ALL_OR_COUNT
     * or PARTIAL; COUNT throws std::invalid_argument), calling |callback| once per released
     * context, in answer order: node ascending, then offset ascending, then orientation
     * (FORWARD, REVERSE, PALINDROMIC) (§5.5, "BOSS edge order, then by offset"; the
     * orientation separates the two contexts an IUPAC pattern can have at one offset).
     *  ALL_OR_COUNT  releases every context iff discovery completed with an EXACT total
     *                <= max_contexts; otherwise none (Extraction::withheld). The release is
     *                buffered: the callback receives all of them, or nothing.
     *  PARTIAL       releases the first max_contexts contexts in answer order among those
     *                discovered, also after a MAX_STEPS or threshold stop (what was built is
     *                delivered, the cut stated), streamed to the callback; none after a TIME
     *                stop in discovery; after a TIME stop in the release, what was called.
     * L > k: nothing is released (PATHS_LATER_INCREMENT) unless anchors are EXACT 0, or
     * request.release_anchors asks for the anchors (Context offset 0).
     * Release runs only while the deadline's work time has not passed, reading the clock every
     * kClockStride edges examined; a stop there is {EXTRACTION, TIME}. Every released context
     * is a valid edge whose k-mer contains the oriented pattern at its offset.
     * Memory: discovery retains at most one leaf range per step charged (ALL_OR_COUNT may drop
     * them once the running lower bound exceeds max_contexts), freed before returning; the
     * contexts themselves are the caller's.
     */
    Result enumerate(const Pattern &pattern, const Request &request, Budget &budget,
                     const std::function<void(const Context&)> &callback) const;

  private:
    const DeBruijnGraph &graph_;
    GraphSupport support_;
    // the succinct graph searched: |graph_| itself, or the PRIMARY graph a CanonicalDBG wraps
    const DBGSuccinct *dbg_succ_ = nullptr;
    // wrapped PRIMARY: the wrapper's id of the reverse complement of stored node y is
    // y + wrapper_offset_ (y itself when y is palindromic); 0 otherwise
    uint64_t wrapper_offset_ = 0;

    // count() and enumerate(): the release runs when |callback| is set
    Result run(const Pattern &pattern, const Request &request, Budget &budget,
               const std::function<void(const Context&)> *callback) const;
};


// ---------------------------------------------------------------- JSON names

inline const char* to_string(Unit unit) {
    switch (unit) {
        case Unit::GRAPH_CONTEXTS: return "graph_contexts";
        case Unit::ANCHORS: return "anchors";
        case Unit::PATHS: return "paths";
        case Unit::PLACED_OCCURRENCES: return "placed_occurrences";
        case Unit::LABELS: return "labels";
    }
    return "unknown";
}

inline const char* to_string(Relation relation) {
    switch (relation) {
        case Relation::EXACT: return "exact";
        case Relation::AT_LEAST: return "at_least";
        case Relation::BOUNDS: return "bounds";
        case Relation::UNKNOWN: return "unknown";
    }
    return "unknown";
}

inline const char* to_string(PatternKind kind) {
    switch (kind) {
        case PatternKind::DNA: return "dna";
        case PatternKind::IUPAC: return "iupac";
    }
    return "unknown";
}

// the request's and the answer's `mode`
inline const char* to_string(Mode mode) {
    switch (mode) {
        case Mode::COUNT: return "count";
        case Mode::ALL_OR_COUNT: return "all_or_count";
        case Mode::PARTIAL: return "partial";
    }
    return "unknown";
}

// the request's and the answer's `scope`
inline const char* to_string(Scope scope) {
    switch (scope) {
        case Scope::SUFFIX: return "suffix";
        case Scope::ANY_OFFSET: return "any_offset";
        case Scope::LONG: return "long";
    }
    return "unknown";
}

// the answer's `absence_scope` (§5.1): over what a zero count would claim absence
inline const char* absence_scope(Scope scope) {
    switch (scope) {
        case Scope::SUFFIX: return "suffix_only";
        case Scope::ANY_OFFSET: return "any_offset";
        case Scope::LONG: return "long";
    }
    return "unknown";
}

// the request's `strands`
inline const char* to_string(Strands strands) {
    switch (strands) {
        case Strands::BOTH: return "both";
        case Strands::FORWARD: return "forward";
        case Strands::REVERSE: return "reverse";
    }
    return "unknown";
}

inline const char* to_string(GraphMode mode) {
    switch (mode) {
        case GraphMode::BASIC: return "basic";
        case GraphMode::CANONICAL: return "canonical";
        case GraphMode::PRIMARY: return "primary";
    }
    return "unknown";
}

// a context's `strand` and the answer's `strands` list (BASIC graphs)
inline const char* strand_symbol(Orientation orientation) {
    switch (orientation) {
        case Orientation::FORWARD: return "+";
        case Orientation::REVERSE: return "-";
        case Orientation::PALINDROMIC: return "=";
    }
    return "unknown";
}

// the key in by_strand (BASIC graphs): a palindromic pattern's contexts under "both" (§3)
inline const char* strand_key(Orientation orientation) {
    switch (orientation) {
        case Orientation::FORWARD: return "+";
        case Orientation::REVERSE: return "-";
        case Orientation::PALINDROMIC: return "both";
    }
    return "unknown";
}

// a context's `orientation`, the answer's `strands` list and the key in by_orientation
// (CANONICAL and PRIMARY graphs)
inline const char* orientation_key(Orientation orientation) {
    switch (orientation) {
        case Orientation::FORWARD: return "forward";
        case Orientation::REVERSE: return "reverse";
        case Orientation::PALINDROMIC: return "palindromic";
    }
    return "unknown";
}

inline const char* to_string(StopReason reason) {
    switch (reason) {
        case StopReason::MAX_STEPS: return "max_steps";
        case StopReason::TIME: return "time";
        case StopReason::MAX_CONTEXTS: return "max_contexts";
        case StopReason::MAX_ANCHORS: return "max_anchors";
    }
    return "unknown";
}

inline const char* to_string(StopPhase phase) {
    switch (phase) {
        case StopPhase::DISCOVERY: return "discovery";
        case StopPhase::MASK_SCAN: return "mask_scan";
        case StopPhase::EXTRACTION: return "extraction";
    }
    return "unknown";
}

inline const char* to_string(Withheld withheld) {
    switch (withheld) {
        case Withheld::COUNT_ABOVE_THRESHOLD: return "count_above_threshold";
        case Withheld::THRESHOLD_CROSSED: return "threshold_crossed";
        case Withheld::DISCOVERY_BUDGET: return "discovery_budget";
        case Withheld::DEADLINE: return "deadline";
        case Withheld::PATHS_LATER_INCREMENT: return "paths_later_increment";
    }
    return "unknown";
}

// the request's `mode`, `scope` (LONG is not requestable) and `strands`; nullopt: not a
// valid value. A valid mode is still refused by the route when its projection is a later
// increment's (§7.1).
inline std::optional<Mode> parse_mode(std::string_view s) {
    if (s == "count") return Mode::COUNT;
    if (s == "all_or_count") return Mode::ALL_OR_COUNT;
    if (s == "partial") return Mode::PARTIAL;
    return std::nullopt;
}

inline std::optional<Scope> parse_scope(std::string_view s) {
    if (s == "suffix") return Scope::SUFFIX;
    if (s == "any_offset") return Scope::ANY_OFFSET;
    return std::nullopt;
}

inline std::optional<Strands> parse_strands(std::string_view s) {
    if (s == "both") return Strands::BOTH;
    if (s == "forward") return Strands::FORWARD;
    if (s == "reverse") return Strands::REVERSE;
    return std::nullopt;
}

} // namespace pattern
} // namespace graph
} // namespace mtg

#endif // __PATTERN_SEARCH_HPP__
