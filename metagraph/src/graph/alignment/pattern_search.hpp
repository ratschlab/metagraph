#ifndef __PATTERN_SEARCH_HPP__
#define __PATTERN_SEARCH_HPP__

/**
 * Pattern search: count, and extract without reading any annotation, every graph context of
 * a short motif or an IUPAC pattern (peptides later) on the succinct graph, and for a pattern
 * longer than k every graph path spelling it (phase 2, the extension).
 *
 * The design is docs/DESIGN-pattern-search.md (v6 + §22); section numbers below refer to it.
 * This header is the contract between the engine (pattern_search.cpp) and the route
 * (src/cli/pattern.cpp) for increments 0-2 and the engine half of increment 4: modes count,
 * all_or_count and partial, the last two only with output.labels "none" (the label-free path,
 * §4.3); no annotation read, no placement, no predicate, single graph. The JSON the route
 * writes from these types is the contract of §7 (src/cli/pattern.cpp writes it); every JSON
 * name comes from the to_string() family at the end of this file.
 *
 * The extension beyond k (§4.2, increment 4) is opt-in: Request::extend_paths. Without it
 * (the default) a pattern longer than k is answered exactly as in increments 1-2 (anchors
 * counted, paths UNKNOWN, nothing released, note paths_later_increment), so the route's
 * answers do not change until it opts in.
 *
 * Wiring the route for `long` (a follow-up; the route is not edited by increment 4's engine
 * work). The paths are OPT-IN per request (owner decision #13 of 2026-10-07; SPEC §12), an
 * addition under contract version 1: a request without the option keeps today's anchor-only
 * answer. In src/cli/pattern.cpp:
 *  1. parse_request: accept the request field "long_search" (reserved: refused by name as
 *     later_increment today, whatever its value): "anchors" (the default) leaves
 *     extend_paths false; "paths" sets req.request.extend_paths = true and admits
 *     "max_paths" (refused by name today; capped like max_anchors, by a new
 *     --pattern-max-paths) into req.request.max_paths;
 *  2. entry_json, the `anchors` branch, with long_search "paths": counts.paths =
 *     count_json(anchors->paths) as today, plus paths["candidates_examined"] =
 *     anchors->candidates_examined, the per-strand split put_orientations(&paths,
 *     anchors->paths_by_orientation, strand_stated), and paths["extension"] =
 *     to_string(anchors->extension); work["extension_edges"] = result->work.extension_edges;
 *     timing["extension_ms"] = result->extension_ms;
 *  3. context_json: a released path (!c.path.empty()) is written with NEW fields, never
 *     `kmer`, which keeps its meaning (the k-mer of a graph context): sequence = c.sequence
 *     (the L spelled bases, §7.2), anchor_kmer = graph.get_node_sequence(c.node) (the k-base
 *     anchor), instance = c.sequence, offset 0, and with output.paths the node ids c.path
 *     (rows: AnnotatedDBG::graph_to_anno_index(search.base_node(n)) per node); the check
 *     `c.offset + length > k` applies to L <= k contexts only;
 *  4. capabilities: advertise the option (long_search values "anchors", "paths"; max_paths
 *     among the caps); long_patterns keeps describing the default answer. For a request
 *     with long_search "paths", L > k, count mode included: counts.paths becomes known, the
 *     note and the withheld reason paths_later_increment do not appear (never produced with
 *     extend_paths), and anchors_above_threshold, stop phase "extension" and reason
 *     "max_paths" can appear. A request without it keeps counts.paths unknown and
 *     paths_later_increment, as today.
 * Everything else (the withheld reasons, the cut, retrieval_complete, the stop and its phase)
 * flows through the existing Extraction and Stop fields with the values added below
 * (Withheld::ANCHORS_ABOVE_THRESHOLD, StopReason::MAX_PATHS, StopPhase::EXTENSION).
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
 *            deferred scans were interrupted: lower <= true <= upper.
 *  UNKNOWN   the phase never ran: a stop came before it started, or not in this increment.
 * A search whose discovery was entered is AT_LEAST even when the stop refused its very first
 * step (AT_LEAST 0): it was interrupted, not skipped (SPEC §7.4). After any discovery stop
 * every per-offset count of the interrupted search is AT_LEAST too: v1 does not track which
 * offsets were completely discovered before the stop.
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
     * the first offending 0-based position. No length cap: a pattern longer than k is
     * charged steps for its anchor windows only (§4.1); parsing it, its information bits,
     * its palindrome test and its low-complexity note cost O(L) time that no step charges
     * and no clock reading interrupts (each a single pass, a few ns per base).
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
     * The extension (§4.2) asks it at every position >= k with |spelled| the whole instance
     * so far: the anchor's k spelled bases (read from the graph) followed by the bases the
     * DFS appended, so that an automaton recomputes its state at the k boundary from the
     * anchor's sequence (§4.1, last bullet) and needs no state carried by the engine.
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
 *              the anchors (the k-mers instantiating positions [0, k)) are counted, and with
 *              Request::extend_paths every admitted anchor is extended along the pattern to
 *              position L (§4.2): the complete paths are the pattern's contexts (offset 0,
 *              n = L - k + 1 k-mers, every one of them retained: §3 "Covered sequence").
 *              Without extend_paths nothing is extended or extracted (increments 1-2).
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
    // (L <= k), anchors against max_anchors (L > k), and with extend_paths the extension
    // stops as soon as more than max_paths paths are complete; the count of that phase is
    // then AT_LEAST, every later phase's UNKNOWN, and the stop names the threshold. Ends only
    // that pattern; the next one starts afresh.
    // Checked in discovery (and the extension) only: the deferred scans (MASK_SCAN) never
    // consult it and run to their end or to a budget stop, so a pattern can answer EXACT
    // above its threshold with no stop (withheld COUNT_ABOVE_THRESHOLD, or cut MAX_CONTEXTS
    // in PARTIAL), or BOUNDS with a budget stop in the scans. The running lower bound lags the
    // count by the discount of masked edges not yet scanned and, on an even-k wrapped
    // PRIMARY graph, by the palindromic k-mers both base searches may find (subtracted as
    // the candidates of the ranges where a palindrome is possible): there the stop can fire
    // late or not at all, and the answer is then the EXACT count with no stop.
    bool stop_at_threshold = false;
    /**
     * Per pattern, all offsets and orientations (§5.3): the stop_at_threshold threshold for
     * L <= k, ALL_OR_COUNT's retrieval threshold (released only when the EXACT total is
     * <= it) and PARTIAL's cap on the contexts returned.
     */
    uint64_t max_contexts = 10'000;
    /**
     * L > k: the stop_at_threshold threshold of the anchors, and with extend_paths the
     * extension's admission (§4.2, §5.2) in every mode, count() included: the anchors are
     * extended only when their count is EXACT and <= max_anchors; above it the paths stay
     * UNKNOWN (Extension::NOT_ADMITTED) and enumerate() withholds ANCHORS_ABOVE_THRESHOLD.
     */
    uint64_t max_anchors = 1'000;
    /**
     * L > k: run phase 2 (§4.2), the depth-first extension of every admitted anchor along
     * outgoing edges to position L, trying only Pattern::allowed() at each position. The
     * complete paths are then the pattern's contexts: AnchorCounts::paths is counted (EXACT
     * when every anchor was extended), enumerate() releases paths (Context::path and
     * Context::sequence) instead of withholding PATHS_LATER_INCREMENT, and the note
     * paths_later_increment is not set. False (the default): increments 1-2, unchanged —
     * anchors counted, paths UNKNOWN (EXACT 0 without anchors), nothing extended or
     * released for L > k, no extension step charged. Not a JSON field: the route will set it
     * for a request with long_search "paths" (opt-in, owner decision #13; see "Wiring the
     * route" at the top of this file); every other request keeps it false.
     */
    bool extend_paths = false;
    /**
     * L > k with extend_paths: the retrieval threshold on completed paths (§4.2, §5.2) —
     * ALL_OR_COUNT releases the paths only when their count is EXACT and <= max_paths,
     * PARTIAL releases the first max_paths paths in answer order, and with
     * stop_at_threshold the extension stops as soon as more than max_paths are complete
     * (stop {EXTENSION, MAX_PATHS}). Retention (§5.2, "Descriptor retention by phase"): a
     * count() keeps no path; ALL_OR_COUNT keeps the paths while their count is <= max_paths
     * and drops them all once it is exceeded; PARTIAL keeps the first max_paths.
     */
    uint64_t max_paths = 1'000;
    /**
     * L > k without extend_paths, enumerate() only: release the anchors (the k-mers
     * instantiating the anchor window, at offset 0) exactly as contexts are released for
     * L <= k, with max_anchors in the place of max_contexts (ALL_OR_COUNT's threshold,
     * PARTIAL's cap). An anchor is not a context of the pattern (its path may not complete,
     * §4.2), so the contract does not publish them and the route leaves this false. Not a
     * JSON field; a test and diagnostic aid.
     * It is not what keeps the anchors for the extension: with extend_paths the anchors
     * are retained in every mode, count() included, through the extension's admission —
     * dropped as soon as their running count exceeds max_anchors (the admission fails) and
     * once every anchor has been listed for the extension (§5.2, "Descriptor retention by
     * phase") — and they are never released as results. enumerate() refuses
     * release_anchors together with extend_paths (std::invalid_argument).
     */
    bool release_anchors = false;
    /**
     * The information floor (§5.3): it gates discovery for every pattern except an exact one
     * (Pattern::is_exact) in SUFFIX scope — IUPAC patterns, ANY_OFFSET, and LONG on the bits
     * of the anchor window of EVERY searched orientation (FORWARD: P[0, k); REVERSE: rc(P)[0,
     * k), the reverse complement of P's last k positions; a palindromic P has one), so that
     * no searched window below it runs; a pattern below it is refused with
     * "information_below_floor" and costs nothing. Server policy
     * (--pattern-min-information-bits), not a request field.
     */
    double min_information_bits = 24;
};

/**
 * Thrown by a Budget whose abort predicate (Budget::set_abort) answers true at a clock reading:
 * the request's caller has left, or the server is stopping. Nothing of the request is
 * answered (the route closes the connection); the engine's state is freed by unwinding.
 */
class Aborted : public std::runtime_error {
  public:
    Aborted() : std::runtime_error("pattern: the request was aborted") {}
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
// the knob an agent can turn. MAX_PATHS: the stop_at_threshold threshold of the extension,
// and PARTIAL's cap on the paths returned
enum class StopReason { MAX_STEPS, TIME, MAX_CONTEXTS, MAX_ANCHORS, MAX_PATHS };

/**
 * Where it stopped (JSON stop.phase):
 *  DISCOVERY   the range DFS over node ranges, the W rule, the flank ranges;
 *  MASK_SCAN   the scans deferred until every range is discovered: a W-rule range's invalid
 *              edges or candidates (§4.1), and on an even-k wrapped PRIMARY graph the check
 *              of which contexts are palindromic k-mers — the only stop that, after complete
 *              discovery, can leave BOUNDS; every other discovery stop leaves AT_LEAST;
 *  EXTRACTION  enumerate()'s release of contexts or paths, stopped by the deadline after
 *              discovery (and extension) completed or stopped (its counts keep their
 *              relation: extraction counts nothing);
 *  EXTENSION   L > k with Request::extend_paths: the listing of the admitted anchors and the
 *              DFS beyond k (§4.2), stopped by max_steps, the deadline, or max_paths with
 *              stop_at_threshold. The anchors stay EXACT, the paths become AT_LEAST.
 */
enum class StopPhase { DISCOVERY, MASK_SCAN, EXTRACTION, EXTENSION };

struct Stop {
    StopPhase phase;
    StopReason reason;
};

/**
 * The request's work budget: max_steps and the deadline, shared by all patterns of the
 * request in request order (§5.3: one share per shard; this increment has one shard), so
 * that the same request on the same index stops at the same step (§5.5). A step is one
 * range evaluation (a tighten_range, successful or not, or a W-rule rank set), one edge
 * examined by a scan, or one outgoing edge examined by the extension (§4.2: every edge
 * DeBruijnGraph::call_outgoing_kmers reports for a node the DFS expands, allowed or not);
 * extraction charges no steps (its work is linear in the contexts it releases, which
 * max_contexts caps, plus the invalid candidates it skips) but reads the clock.
 * Stops are sticky: once a charge is refused, every later charge is refused with the same
 * reason, and the patterns after it answer UNKNOWN counts with that stop.
 *
 * Where the engine reads the clock: at the start of every pattern, at every charge that
 * crosses a multiple of kClockStride steps (discovery, the deferred scans and the extension
 * share one step counter, so the boundary from discovery to the scans has no reading of its
 * own: at most kClockStride - 1 steps pass unread there, as anywhere), before every release,
 * every kClockStride edges the release examines and every kReleaseClockStride contexts it
 * passes on (a released context costs the caller a k-mer spelling, k - 1 BOSS steps, and its
 * result object), and every kReleaseClockStride contexts an ALL_OR_COUNT delivery passes to
 * the caller. Work therefore ends within one such stride after the work time passes.
 */
class Budget {
  public:
    // the clock is read at least every this many steps or examined edges (§5.3), and at
    // every check_time()
    static constexpr uint64_t kClockStride = 4096;
    // ... and every this many contexts released to a caller, whose work per context (a
    // spelling and a result object, microseconds) is far above a step's
    static constexpr uint64_t kReleaseClockStride = 64;

    Budget(uint64_t max_steps, Deadline deadline);

    /**
     * |aborted| is asked at every clock reading (the stride crossings of charge(), every
     * check_time(), the release's and the delivery's readings); when it answers true the
     * reading throws Aborted, and nothing of the request is answered. The route passes
     * "the client left, or the server stops". Unset (the default): never asked.
     */
    void set_abort(std::function<bool()> aborted) { aborted_ = std::move(aborted); }

    /**
     * Admits |n| more steps: true when steps_used() + n <= max_steps and the deadline's
     * work time has not passed at the last clock reading (taken when the charge crosses a
     * multiple of kClockStride). False records the reason (unless a stop was recorded
     * before) — MAX_STEPS when the steps would exceed max_steps (checked first: it does not
     * depend on the machine), else TIME — and does not count the steps; the caller must not
     * do the work.
     */
    bool charge(uint64_t n = 1);
    // reads the clock now (the start of a pattern, the release and its strides, the
    // delivery's strides): false when the work time passed, recording TIME unless a stop was
    // recorded before
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
    std::function<bool()> aborted_;

    void poll_abort() const {
        if (aborted_ && aborted_())
            throw Aborted();
    }
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
 *             palindromic stored k-mer (over A, C, G, T: even k only) is found by both and
 *             counted once. (On a $ACGTN build an odd k-mer with N at its centre can equal its
 *             reverse complement; CanonicalDBG serves such a stored k-mer at two ids, and the
 *             engine, which unites palindromes for even k only, then counts and releases it
 *             twice; a DNA5 build of the engine's tests confirms it. DNA5 is not in CI.) As on
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
    // dummies included, would count as a k-mer, §4). The route narrows it further
    // (cli::route_support): "alphabet_untested" ($ACGTN, not served until a DNA5 build passes
    // the pattern tests) and "mask_invalid" (a mask marking a W = $ edge valid)
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
 * What phase 2 (§4.2) did for a pattern with L > k (JSON, once wired: counts.paths.extension).
 *  NOT_REQUESTED  Request::extend_paths is false: increments 1-2, nothing extended; paths
 *                 UNKNOWN, EXACT 0 when the anchors are EXACT 0.
 *  NO_ANCHORS     the anchors are EXACT 0: nothing to extend; paths EXACT 0 (no path starts
 *                 without an anchor: a derivation, not a promotion).
 *  NOT_STARTED    the anchors are not EXACT (a step, time or anchor-threshold stop in
 *                 discovery or a mask scan, see Result::stop): an anchor set not known
 *                 completely is never extended; paths UNKNOWN (§3: every phase after a stop).
 *  NOT_ADMITTED   the anchors are EXACT and above max_anchors: the extension's admission
 *                 failed (§4.2); paths UNKNOWN; enumerate() withholds ANCHORS_ABOVE_THRESHOLD.
 *  STOPPED        the extension started and stopped (Result::stop, phase EXTENSION: max_steps,
 *                 the deadline, or max_paths with stop_at_threshold): paths AT_LEAST, the
 *                 paths completed before the stop.
 *  COMPLETED      every anchor was extended to L: paths EXACT.
 */
enum class Extension { NOT_REQUESTED, NO_ANCHORS, NOT_STARTED, NOT_ADMITTED, STOPPED, COMPLETED };

/**
 * counts.anchors and counts.paths of a pattern with L > k. An anchor is a k-mer
 * instantiating positions [0, k) of the oriented pattern (§4.2); on a wrapped PRIMARY graph
 * the anchors of Q[0, k) and of rc(Q[0, k)), mapped, united by wrapper node id (§4.1). A
 * path is a walk of n = L - k + 1 k-mers along outgoing edges of the served graph, from an
 * anchor of the oriented pattern Q, spelling an instance of Q; its identity is
 * (orientation, node path) (§3, "Graph context", offset 0). Paths exist only where all n
 * k-mers were retained (§3, "Covered sequence"), and a path is a graph context, not a
 * record occurrence: two records ACG and CGT make the path ACGT (§4.3; per-label support
 * is retrieval's).
 */
struct AnchorCounts {
    Count total = Count::unknown(Unit::ANCHORS);
    std::map<Orientation, Count> by_orientation;
    /**
     * Completed paths (unit PATHS), by Extension: EXACT when COMPLETED (and EXACT 0 with
     * NO_ANCHORS, or without extend_paths when the anchors are EXACT 0), AT_LEAST when
     * STOPPED, UNKNOWN otherwise. Never derived from the anchors (an anchor may have no
     * path: AAAC with AAA present).
     */
    Count paths = Count::unknown(Unit::PATHS);
    /**
     * With extend_paths, one entry per orientation searched (empty without it): EXACT when
     * every anchor of that orientation was extended (EXACT 0 when its anchors are EXACT 0),
     * AT_LEAST when the extension stopped before, UNKNOWN when it did not run. |paths| is
     * their sum with the weakest relation, except that a total whose extension did not run
     * is UNKNOWN (EXACT 0 with NO_ANCHORS) rather than AT_LEAST over its known zeros.
     */
    std::map<Orientation, Count> paths_by_orientation;
    /**
     * The branches the DFS entered (§4.2: candidates_examined beside the paths): every
     * partial path of k + 1 .. L bases the extension formed by appending an allowed
     * outgoing k-mer, complete paths included, anchors not. Work done, not a count of the
     * pattern: exact as such whether or not the extension completed; 0 when it did not run.
     * The outgoing edges examined, allowed or not, are Work::extension_edges.
     */
    uint64_t candidates_examined = 0;
    Extension extension = Extension::NOT_REQUESTED;
};

// JSON work: what the search spent, charged against Budget
struct Work {
    // range evaluations (tighten_range calls and W-rule rank sets), one step each
    uint64_t ranges_visited = 0;
    // ranges whose deferred scan began (MASK_SCAN of StopPhase): a W-rule range's masked
    // edges, zero on BASIC, CANONICAL and odd-k graphs unless masked edges sit among the
    // candidates; on an even-k wrapped PRIMARY graph also the palindrome check of every range
    // at an offset where a palindromic k-mer can hold the pattern (one get_node_sequence, k - 1
    // BOSS steps, per context, one step each): there about one per such range
    uint64_t mask_scans = 0;
    // L > k with extend_paths: the outgoing edges the extension examined, one step each
    // (allowed or not); 0 otherwise
    uint64_t extension_edges = 0;
    // every step this pattern charged: ranges_visited, plus the edges its scans examined,
    // plus extension_edges
    uint64_t steps = 0;
    // diagnostic, not a JSON field of contract version 1: the most range descriptors held at
    // once for the release (24 bytes each; see PatternSearch::enumerate, "Memory")
    uint64_t spans_retained_peak = 0;
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
    // ALL_OR_COUNT: discovery (and for L > k the extension) completed, EXACT total >
    // max_contexts (L <= k), or EXACT paths > max_paths (L > k with extend_paths)
    COUNT_ABOVE_THRESHOLD,
    // ALL_OR_COUNT: stop_at_threshold stopped discovery or the extension (stop reason
    // MAX_CONTEXTS, MAX_ANCHORS or MAX_PATHS)
    THRESHOLD_CROSSED,
    // ALL_OR_COUNT: discovery, a mask scan or the extension stopped at max_steps
    DISCOVERY_BUDGET,
    // ALL_OR_COUNT: discovery, the extension or extraction stopped by the deadline (§5.3)
    DEADLINE,
    // either mode, L > k without extend_paths, anchors not EXACT 0: results of a long
    // pattern are its completed paths (§4.2), not extended in increments 1-2; anchors are
    // never released as results
    PATHS_LATER_INCREMENT,
    // either mode, L > k with extend_paths: the anchors are EXACT and above max_anchors, so
    // the extension was not admitted (§4.2); counts.anchors is exact, the paths UNKNOWN
    ANCHORS_ABOVE_THRESHOLD,
};

/**
 * What enumerate() released: the label-free retrieval (§4.3 "The label-free path first",
 * §5.2 "Projection none"), single shard. Absent after count().
 */
struct Extraction {
    // contexts (L > k with extend_paths: paths) passed to the callback
    uint64_t returned = 0;
    /**
     * Every context of the pattern in its scope and strands was released: discovery completed
     * with an EXACT count and returned equals it (an EXACT 0 included, and for L > k anchors
     * EXACT 0, hence paths EXACT 0; with extend_paths, paths EXACT and all of them
     * returned). JSON retrieval_complete: the one flag that licenses an absence claim over
     * graph contexts (§5.1); it claims nothing about labels, which this increment never reads.
     */
    bool complete = false;
    // set: nothing is published ("results": [], returned 0). The callback received nothing,
    // except when the deadline stopped an ALL_OR_COUNT delivery (DEADLINE, stop {EXTRACTION,
    // TIME}): it then received a prefix of the contexts, which the caller must discard
    std::optional<Withheld> withheld;
    // PARTIAL, neither complete nor withheld: why the list is shorter than the pattern's
    // contexts — the reason of the stop that touched the pattern when there is one, else
    // MAX_CONTEXTS (the first max_contexts in answer order were returned; MAX_ANCHORS for
    // the anchors of Request::release_anchors; MAX_PATHS for paths)
    std::optional<StopReason> cut;
};

/**
 * One graph context, as enumerate() releases it.
 *
 * L <= k (the route writes one JSON result from it: kmer = graph.get_node_sequence(node) on
 * the graph given to PatternSearch, instance = kmer.substr(offset, L), which instantiates P
 * for FORWARD and PALINDROMIC and rc(P) for REVERSE, row =
 * AnnotatedDBG::graph_to_anno_index(base_node)). A context is (orientation, node, offset)
 * (§3): an IUPAC pattern and its reverse complement can both match one instance (NA and TN
 * both match TA), giving two contexts at one (node, offset) that differ in orientation; an
 * exact DNA pattern cannot. |path| and |sequence| are empty.
 *
 * L > k with Request::extend_paths: one path (§4.2), (orientation, path) its identity:
 * offset 0, |node| and |base_node| its anchor (the first k-mer), |path| its n = L - k + 1
 * nodes and |sequence| its L spelled bases, which instantiate P for FORWARD and PALINDROMIC
 * and rc(P) for REVERSE (§7.2, owner decision #13: the result's new field sequence is
 * |sequence|, its new field anchor_kmer the anchor's k-mer; kmer is not used for a path).
 * The anchors of Request::release_anchors are released like L <= k contexts (offset 0,
 * |path| and |sequence| empty).
 */
struct Context {
    Orientation orientation;
    // the oriented pattern's 0-based offset inside the k-mer of |node|, in [0, k - L]; 0 for
    // a path
    uint32_t offset;
    // in the graph given to PatternSearch: the BOSS edge index on a DBGSuccinct, the wrapper
    // id on a wrapped PRIMARY graph (a stored k-mer's id, or that plus the wrapper's offset
    // for a virtual reverse complement); never converted between orientations
    DeBruijnGraph::node_index node;
    // the stored k-mer that carries the annotation row: |node| itself, except for a virtual
    // node on a wrapped PRIMARY graph, where it is CanonicalDBG::get_base_node(node)
    // (PatternSearch::base_node maps any node of |path| the same way)
    DeBruijnGraph::node_index base_node;
    // a path only: its k-mers in reading order, ids as |node|; path.front() == node
    std::vector<DeBruijnGraph::node_index> path = {};
    // a path only: the L bases it spells (its first k are the anchor's k-mer)
    std::string sequence = {};
};

// JSON notes of a pattern (§7.2), the ones this increment can state:
//  low_complexity_pattern    an exact pattern that sdust flags with the seeder's parameters
//                            (T = 20, W = 64, is_low_complexity): why its counts are large
//  strand_unknown_canonical  graph mode CANONICAL or PRIMARY: orientations, not strands
//  paths_later_increment     L > k without Request::extend_paths: anchors counted, paths
//                            neither extended nor extracted (never set with extend_paths)
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
    // the orientations searched, in plan order (FORWARD before REVERSE). The base searches
    // behind them run in the order of their estimated cost, the cheaper first (the plan's
    // order on a tie, as for every exact pattern): a budget stop then leaves the cheap one
    // complete and the expensive one interrupted, rather than the reverse
    std::vector<Orientation> searched;
    double information_bits = 0;
    // L > k: the bits of the anchor window P[0, k), stated separately (§3; its meaning in
    // contract version 1, whatever the strands searched)
    std::optional<double> anchor_information_bits;
    // L > k: the least over the searched orientations' anchor windows (FORWARD and
    // PALINDROMIC: P[0, k); REVERSE: rc(P)[0, k), whose bits are those of P[L - k, L)), the
    // bits the information floor gated on (X-GUARANTEES-01; an addition to the JSON of
    // contract version 1, the owner's decision of 2026-10-07)
    std::optional<double> min_anchor_information_bits;
    // L > k: the bits of each searched orientation's anchor window (an addition to the JSON
    // of contract version 1 when the route publishes it)
    std::map<Orientation, double> anchor_window_bits;

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
    // L > k with extend_paths: the time spent in phase 2 (the listing of the anchors and the
    // DFS; §7.2 timing.extension_ms), part of elapsed_ms; 0 when it did not run
    double extension_ms = 0;
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
 *    steps and counting every valid edge leaving those nodes;
 *  - for L > k with Request::extend_paths (phase 2, §4.2): the anchors above (positions
 *    [0, k) on node ranges and W), listed in answer order once they are counted EXACT and
 *    admitted (<= max_anchors), then extended one by one by a depth-first search along
 *    DeBruijnGraph::call_outgoing_kmers on the graph given to PatternSearch (the wrapper on
 *    a wrapped PRIMARY graph, so a path may pass between stored and virtual k-mers), keeping
 *    only the outgoing k-mers whose last base is in Pattern::allowed(position, spelled), to
 *    position L. Every complete path is a context; the search never scores and stops at the
 *    first disallowed base, so it is complete within its budget. Outgoing edges come from
 *    the valid-edge mask's graph (DBGSuccinct::call_outgoing_kmers skips dummy and pruned
 *    k-mers), so a path exists only where all its k-mers were retained.
 * No dummy edge is relied on and nothing is scored.
 *
 * Cost: a searched window (the oriented pattern, or a long pattern's anchor window) is
 * matched from its first position on the widest ranges, so a run of pattern N there would be
 * branched four ways per position before any specified base narrows a range. On a $ACGT
 * graph, where every base of a valid k-mer is one of A, C, G, T, such a leading run of length
 * r is not searched: the window at offset p is exactly its core at offset p + r, so the core
 * is searched with its offsets shifted (counts, offsets and the release unchanged, only the
 * work). On $ACGTN (pattern N never matches the graph's N) the run is searched as given. An
 * N run inside a window still costs about min(4^run, edges / 4^a) ranges per level, a being
 * the specified bases before it in that orientation: the information bits do not bound it.
 * The base searches run cheapest first by that estimate, so that a budget stop leaves the
 * orientation whose run comes late complete.
 *
 * Answer order (§5.5): contexts by (node, offset, orientation); paths by (anchor node,
 * orientation), then the DFS in symbol order A < C < G < T at every position, i.e. the
 * paths of one anchor by their spelled sequence. Deterministic on the same index.
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
     * The stored k-mer carrying the annotation row of |node| of the graph given to
     * PatternSearch (Context::base_node, for every node of a Context::path): |node| itself,
     * except for a virtual reverse complement on a wrapped PRIMARY graph, where it is
     * CanonicalDBG::get_base_node(node).
     */
    DeBruijnGraph::node_index base_node(DeBruijnGraph::node_index node) const {
        return support_.mode == GraphMode::PRIMARY && node > wrapper_offset_
            ? node - wrapper_offset_
            : node;
    }

    /**
     * Counts |pattern| under |request| (as mode COUNT), charging |budget|. Range descriptors
     * are discarded as they are counted, except those whose scan is deferred until
     * discovery completes (88 bytes each, at most one per step charged): none in the common
     * case, but on an even-k wrapped PRIMARY graph one per range at an offset where a
     * palindromic k-mer can hold the pattern, which for a short or degenerate pattern is
     * nearly every range, so up to max_steps of them; the rest of the memory is the DFS
     * frontier, O(k * alphabet). Refusals (information floor, SUFFIX on a wrapped PRIMARY graph) come
     * back as Result::refusal without charging anything. On a budget already stopped (or
     * whose work time has passed at the pattern's check_time) every count is UNKNOWN and
     * stop is {DISCOVERY, the budget's reason}. Never throws for a parsed pattern, except
     * std::logic_error on a broken internal invariant (never expected: an answer that
     * cannot be stated correctly is not stated at all).
     * L > k with request.extend_paths: the anchors are retained (§5.2: kept through the
     * extension's admission, dropped once their count exceeds max_anchors) and, when
     * admitted, extended to count the paths (AnchorCounts::paths, Extension); no path is
     * retained. The extension's memory is the anchor list (<= max_anchors) and the DFS
     * stack, O(n * alphabet) for n = L - k + 1.
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
     *                buffered (a released set that differs from the EXACT count throws
     *                std::logic_error, nothing being published), then delivered to the
     *                callback under the clock (read every kReleaseClockStride contexts): the
     *                callback receives all of them, or — when the work time passes during
     *                the delivery — a prefix, after which the result says withheld DEADLINE
     *                with returned 0 and stop {EXTRACTION, TIME}, and the caller must discard
     *                what it received.
     *  PARTIAL       releases the first max_contexts contexts in answer order among those
     *                discovered, also after a MAX_STEPS or threshold stop (what was built is
     *                delivered, the cut stated), streamed to the callback; none after a TIME
     *                stop in discovery; after a TIME stop in the release, what was called.
     *                On an even-k wrapped PRIMARY graph a palindromic stored k-mer that only
     *                the reverse-complement search reached before a stop is released too, so
     *                that the list holds every context the counts credit.
     * L > k without request.extend_paths: nothing is released (PATHS_LATER_INCREMENT) unless
     * anchors are EXACT 0, or request.release_anchors asks for the anchors (Context offset 0).
     * L > k with request.extend_paths (release_anchors must then be false, else
     * std::invalid_argument): the discovery and the extension exactly as count() — the same
     * steps, counts and stop — then the paths (Context::path, Context::sequence) in answer
     * order (anchor node, orientation, then the spelled sequence), which is the DFS's own
     * order, so the paths kept during the DFS are released as they are:
     *  ALL_OR_COUNT  every path iff the extension completed (paths EXACT) with paths
     *                <= max_paths, else none: withheld ANCHORS_ABOVE_THRESHOLD (anchors EXACT
     *                > max_anchors), COUNT_ABOVE_THRESHOLD (paths EXACT > max_paths),
     *                THRESHOLD_CROSSED (a MAX_ANCHORS or MAX_PATHS stop), DISCOVERY_BUDGET
     *                (MAX_STEPS in any phase), DEADLINE (TIME in any phase).
     *  PARTIAL       the first max_paths paths in answer order among those completed, also
     *                after a MAX_STEPS or MAX_PATHS stop in the extension (cut = the stop's
     *                reason, else MAX_PATHS when more paths exist); none after a stop in
     *                discovery (an anchor set not known completely is never extended; cut =
     *                that reason), none after a TIME stop (cut TIME), none when the
     *                extension was not admitted (withheld ANCHORS_ABOVE_THRESHOLD).
     *                The release reads the clock once before the first callback.
     *  An EXACT 0 of anchors or of paths is a complete, empty release.
     * Release runs only while the deadline's work time has not passed, reading the clock
     * before it starts, every kClockStride descriptors it prepares, every kClockStride edges
     * it examines and every kReleaseClockStride contexts it passes on (the caller's work per
     * context runs between two readings); a stop there is {EXTRACTION, TIME}. Every released
     * context is a valid edge whose k-mer contains the oriented pattern at its offset; every
     * released path is a walk of valid k-mers whose sequence instantiates the oriented
     * pattern.
     * Memory: discovery retains a 24-byte descriptor per range with contexts, freed before
     * returning. ALL_OR_COUNT keeps them while the running lower bound is <= max_contexts
     * and drops them all once it is above. PARTIAL keeps only those that can hold one of the
     * first max_contexts contexts in answer order: once the descriptors whose node range
     * ends at or before a node T are sure to hold max_contexts contexts, every descriptor
     * whose range starts after T is dropped, at once or at the next compaction (when the
     * retained ones exceed max(4096, twice those kept at the last one)); with max_contexts 0
     * none is kept. What PARTIAL keeps is therefore about max_contexts descriptors plus the
     * ranges that straddle T (O(k) per base search) plus the W-rule ranges whose masked edges
     * leave no context guaranteed, at most 4096 or twice that between compactions, whatever
     * the discovery's size; the release's own index over them (40 bytes each) is as large.
     * The paths retained for the release are at most max_paths, each O(L); the contexts
     * themselves are the caller's.
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
        case StopReason::MAX_PATHS: return "max_paths";
    }
    return "unknown";
}

inline const char* to_string(StopPhase phase) {
    switch (phase) {
        case StopPhase::DISCOVERY: return "discovery";
        case StopPhase::MASK_SCAN: return "mask_scan";
        case StopPhase::EXTRACTION: return "extraction";
        case StopPhase::EXTENSION: return "extension";
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
        case Withheld::ANCHORS_ABOVE_THRESHOLD: return "anchors_above_threshold";
    }
    return "unknown";
}

// counts.paths.extension (once the route is wired): what phase 2 did
inline const char* to_string(Extension extension) {
    switch (extension) {
        case Extension::NOT_REQUESTED: return "not_requested";
        case Extension::NO_ANCHORS: return "no_anchors";
        case Extension::NOT_STARTED: return "not_started";
        case Extension::NOT_ADMITTED: return "not_admitted";
        case Extension::STOPPED: return "stopped";
        case Extension::COMPLETED: return "completed";
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
