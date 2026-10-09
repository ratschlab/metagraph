#ifndef __PATTERN_SEARCH_HPP__
#define __PATTERN_SEARCH_HPP__

/**
 * Pattern search: count, and extract without reading any annotation, every graph context of
 * a short motif, an IUPAC pattern or a peptide (its codon automaton, §6) on the succinct
 * graph, and for a pattern longer than k every graph path spelling it (phase 2, the
 * extension).
 *
 * The design is docs/DESIGN-pattern-search.md; section numbers below refer to it.
 * This header is the contract between the engine (pattern_search.cpp) and the route
 * (src/cli/pattern.cpp): modes count, all_or_count and partial, whose release is label-free
 * (§4.3: the engine reads no annotation, places nothing and evaluates no predicate), single
 * graph. The JSON the route writes from these types is the contract of §7 (src/cli/pattern.cpp
 * writes it); every JSON name comes from the to_string() family at the end of this file.
 *
 * The extension beyond k (§4.2) is opt-in: Request::extend_paths. Without it (the default) a
 * pattern longer than k is answered by its anchors only (anchors counted, paths UNKNOWN,
 * nothing released, note paths_later_increment).
 *
 * How the route serves them (src/cli/pattern.cpp; SPEC-pattern-search.md §12):
 *  - `long`: OPT-IN per request, an addition under contract version 1: the request field
 *    long_search "paths" sets Request::extend_paths and admits max_paths (--pattern-max-paths);
 *    a request without it keeps the anchor-only answer (counts.paths unknown, note and withheld
 *    reason paths_later_increment). With it, counts.paths states the extension's count, its
 *    per-orientation split, candidates_examined and the Extension; a released path is written
 *    with the fields sequence (the L bases), anchor_kmer, nodes and rows, never `kmer`, which
 *    keeps its meaning (a context's k-mer); the withheld reason anchors_above_threshold, stop
 *    phase "extension" and reason "max_paths" can appear. Labels of paths are read by
 *    PatternRetrieval::retrieve_paths (support label_intersection or record_verified).
 *  - `protein` (§6): patterns[i].protein is the third kind; the request field genetic_code (an
 *    integer, default GeneticCode::kStandard) is looked up with GeneticCode::find (an unknown
 *    id: 400 genetic_code_unknown), and Pattern::parse(PatternKind::PROTEIN, text, code)'s
 *    PatternError (bad_alphabet) is the slot's error as for dna and iupac. A peptide is a
 *    Pattern like any other: counts, contexts, anchors, the extension and its paths need no
 *    route code of their own; its entry states kind "protein", pattern = text() (the residues),
 *    length = length() (3m bases: every offset, scope decision, information bit and instance is
 *    in bases), residues (m) and genetic_code. The capabilities list the kind, the residues,
 *    the genetic codes (GeneticCode::ids()) and the default 1. The stop '*' is a residue (a
 *    stop codon of the table; see Pattern::parse): a peptide holding it in a table without an
 *    unconditional stop codon has no instance and is answered with the note kNoteNoStopCodon
 *    (EXACT 0 contexts without a search for L <= k).
 *  - `supported_paths` (not served yet: no route sets a tracker or a sink): the extension is a
 *    DFS over search states (SearchState) with Pattern as its Model; the route is to pass a
 *    SupportTracker and a PathSink through Request::support and Request::sink, and to read
 *    AnchorCounts::walks, supported, branches_pruned and pruned_before_completion. Without them
 *    the extension is the plain one of long_search "paths".
 *  - graphs without the dummy-edge mask: served, see "Graphs without the dummy-edge mask" at
 *    PatternSearch. Their counts that the engine could not resolve are BOUNDS [lower, U]
 *    (Count), the route adding the estimate U x f with f from sample_real_fraction(); the lists
 *    stay exact. A pattern with few unchecked candidates is checked entry by entry and counted
 *    exactly (Request::max_checked_entries; the route passes --pattern-max-checked-entries).
 *
 * For every type here nothing is weakened silently: every count carries its unit and its
 * relation, and a count is never promoted by assumption.
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

#include "graph/alignment/genetic_code.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct;

namespace pattern {

// The version of the /pattern request and answer JSON (capabilities:
// pattern.pattern_contract_version). Raised when a field changes meaning; a field or value that
// is only added does not raise it.
constexpr int kPatternContractVersion = 1;


// ---------------------------------------------------------------- counts

/**
 * What a count counts (§3). Units are never summed into each other: a graph-context count says
 * nothing about placed occurrences, an anchor count nothing about completed paths.
 * PLACED_OCCURRENCES and LABELS need annotation, which the engine never reads: it states them
 * UNKNOWN.
 */
enum class Unit { GRAPH_CONTEXTS, ANCHORS, PATHS, PLACED_OCCURRENCES, LABELS };

/**
 * How a count's value relates to the true number (§3, §4.1).
 *  EXACT     the discovery behind it completed: no step, time or threshold stop touched it.
 *  AT_LEAST  a stop interrupted discovery; the value is the sum of the lower bounds of what
 *            was explored. An undiscovered IUPAC branch, offset or strand has no upper bound,
 *            so this is the relation of every interrupted discovery.
 *  BOUNDS    every range of every branch, offset and orientation was discovered and only the
 *            deferred scans were interrupted: lower <= true <= upper. On a graph without the
 *            dummy-edge mask also (and mostly) a completed discovery whose candidates were not
 *            all checked for source dummies (more of them than Request::max_checked_entries):
 *            upper = U, the candidate entries (source dummies included), lower = the part
 *            known exactly (PatternSearch, "Graphs without the dummy-edge mask").
 *  UNKNOWN   the phase never ran: a stop came before it started, or the request does not run
 *            it.
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

/**
 * The pattern kinds (§3): DNA over A, C, G, T; IUPAC over the 15 codes; PROTEIN a peptide,
 * searched as its codon automaton (§6): 3 positions per residue, the bases allowed at a
 * position depending on the bases already spelled in its codon.
 */
enum class PatternKind { DNA, IUPAC, PROTEIN };

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
 * What may come next along an instance: the automaton the extension (§4.2) asks at every
 * position beyond the anchor window. Pattern is the only Model (DNA and IUPAC: a set per
 * position; a peptide: its codon automaton, §6). A Model is deterministic and of fixed length:
 * from a state, a base leads to at most one state, and an instance is a walk of exactly
 * length() bases.
 *
 * The extension's search state (SearchState) is the oriented node, the position (the bases
 * spelled), the model's state and the support tracker's frame: the DFS of PatternSearch asks
 * start() once per anchor, bases() once per node it expands and next() once per child it
 * enters, so that no model reads the spelled bases again (Pattern::allowed, which does, stays
 * the primitive of the range DFS of phase 1).
 *
 * Open, not designed (nothing below is built):
 *  - bounded mismatches (docs/NOTE-mismatch-search.md): a Model whose State packs (the pattern's
 *    own state, the mismatches used) and whose bases() admits every base while mismatches are
 *    left (next() adds one where the pattern does not admit the base). Deterministic, fixed
 *    length: it plugs in as is, and the trackers and sinks do not change. Its anchor phase (the
 *    BOSS ranges of phase 1) would take the same Model, the range DFS carrying a State per range
 *    in place of the CodonWindow's bases read back from the BOSS.
 *  - a profile model (docs/NOTE-biological-discovery-feature-priorities.md, "Gene family
 *    bundles": profile position, match / insertion / deletion state, coding phase): needs what
 *    this interface leaves out — several next states per base (match and insertion: next()
 *    would append to a short list, and the DFS's Level would hold one child per (base, node,
 *    state) instead of one per base), transitions that consume no base (a deletion: a next()
 *    without a base, the DFS staying on its node), accepting states at several depths (a
 *    variable-length instance: accepting() asked at every depth, the DFS completing there and
 *    going on), and a bound on the score a branch can still reach (a bound(state) the DFS
 *    prunes on). The DFS would then key its work by (node, state); the trackers and sinks keep
 *    their interfaces (a walk is still a walk).
 *  - a selective anchor: a Model that also steps left (prev(), the codon automaton read
 *    backwards), with SearchState::side LEFT along incoming edges. It also needs start(window,
 *    offset) for a window at offset w (Pattern::start derives a peptide's codon phase from the
 *    window's length alone), and a left arm's SearchState: spelled growing at the left,
 *    position = w - depth, base the node's first base (support_step's Arm::LEFT).
 *  - the entry point: only extend_anchor<Pattern> is instantiated, from
 *    PatternRun::extend_listed with the oriented Pattern, and the anchor phase does not ask the
 *    Model. Another Model needs a per-orientation `const Model*` (a Request field or a
 *    PatternRun argument) with extend_anchor<Model> instantiated through the virtual calls, and
 *    the range DFS of phase 1 taking Model::bases() per range (from Model::start of the window)
 *    in place of CodonWindow::allowed.
 */
class Model {
  public:
    // opaque to the engine: the automaton's state after the bases spelled so far
    using State = uint64_t;

    virtual ~Model() = default;

    // the instance's length in bases (L; 3m for a peptide of m residues)
    virtual size_t length() const = 0;
    // the state after the anchor window: the first |anchor_kmer|.size() bases of a walk, as
    // spelled (the anchor's k-mer, read from the graph)
    virtual State start(std::string_view anchor_kmer) const = 0;
    // the bases allowed at |position| (0-based, < length()) in state |s| (after |position|
    // bases); 0: none
    virtual BaseSet bases(State s, uint32_t position) const = 0;
    // the state after |base| (one of bases(s, position)) at |position|
    virtual State next(State s, uint32_t position, char base) const = 0;
    // a walk of |position| bases ending in state |s| is an instance; asked at length() only
    virtual bool accepting(State s, uint32_t position) const {
        (void)s;
        return position == length();
    }

  protected:
    Model() = default;
    Model(const Model&) = default;
    Model& operator=(const Model&) = default;
};

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
 * One oriented pattern: L positions, each a BaseSet (§3); for a peptide, L = 3m positions
 * whose allowed bases follow its codon automaton (§6). The Model of the extension (§4.2): its
 * State is 0 for DNA and IUPAC (a position's set does not depend on what was spelled) and for
 * a peptide the bases already spelled in the current codon (the prefix allowed() reads back
 * from |spelled|), so that bases(s, position) == allowed(position, spelled) for the state s that
 * start() and next() reach along |spelled|.
 */
class Pattern final : public Model {
  public:
    /**
     * Parses |text| (case-insensitive) as |kind|: DNA over A, C, G, T; IUPAC over the 15
     * codes A C G T R Y S W K M B D H V N; PROTEIN as parse(PROTEIN, text,
     * GeneticCode::standard()) below. Throws PatternError with code "bad_alphabet" on
     * an empty text or on any other character (U, '-', '.', whitespace included), naming
     * the first offending 0-based position. No length cap: a pattern longer than k is
     * charged steps for its anchor windows only (§4.1); parsing it, its information bits
     * and its palindrome test cost O(L) time that no step charges and no clock reading
     * interrupts (each a single pass, a few ns per base). Its low-complexity note (sdust,
     * O(L) too, but tens of microseconds per base of a repeat) is read under the clock and
     * stops at the first piece flagged (kNoteLowComplexity).
     */
    static Pattern parse(PatternKind kind, std::string_view text);

    /**
     * As above; for PROTEIN (§6) |code| is the genetic code (the route passes
     * GeneticCode::get(request genetic_code), 1 by default), for DNA and IUPAC it is not read.
     * A peptide is a string over the 20 residues A C D E F G H I K L M N P Q R S T V W Y, the
     * ambiguity codes X (any residue: every codon of |code| that is not a stop), B (D or N), Z
     * (E or Q) and J (I or L), and the stop '*' (a stop codon of |code| at that position, the
     * codons whose residue in NCBI's ncbieaa is '*', GeneticCode::stops()), case-insensitive;
     * it is the pattern of L = 3m positions whose instances are exactly the codon strings c_1
     * .. c_m with c_i a codon of residue i in |code| (§6: no superset; a stop codon only where
     * the peptide has '*', X never matching one). The codons tables 27, 28 and 31 list as a
     * residue that ends translation only in context (GeneticCode::context_stops()) are that
     * residue's and never a '*': those tables have no unconditional stop, so there '*' admits
     * no codon and the peptide has no instance (has_instances() false; the engine answers it
     * with the note kNoteNoStopCodon, see PatternSearch::count).
     * Refused with "bad_alphabet": U (selenocysteine), O (pyrrolysine), '-', digits,
     * whitespace and every other character outside these 24 letters and '*', naming the
     * first such 0-based residue position, and an empty text.
     */
    static Pattern parse(PatternKind kind, std::string_view text, const GeneticCode &code);

    PatternKind kind() const { return kind_; }
    /**
     * The pattern in upper case; for a DNA or IUPAC reverse complement, the IUPAC complement
     * reversed. A peptide's is its residues (the route's `pattern`); a peptide's reverse
     * complement, an automaton with no residue letters of its own, has its residues in
     * reverse order (a diagnostic only, never parsed back).
     */
    const std::string& text() const { return text_; }
    // L, in bases: 3m for a peptide of m residues
    size_t length() const override { return positions_.size(); }
    /**
     * Per position the bases it admits. For DNA and IUPAC exactly the pattern; for a peptide
     * the union over the residue's codons at that codon position, a superset of what the
     * automaton admits after a given prefix (allowed() is the exact primitive; the engine
     * uses these sets only where a superset is safe: its cost estimate and the test whether
     * a palindromic k-mer can hold the pattern, §4.1).
     */
    const std::vector<BaseSet>& positions() const { return positions_; }

    // PROTEIN: the genetic code's NCBI id; 0 for DNA and IUPAC
    int genetic_code() const { return genetic_code_; }
    /**
     * PROTEIN: per residue of this oriented pattern, in reading order, the codons it admits
     * (residue i covers positions [3i, 3i + 3)); for the reverse complement, the original's
     * codon sets in reverse order, each codon reverse-complemented (§6: GCN becomes NGC).
     * Empty for DNA and IUPAC.
     */
    const std::vector<CodonSet>& codon_sets() const { return codons_; }

    /**
     * The engine's one primitive (§3): the bases allowed at |position| after the bases
     * |spelled| at positions [0, position) of this oriented pattern. For DNA and IUPAC it is
     * positions()[position] whatever was spelled. For a peptide it is the codon automaton
     * (§6): with residue i = position / 3 and j = position % 3, exactly the bases b such
     * that some codon of residue i agrees with the last j bases of |spelled| (that codon's
     * first j) and has b at its position j; 0 when none does (a prefix off the automaton,
     * an N included). Only the last j <= 2 bases of |spelled| are read; when |spelled| holds
     * fewer than j, the missing ones are taken as unknown (any base).
     * The extension (§4.2) asks it at every position >= k with |spelled| the whole instance
     * so far: the anchor's k spelled bases (read from the graph) followed by the bases the
     * DFS appended, so that an automaton recomputes its state at the k boundary from the
     * anchor's sequence (§4.1, last bullet) and needs no state carried by the engine. The
     * range DFS of phase 1 asks the same automaton with the bases its range's nodes end with
     * (pattern_search.cpp, CodonWindow).
     */
    BaseSet allowed(size_t position, std::string_view spelled) const;

    /**
     * The Model: allowed() with its state carried instead of read back. start() reads the
     * anchor's last (k mod 3) bases for a peptide (its codon prefix at position k; a base other
     * than A, C, G, T there admits nothing after it, as allowed() answers 0), nothing
     * otherwise; next() appends a base to the codon prefix, emptied at a codon's end. O(1) but
     * bases(), which is allowed()'s codon test for a peptide.
     */
    State start(std::string_view anchor_kmer) const override;
    BaseSet bases(State s, uint32_t position) const override;
    State next(State s, uint32_t position, char base) const override;
    // every walk of L bases the automaton admitted is an instance
    bool accepting(State, uint32_t position) const override { return position == length(); }

    /**
     * The information of the positions [begin, end) (§3): DNA and IUPAC, the sum of
     * log2(4 / |set_i|); PROTEIN, 2 (end - begin) - log2(the number of distinct strings the
     * automaton admits over those positions), computed per residue as the strings its codons
     * spell over the positions of [begin, end) it covers (exact also for a window that cuts
     * a codon: the anchor window [0, k), its reverse P[L - k, L)), never a per-position sum
     * (which overstates a residue whose codons share no position-wise structure). For a
     * whole peptide: the sum of log2(64 / |codons_i|) = 6m - log2(the number of codon
     * strings it admits). A residue admitting no codon (a '*' in a table without a stop
     * codon, has_instances() false) is counted as one exact codon, 2 bits per position it
     * covers, so that the value stays finite (a window holding it admits no string, and its
     * search ends at that position).
     */
    double information_bits(size_t begin, size_t end) const;
    double information_bits() const { return information_bits(0, length()); }

    /**
     * False iff the pattern has no instance at all: a peptide holding a '*' read in a genetic
     * code without an unconditional stop codon (tables 27, 28, 31). Every DNA and IUPAC pattern
     * has instances.
     */
    bool has_instances() const { return has_instances_; }

    // every position is a single base, whatever the kind (the floor is waived for such a
    // pattern in `suffix` scope, §5.3: one range, a few ranks); a peptide is exact when
    // every residue has one codon in its code (M and W in the standard code)
    bool is_exact() const;

    /**
     * rc(P): positions reversed, each set complemented (A<->T, C<->G, R<->Y, K<->M, B<->V,
     * D<->H; S, W, N to themselves, COMPL_TAB of reverse_complement.hpp); the kind is kept.
     * A peptide's is its reverse-complemented automaton (§6): residues in reverse order,
     * each codon reverse-complemented. Searched completely on its own and never converted to
     * P's ids (§4.1, "Orientation").
     */
    Pattern reverse_complement() const;

    /**
     * P == rc(P) (ACGT, RY, NN): searched once, its contexts counted once in every total,
     * strand "=" per context and key "both" in by_strand (§3, "Strand"). DNA and IUPAC
     * position by position; a peptide residue by residue (its instances are a product of
     * codon sets, equal to rc's iff every residue's codons are the reverse complements of
     * those of its mirror residue): only X runs in the tables without a stop codon (27, 28,
     * 31) are palindromic peptides (and there a peptide whose '*' mirror one another, which
     * has no instance: has_instances()). In a table with stop codons the stops are never the
     * reverse complements of a residue's codons, nor of themselves (PatternPeptide,
     * StopIsNeverAMirror).
     */
    bool is_palindromic() const;

  private:
    Pattern(PatternKind kind, std::string text, std::vector<BaseSet> positions)
          : kind_(kind), text_(std::move(text)), positions_(std::move(positions)) {}

    PatternKind kind_;
    std::string text_;
    std::vector<BaseSet> positions_;
    // PROTEIN only
    std::vector<CodonSet> codons_;
    int genetic_code_ = 0;
    bool has_instances_ = true;
};


// ---------------------------------------------------------------- request and budget

/**
 * What a request answers (§5.2). The engine's retrieval in the two retrieval modes is the
 * label-free extraction of enumerate(): it reads no annotation row (the route's
 * PatternRetrieval reads the labels).
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
 *              Without extend_paths nothing is extended or extracted.
 */
enum class Scope { SUFFIX, ANY_OFFSET, LONG };

/**
 * Which oriented patterns a request searches. FORWARD is P as given, REVERSE is rc(P), BOTH
 * both; a palindromic P is searched once whatever this says. On a BASIC graph FORWARD is the
 * + strand and REVERSE the - strand; elsewhere they name orientations of the pattern only.
 */
enum class Strands { BOTH, FORWARD, REVERSE };

// the extension's support tracker and path sink, declared with the results below
class SupportTracker;
class PathSink;

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
    // On a graph without the dummy-edge mask the running count compared is the running UPPER
    // bound U (conservative: the stop can fire while the true count is within the threshold;
    // the note kNoteThresholdUpperBound says when the lower bound did not cross it), and so are
    // ALL_OR_COUNT's retention and the extension's admission. On an even-k wrapped PRIMARY
    // graph the running U leaves out the ranges whose palindrome scan is pending (the scan can
    // only lower their share), so that it never exceeds the final U: there the stop can fire
    // late or not at all, the admission after discovery then comparing the final U. (While the
    // pattern's unchecked candidates are at most max_checked_entries, the retention compares
    // the running lower bound instead, since the check after discovery may make the count EXACT
    // and within the threshold; this stop still compares U.)
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
     * paths_later_increment is not set. False (the default): anchors counted, paths UNKNOWN
     * (EXACT 0 without anchors), nothing extended or released for L > k, no extension step
     * charged. Not a JSON field: the route sets it for a request with long_search "paths"
     * (opt-in; see "How the route serves them" at the top of this file); every other request
     * keeps it false.
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
    /**
     * A graph without the dummy-edge mask only (never read with the mask): when discovery and
     * its deferred scans completed with no stop and the pattern's unchecked candidates (the
     * entries counted into U and not into the lower bound, summed over its base searches:
     * U - lower of the total on a BASIC or CANONICAL graph; on a wrapped PRIMARY graph each
     * entry enters both orientations' counts) number at most this, each is tested with
     * BOSS::node_has_sentinel and only the real k-mers are counted: every count of the pattern
     * (total, suffix, by_offset, by_orientation; the anchors of a long one) is then EXACT, 0
     * when all of them were source dummies. Each entry tested is charged k - 1 steps (the most
     * symbols the test reads), so the check costs at most max_checked_entries x (k - 1) steps;
     * a stop during it (phase MASK_SCAN) leaves the counts BOUNDS. With more unchecked
     * candidates nothing is tested (the counts stay BOUNDS). 0 tests nothing (the default here;
     * the route passes the server's --pattern-max-checked-entries, kDefaultMaxCheckedEntries
     * unless set). Server policy, not a request field.
     */
    uint64_t max_checked_entries = 0;
    /**
     * L > k with extend_paths: the supported-path search, read by the extension only (not
     * served yet: no route sets them). Not JSON fields: the route is to set them per pattern
     * for a request with long_search "supported_paths". Both null (the default): the plain
     * extension — the same steps, clock readings, counts and paths as without these fields.
     *  support  carries each branch's support beside the DFS (SupportTracker: a frame per level,
     *           pushed when a child is entered and popped with its level) and prunes a branch
     *           it declares DEAD; the paths are then the SUPPORTED complete walks
     *           (AnchorCounts::supported), every rule of max_paths (retention, ALL_OR_COUNT's
     *           threshold, PARTIAL's cap, stop_at_threshold) applies to them, and
     *           AnchorCounts::paths becomes a plain count of the complete walks (AT_LEAST once
     *           a branch was pruned before reaching L). Must outlive the call.
     *  sink     receives every path as it completes (with a tracker the supported ones), in
     *           answer order, instead of the engine's own list: the sink is the release, so
     *           count() serves it and enumerate() refuses an extending pattern with a sink
     *           (std::invalid_argument). Must outlive the call.
     */
    SupportTracker *support = nullptr;
    PathSink *sink = nullptr;
};

// the server's default of Request::max_checked_entries (--pattern-max-checked-entries: a block
// of at most ~50 entries, about 1,500 backward steps at k = 31)
constexpr uint64_t kDefaultMaxCheckedEntries = 50;

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

    bool work_expired() const;
    bool respond_expired() const;
    double elapsed_ms() const;
    double time_budget_ms() const { return time_budget_ms_; }
    double finalize_reserve_ms() const { return finalize_reserve_ms_; }

  private:
    Clock::time_point start_;
    double time_budget_ms_;
    double finalize_reserve_ms_;
    std::function<Clock::time_point()> clock_;
};

// Why a pattern's work stopped, or why PARTIAL's list was cut (JSON stop.reason, cut.reason):
// the knob an agent can turn. MAX_PATHS: the stop_at_threshold threshold of the extension, and
// PARTIAL's cap on the paths returned. EXTERNAL (internal): the extension's support tracker or
// path sink (Request::support, Request::sink) answered STOPPED; the route, which owns them,
// writes their own reason (SupportTracker::stop_reason, PathSink::stop_reason:
// max_annotation_work, max_memory), never "external". Never set without a tracker or a sink
enum class StopReason { MAX_STEPS, TIME, MAX_CONTEXTS, MAX_ANCHORS, MAX_PATHS, EXTERNAL };

/**
 * Where it stopped (JSON stop.phase):
 *  DISCOVERY   the range DFS over node ranges, the W rule, the flank ranges;
 *  MASK_SCAN   the scans deferred until every range is discovered: a W-rule range's invalid
 *              edges or candidates (§4.1), on an even-k wrapped PRIMARY graph the check
 *              of which contexts are palindromic k-mers, and on a graph without the mask the
 *              check of a pattern's few unchecked candidates (Request::max_checked_entries,
 *              after the other scans) — the only stop that, after complete discovery, can
 *              leave BOUNDS; every other discovery stop leaves AT_LEAST;
 *  EXTRACTION  enumerate()'s release of contexts or paths, stopped by the deadline after
 *              discovery (and extension) completed or stopped (its counts keep their
 *              relation: extraction counts nothing);
 *  EXTENSION   L > k with Request::extend_paths: the listing of the admitted anchors and the
 *              DFS beyond k (§4.2), stopped by max_steps, the deadline, or max_paths with
 *              stop_at_threshold, and with a support tracker or a path sink by them (EXTERNAL,
 *              or TIME when they found the work time passed through the request's Budget). The
 *              anchors stay EXACT, the paths become AT_LEAST.
 */
enum class StopPhase { DISCOVERY, MASK_SCAN, EXTRACTION, EXTENSION };

struct Stop {
    StopPhase phase;
    StopReason reason;
};

/**
 * The request's work budget: max_steps and the deadline, shared by all patterns of the request
 * in request order (§5.3: one share per shard; a request runs on one shard), so that the same
 * request on the same index stops at the same step (§5.5). A step is one range evaluation (a
 * tighten_range, successful or not, or a W-rule rank set), one edge examined by a scan, or one
 * outgoing edge examined by the extension (§4.2: every edge DeBruijnGraph::call_outgoing_kmers
 * reports for a node the DFS expands, allowed or not); extraction charges no steps (its work is
 * linear in the contexts it releases, which max_contexts caps, plus the invalid candidates it
 * skips) but reads the clock.
 * Stops are sticky: once a charge is refused, every later charge is refused with the same
 * reason, and the patterns after it answer UNKNOWN counts with that stop.
 *
 * Where the engine reads the clock: at the start of every pattern, at every charge that crosses
 * a multiple of kClockStride steps (discovery, the deferred scans and the extension share one
 * step counter, so the boundary from discovery to the scans has no reading of its own: at most
 * kClockStride - 1 steps pass unread there, as anywhere), before every release, every
 * kClockStride edges the release examines and every kReleaseClockStride contexts it passes on
 * (a released context costs the caller a k-mer spelling, k - 1 BOSS steps, and its result
 * object), and every kReleaseClockStride contexts an ALL_OR_COUNT delivery passes to the
 * caller; in the extension also every kReleaseClockStride anchors it lists, before every anchor
 * it extends (spelled first, k - 1 BOSS steps no step charges) and before every
 * kReleaseClockStride-th node its DFS expands; after a search that completed or stopped at a
 * threshold, before every piece but the first of the low-complexity diagnostic
 * (kNoteLowComplexity). Work therefore ends within one such stride after the work time passes.
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
    // not "$ACGT" or "$ACGTN"). The route narrows it further (cli::route_support):
    // "alphabet_untested" ($ACGTN, not served until a DNA5 build passes the pattern tests) and
    // "mask_invalid" (a mask marking a W = $ edge valid). A graph without its mask is served:
    // no reason is "mask_required"
    std::string reason;
    GraphMode mode = GraphMode::BASIC;
    /**
     * The graph has its dummy-edge (valid-edge) mask: every count of a completed discovery is
     * EXACT. Without it the graph is served all the same, its unresolved counts BOUNDS [lower,
     * U] (PatternSearch, "Graphs without the dummy-edge mask"); the mask is derived data of the
     * graph, never a different number for the same count: it only changes which counts are
     * exact.
     */
    bool mask_present = false;
    size_t k = 0;
    // the BOSS alphabet, sentinel first: "$ACGT" (DNA4) or "$ACGTN" (DNA5)
    std::string alphabet;
    // true on BASIC only: FORWARD/REVERSE are the + and - strands
    bool strand_stated = false;
    // the requestable scopes: SUFFIX and ANY_OFFSET, or ANY_OFFSET only (PRIMARY)
    std::vector<Scope> scopes;
};

/**
 * f, the fraction of real k-mers among the entries of a succinct graph that a pattern can
 * count: the BOSS edges whose W is not $ (a sink dummy, W = $, never carries a pattern base),
 * each a k-mer or a source dummy (a k-mer starting with '$', BOSS::node_has_sentinel). On a
 * graph without its dummy-edge mask a count is the bounds [lower, U] and the route states the
 * additive estimate U x f beside it: what the count would be if the source dummies among its U
 * entries were as frequent as in the whole graph (not a bound, never EXACT).
 * Over the whole graph f = real edges / (edges with W != $), where `stats --count-dummy`'s
 * real edges = edges - source dummies - sink dummies (its source dummies include the main
 * dummy edge 1, whose W is $, its sinks do not) and the edges with W != $ number edges - sink
 * dummies - 1. build/mini_refseq: 8,335,760 edges, 375 source and 12 sink dummies, f =
 * 8,335,373 / 8,335,747 = 0.99995513.
 */
struct RealFraction {
    // real / samples; 1 when the graph has no edge with W != $ (nothing can be counted)
    double value = 1;
    // Wilson's score interval at 95% (z = 1.959963984540054) around |value|, clamped to
    // [0, 1] (and holding |value| despite rounding); [0, 1] without samples
    double lower = 0;
    double upper = 1;
    // the entries drawn (W != $) and the real k-mers among them
    uint64_t samples = 0;
    uint64_t real = 0;
    // the edges with W = $ (plain or marked: the sink dummies and the main dummy edge 1),
    // counted exactly by ranks of W, and all edges (BOSS::num_edges)
    uint64_t sentinel_edges = 0;
    uint64_t edges = 0;
    // the seed of the draws: the graph's number of edges
    uint64_t seed = 0;
};

// the entries drawn for f: 10,000
constexpr uint64_t kRealFractionSamples = 10'000;

/**
 * Samples f (RealFraction) from |samples| entries drawn uniformly, with replacement, among
 * the edges whose W is not $: std::mt19937_64 seeded with the graph's number of edges (the
 * generator the standard specifies, so that the same graph gives the same f in every process
 * and on every platform), the edge 1 + rng() % edges (the modulo's bias is below edges / 2^64),
 * a draw with W = $ drawn again; a drawn edge is real iff its source node holds no '$'
 * (BOSS::node_has_sentinel: at most k - 1 symbols read). Reads the BOSS only, never the
 * mask. Cost: |samples| walks of at most k - 1 bwd steps (less with an index of suffix
 * ranges), plus a redraw per W = $ entry drawn: a few tens of milliseconds at k = 31 (see
 * test_pattern_unmasked.cpp, RealFractionCost). Pure: the same graph and |samples| give the
 * same result. The route samples it once per loaded graph.
 */
RealFraction sample_real_fraction(const DBGSuccinct &graph,
                                  uint64_t samples = kRealFractionSamples);


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
 *  NOT_REQUESTED  Request::extend_paths is false: nothing extended; paths
 *                 UNKNOWN, EXACT 0 when the anchors are EXACT 0.
 *  NO_ANCHORS     the anchors are EXACT 0: nothing to extend; paths EXACT 0 (no path starts
 *                 without an anchor: a derivation, not a promotion).
 *  NOT_STARTED    the anchors are not EXACT (a step, time or anchor-threshold stop in
 *                 discovery or a mask scan, see Result::stop): an anchor set not known
 *                 completely is never extended; paths UNKNOWN (§3: every phase after a stop).
 *  NOT_ADMITTED   the anchors are EXACT and above max_anchors: the extension's admission
 *                 failed (§4.2); paths UNKNOWN; enumerate() withholds ANCHORS_ABOVE_THRESHOLD.
 *                 On a graph without the dummy-edge mask also: the anchors are BOUNDS (no stop)
 *                 and their upper bound U is above max_anchors (the admission compares U;
 *                 kNoteThresholdUpperBound when the lower bound is not).
 *                 With U <= max_anchors the anchors are listed first, which drops the source
 *                 dummies and makes their count EXACT (NO_ANCHORS when none is left), and then
 *                 extended.
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
     * With a support tracker (Request::support) a plain count of the complete walks, supported
     * or not: EXACT only when COMPLETED and no branch was pruned before it reached L
     * (pruned_before_completion false), AT_LEAST (the walks completed) otherwise; a branch
     * pruned at its last k-mer is a complete walk and counted.
     */
    Count paths = Count::unknown(Unit::PATHS);
    /**
     * With extend_paths, one entry per orientation searched (empty without it): EXACT when
     * every anchor of that orientation was extended (EXACT 0 when its anchors are EXACT 0),
     * AT_LEAST when the extension stopped before, UNKNOWN when it did not run. |paths| is
     * their sum with the weakest relation, except that a total whose extension did not run
     * is UNKNOWN (EXACT 0 with NO_ANCHORS) rather than AT_LEAST over its known zeros. With a
     * support tracker, AT_LEAST also when a branch of that orientation was pruned before L.
     */
    std::map<Orientation, Count> paths_by_orientation;
    /**
     * The supported-path search (Request::support). |walks|: the complete walks the extension
     * counted, a plain number (paths.value once the extension ran, with or without a tracker; 0
     * when it did not). |supported|: the complete walks the tracker supported (unit PATHS),
     * with |paths|'s relation rule for the extension but without the pruning (a pruned branch
     * holds no supported path): EXACT when COMPLETED, AT_LEAST when STOPPED, EXACT 0 without
     * anchors, UNKNOWN when the extension did not run and always without a tracker; per
     * orientation in |supported_by_orientation| (empty without a tracker), summed as |paths|.
     * |branches_pruned|: the tracker's DEAD verdicts — anchors and branches whose support ran
     * out before L, and complete walks without support (work, exact as such).
     * |pruned_before_completion|: some anchor or branch was pruned before L, so that the walks
     * below it were not counted (|paths| AT_LEAST).
     */
    uint64_t walks = 0;
    Count supported = Count::unknown(Unit::PATHS);
    std::map<Orientation, Count> supported_by_orientation;
    uint64_t branches_pruned = 0;
    bool pruned_before_completion = false;
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
    // BOSS steps, per context, one step each): there about one per such range. On a graph
    // without the dummy-edge mask only these palindrome checks exist (each also tells a source
    // dummy from a k-mer), and the ranges of a pattern whose few unchecked candidates were
    // checked (Request::max_checked_entries: k - 1 steps per candidate)
    uint64_t mask_scans = 0;
    // L > k with extend_paths: the outgoing edges the extension examined, one step each
    // (allowed or not); 0 otherwise
    uint64_t extension_edges = 0;
    // L > k with extend_paths, counters of the extension's work beside its steps (not steps,
    // never charged): the anchors whose extension began, each spelled once (k - 1 BOSS steps no
    // step charges) before its DFS, at most the anchors listed; and the DFS's branchings, the
    // nodes it expanded (anchors included) with two or more outgoing k-mers allowed at the next
    // position. 0 when the extension did not run
    uint64_t extension_anchors = 0;
    uint64_t extension_branches = 0;
    // every step this pattern charged: ranges_visited, plus the edges its scans examined,
    // plus k - 1 per candidate checked (Request::max_checked_entries), plus extension_edges
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
    // ALL_OR_COUNT: discovery (and for L > k the extension) completed, EXACT
    // total > max_contexts (L <= k), or EXACT paths > max_paths (L > k with extend_paths). On a
    // graph without the dummy-edge mask also a BOUNDS total whose upper bound U is above
    // max_contexts (the admission compares U, conservative;
    // kNoteThresholdUpperBound when the lower bound is not above it)
    COUNT_ABOVE_THRESHOLD,
    // ALL_OR_COUNT: stop_at_threshold stopped discovery or the extension (stop reason
    // MAX_CONTEXTS, MAX_ANCHORS or MAX_PATHS)
    THRESHOLD_CROSSED,
    // ALL_OR_COUNT: discovery, a mask scan or the extension stopped at max_steps
    DISCOVERY_BUDGET,
    // ALL_OR_COUNT: discovery, the extension or extraction stopped by the deadline (§5.3)
    DEADLINE,
    // either mode, L > k without extend_paths, anchors not EXACT 0: results of a long pattern
    // are its completed paths (§4.2), which only extend_paths extends; anchors are never
    // released as results
    PATHS_LATER_INCREMENT,
    // either mode, L > k with extend_paths: the anchors are EXACT and above max_anchors, so
    // the extension was not admitted (§4.2); counts.anchors is exact, the paths UNKNOWN
    ANCHORS_ABOVE_THRESHOLD,
    // ALL_OR_COUNT, L > k with a support tracker: the tracker stopped the extension
    // (StopReason::EXTERNAL); the route writes the tracker's reason (annotation_budget), never
    // "external"
    EXTERNAL,
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
     * EXACT 0, hence paths EXACT 0; with extend_paths, paths EXACT and all of them returned).
     * JSON retrieval_complete: the one flag that licenses an absence claim over graph contexts
     * (§5.1); it claims nothing about labels, which the engine never reads.
     */
    bool complete = false;
    // set: nothing is published ("results": [], returned 0). The callback received nothing,
    // except when the deadline stopped an ALL_OR_COUNT delivery (DEADLINE, stop {EXTRACTION,
    // TIME}): it then received a prefix of the contexts, which the caller must discard
    std::optional<Withheld> withheld;
    // PARTIAL, neither complete nor withheld: why the list is shorter than the pattern's
    // contexts — the reason of the stop that touched the pattern when there is one, else
    // MAX_CONTEXTS (the first max_contexts in answer order were returned; MAX_ANCHORS for
    // the anchors of Request::release_anchors; MAX_PATHS for paths; EXTERNAL for the supported
    // paths completed before a support tracker stopped the extension)
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
 * L > k with Request::extend_paths: one path (§4.2), (orientation, path) its identity: offset
 * 0, |node| and |base_node| its anchor (the first k-mer), |path| its n = L - k + 1 nodes and
 * |sequence| its L spelled bases, which instantiate P for FORWARD and PALINDROMIC and rc(P) for
 * REVERSE (§7.2: the result's field sequence is |sequence|, its field anchor_kmer the anchor's
 * k-mer; kmer is not used for a path).
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


// ---------------------------------------------------------------- the search state

/**
 * The arm a step of the extension takes: RIGHT along outgoing edges (every extension: each
 * oriented pattern is read from its first k-window on), LEFT along incoming edges (a selective
 * anchor: not built).
 */
enum class Side : uint8_t { RIGHT, LEFT };

/**
 * One state of the extension's DFS, as a SupportTracker sees it:
 * the oriented node entered, the position in the oriented pattern, the Model's state and the
 * tracker's frame. Built only when a tracker is given (Request::support).
 */
struct SearchState {
    // the k-mer entered (the anchor, or the child a step chose), as the served graph names it
    DeBruijnGraph::node_index node;
    // the stored k-mer that carries its annotation row (PatternSearch::base_node)
    DeBruijnGraph::node_index base_node;
    // the stored k-mer is the reverse complement of the node's k-mer (so that a tracker passes
    // support_step the row's k-mer as rc of the k-mer spelled, without the wrapper's id
    // arithmetic): true on the wrapper of a PRIMARY graph for a node that is not stored as
    // spelled (base_node != node), false on BASIC and CANONICAL graphs. A palindromic k-mer is
    // its own reverse complement: either value names the same row
    bool stored_reverse_complement;
    // the oriented pattern extended: FORWARD and PALINDROMIC read P, REVERSE rc(P)
    Orientation orientation;
    Side side;
    // the base the step spelled (the node's last base on the right arm); '\0' for an anchor
    char base;
    // the bases of the oriented pattern spelled, this step's included: k at an anchor
    uint32_t position;
    // the Model's state after them
    Model::State model;
    // the tracker's frame: 0 at the anchor, one more per step (the level of the DFS)
    uint32_t depth;
    // the bases spelled, in reading order (the anchor's k-mer, then one per step): the node's
    // k-mer is its last k bases. A view into the DFS's own string, valid during the call only
    std::string_view spelled;
};

/**
 * A complete walk as the extension hands it to a tracker and a sink (valid during the call
 * only): nothing is copied before a sink admits it. copy() is the Context enumerate() releases.
 */
struct PathView {
    // the walk's anchor (offset 0, node == path.front())
    const Context &anchor;
    // its n = L - k + 1 k-mers in reading order
    const std::vector<DeBruijnGraph::node_index> &path;
    // the L bases it spells
    std::string_view sequence;

    Context copy() const {
        return Context { anchor.orientation, 0, anchor.node, anchor.base_node, path,
                         std::string(sequence) };
    }
};

/**
 * The support of each branch of the extension: one frame per level of the DFS. Implemented
 * outside the engine (the label and trace trackers over the annotation, on
 * graph/traversal/support_step; not built yet); the engine only asks it and counts its
 * verdicts. With Request::support null the extension asks nothing.
 *
 * The protocol of an admitted extension:
 *  - prepare(anchors): once, before the first anchor is opened, with every anchor the extension
 *    will open (both orientations', in answer order), so that a tracker can read what it needs
 *    of all of them first (the labels any of their walks can be supported by); false stops the
 *    extension before its first anchor, as a STOPPED verdict;
 * then per anchor, in answer order:
 *  - open(anchor): frame 0 from the anchor (depth 0, position k), before its first expansion;
 *  - push(child): frame d + 1 from frame d and the child entered (depth d + 1), right after
 *    the DFS chose it (its edge already examined and charged as a step) and before it is
 *    expanded or, at L, completed;
 *  - complete(walk): the walk of the frames open reached L (the Model accepting); the walk is
 *    a supported path iff it answers ALIVE (then the sink, or the engine's list, receives it);
 *  - pop(): drops the top frame — when the DFS leaves a child it pushed (after its subtree, or
 *    after its completion), and frame 0 when the anchor is done. Every open() and push() that
 *    answered ALIVE is popped exactly once, also when the extension stops (the frames unwind
 *    before extend returns); a DEAD or STOPPED verdict leaves no frame to pop.
 * Verdicts: ALIVE (some support survives: go on), DEAD (no support: the branch is pruned —
 * not expanded, not completed, counted in AnchorCounts::branches_pruned), STOPPED (the
 * tracker's own budget, memory or the deadline: the extension ends with stop {EXTENSION,
 * EXTERNAL}, or {EXTENSION, TIME} when the tracker found the work time passed through the
 * request's Budget::check_time, which then holds TIME; the route writes stop_reason()).
 * The extension's own work is unchanged: every outgoing edge of an expanded node is a step,
 * the clock is read as without a tracker; a pruned branch is not expanded, so its edges are
 * never examined. The tracker charges its own work (rows, merges) to its own budgets and
 * reads the deadline before its own reads (never Budget::charge). pop() must not throw.
 */
class SupportTracker {
  public:
    enum class Verdict : uint8_t { ALIVE, DEAD, STOPPED };

    virtual ~SupportTracker() = default;

    // the default reads nothing
    virtual bool prepare(const std::vector<Context> &anchors) {
        (void)anchors;
        return true;
    }
    virtual Verdict open(const SearchState &anchor) = 0;
    virtual Verdict push(const SearchState &child) = 0;
    virtual void pop() = 0;
    virtual Verdict complete(const PathView &walk) = 0;
    // after a STOPPED verdict: the reason the route writes (max_annotation_work, max_memory,
    // time); nullptr before
    virtual const char* stop_reason() const = 0;
};

/**
 * What a complete walk becomes: the route's list of supported paths, an alignment projection or
 * a ranked list (none built yet). Given (Request::sink), it receives every path the extension
 * completes — with a tracker the supported ones, the tracker's frames still open so that it can
 * read their support — in answer order, through count(): the engine then keeps no path of its
 * own, and enumerate() refuses an extending pattern with a sink. Null: the engine's own list
 * (enumerate() only), with the retention rule of Request::max_paths. The engine applies
 * stop_at_threshold (more than max_paths paths complete: stop {EXTENSION, MAX_PATHS}) after the
 * sink's call either way; the retention of a given sink is its own (it admits a path's memory
 * before it copies it: PathView::copy).
 * The engine's threshold counts the paths it handed over (with a tracker the supported ones),
 * not what the sink kept. A sink that selects among them (a predicate over supported paths: the
 * threshold is on the selected paths) owns max_paths itself: the route then leaves
 * Request::stop_at_threshold off, and the sink answers false once more than max_paths paths are
 * selected, with stop_reason() "max_paths" (the route writes stop {extension, max_paths}).
 */
class PathSink {
  public:
    virtual ~PathSink() = default;

    // false: the sink stopped the extension (stop {EXTENSION, EXTERNAL}, or TIME as for a
    // tracker), after which it is asked nothing more
    virtual bool accept(const PathView &path, const SupportTracker *support) = 0;
    // after a false accept(): the reason the route writes; nullptr before
    virtual const char* stop_reason() const = 0;
    // after a false accept(): whether that stop is a threshold of the request's
    // (stop_at_threshold, which a selecting sink owns: the selected paths above max_paths),
    // not a budget; the engine treats it as its own threshold stops (kNoteLowComplexity)
    virtual bool stopped_at_threshold() const { return false; }
};

// JSON notes of a pattern (§7.2), the ones the engine states:
//  low_complexity_pattern    an exact pattern that sdust flags with the seeder's parameters
//                            (T = 20, W = 64, is_low_complexity): why its counts are large. An
//                            optional diagnostic: never stated beside a budget stop (max_steps,
//                            time, or a support tracker's or path sink's own: any phase, the
//                            request-wide stop of an earlier pattern included), which leaves no
//                            room for optional work; beside a threshold stop (max_contexts,
//                            max_anchors, max_paths, or a path sink's threshold:
//                            PathSink::stopped_at_threshold) it is diagnosed, those being the
//                            answers whose counts are large. Not stated either when the work time passed
//                            before sdust had read the pattern (read in pieces of 128 bases, the
//                            clock before each but the first, so only a pattern longer than 191
//                            bases can lose it so; Result::time_limited is then set, the counts
//                            as they were). Otherwise stated iff sdust flags the whole pattern
//  strand_unknown_canonical  graph mode CANONICAL or PRIMARY: orientations, not strands
//  paths_later_increment     L > k without Request::extend_paths: anchors counted, paths
//                            neither extended nor extracted (never set with extend_paths)
//  threshold_upper_bound     a graph without the dummy-edge mask: a
//                            threshold decision went against the request on the count's
//                            upper bound U while its lower bound did not cross the threshold
//                            (ALL_OR_COUNT withheld COUNT_ABOVE_THRESHOLD, stop_at_threshold
//                            stopped, the extension NOT_ADMITTED): the true count may be
//                            within the threshold (PARTIAL lists the contexts regardless)
//  no_stop_codon             a peptide holding '*' read in a genetic code without an
//                            unconditional stop codon (tables 27, 28, 31):
//                            '*' matches nothing there, so the pattern has no instance: its
//                            contexts (L <= k) and paths (L > k) are 0 for that reason (never a
//                            silent 0); a long one's anchors are the anchor window's, which may
//                            not reach the '*'. Set on every answered such pattern
constexpr const char kNoteLowComplexity[] = "low_complexity_pattern";
constexpr const char kNoteStrandUnknown[] = "strand_unknown_canonical";
constexpr const char kNotePathsLater[] = "paths_later_increment";
constexpr const char kNoteThresholdUpperBound[] = "threshold_upper_bound";
constexpr const char kNoteNoStopCodon[] = "no_stop_codon";

/**
 * The answer for one pattern. With |refusal| set nothing was searched and only the pattern
 * description (scope, palindromic, information bits) is meaningful.
 */
struct Result {
    std::optional<Refusal> refusal;

    // SUFFIX | ANY_OFFSET as requested, LONG when L > k
    Scope scope = Scope::ANY_OFFSET;
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
    // L > k: the least over the searched orientations' anchor windows (FORWARD and PALINDROMIC:
    // P[0, k); REVERSE: rc(P)[0, k), whose bits are those of P[L - k, L)), the bits the
    // information floor gated on (an addition to the JSON of contract version 1)
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
    // determinism "time_limited", else "full"; §5.5). Also set when the work time passed
    // during the low-complexity diagnostic of a pattern of more than 191 bases (its search
    // completed or stopped at a threshold), which then leaves its note out
    // (kNoteLowComplexity): only that note depends on the machine there
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
 *    k-mers), so a path exists only where all its k-mers were retained; without the mask the
 *    only dummy an outgoing edge of a k-mer can reach is a sink (base '$'), which no pattern
 *    position allows.
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
 * orientation whose run comes late complete. A peptide's leading X residues are searched as
 * given on every graph (an X codon is not every 3-mer: its stops are excluded), at the cost of
 * a codon's worth of ranges each, as an IUPAC N^3 inside a window.
 *
 * Peptides (§6): a window of a peptide is searched with its codon automaton: at every depth of
 * the range DFS the bases tried are Pattern::allowed()'s after the codon bases the range's
 * nodes already end with, read from the BOSS (the last by its F array, the one before through
 * bwd), so no prefix state travels outside the range; the W rule and the extension ask the
 * same automaton. On a wrapped PRIMARY graph the reverse complement of a long peptide's anchor
 * window, rc(Q[0, k)) = rc(Q)[L - k, L), may start inside a codon: its first codon is matched
 * on the bases inside the window only (the codon's earlier bases are free).
 *
 * Graphs without the dummy-edge mask (§4.4). The mask only says which BOSS edges are dummies;
 * without it the same ranges are discovered with the same steps, and every place that read the
 * mask counts the candidates instead:
 *  - a flank range counts its non-sink edges (W != $, DBGSuccinct::count_non_sink_edges_in_
 *    range) instead of its valid ones, a W-rule leaf its candidates (W in {c, c + alph_size})
 *    instead of the valid ones among them: no INVALID scan exists, a sink never counts;
 *  - a range whose nodes the search spelled wholly (depth k - 1: every node symbol a pattern
 *    base or a flank symbol, so no room for '$') holds no source dummy: its count is exact.
 *    Any other range's candidates may include source dummies (k-mers starting with '$'):
 *    they are UNCHECKED, counted into the upper bound U and not into the lower bound;
 *  - the palindrome scans of an even-k wrapped PRIMARY graph spell every candidate anyway,
 *    which checks it too (a k-mer holding '$' is a dummy, and never a palindrome).
 *  - a pattern whose discovery and scans completed with at most Request::max_checked_entries
 *    unchecked candidates in all has each of them tested after the scans
 *    (BOSS::node_has_sentinel, k - 1 steps each, phase MASK_SCAN): the real ones are counted
 *    exactly and the rest dropped, so every count of the pattern is EXACT. Its unchecked ranges
 *    are kept for that while their candidates number at most the limit (at most that many
 *    ranges), and freed once there are more.
 * A count is then EXACT when nothing of it is unchecked (an empty block, U = 0, included:
 * absence holds), else BOUNDS {lower, U}; stops as on a masked graph (AT_LEAST, UNKNOWN). The
 * admission decisions compare U (Request::stop_at_threshold, ALL_OR_COUNT's threshold, the
 * extension's admission): conservative, stated with kNoteThresholdUpperBound when the lower
 * bound did not cross (after a check the count is EXACT, and they compare it as with the
 * mask; while the unchecked candidates are at most the limit, the retention for ALL_OR_COUNT
 * and the extension compares the lower bound, so that what the check admits is still there to
 * release). The lists stay exact: the release tests every unchecked candidate with
 * BOSS::node_has_sentinel (at most k - 1 symbols read, about the cost of the spelling the
 * route does per context anyway) and never releases a source dummy, and a release that
 * enumerated every candidate makes the counts EXACT (enumerate()). Memory: PARTIAL's
 * retention bound counts only the contexts it is sure of, which an unchecked range is not:
 * on such a graph PARTIAL keeps every discovered range whose nodes are not wholly spelled (24
 * bytes each, at most one per step charged) rather than about max_contexts of them.
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
     * frontier, O(k * alphabet); on a graph without the mask also the unchecked ranges kept for
     * the check of Request::max_checked_entries (at most that many). Refusals (information
     * floor, SUFFIX on a wrapped PRIMARY graph) come
     * back as Result::refusal without charging anything. On a budget already stopped (or
     * whose work time has passed at the pattern's check_time) every count is UNKNOWN and
     * stop is {DISCOVERY, the budget's reason}. A pattern without instances
     * (Pattern::has_instances: '*' in a table without a stop codon) of L <= k is answered
     * EXACT 0 in every count, before the information floor and without reading the budget
     * (nothing is searched, nothing charged); one of L > k is searched as any other (its
     * anchors instantiate the anchor window only, which may not reach the '*'; the extension
     * finds no path); both carry the note kNoteNoStopCodon. Never throws for a parsed
     * pattern, except
     * std::logic_error on a broken internal invariant (never expected: an answer that
     * cannot be stated correctly is not stated at all).
     * L > k with request.extend_paths: the anchors are retained (§5.2: kept through the
     * extension's admission, dropped once their count exceeds max_anchors) and, when
     * admitted, extended to count the paths (AnchorCounts::paths, Extension); no path is
     * retained. The extension's memory is the anchor list (<= max_anchors) and the DFS
     * stack, O(n * alphabet) for n = L - k + 1.
     * With Request::support and Request::sink the extension asks the tracker per level and
     * hands every (supported) complete walk to the sink as it completes; the counts are those
     * of AnchorCounts (walks, supported, branches_pruned).
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
     *  With a support tracker (Request::support) the paths are the supported ones and every
     *  rule above reads their count (AnchorCounts::supported, EXACT when the extension
     *  COMPLETED whatever it pruned); a tracker's stop is withheld EXTERNAL (ALL_OR_COUNT) or
     *  cut EXTERNAL after the prefix completed before it (PARTIAL). With a path sink
     *  (Request::sink) the sink is the release: an extending pattern is refused
     *  (std::invalid_argument; count() serves it).
     * On a graph without the dummy-edge mask (see the class): ALL_OR_COUNT releases iff
     * discovery completed and the count is EXACT, or BOUNDS with U <= the threshold; the
     * release then enumerates every candidate, drops the source dummies, and the counts
     * (total, suffix, by_offset, by_orientation; anchors with release_anchors) become EXACT,
     * the number of contexts released (a number outside [lower, U] is std::logic_error).
     * PARTIAL: a release that enumerated every candidate of a completed discovery (nothing
     * pruned) makes the counts EXACT likewise; one cut by its cap or after a stop raises each
     * count's lower bound to the contexts it released at that orientation and offset. With
     * extend_paths, anchors BOUNDS and U <= max_anchors, the anchors are listed (dummies
     * dropped, their counts EXACT) and extended. count() never enumerates: its counts stay
     * BOUNDS, and the steps charged are the same in both. The check of a pattern with at most
     * Request::max_checked_entries unchecked candidates runs in both alike (before any
     * release): its counts are then EXACT in count() as in enumerate(), the same steps
     * charged.
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

    // a pattern with no instance (Pattern::has_instances) and L <= k:
    // every count EXACT 0, nothing charged
    void answer_no_instance(const Pattern &pattern, const Request &request, bool enumerating,
                            Result *result) const;
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

// the answer's `kind` and the request's pattern key
inline const char* to_string(PatternKind kind) {
    switch (kind) {
        case PatternKind::DNA: return "dna";
        case PatternKind::IUPAC: return "iupac";
        case PatternKind::PROTEIN: return "protein";
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
        // never written by the route, which names the tracker's or the sink's reason
        case StopReason::EXTERNAL: return "external";
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
        // never written by the route (see Withheld::EXTERNAL)
        case Withheld::EXTERNAL: return "external";
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

// the request's `mode`, `scope` (LONG is not requestable) and `strands`; nullopt: not a valid
// value. A valid mode is still refused by the route when its projection is not served (§7.1).
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
