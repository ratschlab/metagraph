#include "pattern_search.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <queue>
#include <random>
#include <sstream>
#include <tuple>
#include <type_traits>

#include <sdust.h>

#include "aligner_seeder_methods.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {
namespace pattern {

using boss::BOSS;
using node_index = DeBruijnGraph::node_index;
using edge_index = BOSS::edge_index;
using TAlphabet = BOSS::TAlphabet;


// ---------------------------------------------------------------- patterns

namespace {

BaseSet iupac_set(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return kBaseA;
        case 'C': return kBaseC;
        case 'G': return kBaseG;
        case 'T': return kBaseT;
        case 'R': return kBaseA | kBaseG;
        case 'Y': return kBaseC | kBaseT;
        case 'S': return kBaseC | kBaseG;
        case 'W': return kBaseA | kBaseT;
        case 'K': return kBaseG | kBaseT;
        case 'M': return kBaseA | kBaseC;
        case 'B': return kBaseC | kBaseG | kBaseT;
        case 'D': return kBaseA | kBaseG | kBaseT;
        case 'H': return kBaseA | kBaseC | kBaseT;
        case 'V': return kBaseA | kBaseC | kBaseG;
        case 'N': return kAllBases;
        default: return 0;
    }
}

// the IUPAC letter of a non-empty set (index: the set's bits A=1, C=2, G=4, T=8)
char iupac_letter(BaseSet set) {
    static constexpr char kLetters[] = "?ACMGRSVTWYHKDBN";
    assert(set && set <= kAllBases);
    return kLetters[set];
}

// A <-> T, C <-> G, bit by bit, so that every IUPAC set maps to its complement's set
// (R <-> Y, K <-> M, B <-> V, D <-> H; S, W, N to themselves), as COMPL_TAB does for letters
BaseSet complement_set(BaseSet set) {
    return ((set & kBaseA) << 3) | ((set & kBaseT) >> 3)
         | ((set & kBaseC) << 1) | ((set & kBaseG) >> 1);
}

uint32_t set_size(BaseSet set) {
    return __builtin_popcount(set);
}

// the bases a codon set admits at codon position |j| (0, 1, 2), whatever the others
BaseSet codon_position_set(CodonSet set, size_t j) {
    return next_bases(set, j);
}

// a peptide's residue (upper case) as a codon set of |code|: the 20 residues, X (every codon
// that is not a stop), B (D or N), Z (E or Q), J (I or L); 0 for anything else
CodonSet residue_codons(char residue, const GeneticCode &code) {
    switch (residue) {
        case 'X': return kAllCodons & ~code.stops();
        case 'B': return code.codons('D') | code.codons('N');
        case 'Z': return code.codons('E') | code.codons('Q');
        case 'J': return code.codons('I') | code.codons('L');
        default:
            return residue && std::strchr("ACDEFGHIKLMNPQRSTVWY", residue)
                ? code.codons(residue)
                : 0;
    }
}

} // namespace

Pattern Pattern::parse(PatternKind kind, std::string_view text) {
    if (kind == PatternKind::PROTEIN)
        return parse(kind, text, GeneticCode::standard());

    if (text.empty())
        throw PatternError("bad_alphabet", "pattern: empty");

    std::string upper;
    upper.reserve(text.size());
    std::vector<BaseSet> positions;
    positions.reserve(text.size());

    for (size_t i = 0; i < text.size(); ++i) {
        BaseSet set = iupac_set(text[i]);
        // a DNA pattern is a string over A, C, G, T only; IUPAC admits the 15 codes
        if (!set || (kind == PatternKind::DNA && set_size(set) != 1)) {
            std::ostringstream msg;
            msg << "pattern: character '" << text[i] << "' at position " << i
                << " is not in the " << to_string(kind) << " alphabet ("
                << (kind == PatternKind::DNA ? "A, C, G, T" : "A C G T R Y S W K M B D H V N")
                << ")";
            throw PatternError("bad_alphabet", msg.str());
        }
        upper.push_back(std::toupper(static_cast<unsigned char>(text[i])));
        positions.push_back(set);
    }

    return Pattern(kind, std::move(upper), std::move(positions));
}

Pattern Pattern::parse(PatternKind kind, std::string_view text, const GeneticCode &code) {
    if (kind != PatternKind::PROTEIN)
        return parse(kind, text);

    if (text.empty())
        throw PatternError("bad_alphabet", "pattern: empty");

    std::string upper;
    upper.reserve(text.size());
    std::vector<CodonSet> codons;
    codons.reserve(text.size());
    std::vector<BaseSet> positions;
    positions.reserve(3 * text.size());

    bool has_instances = true;
    for (size_t i = 0; i < text.size(); ++i) {
        const char residue = std::toupper(static_cast<unsigned char>(text[i]));
        // the stop '*': the table's stop codons, none in a table without an unconditional stop
        // (27, 28, 31), whose context stops code their residue: the peptide then has no
        // instance
        const CodonSet set = residue == '*' ? code.stops() : residue_codons(residue, code);
        if (!set && residue != '*') {
            std::ostringstream msg;
            // ('*', a residue too, is not listed: the message is part of the answers)
            msg << "pattern: character '" << text[i] << "' at position " << i
                << " is not in the protein alphabet (A C D E F G H I K L M N P Q R S T V W Y, "
                   "and X B Z J)";
            throw PatternError("bad_alphabet", msg.str());
        }
        has_instances &= set != 0;
        upper.push_back(residue);
        codons.push_back(set);
        for (size_t j = 0; j < 3; ++j) {
            positions.push_back(codon_position_set(set, j));
        }
    }

    Pattern pattern(kind, std::move(upper), std::move(positions));
    pattern.codons_ = std::move(codons);
    pattern.genetic_code_ = code.id();
    pattern.has_instances_ = has_instances;
    return pattern;
}

BaseSet Pattern::allowed(size_t position, std::string_view spelled) const {
    assert(position < positions_.size());
    if (kind_ != PatternKind::PROTEIN)
        return positions_[position];

    // the codon automaton (§6): the residue's codons agreeing with the bases spelled so far
    // in this codon
    const size_t j = position % 3;
    int known[2] = { -1, -1 };
    for (size_t i = 0; i < j && i < spelled.size(); ++i) {
        const int b = base_index(spelled[spelled.size() - 1 - i]);
        if (b < 0)
            return 0;
        known[j - 1 - i] = b;
    }
    return next_bases(codons_[position / 3], j, known[0], known[1]);
}

namespace {

// a peptide's Model state (Pattern::start, next): the bases spelled so far in the current codon,
// codon position i in bits [3i, 3i + 3) as its base index + 1 (0: not known, as allowed() takes
// a base before the spelled ones), and a flag for a base other than A, C, G, T among them, after
// which allowed() admits nothing
constexpr Model::State kCodonDead = Model::State(1) << 6;

int codon_prefix_base(Model::State s, size_t i) {
    return static_cast<int>((s >> (3 * i)) & 7) - 1;
}

} // namespace

Model::State Pattern::start(std::string_view anchor_kmer) const {
    if (kind_ != PatternKind::PROTEIN)
        return 0;

    // the last j bases of the anchor, as allowed(|anchor_kmer|.size(), anchor_kmer) reads them
    const size_t j = anchor_kmer.size() % 3;
    State s = 0;
    for (size_t i = 0; i < j && i < anchor_kmer.size(); ++i) {
        const int b = base_index(anchor_kmer[anchor_kmer.size() - 1 - i]);
        if (b < 0)
            return kCodonDead;
        s |= State(b + 1) << (3 * (j - 1 - i));
    }
    return s;
}

BaseSet Pattern::bases(State s, uint32_t position) const {
    assert(position < positions_.size());
    if (kind_ != PatternKind::PROTEIN)
        return positions_[position];
    if (s & kCodonDead)
        return 0;

    return next_bases(codons_[position / 3], position % 3, codon_prefix_base(s, 0),
                      codon_prefix_base(s, 1));
}

Model::State Pattern::next(State s, uint32_t position, char base) const {
    assert(position < positions_.size());
    // |base| is one the state admits (the precondition of Model::next)
    assert(bases(s, position) & iupac_set(base));
    if (kind_ != PatternKind::PROTEIN)
        return 0;

    const size_t j = position % 3;
    // the codon's last base: the next codon starts with nothing spelled
    if (j == 2)
        return 0;
    // base_index()'s codes, inline (once per child the extension enters); |base| is A, C, G or T
    State b = 0;
    switch (base) {
        case 'A': b = 0; break;
        case 'C': b = 1; break;
        case 'G': b = 2; break;
        case 'T': b = 3; break;
        default: return kCodonDead;
    }
    return (s & ~(State(7) << (3 * j))) | (b + 1) << (3 * j);
}

double Pattern::information_bits(size_t begin, size_t end) const {
    assert(begin <= end && end <= positions_.size());
    if (kind_ == PatternKind::PROTEIN) {
        // per residue covered: the positions of [begin, end) inside its codon, and the
        // number of distinct strings its codons spell there
        double bits = 0;
        for (size_t i = begin / 3; 3 * i < end; ++i) {
            const size_t first = std::max(begin, 3 * i) - 3 * i;
            const size_t last = std::min(end, 3 * i + 3) - 3 * i;
            // a residue admitting no codon (has_instances() false) as one exact codon: finite
            bits += 2.0 * (last - first)
                    - std::log2(static_cast<double>(std::max<uint32_t>(
                            1, distinct_projections(codons_[i], first, last))));
        }
        return bits;
    }
    // log2(4 / |set|) per set size, each computed once by the expression a per-position sum
    // would use, so that the sums are the same doubles bit for bit (one table read per base
    // instead of a log2: the O(L) pass of a long pattern)
    static const std::array<double, 5> kBits {
        0, std::log2(4.0 / 1u), std::log2(4.0 / 2u), std::log2(4.0 / 3u), std::log2(4.0 / 4u)
    };
    double bits = 0;
    for (size_t i = begin; i < end; ++i) {
        bits += kBits[set_size(positions_[i])];
    }
    return bits;
}

bool Pattern::is_exact() const {
    return std::all_of(positions_.begin(), positions_.end(),
                       [](BaseSet set) { return set_size(set) == 1; });
}

Pattern Pattern::reverse_complement() const {
    std::vector<BaseSet> positions(positions_.rbegin(), positions_.rend());
    if (kind_ == PatternKind::PROTEIN) {
        // the reverse-complemented automaton (§6): residues reversed, each codon set
        // reverse-complemented; its per-position unions are those of P complemented and
        // reversed
        for (BaseSet &set : positions) {
            set = complement_set(set);
        }
        Pattern rc(kind_, std::string(text_.rbegin(), text_.rend()), std::move(positions));
        rc.codons_.reserve(codons_.size());
        for (auto it = codons_.rbegin(); it != codons_.rend(); ++it) {
            rc.codons_.push_back(reverse_complement_codons(*it));
        }
        rc.genetic_code_ = genetic_code_;
        rc.has_instances_ = has_instances_;
        return rc;
    }
    std::string text;
    text.reserve(positions.size());
    for (BaseSet &set : positions) {
        set = complement_set(set);
        text.push_back(iupac_letter(set));
    }
    return Pattern(kind_, std::move(text), std::move(positions));
}

bool Pattern::is_palindromic() const {
    if (kind_ == PatternKind::PROTEIN) {
        // the instances are the product of the residues' codon sets: equal to rc's product
        // iff residue by residue (fixed-length factors)
        for (size_t i = 0, j = codons_.size(); i < codons_.size(); ++i) {
            if (codons_[i] != reverse_complement_codons(codons_[--j]))
                return false;
        }
        return true;
    }
    for (size_t i = 0, j = positions_.size(); i < positions_.size(); ++i) {
        if (positions_[i] != complement_set(positions_[--j]))
            return false;
    }
    return true;
}


// ---------------------------------------------------------------- deadline and budget

Deadline::Deadline(Clock::time_point start, double time_budget_ms, double finalize_reserve_ms,
                   std::function<Clock::time_point()> clock)
      : start_(start),
        time_budget_ms_(time_budget_ms),
        finalize_reserve_ms_(finalize_reserve_ms),
        clock_(clock ? std::move(clock) : std::function<Clock::time_point()>(&Clock::now)) {}

Deadline Deadline::unbounded() {
    Deadline deadline(Clock::now(), std::numeric_limits<double>::infinity(), 0);
    deadline.unbounded_ = true;
    return deadline;
}

bool Deadline::work_expired() const {
    return !unbounded_ && elapsed_ms() >= time_budget_ms_ - finalize_reserve_ms_;
}

bool Deadline::respond_expired() const {
    return !unbounded_ && elapsed_ms() >= time_budget_ms_;
}

double Deadline::elapsed_ms() const {
    return std::chrono::duration<double, std::milli>(clock_() - start_).count();
}

Budget::Budget(uint64_t max_steps, Deadline deadline)
      : max_steps_(max_steps), deadline_(std::move(deadline)) {}

bool Budget::charge(uint64_t n) {
    if (stopped_)
        return false;

    // the step cap first: it does not depend on the machine, so a request stops at the
    // same step on every run (§5.5)
    if (n > max_steps_ - steps_used_) {
        stopped_ = StopReason::MAX_STEPS;
        return false;
    }

    // the clock is read whenever the charge crosses a multiple of the stride (§5.3)
    if (steps_used_ / kClockStride != (steps_used_ + n) / kClockStride) {
        poll_abort();
        if (deadline_.work_expired()) {
            stopped_ = StopReason::TIME;
            return false;
        }
    }

    steps_used_ += n;
    return true;
}

bool Budget::check_time() {
    poll_abort();
    if (!deadline_.work_expired())
        return true;

    if (!stopped_)
        stopped_ = StopReason::TIME;

    return false;
}



// ---------------------------------------------------------------- the engine

namespace {

// a BOSS edge range: first, last (inclusive, whole node groups) and the length of the
// suffix its source nodes share (the DFS depth); the BOSSEdgeRange of suffix_to_prefix
typedef std::tuple<edge_index, edge_index, size_t> Range;

constexpr uint64_t kNoLimit = std::numeric_limits<uint64_t>::max();

// the seeder's sdust parameters (aligner_seeder_methods.cpp, T = 20, W = 64), and the bases
// is_low_complexity reads between two clock readings
constexpr int kSdustThreshold = 20;
constexpr size_t kSdustWindow = 64;
constexpr size_t kSdustPiece = 128;

/**
 * Whether sdust flags |s| anywhere with the seeder's parameters: the seeder's filter, repeated
 * here because that function is file-local to the aligner, which the pattern search does not
 * edit (§11). An optional diagnostic of a completed search (kNoteLowComplexity), run under the
 * deadline (one sdust over a 30,000-base repeat takes 0.7 s). sdust decides at each base from
 * the window of W bases ending there (its triplet counts, and the longest suffix whose counts
 * stay within T / 5), and it flags |s| iff some window holds a perfect interval. So |s| is read
 * in pieces of kSdustPiece + W - 1 bases overlapping by W - 1, which hold every window of |s|;
 * a window cut at a piece's start is a suffix of the window of |s| at that base and holds a
 * perfect interval only if that one does: some piece is flagged iff |s| is. Stops at the first
 * piece flagged, so a repeat costs one piece (its perfect intervals are what makes sdust slow:
 * about 6 ms for 191 bases of ATG), and reads the clock before every piece but the first:
 * nullopt when the work time passed first. The first piece is read whatever the clock says, so
 * that a pattern of at most kSdustPiece + W - 1 bases is always diagnosed (its answer never
 * depends on the machine for it), at the cost of one piece past the work time.
 */
std::optional<bool> is_low_complexity(std::string_view s, Budget &budget) {
    for (size_t begin = 0; begin < s.size(); begin += kSdustPiece) {
        if (begin && !budget.check_time())
            return std::nullopt;
        const std::string_view piece = s.substr(begin, kSdustPiece + kSdustWindow - 1);
        int n = 0;
        std::unique_ptr<uint64_t, decltype(std::free)*> r {
            sdust(0, reinterpret_cast<const uint8_t*>(piece.data()), piece.size(),
                  kSdustThreshold, kSdustWindow, &n),
            std::free
        };
        if (n > 0)
            return true;
        // the last piece reached the end of |s|
        if (begin + piece.size() == s.size())
            break;
    }
    return false;
}

std::vector<BaseSet> reverse_complement_sets(std::vector<BaseSet> q) {
    std::reverse(q.begin(), q.end());
    for (BaseSet &set : q) {
        set = complement_set(set);
    }
    return q;
}

// the orientations a request searches, in search order: a palindromic P once (§3)
std::vector<Orientation> searched_orientations(const Pattern &pattern, Strands strands) {
    if (pattern.is_palindromic())
        return { Orientation::PALINDROMIC };
    if (strands == Strands::FORWARD)
        return { Orientation::FORWARD };
    if (strands == Strands::REVERSE)
        return { Orientation::REVERSE };
    return { Orientation::FORWARD, Orientation::REVERSE };
}

/**
 * Whether a k-mer x equal to its reverse complement can contain the oriented pattern |q| at
 * offset |p|. Such an x exists only for even k and satisfies x[i] = complement(x[k-1-i]), so
 * where the pattern covers both i and k-1-i the two sets must admit complementary bases.
 * Skips the palindrome checks of a wrapped PRIMARY graph where no palindrome can match
 * (always at odd k, and for most exact patterns longer than k/2).
 */
bool may_be_palindromic(const std::vector<BaseSet> &q, size_t p, size_t k) {
    if (k % 2)
        return false;

    for (size_t i = 0; i < q.size(); ++i) {
        size_t mirror = k - 1 - (p + i);
        if (mirror >= p && mirror < p + q.size()
                && !(q[i] & complement_set(q[mirror - p])))
            return false;
    }
    return true;
}

/**
 * A range whose edges are contexts of one base search at one offset: a flank range (c == 0:
 * every valid edge leaving its nodes) or a W-rule leaf (the valid edges with W in
 * {c, c + alph_size}). Kept for the deferred scans and for enumerate()'s release.
 */
struct Item {
    edge_index first = 0;
    edge_index last = 0;
    uint32_t offset = 0;
    TAlphabet c = 0;
    // W rule: the candidate edges (plain + marked, R); flank: the valid edges (exact)
    uint64_t candidates = 0;
    // W rule: invalid edges whose W is not the sentinel (J): the most candidates can lose
    uint64_t invalid_ns = 0;
    // W rule: all invalid edges of the range (I), to choose the cheaper scan
    uint64_t invalid = 0;

    // a scan resolves what ranks cannot: which candidates are valid (a W rule with J > 0),
    // and on an even-k wrapped PRIMARY graph which contexts are palindromic k-mers
    enum class Scan : uint8_t {
        NONE,
        // the range's invalid edges through the mask (select0), testing W (§4.1)
        INVALID,
        // the candidates (W rule) or the valid edges (flank), testing validity and, with
        // count_palindromes, whether the k-mer is its own reverse complement
        CANDIDATES,
    };
    Scan scan = Scan::NONE;
    bool count_palindromes = false;
    bool scan_done = false;
    // a graph without the dummy-edge mask, a range whose nodes were not wholly spelled: its
    // candidates may include source dummies (k-mers starting with '$'), which nothing has
    // checked. Never set on a masked graph
    bool unchecked = false;
    uint64_t examined = 0;
    // INVALID: the examined invalid edges whose W is not the sentinel (s)
    uint64_t examined_ns = 0;
    // INVALID: the examined invalid edges with W in {c, c + alph_size} (t);
    // CANDIDATES on a W rule: the valid candidates found; on a graph without the mask, the
    // candidates found real (any range: the scan spells each)
    uint64_t hits = 0;
    uint64_t palindromes = 0;

    bool flank() const { return !c; }

    bool pending() const { return scan != Scan::NONE; }

    /**
     * The bounds of the range's valid contexts. §4.1: an interrupted scan of the invalid
     * edges leaves lower = max(0, R - t - (I - s)) and upper = R - t; I and s are counted
     * here over non-sentinel invalid edges (J), since an edge with W = $ never carries c:
     * the design's bound, tightened, never loosened. Unchecked candidates (no mask): the
     * real ones a scan found below, every candidate not found a dummy above.
     */
    uint64_t lower() const {
        if (unchecked)
            return scan == Scan::CANDIDATES ? hits : 0;
        if (flank() || scan == Scan::NONE)
            return candidates;

        if (scan == Scan::INVALID) {
            if (scan_done)
                return candidates - hits;
            uint64_t lost = hits + (invalid_ns - examined_ns);
            return candidates > lost ? candidates - lost : 0;
        }

        // CANDIDATES on a W rule: the unexamined candidates can still lose the invalid
        // non-sentinel edges not yet seen among the examined ones
        if (scan_done)
            return hits;
        uint64_t rest = candidates - examined;
        uint64_t seen_invalid = examined - hits;
        uint64_t rest_invalid = invalid_ns > seen_invalid ? invalid_ns - seen_invalid : 0;
        return hits + (rest > rest_invalid ? rest - rest_invalid : 0);
    }

    uint64_t upper() const {
        if (unchecked)
            return scan == Scan::CANDIDATES ? hits + (candidates - examined) : candidates;
        if (flank() || scan == Scan::NONE)
            return candidates;
        if (scan == Scan::INVALID)
            return candidates - hits;
        return hits + (candidates - examined);
    }

    // the count is known exactly: by ranks, or by a completed scan (never by bounds that
    // happen to meet); unchecked candidates only by a completed scan
    bool count_exact() const {
        if (unchecked)
            return scan == Scan::CANDIDATES && scan_done;
        return flank() || scan == Scan::NONE || scan_done;
    }

    uint64_t palindromes_upper() const {
        return count_palindromes ? palindromes + (candidates - examined) : 0;
    }
};

/**
 * A range kept for enumerate()'s release: only what the release reads (24 bytes; the scan
 * state of Item is not needed there), so that the retained ranges stay small.
 */
struct Span {
    edge_index first;
    edge_index last;
    // the contexts the release is sure to find in the range (Item::lower() when it was
    // counted, capped to 32 bits, still a lower bound): what PARTIAL's retention bound sums
    uint32_t guaranteed;
    uint16_t offset;
    // the W-rule symbol; 0 for a flank range
    TAlphabet c;
    // Item::unchecked: the release tests each candidate for a source dummy
    bool unchecked;

    bool flank() const { return !c; }
};
static_assert(sizeof(Span) == 24, "a retained range stays 24 bytes");

/**
 * A graph without the dummy-edge mask: an unchecked range kept for the check after discovery,
 * while the pattern's unchecked candidates number at most Request::max_checked_entries (so at
 * most that many of these). Its candidates: a flank's non-sink edges (c == 0), a W rule's edges
 * with W in {c, c + alph_size}.
 */
struct CheckRange {
    edge_index first;
    edge_index last;
    uint32_t offset;
    TAlphabet c;
};

/**
 * What is known of one count: started or not, its offset discovered completely or not,
 * every scan behind it finished or not, and its bounds.
 */
struct Estimate {
    bool started = false;
    bool discovered = false;
    bool exact = false;
    uint64_t lower = 0;
    uint64_t upper = 0;

    static Estimate zero() { return { true, true, true, 0, 0 }; }
};

Count to_count(const Estimate &e, Unit unit) {
    if (!e.started)
        return Count::unknown(unit);
    if (!e.discovered)
        return Count::at_least(unit, e.lower);
    if (!e.exact)
        return Count::bounds(unit, e.lower, e.upper);
    return Count::exact(unit, e.lower);
}

/**
 * A window [begin, end) of an oriented peptide's codon automaton (§6), the language one
 * BaseSearch matches: the codon sets of the residues it overlaps, the first and the last cut
 * to the codon positions inside the window as cylinders (the positions outside free), so that
 * two windows admitting the same strings compare equal. |phase| is the place of the window's
 * first position in its codon: 0 for a window starting an oriented pattern (position 0 is a
 * codon start), (L - k) % 3 for the reverse complement of a long peptide's anchor window
 * searched on a wrapped PRIMARY graph (rc(Q[0, k)) = rc(Q)[L - k, L), §4.1).
 */
struct CodonWindow {
    std::vector<CodonSet> sets;
    size_t phase = 0;

    bool operator==(const CodonWindow &other) const {
        return phase == other.phase && sets == other.sets;
    }

    // how many bases of window position j's codon before it lie inside the window: the ones
    // the automaton reads (the earlier ones are free in the cylinder)
    size_t known(size_t j) const { return std::min((phase + j) % 3, j); }

    // the bases allowed at window position |j| after the codon's bases b0 and b1 at codon
    // positions 0 and 1 (-1: not known, before the window)
    BaseSet allowed(size_t j, int b0, int b1) const {
        const size_t at = phase + j;
        return next_bases(sets[at / 3], at % 3, b0, b1);
    }

    /**
     * The window [begin, end) of the oriented peptide whose residue i has the codons
     * codons[i] (|reversed| false: P itself) or rc(codons[m - 1 - i]) (|reversed| true: rc(P)),
     * m = codons.size(). Builds only the residues the window overlaps (O(k)).
     */
    static CodonWindow of(const std::vector<CodonSet> &codons, bool reversed,
                          size_t begin, size_t end) {
        assert(begin < end && end <= 3 * codons.size());
        CodonWindow w;
        w.phase = begin % 3;
        for (size_t i = begin / 3; 3 * i < end; ++i) {
            const CodonSet set = reversed
                ? reverse_complement_codons(codons[codons.size() - 1 - i])
                : codons[i];
            const size_t first = std::max(begin, 3 * i) - 3 * i;
            const size_t last = std::min(end, 3 * i + 3) - 3 * i;
            w.sets.push_back(cylinder(set, first, last));
        }
        return w;
    }
};

/**
 * One search of one oriented pattern string q (|q| <= k) on the base BOSS, over the offsets
 * [0, k - |q|] (any_offset) or k - |q| only (suffix, and a long pattern's anchor window,
 * where |q| = k). §4.1: positions 0 .. |q|-2 on node ranges, the last on W, then the flank.
 * q is a list of base sets; for a peptide window (|automaton|) these are the per-position
 * unions, a superset used only for the cost estimate and the palindrome test, and the bases
 * tried at each depth are the automaton's after the bases the range's nodes end with
 * (PatternRun::symbols_at): a range of depth d is the set of nodes ending with the window's
 * first d spelled bases, so the spelled prefix travels with the range itself.
 *
 * A leading run of pattern N on a $ACGT graph is not searched (|lead|): every base of a valid
 * k-mer there is one of A, C, G, T, so q at offset p is exactly its core q[lead, |q|) at offset
 * p + lead, and the DFS runs on the core with its flank stopped where the offset of q reaches 0
 * (the run would otherwise be branched four ways per position on the widest ranges). Every
 * offset stored here is q's.
 */
struct BaseSearch {
    std::vector<BaseSet> q;
    // a peptide window: its codon automaton (§6); unset for DNA and IUPAC
    std::optional<CodonWindow> automaton;
    // the leading positions of q not searched (a pattern-N run, $ACGT graphs only; never in
    // a peptide window, whose first codon constrains the positions after the run)
    size_t lead = 0;
    // the BOSS codes of each core position's bases (q[lead + i]), in A, C, G, T order
    std::vector<std::vector<TAlphabet>> codes;
    bool any_offset = false;
    // count the palindromic contexts (even k, wrapped PRIMARY: what the union subtracts)
    bool count_palindromes = false;
    // some orientation reads this search's contexts directly (not only mapped to the
    // wrapper's reverse complements): its ranges' contexts start at their own first node
    bool read_directly = false;

    bool started = false;
    // discovery was interrupted (budget or threshold): every offset is then undiscovered
    // (AT_LEAST). v1 does not track which offsets were completely discovered before the stop
    bool halted = false;

    // per offset: the contexts known exactly when their range was discovered
    std::vector<uint64_t> exact;
    // per offset, a graph without the mask: the unchecked candidates of the ranges no scan
    // reads (Item::unchecked), counted into the upper bound only
    std::vector<uint64_t> unchecked;
    // a graph without the mask: those ranges, kept for the check of
    // Request::max_checked_entries while the pattern's unchecked candidates are within it
    std::vector<CheckRange> checks;
    // the ranges that need a scan, in discovery order
    std::vector<Item> pending;
    // enumerate(): the ranges with candidates the release may need, in discovery order
    std::vector<Span> release;
    // the lower bound of everything counted so far (stop_at_threshold, retention)
    uint64_t running_lower = 0;
    // the upper bound of everything counted so far: what the thresholds compare on a graph
    // without the mask
    uint64_t running_upper = 0;
    // the candidates of the items whose palindromes are counted (count_palindromes): at least
    // the palindromic contexts among them, which a wrapped PRIMARY union finds twice
    uint64_t running_palindrome_candidates = 0;

    // per offset, filled by tally() after the scans
    std::vector<Estimate> contexts;
    std::vector<Estimate> palindromes;

    uint32_t max_offset(size_t k) const { return static_cast<uint32_t>(k - q.size()); }
    // the positions the DFS matches
    size_t core_length() const { return q.size() - lead; }
    // the deepest range the DFS creates: a flank range there carries q's offset 0
    size_t max_depth(size_t k) const { return k - 1 - lead; }

    void tally() {
        // not started: nothing known, palindromes included
        contexts.assign(exact.size(), Estimate());
        palindromes.assign(exact.size(), Estimate());
        if (!started)
            return;

        const bool discovered = !halted;
        for (uint32_t p = 0; p < exact.size(); ++p) {
            contexts[p] = { true, discovered, discovered && !unchecked[p], exact[p],
                            exact[p] + unchecked[p] };
            palindromes[p] = { true, discovered, discovered, 0, 0 };
        }
        for (const Item &item : pending) {
            Estimate &e = contexts[item.offset];
            e.lower += item.lower();
            e.upper += item.upper();
            e.exact &= item.count_exact();
            if (item.count_palindromes) {
                Estimate &pal = palindromes[item.offset];
                pal.lower += item.palindromes;
                pal.upper += item.palindromes_upper();
                pal.exact &= item.scan_done;
            }
        }
    }
};

// how an orientation's contexts are made of base searches (§4.1, "Three graph modes")
struct OrientationPlan {
    Orientation orientation;
    // the base search of the oriented pattern (or of its anchor window) itself
    size_t direct;
    // wrapped PRIMARY only: the base search of the reverse complement, whose context
    // (y, p) is the wrapper's (y + offset, max_offset - p)
    std::optional<size_t> mapped;
    // wrapped PRIMARY, even k: the search counting the palindromic stored k-mers, which
    // both searches find and the union counts once, and whether its offsets are the
    // mirror of this orientation's (a pair already counting them is reused)
    std::optional<size_t> palindromes;
    bool palindromes_mirrored = false;
};

// one context in the release, ordered by (node, offset, orientation): §5.5
struct Key {
    node_index node = 0;
    uint32_t offset = 0;
    Orientation orientation = Orientation::FORWARD;

    bool operator<(const Key &other) const {
        return std::tie(node, offset, orientation)
                < std::tie(other.node, other.offset, other.orientation);
    }
    bool operator==(const Key &other) const {
        return std::tie(node, offset, orientation)
                == std::tie(other.node, other.offset, other.orientation);
    }
};

// the BaseSet bit of a spelled base; 0 for N and $, which no pattern position allows (§3)
BaseSet base_of_char(char c) {
    switch (c) {
        case 'A': return kBaseA;
        case 'C': return kBaseC;
        case 'G': return kBaseG;
        case 'T': return kBaseT;
        default: return 0;
    }
}

// the bases of an exact pattern (every position one base)
std::string exact_bases(const Pattern &pattern) {
    std::string bases;
    bases.reserve(pattern.length());
    for (BaseSet set : pattern.positions()) {
        bases.push_back(iupac_letter(set));
    }
    return bases;
}

std::string format_bits(double bits) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(1) << bits;
    return out.str();
}

/**
 * Counts, and for enumerate() retains and releases, the contexts of one pattern.
 * One object per call; all state is its own, so calls run concurrently.
 */
class PatternRun {
  public:
    PatternRun(const DeBruijnGraph &graph,
               const DBGSuccinct &dbg_succ,
               const GraphSupport &support,
               uint64_t wrapper_offset,
               const Pattern &pattern,
               const Request &request,
               Budget &budget,
               bool retain)
          : graph_(graph),
            canonical_(dynamic_cast<const CanonicalDBG*>(&graph)),
            dbg_succ_(dbg_succ),
            boss_(dbg_succ.get_boss()),
            support_(support),
            k_(support.k),
            wrapper_offset_(wrapper_offset),
            pattern_(pattern),
            request_(request),
            budget_(budget),
            retain_(retain),
            long_(pattern.length() > support.k),
            extending_(long_ && request.extend_paths),
            cap_(long_ ? request.max_anchors : request.max_contexts),
            // PARTIAL releases at most cap_ contexts: keep only the ranges that can hold one of
            // the first cap_ in answer order
            prune_(retain && !extending_ && request.mode == Mode::PARTIAL),
            masked_(support.mask_present),
            // the check of a few unchecked candidates: never with the mask
            check_limit_(support.mask_present ? 0 : request.max_checked_entries),
            checkable_(check_limit_ > 0) {
        for (char base : { 'A', 'C', 'G', 'T' }) {
            codes_.push_back(boss_.encode(base));
        }
        for (BaseSet set = 0; set <= kAllBases; ++set) {
            for (size_t b = 0; b < 4; ++b) {
                if (set & (1 << b))
                    codes_of_set_[set].push_back(codes_[b]);
            }
        }
        plan();
        // nothing to release: nothing to keep
        if (prune_ && !cap_)
            retain_ = false;
    }

    const std::vector<Orientation>& searched() const { return searched_; }
    const std::vector<OrientationPlan>& plans() const { return plans_; }
    const std::optional<Stop>& stop() const { return stop_; }
    bool time_limited() const { return time_limited_; }

    // kNoteThresholdUpperBound: a threshold decision went against the request on an upper
    // bound (no mask) that the lower bound did not cross
    bool upper_bound_decision() const { return upper_bound_decision_; }
    void note_upper_bound_decision() { upper_bound_decision_ = true; }

    // every discovered range with candidates is retained for the release (none dropped by
    // ALL_OR_COUNT's threshold, the extension's admission, or PARTIAL's retention bound): a
    // release that drains them enumerates every candidate
    bool retained_all() const {
        return retain_ && prune_bound_ == std::numeric_limits<node_index>::max();
    }
    Work work() const {
        Work work = work_;
        work.spans_retained_peak = retained_peak_;
        return work;
    }

    // the offsets of the pattern's scope, ascending; 0 for a long pattern's anchors
    uint32_t max_offset() const {
        return long_ ? 0 : static_cast<uint32_t>(k_ - pattern_.length());
    }

    std::vector<uint32_t> offsets() const {
        if (long_ || request_.scope == Scope::SUFFIX)
            return { max_offset() };
        std::vector<uint32_t> result;
        for (uint32_t p = 0; p <= max_offset(); ++p) {
            result.push_back(p);
        }
        return result;
    }

    /**
     * Discovery of every base search, then the deferred scans, so that a stop in a scan
     * leaves every range discovered and the counts BOUNDS (§4.1), then the tallies. The base
     * searches run in the order of their estimated cost (the cheaper first, the plan's order
     * on a tie: always for an exact pattern), so that a budget stop leaves the cheap ones
     * complete rather than not started.
     */
    void run() {
        steps_before_ = budget_.steps_used();

        if (budget_.stopped() || !budget_.check_time()) {
            // a request-wide stop before this pattern: nothing runs, every count UNKNOWN
            record_stop(StopPhase::DISCOVERY, *budget_.stopped());
        } else {
            std::vector<size_t> order(searches_.size());
            std::vector<double> cost(searches_.size());
            for (size_t i = 0; i < searches_.size(); ++i) {
                order[i] = i;
                cost[i] = estimated_cost(searches_[i]);
            }
            std::stable_sort(order.begin(), order.end(),
                             [&](size_t a, size_t b) { return cost[a] < cost[b]; });
            for (size_t i : order) {
                discover(searches_[i]);
                if (stop_)
                    break;
            }
            if (!stop_)
                scan_all(order);
            // a graph without the mask, every range discovered and scanned: a few unchecked
            // candidates are tested one by one
            if (!stop_)
                check_unchecked(order);
        }

        work_.steps = budget_.steps_used() - steps_before_;
        for (BaseSearch &search : searches_) {
            search.tally();
        }
    }

    // the contexts (or anchors) of one orientation at one offset of the served graph
    Count count(const OrientationPlan &o, uint32_t p, Unit unit) const {
        const Estimate &a = searches_[o.direct].contexts[p];
        if (!o.mapped)
            return to_count(a, unit);

        // wrapped PRIMARY: U = A + B - P, with B the reverse-complement search at the mirrored
        // offset and P the palindromic stored k-mers both find (§4.1: the probes are united
        // by wrapper node id before anything counts them)
        uint32_t mirror = max_offset() - p;
        const Estimate &b = searches_[*o.mapped].contexts[mirror];
        Estimate pal = Estimate::zero();
        if (o.palindromes) {
            uint32_t q = o.palindromes_mirrored ? mirror : p;
            const BaseSearch &source = searches_[*o.palindromes];
            if (may_be_palindromic(source.q, q, k_))
                pal = source.palindromes[q];
        }

        if (!a.started && !b.started)
            return Count::unknown(unit);

        if (!a.discovered || !b.discovered) {
            // an undiscovered part has no upper bound: AT_LEAST over what is known. The two
            // parts are disjoint when no palindrome can be among them; otherwise the union
            // is at least its larger part
            bool disjoint = pal.discovered && pal.exact && pal.upper == 0;
            return Count::at_least(unit, disjoint ? a.lower + b.lower
                                                  : std::max(a.lower, b.lower));
        }

        uint64_t sum_lower = a.lower + b.lower;
        uint64_t sum_upper = a.upper + b.upper;
        // the union is at least each part, and at least the sum less every possible duplicate
        uint64_t lower = std::max({ a.lower, b.lower,
                                    sum_lower > pal.upper ? sum_lower - pal.upper : 0 });
        uint64_t upper = sum_upper - std::min(pal.lower, sum_upper);
        if (a.exact && b.exact && pal.exact) {
            assert(lower == upper);
            return Count::exact(unit, upper);
        }
        return Count::bounds(unit, std::min(lower, upper), upper);
    }

    /**
     * The release (§4.3, label-free): the retained contexts of every orientation in
     * (node, offset, orientation) order, at most |limit| of them, each passed to |emit|.
     * Charges no steps; reads the clock before it starts, every kClockStride cursors it
     * prepares and every kClockStride edges it examines, and, when |emit| is the caller's
     * (|callers_emit|: a spelling and a result object per context, not a buffer's
     * push_back; also the listing of the anchors to extend, list_anchors), every
     * kReleaseClockStride contexts it emits. Returns false when the
     * deadline stopped it ({|phase|, TIME} recorded: EXTRACTION for a release, EXTENSION for
     * the listing of the anchors to extend). |exhausted|, when given and true is returned:
     * every retained range was drained (no context left in them; a cut at |limit| with
     * contexts left says false).
     */
    bool release(uint64_t limit, const std::function<void(const Context&)> &emit,
                 bool callers_emit, StopPhase phase = StopPhase::EXTRACTION,
                 bool *exhausted = nullptr);

    /**
     * ALL_OR_COUNT's delivery of its buffered release to |callback|, whose work per context
     * (a spelling, a result object) is the caller's: the clock is read before every
     * kReleaseClockStride-th context. False when the work time passed before the last one
     * ({EXTRACTION, TIME} recorded): |callback| then received a prefix, to be discarded.
     */
    bool deliver(const std::vector<Context> &buffer,
                 const std::function<void(const Context&)> &callback) {
        for (size_t i = 0; i < buffer.size(); ++i) {
            if (i && !(i % Budget::kReleaseClockStride) && !budget_.check_time()) {
                record_stop(StopPhase::EXTRACTION, StopReason::TIME);
                return false;
            }
            callback(buffer[i]);
        }
        return true;
    }

    // ---------------------------------------------------------------- extension (§4.2)

    /**
     * Phase 2 for an admitted anchor set (EXACT, <= max_anchors): lists the anchors in
     * answer order through release() (the anchor ranges are then discarded, §5.2), and
     * extends each by extend_anchor(). |keep|: retain complete paths for the release under
     * the rule of Request::max_paths (enumerate()); count() keeps none. |anchors| is their
     * EXACT count, which the listing must reproduce (std::logic_error otherwise: an anchor
     * dropped before its extension would make a path count look exact). True when every
     * anchor was extended to L; false when a stop ended it (recorded, phase EXTENSION).
     */
    bool extend(bool keep, uint64_t anchors);

    /**
     * The two halves of extend(): the listing of the anchors in answer order (the anchor ranges
     * then discarded; false when the deadline stopped it, {EXTENSION, TIME} recorded), and the
     * extension of a listed set (false when a stop ended it). Without the mask the listing
     * drops the source dummies, so that the caller learns the exact anchor count from it before
     * extending.
     */
    bool list_anchors(std::vector<Context> *anchors);
    bool extend_listed(bool keep, const std::vector<Context> &anchors);

    // the anchors' ranges are no longer needed: the extension's admission failed (§5.2)
    void drop_anchors() {
        retain_ = false;
        retained_ = 0;
        for (BaseSearch &search : searches_) {
            std::vector<Span>().swap(search.release);
        }
    }

    // the complete paths of anchors of |orientation| (all of them, retained or not)
    uint64_t paths_found(Orientation orientation) const {
        auto it = found_.find(orientation);
        return it == found_.end() ? 0 : it->second;
    }

    // every anchor of |orientation| was extended to L: its path count is exact
    bool orientation_extended(Orientation orientation) const {
        auto it = anchors_left_.find(orientation);
        return listed_ && (it == anchors_left_.end() || !it->second);
    }

    uint64_t candidates_examined() const { return candidates_; }

    // with a support tracker (Request::support): the supported paths of anchors of
    // |orientation|, the tracker's DEAD verdicts, and whether one came before L (in
    // |orientation|, or in any)
    uint64_t supported_found(Orientation orientation) const {
        auto it = supported_.find(orientation);
        return it == supported_.end() ? 0 : it->second;
    }
    uint64_t branches_pruned() const { return pruned_; }
    // the complete walks, supported or not (AnchorCounts::walks)
    uint64_t walks() const { return found_total_; }
    bool pruned_before_completion(Orientation orientation) const {
        return pruned_early_.count(orientation) > 0;
    }
    bool pruned_before_completion() const { return !pruned_early_.empty(); }

    /**
     * The retained paths, in answer order (the DFS's own), each passed to |emit|. Reads the
     * clock before the first and before every kReleaseClockStride-th (the caller's work per
     * path runs between two readings): false when the work time has passed ({EXTRACTION,
     * TIME} recorded), |emit| having received the paths before that reading.
     */
    bool release_paths(const std::function<void(const Context&)> &emit) {
        auto late = [&]() {
            if (budget_.check_time())
                return false;
            record_stop(StopPhase::EXTRACTION, StopReason::TIME);
            return true;
        };
        if (late())
            return false;
        for (size_t i = 0; i < paths_.size(); ++i) {
            if (i && !(i % Budget::kReleaseClockStride) && late())
                return false;
            emit(paths_[i]);
        }
        return true;
    }

    uint64_t paths_retained() const { return paths_.size(); }

  private:
    // the served graph: the DBGSuccinct itself, or the CanonicalDBG wrapping a PRIMARY one
    const DeBruijnGraph &graph_;
    const CanonicalDBG *canonical_;
    const DBGSuccinct &dbg_succ_;
    const BOSS &boss_;
    const GraphSupport &support_;
    const size_t k_;
    const uint64_t wrapper_offset_;
    const Pattern &pattern_;
    const Request &request_;
    Budget &budget_;
    // enumerate(): keep the ranges for the release; dropped in ALL_OR_COUNT once the
    // running lower bound shows that nothing will be released
    bool retain_;
    const bool long_;
    // L > k with Request::extend_paths: the anchors are retained for the extension in every
    // mode, and dropped once their count exceeds max_anchors (§5.2)
    const bool extending_;
    // ALL_OR_COUNT's threshold and PARTIAL's cap on what is released
    const uint64_t cap_;
    // PARTIAL: the retained ranges are bounded by what the release can use
    const bool prune_;
    // the graph has its dummy-edge mask; without it the candidates of a range not wholly
    // spelled are unchecked, the thresholds compare upper bounds, and the release tests every
    // unchecked candidate for a source dummy
    const bool masked_;
    // a graph without the mask: the most unchecked candidates the check after discovery tests
    // (Request::max_checked_entries; 0 with the mask or when disabled), the unchecked
    // candidates counted so far over every base search, and whether they are still within the
    // limit (their ranges kept in BaseSearch::checks)
    const uint64_t check_limit_;
    uint64_t unchecked_entries_ = 0;
    bool checkable_;
    // a threshold decision went against the request on an upper bound its lower bound did
    // not cross (kNoteThresholdUpperBound)
    bool upper_bound_decision_ = false;
    // the ranges retained now, and the most at once (Work::spans_retained_peak)
    uint64_t retained_ = 0;
    uint64_t retained_peak_ = 0;
    // PARTIAL: a range whose contexts all lie at nodes above this bound cannot hold one of the
    // first cap_ contexts in answer order, which lie at or below it
    node_index prune_bound_ = std::numeric_limits<node_index>::max();
    // PARTIAL: the next compaction when more ranges than this are retained
    uint64_t next_compaction_ = kMinCompaction;
    static constexpr uint64_t kMinCompaction = 4096;

    std::vector<TAlphabet> codes_;
    // the BOSS codes of every BaseSet's bases, in A, C, G, T order (a peptide window's
    // symbols at a range, symbols_at)
    std::array<std::vector<TAlphabet>, kAllBases + 1> codes_of_set_;
    std::vector<BaseSearch> searches_;
    std::vector<OrientationPlan> plans_;
    std::vector<Orientation> searched_;

    std::optional<Stop> stop_;
    bool time_limited_ = false;
    Work work_;
    uint64_t steps_before_ = 0;
    // the edges examined by the release, for the clock stride
    uint64_t examined_ = 0;
    // the phase a clock stop in release() is recorded with
    StopPhase release_phase_ = StopPhase::EXTRACTION;

    // the extension: the anchors listed, those not yet extended to L per orientation, the
    // complete paths per orientation and in all, the branches entered, the retained paths
    bool listed_ = false;
    std::map<Orientation, uint64_t> anchors_left_;
    std::map<Orientation, uint64_t> found_;
    uint64_t found_total_ = 0;
    uint64_t candidates_ = 0;
    // the nodes the DFS expanded, for its clock stride
    uint64_t expanded_ = 0;
    bool keep_paths_ = false;
    std::vector<Context> paths_;
    // with a support tracker: the supported paths per orientation and in all, the DEAD
    // verdicts, and the orientations with one before L
    std::map<Orientation, uint64_t> supported_;
    uint64_t supported_total_ = 0;
    uint64_t pruned_ = 0;
    std::map<Orientation, bool> pruned_early_;

    // the first stop is the pattern's; a later one only adds whether time touched it
    void record_stop(StopPhase phase, StopReason reason) {
        if (!stop_)
            stop_ = Stop { phase, reason };
        time_limited_ |= reason == StopReason::TIME;
    }

    // PatternSearch::base_node: the stored k-mer of a node of the served graph
    node_index base_of(node_index node) const {
        return support_.mode == GraphMode::PRIMARY && node > wrapper_offset_
            ? node - wrapper_offset_
            : node;
    }

    size_t add_search(const std::vector<BaseSet> &q, bool count_palindromes,
                      const std::optional<CodonWindow> &automaton = std::nullopt) {
        for (size_t i = 0; i < searches_.size(); ++i) {
            if (searches_[i].q == q && searches_[i].automaton == automaton) {
                searches_[i].count_palindromes |= count_palindromes;
                return i;
            }
        }
        BaseSearch search;
        search.q = q;
        search.automaton = automaton;
        search.any_offset = !long_ && request_.scope == Scope::ANY_OFFSET;
        search.count_palindromes = count_palindromes;
        // a leading pattern-N run is not searched where the graph has no N symbol: every
        // base of a valid k-mer is then one of A, C, G, T, which N admits (at least one
        // position is kept: an all-N window is searched as its last N). Never for a peptide
        // window: an X codon's first bases admit every base, but not every codon
        if (support_.alphabet == "$ACGT" && !automaton) {
            while (search.lead + 1 < q.size() && q[search.lead] == kAllBases) {
                ++search.lead;
            }
        }
        for (size_t i = search.lead; i < q.size(); ++i) {
            search.codes.emplace_back();
            for (size_t b = 0; b < 4; ++b) {
                if (q[i] & (1 << b))
                    search.codes.back().push_back(codes_[b]);
            }
        }
        search.exact.assign(search.max_offset(k_) + 1, 0);
        search.unchecked.assign(search.max_offset(k_) + 1, 0);
        searches_.push_back(std::move(search));
        return searches_.size() - 1;
    }

    /**
     * The orientations and the base searches behind them. A palindromic P is searched once
     * (§3). BASIC and native CANONICAL: one base search per orientation. Wrapped PRIMARY
     * (§4.1): the wrapper's contexts of Q are the stored k-mers containing Q and the virtual
     * reverse complements of the stored ones containing rc(Q); for a long pattern, the
     * anchor window Q[0, k) and rc(Q[0, k)) (not rc(Q)[0, k), which is Q's last k-mer).
     * Identical base searches are run once (P and rc(P) serve both orientations).
     * Only the windows are built (O(k)), never rc of the whole pattern: REVERSE's window
     * rc(P)[0, k) is rc(P[L - k, L)). A peptide's windows carry their codon automaton (§6):
     * Q[0, m) of the oriented pattern, and on a wrapped PRIMARY graph rc(Q[0, m)), the window
     * [L - m, L) of the other orientation's automaton.
     */
    void plan() {
        searched_ = searched_orientations(pattern_, request_.strands);

        const bool primary = support_.mode == GraphMode::PRIMARY;
        const bool even_primary = primary && !(k_ % 2);
        const auto &positions = pattern_.positions();
        const size_t L = positions.size();
        const size_t m = std::min(L, k_);
        const bool protein = pattern_.kind() == PatternKind::PROTEIN;

        for (Orientation orientation : searched_) {
            std::vector<BaseSet> window = orientation == Orientation::REVERSE
                ? reverse_complement_sets(std::vector<BaseSet>(positions.end() - m,
                                                               positions.end()))
                : std::vector<BaseSet>(positions.begin(), positions.begin() + m);
            const bool reverse = orientation == Orientation::REVERSE;
            std::optional<CodonWindow> automaton;
            if (protein)
                automaton = CodonWindow::of(pattern_.codon_sets(), reverse, 0, m);

            OrientationPlan o { orientation, 0, std::nullopt, std::nullopt, false };
            if (!primary) {
                o.direct = add_search(window, false, automaton);
                searches_[o.direct].read_directly = true;
                plans_.push_back(o);
                continue;
            }

            std::vector<BaseSet> window_rc = reverse_complement_sets(window);
            // rc(Q[0, m)): the positions [L - m, L) of the other orientation's automaton
            std::optional<CodonWindow> automaton_rc;
            if (protein)
                automaton_rc = CodonWindow::of(pattern_.codon_sets(), !reverse, L - m, L);
            // the palindromic k-mers at offset p of the window's search are those at the
            // mirrored offset of its reverse complement's: count them on one of the two
            bool mirrored = false;
            if (even_primary) {
                for (const BaseSearch &search : searches_) {
                    mirrored |= search.q == window_rc && search.automaton == automaton_rc
                                    && search.count_palindromes;
                }
            }
            o.direct = add_search(window, even_primary && !mirrored, automaton);
            searches_[o.direct].read_directly = true;
            o.mapped = add_search(window_rc, false, automaton_rc);
            if (even_primary) {
                o.palindromes = mirrored ? *o.mapped : o.direct;
                o.palindromes_mirrored = mirrored;
            }
            plans_.push_back(o);
        }
    }

    // ---------------------------------------------------------------- discovery

    // the running lower bound of the whole pattern (§5.2 stop_at_threshold): the sum over
    // orientations. A wrapped PRIMARY union of a + b discovered contexts counts both parts when
    // k is odd (no k-mer of odd length over A, C, G, T is its own reverse complement); for even
    // k it can count a palindromic stored k-mer twice, at most once per candidate of the ranges
    // where the counting search may meet one (pc), so it holds at least a + b - min(pc, a, b)
    // (max(a, b) would be about half the count)
    uint64_t running_lower() const {
        uint64_t total = 0;
        for (const OrientationPlan &o : plans_) {
            uint64_t a = searches_[o.direct].running_lower;
            if (!o.mapped) {
                total += a;
            } else {
                uint64_t b = searches_[*o.mapped].running_lower;
                uint64_t twice = 0;
                if (o.palindromes) {
                    twice = std::min({ searches_[*o.palindromes].running_palindrome_candidates,
                                       a, b });
                }
                total += a + b - twice;
            }
        }
        return total;
    }

    // the running upper bound of the whole pattern, as the thresholds compare it without the
    // mask: every candidate counted so far, of both parts of a wrapped PRIMARY union, but for
    // the ranges whose palindrome scan is pending (add_item): at most the final U, which it
    // equals on BASIC, CANONICAL and odd-k graphs
    uint64_t running_upper() const {
        uint64_t total = 0;
        for (const OrientationPlan &o : plans_) {
            total += searches_[o.direct].running_upper;
            if (o.mapped)
                total += searches_[*o.mapped].running_upper;
        }
        return total;
    }

    // after every count: the release's retention (ALL_OR_COUNT releases nothing once the count
    // exceeds max_contexts, §5.2) and the threshold stop of stop_at_threshold. The count
    // compared is the running lower bound with the mask, and the running upper bound without it
    // (conservative, as the admissions after discovery)
    bool threshold_crossed() {
        const uint64_t threshold = long_ ? request_.max_anchors : request_.max_contexts;
        const uint64_t compared = masked_ ? running_lower() : running_upper();
        // ALL_OR_COUNT releases nothing above its threshold, and the extension is not admitted
        // above max_anchors in any mode (§5.2: anchors kept through the admission). Without the
        // mask, while the unchecked candidates are few enough for the check after discovery,
        // the count may still become EXACT and within the threshold: the ranges are dropped
        // only once the lower bound is above it
        const uint64_t retained_on = !masked_ && checkable_ ? running_lower() : compared;
        if (retain_ && (request_.mode == Mode::ALL_OR_COUNT || extending_)
                && retained_on > threshold) {
            drop_anchors();
        }
        if (request_.stop_at_threshold && compared > threshold) {
            if (!masked_ && running_lower() <= threshold)
                upper_bound_decision_ = true;
            record_stop(StopPhase::DISCOVERY,
                        long_ ? StopReason::MAX_ANCHORS : StopReason::MAX_CONTEXTS);
            return true;
        }
        return false;
    }

    // one range evaluation (a tighten_range, or the rank set of the W rule): one step
    bool charge_range() {
        if (!budget_.charge(1)) {
            record_stop(StopPhase::DISCOVERY, *budget_.stopped());
            return false;
        }
        ++work_.ranges_visited;
        return true;
    }

    static void halt(BaseSearch &search) {
        search.halted = true;
    }

    // whether the release may read a range of |search| at |offset| through a cursor whose
    // contexts start at its own first node (the direct read, or the palindromes of a mapped
    // one), rather than only shifted by the wrapper's offset
    bool read_at_own_node(const BaseSearch &search, uint32_t offset) const {
        return search.read_directly || may_be_palindromic(search.q, offset, k_);
    }

    void add_item(BaseSearch &search, const Item &item) {
        if (retain_) {
            const uint64_t lower = item.lower();
            const node_index lowest = read_at_own_node(search, item.offset)
                ? item.first
                : item.first + wrapper_offset_;
            // PARTIAL: a range wholly above the bound holds none of the first cap_ contexts
            if (!prune_ || lowest <= prune_bound_) {
                assert(item.offset <= std::numeric_limits<uint16_t>::max());
                search.release.push_back(Span {
                    item.first, item.last,
                    static_cast<uint32_t>(std::min<uint64_t>(
                            lower, std::numeric_limits<uint32_t>::max())),
                    static_cast<uint16_t>(item.offset), item.c, item.unchecked
                });
                retained_peak_ = std::max(retained_peak_, ++retained_);
                if (prune_ && retained_ > next_compaction_)
                    compact_release();
            }
        }
        if (item.pending()) {
            search.pending.push_back(item);
        } else if (item.unchecked) {
            search.unchecked[item.offset] += item.candidates;
            keep_for_check(search, item);
        } else {
            search.exact[item.offset] += item.candidates;
        }
        search.running_lower += item.lower();
        // a range whose scan is pending (the palindromes of an even-k wrapped PRIMARY graph)
        // can only lose candidates to it, and its palindromes come off the union: counted
        // when scanned, so that the running count never exceeds the final U and a retention
        // dropped or a stop taken on it agrees with the admission after discovery. Without
        // a pending scan (every range on BASIC and CANONICAL graphs) it is the final U
        if (!item.pending())
            search.running_upper += item.upper();
        if (item.count_palindromes)
            search.running_palindrome_candidates += item.candidates;
    }

    /**
     * A graph without the mask: an unchecked range no scan reads, kept for the check after
     * discovery while the pattern's unchecked candidates number at most check_limit_; once they
     * are more, nothing will be checked and every kept range is freed.
     */
    void keep_for_check(BaseSearch &search, const Item &item) {
        if (!checkable_)
            return;
        unchecked_entries_ += item.candidates;
        if (unchecked_entries_ <= check_limit_) {
            search.checks.push_back(CheckRange { item.first, item.last, item.offset, item.c });
            return;
        }
        checkable_ = false;
        for (BaseSearch &s : searches_) {
            std::vector<CheckRange>().swap(s.checks);
        }
    }

    /**
     * PARTIAL's retention bound: every cursor the release would build over the retained ranges,
     * as (lowest node, highest node, contexts it is sure to emit), where distinct cursors never
     * emit one context twice except a palindrome cursor, which is credited none. If the cursors
     * that end at or before a node T are sure to emit cap_ contexts, the first cap_ contexts in
     * answer order (node first) lie at or below T, so a range none of whose cursors starts at
     * or below T can be dropped. T found over the ranges retained so far can only fall as more
     * arrive, so the bound holds for every later range too (add_item drops those at once).
     */
    void compact_release() {
        struct Bound {
            node_index high;
            uint64_t sure;
        };
        std::vector<Bound> bounds;
        bounds.reserve(retained_ + retained_ / 2);
        for_each_use([&](const OrientationPlan &, const Span &span, Use use) {
            bounds.push_back(Bound { use_high(span, use), use_sure(span, use) });
        });
        std::sort(bounds.begin(), bounds.end(),
                  [](const Bound &a, const Bound &b) { return a.high < b.high; });
        uint64_t sure = 0;
        for (const Bound &bound : bounds) {
            sure += bound.sure;
            if (sure >= cap_) {
                prune_bound_ = std::min(prune_bound_, bound.high);
                break;
            }
        }

        retained_ = 0;
        for (size_t s = 0; s < searches_.size(); ++s) {
            BaseSearch &search = searches_[s];
            auto end = std::remove_if(search.release.begin(), search.release.end(),
                                      [&](const Span &span) {
                const node_index lowest = read_at_own_node(search, span.offset)
                    ? span.first
                    : span.first + wrapper_offset_;
                return lowest > prune_bound_;
            });
            search.release.erase(end, search.release.end());
            if (search.release.capacity() > 2 * search.release.size() + kMinCompaction)
                search.release.shrink_to_fit();
            retained_ += search.release.size();
        }
        next_compaction_ = std::max(kMinCompaction, 2 * retained_);
    }

    /**
     * How the release reads a retained range (one cursor each):
     *  DIRECT           the contexts of the orientation's own base search, at their nodes;
     *  MAPPED           wrapped PRIMARY: the reverse-complement search's, at the wrapper's
     *                   virtual node (node + wrapper offset) and the mirrored offset;
     *  MAPPED_SKIPPING  the same at an offset where a palindromic k-mer is possible (even k):
     *                   palindromic k-mers are skipped (CanonicalDBG serves them once, at
     *                   their own id);
     *  PALINDROMES      beside MAPPED_SKIPPING: only those palindromic k-mers, at their own id
     *                   and the mirrored offset, which the direct search may not have reached
     *                   before a stop; the merge drops the duplicate when it has.
     */
    enum class Use : uint8_t { DIRECT, MAPPED, MAPPED_SKIPPING, PALINDROMES };

    node_index use_low(const Span &span, Use use) const {
        return use == Use::MAPPED || use == Use::MAPPED_SKIPPING
            ? span.first + wrapper_offset_
            : span.first;
    }
    node_index use_high(const Span &span, Use use) const {
        return use == Use::MAPPED || use == Use::MAPPED_SKIPPING
            ? span.last + wrapper_offset_
            : span.last;
    }
    // the contexts a cursor is sure to emit that no other cursor emits
    static uint64_t use_sure(const Span &span, Use use) {
        return use == Use::DIRECT || use == Use::MAPPED ? span.guaranteed : 0;
    }

    /**
     * Calls f(plan, span, use) for every cursor the release builds, in the order it builds
     * them: per orientation, the direct search's ranges, then (wrapped PRIMARY) the mapped
     * search's, each mapped one at a palindrome-capable offset followed by its PALINDROMES.
     */
    template <class F>
    void for_each_use(F f) const {
        for (const OrientationPlan &o : plans_) {
            for (const Span &span : searches_[o.direct].release) {
                f(o, span, Use::DIRECT);
            }
            if (!o.mapped)
                continue;
            const BaseSearch &mapped = searches_[*o.mapped];
            for (const Span &span : mapped.release) {
                if (may_be_palindromic(mapped.q, span.offset, k_)) {
                    f(o, span, Use::MAPPED_SKIPPING);
                    f(o, span, Use::PALINDROMES);
                } else {
                    f(o, span, Use::MAPPED);
                }
            }
        }
    }

    // a flank range at depth d >= L: every valid edge leaving its nodes has the pattern at
    // offset k - 1 - d (§4.1, all offsets; the core's offset less the lead for q's)
    void count_flank(BaseSearch &search, const Range &range) {
        const auto &[first, last, depth] = range;
        Item item;
        item.first = first;
        item.last = last;
        item.offset = static_cast<uint32_t>(k_ - 1 - depth - search.lead);
        if (masked_) {
            item.candidates = dbg_succ_.count_valid_edges_in_range(first, last);
        } else {
            // no mask: every edge that can carry a base, the source dummies among them unless
            // the nodes are wholly spelled (depth k - 1)
            item.candidates = dbg_succ_.count_non_sink_edges_in_range(first, last);
            item.unchecked = depth < boss_.get_k();
        }
        if (!item.candidates)
            return;
        if (search.count_palindromes && may_be_palindromic(search.q, item.offset, k_)) {
            item.scan = Item::Scan::CANDIDATES;
            item.count_palindromes = true;
        }
        add_item(search, item);
    }

    // the W rule at a leaf (nodes ending with q[0, L-1)): the k-mers ending with q[0, L)
    // are the valid edges with W in {c, c + alph_size} (§4.1, "Counting")
    void count_w_rule(BaseSearch &search, const Range &range, TAlphabet c) {
        const auto &[first, last, depth] = range;
        if (!masked_) {
            // no mask: the candidates, unchecked unless the nodes are wholly spelled (depth
            // k - 1: no room for '$'); no INVALID scan exists
            Item item;
            item.first = first;
            item.last = last;
            item.offset = search.max_offset(k_);
            item.c = c;
            item.candidates = dbg_succ_.count_edges_with_symbol(first, last, c);
            if (!item.candidates)
                return;
            item.unchecked = depth < boss_.get_k();
            if (search.count_palindromes && may_be_palindromic(search.q, item.offset, k_)) {
                item.scan = Item::Scan::CANDIDATES;
                item.count_palindromes = true;
            }
            add_item(search, item);
            return;
        }
        auto edges = dbg_succ_.count_edges_with_last_symbol(first, last, c);
        if (!edges.candidates)
            return;

        Item item;
        item.first = first;
        item.last = last;
        item.offset = search.max_offset(k_);
        item.c = c;
        item.candidates = edges.candidates;
        item.invalid_ns = edges.invalid_non_sentinel;
        item.invalid = edges.invalid;
        if (search.count_palindromes && may_be_palindromic(search.q, item.offset, k_)) {
            item.scan = Item::Scan::CANDIDATES;
            item.count_palindromes = true;
        } else if (item.invalid_ns) {
            // whichever is shorter: the range's invalid edges (the design's select0 scan)
            // or the candidates themselves; either resolves the count exactly
            item.scan = item.candidates < item.invalid ? Item::Scan::CANDIDATES
                                                       : Item::Scan::INVALID;
        }
        add_item(search, item);
    }

    /**
     * Expands one range at its depth d (the length of the suffix its nodes share), L being the
     * core's length (the window less its unsearched lead):
     *  - d >= L: a flank range, counted;
     *  - d == L - 1: a leaf: the W rule for each base of the last position, then (any_offset)
     *    the ranges of the whole core;
     *  - children: the core's bases at d < L, every non-sentinel symbol in the flank, down to
     *    the search's max_depth (k - 1 less the lead, where the window's offset reaches 0).
     * A child at depth 1 starts a suffix_to_prefix DFS; a child at depth k - 1 (whole nodes)
     * is expanded here, since suffix_to_prefix would call its edges one by one; the others
     * go through |push|, the try_symbol of the suffix_to_prefix DFS that popped |range|.
     * After a halt the ranges still created are counted (they were charged) but not expanded.
     */
    void expand(BaseSearch &search, const Range &range,
                const std::function<void(TAlphabet)> *push) {
        const size_t d = std::get<2>(range);
        const size_t L = search.core_length();
        const size_t max_depth = search.max_depth(k_);

        if (d >= L) {
            count_flank(search, range);
            if (threshold_crossed())
                halt(search);
        }
        if (search.halted)
            return;

        if (d + 1 == L) {
            const std::vector<TAlphabet> &symbols = symbols_at(search, range, L - 1);
            for (TAlphabet c : symbols) {
                if (!charge_range()) {
                    halt(search);
                    return;
                }
                count_w_rule(search, range, c);
                if (threshold_crossed()) {
                    halt(search);
                    return;
                }
            }
            if (search.any_offset && d < max_depth)
                children(search, range, symbols, push);

        } else if (d + 1 < L) {
            children(search, range, symbols_at(search, range, d), push);

        } else if (d < max_depth) {
            // the flank admits every symbol of the graph's alphabet but $, N on a DNA5 build
            // included (§4.1): the seeder's own symbol set
            std::vector<TAlphabet> symbols;
            align::NonSentinelSymbols()(boss_, range,
                                        [&](TAlphabet s) { symbols.push_back(s); });
            children(search, range, symbols, push);
        }
    }

    /**
     * The BOSS codes of the bases |search|'s window allows at its core position |x| for the
     * nodes of |range|, whose depth is x (they end with the window's first x bases): for DNA
     * and IUPAC the position's own (search.codes[x], whatever the range); for a peptide
     * window the automaton's after the codon's bases already spelled (§6), read from the
     * range's nodes — the last one by get_node_last_value, the one before through bwd (both
     * O(1) BOSS reads per range evaluated, which is charged its step) — in A, C, G, T order.
     */
    const std::vector<TAlphabet>& symbols_at(const BaseSearch &search, const Range &range,
                                             size_t x) const {
        if (!search.automaton)
            return search.codes[x];

        const CodonWindow &window = *search.automaton;
        assert(!search.lead && std::get<2>(range) == x);
        const size_t j = (window.phase + x) % 3;
        int known[2] = { -1, -1 };
        edge_index e = std::get<0>(range);
        for (size_t i = 0; i < window.known(x); ++i) {
            if (i)
                e = boss_.bwd(e);
            const TAlphabet s = boss_.get_node_last_value(e);
            const auto it = std::find(codes_.begin(), codes_.end(), s);
            // the window admits A, C, G, T only, so its spelled bases are among them
            assert(it != codes_.end());
            if (it == codes_.end())
                return codes_of_set_[0];
            known[j - 1 - i] = static_cast<int>(it - codes_.begin());
        }
        return codes_of_set_[window.allowed(x, known[0], known[1])];
    }

    void children(BaseSearch &search, const Range &range,
                  const std::vector<TAlphabet> &symbols,
                  const std::function<void(TAlphabet)> *push) {
        const size_t d = std::get<2>(range);
        const size_t K1 = boss_.get_k();

        for (TAlphabet s : symbols) {
            if (!charge_range()) {
                halt(search);
                return;
            }
            if (d && d + 1 < K1) {
                assert(push);
                (*push)(s);
                continue;
            }

            Range child = range;
            auto &[first, last, depth] = child;
            ++depth;
            if (boss_.tighten_range(&first, &last, s)) {
                if (depth == K1) {
                    expand(search, child, nullptr);
                } else {
                    dfs(search, child);
                }
            }
            if (search.halted)
                return;
        }
    }

    // the seeder's range DFS (suffix_to_prefix) from a range at depth 1 <= d < k - 1, with
    // the symbols chosen per depth by expand() instead of every symbol at every depth
    void dfs(BaseSearch &search, const Range &start) {
        assert(std::get<2>(start) >= 1 && std::get<2>(start) < boss_.get_k());
        align::suffix_to_prefix(
            dbg_succ_, start,
            [](node_index) {
                // never reached: whole-node ranges are expanded by expand() itself
                assert(false);
            },
            [&](const BOSS&, const Range &incremented, const auto &try_symbol) {
                // suffix_to_prefix hands over the popped range with its length already
                // incremented to that of its children
                Range range = incremented;
                --std::get<2>(range);
                std::function<void(TAlphabet)> push = [&](TAlphabet s) { try_symbol(s); };
                expand(search, range, &push);
            }
        );
    }

    /**
     * The ranges a base search's DFS can create over its core, level by level: the core's
     * combinations so far, of which at most the graph's share (edges / 4^d) can be non-empty. A
     * planning estimate (an N run early in a window is wide, late in it narrow), used only to
     * order the searches; the flank is the same for every search.
     */
    double estimated_cost(const BaseSearch &search) const {
        const double edges = static_cast<double>(boss_.num_edges());
        double combinations = 1;
        double strings = 1;
        double cost = 0;
        for (const auto &symbols : search.codes) {
            combinations *= symbols.size();
            strings *= 4;
            cost += combinations * std::min(1.0, edges / strings);
        }
        return cost;
    }

    void discover(BaseSearch &search) {
        if (search.started)
            return;
        search.started = true;
        // the empty suffix: every edge, the ranges of whole node groups from the first
        Range all { 1, boss_.num_edges(), 0 };
        expand(search, all, nullptr);
    }

    // ---------------------------------------------------------------- extension

    // one outgoing edge examined by the extension: one step
    bool charge_edge() {
        if (!budget_.charge(1)) {
            record_stop(StopPhase::EXTENSION, *budget_.stopped());
            return false;
        }
        ++work_.extension_edges;
        return true;
    }

    // the outgoing k-mers of |node|, whose sequence is the last k bases of |spelled|; on the
    // wrapper with that sequence as the spelling hint (it is the node's sequence, read once
    // for the anchor and then extended base by base)
    void call_outgoing(node_index node, const std::string &spelled,
                       const DeBruijnGraph::OutgoingEdgeCallback &callback) const {
        if (canonical_) {
            canonical_->call_outgoing_kmers(node, spelled.substr(spelled.size() - k_), callback);
        } else {
            graph_.call_outgoing_kmers(node, callback);
        }
    }

    // a STOPPED verdict of the support tracker or a false accept() of the path sink: TIME when
    // they found the work time passed through the request's Budget (check_time, which then
    // holds TIME: nothing else stops it while the extension runs), else EXTERNAL, whose reason
    // the route reads from them
    void external_stop() {
        record_stop(StopPhase::EXTENSION, budget_.stopped().value_or(StopReason::EXTERNAL));
    }

    // a DEAD verdict of the support tracker: the branch of |orientation| is pruned; before
    // L, the walks below it are not counted
    void prune(Orientation orientation, bool before_completion) {
        ++pruned_;
        if (before_completion)
            pruned_early_[orientation] = true;
    }

    /**
     * A complete walk of |anchor| (the Model accepting at L): counted; with a support tracker
     * (Request::support) supported or not by its complete() (a DEAD answer prunes it, a
     * complete walk all the same); then a path: handed to the sink (Request::sink) or kept
     * under the retention rule of Request::max_paths (§5.2), and checked against max_paths
     * for stop_at_threshold — under a tracker both count the supported paths. False when a
     * stop ended the extension (the threshold, the tracker's or the sink's).
     */
    bool complete_path(const Context &anchor, const std::vector<node_index> &path,
                       const std::string &spelled) {
        ++found_[anchor.orientation];
        ++found_total_;
        SupportTracker *const tracker = request_.support;
        const PathView view { anchor, path, spelled };
        if (tracker) {
            switch (tracker->complete(view)) {
                case SupportTracker::Verdict::ALIVE:
                    break;
                case SupportTracker::Verdict::DEAD:
                    prune(anchor.orientation, false);
                    return true;
                case SupportTracker::Verdict::STOPPED:
                    external_stop();
                    return false;
            }
            ++supported_[anchor.orientation];
            ++supported_total_;
        }
        // the paths the rules of max_paths read: the supported ones under a tracker
        const uint64_t released = tracker ? supported_total_ : found_total_;
        if (request_.sink) {
            if (!request_.sink->accept(view, tracker)) {
                external_stop();
                return false;
            }
        } else if (keep_paths_) {
            const bool keep = request_.mode == Mode::PARTIAL
                ? paths_.size() < request_.max_paths
                : released <= request_.max_paths;
            if (keep) {
                paths_.push_back(Context { anchor.orientation, 0, anchor.node, anchor.base_node,
                                           path, spelled });
            } else if (request_.mode != Mode::PARTIAL) {
                // ALL_OR_COUNT releases all paths or none: none, once they are too many
                keep_paths_ = false;
                std::vector<Context>().swap(paths_);
            }
        }
        if (request_.stop_at_threshold && released > request_.max_paths) {
            record_stop(StopPhase::EXTENSION, StopReason::MAX_PATHS);
            return false;
        }
        return true;
    }

    /**
     * The depth-first search of one anchor (§4.2) over search states (SearchState: the node,
     * the position, the Model's state, the tracker's frame): from the anchor's k spelled bases,
     * at every position the outgoing k-mers of the path's last node whose base the |model|
     * allows in its state there (Model::bases; the state after each base entered, Model::next),
     * in symbol order (A, C, G, T), to position L = model.length() of the oriented pattern. One
     * step per outgoing edge examined, allowed or not; the clock before every
     * kReleaseClockStride-th node expanded, counted over the whole extension. With a support
     * tracker (Request::support) one frame per level: opened at the anchor, pushed when a child
     * is entered (before it is expanded or completed), popped when the DFS leaves it; a DEAD
     * verdict prunes the branch (not expanded), a STOPPED one ends the extension. Without one,
     * the plain DFS of long_search "paths". False when a stop ended it.
     * A template on the Model's type: with a Pattern (final, its Model methods defined in this
     * file) every call is resolved and inlined, so that the plain DFS costs no more than one
     * written for Pattern alone (per extension edge, within the noise of an A/A timing);
     * another Model instantiates it as itself, or as Model through the virtual calls.
     */
    template <class M>
    bool extend_anchor(const Context &anchor, const M &model) {
        static_assert(std::is_base_of_v<Model, M>, "the extension's automaton is a Model");
        const size_t L = model.length();
        assert(L > k_);

        std::string spelled = graph_.get_node_sequence(anchor.node);
        assert(spelled.size() == k_);
        std::vector<node_index> path { anchor.node };

        // the allowed outgoing k-mers of one node of the path, in symbol order, the next one
        // to enter, and the Model's state at that node (after its bases). At most four: a node
        // has one outgoing k-mer per last base, and a pattern position allows only A, C, G, T
        // (never N or $)
        struct Level {
            std::array<std::pair<char, node_index>, 4> children;
            Model::State state = 0;
            uint8_t size = 0;
            uint8_t next = 0;
        };
        std::vector<Level> levels;
        levels.reserve(L - k_);

        // the support tracker's frames open, popped on every way out (SupportTracker's
        // protocol: each ALIVE open() and push() popped exactly once)
        SupportTracker *const tracker = request_.support;
        struct Frames {
            SupportTracker *tracker;
            uint32_t open = 0;
            void pop() {
                assert(open);
                tracker->pop();
                --open;
            }
            ~Frames() {
                while (open) {
                    pop();
                }
            }
        } frames { tracker };
        // the search state of the node just entered (built only for a tracker)
        auto state_of = [&](node_index node, char base, Model::State s) {
            const node_index stored = base_of(node);
            return SearchState { node, stored, stored != node, anchor.orientation, Side::RIGHT,
                                 base, static_cast<uint32_t>(spelled.size()), s,
                                 static_cast<uint32_t>(path.size() - 1), spelled };
        };

        auto expand = [&](Model::State state) {
            // the clock before every kReleaseClockStride-th node expanded: a node's outgoing
            // k-mers cost a few BOSS steps each (more on the wrapper of a PRIMARY graph) and
            // one step is charged per edge, so that a stride of steps can take long on a cold
            // index
            if (!(++expanded_ % Budget::kReleaseClockStride) && !budget_.check_time()) {
                record_stop(StopPhase::EXTENSION, StopReason::TIME);
                return false;
            }
            const BaseSet allowed = model.bases(state, static_cast<uint32_t>(spelled.size()));
            Level level;
            level.state = state;
            bool charged = true;
            call_outgoing(path.back(), spelled, [&](node_index next, char c) {
                if (!charged || !(charged = charge_edge()))
                    return;
                if (!(allowed & base_of_char(c)))
                    return;
                assert(level.size < level.children.size());
                if (level.size < level.children.size())
                    level.children[level.size++] = { c, next };
            });
            if (!charged)
                return false;
            if (level.size > 1)
                ++work_.extension_branches;
            // symbol order, by an insertion sort of at most four (std::sort's path for more
            // than 16 elements makes GCC 13 -O3 report -Warray-bounds on this array)
            for (uint8_t i = 1; i < level.size; ++i) {
                for (uint8_t j = i;
                        j && level.children[j].first < level.children[j - 1].first; --j) {
                    std::swap(level.children[j], level.children[j - 1]);
                }
            }
            levels.push_back(level);
            return true;
        };

        const Model::State start = model.start(spelled);
        if (tracker) {
            // frame 0: the anchor's support, before its first expansion
            switch (tracker->open(state_of(anchor.node, '\0', start))) {
                case SupportTracker::Verdict::ALIVE:
                    ++frames.open;
                    break;
                case SupportTracker::Verdict::DEAD:
                    // nothing supports the anchor: none of its walks is followed
                    prune(anchor.orientation, true);
                    return true;
                case SupportTracker::Verdict::STOPPED:
                    external_stop();
                    return false;
            }
        }

        if (!expand(start))
            return false;

        while (levels.size()) {
            Level &top = levels.back();
            if (top.next == top.size) {
                levels.pop_back();
                // the anchor itself is never popped (its frame 0 after the loop)
                if (levels.size()) {
                    path.pop_back();
                    spelled.pop_back();
                    if (tracker)
                        frames.pop();
                }
                continue;
            }
            const auto [c, next] = top.children[top.next++];
            const Model::State state = model.next(top.state, static_cast<uint32_t>(spelled.size()),
                                                  c);
            path.push_back(next);
            spelled.push_back(c);
            ++candidates_;
            const bool at_end = spelled.size() == L;
            // a fixed-length Model accepts every walk that reaches L (a Pattern always): a
            // walk it does not accept is no instance, neither counted nor supported
            if (at_end && !model.accepting(state, static_cast<uint32_t>(L))) {
                path.pop_back();
                spelled.pop_back();
                continue;
            }
            if (tracker) {
                switch (tracker->push(state_of(next, c, state))) {
                    case SupportTracker::Verdict::ALIVE:
                        ++frames.open;
                        break;
                    case SupportTracker::Verdict::DEAD:
                        // pruned: not expanded; at L a complete walk all the same
                        if (at_end) {
                            ++found_[anchor.orientation];
                            ++found_total_;
                        }
                        prune(anchor.orientation, !at_end);
                        path.pop_back();
                        spelled.pop_back();
                        continue;
                    case SupportTracker::Verdict::STOPPED:
                        external_stop();
                        return false;
                }
            }
            if (!at_end) {
                if (!expand(state))
                    return false;
                continue;
            }
            if (!complete_path(anchor, path, spelled))
                return false;
            if (tracker)
                frames.pop();
            path.pop_back();
            spelled.pop_back();
        }
        if (tracker)
            frames.pop();
        return true;
    }

    // ---------------------------------------------------------------- scans

    static bool is_palindrome(const std::string &kmer) {
        std::string rc = kmer;
        ::reverse_complement(rc.begin(), rc.end());
        return kmer == rc;
    }

    bool is_palindrome(edge_index edge) const {
        return is_palindrome(dbg_succ_.get_node_sequence(edge));
    }

    // one examined edge of a scan: one step
    bool charge_scan() {
        if (!budget_.charge(1)) {
            record_stop(StopPhase::MASK_SCAN, *budget_.stopped());
            return false;
        }
        return true;
    }

    // false when the budget stopped it
    bool scan(Item &item) {
        ++work_.mask_scans;
        if (!masked_) {
            // no mask: only the palindrome scans of an even-k wrapped PRIMARY graph exist
            // (CANDIDATES). Each candidate is spelled for its palindrome test, which shows a
            // source dummy too (its k-mer starts with '$'; never a palindrome, its last base
            // being one): |hits| counts the real ones
            assert(item.scan == Item::Scan::CANDIDATES && item.count_palindromes);
            // the candidates: a flank's non-sink edges, a W rule's edges with W in
            // {c, c + alph_size}
            auto next = [&](edge_index from) {
                return item.flank()
                    ? dbg_succ_.next_non_sink_edge(from, item.last)
                    : dbg_succ_.next_edge_with_last_symbol(from, item.last, item.c);
            };
            for (edge_index e = next(item.first); e; e = next(e + 1)) {
                if (!charge_scan())
                    return false;
                ++item.examined;
                const std::string kmer = dbg_succ_.get_node_sequence(e);
                // the sentinel '$' (BOSS::kSentinel), which only a source dummy starts with
                if (kmer.front() == '$')
                    continue;
                ++item.hits;
                if (is_palindrome(kmer))
                    ++item.palindromes;
            }
            item.scan_done = true;
            return true;
        }
        if (item.scan == Item::Scan::INVALID) {
            for (edge_index e = dbg_succ_.next_invalid_edge(item.first, item.last); e;
                    e = dbg_succ_.next_invalid_edge(e + 1, item.last)) {
                if (!charge_scan())
                    return false;
                ++item.examined;
                TAlphabet w = boss_.get_W(e) % boss_.alph_size;
                if (w)
                    ++item.examined_ns;
                if (w == item.c)
                    ++item.hits;
            }
        } else if (item.flank()) {
            for (edge_index e = dbg_succ_.next_valid_edge(item.first, item.last); e;
                    e = dbg_succ_.next_valid_edge(e + 1, item.last)) {
                if (!charge_scan())
                    return false;
                ++item.examined;
                if (is_palindrome(e))
                    ++item.palindromes;
            }
        } else {
            for (edge_index e = dbg_succ_.next_edge_with_last_symbol(item.first, item.last,
                                                                     item.c);
                    e; e = dbg_succ_.next_edge_with_last_symbol(e + 1, item.last, item.c)) {
                if (!charge_scan())
                    return false;
                ++item.examined;
                if (dbg_succ_.in_graph(e)) {
                    ++item.hits;
                    if (item.count_palindromes && is_palindrome(e))
                        ++item.palindromes;
                }
            }
        }
        item.scan_done = true;
        return true;
    }

    void scan_all(const std::vector<size_t> &order) {
        for (size_t i : order) {
            for (Item &item : searches_[i].pending) {
                if (!scan(item))
                    return;
            }
        }
    }

    // one candidate tested by the check: k - 1 steps, the most symbols BOSS::node_has_sentinel
    // reads
    bool charge_check() {
        if (!budget_.charge(k_ - 1)) {
            record_stop(StopPhase::MASK_SCAN, *budget_.stopped());
            return false;
        }
        return true;
    }

    /**
     * A graph without the mask, discovery and the deferred scans complete: when the pattern's
     * unchecked candidates number at most check_limit_, each is tested with
     * BOSS::node_has_sentinel, in the order of the base searches (as the scans) and of their
     * ranges' discovery, edge by edge: the real k-mers become exact contexts and the source
     * dummies drop out, so that every count of the pattern is EXACT. Nothing is applied before
     * the last candidate is tested: a stop on the way (phase MASK_SCAN) leaves every count as
     * discovery left it, BOUNDS.
     */
    void check_unchecked(const std::vector<size_t> &order) {
        if (!checkable_ || !unchecked_entries_)
            return;

        std::vector<std::vector<uint64_t>> real(searches_.size());
        uint64_t tested = 0;
        for (size_t i : order) {
            BaseSearch &search = searches_[i];
            real[i].assign(search.unchecked.size(), 0);
            for (const CheckRange &range : search.checks) {
                ++work_.mask_scans;
                auto next = [&](edge_index from) {
                    return !range.c
                        ? dbg_succ_.next_non_sink_edge(from, range.last)
                        : dbg_succ_.next_edge_with_last_symbol(from, range.last, range.c);
                };
                for (edge_index e = next(range.first); e; e = next(e + 1)) {
                    if (!charge_check())
                        return;
                    ++tested;
                    real[i][range.offset] += !boss_.node_has_sentinel(e);
                }
            }
        }
        // every unchecked candidate counted is among those tested (a broken invariant is never
        // published as an exact count)
        if (tested != unchecked_entries_) {
            throw std::logic_error("pattern: " + std::to_string(tested)
                                   + " candidates checked of "
                                   + std::to_string(unchecked_entries_) + " unchecked");
        }
        for (size_t i = 0; i < searches_.size(); ++i) {
            BaseSearch &search = searches_[i];
            for (size_t p = 0; p < search.unchecked.size(); ++p) {
                assert(real[i].size() == search.unchecked.size());
                assert(real[i][p] <= search.unchecked[p]);
                search.exact[p] += real[i][p];
                search.unchecked[p] = 0;
            }
            std::vector<CheckRange>().swap(search.checks);
        }
    }

    // ---------------------------------------------------------------- release

    struct Cursor {
        const Span *item;
        Orientation orientation;
        Use use;
        edge_index next;
        Key key;
    };

    // one examined edge of the release; the clock every kClockStride of them
    bool tick() {
        if (++examined_ % Budget::kClockStride || budget_.check_time())
            return true;
        record_stop(release_phase_, StopReason::TIME);
        return false;
    }

    // moves |cursor| to its next context; false when it is exhausted or the clock stopped
    // the release (then *stopped is set)
    bool advance(Cursor &cursor, bool *stopped) {
        const Span &item = *cursor.item;
        while (cursor.next && cursor.next <= item.last) {
            edge_index e = item.flank()
                ? (masked_ ? dbg_succ_.next_valid_edge(cursor.next, item.last)
                           : dbg_succ_.next_non_sink_edge(cursor.next, item.last))
                : dbg_succ_.next_edge_with_last_symbol(cursor.next, item.last, item.c);
            if (!e)
                break;
            cursor.next = e + 1;
            if (!tick()) {
                *stopped = true;
                return false;
            }
            if (!item.flank() && !dbg_succ_.in_graph(e))
                continue;
            // no mask: a source dummy is never released
            if (item.unchecked && boss_.node_has_sentinel(e))
                continue;
            switch (cursor.use) {
                case Use::DIRECT:
                    cursor.key = Key { e, item.offset, cursor.orientation };
                    return true;
                case Use::MAPPED_SKIPPING:
                    // a palindromic stored k-mer is its own reverse complement in the wrapper
                    // (CanonicalDBG::reverse_complement maps it to itself): the PALINDROMES
                    // cursor beside this one releases it at its own id (§4.1)
                    if (is_palindrome(e))
                        continue;
                    [[fallthrough]];
                case Use::MAPPED:
                    cursor.key = Key { e + wrapper_offset_, max_offset() - item.offset,
                                       cursor.orientation };
                    return true;
                case Use::PALINDROMES:
                    if (!is_palindrome(e))
                        continue;
                    cursor.key = Key { e, max_offset() - item.offset, cursor.orientation };
                    return true;
            }
        }
        cursor.next = 0;
        return false;
    }
};

bool PatternRun::release(uint64_t limit, const std::function<void(const Context&)> &emit,
                         bool callers_emit, StopPhase phase, bool *exhausted) {
    release_phase_ = phase;
    if (exhausted)
        *exhausted = false;
    if (!budget_.check_time()) {
        record_stop(phase, StopReason::TIME);
        return false;
    }

    // the cursors over the retained ranges (PARTIAL keeps only those that can hold one of the
    // first cap_ contexts, so this index is small); the clock every kClockStride of them
    std::vector<Cursor> cursors;
    cursors.reserve(retained_ + retained_ / 4);
    bool late = false;
    for_each_use([&](const OrientationPlan &o, const Span &span, Use use) {
        if (late)
            return;
        cursors.push_back(Cursor { &span, o.orientation, use, span.first, {} });
        late = !(cursors.size() % Budget::kClockStride) && !budget_.check_time();
    });
    if (late) {
        record_stop(phase, StopReason::TIME);
        return false;
    }

    // a cursor's contexts are at nodes >= this bound: cursors join the merge only when the
    // merge reaches their bound, so that a small limit examines few ranges
    auto bound = [&](const Cursor &c) { return use_low(*c.item, c.use); };
    std::stable_sort(cursors.begin(), cursors.end(), [&](const Cursor &a, const Cursor &b) {
        return bound(a) < bound(b);
    });
    if (cursors.size() >= Budget::kClockStride && !budget_.check_time()) {
        record_stop(phase, StopReason::TIME);
        return false;
    }

    auto later = [&](size_t a, size_t b) { return cursors[b].key < cursors[a].key; };
    std::priority_queue<size_t, std::vector<size_t>, decltype(later)> heap(later);

    size_t joined = 0;
    uint64_t released = 0;
    bool stopped = false;
    // the last context emitted: a PALINDROMES cursor can meet one the direct search also
    // found, and keys come out of the merge in order, so a duplicate follows its original
    std::optional<Key> last;
    while (released < limit) {
        if (heap.empty() && joined == cursors.size())
            break;

        node_index next_bound = heap.empty() ? bound(cursors[joined])
                                             : cursors[heap.top()].key.node;
        bool any_joined = false;
        while (joined < cursors.size() && bound(cursors[joined]) <= next_bound) {
            if (advance(cursors[joined], &stopped))
                heap.push(joined);
            if (stopped)
                return false;
            ++joined;
            any_joined = true;
        }
        if (any_joined)
            continue;

        size_t top = heap.top();
        heap.pop();
        const Key key = cursors[top].key;
        if (!last || !(key == *last)) {
            // the caller's work per context runs between two clock readings
            if (callers_emit && released && !(released % Budget::kReleaseClockStride)
                    && !budget_.check_time()) {
                record_stop(phase, StopReason::TIME);
                return false;
            }
            emit(Context { key.orientation, key.offset, key.node, base_of(key.node), {}, {} });
            ++released;
            last = key;
        }

        if (advance(cursors[top], &stopped))
            heap.push(top);
        if (stopped)
            return false;
    }
    if (exhausted)
        *exhausted = heap.empty() && joined == cursors.size();
    return true;
}

bool PatternRun::list_anchors(std::vector<Context> *anchors) {
    assert(extending_);
    // the anchors in answer order (§5.5); their ranges are discarded once they are listed: the
    // anchors are kept through the extension, not beyond (§5.2). The clock every
    // kReleaseClockStride anchors listed, as for a caller's emit: on a cold index each costs a
    // few page reads that no step charges
    listed_ = release(kNoLimit, [&](const Context &c) { anchors->push_back(c); }, true,
                      StopPhase::EXTENSION);
    drop_anchors();
    work_.steps = budget_.steps_used() - steps_before_;
    return listed_;
}

bool PatternRun::extend_listed(bool keep, const std::vector<Context> &anchors) {
    assert(extending_ && listed_);
    keep_paths_ = keep;

    bool done = true;
    for (const Context &anchor : anchors) {
        ++anchors_left_[anchor.orientation];
    }
    const Pattern rc = pattern_.reverse_complement();
    for (const Context &anchor : anchors) {
        // the clock before every anchor: its spelling (k - 1 BOSS steps) is work no step
        // charges, and the anchors' DFS charges too few steps to cross a stride (674 anchors
        // spelled on a cold index would run 3 s past the work time unread)
        if (!budget_.check_time()) {
            record_stop(StopPhase::EXTENSION, StopReason::TIME);
            done = false;
            break;
        }
        ++work_.extension_anchors;
        // each oriented pattern is extended in its own reading direction from its own
        // first k-mer (§4.1, "Orientation")
        if (!extend_anchor(anchor, anchor.orientation == Orientation::REVERSE ? rc
                                                                              : pattern_)) {
            done = false;
            break;
        }
        --anchors_left_[anchor.orientation];
    }

    work_.steps = budget_.steps_used() - steps_before_;
    return done;
}

bool PatternRun::extend(bool keep, uint64_t num_anchors) {
    assert(extending_);
    keep_paths_ = keep;

    std::vector<Context> anchors;
    if (!list_anchors(&anchors))
        return false;
    if (anchors.size() != num_anchors) {
        throw std::logic_error("pattern: " + std::to_string(anchors.size())
                               + " anchors listed for the extension of "
                               + std::to_string(num_anchors) + " counted");
    }
    return extend_listed(keep, anchors);
}

// why ALL_OR_COUNT withholds the results of a pattern its stop touched (§5.2)
Withheld withheld_for(StopReason reason) {
    switch (reason) {
        case StopReason::MAX_CONTEXTS:
        case StopReason::MAX_ANCHORS:
        case StopReason::MAX_PATHS:
            return Withheld::THRESHOLD_CROSSED;
        case StopReason::TIME:
            return Withheld::DEADLINE;
        case StopReason::MAX_STEPS:
            return Withheld::DISCOVERY_BUDGET;
        case StopReason::EXTERNAL:
            return Withheld::EXTERNAL;
    }
    return Withheld::DISCOVERY_BUDGET;
}

/**
 * The release of a long pattern's paths (Request::extend_paths), from what the extension
 * counted and kept: see PatternSearch::enumerate for the rules per mode. |released| is the count
 * of the paths the rules read: AnchorCounts::paths, or with a support tracker
 * AnchorCounts::supported (EXACT when the extension COMPLETED, whatever it pruned).
 */
Extraction extract_paths(const AnchorCounts &anchors, const Count &released,
                         const Request &request, PatternRun &engine,
                         const std::function<void(const Context&)> &callback) {
    Extraction extraction;
    const bool all_or_count = request.mode == Mode::ALL_OR_COUNT;
    auto release = [&]() {
        return engine.release_paths([&](const Context &path) {
            callback(path);
            ++extraction.returned;
        });
    };

    switch (anchors.extension) {
        case Extension::NOT_REQUESTED:
            throw std::logic_error("pattern: paths released without extension");

        case Extension::NO_ANCHORS:
            // no anchor, no path: the empty answer is complete
            extraction.complete = true;
            break;

        case Extension::NOT_ADMITTED:
            extraction.withheld = Withheld::ANCHORS_ABOVE_THRESHOLD;
            break;

        case Extension::NOT_STARTED: {
            // the anchors stopped (discovery, a mask scan, or their threshold): nothing was
            // extended, so nothing is delivered in either mode
            assert(engine.stop());
            const StopReason reason = engine.stop() ? engine.stop()->reason
                                                    : StopReason::MAX_STEPS;
            if (all_or_count) {
                extraction.withheld = withheld_for(reason);
            } else {
                extraction.cut = reason;
            }
            break;
        }

        case Extension::STOPPED: {
            assert(engine.stop());
            const StopReason reason = engine.stop() ? engine.stop()->reason
                                                    : StopReason::MAX_STEPS;
            if (all_or_count) {
                extraction.withheld = withheld_for(reason);
            } else if (reason == StopReason::TIME) {
                // membership depends on the machine: nothing (as for L <= k)
                extraction.cut = StopReason::TIME;
            } else {
                // the paths completed before the stop are a prefix of the answer order
                extraction.cut = release() ? reason : StopReason::TIME;
            }
            break;
        }

        case Extension::COMPLETED: {
            assert(released.relation == Relation::EXACT);
            const uint64_t paths = released.value;
            if (all_or_count) {
                if (paths > request.max_paths) {
                    extraction.withheld = Withheld::COUNT_ABOVE_THRESHOLD;
                } else if (engine.paths_retained() != paths) {
                    throw std::logic_error("pattern: the extension kept "
                                           + std::to_string(engine.paths_retained())
                                           + " of " + std::to_string(paths) + " paths");
                } else if (release()) {
                    extraction.complete = true;
                } else {
                    // a deadline during the delivery: the prefix delivered is the caller's to
                    // discard (all or nothing)
                    extraction.withheld = Withheld::DEADLINE;
                    extraction.returned = 0;
                }
            } else if (!release()) {
                extraction.cut = StopReason::TIME;
            } else if (extraction.returned == paths) {
                extraction.complete = true;
            } else {
                extraction.cut = StopReason::MAX_PATHS;
            }
            break;
        }
    }
    return extraction;
}

} // namespace


namespace {

// Wilson's score interval at 95% around real / samples, clamped to [0, 1]; it holds the
// point estimate (min and max only undo a rounding at 0 and 1)
void set_wilson_interval(RealFraction *f) {
    if (!f->samples)
        return;
    const double n = static_cast<double>(f->samples);
    const double p = static_cast<double>(f->real) / n;
    const double z = 1.959963984540054;
    const double denominator = 1 + z * z / n;
    const double centre = (p + z * z / (2 * n)) / denominator;
    const double half = z / denominator * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n));
    f->value = p;
    // (the interval holds p; min and max only undo the rounding at p = 0 or 1)
    f->lower = std::min(p, std::max(0.0, centre - half));
    f->upper = std::max(p, std::min(1.0, centre + half));
}

// the entries no pattern base matches: W = $, plain or marked (the sinks and edge 1)
uint64_t sentinel_edges(const BOSS &boss) {
    const auto sentinel = static_cast<TAlphabet>(BOSS::kSentinelCode);
    return boss.rank_W(boss.num_edges(), sentinel)
         + boss.rank_W(boss.num_edges(), sentinel + boss.alph_size);
}

} // namespace

RealFraction sample_real_fraction(const DBGSuccinct &graph, uint64_t samples) {
    const BOSS &boss = graph.get_boss();
    RealFraction f;
    f.edges = boss.num_edges();
    f.sentinel_edges = sentinel_edges(boss);
    f.seed = f.edges;
    if (f.sentinel_edges >= f.edges || !samples)
        return f;

    // the generator the standard specifies (not a distribution of the library, whose draws
    // differ between implementations): the same draws on every platform
    std::mt19937_64 rng(f.seed);
    while (f.samples < samples) {
        // uniform over [1, edges]; an entry with W = $ (W modulo alph_size 0, plain or
        // marked) is drawn again, so that the draws are uniform over the entries a pattern
        // can count
        const uint64_t edge = 1 + rng() % f.edges;
        if (!(boss.get_W(edge) % boss.alph_size))
            continue;
        ++f.samples;
        f.real += !boss.node_has_sentinel(edge);
    }
    set_wilson_interval(&f);
    return f;
}

RealFraction exact_real_fraction(const DBGSuccinct &graph) {
    const BOSS &boss = graph.get_boss();
    RealFraction f;
    f.edges = boss.num_edges();
    f.sentinel_edges = sentinel_edges(boss);
    f.seed = f.edges;
    f.exact = true;
    f.samples = f.edges - f.sentinel_edges;
    // the source dummies by BOSS's own traversal of the dummy tree (`stats --count-dummy`),
    // not by the walk sample_real_fraction tests with: those with W != $ are the dummies
    // among the entries
    sdsl::bit_vector source_dummies(boss.get_W().size(), false);
    boss.mark_source_dummy_edges(&source_dummies, 1);
    uint64_t dummies = 0;
    for (uint64_t e = 1; e < source_dummies.size(); ++e) {
        // (W modulo alph_size is 0 for $, plain or marked: BOSS::kSentinelCode)
        if (source_dummies[e] && boss.get_W(e) % boss.alph_size)
            ++dummies;
    }
    assert(dummies <= f.samples);
    f.real = f.samples - dummies;
    if (f.samples) {
        f.value = static_cast<double>(f.real) / static_cast<double>(f.samples);
        f.lower = f.upper = f.value;
    }
    return f;
}

GraphSupport PatternSearch::support(const DeBruijnGraph &graph) {
    GraphSupport result;
    result.k = graph.get_k();

    const DBGSuccinct *dbg_succ = nullptr;
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph)) {
        dbg_succ = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());
        if (!dbg_succ || dbg_succ->get_mode() != DeBruijnGraph::PRIMARY) {
            result.reason = "representation_unsupported";
            return result;
        }
        result.mode = GraphMode::PRIMARY;
    } else if ((dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph))) {
        switch (dbg_succ->get_mode()) {
            case DeBruijnGraph::BASIC:
                result.mode = GraphMode::BASIC;
                break;
            case DeBruijnGraph::CANONICAL:
                result.mode = GraphMode::CANONICAL;
                break;
            case DeBruijnGraph::PRIMARY:
                // its contexts would be the stored orientation's only (§4.1)
                result.mode = GraphMode::PRIMARY;
                result.reason = "primary_unwrapped";
                return result;
        }
    } else {
        result.reason = "representation_unsupported";
        return result;
    }

    result.alphabet = dbg_succ->get_boss().alphabet;
    result.mask_present = dbg_succ->get_mask() != nullptr;
    result.strand_stated = result.mode == GraphMode::BASIC;
    if (result.mode == GraphMode::PRIMARY) {
        result.scopes = { Scope::ANY_OFFSET };
    } else {
        result.scopes = { Scope::SUFFIX, Scope::ANY_OFFSET };
    }

    if (result.alphabet != "$ACGT" && result.alphabet != "$ACGTN") {
        result.reason = "alphabet_unsupported";
        return result;
    }
    if (result.k < 2) {
        // a BOSS node is a (k-1)-mer: the range narrowing needs one symbol at least
        result.reason = "representation_unsupported";
        return result;
    }
    // without the mask the graph is served all the same: its unresolved counts are upper bounds
    // (PatternSearch, "Graphs without the dummy-edge mask"); no reason is mask_required

    result.supported = true;
    return result;
}

PatternSearch::PatternSearch(const DeBruijnGraph &graph)
      : graph_(graph), support_(support(graph)) {
    if (!support_.supported)
        throw std::invalid_argument("pattern: graph not supported: " + support_.reason);

    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph_)) {
        dbg_succ_ = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());
        // CanonicalDBG numbers the reverse complement of a stored node y as y plus the
        // stored graph's max_index (canonical_dbg.cpp, offset_), and a palindromic y as y
        wrapper_offset_ = canonical->get_graph().max_index();
    } else {
        dbg_succ_ = dynamic_cast<const DBGSuccinct*>(&graph_);
        wrapper_offset_ = 0;
    }
    assert(dbg_succ_);
}

Result PatternSearch::count(const Pattern &pattern, const Request &request,
                            Budget &budget) const {
    return run(pattern, request, budget, nullptr);
}

Result PatternSearch::enumerate(const Pattern &pattern, const Request &request,
                                Budget &budget,
                                const std::function<void(const Context&)> &callback) const {
    if (request.mode == Mode::COUNT)
        throw std::invalid_argument("pattern: enumerate() needs mode all_or_count or partial");
    if (request.extend_paths && request.release_anchors) {
        throw std::invalid_argument("pattern: release_anchors and extend_paths exclude each "
                                    "other (the results of a long pattern are its paths)");
    }
    if (request.sink && request.extend_paths && pattern.length() > support_.k) {
        throw std::invalid_argument("pattern: the paths of the extension go to the path sink, "
                                    "which is their release: use count()");
    }

    return run(pattern, request, budget, &callback);
}

Result PatternSearch::run(const Pattern &pattern, const Request &request, Budget &budget,
                          const std::function<void(const Context&)> *callback) const {
    auto start = std::chrono::steady_clock::now();

    const size_t k = support_.k;
    const size_t L = pattern.length();
    const bool is_long = L > k;

    Result result;
    result.scope = is_long ? Scope::LONG : request.scope;
    result.graph_mode = support_.mode;
    result.palindromic = pattern.is_palindromic();
    result.information_bits = pattern.information_bits();
    // the bits the information floor reads: for a long pattern the least over the searched
    // orientations' anchor windows (FORWARD P[0, k), REVERSE rc(P)[0, k), whose bits are
    // P[L - k, L)'s), so that no searched window below the floor runs (kept in a plain double
    // rather than read back from the optional, which GCC 13 -O3 flags as maybe uninitialized
    // under -Werror). anchor_information_bits keeps its contract-version-1 meaning, the bits of
    // P[0, k) whatever the strands; the least is stated apart, in min_anchor_information_bits
    double floor_bits = result.information_bits;
    // the orientation whose window is the least informative (FORWARD before REVERSE on a tie)
    Orientation floor_window = Orientation::FORWARD;
    if (is_long) {
        bool first = true;
        for (Orientation o : searched_orientations(pattern, request.strands)) {
            const double bits = o == Orientation::REVERSE
                ? pattern.information_bits(L - k, L)
                : pattern.information_bits(0, k);
            result.anchor_window_bits[o] = bits;
            if (first || bits < floor_bits) {
                floor_bits = bits;
                floor_window = o;
            }
            first = false;
        }
        result.anchor_information_bits = pattern.information_bits(0, k);
        result.min_anchor_information_bits = floor_bits;
    }

    auto finish = [&]() {
        result.elapsed_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        return result;
    };

    // per-pattern refusals (§7.2), decided before anything is charged
    if (!is_long && request.scope == Scope::SUFFIX && support_.mode == GraphMode::PRIMARY) {
        result.refusal = Refusal {
            "scope_unsupported",
            "pattern: scope suffix is not served on a primary graph (the suffix of a virtual "
            "k-mer is the prefix of a stored one); use scope any_offset"
        };
        return finish();
    }
    // a pattern without instances ('*' read in a table without an unconditional stop codon) and
    // L <= k, where a context instantiates the whole pattern: answered EXACT 0 with nothing
    // searched, so that there is nothing for the information floor to gate. A longer one is
    // searched as any other: its anchors instantiate the anchor window only, which may not
    // reach the '*', and its paths (none) are the extension's
    if (!pattern.has_instances() && !is_long) {
        answer_no_instance(pattern, request, callback != nullptr, &result);
        return finish();
    }
    // the information floor gates discovery, except for an exact pattern in suffix scope:
    // one range, a few ranks (§5.3)
    if (!(result.scope == Scope::SUFFIX && pattern.is_exact())) {
        const double bits = floor_bits;
        if (bits < request.min_information_bits) {
            result.refusal = Refusal {
                "information_below_floor",
                "pattern: " + format_bits(bits) + " information bits"
                    + (!is_long ? ""
                        : floor_window == Orientation::REVERSE
                            ? " in the anchor window of the reverse orientation (the reverse "
                              "complement of the pattern's last k positions)"
                            : " in the anchor window")
                    + ", below the floor of " + format_bits(request.min_information_bits)
                    + " for scope " + to_string(result.scope)
            };
            return finish();
        }
    }

    // phase 2 (§4.2): a long pattern's anchors extended to paths, on request
    const bool extending = is_long && request.extend_paths;
    // the release of contexts through PatternRun::release: the contexts of L <= k, or the
    // anchors of a long pattern on request (Request::release_anchors, never with extension)
    const bool releasing = callback && !extending && (!is_long || request.release_anchors);
    // the threshold of ALL_OR_COUNT and the cap of PARTIAL
    const uint64_t max_released = is_long ? request.max_anchors : request.max_contexts;
    // the anchors are retained for the extension in every mode (§5.2)
    PatternRun engine(graph_, *dbg_succ_, support_, wrapper_offset_, pattern, request, budget,
                      releasing || extending);
    engine.run();

    result.searched = engine.searched();

    // the count of every (orientation, offset) of the served graph, in plan order: what the
    // totals sum (§3)
    const Unit unit = is_long ? Unit::ANCHORS : Unit::GRAPH_CONTEXTS;
    struct Cell {
        Orientation orientation;
        uint32_t offset;
        Count count;
    };
    std::vector<Cell> cells;
    for (const OrientationPlan &o : engine.plans()) {
        for (uint32_t p : engine.offsets()) {
            cells.push_back(Cell { o.orientation, p, engine.count(o, p, unit) });
        }
    }

    // counts by orientation and offset, summed with the weakest relation (§3)
    std::map<uint32_t, Count> by_offset;
    std::map<Orientation, Count> by_orientation;
    std::optional<Count> total;
    auto aggregate = [&]() {
        by_offset.clear();
        by_orientation.clear();
        total.reset();
        // one plan's cells are consecutive, one plan per orientation
        for (size_t i = 0; i < cells.size(); ) {
            const Orientation o = cells[i].orientation;
            std::optional<Count> orientation_total;
            for (; i < cells.size() && cells[i].orientation == o; ++i) {
                const Count &c = cells[i].count;
                auto [it, inserted] = by_offset.emplace(cells[i].offset, c);
                if (!inserted)
                    it->second += c;
                if (orientation_total) {
                    *orientation_total += c;
                } else {
                    orientation_total = c;
                }
            }
            by_orientation.emplace(o, *orientation_total);
            if (total) {
                *total += *orientation_total;
            } else {
                total = orientation_total;
            }
        }
        assert(total);
    };
    aggregate();

    /**
     * A graph without the dummy-edge mask: the contexts (or anchors) a release enumerated, per
     * orientation and offset. A release that enumerated every candidate of a completed
     * discovery makes each count the number it released (exactify); one cut short raises each
     * count's lower bound to it (raise_lower). Real contexts only (the release drops the source
     * dummies), so a number outside a count's bounds is a broken invariant, never published.
     */
    typedef std::map<std::pair<Orientation, uint32_t>, uint64_t> Released;
    auto released_at = [](const Released &released, const Cell &cell) -> uint64_t {
        auto it = released.find({ cell.orientation, cell.offset });
        return it == released.end() ? 0 : it->second;
    };
    auto broken = [](const Cell &cell, uint64_t n) {
        const Count &c = cell.count;
        return std::logic_error(
            "pattern: " + std::to_string(n) + " released at offset "
            + std::to_string(cell.offset) + " " + orientation_key(cell.orientation)
            + " of a count " + to_string(c.relation) + " " + std::to_string(c.value) + ".."
            + std::to_string(c.relation == Relation::BOUNDS ? c.upper : c.value));
    };
    auto exactify = [&](const Released &released) {
        for (Cell &cell : cells) {
            const uint64_t n = released_at(released, cell);
            const Count &c = cell.count;
            const bool within = c.relation == Relation::EXACT
                ? n == c.value
                : c.relation == Relation::BOUNDS && c.lower <= n && n <= c.upper;
            if (!within)
                throw broken(cell, n);
            cell.count = Count::exact(unit, n);
        }
        aggregate();
    };
    auto raise_lower = [&](const Released &released) {
        for (Cell &cell : cells) {
            const uint64_t n = released_at(released, cell);
            Count &c = cell.count;
            switch (c.relation) {
                case Relation::EXACT:
                    if (n > c.value)
                        throw broken(cell, n);
                    break;
                case Relation::BOUNDS:
                    if (n > c.upper)
                        throw broken(cell, n);
                    if (n > c.lower)
                        c = Count::bounds(unit, n, c.upper);
                    break;
                case Relation::AT_LEAST:
                    if (n > c.value)
                        c = Count::at_least(unit, n);
                    break;
                case Relation::UNKNOWN:
                    if (n)
                        c = Count::at_least(unit, n);
                    break;
            }
        }
        aggregate();
    };
    // a graph without the mask, discovery completed: a count BOUNDS for its unchecked
    // candidates only, which a release resolves
    auto unresolved = [&]() {
        return !support_.mask_present && !engine.stop()
                && total->relation == Relation::BOUNDS;
    };

    AnchorCounts anchors;
    if (is_long) {
        // no path starts without an anchor (a derivation, not a promotion); otherwise the
        // paths need the extension (§4.2)
        const bool no_anchor = total->relation == Relation::EXACT && total->value == 0;
        if (no_anchor)
            anchors.paths = Count::exact(Unit::PATHS, 0);

        if (extending) {
            bool ran = false;
            if (no_anchor) {
                anchors.extension = Extension::NO_ANCHORS;
            } else if (unresolved()) {
                // no mask: the admission compares the upper bound
                if (total->upper > request.max_anchors) {
                    anchors.extension = Extension::NOT_ADMITTED;
                    engine.drop_anchors();
                    if (total->lower <= request.max_anchors)
                        engine.note_upper_bound_decision();
                } else {
                    // the listing drops the source dummies: the anchors it lists are the
                    // exact set, extended as on a masked graph (the running count the
                    // retention compares never exceeds the final U: nothing was dropped)
                    if (!engine.retained_all())
                        throw std::logic_error("pattern: anchors dropped for an admitted "
                                               "extension");
                    auto extension_start = std::chrono::steady_clock::now();
                    std::vector<Context> listed;
                    if (!engine.list_anchors(&listed)) {
                        anchors.extension = Extension::STOPPED;
                        ran = true;
                    } else {
                        Released counted;
                        for (const Context &c : listed) {
                            ++counted[{ c.orientation, c.offset }];
                        }
                        exactify(counted);
                        if (listed.empty()) {
                            anchors.extension = Extension::NO_ANCHORS;
                            anchors.paths = Count::exact(Unit::PATHS, 0);
                        } else {
                            const bool keep = callback != nullptr;
                            anchors.extension = engine.extend_listed(keep, listed)
                                    ? Extension::COMPLETED
                                    : Extension::STOPPED;
                            ran = true;
                        }
                    }
                    result.extension_ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - extension_start).count();
                    anchors.candidates_examined = engine.candidates_examined();
                }
            } else if (total->relation != Relation::EXACT || engine.stop()) {
                // an anchor set not known completely is never extended (§3)
                anchors.extension = Extension::NOT_STARTED;
            } else if (total->value > request.max_anchors) {
                // the extension's admission (§4.2); the anchors are discarded (§5.2)
                anchors.extension = Extension::NOT_ADMITTED;
                engine.drop_anchors();
            } else {
                auto extension_start = std::chrono::steady_clock::now();
                anchors.extension = engine.extend(callback != nullptr, total->value)
                        ? Extension::COMPLETED
                        : Extension::STOPPED;
                result.extension_ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - extension_start).count();
                anchors.candidates_examined = engine.candidates_examined();
                ran = true;
            }

            // per orientation: EXACT 0 without anchors; EXACT when every anchor of it was
            // extended, AT_LEAST when the extension stopped before; UNKNOWN when it did not
            // run. With a support tracker the walks are a plain count, AT_LEAST also when a
            // branch was pruned before L, and the supported paths follow the rule of the
            // extension alone
            const bool tracked = request.support != nullptr;
            std::optional<Count> sum;
            std::optional<Count> supported_sum;
            for (Orientation o : engine.searched()) {
                const Count &a = by_orientation.at(o);
                Count paths = Count::unknown(Unit::PATHS);
                Count supported = Count::unknown(Unit::PATHS);
                if (a.relation == Relation::EXACT && a.value == 0) {
                    paths = Count::exact(Unit::PATHS, 0);
                    supported = Count::exact(Unit::PATHS, 0);
                } else if (ran) {
                    const bool extended = engine.orientation_extended(o);
                    paths = extended && !engine.pruned_before_completion(o)
                        ? Count::exact(Unit::PATHS, engine.paths_found(o))
                        : Count::at_least(Unit::PATHS, engine.paths_found(o));
                    supported = extended
                        ? Count::exact(Unit::PATHS, engine.supported_found(o))
                        : Count::at_least(Unit::PATHS, engine.supported_found(o));
                }
                anchors.paths_by_orientation.emplace(o, paths);
                if (sum) {
                    *sum += paths;
                } else {
                    sum = paths;
                }
                if (!tracked)
                    continue;
                anchors.supported_by_orientation.emplace(o, supported);
                if (supported_sum) {
                    *supported_sum += supported;
                } else {
                    supported_sum = supported;
                }
            }
            if (ran) {
                assert(sum);
                anchors.paths = *sum;
                anchors.walks = engine.walks();
                assert((anchors.paths.relation == Relation::EXACT)
                            == (anchors.extension == Extension::COMPLETED
                                    && !engine.pruned_before_completion()));
            }
            if (tracked) {
                assert(supported_sum);
                if (ran) {
                    anchors.supported = *supported_sum;
                } else if (anchors.extension == Extension::NO_ANCHORS) {
                    anchors.supported = Count::exact(Unit::PATHS, 0);
                }
                anchors.branches_pruned = engine.branches_pruned();
                anchors.pruned_before_completion = engine.pruned_before_completion();
            }
        }
        anchors.total = *total;
        anchors.by_orientation = by_orientation;
    }

    // the label-free release (§4.3, §5.2 "Projection none")
    if (callback) {
        Extraction extraction;
        const bool exact = total->relation == Relation::EXACT && !engine.stop();
        // no mask: discovery complete, the count BOUNDS for its unchecked candidates only
        const bool resolvable = unresolved();
        std::optional<StopReason> reason;
        if (engine.stop())
            reason = engine.stop()->reason;

        if (extending) {
            extraction = extract_paths(anchors, request.support ? anchors.supported
                                                                : anchors.paths,
                                       request, engine, *callback);

        } else if (!releasing) {
            // the results of a long pattern are its paths, not extended without
            // extend_paths; with no anchor there is no path, so the empty answer is complete
            if (exact && total->value == 0) {
                extraction.complete = true;
            } else {
                extraction.withheld = Withheld::PATHS_LATER_INCREMENT;
            }

        } else if (request.mode == Mode::ALL_OR_COUNT) {
            // without the mask the threshold is compared with the upper bound (conservative),
            // and the release resolves the count
            if (resolvable && total->upper <= max_released && !engine.retained_all()) {
                // the running count the retention compares never exceeds the final U
                throw std::logic_error("pattern: ranges dropped for an admitted release");
            }
            if ((exact && total->value <= max_released)
                    || (resolvable && total->upper <= max_released)) {
                // buffered, so that a deadline in the release withholds everything (§5.2),
                // then delivered under the clock: the caller's work per context (a spelling,
                // a result object) is work too, and a deadline during it withholds everything
                // as well, the prefix delivered being the caller's to discard
                std::vector<Context> buffer;
                if (engine.release(kNoLimit, [&](const Context &c) { buffer.push_back(c); },
                                   false)) {
                    // the absence licence rests on this: checked in every build
                    if (exact && buffer.size() != total->value) {
                        throw std::logic_error("pattern: " + std::to_string(buffer.size())
                                               + " contexts released of "
                                               + std::to_string(total->value)
                                               + " counted exactly");
                    }
                    if (resolvable) {
                        // every candidate enumerated, the source dummies dropped: the
                        // counts are what was released
                        Released released;
                        for (const Context &c : buffer) {
                            ++released[{ c.orientation, c.offset }];
                        }
                        exactify(released);
                    }
                    if (engine.deliver(buffer, *callback)) {
                        extraction.returned = buffer.size();
                        extraction.complete = true;
                    } else {
                        extraction.withheld = Withheld::DEADLINE;
                    }
                } else {
                    extraction.withheld = Withheld::DEADLINE;
                }
            } else if (exact || resolvable) {
                extraction.withheld = Withheld::COUNT_ABOVE_THRESHOLD;
                if (resolvable && total->lower <= max_released)
                    engine.note_upper_bound_decision();
            } else if (reason == StopReason::MAX_CONTEXTS || reason == StopReason::MAX_ANCHORS) {
                extraction.withheld = Withheld::THRESHOLD_CROSSED;
            } else if (reason == StopReason::TIME) {
                extraction.withheld = Withheld::DEADLINE;
            } else {
                extraction.withheld = Withheld::DISCOVERY_BUDGET;
            }

        } else {
            // PARTIAL: what was discovered, in answer order, up to max_contexts, every cut
            // stated; nothing after a time stop, whose membership depends on the machine
            if (reason == StopReason::TIME) {
                extraction.cut = StopReason::TIME;
            } else {
                const bool masked = support_.mask_present;
                Released released;
                bool exhausted = false;
                const bool finished = engine.release(max_released, [&](const Context &c) {
                    (*callback)(c);
                    ++extraction.returned;
                    if (!masked)
                        ++released[{ c.orientation, c.offset }];
                }, true, StopPhase::EXTRACTION, &exhausted);
                bool now_exact = exact;
                if (!masked && !exact) {
                    // no mask: a release that drained every retained range of a completed
                    // discovery, none dropped, enumerated every candidate
                    if (finished && resolvable && exhausted && engine.retained_all()) {
                        exactify(released);
                        now_exact = true;
                    } else {
                        raise_lower(released);
                    }
                }
                if (finished) {
                    // an exact count and a release run to its end must agree, as in
                    // all_or_count: min(count, cap) contexts, else the list would be stated
                    // complete short of the count, or cut by a cap it did not reach
                    if (now_exact && (extraction.returned > total->value
                                        || (extraction.returned < total->value
                                                && extraction.returned < max_released))) {
                        throw std::logic_error("pattern: "
                                               + std::to_string(extraction.returned)
                                               + " contexts released of "
                                               + std::to_string(total->value)
                                               + " counted exactly (cap "
                                               + std::to_string(max_released) + ")");
                    }
                    if (now_exact && extraction.returned == total->value) {
                        extraction.complete = true;
                    } else {
                        // the cap that cut the list: max_anchors when anchors are released
                        extraction.cut = reason ? *reason
                                                : is_long ? StopReason::MAX_ANCHORS
                                                          : StopReason::MAX_CONTEXTS;
                    }
                } else {
                    extraction.cut = StopReason::TIME;
                }
            }
        }
        result.extraction = extraction;
    }

    if (is_long) {
        // (the counts as the release left them: EXACT once it enumerated every candidate)
        anchors.total = *total;
        anchors.by_orientation = std::move(by_orientation);
        result.anchors = std::move(anchors);
    } else {
        ContextCounts contexts;
        contexts.total = *total;
        contexts.suffix = by_offset.at(static_cast<uint32_t>(k - L));
        contexts.by_offset = std::move(by_offset);
        contexts.by_orientation = std::move(by_orientation);
        result.contexts = std::move(contexts);
    }

    result.work = engine.work();
    result.stop = engine.stop();
    result.time_limited = engine.time_limited();

    // the bases of an exact pattern: its text, or for a peptide (every residue one codon) the
    // codons it spells. An optional diagnostic: not run after any stop, and left out when the
    // work time passes before sdust has its answer, a time stop the answer states as
    // time_limited only (its counts complete, its stop none)
    bool low_complexity = false;
    if (pattern.is_exact() && !result.stop) {
        const std::optional<bool> flagged = is_low_complexity(
                pattern.kind() == PatternKind::PROTEIN ? exact_bases(pattern) : pattern.text(),
                budget);
        low_complexity = flagged.value_or(false);
        result.time_limited |= !flagged;
    }
    if (low_complexity)
        result.notes.push_back(kNoteLowComplexity);
    if (support_.mode != GraphMode::BASIC)
        result.notes.push_back(kNoteStrandUnknown);
    if (is_long && !extending)
        result.notes.push_back(kNotePathsLater);
    if (engine.upper_bound_decision())
        result.notes.push_back(kNoteThresholdUpperBound);
    if (!pattern.has_instances())
        result.notes.push_back(kNoteNoStopCodon);

    return finish();
}

void PatternSearch::answer_no_instance(const Pattern &pattern, const Request &request,
                                       bool enumerating, Result *result) const {
    const size_t k = support_.k;
    const size_t L = pattern.length();
    assert(!pattern.has_instances() && L <= k);
    // the orientations a search would have had, each with every count EXACT 0: an absence
    // derived from the pattern, nothing searched or charged
    result->searched = searched_orientations(pattern, request.strands);
    ContextCounts contexts;
    contexts.total = Count::exact(Unit::GRAPH_CONTEXTS, 0);
    contexts.suffix = Count::exact(Unit::GRAPH_CONTEXTS, 0);
    if (request.scope == Scope::SUFFIX) {
        contexts.by_offset.emplace(static_cast<uint32_t>(k - L),
                                   Count::exact(Unit::GRAPH_CONTEXTS, 0));
    } else {
        for (uint32_t p = 0; p + L <= k; ++p) {
            contexts.by_offset.emplace(p, Count::exact(Unit::GRAPH_CONTEXTS, 0));
        }
    }
    for (Orientation o : result->searched) {
        contexts.by_orientation.emplace(o, Count::exact(Unit::GRAPH_CONTEXTS, 0));
    }
    result->contexts = std::move(contexts);
    if (enumerating) {
        // nothing to release, and that is all of it
        Extraction extraction;
        extraction.complete = true;
        result->extraction = extraction;
    }
    if (support_.mode != GraphMode::BASIC)
        result->notes.push_back(kNoteStrandUnknown);
    result->notes.push_back(kNoteNoStopCodon);
}

} // namespace pattern
} // namespace graph
} // namespace mtg
