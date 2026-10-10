#ifndef __TRAVERSAL_SUPPORT_STEP_HPP__
#define __TRAVERSAL_SUPPORT_STEP_HPP__

#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "traversal_types.hpp"


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The support step: how the support a label gives a walk narrows when the walk takes one step
 * (the walker's step rule, factored out as pure linear merges so that the pattern's
 * supported-path search and the walker can share it). Nothing here reads the annotation, the
 * graph or a clock of its own: the caller reads the rows, says which k-mer each row belongs
 * to, charges the units a call returns and owns the deadline.
 *
 * The two support levels are the walker's (traversal_types.hpp):
 *  - Support::KMER (`label_intersection`): a label supports a walk when it annotates every
 *    k-mer of it. A step intersects the walk's label set with the next k-mer's row.
 *  - Support::TRACE (`record_verified`): a label supports a walk when one of its records holds
 *    the whole walk, i.e. a column coordinate c of the first k-mer with c + i a coordinate of
 *    the i-th, all in one record. A step continues the CHAINS: the walker keeps a coordinate x
 *    of the next k-mer when x - 1 is live on its right arm and x + 1 on its left arm
 *    (walker.cpp, process_item and validate_seed: a binary search per coordinate). Here the
 *    live coordinates are runs of consecutive coordinates (ChainRun), shifted by one and
 *    intersected with the row's coordinates turned into runs once per row (CoordRun): linear
 *    in runs, so a homopolymer's thousands of coordinates are one run and one comparison per
 *    step. Each run also carries the end of its
 *    record in the direction it moves, and a chain is not continued past it: consecutive column
 *    coordinates in two records of one column are no occurrence (the walker's column-label
 *    trace does not see this, `trace_record_boundaries`). At this level a frame also keeps the
 *    KMER membership beside the chains, so that a label carried on every k-mer but held whole
 *    by no record can still be listed with `label_intersection` support (DESIGN-pattern-search.md §4.3).
 *    Only the anchor frame (depth 0) stores the chain runs: a chain's coordinate at depth d is
 *    its anchor coordinate plus d (UP) or minus d (DOWN) and its record end never changes, so
 *    a frame below the anchor keeps one bit per chain of the anchor (set while the chain
 *    carries the walk) and reads the runs through the anchor's chains, which the frames of a
 *    walk share (AnchorChains). A step scans a label's set bits as runs of consecutive
 *    surviving chains of one anchor run — the same runs the merge would have stored — and
 *    merges them with the row's coordinate runs as before; a walk of n k-mers over a k-mer
 *    carried by tens of thousands of records costs n bitmaps of as many bits, not n copies of
 *    its chain runs (the memory measurement of 2026-10-10, the owner's decision on it).
 *
 * Strands (one orientation per walk: the strand is consistent for one search direction and
 * never flips half-way): a label supports a walk in ONE orientation as a whole —
 * every k-mer of the walk as spelled, or every k-mer of its reverse-complement walk, never
 * some on one strand and the rest on the other. A Frame therefore has one Orientation, fixed
 * when it is opened, and every row it is given must be the row of the k-mer that orientation
 * asks for at that step: the frame spells the walk itself (it holds the walk's current k-mer
 * and extends it by the base of each step) and refuses (std::invalid_argument) a row whose
 * k-mer is not the spelled one (SPELLED) or its reverse complement (REVERSE_COMPLEMENT). A
 * label set or a chain of one orientation is never continued with the other's row, and
 * nothing here unites the rows of the two strands at a k-mer: the two orientations meet only
 * at the end of a walk, label by label (combine()), each label with the orientation that
 * supported it.
 *
 * How the reverse-complement walk is checked: walking w = x_0 .. x_{n-1} to the right, the
 * REVERSE_COMPLEMENT frame reads at step i the row of rc(x_i). Its label set is the
 * intersection of those rows (rc(w)'s k-mers, in the reverse order: the same set). A record
 * holding rc(w) holds rc(x_{n-1}) at some coordinate d and rc(x_i) at d + n - 1 - i, so the
 * frame opens its chains at rc(x_0)'s coordinates (d + n - 1) and moves them DOWN by one per
 * step, each bounded by its record's FIRST coordinate; after the last step a surviving chain
 * stands at d, where rc(w) starts. On the left arm everything is mirrored (direction()). On
 * CANONICAL and PRIMARY graphs a k-mer and its reverse complement share one row and the
 * coordinates carry no strand: only SPELLED frames at the KMER level are meaningful there.
 *
 * Units, charged by the caller: 1 per label of either list a label merge
 * compares, 1 per chain run and per coordinate run of a label a chain merge compares (a chain
 * run: a maximal run of consecutive chains of one record that still carry the walk — the
 * anchor's runs as opened, below it the runs of set bits inside one anchor run), 1 per
 * coordinate turned into runs (a row's, once), 1 per label, per run and per record lookup when
 * a frame is opened. Every call returns its units; Frame::price() gives, before a step, an upper bound of
 * its units and of the bytes of the frame it builds, so that the caller can gate the work and
 * admit the memory before the step and charge the exact units after it (before its next call).
 * A long call reads the caller's Clock before every kPollIterations rounds of its loops (a
 * round compares or maps one or two elements), never less often, and stops when it is expired
 * (Clock::stopped; its output is then incomplete).
 *
 * Memory (the walker's CostModel convention: an element of a vector that grows is charged
 * twice its size; a bitmap assigned to its size once): Frame::bytes() and RowRuns::bytes(),
 * and their bounds before they are built. The anchor's chains count in the anchor frame's
 * bytes only; a frame below it is charged its labels and its bits.
 */
namespace support_step {

// The walk whose rows a frame reads (the strand rule, above)
enum class Orientation : uint8_t {
    SPELLED = 0,             // the rows of the walk's k-mers as it spells them
    REVERSE_COMPLEMENT = 1,  // the rows of their reverse complements: the reverse-complement walk
};
const char* to_string(Orientation orientation);

// How the coordinates of a frame's chains move when its walk takes a step
enum class Direction : uint8_t { DOWN = 0, UP = 1 };
const char* to_string(Direction direction);
// SPELLED: RIGHT moves UP (the walker's right arm: x continues x - 1), LEFT moves DOWN (x
// continues x + 1); REVERSE_COMPLEMENT the other way round
Direction direction(Arm arm, Orientation orientation);

// The record end of a chain that no record bounds (no record mapping): the last coordinate
// a chain moving in |direction| could reach
constexpr Coord unbounded(Direction direction) {
    return direction == Direction::UP ? std::numeric_limits<Coord>::max() : 0;
}

// Consecutive column coordinates [first, last] of one label in one row
struct CoordRun {
    Coord first = 0;
    Coord last = 0;
    bool operator==(const CoordRun &other) const {
        return first == other.first && last == other.last;
    }
};

// Consecutive chains of one label, all in one record: their current coordinates are first,
// first + 1, .., last, and they may move on to |record_end| (inclusive) — the record's last
// coordinate for chains moving UP, its first for chains moving DOWN
struct ChainRun {
    Coord first = 0;
    Coord last = 0;
    Coord record_end = 0;
    bool operator==(const ChainRun &other) const {
        return first == other.first && last == other.last && record_end == other.record_end;
    }
};

// The first and last column coordinates of a record (CoordToHeader::sequence_range)
struct RecordRange {
    Coord first = 0;
    Coord last = 0;
};
// The record of a label's coordinate. An empty function: no record mapping, chains unbounded
using RecordOf = std::function<RecordRange(LabelId label, Coord coord)>;

static constexpr uint64_t kPollIterations = 4096;

/**
 * The caller's deadline for long calls: |expired| is read before every kPollIterations rounds
 * of their loops (counted across calls in |since|), and a call that finds it true stops where
 * it is, its output incomplete. Stopped stays set (every later call stops at once) until the
 * caller clears it. Without |expired| the rounds are counted and nothing stops.
 */
struct Clock {
    std::function<bool()> expired;
    uint64_t since = 0;         // rounds since the last reading
    uint64_t readings = 0;
    bool stopped = false;

    // before one round of work: false when it must not be done
    bool tick() {
        if (stopped)
            return false;
        if (++since < kPollIterations)
            return true;
        since = 0;
        ++readings;
        stopped = expired && expired();
        return !stopped;
    }
};

// Each appends to |out| and returns its units (on a stop by |clock|: of the work done).

// a ∩ b of two ascending lists of distinct label ids; units |a| + |b|
uint64_t intersect_labels(const LabelId *a, size_t na, const LabelId *b, size_t nb,
                          std::vector<LabelId> *out, Clock *clock = nullptr);

// the maximal runs of n ascending, distinct coordinates (std::invalid_argument otherwise);
// units n
uint64_t coordinate_runs(const Coord *coords, size_t n, std::vector<CoordRun> *out);

// the chains of |label| that start at the coordinates of runs[0, n), split at the records
// |record_of| says they lie in, moving in |direction|; units n + the record lookups
uint64_t open_chains(Direction direction, LabelId label, const CoordRun *runs, size_t n,
                     const RecordOf &record_of, std::vector<ChainRun> *out,
                     Clock *clock = nullptr);

// chains[0, n) after one step in |direction|: each chain moved by one, unless that passes its
// record's end, and kept where the next k-mer's row runs[0, m) holds the coordinate it moved
// to (the walker's rule, one merge of the two run lists); units n + m
uint64_t continue_chains(Direction direction, const ChainRun *chains, size_t n,
                         const CoordRun *row, size_t m, std::vector<ChainRun> *out,
                         Clock *clock = nullptr);

/**
 * One k-mer's row as the steps read it: its labels (ids ascending, distinct) and, in a row
 * with coordinates, each label's coordinates as maximal runs. Built once per row and kept with
 * it, so that every step through the k-mer reuses the runs. It does not know which k-mer it
 * belongs to or in which orientation it is used: the frame given it is told (Frame::step).
 */
class RowRuns {
  public:
    explicit RowRuns(bool with_coords = false) : with_coords_(with_coords) {}

    void clear(bool with_coords);
    // appends |label| (above every label added before) with its |coords| (ascending,
    // distinct, at least one; read only in a row with coordinates), std::invalid_argument
    // otherwise; units 1 + n
    uint64_t add(LabelId label, const Coord *coords = nullptr, size_t n = 0);
    // the row of hits sorted by label, each with |label| and ascending |coords| (as
    // LabelQuery::NodeHits are); units 1 per hit and per coordinate
    template <class Hits>
    uint64_t assign(const Hits &hits, bool with_coords) {
        clear(with_coords);
        uint64_t units = 0;
        for (const auto &h : hits) {
            units += with_coords ? add(h.label, h.coords.data(), h.coords.size())
                                 : add(h.label);
        }
        return units;
    }

    bool with_coords() const { return with_coords_; }
    const std::vector<LabelId>& labels() const { return labels_; }
    size_t num_labels() const { return labels_.size(); }
    // the coordinate runs of labels()[i] (a row with coordinates)
    const CoordRun* runs(size_t i) const { return runs_.data() + run_begin_[i]; }
    size_t num_runs(size_t i) const { return run_begin_[i + 1] - run_begin_[i]; }
    size_t total_runs() const { return runs_.size(); }

    uint64_t bytes() const { return model_bytes(with_coords_, labels_.size(), runs_.size()); }
    // what a row of |labels| labels with |runs| runs holds (a bound before it is built: its
    // coordinates for |runs|)
    static uint64_t model_bytes(bool with_coords, size_t labels, size_t runs);

  private:
    bool with_coords_;
    std::vector<LabelId> labels_;
    std::vector<uint32_t> run_begin_ = { 0 };   // with coordinates: labels_.size() + 1 offsets
    std::vector<CoordRun> runs_;
};

// What a step will cost at most: its units, and the bytes of the frame it builds
struct Price {
    uint64_t units = 0;
    uint64_t bytes = 0;
};

/**
 * The chains of a TRACE anchor: its labels' chain runs as opened, which the frames stepped
 * from the anchor share, each holding one bit per chain. The chains are numbered across the
 * anchor, label by label, run by run, coordinate by coordinate: label i's chains are
 * chain_begin[i] .. chain_begin[i + 1] - 1, its runs runs[run_begin[i] .. run_begin[i + 1]).
 * Owned jointly (std::shared_ptr) by the anchor frame and the frames below it, so that a
 * frame may be moved, swapped or stepped into while other frames stand on its anchor (a
 * stack that grows as a std::vector, two buffers used in turn).
 */
struct AnchorChains {
    std::vector<uint32_t> run_begin;        // labels + 1 offsets into |runs|
    std::vector<uint64_t> chain_begin;      // labels + 1: the number of each label's first chain
    std::vector<ChainRun> runs;

    uint64_t num_chains() const { return chain_begin.empty() ? 0 : chain_begin.back(); }
};

/**
 * The support of one walk in one orientation at one level, as it stands after the walk's
 * last step: the labels carrying every k-mer so far (ascending ids) and, at the TRACE level,
 * each label's chains (none when no record holds the walk whole: the label is then carried,
 * not verified). Opened at the walk's first k-mer (its anchor) and stepped one base at a time
 * on the frame's arm into another frame (the caller's stack reuses their buffers). The anchor
 * holds the chain runs; a frame below it holds which of the anchor's chains survive, one bit
 * each, and per label its index among the anchor's labels and its number of chain runs.
 */
class Frame {
  public:
    /**
     * The frame of the walk spelling |kmer|, on |arm|, in |orientation|, at |level|. |row| is
     * the row of |row_kmer|, which must be |kmer| (SPELLED) or its reverse complement
     * (REVERSE_COMPLEMENT), else std::invalid_argument. TRACE needs a row with coordinates;
     * its chains start at the row's coordinates, each bounded by its record (|record_of|).
     * Returns the units. Stopped by |clock|: the frame is incomplete (complete() false).
     */
    uint64_t open(Orientation orientation, Arm arm, Support level, std::string_view kmer,
                  std::string_view row_kmer, const RowRuns &row,
                  const RecordOf &record_of = RecordOf(), Clock *clock = nullptr);

    // An upper bound of the units of step(.., |row|, ..) and of the bytes of the frame it
    // builds, known before it (sizes only)
    Price price(const RowRuns &row) const;

    /**
     * Into |next|: the walk extended by |base| on the frame's arm (its new k-mer: kmer()
     * without its first base followed by |base| on the right arm, |base| followed by kmer()
     * without its last base on the left), narrowed by |row|, the row of |row_kmer|, which
     * must be that new k-mer (SPELLED) or its reverse complement (REVERSE_COMPLEMENT), else
     * std::invalid_argument. Returns the units. Stopped by |clock|: |next| is incomplete.
     */
    uint64_t step(char base, std::string_view row_kmer, const RowRuns &row, Frame *next,
                  Clock *clock = nullptr) const;

    Orientation orientation() const { return orientation_; }
    Arm arm() const { return arm_; }
    Support level() const { return level_; }
    Direction direction() const { return support_step::direction(arm_, orientation_); }
    // the walk's current k-mer as spelled: its last on the right arm, its first on the left
    const std::string& kmer() const { return kmer_; }
    // the steps since the anchor: the walk spells depth() + 1 k-mers
    uint32_t depth() const { return depth_; }
    // false after a step or an opening stopped by its clock: not a frame to use or step
    bool complete() const { return complete_; }

    const std::vector<LabelId>& labels() const { return labels_; }
    // TRACE: the chain runs of labels()[i] (the maximal runs of consecutive chains of one
    // record that carry the walk), and all labels' together; 0 at the KMER level
    size_t num_chain_runs(size_t i) const;
    size_t total_chain_runs() const { return runs_; }
    // TRACE: the chains of labels()[i] that carry the walk (below the anchor: its set bits,
    // counted word by word)
    uint64_t num_alive(size_t i) const;
    // TRACE: the anchor's chains, shared by the frames stepped from it (null at the KMER level)
    const AnchorChains* anchor() const { return anchor_.get(); }
    // whether some label supports the walk at the frame's level (TRACE: some chain survives)
    bool supported() const;
    // TRACE: the first coordinates of the record occurrences of the frame's walk (SPELLED) or
    // of its reverse complement (REVERSE_COMPLEMENT) in labels()[i]'s records, as runs: each
    // chain's lowest coordinate along the walk; appended to |out|, units its runs
    uint64_t starts(size_t i, std::vector<CoordRun> *out) const;

    uint64_t bytes() const;
    // what an anchor frame of |labels| labels and |runs| chain runs holds, its k-mer of k
    // bases with it
    static uint64_t anchor_bytes(Support level, size_t labels, size_t runs, size_t k);
    // what a frame below an anchor of |chains| chains holds with |labels| labels
    static uint64_t step_bytes(Support level, size_t labels, uint64_t chains, size_t k);

  private:
    // a label of a frame below the anchor: the same label's index among the anchor's labels,
    // and the frame's chain runs of it
    struct LabelRef {
        uint32_t anchor = 0;
        uint32_t runs = 0;
    };

    Orientation orientation_ = Orientation::SPELLED;
    Arm arm_ = Arm::RIGHT;
    Support level_ = Support::KMER;
    bool complete_ = false;
    uint32_t depth_ = 0;
    std::string kmer_;
    std::vector<LabelId> labels_;
    uint64_t runs_ = 0;                         // TRACE: the chain runs of all labels
    std::shared_ptr<AnchorChains> anchor_;      // TRACE: built by open(), shared by step()
    std::vector<LabelRef> refs_;                // TRACE below the anchor: one per label
    std::vector<uint64_t> bits_;                // TRACE below the anchor: one per anchor chain
};

// How a label supports a walk in one orientation
enum class Held : uint8_t {
    NONE = 0,
    LABEL_INTERSECTION = 1,     // on every k-mer of the walk in that orientation
    RECORD_VERIFIED = 2,        // a record holds the walk in that orientation whole (TRACE)
};
const char* to_string(Held held);

struct LabelSupport {
    LabelId label = 0;
    Held spelled = Held::NONE;
    Held reverse_complement = Held::NONE;
    bool operator==(const LabelSupport &other) const {
        return label == other.label && spelled == other.spelled
            && reverse_complement == other.reverse_complement;
    }
};

/**
 * The one place where the two orientations of a walk meet, at its end: every label of either
 * frame, ascending, with how it supports the walk as spelled and as its reverse complement —
 * a label is never supported by some k-mers of one strand and the rest of the other. The
 * frames must be complete and of one walk (the same k-mer, arm, level and depth), |spelled|
 * SPELLED and |reverse_complement| REVERSE_COMPLEMENT, else std::invalid_argument. Appended to
 * |out|; units the labels of both.
 */
uint64_t combine(const Frame &spelled, const Frame &reverse_complement,
                 std::vector<LabelSupport> *out);

} // namespace support_step

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_SUPPORT_STEP_HPP__
