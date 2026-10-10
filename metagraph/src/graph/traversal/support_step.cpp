#include "support_step.hpp"

#include <algorithm>
#include <stdexcept>

#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {
namespace traversal {
namespace support_step {

const char* to_string(Orientation orientation) {
    return orientation == Orientation::SPELLED ? "spelled" : "reverse_complement";
}

const char* to_string(Direction direction) {
    return direction == Direction::UP ? "up" : "down";
}

const char* to_string(Held held) {
    switch (held) {
        case Held::NONE: return "none";
        case Held::LABEL_INTERSECTION: return "label_intersection";
        case Held::RECORD_VERIFIED: return "record_verified";
    }
    return "unknown";
}

Direction direction(Arm arm, Orientation orientation) {
    // walking w to the right walks rc(w) to the left: its coordinates go down
    const bool up = (arm == Arm::RIGHT) == (orientation == Orientation::SPELLED);
    return up ? Direction::UP : Direction::DOWN;
}

uint64_t intersect_labels(const LabelId *a, size_t na, const LabelId *b, size_t nb,
                          std::vector<LabelId> *out, Clock *clock) {
    size_t i = 0, j = 0;
    while (i < na && j < nb) {
        if (clock && !clock->tick())
            return i + j;
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            out->push_back(a[i]);
            ++i;
            ++j;
        }
    }
    // the rest of the longer list is not compared, but counted: the units are the sizes
    return na + nb;
}

uint64_t coordinate_runs(const Coord *coords, size_t n, std::vector<CoordRun> *out) {
    for (size_t i = 0; i < n; ) {
        CoordRun run { coords[i], coords[i] };
        size_t j = i + 1;
        for (; j < n && coords[j] > run.last; ++j) {
            if (coords[j] != run.last + 1)
                break;
            run.last = coords[j];
        }
        if (j < n && coords[j] <= run.last)
            throw std::invalid_argument("support_step: the coordinates of a row are not ascending");
        out->push_back(run);
        i = j;
    }
    return n;
}

uint64_t open_chains(Direction direction, LabelId label, const CoordRun *runs, size_t n,
                     const RecordOf &record_of, std::vector<ChainRun> *out, Clock *clock) {
    uint64_t units = 0;
    for (size_t i = 0; i < n; ++i) {
        if (clock && !clock->tick())
            return units;
        ++units;
        if (!record_of) {
            out->push_back(ChainRun { runs[i].first, runs[i].last, unbounded(direction) });
            continue;
        }
        // a run of coordinates may span records (the last k-mer of one record and the first
        // of the next are consecutive coordinates of the column): one chain run per record
        for (Coord c = runs[i].first; ; ) {
            if (clock && !clock->tick())
                return units;
            ++units;
            const RecordRange record = record_of(label, c);
            if (record.first > c || record.last < c) {
                throw std::runtime_error("support_step: the record mapping does not hold a "
                                         "coordinate of its label");
            }
            const Coord last = std::min(runs[i].last, record.last);
            out->push_back(ChainRun { c, last, direction == Direction::UP ? record.last
                                                                        : record.first });
            if (last == runs[i].last)
                break;
            c = last + 1;
        }
    }
    return units;
}

namespace {

// The coordinates the chains of |run| move to, [*lo, *hi]: false when every chain of it
// stands at its record's end
inline bool moved(Direction direction, const ChainRun &run, Coord *lo, Coord *hi) {
    if (direction == Direction::UP) {
        if (run.first >= run.record_end)
            return false;
        *lo = run.first + 1;
        *hi = std::min(run.last, run.record_end - 1) + 1;
    } else {
        if (run.last <= run.record_end)
            return false;
        *lo = std::max(run.first, run.record_end + 1) - 1;
        *hi = run.last - 1;
    }
    return true;
}

// ----- bitmaps: one bit per chain, bit c in word c / 64

constexpr uint64_t kAll = ~uint64_t(0);

inline size_t words_for(uint64_t bits) {
    return (bits + 63) / 64;
}

// the first set bit of |bits| in [from, end), or end
inline uint64_t next_set(const uint64_t *bits, uint64_t from, uint64_t end) {
    if (from >= end)
        return end;
    uint64_t w = from >> 6;
    uint64_t word = bits[w] & (kAll << (from & 63));
    const uint64_t last = (end - 1) >> 6;
    while (!word) {
        if (++w > last)
            return end;
        word = bits[w];
    }
    return std::min(end, (w << 6) + static_cast<uint64_t>(__builtin_ctzll(word)));
}

// the first clear bit of |bits| in [from, end), or end
inline uint64_t next_clear(const uint64_t *bits, uint64_t from, uint64_t end) {
    if (from >= end)
        return end;
    uint64_t w = from >> 6;
    uint64_t word = ~bits[w] & (kAll << (from & 63));
    const uint64_t last = (end - 1) >> 6;
    while (!word) {
        if (++w > last)
            return end;
        word = ~bits[w];
    }
    return std::min(end, (w << 6) + static_cast<uint64_t>(__builtin_ctzll(word)));
}

// sets the bits [first, last] of |bits|
inline void set_bits(uint64_t *bits, uint64_t first, uint64_t last) {
    const uint64_t wa = first >> 6, wb = last >> 6;
    const uint64_t ma = kAll << (first & 63);
    const uint64_t mb = kAll >> (63 - (last & 63));
    if (wa == wb) {
        bits[wa] |= ma & mb;
        return;
    }
    bits[wa] |= ma;
    for (uint64_t w = wa + 1; w < wb; ++w) {
        bits[w] = kAll;
    }
    bits[wb] |= mb;
}

// the set bits of |bits| in [from, end)
inline uint64_t count_set(const uint64_t *bits, uint64_t from, uint64_t end) {
    uint64_t count = 0;
    for (uint64_t p = from; p < end; ) {
        uint64_t word = bits[p >> 6] & (kAll << (p & 63));
        const uint64_t word_end = (p | 63) + 1;
        if (end < word_end)
            word &= kAll >> (word_end - end);
        count += __builtin_popcountll(word);
        p = word_end;
    }
    return count;
}

/**
 * The chain runs of one label of a frame, read in order, each in the ANCHOR's coordinates
 * (where its chains stood when the walk began) with the number of its first chain: the
 * anchor's own runs when |bits| is null, else the maximal runs of set bits inside each
 * anchor run — the chains that still carry the walk, consecutive in one record.
 */
class ChainRuns {
  public:
    ChainRuns(const AnchorChains &anchor, size_t label, const uint64_t *bits)
          : run_(anchor.runs.data() + anchor.run_begin[label]),
            end_(anchor.runs.data() + anchor.run_begin[label + 1]),
            chain_(anchor.chain_begin[label]), pos_(chain_), bits_(bits) {}

    bool next(ChainRun *run, uint64_t *chain) {
        while (run_ != end_) {
            const uint64_t end = chain_ + (run_->last - run_->first + 1);
            if (!bits_) {
                *run = *run_;
                *chain = chain_;
                ++run_;
                chain_ = end;
                return true;
            }
            const uint64_t p = next_set(bits_, pos_, end);
            if (p == end) {
                ++run_;
                chain_ = pos_ = end;
                continue;
            }
            const uint64_t q = next_clear(bits_, p + 1, end);
            *run = ChainRun { run_->first + (p - chain_), run_->first + (q - 1 - chain_),
                              run_->record_end };
            *chain = p;
            pos_ = q;
            if (q == end) {
                ++run_;
                chain_ = end;
            }
            return true;
        }
        return false;
    }

  private:
    const ChainRun *run_;
    const ChainRun *end_;
    uint64_t chain_;            // the number of run_'s first chain
    uint64_t pos_;              // the next bit to read, in run_
    const uint64_t *bits_;
};

/**
 * One step of n chain runs, read from |chains| (next(ChainRun*, uint64_t*): the run at its
 * current coordinates and the number of its first chain), against the row runs row[0, m):
 * each chain moved by one in |direction|, unless that passes its record's end, and kept where
 * the row holds the coordinate it moved to. Both lists are ascending and disjoint, and moving
 * every chain by the same one keeps them so (clipping at a record's end only shortens a run),
 * so the intersection is one merge: each round consumes a chain run or a row run. Every
 * surviving piece of a chain run goes to |out| (piece(first, last, run, chain): the next
 * coordinates [first, last], the run it comes from and the number of its first chain).
 * Returns the units: n + m, or on a stop by |clock| the runs consumed so far.
 */
template <class Source, class Sink>
uint64_t merge_chains(Direction direction, size_t n, Source &chains,
                      const CoordRun *row, size_t m, Sink &out, Clock *clock) {
    size_t i = 0, j = 0;
    ChainRun run;
    uint64_t chain = 0;
    Coord lo = 0, hi = 0;
    bool have = false;
    while (i < n && j < m) {
        if (clock && !clock->tick())
            return i + j;
        if (!have) {
            if (!chains.next(&run, &chain))
                throw std::logic_error("support_step: a frame holds fewer chain runs than counted");
            if (!moved(direction, run, &lo, &hi)) {
                ++i;
                continue;
            }
            have = true;
        }
        const CoordRun &r = row[j];
        if (r.last < lo) {
            ++j;
        } else if (hi < r.first) {
            ++i;
            have = false;
        } else {
            out.piece(std::max(lo, r.first), std::min(hi, r.last), run, chain);
            if (hi <= r.last) {
                ++i;
                have = false;
            } else {
                ++j;
            }
        }
    }
    return n + m;
}

// chain runs from an array, at their current coordinates
struct ArrayChains {
    const ChainRun *chains;
    bool next(ChainRun *run, uint64_t *chain) {
        *run = *chains++;
        *chain = 0;
        return true;
    }
};

// the pieces as chain runs: a part of one chain run, its record's end staying with it
struct RunSink {
    std::vector<ChainRun> *out;
    void piece(Coord first, Coord last, const ChainRun &run, uint64_t) {
        out->push_back(ChainRun { first, last, run.record_end });
    }
};

} // namespace

uint64_t continue_chains(Direction direction, const ChainRun *chains, size_t n,
                         const CoordRun *row, size_t m, std::vector<ChainRun> *out,
                         Clock *clock) {
    ArrayChains source { chains };
    RunSink sink { out };
    return merge_chains(direction, n, source, row, m, sink, clock);
}

// ------------------------------------------------------------------ RowRuns

void RowRuns::clear(bool with_coords) {
    with_coords_ = with_coords;
    labels_.clear();
    run_begin_.assign(1, 0);
    runs_.clear();
}

uint64_t RowRuns::add(LabelId label, const Coord *coords, size_t n) {
    if (!labels_.empty() && label <= labels_.back())
        throw std::invalid_argument("support_step: the labels of a row are not ascending");
    if (with_coords_ && !n) {
        // a coordinate annotation lists no label without coordinates: a row read without
        // them (a LabelQuery built without coordinates) would verify nothing, silently
        throw std::invalid_argument("support_step: a row with coordinates lists a label "
                                    "without coordinates");
    }
    labels_.push_back(label);
    if (!with_coords_) {
        run_begin_.push_back(0);
        return 1;
    }
    coordinate_runs(coords, n, &runs_);
    if (runs_.size() > std::numeric_limits<uint32_t>::max())
        throw std::length_error("support_step: a row of more than 2^32 coordinate runs");
    run_begin_.push_back(runs_.size());
    return 1 + n;
}

uint64_t RowRuns::model_bytes(bool with_coords, size_t labels, size_t runs) {
    uint64_t bytes = sizeof(RowRuns) + 2 * labels * sizeof(LabelId)
                        + 2 * (labels + 1) * sizeof(uint32_t);
    if (with_coords)
        bytes += 2 * runs * sizeof(CoordRun);
    return bytes;
}

// ------------------------------------------------------------------ Frame

namespace {

// |row_kmer| is the k-mer whose row a frame of |orientation| must read where its walk spells
// |kmer|: the k-mer itself, or its reverse complement
void check_row(Orientation orientation, std::string_view kmer, std::string_view row_kmer) {
    bool ok = row_kmer.size() == kmer.size();
    if (ok && orientation == Orientation::SPELLED) {
        ok = row_kmer == kmer;
    } else if (ok) {
        const size_t k = kmer.size();
        for (size_t i = 0; ok && i < k; ++i) {
            ok = row_kmer[i] == static_cast<char>(
                    COMPL_TAB[static_cast<unsigned char>(kmer[k - 1 - i])]);
        }
    }
    if (!ok) {
        throw std::invalid_argument("support_step: a frame of orientation "
                                    + std::string(to_string(orientation)) + " at "
                                    + std::string(kmer) + " was given the row of "
                                    + std::string(row_kmer));
    }
}

} // namespace

namespace {

// |run|, read in the anchor's coordinates, where its chains stand after |depth| steps in
// |direction|
inline ChainRun at_depth(Direction direction, uint32_t depth, ChainRun run) {
    if (direction == Direction::UP) {
        run.first += depth;
        run.last += depth;
    } else {
        run.first -= depth;
        run.last -= depth;
    }
    return run;
}

// the chain runs of one label of a frame at the coordinates they stand at after |depth| steps
struct CurrentChains {
    ChainRuns runs;
    Direction direction;
    uint32_t depth;

    bool next(ChainRun *run, uint64_t *chain) {
        if (!runs.next(run, chain))
            return false;
        *run = at_depth(direction, depth, *run);
        return true;
    }
};

// the pieces of a step as bits of the next frame: the chain that moved to coordinate x stood
// at x - 1 (UP) or x + 1 (DOWN), |run.first| is where the run's first chain stood, and the
// chains of a run are numbered in the order of their coordinates
struct BitSink {
    Direction direction;
    uint64_t *bits;
    uint32_t runs = 0;

    void piece(Coord first, Coord last, const ChainRun &run, uint64_t chain) {
        const Coord from = direction == Direction::UP ? first - 1 : first + 1;
        const uint64_t offset = chain + (from - run.first);
        set_bits(bits, offset, offset + (last - first));
        ++runs;
    }
};

} // namespace

uint64_t Frame::open(Orientation orientation, Arm arm, Support level, std::string_view kmer,
                     std::string_view row_kmer, const RowRuns &row,
                     const RecordOf &record_of, Clock *clock) {
    if (kmer.empty())
        throw std::invalid_argument("support_step: an empty k-mer");
    check_row(orientation, kmer, row_kmer);
    if (level == Support::TRACE && !row.with_coords())
        throw std::invalid_argument("support_step: a trace frame needs a row with coordinates");
    orientation_ = orientation;
    arm_ = arm;
    level_ = level;
    complete_ = false;
    depth_ = 0;
    kmer_.assign(kmer.data(), kmer.size());
    labels_.assign(row.labels().begin(), row.labels().end());
    runs_ = 0;
    refs_.clear();
    bits_.clear();
    uint64_t units = labels_.size();
    if (level_ != Support::TRACE) {
        anchor_.reset();
        complete_ = true;
        return units;
    }
    // the block is reused when no frame stepped from an earlier anchor still shares it
    if (!anchor_ || anchor_.use_count() != 1)
        anchor_ = std::make_shared<AnchorChains>();
    AnchorChains &a = *anchor_;
    a.run_begin.clear();
    a.chain_begin.clear();
    a.runs.clear();
    const Direction dir = direction();
    a.run_begin.reserve(labels_.size() + 1);
    a.chain_begin.reserve(labels_.size() + 1);
    uint64_t chains = 0;
    for (size_t i = 0; i < labels_.size(); ++i) {
        a.run_begin.push_back(a.runs.size());
        a.chain_begin.push_back(chains);
        const size_t from = a.runs.size();
        units += open_chains(dir, labels_[i], row.runs(i), row.num_runs(i), record_of,
                             &a.runs, clock);
        if (clock && clock->stopped)
            return units;
        for (size_t r = from; r < a.runs.size(); ++r) {
            chains += a.runs[r].last - a.runs[r].first + 1;
        }
    }
    if (a.runs.size() > std::numeric_limits<uint32_t>::max())
        throw std::length_error("support_step: a frame of more than 2^32 chain runs");
    a.run_begin.push_back(a.runs.size());
    a.chain_begin.push_back(chains);
    runs_ = a.runs.size();
    complete_ = true;
    return units;
}

Price Frame::price(const RowRuns &row) const {
    Price price;
    price.units = labels_.size() + row.num_labels();
    if (level_ == Support::TRACE) {
        // a chain merge compares the runs of a label on both lists: at most both lists' runs
        price.units += runs_ + row.total_runs();
    }
    price.bytes = step_bytes(level_, std::min(labels_.size(), row.num_labels()),
                             anchor_ ? anchor_->num_chains() : 0, kmer_.size());
    return price;
}

uint64_t Frame::step(char base, std::string_view row_kmer, const RowRuns &row, Frame *next,
                     Clock *clock) const {
    if (!complete_)
        throw std::logic_error("support_step: an incomplete frame was stepped");
    if (next == this)
        throw std::invalid_argument("support_step: a frame stepped into itself");
    if (level_ == Support::TRACE && !row.with_coords())
        throw std::invalid_argument("support_step: a trace frame needs a row with coordinates");
    // incomplete until the merge is done: a refused row or a stop leaves nothing usable
    next->complete_ = false;
    const size_t k = kmer_.size();
    std::string &spelled = next->kmer_;
    if (arm_ == Arm::RIGHT) {
        spelled.assign(kmer_, 1, k - 1);
        spelled.push_back(base);
    } else {
        spelled.assign(1, base);
        spelled.append(kmer_, 0, k - 1);
    }
    check_row(orientation_, spelled, row_kmer);
    next->orientation_ = orientation_;
    next->arm_ = arm_;
    next->level_ = level_;
    next->depth_ = depth_ + 1;
    next->labels_.clear();
    next->runs_ = 0;
    next->refs_.clear();
    const bool trace = level_ == Support::TRACE;
    if (trace) {
        next->anchor_ = anchor_;
        next->bits_.assign(words_for(anchor_->num_chains()), uint64_t(0));
    } else {
        next->anchor_.reset();
        next->bits_.clear();
    }

    // the label merge, and for a label on both lists (TRACE) the merge of its chains with its
    // coordinate runs
    const Direction dir = direction();
    const std::vector<LabelId> &a = labels_;
    const std::vector<LabelId> &b = row.labels();
    uint64_t chain_units = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (clock && !clock->tick())
            return i + j + chain_units;
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            next->labels_.push_back(a[i]);
            if (trace) {
                const uint32_t ai = depth_ ? refs_[i].anchor : static_cast<uint32_t>(i);
                CurrentChains source { ChainRuns(*anchor_, ai, depth_ ? bits_.data() : nullptr),
                                       dir, depth_ };
                BitSink sink { dir, next->bits_.data() };
                chain_units += merge_chains(dir, num_chain_runs(i), source, row.runs(j),
                                            row.num_runs(j), sink, clock);
                next->refs_.push_back(LabelRef { ai, sink.runs });
                next->runs_ += sink.runs;
                if (clock && clock->stopped)
                    return i + j + chain_units;
            }
            ++i;
            ++j;
        }
    }
    next->complete_ = true;
    return a.size() + b.size() + chain_units;
}

size_t Frame::num_chain_runs(size_t i) const {
    if (level_ != Support::TRACE)
        return 0;
    if (depth_)
        return refs_[i].runs;
    return anchor_->run_begin[i + 1] - anchor_->run_begin[i];
}

uint64_t Frame::num_alive(size_t i) const {
    if (level_ != Support::TRACE)
        return 0;
    const size_t ai = depth_ ? refs_[i].anchor : i;
    const uint64_t from = anchor_->chain_begin[ai], end = anchor_->chain_begin[ai + 1];
    return depth_ ? count_set(bits_.data(), from, end) : end - from;
}

bool Frame::supported() const {
    return level_ == Support::TRACE ? runs_ > 0 : !labels_.empty();
}

uint64_t Frame::starts(size_t i, std::vector<CoordRun> *out) const {
    if (level_ != Support::TRACE)
        throw std::logic_error("support_step: a label-level frame has no chains");
    // the runs stand in the anchor's coordinates: a chain moving UP began at the walk's lowest
    // coordinate, one moving DOWN stands depth() below where it began; every coordinate of
    // a chain lies in its record
    const Coord back = direction() == Direction::UP ? 0 : depth_;
    ChainRuns runs(*anchor_, depth_ ? refs_[i].anchor : i, depth_ ? bits_.data() : nullptr);
    ChainRun run;
    uint64_t chain = 0, n = 0;
    while (runs.next(&run, &chain)) {
        out->push_back(CoordRun { run.first - back, run.last - back });
        ++n;
    }
    return n;
}

uint64_t Frame::bytes() const {
    if (!depth_)
        return anchor_bytes(level_, labels_.size(), runs_, kmer_.size());
    return step_bytes(level_, labels_.size(), anchor_ ? anchor_->num_chains() : 0, kmer_.size());
}

uint64_t Frame::anchor_bytes(Support level, size_t labels, size_t runs, size_t k) {
    uint64_t bytes = sizeof(Frame) + 2 * k + 2 * labels * sizeof(LabelId);
    if (level == Support::TRACE) {
        bytes += sizeof(AnchorChains) + 2 * (labels + 1) * (sizeof(uint32_t) + sizeof(uint64_t))
                    + 2 * runs * sizeof(ChainRun);
    }
    return bytes;
}

uint64_t Frame::step_bytes(Support level, size_t labels, uint64_t chains, size_t k) {
    uint64_t bytes = sizeof(Frame) + 2 * k + 2 * labels * sizeof(LabelId);
    if (level == Support::TRACE)
        bytes += 2 * labels * sizeof(LabelRef) + words_for(chains) * sizeof(uint64_t);
    return bytes;
}

// ------------------------------------------------------------------ the two orientations

uint64_t combine(const Frame &spelled, const Frame &reverse_complement,
                 std::vector<LabelSupport> *out) {
    if (spelled.orientation() != Orientation::SPELLED
            || reverse_complement.orientation() != Orientation::REVERSE_COMPLEMENT) {
        throw std::invalid_argument("support_step: combine takes a spelled and a "
                                    "reverse-complement frame");
    }
    if (!spelled.complete() || !reverse_complement.complete()
            || spelled.kmer() != reverse_complement.kmer()
            || spelled.arm() != reverse_complement.arm()
            || spelled.level() != reverse_complement.level()
            || spelled.depth() != reverse_complement.depth()) {
        throw std::invalid_argument("support_step: combine takes two complete frames of one walk");
    }
    auto held = [](const Frame &f, size_t i) {
        return f.level() == Support::TRACE && f.num_chain_runs(i) ? Held::RECORD_VERIFIED
                                                                  : Held::LABEL_INTERSECTION;
    };
    const std::vector<LabelId> &a = spelled.labels();
    const std::vector<LabelId> &b = reverse_complement.labels();
    size_t i = 0, j = 0;
    while (i < a.size() || j < b.size()) {
        LabelSupport s;
        if (j == b.size() || (i < a.size() && a[i] < b[j])) {
            s.label = a[i];
            s.spelled = held(spelled, i++);
        } else if (i == a.size() || b[j] < a[i]) {
            s.label = b[j];
            s.reverse_complement = held(reverse_complement, j++);
        } else {
            s.label = a[i];
            s.spelled = held(spelled, i++);
            s.reverse_complement = held(reverse_complement, j++);
        }
        out->push_back(s);
    }
    return a.size() + b.size();
}

} // namespace support_step
} // namespace traversal
} // namespace graph
} // namespace mtg
