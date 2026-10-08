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

} // namespace

uint64_t continue_chains(Direction direction, const ChainRun *chains, size_t n,
                         const CoordRun *row, size_t m, std::vector<ChainRun> *out,
                         Clock *clock) {
    // Both lists are ascending and disjoint, and moving every chain by the same one keeps
    // them so (clipping at a record's end only shortens a run), so the intersection is one
    // merge: each round consumes a chain run or a row run
    size_t i = 0, j = 0;
    Coord lo = 0, hi = 0;
    bool have = false;
    while (i < n && j < m) {
        if (clock && !clock->tick())
            return i + j;
        if (!have) {
            if (!moved(direction, chains[i], &lo, &hi)) {
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
            // a part of one chain run: its record's end stays with it
            out->push_back(ChainRun { std::max(lo, r.first), std::min(hi, r.last),
                                      chains[i].record_end });
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
    run_begin_.clear();
    runs_.clear();
    uint64_t units = labels_.size();
    if (level_ == Support::TRACE) {
        const Direction dir = direction();
        run_begin_.reserve(labels_.size() + 1);
        for (size_t i = 0; i < labels_.size(); ++i) {
            run_begin_.push_back(runs_.size());
            units += open_chains(dir, labels_[i], row.runs(i), row.num_runs(i), record_of,
                                 &runs_, clock);
            if (clock && clock->stopped)
                return units;
        }
        if (runs_.size() > std::numeric_limits<uint32_t>::max())
            throw std::length_error("support_step: a frame of more than 2^32 chain runs");
        run_begin_.push_back(runs_.size());
    }
    complete_ = true;
    return units;
}

Price Frame::price(const RowRuns &row) const {
    Price price;
    price.units = labels_.size() + row.num_labels();
    size_t runs = 0;
    if (level_ == Support::TRACE) {
        // a chain merge compares the runs of a label on both lists, and each run it outputs
        // is a part of one of them: at most both lists' runs
        runs = runs_.size() + row.total_runs();
        price.units += runs;
    }
    price.bytes = model_bytes(level_, std::min(labels_.size(), row.num_labels()), runs,
                              kmer_.size());
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
    next->run_begin_.clear();
    next->runs_.clear();

    // the label merge, and for a label on both lists (TRACE) the merge of its chains with its
    // coordinate runs
    const bool trace = level_ == Support::TRACE;
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
                next->run_begin_.push_back(next->runs_.size());
                chain_units += continue_chains(dir, chains(i), num_chains(i), row.runs(j),
                                               row.num_runs(j), &next->runs_, clock);
                if (clock && clock->stopped)
                    return i + j + chain_units;
            }
            ++i;
            ++j;
        }
    }
    if (trace) {
        if (next->runs_.size() > std::numeric_limits<uint32_t>::max())
            throw std::length_error("support_step: a frame of more than 2^32 chain runs");
        next->run_begin_.push_back(next->runs_.size());
    }
    next->complete_ = true;
    return a.size() + b.size() + chain_units;
}

bool Frame::supported() const {
    return level_ == Support::TRACE ? !runs_.empty() : !labels_.empty();
}

uint64_t Frame::starts(size_t i, std::vector<CoordRun> *out) const {
    if (level_ != Support::TRACE)
        throw std::logic_error("support_step: a label-level frame has no chains");
    // a chain moving UP stands at the walk's highest coordinate, one moving DOWN at its
    // lowest; every coordinate of a chain lies in its record
    const Coord back = direction() == Direction::UP ? depth_ : 0;
    const size_t n = num_chains(i);
    const ChainRun *runs = chains(i);
    for (size_t r = 0; r < n; ++r) {
        out->push_back(CoordRun { runs[r].first - back, runs[r].last - back });
    }
    return n;
}

uint64_t Frame::model_bytes(Support level, size_t labels, size_t runs, size_t k) {
    uint64_t bytes = sizeof(Frame) + 2 * k + 2 * labels * sizeof(LabelId);
    if (level == Support::TRACE)
        bytes += 2 * (labels + 1) * sizeof(uint32_t) + 2 * runs * sizeof(ChainRun);
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
        return f.level() == Support::TRACE && f.num_chains(i) ? Held::RECORD_VERIFIED
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
