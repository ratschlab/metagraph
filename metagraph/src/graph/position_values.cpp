#include "position_values.hpp"

#include <limits>

#include <sdsl/util.hpp>
#if _USE_FOLLY
#include <folly/memory/Malloc.h>
#endif

#include "annotation/binary_matrix/base/decode_budget.hpp"


namespace mtg {
namespace graph {

using mtg::annot::matrix::alloc_bytes;
using mtg::annot::matrix::buffer_bytes;

namespace {

// a run of 32-bit positions ends at most here: a larger query is held densely from the start
constexpr size_t kMaxRunsSize = std::numeric_limits<uint32_t>::max();

// Make room for one more element in |v| within |budget| bytes: double its capacity,
// up to what |budget| holds. Return false if |budget| holds no more elements.
// The growth is charged to |account| (null: none) before it is made: the new buffer less
// the old one
template <typename T>
bool grow_within(std::vector<T> &v, size_t budget, SearchAccount *account) {
    const size_t max_capacity = budget / sizeof(T);
    if (v.size() >= max_capacity)
        return false;
    const size_t capacity = std::min(max_capacity, std::max<size_t>(2, 2 * v.capacity()));
    if (account && capacity > v.capacity()) {
        account->charge(buffer_bytes(capacity, sizeof(T))
                        - buffer_bytes(v.capacity(), sizeof(T)));
    }
    v.reserve(capacity);
    return true;
}

// the net bytes of a buffer of |to| bytes that replaces one of |from| (none if not more)
uint64_t net_bytes(uint64_t to, uint64_t from) {
    return to > from ? to - from : 0;
}

} // namespace

PresenceRuns::PresenceRuns(size_t size, SearchAccount *account)
      : size_(size), account_(account) {
    if (size_ > kMaxRunsSize)
        to_dense_form();
}

PresenceRuns::PresenceRuns(const sdsl::bit_vector &mask) : PresenceRuns(mask.size()) {
    call_runs_of_ones(mask, [&](size_t begin, size_t end) { add_run(begin, end); });
}

bool PresenceRuns::grow() {
    return grow_within(runs_, dense_bytes(size_), account_);
}

void PresenceRuns::to_dense_form() {
    // the mask's words (and sdsl's one word beyond them), the runs freed after
    if (account_) {
        account_->charge(net_bytes(alloc_bytes(dense_bytes(size_) + sizeof(uint64_t)),
                                   buffer_bytes(runs_.capacity(), sizeof(Run))));
    }
    dense_ = sdsl::bit_vector(size_, false);
    dense_form_ = true;
    for (const Run &run : runs_) {
        mark(run.begin, run.end);
    }
    std::vector<Run>().swap(runs_);
}

void PresenceRuns::mark(size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
        dense_[i] = true;
    }
}

size_t PresenceRuns::count() const {
    if (dense_form_)
        return sdsl::util::cnt_one_bits(dense_);
    size_t count = 0;
    for (const Run &run : runs_) {
        count += run.end - run.begin;
    }
    return count;
}

sdsl::bit_vector PresenceRuns::to_dense() const {
    if (dense_form_)
        return dense_;
    sdsl::bit_vector mask(size_, false);
    for (const Run &run : runs_) {
        for (size_t i = run.begin; i < run.end; ++i) {
            mask[i] = true;
        }
    }
    return mask;
}

size_t PresenceRuns::heap_bytes() const {
    return dense_form_ ? dense_bytes(size_) : runs_.capacity() * sizeof(Run);
}

bool PresenceRuns::operator==(const PresenceRuns &other) const {
    return size_ == other.size_ && to_dense() == other.to_dense();
}

CountRuns::CountRuns(size_t size, SearchAccount *account) : size_(size), account_(account) {
    if (size_ > kMaxRunsSize)
        to_dense_form();
}

CountRuns::CountRuns(const std::vector<size_t> &values) : CountRuns(values.size()) {
    for (size_t i = 0; i < values.size(); ++i) {
        add(i, values[i]);
    }
}

bool CountRuns::grow() {
    return grow_within(runs_, dense_bytes(size_), account_);
}

void CountRuns::to_dense_form() {
    // the counts, the runs freed after
    if (account_) {
        account_->charge(net_bytes(buffer_bytes(size_, sizeof(size_t)),
                                   buffer_bytes(runs_.capacity(), sizeof(Run))));
    }
    dense_.assign(size_, 0);
    dense_form_ = true;
    for (const Run &run : runs_) {
        std::fill(dense_.begin() + run.begin, dense_.begin() + run.end, run.value);
    }
    std::vector<Run>().swap(runs_);
}

size_t CountRuns::count() const {
    if (dense_form_)
        return size_ - std::count(dense_.begin(), dense_.end(), size_t(0));
    size_t count = 0;
    for (const Run &run : runs_) {
        count += run.end - run.begin;
    }
    return count;
}

std::vector<size_t> CountRuns::to_dense() const {
    if (dense_form_)
        return dense_;
    std::vector<size_t> values(size_, 0);
    for (const Run &run : runs_) {
        std::fill(values.begin() + run.begin, values.begin() + run.end, run.value);
    }
    return values;
}

size_t CountRuns::heap_bytes() const {
    return dense_form_ ? dense_.capacity() * sizeof(size_t) : runs_.capacity() * sizeof(Run);
}

bool CountRuns::operator==(const CountRuns &other) const {
    return size_ == other.size_ && to_dense() == other.to_dense();
}

PositionCoords::PositionCoords(size_t size, SearchAccount *account)
      : size_(size), account_(account) {
    if (size_ > kMaxRunsSize)
        to_dense_form();
}

PositionCoords::PositionCoords(const std::vector<Coords> &values) : PositionCoords(values.size()) {
    for (size_t i = 0; i < values.size(); ++i) {
        add(i, Coords(values[i]));
    }
}

bool PositionCoords::open(size_t position) {
    const bool new_run = runs_.empty() || runs_.back().end != position;
    // the runs and the coordinate lists share the budget: each grows within what
    // the other's capacity leaves
    const size_t budget = dense_bytes(size_);
    if ((coords_.size() == coords_.capacity()
            && !grow_within(coords_, budget - runs_.capacity() * sizeof(Run), account_))
        || (new_run && runs_.size() == runs_.capacity()
            && !grow_within(runs_, budget - coords_.capacity() * sizeof(Coords), account_))) {
        to_dense_form();
        return false;
    }
    if (new_run) {
        runs_.push_back({ static_cast<uint32_t>(position), static_cast<uint32_t>(position + 1) });
    } else {
        ++runs_.back().end;
    }
    return true;
}

void PositionCoords::to_dense_form() {
    // a list per position, the runs and the lists' vector freed after (the lists' own
    // buffers move)
    if (account_) {
        account_->charge(net_bytes(buffer_bytes(size_, sizeof(Coords)),
                                   buffer_bytes(runs_.capacity(), sizeof(Run))
                                   + buffer_bytes(coords_.capacity(), sizeof(Coords))));
    }
    dense_.resize(size_);
    dense_form_ = true;
    size_t c = 0;
    for (const Run &run : runs_) {
        for (size_t i = run.begin; i < run.end; ++i) {
            dense_[i] = std::move(coords_[c++]);
        }
    }
    assert(c == coords_.size());
    std::vector<Run>().swap(runs_);
    std::vector<Coords>().swap(coords_);
}

size_t PositionCoords::count() const {
    if (!dense_form_)
        return coords_.size();
    return std::count_if(dense_.begin(), dense_.end(),
                         [](const Coords &coords) { return !coords.empty(); });
}

std::vector<PositionCoords::Coords> PositionCoords::to_dense() const {
    if (dense_form_)
        return dense_;
    std::vector<Coords> values(size_);
    call_positions([&](size_t i, const Coords &coords) { values[i] = coords; });
    return values;
}

size_t PositionCoords::heap_bytes() const {
    return dense_form_ ? dense_.capacity() * sizeof(Coords)
                       : runs_.capacity() * sizeof(Run) + coords_.capacity() * sizeof(Coords);
}

bool PositionCoords::operator==(const PositionCoords &other) const {
    return size_ == other.size_ && to_dense() == other.to_dense();
}

// folly's small_vector<uint64_t, 2, uint32_t> (the deployed build's, folly as pinned): no
// buffer while its two elements are inline, else its capacity's (the capacity is stored
// beside the pointer, in the inline storage: no bytes in front of the buffer); it grows to the
// larger of one more and 3/2 its capacity plus one, its capacity filling the allocation
// folly::goodMallocSize gives (with jemalloc the buffer's size class, which alloc_bytes
// models). A std::vector without folly doubles
uint64_t PositionCoords::list_bytes(const Coords &coords) {
#if _USE_FOLLY
    return coords.capacity() > 2 ? alloc_bytes(coords.capacity() * sizeof(uint64_t)) : 0;
#else
    return buffer_bytes(coords.capacity(), sizeof(uint64_t));
#endif
}

uint64_t PositionCoords::push_bytes(const Coords &coords) {
    if (coords.size() < coords.capacity())
        return 0;
#if _USE_FOLLY
    const uint64_t size = std::max<uint64_t>(coords.size() + 1,
                                             3 * uint64_t(coords.capacity()) / 2 + 1);
    const uint64_t capacity = folly::goodMallocSize(size * sizeof(uint64_t)) / sizeof(uint64_t);
    return net_bytes(alloc_bytes(capacity * sizeof(uint64_t)), list_bytes(coords));
#else
    const uint64_t capacity = std::max<uint64_t>(1, 2 * uint64_t(coords.capacity()));
    return net_bytes(buffer_bytes(capacity, sizeof(uint64_t)), list_bytes(coords));
#endif
}

} // namespace graph
} // namespace mtg
