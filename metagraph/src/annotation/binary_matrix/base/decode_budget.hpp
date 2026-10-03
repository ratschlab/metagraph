#ifndef __DECODE_BUDGET_HPP__
#define __DECODE_BUDGET_HPP__

#include <cassert>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <limits>
#include <vector>

#include "common/vector.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * The budget-aware decode path of the annotation library (DESIGN-traverse-graphlet.md
 * §14, "row decoding is budget-aware inside the decoder"; stage 3 of §14.1).
 *
 * A caller that opts in passes a DecodeBudget to a decode call (IRowDiff::decode_rows,
 * IRowDiff::decode_row_tuples). The decoder charges every heap buffer it allocates to the
 * budget BEFORE allocating it — the row-diff trace, each dependency row, each coordinate
 * tuple (whose size the delimiters give before it is copied), the reconstruction's
 * buffers — and releases it when the buffer is freed; a charge that does not fit refuses
 * the whole call, which then returns REFUSED with no side effect: its outputs untouched
 * and the budget's held bytes as they were at entry. The default decode functions are not
 * touched by any of this: the path is separate, single-threaded and taken only when a
 * caller asks for it.
 *
 * The charged bytes are a deterministic MODEL of the heap, never a measurement, so that
 * whether a read fits is reproducible: every buffer of n bytes at jemalloc's size class
 * (alloc_bytes()), with the containers' growth policies modelled explicitly. Tests check
 * the model against jemalloc's own peak (tests/annotation/row_diff).
 */

// Bytes a heap buffer of |n| bytes occupies: jemalloc's size classes — 16 B minimum, 16 B
// steps up to 128 B, then four classes per doubling. A pure function, so that a charge
// never depends on the allocator's state; for another allocator it is the model, not the
// truth (the tests that compare it with the heap need jemalloc).
inline uint64_t alloc_bytes(uint64_t n) {
    if (!n)
        return 0;
    if (n <= 16)
        return 16;
    if (n <= 128)
        return (n + 15) / 16 * 16;
    const uint64_t p = 64 - __builtin_clzll(n - 1);      // 2^p >= n > 2^(p-1)
    const uint64_t step = (uint64_t(1) << p) / 8;        // 4 classes in (2^(p-1), 2^p]
    return (n + step - 1) / step * step;
}

// The heap of a std::vector / folly::fbvector buffer for |n| elements of |size| bytes
inline uint64_t buffer_bytes(uint64_t n, uint64_t size) {
    return n ? alloc_bytes(n * size) : 0;
}

// The heap of a SmallVector (folly::small_vector<T, 2, uint32_t>) constructed, assigned
// or reserved for |n| elements from an empty one: none up to its two inline elements,
// otherwise at least four elements (folly grows past the inline capacity to 3 * 2 / 2 + 1)
// and the capacity it stores in front of the buffer (at most 8 B). Without folly a
// SmallVector is a std::vector.
inline uint64_t small_vector_bytes(uint64_t n, uint64_t size) {
#if _USE_FOLLY
    return n > 2 ? alloc_bytes(std::max<uint64_t>(n, 4) * size + 8) : 0;
#else
    return buffer_bytes(n, size);
#endif
}

enum class DecodeStatus {
    OK,
    // a charge did not fit (or the test hook denied it): nothing was returned, and the
    // budget holds what it held at entry
    REFUSED,
    // this annotation has no budget-aware decode path; the caller decodes by the default
    // path, whose buffers are not charged
    UNSUPPORTED,
};

/**
 * The bytes a decode call may hold at once, owned by the caller — one per request, used by
 * one thread (thread safety is the caller's, as for the rest of a request's state). It
 * accumulates: what a successful call returns stays charged until the caller releases it,
 * so consecutive calls share one allowance.
 */
class DecodeBudget {
  public:
    explicit DecodeBudget(uint64_t max_bytes = std::numeric_limits<uint64_t>::max())
          : max_(max_bytes) {}

    // Reserve |bytes| before they are allocated. false: nothing was reserved — held + bytes
    // would exceed the maximum, or the test hook denied this charge.
    bool charge(uint64_t bytes) {
        const uint64_t ordinal = charges_++;
        if ((deny && deny(ordinal)) || bytes > max_ - held_) {
            refused_need_ = bytes > std::numeric_limits<uint64_t>::max() - held_
                ? std::numeric_limits<uint64_t>::max() : held_ + bytes;
            return false;
        }
        held_ += bytes;
        peak_ = std::max(peak_, held_);
        return true;
    }
    void release(uint64_t bytes) {
        assert(bytes <= held_);
        held_ -= bytes;
    }
    // what a refused call does: back to what was held when it began
    void restore(uint64_t held_at_entry) {
        assert(held_at_entry <= held_);
        held_ = held_at_entry;
    }

    uint64_t max_bytes() const { return max_; }
    uint64_t held() const { return held_; }
    uint64_t peak() const { return peak_; }
    // the ordinal the next charge will have (the denial sweeps of the tests count them)
    uint64_t charges() const { return charges_; }
    // held + the charge that did not fit, at the last refusal
    uint64_t refused_need() const { return refused_need_; }
    // start measuring a call's own peak: returns the peak so far and restarts it at held
    // (end_peak() combines them again)
    uint64_t begin_peak() {
        const uint64_t before = peak_;
        peak_ = held_;
        return before;
    }
    void end_peak(uint64_t before) { peak_ = std::max(peak_, before); }

    // Test hook: asked at every charge with its ordinal; true refuses that charge as if it
    // did not fit. Empty by default.
    std::function<bool(uint64_t ordinal)> deny;

  private:
    uint64_t max_;
    uint64_t held_ = 0;
    uint64_t peak_ = 0;
    uint64_t charges_ = 0;
    uint64_t refused_need_ = 0;
};

/**
 * Per requested row, over its WHOLE row-diff path to the anchor: properties of the row,
 * the same in every batch it is decoded in (the path stops only at an anchor or at a row
 * decoded earlier in the same call, whose costs were computed over its own whole path).
 */
struct RowCost {
    // rows on the path besides the row itself (0: the row is an anchor)
    uint64_t dependency_rows = 0;
    // the entries those rows store: columns, and coordinates for a coordinate annotation
    uint64_t dependency_entries = 0;
    // an upper bound of the bytes decoding this row ALONE holds at its peak (the trace,
    // the stored rows of its path, the scratch, the reconstruction and the returned row)
    uint64_t demand = 0;
};

/**
 * A buffer reused by the rows of one decode call that grows only through grow(): each
 * growth is charged before the larger buffer is allocated and the smaller one released
 * after it is freed, so that the transient of the copy (both buffers) is charged too. The
 * growth policy is fixed (16 elements, then doubling), so that what a row of n entries
 * makes it hold alone is a function of n (peak_for()), which a row's demand states.
 */
template <class T>
class ChargedBuffer {
  public:
    Vector<T> data;

    // room for one more element; false: the growth did not fit (nothing changed)
    bool reserve_one(DecodeBudget &budget) {
        if (data.size() < capacity_)
            return true;
        const uint64_t next = std::max<uint64_t>(16, 2 * capacity_);
        const uint64_t bytes = buffer_bytes(next, sizeof(T));
        if (!budget.charge(bytes))
            return false;
        data.reserve(next);
        budget.release(charged_);
        charged_ = bytes;
        capacity_ = next;
        return true;
    }
    // the bytes charged for the buffer now (the caller releases them when it frees it)
    uint64_t charged() const { return charged_; }

    // the most growing an empty buffer to |n| elements holds at once: its last buffer and
    // the one before it, both alive while the elements are copied
    static uint64_t peak_for(uint64_t n) {
        uint64_t prev = 0, cap = 0;
        while (cap < n) {
            prev = cap;
            cap = std::max<uint64_t>(16, 2 * cap);
        }
        return buffer_bytes(cap, sizeof(T)) + buffer_bytes(prev, sizeof(T));
    }

  private:
    uint64_t capacity_ = 0;     // the modelled capacity (folly may round the real one up)
    uint64_t charged_ = 0;
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __DECODE_BUDGET_HPP__
