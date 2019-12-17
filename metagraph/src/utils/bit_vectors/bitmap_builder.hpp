#ifndef __BITMAP_BUILDER_HPP__
#define __BITMAP_BUILDER_HPP__

#include <cassert>
#include <functional>

#include "bitmap.hpp"
#include "common/sorted_set.hpp"


/**
 * A class for building a bitmap from positions of its set bits (ones).
 * The positions may be not distinct and can be passed in arbitrary order.
 */
class bitmap_builder_set : public bitmap_builder {
  public:
    bitmap_builder_set(uint64_t size, size_t num_threads = 0)
          : size_(size), set_bit_positions_([](const auto &) {}, num_threads) {}

    virtual void add_one(uint64_t pos) { set_bit_positions_.insert(&pos, &pos + 1); }
    virtual void add_ones(const uint64_t *begin, const uint64_t *end) {
        set_bit_positions_.insert(begin, end);
    }

  private:
    virtual uint64_t size() const { return size_; }
    virtual uint64_t num_set_bits() const {
        return const_cast<SortedSet<uint64_t>&>(set_bit_positions_).data().size();
    }
    virtual void call_ones(const VoidCall<uint64_t> &callback) const {
        for (uint64_t pos : const_cast<SortedSet<uint64_t>&>(set_bit_positions_).data()) {
            assert(pos < size_ && "Indexes cannot be greater than bitmap's size");
            callback(pos);
        }
    }

    const uint64_t size_;
    SortedSet<uint64_t> set_bit_positions_;
};


#endif // __BITMAP_BUILDER_HPP__
