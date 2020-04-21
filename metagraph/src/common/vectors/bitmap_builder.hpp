#ifndef __BITMAP_BUILDER_HPP__
#define __BITMAP_BUILDER_HPP__

#include <cassert>
#include <functional>

#include "common/sorted_set.hpp"


/**
 * An abstract interface for a bitmap constructor.
 */
class bitmap_builder {
  public:
    virtual ~bitmap_builder() {}

    virtual void add_one(uint64_t pos) = 0;
    virtual void add_ones(const uint64_t *begin, const uint64_t *end) {
        std::for_each(begin, end, [&](uint64_t pos) { add_one(pos); });
    }

    // Data for initializing a bitmap
    struct InitializationData {
        const uint64_t size;
        const uint64_t num_set_bits;
        const std::function<void(const std::function<void(uint64_t)> &callback)> call_ones;
    };

    virtual InitializationData get_initialization_data() const = 0;
};


/**
 * A class for building a bitmap from positions of its set bits (ones).
 * The positions may be not distinct and can be passed in arbitrary order.
 */
class bitmap_builder_set : public bitmap_builder {
  public:
    /**
     * @brief      Constructs a new bitmap builder.
     *
     * @param[in]  size         bitmap size
     * @param[in]  num_threads  threads used by SortedSet for parallel processing
     * @param[in]  reserved     number of elements in buffer (grows automatically)
     */
    bitmap_builder_set(uint64_t size, size_t num_threads = 0, uint64_t reserved = 0)
          : size_(size), set_bit_positions_([](const auto &) {}, num_threads) {
        set_bit_positions_.reserve(reserved);
    }

    virtual void add_one(uint64_t pos) { set_bit_positions_.insert(&pos, &pos + 1); }
    virtual void add_ones(const uint64_t *begin, const uint64_t *end) {
        set_bit_positions_.insert(begin, end);
    }

    virtual InitializationData get_initialization_data() const {
        return { size(), num_set_bits(),
                 [&](auto callback) { call_ones(callback); } };
    }

  private:
    virtual uint64_t size() const { return size_; }
    virtual uint64_t num_set_bits() const {
        return const_cast<mg::common::SortedSet<uint64_t>&>(set_bit_positions_).data().size();
    }
    virtual void call_ones(const std::function<void(uint64_t)> &callback) const {
        for (uint64_t pos : const_cast<mg::common::SortedSet<uint64_t>&>(set_bit_positions_).data()) {
            assert(pos < size_ && "Indexes cannot be greater than bitmap's size");
            callback(pos);
        }
    }

    const uint64_t size_;
    mg::common::SortedSet<uint64_t> set_bit_positions_;
};


#endif // __BITMAP_BUILDER_HPP__
