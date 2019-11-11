#ifndef __ROLLING_HASHER_HPP__
#define __ROLLING_HASHER_HPP__

#include <vector>

#include <cyclichash.h>

#include "utils.hpp"


template <typename TAlphabet = uint8_t>
class RollingHasher {
  public:
    // Note: this constructor is expensive. Try to construct it once and
    // make copies of the object.
    explicit RollingHasher(size_t k, uint32_t seed1 = 0, uint32_t seed2 = 1)
          : hash_(k, seed1, seed2, 64),
            ring_buffer_(k) {
        for (int i = 0; i < hash_.n; ++i) {
            hash_.eat(0);
        }
    }

    void reset(const TAlphabet *it) {
        hash_.reset();
        ring_buffer_.reset(it);
        for (int i = 0; i < hash_.n; ++i) {
            hash_.eat(*it);
            ++it;
        }
    }

    void next(TAlphabet next_char) {
        hash_.update(ring_buffer_.front(), next_char);
        ring_buffer_.push_back(next_char);
        assert(next_char == ring_buffer_.back());
    }

    void prev(TAlphabet prev_char) {
        hash_.reverse_update(prev_char, ring_buffer_.back());
        ring_buffer_.push_front(prev_char);
        assert(prev_char == ring_buffer_.front());
    }

    bool operator<(const RollingHasher &other) const {
        return hash_.hashvalue < other.hash_.hashvalue;
    }

    bool operator>(const RollingHasher &other) const {
        return hash_.hashvalue > other.hash_.hashvalue;
    }

    bool operator==(const RollingHasher &other) const {
        return hash_.hashvalue == other.hash_.hashvalue;
    }

    size_t get_k() const { return hash_.n; }

    uint64_t get_hash() const { return hash_.hashvalue; }

  private:
    CyclicHash<uint64_t, TAlphabet> hash_;
    utils::RingBuffer<TAlphabet> ring_buffer_;
};


template <int h,
          typename TAlphabet = uint8_t,
          class RollingHasher = ::RollingHasher<TAlphabet>>
class RollingMultiHasher {
  public:
    explicit RollingMultiHasher(size_t k) {
        static_assert(h);

        hashers_.reserve(h);
        for (size_t i = 0; i < h; ++i) {
            hashers_.emplace_back(k, i * 2, i * 2 + 1);
        }

        assert(hashers_.size() == h);
    }

    void reset(const TAlphabet *it) {
        assert(hashers_.size() == h);
        for (size_t i = 0; i < h; ++i) {
            hashers_[i].reset(it);
        }
    }

    void next(TAlphabet next_char) {
        assert(hashers_.size() == h);
        for (size_t i = 0; i < h; ++i) {
            hashers_[i].next(next_char);
        }
    }

    void prev(TAlphabet prev_char) {
        assert(hashers_.size() == h);
        for (size_t i = 0; i < h; ++i) {
            hashers_[i].prev(prev_char);
        }
    }

    bool operator<(const RollingMultiHasher &other) const {
        assert(hashers_.size() == h);
        return hashers_ < other.hashers_;
    }

    bool operator>(const RollingMultiHasher &other) const {
        assert(hashers_.size() == h);
        return hashers_ > other.hashers_;
    }

    bool operator==(const RollingMultiHasher &other) const {
        assert(hashers_.size() == h);
        return hashers_ == other.hashers_;
    }

    size_t get_k() const {
        assert(hashers_.size() == h);
        return hashers_.at(0).get_k();
    }

    template <int j>
    uint64_t get_hash() const {
        static_assert(j < h);
        assert(hashers_.size() == h);
        return hashers_.at(j).get_hash();
    }

  private:
    std::vector<RollingHasher> hashers_;
};


#endif // __ROLLING_HASHER_HPP__
