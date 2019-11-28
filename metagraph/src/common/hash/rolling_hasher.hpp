#ifndef __ROLLING_HASHER_HPP__
#define __ROLLING_HASHER_HPP__

#include <vector>

#include <cyclichash.h>

#include "common/ring_buffer.hpp"


template <typename TAlphabet = uint8_t>
class RollingHash {
  public:
    // Note: this constructor is expensive. Try to construct it once and
    // make copies of the object.
    explicit RollingHash(size_t k, uint32_t seed1 = 0, uint32_t seed2 = 1)
          : hash_(k, seed1, seed2, 64), ring_buffer_(k) {}

    // TODO: replace with `initialize(Iterator kmer_begin)`
    void reset(const TAlphabet *kmer_begin) {
        hash_.reset();
        ring_buffer_.reset(kmer_begin);
        for (int i = 0; i < hash_.n; ++i) {
            hash_.eat(*kmer_begin);
            ++kmer_begin;
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

    bool operator<(const RollingHash &other) const {
        return hash_.hashvalue < other.hash_.hashvalue;
    }

    bool operator>(const RollingHash &other) const {
        return hash_.hashvalue > other.hash_.hashvalue;
    }

    bool operator==(const RollingHash &other) const {
        return hash_.hashvalue == other.hash_.hashvalue;
    }

    size_t get_k() const { return hash_.n; }

    operator uint64_t() const { return hash_.hashvalue; }

  private:
    CyclicHash<uint64_t, TAlphabet> hash_;
    RingBuffer<TAlphabet> ring_buffer_;
};


template <int num_hashes,
          typename TAlphabet = uint8_t,
          class RollingHash = ::RollingHash<TAlphabet>>
class RollingMultiHash {
  public:
    explicit RollingMultiHash(size_t k) {
        static_assert(num_hashes);

        hashers_.reserve(num_hashes);
        for (size_t i = 0; i < num_hashes; ++i) {
            hashers_.emplace_back(k, i * 2, i * 2 + 1);
        }

        assert(hashers_.size() == num_hashes);
    }

    void reset(const TAlphabet *it) {
        for (size_t i = 0; i < num_hashes; ++i) {
            hashers_[i].reset(it);
        }
    }

    void next(TAlphabet next_char) {
        for (size_t i = 0; i < num_hashes; ++i) {
            hashers_[i].next(next_char);
        }
    }

    void prev(TAlphabet prev_char) {
        for (size_t i = 0; i < num_hashes; ++i) {
            hashers_[i].prev(prev_char);
        }
    }

    bool operator<(const RollingMultiHash &other) const {
        return hashers_ < other.hashers_;
    }

    bool operator>(const RollingMultiHash &other) const {
        return hashers_ > other.hashers_;
    }

    bool operator==(const RollingMultiHash &other) const {
        return hashers_ == other.hashers_;
    }

    size_t get_k() const { return hashers_[0].get_k(); }

    template <int j>
    uint64_t get_hash() const {
        static_assert(j < num_hashes);
        return hashers_[j];
    }

  private:
    std::vector<RollingHash> hashers_;
};


#endif // __ROLLING_HASHER_HPP__
