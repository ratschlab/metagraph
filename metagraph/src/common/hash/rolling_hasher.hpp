#ifndef __ROLLING_HASHER_HPP__
#define __ROLLING_HASHER_HPP__

#include <memory>

#include <cyclichash.h>

#include "common/ring_buffer.hpp"


template <typename TAlphabet = uint8_t>
class RollingHash {
  public:
    using Hasher = CyclicHash<uint64_t, TAlphabet>;

    // Note: this constructor is expensive. Try to construct it once and
    // make copies of the object.
    explicit RollingHash(size_t k, uint32_t seed1 = 0, uint32_t seed2 = 1)
          : hash_(new Hasher(k, seed1, seed2, 64)),
            ring_buffer_(k) { assert(hash_); }

    RollingHash(const RollingHash &other)
          : hash_(new Hasher(*other.hash_)),
            ring_buffer_(other.ring_buffer_) {}

    RollingHash& operator=(const RollingHash &other) {
        hash_ = std::make_unique<Hasher>(*other.hash_);
        ring_buffer_ = other.ring_buffer_;

        return *this;
    }

    template <typename Iterator>
    void reset(Iterator it) {
        hash_->reset();
        ring_buffer_.reset(it);
        for (int i = 0; i < hash_->n; ++i) {
            hash_->eat(*it);
            ++it;
        }
    }

    void next(TAlphabet next_char) {
        hash_->update(ring_buffer_.front(), next_char);
        ring_buffer_.push_back(next_char);
    }

    void prev(TAlphabet prev_char) {
        hash_->reverse_update(prev_char, ring_buffer_.back());
        ring_buffer_.push_front(prev_char);
    }

    inline operator uint64_t() const { return hash_->hashvalue; }

    size_t get_k() const { return hash_->n; }

  private:
    std::unique_ptr<Hasher> hash_;
    RingBuffer<TAlphabet> ring_buffer_;
};


#endif // __ROLLING_HASHER_HPP__
