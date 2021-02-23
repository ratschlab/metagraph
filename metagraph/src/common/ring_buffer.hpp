#ifndef __RING_BUFFER_HPP__
#define __RING_BUFFER_HPP__

#include <vector>

#include <sdsl/bits.hpp>


template <typename T, class Storage = std::vector<T>>
class RingBuffer {
  public:
    explicit RingBuffer(size_t size)
          : ring_buffer_(size ? (1llu << (sdsl::bits::hi(size - 1) + 1)) : 0),
            size_(size),
            buffer_it_mask_(ring_buffer_.size() - 1),
            buffer_back_it_(buffer_it_mask_) {
        assert(ring_buffer_.size() >= size);
    }

    void reset() { buffer_back_it_ = buffer_it_mask_; }

    template <typename Iterator>
    void reset(Iterator begin) {
        buffer_back_it_ = size_ - 1;
        std::copy(begin, begin + size_, ring_buffer_.begin());
    }

    void push_back(const T &obj) {
        buffer_back_it_ = (buffer_back_it_ + 1) & buffer_it_mask_;
        assert(buffer_back_it_ < ring_buffer_.size());

        ring_buffer_[buffer_back_it_] = obj;
    }

    void push_front(const T &obj) {
        if (buffer_back_it_) {
            --buffer_back_it_;
        } else {
            buffer_back_it_ = buffer_it_mask_;
        }
        assert(buffer_back_it_ < ring_buffer_.size());

        ring_buffer_[get_front_index()] = obj;
    }

    const T& back() const { return ring_buffer_[buffer_back_it_]; }
    T& back() { return ring_buffer_[buffer_back_it_]; }

    const T& front() const { return ring_buffer_[get_front_index()]; }
    T& front() { return ring_buffer_[get_front_index()]; }

    size_t capacity() const { return size_; }

  private:
    size_t get_front_index() const {
        return (buffer_back_it_ + buffer_it_mask_ + 2 - size_) & buffer_it_mask_;
    }

    Storage ring_buffer_;
    size_t size_;
    size_t buffer_it_mask_;
    size_t buffer_back_it_;
};


#endif // __RING_BUFFER_HPP__
