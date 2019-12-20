#pragma once

#include <memory>

namespace mg {
namespace common {

/**
 * A circular buffer of given size. When full, the oldest element gets overwritten.
 *
 * This class is NOT thread-safe.
 * For performance-critical code, consider using #RollingWindow.
 */
 //Implementation note: only the minimal operations needed by the current clients were
 // implemented. Feel free to extend if you need more features.
template <typename T>
class CircularBuffer {
  public:
    class ReverseIterator;
    using reverse_iterator = ReverseIterator;
    explicit CircularBuffer(size_t size)
        : buf_(std::unique_ptr<T[]>(new T[size])), size_(size) {}

    void push_back(T item) {
        buf_[end_] = item;

        if (full_) { // throw away oldest element
            front_ = (front_ + 1) % size_;
        }

        end_ = (end_ + 1) % size_;

        full_ = end_ == front_;
    }

    T pop_front() {
        assert(!empty() && "Attempting to pop empty buffer");
        T val = buf_[front_];
        full_ = false;
        front_ = (front_ + 1) % size_;

        return val;
    }

    void reset() {
        end_ = front_;
        full_ = false;
    }

    bool empty() const { return !full_ && (end_ == front_); }

    bool full() const { return full_; }

    size_t capacity() const { return size_; }

    size_t size() const {
        if (full_)
            return size_;

        if (end_ >= front_)
            return end_ - front_;

        return size_ + end_ - front_;
    }

    ReverseIterator rbegin() { return ReverseIterator(this, (size_ + end_ - 1) % size_); }

    friend
    std::ostream& operator<<(std::ostream& o, const CircularBuffer<T> &buf);

  private:
    std::unique_ptr<T[]> buf_;
    size_t end_ = 0;
    size_t front_ = 0;
    const size_t size_;
    bool full_ = 0;
};

template <typename T>
class CircularBuffer<T>::ReverseIterator {
  public:
    explicit ReverseIterator(CircularBuffer<T> *parent, size_t idx)
        : parent_(parent), idx_(idx) {}

    T& operator*()  { return parent_->buf_[idx_]; }

    ReverseIterator & operator++() {
        idx_ = idx_ > 0 ? idx_ - 1 : parent_->size() - 1;
        return *this;
    }

    bool at_begin() {
        return idx_ == parent_->front_;
    }

  private:
    CircularBuffer<T> *parent_;
    size_t idx_;
};

} // namespace common
} // namespace mg
