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
    class Iterator;
    explicit CircularBuffer(size_t size)
        : buf_(std::unique_ptr<T[]>(new T[size])), size_(size) {}

    void push_back(T item) {
        buf_[back_] = item;

        if (full_) { // throw away oldest element
            front_ = (front_ + 1) % size_;
        }

        back_ = (back_ + 1) % size_;

        full_ = back_ == front_;
    }

    T pop_front() {
        if (empty()) {
            return T();
        }

        T val = buf_[front_];
        full_ = false;
        front_ = (front_ + 1) % size_;

        return val;
    }

    void reset() {
        back_ = front_;
        full_ = false;
    }

    bool empty() const { return (!full_ && (back_ == front_)); }

    bool full() const { return full_; }

    size_t capacity() const { return size_; }

    size_t size() const {
        size_t size = size_;

        if (!full_) {
            if (back_ >= front_) {
                size = back_ - front_;
            } else {
                size = size_ + back_ - front_;
            }
        }

        return size;
    }

    Iterator rbegin() { return Iterator(this, (size_ + back_ - 1) % size_); }

    friend
    std::ostream& operator<<(std::ostream& o, const CircularBuffer<T> &buf);

  private:
    std::unique_ptr<T[]> buf_;
    size_t back_ = 0;
    size_t front_ = 0;
    const size_t size_;
    bool full_ = 0;
};

template <typename T>
class CircularBuffer<T>::Iterator {
  public:
    explicit Iterator(CircularBuffer<T> *parent, size_t idx)
        : parent_(parent), idx_(idx) {}

    T& operator*()  { return parent_->buf_[idx_]; }

    Iterator& operator--() {
        idx_ = idx_ > 0 ? idx_ - 1 : parent_->size() - 1;
        return *this;
    }

  private:
    CircularBuffer<T> *parent_;
    size_t idx_;
};

} // namespace common
} // namespace mg
