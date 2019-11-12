#ifndef __DEQUE_STORAGE_HPP__
#define __DEQUE_STORAGE_HPP__

#include <deque>


/**
 * Implements interface of vector using std::deque container
 */
template <typename T>
class DequeStorage {
  public:
    typedef T value_type;
    using iterator = typename std::deque<T>::iterator;
    using const_iterator = typename std::deque<T>::const_iterator;

    template <class... Args>
    DequeStorage(Args&&... args) : deque_(std::forward<Args>(args)...),
                                   size_(deque_.size()) {}
    ~DequeStorage() {}

    DequeStorage& operator=(const DequeStorage &other) = default;
    DequeStorage& operator=(DequeStorage&& other) = default;

    inline void push_back(T&& value) {
        if (size_ == deque_.size())
            try_reserve(size_ * growth_factor, size_ + 1);

        deque_[size_++] = std::move(value);
    }

    inline void push_back(const T &value) {
        if (size_ == deque_.size())
            try_reserve(size_ * growth_factor, size_ + 1);

        deque_[size_++] = value;
    }

    template <class... Args>
    inline void emplace_back(Args&&... args) {
        if (size_ == deque_.size())
            try_reserve(size_ * growth_factor, size_ + 1);

        deque_[size_++] = T(std::forward<Args>(args)...);
    }

    inline iterator erase(iterator first, iterator last) {
        return __erase(first, last);
    }

    inline const_iterator erase(const_iterator first, const_iterator last) {
        return __erase(first, last);
    }

    inline void clear() {
        size_ = 0;
        deque_.clear();
    }

    inline void resize(size_t size) {
        reserve(size);
        size_ = size;
    }

    inline void reserve(size_t size) {
        if (size > deque_.size())
            deque_.resize(size);
    }

    inline void try_reserve(size_t size, size_t min_size = 0) {
        size = std::max(size, min_size);

        while (size > min_size) {
            try {
                reserve(size);
                return;
            } catch (const std::bad_alloc &exception) {
                size = min_size + (size - min_size) * 2 / 3;
            }
        }
        reserve(min_size);
    }

    inline size_t size() const { return size_; }

    inline size_t capacity() const { return deque_.size(); }

    inline T& operator[](size_t pos) { return deque_[pos]; }
    inline const T& operator[](size_t pos) const { return deque_[pos]; }

    inline T& at(size_t pos) {
        if (pos > size_)
            throw std::out_of_range("Out of range error");
        return deque_[pos];
    }

    inline const T& at(size_t pos) const {
        if (pos > size_)
            throw std::out_of_range("Out of range error");
        return deque_[pos];
    }

    inline T& front() { return deque_.front(); }
    inline const T& front() const { return deque_.front(); }

    inline T& back() { return deque_[size_ - 1]; }
    inline const T& back() const { return deque_[size_ - 1]; }

    inline void shrink_to_fit() {
        deque_.resize(size_);
        deque_.shrink_to_fit();
    }

    inline auto begin() { return deque_.begin(); }
    inline auto end() { return deque_.begin() + size_; }

    inline auto begin() const { return deque_.begin(); }
    inline auto end() const { return deque_.begin() + size_; }

  private:
    template <typename Iterator>
    inline Iterator __erase(Iterator first, Iterator last) {
        assert(first <= last);

        if (last == end()) {
            size_ -= (last - first);
            return end();
        }

        size_t old_size = deque_.size();
        auto it = deque_.erase(first, last);
        size_ -= old_size - deque_.size();
        return it;
    }

    static constexpr double growth_factor = 3. / 2;
    std::deque<T> deque_;
    size_t size_;
};

#endif // __DEQUE_STORAGE_HPP__
