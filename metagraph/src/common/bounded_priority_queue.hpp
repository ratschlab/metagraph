#ifndef __BOUNDED_PRIORITY_QUEUE_HPP__
#define __BOUNDED_PRIORITY_QUEUE_HPP__

#include <vector>
#include <functional>
#include <algorithm>
#include <cassert>

#include <priority_deque.hpp>


template <class T,
          class Container = std::vector<T>,
          class Compare = std::less<typename Container::value_type>>
class BoundedPriorityQueue {
    typedef boost::container::priority_deque<T, Container, Compare> MinMaxHeap;

  public:
    typedef typename MinMaxHeap::const_iterator const_iterator;

    BoundedPriorityQueue(size_t size = std::numeric_limits<size_t>::max())
          : max_size_(size), minmaxheap_() { }

    template <typename... Args>
    void emplace(Args&&... args) {
        if (minmaxheap_.size() < max_size_) {
            minmaxheap_.emplace(std::forward<Args>(args)...);
            return;
        }

        T value(std::forward<Args>(args)...);
        if (compare_(minmaxheap_.minimum(), value))
            minmaxheap_.update(minmaxheap_.begin(), value);
    }

    void push(const T &value) { emplace(value); }
    void push(T&& value) { emplace(std::move(value)); }

    void pop() { minmaxheap_.pop_maximum(); }

    const T& top() const { return minmaxheap_.maximum(); }

    // based on the implementation in priority_deque.hpp
    const_iterator top_it() const {
        auto it = minmaxheap_.begin() + 1;
        return it == minmaxheap_.end() ? minmaxheap_.begin() : it;
    }

    T pop_top() {
        T value = top();
        pop();
        return value;
    }

    const T& bottom() const { return minmaxheap_.minimum(); }
    const_iterator bottom_it() const { return minmaxheap_.begin(); }

    void update(const_iterator it, T&& value) {
        minmaxheap_.update(it, std::move(value));
    }

    size_t size() const { return minmaxheap_.size(); }
    bool empty() const { return minmaxheap_.empty(); }
    void clear() { minmaxheap_.clear(); }

  private:
    size_t max_size_;
    MinMaxHeap minmaxheap_;
    static Compare compare_;
};

#endif //  __BOUNDED_PRIORITY_QUEUE_HPP__
