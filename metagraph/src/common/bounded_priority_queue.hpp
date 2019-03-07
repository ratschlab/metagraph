#ifndef __BOUNDED_PRIORITY_QUEUE_HPP__
#define __BOUNDED_PRIORITY_QUEUE_HPP__

#include <queue>
#include <vector>
#include <functional>
#include <iterator>

template <class T, class Container = std::vector<T>,
          class Compare = std::less<typename Container::value_type>>
class BoundedPriorityQueue : public std::priority_queue<T, Container, Compare> {
  public:
    BoundedPriorityQueue() = delete;
    BoundedPriorityQueue(uint64_t size) : max_size(size) {}

    // Push the value in the queue and remove the element with lower priority
    // in case the size of queue was about to exceed max_size.
    void push(T&& value) {
        if (this->size() < max_size) {
            this->std::priority_queue<T, Container, Compare>::push(std::move(value));
        }
        else {
            // Search for the lowes priority element among the leaves of the heap.
            auto mid = std::begin(this->c) + this->size() / 2;
            auto end = std::end(this->c);
            auto min_element = std::min_element(mid, end);
            // If the new value has a value greater than the smallest value in the heap,
            // replace it and heapify. Otherwise, the new value is not among the top priority
            // elements. So ignore it.
            if (*min_element < value) {
                *min_element = std::move(value);
                std::push_heap(std::begin(this->c), min_element + 1, this->comp);
            }
        }
    }
    // Push the value in the queue and remove the element with lower priority
    // in case the size of queue was about to exceed max_size.
    void push(const T& value) {
        if (this->size() < max_size) {
            this->std::priority_queue<T, Container, Compare>::push(value);
        }
        else {
            // Search for the lowes priority element among the leaves of the heap.
            auto mid = std::begin(this->c) + this->size() / 2;
            auto end = std::end(this->c);
            auto min_element = std::min_element(mid, end);
            // If the new value has a value greater than the smallest value in the heap,
            // replace it and heapify. Otherwise, the new value is not among the top priority
            // elements. So ignore it.
            if (*min_element < value) {
                *min_element = value;
                std::push_heap(std::begin(this->c), min_element + 1, this->comp);
            }
        }
    }

  private:
    uint64_t max_size;
};

#endif //  __BOUNDED_PRIORITY_QUEUE_HPP__
