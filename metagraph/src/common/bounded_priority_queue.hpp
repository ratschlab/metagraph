#ifndef __BOUNDED_PRIORITY_QUEUE_HPP__
#define __BOUNDED_PRIORITY_QUEUE_HPP__

#include <queue>
#include <vector>
#include <functional>
#include <iterator>
#include <cmath>
#include <cassert>

template <class T, class Container = std::vector<T>,
          class Compare = std::less<typename Container::value_type>>
class BoundedPriorityQueue : public std::priority_queue<T, Container, Compare> {
  public:
    BoundedPriorityQueue() = delete;
    BoundedPriorityQueue(size_t size) : max_size(size) { assert(size > 0); }

    // Push the value in the queue and remove the element with lower priority
    // in case the size of queue was about to exceed max_size.
    void push(T&& value, std::function<void(T&)>callback = [](T&){}) {
        if (this->size() < max_size) {
            this->std::priority_queue<T, Container, Compare>::push(std::move(value));
        }
        else {
            // Search for the lowes priority element among the leaves of the heap.
            auto mid = std::begin(this->c) + std::floor((this->size() - 1) / 2.0);
            auto end = std::end(this->c);
            auto min_element = std::min_element(mid, end);
            // Update the true score of the least scoring element.
            while (*min_element < value) {
                callback(*min_element);
                if (*min_element < value)
                    break;
                std::push_heap(std::begin(this->c), min_element + 1, this->comp);
                auto update_min_element = std::min_element(mid, end);
                if (update_min_element == min_element)
                    break;
                min_element = update_min_element;

            }
            if (value <= *min_element)
                callback(value);
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
    void push(T& value, std::function<void(T&)>callback = [](T&){}) {
        if (this->size() < max_size) {
            this->std::priority_queue<T, Container, Compare>::push(value);
        }
        else {
            // Search for the lowes priority element among the leaves of the heap.
            auto mid = std::begin(this->c) + std::floor((this->size() - 1) / 2.0);
            auto end = std::end(this->c);
            auto min_element = std::min_element(mid, end);
            // Update the true score of the least scoring element.
            while (*min_element < value) {
                callback(*min_element);
                if (*min_element < value)
                    break;
                std::push_heap(std::begin(this->c), min_element + 1, this->comp);
                auto update_min_element = std::min_element(mid, end);
                if (update_min_element == min_element)
                    break;
                min_element = update_min_element;

            }
            if (value <= *min_element)
                callback(value);
            // If the new value has a value greater than the smallest value in the heap,
            // replace it and heapify. Otherwise, the new value is not among the top priority
            // elements. So ignore it.
            if (*min_element < value) {
                *min_element = value;
                std::push_heap(std::begin(this->c), min_element + 1, this->comp);
            }
        }
    }

    T back() {
        assert(this->size() > 0);
        auto mid = std::begin(this->c) + std::floor((this->size() - 1) / 2.0);
        auto end = std::end(this->c);
        auto min_element = std::min_element(mid, end);
        return *(min_element);
    }

  private:
    size_t max_size;
};

#endif //  __BOUNDED_PRIORITY_QUEUE_HPP__
