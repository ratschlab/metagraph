#include "bounded_priority_queue.hpp"

#include <gtest/gtest.h>
#include <algorithm>

const std::vector<uint64_t> values = {16, 4, 32, 8, 0, 128, 2, 3, 4, 5, 128, 0, 10};

TEST(BoundedPriorityQueue, push_all_range) {
    auto sorted_values = values;
    std::sort(std::begin(sorted_values), std::end(sorted_values), std::greater<uint64_t>());
    for(size_t limit = 1; limit < values.size() * 2; ++limit) {
       BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values)
            queue.push(std::move(value));
        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(limit);

        for (auto value: bounded_sorted_values) {
            EXPECT_EQ(value, queue.top()) << "Limit on priority queue: " << limit << std::endl;;
            queue.pop();
        }
    }
}

TEST(BoundedPriorityQueue, push_rvalue_all_range) {
    auto sorted_values = values;
    std::sort(std::begin(sorted_values), std::end(sorted_values), std::greater<uint64_t>());
    for(size_t limit = 1; limit < values.size() * 2; ++limit) {
       BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values)
            queue.push(value);
        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(limit);

        for (auto value: bounded_sorted_values) {
            EXPECT_EQ(value, queue.top()) << "Limit on priority queue: " << limit << std::endl;;
            queue.pop();
        }
    }
}
