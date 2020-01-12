#include "common/bounded_priority_queue.hpp"

#include <gtest/gtest.h>
#include <algorithm>


//TODO: unit tests for different priority functions

const std::vector<uint64_t> values = {16, 4, 32, 8, 0, 128, 2, 3, 4, 5, 128, 0, 10};

TEST(BoundedPriorityQueue, push_all_range) {
    auto sorted_values = values;
    std::sort(sorted_values.begin(), sorted_values.end(), std::greater<uint64_t>());

    for (size_t limit = 1; limit < values.size() * 2; ++limit) {
        BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values) {
            queue.push(std::move(value));
        }

        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(std::min(limit, values.size()));
        ASSERT_EQ(bounded_sorted_values.size(), queue.size());

        for (auto value : bounded_sorted_values) {
            ASSERT_FALSE(queue.empty());
            EXPECT_EQ(value, queue.top()) << "Limit on priority queue: " << limit << std::endl;
            queue.pop();
        }
    }
}

TEST(BoundedPriorityQueue, push_all_range_pop) {
    auto sorted_values = values;
    std::sort(sorted_values.begin(), sorted_values.end(), std::greater<uint64_t>());

    for (size_t limit = 1; limit < values.size() * 2; ++limit) {
        BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values) {
            queue.push(std::move(value));
        }

        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(std::min(limit, values.size()));
        ASSERT_EQ(bounded_sorted_values.size(), queue.size());

        for (auto value : bounded_sorted_values) {
            ASSERT_FALSE(queue.empty());
            EXPECT_EQ(value, queue.pop_top()) << "Limit on priority queue: " << limit << std::endl;
        }
    }
}

TEST(BoundedPriorityQueue, push_rvalue_all_range) {
    auto sorted_values = values;
    std::sort(sorted_values.begin(), sorted_values.end(), std::greater<uint64_t>());

    for (size_t limit = 1; limit < values.size() * 2; ++limit) {
        BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values) {
            queue.push(value);
        }

        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(std::min(limit, values.size()));
        ASSERT_EQ(bounded_sorted_values.size(), queue.size());

        for (auto value : bounded_sorted_values) {
            ASSERT_FALSE(queue.empty());
            EXPECT_EQ(value, queue.top()) << "Limit on priority queue: " << limit << std::endl;
            queue.pop();
        }
    }
}

TEST(BoundedPriorityQueue, push_rvalue_all_range_pop) {
    auto sorted_values = values;
    std::sort(sorted_values.begin(), sorted_values.end(), std::greater<uint64_t>());

    for (size_t limit = 1; limit < values.size() * 2; ++limit) {
        BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values) {
            queue.push(value);
        }

        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(std::min(limit, values.size()));
        ASSERT_EQ(bounded_sorted_values.size(), queue.size());

        for (auto value : bounded_sorted_values) {
            ASSERT_FALSE(queue.empty());
            EXPECT_EQ(value, queue.pop_top()) << "Limit on priority queue: " << limit << std::endl;
        }
    }
}

TEST(BoundedPriorityQueue, back) {
    std::vector<uint64_t> values =                {16, 4, 32, 8, 0, 128, 2, 3, 4, 5, 128, 0, 10, 11};
    const std::vector<uint64_t> expected_values = {16, 4, 4,  4, 0, 0,   0, 2, 3, 4, 4,   4, 5,  8};
    BoundedPriorityQueue<uint64_t> queue(values.size() / 2);

    for (size_t i = 0; i < values.size(); ++i) {
        queue.push(values[i]);
        EXPECT_EQ(expected_values[i], queue.bottom()) << " i: " << i << std::endl;
    }
}
