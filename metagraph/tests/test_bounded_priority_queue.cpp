#include "bounded_priority_queue.hpp"
#include "path.hpp"

#include <gtest/gtest.h>
#include <algorithm>


const std::vector<uint64_t> values = {16, 4, 32, 8, 0, 128, 2, 3, 4, 5, 128, 0, 10};

TEST(BoundedPriorityQueue, push_all_range) {
    auto sorted_values = values;
    std::sort(std::begin(sorted_values), std::end(sorted_values), std::greater<uint64_t>());
    for (size_t limit = 1; limit < values.size() * 2; ++limit) {
       BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values)
            queue.push(std::move(value));
        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(limit);

        for (auto value: bounded_sorted_values) {
            EXPECT_EQ(value, queue.top()) << "Limit on priority queue: " << limit << std::endl;
            queue.pop();
        }
    }
}

TEST(BoundedPriorityQueue, push_rvalue_all_range) {
    auto sorted_values = values;
    std::sort(std::begin(sorted_values), std::end(sorted_values), std::greater<uint64_t>());
    for (size_t limit = 1; limit < values.size() * 2; ++limit) {
        BoundedPriorityQueue<uint64_t> queue(limit);
        for (auto value : values)
            queue.push(value);
        auto bounded_sorted_values = sorted_values;
        bounded_sorted_values.resize(limit);

        for (auto value: bounded_sorted_values) {
            EXPECT_EQ(value, queue.top()) << "Limit on priority queue: " << limit << std::endl;
            queue.pop();
        }
    }
}

TEST(BoundedPriorityQueue, back) {
    const std::vector<uint64_t> values =          {16, 4, 32, 8, 0, 128, 2, 3, 4, 5, 128, 0, 10, 11};
    const std::vector<uint64_t> expected_values = {16, 4, 4,  4, 0, 0,   0, 2, 3, 4, 4,   4, 5,  8};
    BoundedPriorityQueue<uint64_t> queue(values.size()/2);
    for (size_t i = 0; i < values.size(); ++i) {
        queue.push(values[i]);
        EXPECT_EQ(expected_values[i], queue.back()) << " i: " << i << std::endl;
    }
}

TEST(Path, cigar_update_insertion) {
    size_t k = 4;
    std::string query = "GGATTTTCAAA";
    Path<uint64_t, std::vector<std::string>> path(k, std::begin(query), std::begin(query));
    path.seed(1, {}, "GGAT", 8);
    path.extend(2, {}, 'T', 2);
    path.extend(3, {}, 'C', -2);
    path.extend(4, {}, 'A', -2);
    path.extend(5, {}, 'A', -2);
    path.extend(6, {}, 'A', 2);

    uint64_t expected_query_size = path.get_query_it() - path.get_query_begin_it() + k - 1;
    EXPECT_EQ(9ul, expected_query_size);

    path.set_cigar("5=2I4=");

    expected_query_size = path.get_query_it() - path.get_query_begin_it() + k - 1;
    EXPECT_EQ(query.size(), expected_query_size);
    EXPECT_EQ("GGATTCAAA", path.get_sequence());
}

TEST(Path, cigar_update_deletion) {
    size_t k = 4;
    std::string query = "GGATCAAA";
    Path<uint64_t, std::vector<std::string>> path(k, std::begin(query), std::begin(query));
    path.seed(1, {}, "GGAT", 8);
    path.extend(2, {}, 'T', -2);
    path.extend(2, {}, 'T', -2);
    path.extend(3, {}, 'C', -2);
    path.extend(4, {}, 'A', 2);
    path.extend(5, {}, 'A', 2);
    path.extend(6, {}, 'A', 2);

    uint64_t expected_query_size = path.get_query_it() - path.get_query_begin_it() + k - 1;
    EXPECT_EQ(10ul, expected_query_size);

    path.set_cigar("4=2D4=");

    expected_query_size = path.get_query_it() - path.get_query_begin_it() + k - 1;
    EXPECT_EQ(query.size(), expected_query_size);
    EXPECT_EQ("GGATTTCAAA", path.get_sequence());
}


TEST(Path, cigar_update_clipped) {
    size_t k = 4;
    std::string reference = "GGATTTCA";
    std::string query     = "GGATTTCAAA";
    Path<uint64_t, std::vector<std::string>> path(k, std::begin(query), std::begin(query));
    path.seed(1, {}, "GGAT", 8);
    path.extend(2, {}, 'T', -2);
    path.extend(2, {}, 'T', -2);
    path.extend(3, {}, 'C', -2);
    path.extend(6, {}, 'A', 2);

    uint64_t expected_query_size = path.get_query_it() - path.get_query_begin_it() + k - 1;
    EXPECT_EQ(8ul, expected_query_size);

    path.set_cigar("8=2S");

    expected_query_size = path.get_query_it() - path.get_query_begin_it() + k - 1;
    EXPECT_EQ(reference.size(), expected_query_size);
    EXPECT_EQ(reference, path.get_sequence());
}
