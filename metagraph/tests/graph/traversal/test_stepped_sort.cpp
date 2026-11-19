#include "gtest/gtest.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "graph/traversal/stepped_sort.hpp"


// The sorts in steps of a /resolve's loops after the work (stepped_sort.hpp; lane T8-RL): the
// order std::sort and std::partial_sort leave by a total order, every element kept by an order
// with ties, the steps told about every kSortStepElements elements (none of a larger range left
// untold: the readings fall inside a sort of millions of elements), a throwing step ending the
// sort with no element lost; both partitions (by copies for trivially copyable elements, Hoare's
// scheme otherwise) and the heap sort past the partitions' depth

namespace {

using namespace mtg::graph::traversal;

// what the steps told: their sum and the largest one
struct Told {
    uint64_t sum = 0;
    uint64_t largest = 0;
    size_t calls = 0;
    std::function<void(uint64_t)> step() {
        return [this](uint64_t n) {
            sum += n;
            largest = std::max(largest, n);
            ++calls;
        };
    }
};

// the shapes of the round's starts and ends and of their patterns: random, sorted, reversed,
// three values, one value, organ pipe
std::vector<uint64_t> numbers(size_t n, int shape, std::mt19937_64 &rng) {
    std::vector<uint64_t> v(n);
    for (size_t i = 0; i < n; ++i) {
        switch (shape) {
            case 0: v[i] = rng(); break;
            case 1: v[i] = i; break;
            case 2: v[i] = n - i; break;
            case 3: v[i] = rng() % 3; break;
            case 4: v[i] = 7; break;
            default: v[i] = i % 2 ? i : n - i; break;
        }
    }
    return v;
}

TEST(SteppedSort, NumbersAreSortedAsByStdSort) {
    std::mt19937_64 rng(1);
    for (size_t n : { 0, 1, 2, 3, 10, 100, 16384, 16385, 20000, 100000, 1000003 }) {
        for (int shape = 0; shape < 6; ++shape) {
            const std::vector<uint64_t> v = numbers(n, shape, rng);
            std::vector<uint64_t> expected = v;
            std::sort(expected.begin(), expected.end());
            std::vector<uint64_t> got = v;
            Told told;
            sort_in_steps(got.begin(), got.end(), std::less<uint64_t>(), told.step());
            ASSERT_EQ(expected, got) << n << " " << shape;
            if (n > 1) {
                EXPECT_LE(n, told.sum) << n << " " << shape;
            }
            // told at least once per kSortStepElements of the elements, none of a range past a
            // piece left untold
            EXPECT_LE(told.largest, kSortPieceElements + kSortStepElements) << n << " " << shape;
            EXPECT_LE(n / (kSortPieceElements + kSortStepElements), told.calls) << n;
        }
    }
}

TEST(SteppedSort, PairsAreSortedAsByStdSort) {
    struct Run {
        uint64_t begin;
        uint64_t end;
    };
    auto by_begin = [](const Run &a, const Run &b) {
        return a.begin != b.begin ? a.begin < b.begin : a.end < b.end;
    };
    std::mt19937_64 rng(2);
    for (size_t n : { 16385, 100000, 1000003 }) {
        // random; 1,025 pieces of many labels shuffled and in label order (a max_support
        // round's); descending thirds; organ pipe; one value
        for (int shape = 0; shape < 6; ++shape) {
            std::vector<Run> v(n);
            for (size_t i = 0; i < n; ++i) {
                const uint64_t p = i % 1025;
                switch (shape) {
                    case 0: v[i] = { rng() % 10000, rng() % 500 }; break;
                    case 1:
                    case 2: v[i] = { p * 10, p * 10 + 9 }; break;
                    case 3: v[i] = { (n - i) / 3, 1 }; break;
                    case 4: v[i] = { i % 2 ? i : n - i, 0 }; break;
                    default: v[i] = { 5, 5 }; break;
                }
            }
            if (shape == 1)
                std::shuffle(v.begin(), v.end(), rng);
            std::vector<Run> expected = v;
            std::sort(expected.begin(), expected.end(), by_begin);
            Told told;
            sort_in_steps(v.begin(), v.end(), by_begin, told.step());
            for (size_t i = 0; i < n; ++i) {
                ASSERT_EQ(expected[i].begin, v[i].begin) << n << " " << shape << " " << i;
                ASSERT_EQ(expected[i].end, v[i].end) << n << " " << shape << " " << i;
            }
            EXPECT_LE(n, told.sum);
            EXPECT_LE(told.largest, kSortPieceElements + kSortStepElements);
        }
    }
}

// elements that are not trivially copyable (the grouping's candidates: Hoare's scheme)
TEST(SteppedSort, ObjectsAreSortedAsByStdSort) {
    struct Candidate {
        uint64_t begin;
        std::vector<uint32_t> labels;
    };
    auto total = [](const Candidate &a, const Candidate &b) {
        return a.begin != b.begin ? a.begin < b.begin : a.labels < b.labels;
    };
    auto ties = [](const Candidate &a, const Candidate &b) { return a.begin < b.begin; };
    std::mt19937_64 rng(3);
    for (size_t n : { 0, 1, 5, 16385, 50000, 300007 }) {
        for (int shape = 0; shape < 4; ++shape) {
            std::vector<Candidate> v(n);
            for (size_t i = 0; i < n; ++i) {
                v[i].begin = shape == 0 ? rng() : shape == 1 ? i : shape == 2 ? rng() % 5 : 3;
                v[i].labels.assign(1 + rng() % 3, static_cast<uint32_t>(i));
            }
            std::vector<Candidate> expected = v;
            std::sort(expected.begin(), expected.end(), total);
            std::vector<Candidate> got = v;
            Told told;
            sort_in_steps(got.begin(), got.end(), total, told.step());
            for (size_t i = 0; i < n; ++i) {
                ASSERT_EQ(expected[i].begin, got[i].begin) << n << " " << shape << " " << i;
                ASSERT_EQ(expected[i].labels, got[i].labels) << n << " " << shape << " " << i;
            }
            if (n > 1) {
                EXPECT_LE(n, told.sum);
            }
            // by an order with ties: sorted by it, every element kept
            std::vector<Candidate> tied = v;
            sort_in_steps(tied.begin(), tied.end(), ties, [](uint64_t) {});
            ASSERT_TRUE(std::is_sorted(tied.begin(), tied.end(), ties)) << n << " " << shape;
            std::sort(tied.begin(), tied.end(), total);
            for (size_t i = 0; i < n; ++i) {
                ASSERT_EQ(expected[i].labels, tied[i].labels) << n << " " << shape << " " << i;
            }
        }
    }
}

// a range partitioned past its depth is heap-sorted (std::sort's guard against a quadratic sort)
TEST(SteppedSort, TheHeapSortPastTheDepth) {
    std::mt19937_64 rng(4);
    for (int depth : { 0, 1, 2 }) {
        for (size_t n : { 16385, 50000, 300001 }) {
            std::vector<uint64_t> v(n);
            for (uint64_t &x : v) {
                x = rng() % 1000;
            }
            std::vector<uint64_t> expected = v;
            std::sort(expected.begin(), expected.end());
            auto comp = std::less<uint64_t>();
            Told told;
            stepped_sort_detail::sort_ranges(v.begin(), v.end(), comp, told.step(), depth);
            EXPECT_EQ(expected, v) << depth << " " << n;
            EXPECT_LE(n, told.sum);
            EXPECT_LE(told.largest, kSortPieceElements + kSortStepElements);

            struct Object {
                uint64_t x;
                std::string name;
            };
            std::vector<Object> objects;
            for (size_t i = 0; i < n; i += 7) {
                objects.push_back({ rng() % 1000, std::to_string(i) });
            }
            auto by_x = [](const Object &a, const Object &b) {
                return a.x != b.x ? a.x < b.x : a.name < b.name;
            };
            std::vector<Object> sorted = objects;
            std::sort(sorted.begin(), sorted.end(), by_x);
            stepped_sort_detail::sort_ranges(objects.begin(), objects.end(), by_x,
                                             [](uint64_t) {}, depth);
            for (size_t i = 0; i < objects.size(); ++i) {
                ASSERT_EQ(sorted[i].name, objects[i].name) << depth << " " << n << " " << i;
            }
        }
    }
}

// a step that throws ends the sort where it is: the elements permuted, none lost
TEST(SteppedSort, AThrowingStepLosesNothing) {
    std::mt19937_64 rng(5);
    std::vector<uint64_t> v(200000);
    for (uint64_t &x : v) {
        x = rng();
    }
    std::vector<uint64_t> got = v;
    int calls = 0;
    EXPECT_THROW(sort_in_steps(got.begin(), got.end(), std::less<uint64_t>(), [&](uint64_t) {
        if (++calls == 3)
            throw std::runtime_error("stop");
    }), std::runtime_error);
    EXPECT_EQ(3, calls);
    std::sort(v.begin(), v.end());
    std::sort(got.begin(), got.end());
    EXPECT_EQ(v, got);
}

TEST(SteppedSort, PartialSortIsStdPartialSorts) {
    std::mt19937_64 rng(6);
    for (size_t n : { 11, 100, 16385, 100000, 1000003 }) {
        for (int shape = 0; shape < 6; ++shape) {
            const std::vector<uint64_t> v = numbers(n, shape, rng);
            const size_t m = n / 7;
            std::vector<uint64_t> expected = v;
            std::partial_sort(expected.begin(), expected.begin() + m, expected.end());
            std::vector<uint64_t> got = v;
            Told told;
            partial_sort_in_steps(got.begin(), got.begin() + m, got.end(),
                                  std::less<uint64_t>(), told.step());
            ASSERT_TRUE(std::equal(expected.begin(), expected.begin() + m, got.begin()))
                    << n << " " << shape;
            std::sort(expected.begin() + m, expected.end());
            std::sort(got.begin() + m, got.end());
            ASSERT_EQ(expected, got) << n << " " << shape;
            // (each element offered to the heap, and each sift)
            EXPECT_LE(n - m, told.sum);
            EXPECT_LE(told.largest, kSortStepElements);
        }
    }
}

TEST(SteppedSort, TheVectorPassesInSteps) {
    std::mt19937_64 rng(7);
    for (size_t n : { 0, 1, 4096, 4097, 100000 }) {
        const std::vector<uint64_t> v = numbers(n, 3, rng);
        std::vector<uint64_t> sorted = v;
        std::sort(sorted.begin(), sorted.end());

        std::vector<uint64_t> expected = sorted;
        expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
        std::vector<uint64_t> got = sorted;
        Told unique_told;
        unique_in_steps(got, unique_told.step());
        EXPECT_EQ(expected, got) << n;
        EXPECT_EQ(n, unique_told.sum);
        EXPECT_LE(unique_told.largest, kSortStepElements);

        std::vector<uint64_t> grown = v;
        Told grown_told;
        reserve_in_steps(grown, 2 * n + 5, grown_told.step());
        EXPECT_EQ(v, grown) << n;
        EXPECT_LE(2 * n + 5, grown.capacity());
        EXPECT_EQ(n, grown_told.sum);
        EXPECT_LE(grown_told.largest, kSortStepElements);
        // a throwing step leaves the vector as it was
        if (n > kSortStepElements) {
            std::vector<uint64_t> kept = v;
            const size_t capacity = kept.capacity();
            EXPECT_THROW(reserve_in_steps(kept, 4 * n, [](uint64_t) {
                throw std::runtime_error("stop");
            }), std::runtime_error);
            EXPECT_EQ(v, kept);
            EXPECT_EQ(capacity, kept.capacity());
        }

        std::vector<uint32_t> zeros = { 1, 2, 3 };
        Told zeros_told;
        assign_zeros_in_steps(zeros, n, zeros_told.step());
        EXPECT_EQ(std::vector<uint32_t>(n, 0), zeros);
        EXPECT_EQ(n, zeros_told.sum);
        EXPECT_LE(zeros_told.largest, kSortStepElements);
    }
}

} // namespace
