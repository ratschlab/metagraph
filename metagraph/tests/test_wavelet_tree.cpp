#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include "common/vectors/wavelet_tree.hpp"
#include "common/threads/threading.hpp"


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bit_vector_dump_test";


template <typename Bitmap>
class WaveletTreeTest : public ::testing::Test { };

typedef ::testing::Types<wavelet_tree_stat,
                         wavelet_tree_fast,
                         wavelet_tree_dyn,
                         wavelet_tree_small>
        WaveletTreeTypes;

TYPED_TEST_SUITE(WaveletTreeTest, WaveletTreeTypes);


void test_next(const wavelet_tree &vector) {
    ASSERT_DEBUG_DEATH(vector.next(0, 1 << vector.logsigma()), "");
    ASSERT_DEBUG_DEATH(vector.next(0, 10 << vector.logsigma()), "");

    if (vector.size() == 0)
        return;

    EXPECT_EQ(0u, vector.next(0, vector[0]));
    EXPECT_EQ(vector.size() / 5, vector.next(vector.size() / 5, vector[vector.size() / 5]));
    EXPECT_EQ(vector.size() / 2, vector.next(vector.size() / 2, vector[vector.size() / 2]));
    EXPECT_EQ(vector.size()*2/3, vector.next(vector.size()*2/3, vector[vector.size()*2/3]));
    EXPECT_EQ(vector.size() - 1, vector.next(vector.size() - 1, vector[vector.size() - 1]));

    for (uint64_t c = 0; c < (1llu << vector.logsigma()); ++c) {
        for (uint64_t i : { uint64_t(0),
                            vector.size() / 5,
                            vector.size() / 2,
                            vector.size() * 2 / 3,
                            vector.size() - 1 }) {
            auto next = vector.next(i, c);

            EXPECT_TRUE(next >= i);

            if (next == vector.size()) {
                EXPECT_TRUE(vector[i] != c);
            } else {
                if (vector[i] == c) {
                    ASSERT_EQ(i, next);
                }
                EXPECT_TRUE(vector[next] == c);
            }
        }
    }
}

void test_prev(const wavelet_tree &vector) {
    ASSERT_DEBUG_DEATH(vector.prev(vector.size(), 1 << vector.logsigma()), "");
    ASSERT_DEBUG_DEATH(vector.prev(vector.size(), 10 << vector.logsigma()), "");

    if (vector.size() == 0)
        return;

    EXPECT_EQ(0u, vector.prev(0, vector[0]));
    EXPECT_EQ(vector.size() / 5, vector.prev(vector.size() / 5, vector[vector.size() / 5]));
    EXPECT_EQ(vector.size() / 2, vector.prev(vector.size() / 2, vector[vector.size() / 2]));
    EXPECT_EQ(vector.size() * 2 / 3, vector.prev(vector.size() * 2 / 3, vector[vector.size() * 2 / 3]));
    EXPECT_EQ(vector.size() - 1, vector.prev(vector.size() - 1, vector[vector.size() - 1]));

    for (uint64_t c = 0; c < (1llu << vector.logsigma()); ++c) {
        for (uint64_t i : { uint64_t(0),
                            vector.size() / 5,
                            vector.size() / 2,
                            vector.size() * 2 / 3,
                            vector.size() - 1 }) {
            auto prev = vector.prev(i, c);


            if (prev == vector.size()) {
                EXPECT_TRUE(vector[i] != c);
            } else {
                EXPECT_TRUE(prev <= i);
                if (vector[i] == c) {
                    ASSERT_EQ(i, prev);
                }
                EXPECT_TRUE(vector[prev] == c);
            }
        }
    }
}

void reference_based_test(const wavelet_tree &vector,
                          const std::vector<uint64_t> &reference) {
    auto int_vector = vector.to_vector();
    ASSERT_TRUE(std::equal(int_vector.begin(), int_vector.end(), reference.begin()));

    for (uint64_t c = 0; c < (1ull << vector.logsigma()); ++c) {
        uint64_t max_rank = std::count(reference.begin(), reference.end(), c);
        ASSERT_EQ(max_rank, vector.count(c));

        ASSERT_DEBUG_DEATH(vector.select(c, 0), "");
        ASSERT_DEBUG_DEATH(vector.select(c, max_rank + 1), "");

        for (size_t i : {1, 2, 10, 100, 1000}) {
            EXPECT_EQ(max_rank, vector.rank(c, vector.size() - 1 + i - 1));
        }

        for (size_t i = 1; i <= max_rank; ++i) {
            EXPECT_EQ(i, vector.rank(c, vector.select(c, i)));
        }

        EXPECT_EQ(vector[0] == c, vector.rank(c, 0) == 1);
        EXPECT_EQ(vector[0], reference[0]);

        for (size_t i = 1; i < vector.size(); ++i) {
            EXPECT_EQ(vector[i] == c, vector.rank(c, i) != vector.rank(c, i - 1));
            EXPECT_EQ(vector[i], reference[i]);
        }
    }

    for (uint64_t i = 0; i < vector.size(); ++i) {
        EXPECT_EQ(std::make_pair(vector.rank(reference[i], i), reference[i]),
                  vector.inverse_select(i));
    }

    test_next(vector);
    test_prev(vector);
}

TYPED_TEST(WaveletTreeTest, Queries) {
    wavelet_tree *vector = new TypeParam(4);
    ASSERT_TRUE(vector);
    delete vector;

    vector = new TypeParam(4, std::vector<int>(10, 0));
    ASSERT_TRUE(vector);

    EXPECT_EQ(10u, vector->size());
    ASSERT_DEBUG_DEATH((*vector)[vector->size()], "");
    ASSERT_DEBUG_DEATH((*vector)[vector->size() + 1], "");

    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(0u, (*vector)[i]);
        EXPECT_EQ(i + 1, vector->rank(0, i));
        EXPECT_EQ(0u, vector->rank(1, i));
        EXPECT_EQ(0u, vector->rank(2, i));
        EXPECT_EQ(i, vector->select(0, i + 1));
        ASSERT_DEBUG_DEATH(vector->select(1, i + 1), "");
    }
    EXPECT_EQ(0u, vector->rank(2, 0));
    EXPECT_EQ(0u, vector->rank(2, 1'000));
    EXPECT_EQ(10u, vector->rank(0, 1'000));
    ASSERT_DEBUG_DEATH(vector->select(1, 1'000), "");
    ASSERT_DEBUG_DEATH(vector->select(1, 0), "");
    delete vector;

    vector = new TypeParam(4, std::vector<int>(10, 2));
    ASSERT_TRUE(vector);

    EXPECT_EQ(10u, vector->size());
    ASSERT_DEBUG_DEATH((*vector)[vector->size()], "");
    ASSERT_DEBUG_DEATH((*vector)[vector->size() + 1], "");

    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(2u, (*vector)[i]);
        EXPECT_EQ(i + 1, vector->rank(2, i));
        EXPECT_EQ(0u, vector->rank(1, i));
        EXPECT_EQ(0u, vector->rank(3, i));
        EXPECT_EQ(i, vector->select(2, i + 1));
        ASSERT_DEBUG_DEATH(vector->select(1, i + 1), "");
    }
    EXPECT_EQ(0u, vector->rank(0, 0));
    EXPECT_EQ(0u, vector->rank(0, 1'000));
    EXPECT_EQ(10u, vector->rank(2, 1'000));
    ASSERT_DEBUG_DEATH(vector->select(3, 1'000), "");
    ASSERT_DEBUG_DEATH(vector->select(1, 0), "");
    delete vector;

    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    vector = new TypeParam(4, numbers);
    ASSERT_TRUE(vector);
    reference_based_test(*vector, numbers);

    delete vector;
}


TEST(wavelet_tree_dyn, Set) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree_dyn vector(4, numbers);

    reference_based_test(vector, numbers);

    for (size_t i = 0; i < numbers.size(); ++i) {
        uint64_t value = numbers.at(i);

        numbers.at(i) = 1;
        vector.set(i, 1);
        reference_based_test(vector, numbers);

        numbers.at(i) = 0;
        vector.set(i, 0);
        reference_based_test(vector, numbers);

        numbers.at(i) = value;
        vector.set(i, value);
    }
}


void test_wavelet_tree_ins_del(wavelet_tree_dyn *vector,
                               std::vector<uint64_t> *numbers) {
    reference_based_test(*vector, *numbers);

    for (size_t i = 0; i < numbers->size(); ++i) {
        numbers->insert(numbers->begin() + i, 1);
        vector->insert(i, 1);
        reference_based_test(*vector, *numbers);
        numbers->erase(numbers->begin() + i, numbers->begin() + i + 1);
        vector->remove(i);

        numbers->insert(numbers->begin() + i, 0);
        vector->insert(i, 0);
        reference_based_test(*vector, *numbers);
        numbers->erase(numbers->begin() + i, numbers->begin() + i + 1);
        vector->remove(i);
    }
}

void test_wavelet_tree_batch_ins_del(wavelet_tree_dyn *vector,
                                     const std::vector<uint64_t> &numbers) {
    assert(vector);

    vector->clear();

    // insert numbers in the given order
    for (size_t i = 0; i < numbers.size(); ++i) {
        ASSERT_EQ(i, vector->size());
        vector->insert(i, numbers[i]);
        ASSERT_EQ(i + 1, vector->size());
    }
    reference_based_test(*vector, numbers);

    vector->clear();

    // insert numbers in reversed order
    for (auto it = numbers.rbegin(); it != numbers.rend(); ++it) {
        vector->insert(0, *it);
    }
    reference_based_test(*vector, numbers);

    vector->clear();

    // insert reversed numbers, push_back
    for (auto it = numbers.rbegin(); it != numbers.rend(); ++it) {
        vector->insert(vector->size(), *it);
    }
    // erase all, pop_back
    while (vector->size()) {
        vector->remove(vector->size() - 1);
    }
    // insert all numbers from the batch, push_back
    for (auto it = numbers.begin(); it != numbers.end(); ++it) {
        vector->insert(vector->size(), *it);
    }
    reference_based_test(*vector, numbers);

    // insert reversed numbers, push_front
    for (auto it = numbers.begin(); it != numbers.end(); ++it) {
        vector->insert(0, *it);
    }
    // erase all, pop_back
    while (vector->size()) {
        vector->remove(vector->size() - 1);
    }
    // insert all numbers from the batch, push_back
    for (auto it = numbers.begin(); it != numbers.end(); ++it) {
        vector->insert(vector->size(), *it);
    }
    reference_based_test(*vector, numbers);

    // insert reversed numbers, push_back
    for (auto it = numbers.rbegin(); it != numbers.rend(); ++it) {
        vector->insert(vector->size(), *it);
    }
    // erase all, pop_front
    while (vector->size()) {
        vector->remove(0);
    }
    // insert all numbers from the batch, push_back
    for (auto it = numbers.begin(); it != numbers.end(); ++it) {
        vector->insert(vector->size(), *it);
    }
    reference_based_test(*vector, numbers);

    // insert reversed numbers, push_back
    for (auto it = numbers.rbegin(); it != numbers.rend(); ++it) {
        vector->insert(vector->size(), *it);
    }
    // erase all, pop_back
    while (vector->size()) {
        vector->remove(vector->size() - 1);
    }
    // insert all numbers from the batch, push_front
    for (auto it = numbers.rbegin(); it != numbers.rend(); ++it) {
        vector->insert(0, *it);
    }
    reference_based_test(*vector, numbers);
}

TEST(wavelet_tree_dyn, InsertDelete) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree_dyn vector(4, numbers);

    test_wavelet_tree_ins_del(&vector, &numbers);
}

TEST(wavelet_tree_dyn, InsertDeleteBatch) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1,
                                      3, 1, 0, 2, 3, 2, 2, 1,
                                      3, 2, 1, 2, 3, 3, 3, 0 };
    wavelet_tree_dyn vector(4);

    test_wavelet_tree_batch_ins_del(&vector, numbers);
}

TEST(wavelet_tree_dyn, Insert) {
    std::vector<uint64_t> numbers(1024, 0);
    auto vector = std::make_unique<wavelet_tree_dyn>(4, numbers);
    ASSERT_EQ(1024u, vector->size());
    ASSERT_EQ(0u, vector->operator[](0));
    EXPECT_EQ(vector->size(), vector->rank(0, vector->size() - 1));
    EXPECT_EQ(0u, vector->rank(3u, vector->size() - 1));

    vector->insert(965, 0);

    EXPECT_EQ(1025u, vector->size());
    EXPECT_EQ(0u, vector->operator[](0));
    EXPECT_EQ(vector->size(), vector->rank(0, vector->size() - 1));
    EXPECT_EQ(0u, vector->rank(3u, vector->size() - 1));
}


std::vector<uint64_t> to_std_vector(const sdsl::int_vector<> &int_vector) {
    return std::vector<uint64_t>(int_vector.begin(), int_vector.end());
}

TYPED_TEST(WaveletTreeTest, ToVector) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new TypeParam(2, numbers);
    ASSERT_TRUE(vector);
    EXPECT_EQ(numbers, to_std_vector(vector->to_vector()));
    delete vector;
}

TYPED_TEST(WaveletTreeTest, Serialization) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new TypeParam(4, numbers);
    ASSERT_TRUE(vector);
    std::ofstream outstream(test_dump_basename);
    vector->serialize(outstream);
    outstream.close();
    delete vector;

    vector = new TypeParam(1);
    ASSERT_TRUE(vector);
    std::ifstream instream(test_dump_basename);
    ASSERT_TRUE(vector->load(instream));

    reference_based_test(*vector, numbers);

    delete vector;
}

TYPED_TEST(WaveletTreeTest, MoveConstructor) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size());
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    TypeParam first(2, std::move(int_vector));
    TypeParam second(std::move(first));

    reference_based_test(second, numbers);
}

TYPED_TEST(WaveletTreeTest, MoveAssignment) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size());
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    TypeParam first(2, std::move(int_vector));
    TypeParam second(1);
    second = std::move(first);

    reference_based_test(second, numbers);
}

template <class FirstType, class SecondType>
void test_convert_to(const sdsl::int_vector<> &vector) {
    FirstType first(vector.width(), vector);

    ASSERT_EQ(vector.size(), first.size());
    ASSERT_EQ(vector.width(), first.logsigma());
    ASSERT_EQ(vector, first.to_vector());

    SecondType second = first.template convert_to<SecondType>();

    EXPECT_EQ(vector.size(), second.size());
    EXPECT_EQ(vector.width(), second.logsigma());
    EXPECT_EQ(vector, second.to_vector());
}

TYPED_TEST(WaveletTreeTest, Convert) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size(), 0, 2);
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    test_convert_to<TypeParam, wavelet_tree_stat>(int_vector);
    test_convert_to<TypeParam, wavelet_tree_fast>(int_vector);
    test_convert_to<TypeParam, wavelet_tree_dyn>(int_vector);
    test_convert_to<TypeParam, wavelet_tree_small>(int_vector);
}

TYPED_TEST(WaveletTreeTest, BeyondTheDNA) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    TypeParam vector(7, numbers);
    EXPECT_EQ(numbers, to_std_vector(vector.to_vector()));
}

TYPED_TEST(WaveletTreeTest, initialize_large) {
    vector<int> initial_content(1'000'000, 6);
    TypeParam vector(3, initial_content);
    EXPECT_EQ(initial_content.size(), vector.size());
}

TYPED_TEST(WaveletTreeTest, operator_eq) {
    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 10000 }) {
        for (int value : { 0, 1, 2, 5 }) {
            const auto vector = sdsl::int_vector<>(size, value);
            const TypeParam wt(3, vector);

            EXPECT_EQ(wt, wavelet_tree_stat(3, vector));
            EXPECT_EQ(wt, wavelet_tree_fast(3, vector));
            EXPECT_EQ(wt, wavelet_tree_dyn(3, vector));
            EXPECT_EQ(wt, wavelet_tree_small(3, vector));
        }
    }
}

TYPED_TEST(WaveletTreeTest, operator_neq) {
    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 10000 }) {
        for (int value : { 0, 1, 2, 5 }) {
            const auto vector = sdsl::int_vector<>(size, value);
            const TypeParam wt(3, sdsl::int_vector<>(size + 1, value));

            EXPECT_NE(wt, wavelet_tree_stat(3, vector));
            EXPECT_NE(wt, wavelet_tree_fast(3, vector));
            EXPECT_NE(wt, wavelet_tree_dyn(3, vector));
            EXPECT_NE(wt, wavelet_tree_small(3, vector));
        }
    }

    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 10000 }) {
        for (int value : { 0, 1, 2, 5 }) {
            const auto vector = sdsl::int_vector<>(size, value);
            const TypeParam wt(3, sdsl::int_vector<>(size + 1, value + 1));

            EXPECT_NE(wt, wavelet_tree_stat(3, vector));
            EXPECT_NE(wt, wavelet_tree_fast(3, vector));
            EXPECT_NE(wt, wavelet_tree_dyn(3, vector));
            EXPECT_NE(wt, wavelet_tree_small(3, vector));
        }
    }

    for (uint64_t size : { 10, 64, 120, 128, 1000, 10000, 10000 }) {
        for (int value : { 0, 1, 2, 5 }) {
            const auto vector = sdsl::int_vector<>(size, value);
            const TypeParam wt(3, sdsl::int_vector<>(size, value + 1));

            EXPECT_NE(wt, wavelet_tree_stat(3, vector));
            EXPECT_NE(wt, wavelet_tree_fast(3, vector));
            EXPECT_NE(wt, wavelet_tree_dyn(3, vector));
            EXPECT_NE(wt, wavelet_tree_small(3, vector));
        }
    }
}

TYPED_TEST(WaveletTreeTest, width) {
    for (uint8_t width : { 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 18 }) {
        const TypeParam wt(width, sdsl::int_vector<>(100, 1llu << (width - 1)));

        EXPECT_EQ(width, wt.logsigma());
    }
}

} // namespace
