#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include "common/vectors/bit_vector.hpp"
#include "common/vectors/int_vector_algorithm.hpp"
#include "common/threads/threading.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bit_vector_dump_test";


template <typename Bitmap>
class BitVectorTest : public ::testing::Test { };

typedef ::testing::Types<bit_vector_stat,
                         bit_vector_dyn,
                         bit_vector_sd,
                         bit_vector_rrr<>,
                         bit_vector_il<>,
                         bit_vector_hyb<>,
                         bit_vector_small,
                         bit_vector_smart>
        BitVectorTypes;

TYPED_TEST_SUITE(BitVectorTest, BitVectorTypes);

template <typename Bitmap>
class BitVectorTestSelect0 : public ::testing::Test { };

typedef ::testing::Types<bit_vector_dyn,
                         bit_vector_sd,
                         bit_vector_rrr<>,
                         bit_vector_il<>,
                         bit_vector_hyb<>>
        BitVectorTypesSelect0;

TYPED_TEST_SUITE(BitVectorTestSelect0, BitVectorTypesSelect0);


void test_next_subvector(const bit_vector &vector, uint64_t idx) {
    if (vector.size() == 0)
        return;

    uint64_t count = idx ? vector.rank1(idx - 1) : 0;

    for (size_t t = 0; t < 1100 && idx < vector.size(); ++t) {
        auto next = vector.next1(idx);
        ASSERT_TRUE(next <= vector.size());

        if (next == vector.size()) {
            ASSERT_EQ(count, vector.num_set_bits());
            break;
        }

        EXPECT_EQ(count + 1, vector.rank1(next));

        count++;
        idx = next + 1;
    }
}

void test_next(const bit_vector &vector) {
    ASSERT_DEATH(vector.next1(vector.size()), "");
    ASSERT_DEATH(vector.next1(vector.size() + 1), "");
    ASSERT_DEATH(vector.next1(vector.size() * 2), "");

    test_next_subvector(vector, 0);
    test_next_subvector(vector, 0);
    if (vector.size() >= 64)
        test_next_subvector(vector, 63);
    if (vector.size() >= 128)
        test_next_subvector(vector, 127);
    test_next_subvector(vector, vector.size() / 5);
    test_next_subvector(vector, vector.size() / 2);
    test_next_subvector(vector, vector.size() * 2 / 3);
    test_next_subvector(vector, vector.size() - 1);
}

void test_prev_subvector(const bit_vector &vector, uint64_t idx) {
    if (vector.size() == 0)
        return;

    uint64_t count = vector.rank1(idx);

    for (size_t t = 0; t < 1100; ++t) {
        auto prev = vector.prev1(idx);
        ASSERT_TRUE(prev <= vector.size());

        if (prev == vector.size()) {
            ASSERT_TRUE(count <= 1);
            break;
        }

        EXPECT_EQ(count, vector.rank1(prev));

        if (!prev)
            break;

        count--;
        idx = prev - 1;
    }
}

void test_prev(const bit_vector &vector) {
    ASSERT_DEATH(vector.prev1(vector.size()), "");
    ASSERT_DEATH(vector.prev1(vector.size() + 1), "");
    ASSERT_DEATH(vector.prev1(vector.size() * 2), "");

    test_prev_subvector(vector, 0);
    if (vector.size() >= 64)
        test_prev_subvector(vector, 63);
    if (vector.size() >= 128)
        test_prev_subvector(vector, 127);
    test_prev_subvector(vector, vector.size() / 5);
    test_prev_subvector(vector, vector.size() / 2);
    test_prev_subvector(vector, vector.size() * 2 / 3);
    test_prev_subvector(vector, vector.size() - 1);
}

void reference_based_test(const bit_vector &vector,
                          const sdsl::bit_vector &reference) {
    EXPECT_EQ(reference.size(), vector.size());

    size_t max_rank = std::accumulate(reference.begin(), reference.end(), 0u);

    ASSERT_DEATH(vector.select1(0), "");

    for (size_t i : {1, 2, 10, 100, 1000}) {
        ASSERT_DEATH(vector.select1(max_rank + i), "");
        EXPECT_EQ(max_rank, vector.rank1(vector.size() + i - 2))
            << bit_vector_stat(reference);
    }
    ASSERT_DEATH(vector.select1(vector.size() + 1), "");
    ASSERT_DEATH(vector[vector.size()], "");
    ASSERT_DEATH(vector[vector.size() + 1], "");

    for (size_t i = 1; i <= max_rank; ++i) {
        EXPECT_TRUE(vector[vector.select1(i)]);
        EXPECT_EQ(i, vector.rank1(vector.select1(i)));
    }

    if (vector.size()) {
        EXPECT_EQ(vector[0], vector.rank1(0));
        EXPECT_EQ(vector[0], reference[0]);
    }

    for (size_t i = 1; i < vector.size(); ++i) {
        EXPECT_EQ(vector[i], vector.rank1(i) - vector.rank1(i - 1));
        EXPECT_EQ(vector[i], reference[i]);
    }

    EXPECT_EQ(reference, vector.to_vector());

    test_next(vector);
    test_prev(vector);
}


template <class T>
void test_bit_vector_queries() {
    std::unique_ptr<bit_vector> vector { new T() };
    ASSERT_TRUE(vector);
    reference_based_test(*vector, sdsl::bit_vector());

    vector.reset(new T(0, 1));
    ASSERT_TRUE(vector);
    reference_based_test(*vector, sdsl::bit_vector());

    vector.reset(new T(10, 0));
    ASSERT_TRUE(vector);
    reference_based_test(*vector, sdsl::bit_vector(10, 0));
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(0, (*vector)[i]);
        EXPECT_EQ(i + 1, vector->rank0(i));
        EXPECT_EQ(0u, vector->rank1(i));
        ASSERT_DEATH(vector->select1(i), "");
    }
    EXPECT_EQ(0u, vector->rank1(0));
    EXPECT_EQ(0u, vector->rank1(1'000));
    EXPECT_EQ(10u, vector->rank0(1'000));

    vector.reset(new T(10, 1));
    ASSERT_TRUE(vector);
    reference_based_test(*vector, sdsl::bit_vector(10, 1));
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
        EXPECT_EQ(0u, vector->rank0(i));
        EXPECT_TRUE((*vector)[vector->select1(i + 1)]);
        EXPECT_EQ(i, vector->select1(i + 1));
        EXPECT_EQ(i + 1, vector->rank1(i));
        EXPECT_EQ(i + 1, vector->rank1(vector->select1(i + 1)));
        EXPECT_EQ(i, vector->select1(vector->rank1(i)));
    }
    EXPECT_EQ(10u, vector->rank1(1'000));
    EXPECT_EQ(0u, vector->rank0(1'000));
    EXPECT_EQ(0u, vector->rank0(0));
    ASSERT_DEATH(vector->select1(1'000), "");
    ASSERT_DEATH(vector->select1(0), "");

    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    vector.reset(new T(numbers));
    ASSERT_TRUE(vector);
    reference_based_test(*vector, numbers);

    // set to 2 * max_step_size of bit_vector_stat
    sdsl::bit_vector large(2000, 0);
    large[0] = 1;
    large[1002] = 1;
    large[1999] = 1;
    vector.reset(new T(large));

    ASSERT_TRUE(vector);
    reference_based_test(*vector, large);
    EXPECT_EQ(2000u, vector->size());
}


TYPED_TEST(BitVectorTest, queries) {
    test_bit_vector_queries<TypeParam>();
}

TYPED_TEST(BitVectorTest, select1) {
    for (size_t size : { 1, 2, 3, 4, 5, 50, 51, 52 }) {
        for (size_t i = 0; i < size; ++i) {
            sdsl::bit_vector bv(size, 0);
            bv[i] = 1;
            TypeParam bit_vector(std::move(bv));
            EXPECT_EQ(i, bit_vector.select1(1));
        }
    }
}

TYPED_TEST(BitVectorTestSelect0, select0) {
    // Mainly test select0.
    auto vector = std::make_unique<TypeParam>(10, 1);
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
        ASSERT_DEATH(vector->select0(i), "");
    }

    vector.reset(new TypeParam(10, 0));
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(0, (*vector)[i]);
        EXPECT_EQ(i, vector->select0(i + 1));
        EXPECT_EQ(i + 1, vector->rank0(vector->select0(i + 1)));
        EXPECT_EQ(i, vector->select0(vector->rank0(i)));
    }

    vector.reset(new TypeParam({ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 }));
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    EXPECT_EQ(0u, vector->select0(1));
    EXPECT_EQ(1u, vector->select0(2));

    for (size_t size : { 1, 2, 3, 4, 5, 50, 51, 52 }) {
        for (size_t i = 0; i < size; ++i) {
            sdsl::bit_vector bv(size, 1);
            bv[i] = 0;
            TypeParam bit_vector(std::move(bv));
            EXPECT_EQ(i, bit_vector.select0(1));
        }
    }
}

void test_bit_vector_set(bit_vector *vector, sdsl::bit_vector *numbers) {
    reference_based_test(*vector, *numbers);

    for (size_t i = 0; i < numbers->size(); ++i) {
        bool value = (*numbers)[i];

        (*numbers)[i] = 1;
        vector->set(i, 1);
        reference_based_test(*vector, *numbers);

        (*numbers)[i] = 0;
        vector->set(i, 0);
        reference_based_test(*vector, *numbers);

        (*numbers)[i] = 1;
        vector->set(i, 1);
        reference_based_test(*vector, *numbers);

        (*numbers)[i] = 0;
        vector->set(i, 0);
        reference_based_test(*vector, *numbers);

        (*numbers)[i] = value;
        vector->set(i, value);
    }
}


TEST(bit_vector_dyn, set) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_dyn(numbers);
    ASSERT_TRUE(vector);

    test_bit_vector_set(vector, &numbers);

    delete vector;
}

TEST(bit_vector_stat, set) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_stat(numbers);
    ASSERT_TRUE(vector);

    test_bit_vector_set(vector, &numbers);

    delete vector;
}

TEST(bit_vector_sd, setException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_sd(numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}

TEST(bit_vector_rrr, setException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_rrr<>(numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}

TEST(bit_vector_small, setException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_small(numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}


void test_bit_vector_ins_del(bit_vector *vector,
                             const sdsl::bit_vector &reference) {
    reference_based_test(*vector, reference);

    std::vector<bool> numbers(reference.begin(), reference.end());

    for (size_t i = 0; i < numbers.size(); ++i) {
        numbers.insert(numbers.begin() + i, i % 16);
        vector->insert_bit(i, i % 16);
        reference_based_test(*vector, to_sdsl(numbers));
        numbers.erase(numbers.begin() + i, numbers.begin() + i + 1);
        vector->delete_bit(i);

        numbers.insert(numbers.begin() + i, 0);
        vector->insert_bit(i, 0);
        reference_based_test(*vector, to_sdsl(numbers));
        numbers.erase(numbers.begin() + i, numbers.begin() + i + 1);
        vector->delete_bit(i);
    }
}


TEST(bit_vector_dyn, InsertDelete) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_dyn(numbers);
    ASSERT_TRUE(vector);

    test_bit_vector_ins_del(vector, numbers);

    delete vector;
}

TEST(bit_vector_stat, InsertDelete) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_stat(numbers);
    ASSERT_TRUE(vector);

    test_bit_vector_ins_del(vector, numbers);

    delete vector;
}

TEST(bit_vector_sd, InsertDeleteException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_sd(numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_ins_del(vector, numbers), "");

    delete vector;
}

TEST(bit_vector_rrr, InsertDeleteException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_rrr<>(numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_ins_del(vector, numbers), "");

    delete vector;
}

TEST(bit_vector_small, InsertDeleteException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector *vector = new bit_vector_small(numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_ins_del(vector, numbers), "");

    delete vector;
}


TEST(bit_vector_dyn, Serialization) {
    std::vector<std::initializer_list<bool>> init_lists = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    for (auto init_list : init_lists) {
        sdsl::bit_vector numbers(init_list);
        std::unique_ptr<bit_vector> vector { new bit_vector_dyn(numbers) };
        ASSERT_TRUE(vector);
        std::ofstream outstream(test_dump_basename, std::ios::binary);
        vector->serialize(outstream);
        outstream.close();

        vector.reset(new bit_vector_dyn());
        ASSERT_TRUE(vector);
        std::ifstream instream(test_dump_basename, std::ios::binary);
        ASSERT_TRUE(vector->load(instream));

        reference_based_test(*vector, numbers);
    }
}

TYPED_TEST(BitVectorTest, Serialization) {
    std::vector<std::initializer_list<bool>> init_lists = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    for (auto init_list : init_lists) {
        sdsl::bit_vector numbers(init_list);
        std::unique_ptr<bit_vector> vector { new TypeParam(numbers) };
        ASSERT_TRUE(vector);
        std::ofstream outstream(test_dump_basename, std::ios::binary);
        vector->serialize(outstream);
        outstream.close();

        vector.reset(new TypeParam());
        ASSERT_TRUE(vector);
        std::ifstream instream(test_dump_basename, std::ios::binary);
        ASSERT_TRUE(vector->load(instream));

        reference_based_test(*vector, numbers);
    }
}

TEST(bit_vector_sd, SerializationCatchErrorWhenLoadingSdVector) {
    std::vector<std::initializer_list<bool>> init_lists = {
        { },
        { 1, },
        { 0, },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    for (auto init_list : init_lists) {
        sdsl::bit_vector bv(init_list);
        {
            std::ofstream outstream(test_dump_basename, std::ios::binary);
            sdsl::sd_vector<>(bv).serialize(outstream);
            outstream.close();
            bit_vector_sd vector;
            std::ifstream instream(test_dump_basename, std::ios::binary);
            ASSERT_FALSE(vector.load(instream));
        }
        {
            std::ofstream outstream(test_dump_basename, std::ios::binary);
            bit_vector_sd(bv).serialize(outstream);
            outstream.close();
            bit_vector_sd vector;
            std::ifstream instream(test_dump_basename, std::ios::binary);
            ASSERT_TRUE(vector.load(instream));
        }
    }
}


TYPED_TEST(BitVectorTest, MoveConstructor) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    TypeParam first(numbers);
    TypeParam second(std::move(first));
    reference_based_test(second, numbers);
}


TYPED_TEST(BitVectorTest, MoveAssignment) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    TypeParam first(numbers);
    TypeParam second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, MoveAssignmentSparse) {
    std::initializer_list<bool> init_list = { 0, 0, 0, 1, 0, 0, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_sd first(numbers);
    ASSERT_FALSE(first.is_inverted());
    bit_vector_sd second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, MoveAssignmentDense) {
    std::initializer_list<bool> init_list = { 1, 1, 0, 1, 0, 0, 1, 0,
                                              1, 1, 0, 1, 1, 1, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_sd first(numbers);
    ASSERT_TRUE(first.is_inverted());
    bit_vector_sd second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_stat, InitializeByBitsSparse) {
    std::vector<uint64_t> set_bits = { 3, 6, 9, 14, 15 };
    std::initializer_list<bool> init_list = { 0, 0, 0, 1, 0, 0, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_stat first(numbers);
    bit_vector_stat second(
        [&](const auto &callback) {
            for (uint64_t pos : set_bits) {
                callback(pos);
            }
        },
        first.size()
    );
    reference_based_test(second, numbers);
}

TEST(bit_vector_stat, InitializeByBitsDense) {
    std::vector<uint64_t> set_bits = { 0, 1, 3, 6, 8, 9, 11, 12, 13, 14, 15 };
    std::initializer_list<bool> init_list = { 1, 1, 0, 1, 0, 0, 1, 0,
                                              1, 1, 0, 1, 1, 1, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_stat first(numbers);
    bit_vector_stat second(
        [&](const auto &callback) {
            for (uint64_t pos : set_bits) {
                callback(pos);
            }
        },
        first.size()
    );
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, InitializeByBitsSparse) {
    std::vector<uint64_t> set_bits = { 3, 6, 9, 14, 15 };
    std::initializer_list<bool> init_list = { 0, 0, 0, 1, 0, 0, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_sd first(numbers);
    ASSERT_FALSE(first.is_inverted());
    bit_vector_sd second(
        [&](const auto &callback) {
            for (uint64_t pos : set_bits) {
                callback(pos);
            }
        },
        first.size(), set_bits.size());
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, InitializeByBitsDense) {
    std::vector<uint64_t> set_bits = { 0, 1, 3, 6, 8, 9, 11, 12, 13, 14, 15 };
    std::initializer_list<bool> init_list = { 1, 1, 0, 1, 0, 0, 1, 0,
                                              1, 1, 0, 1, 1, 1, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_sd first(numbers);
    ASSERT_TRUE(first.is_inverted());
    bit_vector_sd second(
        [&](const auto &callback) {
            for (uint64_t pos : set_bits) {
                callback(pos);
            }
        },
        first.size(), set_bits.size());
    reference_based_test(second, numbers);
}

TEST(bit_vector_stat, ConcurrentReadingAfterWriting) {
    ThreadPool thread_pool(3);
    bit_vector_stat vector;

    std::vector<bool> bits;
    std::vector<uint64_t> ranks = { 0 };

    for (size_t i = 0; i < 10'000'000; ++i) {
        bits.push_back((i + (i * i) % 31) % 2);
        if (bits.back()) {
            ranks.back()++;
        }
        ranks.push_back(ranks.back());
    }

    for (size_t i = 0; i < bits.size(); ++i) {
        ASSERT_EQ(0u, vector.rank1(i));
    }
    for (size_t i = 0; i < bits.size(); ++i) {
        vector.insert_bit(i, bits[i]);
    }
    for (size_t t = 0; t < 5; ++t) {
        thread_pool.enqueue([&]() {
            for (size_t i = 0; i < bits.size(); ++i) {
                ASSERT_EQ(ranks[i], vector.rank1(i));
            }
        });
    }
    thread_pool.join();
}

TEST(bit_vector_sd, CheckIfInverts) {
    sdsl::bit_vector vector(10);
    for (uint64_t i = 0; i < 1024; ++i) {
        vector.set_int(0, i);
        bit_vector_sd bvs(vector);
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(__builtin_popcountll(i) > 5, bvs.is_inverted());

        bit_vector_sd bvs_bits(
            [&](const auto &callback) {
                bvs.call_ones([&](uint64_t pos) {
                    callback(pos);
                });
            }, bvs.size(), bvs.num_set_bits());
        ASSERT_EQ(bvs.is_inverted(), bvs_bits.is_inverted());
    }
}

TYPED_TEST(BitVectorTest, to_vector_dense) {
    sdsl::bit_vector vector(10, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        TypeParam bvs(vector);
        if constexpr(std::is_same_v<TypeParam, bit_vector_sd>) {
            ASSERT_TRUE(bvs.is_inverted());
        }
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitVectorTest, to_vector_sparse) {
    sdsl::bit_vector vector(10);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = (i % 3) == 0;
        TypeParam bvs(vector);
        if constexpr(std::is_same_v<TypeParam, bit_vector_sd>) {
            ASSERT_FALSE(bvs.is_inverted());
        }
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitVectorTest, call_ones_dense) {
    sdsl::bit_vector vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        TypeParam bvs(vector);
        if constexpr(std::is_same_v<TypeParam, bit_vector_sd>) {
            ASSERT_TRUE(bvs.is_inverted());
        }
        sdsl::bit_vector copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitVectorTest, call_ones_sparse) {
    sdsl::bit_vector vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        TypeParam bvs(vector);
        if constexpr(std::is_same_v<TypeParam, bit_vector_sd>) {
            ASSERT_FALSE(bvs.is_inverted());
        }
        sdsl::bit_vector copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitVectorTest, copy) {
    auto bv_ptr = std::make_unique<TypeParam>(sdsl::bit_vector(10, 1));

    auto bv_copy = bv_ptr->copy();

    ASSERT_TRUE(dynamic_cast<TypeParam*>(bv_copy.get()));
    ASSERT_EQ(10u, bv_copy->rank1(9));
    ASSERT_EQ(0u, bv_copy->rank0(9));
}

template <class BitVectorFrom, class BitVectorTo>
void test_copy_convert_to() {
    std::unique_ptr<bit_vector> bv_ptr
        = std::make_unique<BitVectorFrom>(sdsl::bit_vector(10, 1));

    ASSERT_TRUE(dynamic_cast<BitVectorFrom*>(bv_ptr.get()));

    {
        auto bv_copy = bv_ptr->copy_to<BitVectorTo>();

        ASSERT_EQ(10u, bv_copy.rank1(9));
        ASSERT_EQ(0u, bv_copy.rank0(9));
    }

    {
        auto bv_copy = bv_ptr->convert_to<BitVectorTo>();

        ASSERT_EQ(10u, bv_copy.rank1(9));
        ASSERT_EQ(0u, bv_copy.rank0(9));
    }
}

TYPED_TEST(BitVectorTest, copy_to) {
    test_copy_convert_to< TypeParam, bit_vector_stat >();
    test_copy_convert_to< TypeParam, bit_vector_dyn >();
    test_copy_convert_to< TypeParam, bit_vector_sd >();
    test_copy_convert_to< TypeParam, bit_vector_rrr<> >();
    test_copy_convert_to< TypeParam, bit_vector_il<> >();
    test_copy_convert_to< TypeParam, bit_vector_hyb<> >();
    test_copy_convert_to< TypeParam, bit_vector_small >();
    test_copy_convert_to< TypeParam, bit_vector_smart >();
}


TEST(bit_vector, count_ones_empty) {
    sdsl::bit_vector v;
    EXPECT_EQ(0u, count_ones(v, 0, v.size()));
}

TEST(bit_vector, count_ones_all_zero) {
    {
        sdsl::bit_vector v(1, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(10, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999999, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(99999999, 0);
        EXPECT_EQ(0u, count_ones(v, 0, v.size()));
    }
}

TEST(bit_vector, count_ones_all_ones) {
    {
        sdsl::bit_vector v(1, 1);
        EXPECT_EQ(1u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(10, 1);
        EXPECT_EQ(10u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999, 1);
        EXPECT_EQ(999u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(999999, 1);
        EXPECT_EQ(999999u, count_ones(v, 0, v.size()));
    }
    {
        sdsl::bit_vector v(99999999, 1);
        EXPECT_EQ(99999999u, count_ones(v, 0, v.size()));
    }
}

TEST(bit_vector, count_ones) {
    {
        sdsl::bit_vector v(1);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(10);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(999);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(999999);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i + 100000 <= v.size(); i += 100000) {
            for (size_t j = i; j + 100000 <= v.size(); j += 100000) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    {
        sdsl::bit_vector v(99999999);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i + 10000000 <= v.size(); i += 10000000) {
            for (size_t j = i; j + 10000000 <= v.size(); j += 10000000) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    for (size_t size = 0; size < 500; ++size) {
        sdsl::bit_vector v(size);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i < v.size(); ++i) {
            for (size_t j = i; j < v.size(); ++j) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
    for (size_t size = 99999999; size < 9999999 + 200; ++size) {
        sdsl::bit_vector v(size);
        for (size_t i = 0; i < (v.size() + 63) >> 6; ++i) {
            v.data()[i] = 18446744073709551557llu + i * 32416189321llu;
        }
        sdsl::rank_support_v5<> rk(&v);
        for (size_t i = 0; i + 10000000 <= v.size(); i += 10000000) {
            for (size_t j = i; j + 10000000 <= v.size(); j += 10000000) {
                size_t first = count_ones(v, 0, i);
                EXPECT_EQ(rk(i), first);
                size_t second = count_ones(v, i, j);
                EXPECT_EQ(rk(j) - first, second);
                size_t third = count_ones(v, j, v.size());
                EXPECT_EQ(sdsl::util::cnt_one_bits(v), first + second + third);
            }
        }
    }
}

TEST(bit_vector, inner_prod_empty) {
    sdsl::bit_vector first;
    sdsl::bit_vector second;
    EXPECT_EQ(0u, inner_prod(first, second));
}

TEST(bit_vector, inner_prod_all_zero) {
    {
        sdsl::bit_vector first(1, 0);
        sdsl::bit_vector second(1, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10, 0);
        sdsl::bit_vector second(10, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999, 0);
        sdsl::bit_vector second(999, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999, 0);
        sdsl::bit_vector second(999999, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999, 0);
        sdsl::bit_vector second(99999999, 0);
        EXPECT_EQ(0u, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_all_ones) {
    {
        sdsl::bit_vector first(1, 1);
        sdsl::bit_vector second(1, 1);
        EXPECT_EQ(1u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10, 1);
        sdsl::bit_vector second(10, 1);
        EXPECT_EQ(10u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999, 1);
        sdsl::bit_vector second(999, 1);
        EXPECT_EQ(999u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999, 1);
        sdsl::bit_vector second(999999, 1);
        EXPECT_EQ(999999u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999, 1);
        sdsl::bit_vector second(99999999, 1);
        EXPECT_EQ(99999999u, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_all_ones_offset) {
    for (size_t size = 0; size < 200; ++size) {
        sdsl::bit_vector first(size, 1);
        sdsl::bit_vector second(size, 1);
        EXPECT_EQ(size, inner_prod(first, second));
    }
    for (size_t size = 99999999; size < 99999999 + 200; ++size) {
        sdsl::bit_vector first(size, 1);
        sdsl::bit_vector second(size, 1);
        EXPECT_EQ(size, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_same) {
    {
        sdsl::bit_vector first(1);
        sdsl::bit_vector second(1);
        uint64_t prod = 0;
        uint64_t bits = 1;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10);
        sdsl::bit_vector second(10);
        uint64_t prod = 0;
        uint64_t bits = 10;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999);
        sdsl::bit_vector second(999);
        uint64_t prod = 0;
        uint64_t bits = 999;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999);
        sdsl::bit_vector second(999999);
        uint64_t prod = 0;
        uint64_t bits = 999999;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999);
        sdsl::bit_vector second(99999999);
        uint64_t prod = 0;
        uint64_t bits = 99999999;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    for (size_t size = 0; size < 500; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        uint64_t prod = 0;
        uint64_t bits = size;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
    for (size_t size = 99999999; size < 9999999 + 200; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        uint64_t prod = 0;
        uint64_t bits = size;
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = second.data()[i] = val;
            for (size_t k = 0; k < 64 && bits--; ++k) {
                prod += val & 1;
                val >>= 1;
            }
        }
        EXPECT_EQ(prod, inner_prod(first, second));
    }
}

TEST(bit_vector, inner_prod_disjoint) {
    {
        sdsl::bit_vector first(1);
        sdsl::bit_vector second(1);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(10);
        sdsl::bit_vector second(10);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999);
        sdsl::bit_vector second(999);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(999999);
        sdsl::bit_vector second(999999);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    {
        sdsl::bit_vector first(99999999);
        sdsl::bit_vector second(99999999);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    for (size_t size = 0; size < 500; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
    for (size_t size = 99999999; size < 99999999 + 80; ++size) {
        sdsl::bit_vector first(size);
        sdsl::bit_vector second(size);
        for (size_t i = 0; i < (first.size() + 63) >> 6; ++i) {
            uint64_t val = 18446744073709551557llu + i * 32416189321llu;
            first.data()[i] = val;
            second.data()[i] = ~val;
        }
        EXPECT_EQ(0u, inner_prod(first, second));
    }
}

TYPED_TEST(BitVectorTest, operator_eq) {
    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 100000 }) {
        for (bool value : { false, true }) {
            const TypeParam bit_vector(size, value);

            EXPECT_EQ(bit_vector, bitmap_vector(size, value));
            EXPECT_EQ(bit_vector, bitmap_set(size, value));
            EXPECT_EQ(bit_vector, bitmap_adaptive(size, value));
            EXPECT_EQ(bit_vector, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
            EXPECT_EQ(bit_vector, bit_vector_stat(size, value));
            EXPECT_EQ(bit_vector, bit_vector_dyn(size, value));
            EXPECT_EQ(bit_vector, bit_vector_sd(size, value));
            EXPECT_EQ(bit_vector, bit_vector_rrr<>(size, value));
            EXPECT_EQ(bit_vector, bit_vector_il<>(size, value));
            EXPECT_EQ(bit_vector, bit_vector_hyb<>(size, value));
            EXPECT_EQ(bit_vector, bit_vector_small(size, value));
            EXPECT_EQ(bit_vector, bit_vector_smart(size, value));
        }
    }
}

TYPED_TEST(BitVectorTest, operator_neq) {
    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 100000 }) {
        bool value = false;

        const TypeParam bit_vector(size + 1, false);

        EXPECT_NE(bit_vector, bitmap_vector(size, value));
        EXPECT_NE(bit_vector, bitmap_set(size, value));
        EXPECT_NE(bit_vector, bitmap_adaptive(size, value));
        EXPECT_NE(bit_vector, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
        EXPECT_NE(bit_vector, bit_vector_stat(size, value));
        EXPECT_NE(bit_vector, bit_vector_dyn(size, value));
        EXPECT_NE(bit_vector, bit_vector_sd(size, value));
        EXPECT_NE(bit_vector, bit_vector_rrr<>(size, value));
        EXPECT_NE(bit_vector, bit_vector_il<>(size, value));
        EXPECT_NE(bit_vector, bit_vector_hyb<>(size, value));
        EXPECT_NE(bit_vector, bit_vector_small(size, value));
        EXPECT_NE(bit_vector, bit_vector_smart(size, value));
    }

    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 100000 }) {
        bool value = false;

        const TypeParam bit_vector(size + 1, true);

        EXPECT_NE(bit_vector, bitmap_vector(size, value));
        EXPECT_NE(bit_vector, bitmap_set(size, value));
        EXPECT_NE(bit_vector, bitmap_adaptive(size, value));
        EXPECT_NE(bit_vector, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
        EXPECT_NE(bit_vector, bit_vector_stat(size, value));
        EXPECT_NE(bit_vector, bit_vector_dyn(size, value));
        EXPECT_NE(bit_vector, bit_vector_sd(size, value));
        EXPECT_NE(bit_vector, bit_vector_rrr<>(size, value));
        EXPECT_NE(bit_vector, bit_vector_il<>(size, value));
        EXPECT_NE(bit_vector, bit_vector_hyb<>(size, value));
        EXPECT_NE(bit_vector, bit_vector_small(size, value));
        EXPECT_NE(bit_vector, bit_vector_smart(size, value));
    }

    for (uint64_t size : { 10, 64, 120, 128, 1000, 10000, 100000 }) {
        for (bool value : { false, true }) {
            const TypeParam bit_vector(size, !value);

            EXPECT_NE(bit_vector, bitmap_vector(size, value));
            EXPECT_NE(bit_vector, bitmap_set(size, value));
            EXPECT_NE(bit_vector, bitmap_adaptive(size, value));
            EXPECT_NE(bit_vector, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
            EXPECT_NE(bit_vector, bit_vector_stat(size, value));
            EXPECT_NE(bit_vector, bit_vector_dyn(size, value));
            EXPECT_NE(bit_vector, bit_vector_sd(size, value));
            EXPECT_NE(bit_vector, bit_vector_rrr<>(size, value));
            EXPECT_NE(bit_vector, bit_vector_il<>(size, value));
            EXPECT_NE(bit_vector, bit_vector_hyb<>(size, value));
            EXPECT_NE(bit_vector, bit_vector_small(size, value));
            EXPECT_NE(bit_vector, bit_vector_smart(size, value));
        }
    }
}

TEST(bit_vector, to_sdsl_empty) {
    EXPECT_EQ(0u, to_sdsl(std::vector<bool>()).size());
    EXPECT_EQ(0u, to_sdsl(std::vector<uint8_t>()).size());
}

TEST(bit_vector, to_sdsl_zeros) {
    for (size_t i = 0; i < 1025; ++i) {
        auto vector_bool = to_sdsl(std::vector<bool>(i, false));
        EXPECT_EQ(i, vector_bool.size());
        EXPECT_EQ(0u, sdsl::util::cnt_one_bits(vector_bool));

        auto vector_uint8_t = to_sdsl(std::vector<uint8_t>(i, false));
        EXPECT_EQ(i, vector_uint8_t.size());
        EXPECT_EQ(0u, sdsl::util::cnt_one_bits(vector_uint8_t));
    }
}

TEST(bit_vector, to_sdsl_ones) {
    for (size_t i = 0; i < 1025; ++i) {
        auto vector_bool = to_sdsl(std::vector<bool>(i, true));
        EXPECT_EQ(i, vector_bool.size());
        EXPECT_EQ(i, sdsl::util::cnt_one_bits(vector_bool));

        auto vector_uint8_t = to_sdsl(std::vector<uint8_t>(i, true));
        EXPECT_EQ(i, vector_uint8_t.size());
        EXPECT_EQ(i, sdsl::util::cnt_one_bits(vector_uint8_t));
    }
}

TEST(bit_vector, to_sdsl) {
    std::vector<bool> bits_bool;
    std::vector<uint8_t> bits_uint8_t_a;
    std::vector<uint8_t> bits_uint8_t_b;
    std::vector<uint64_t> ranks = { 0 };

    for (size_t i = 0; i < 10'000'000; ++i) {
        bits_bool.push_back((i + (i * i) % 31) % 2);
        bits_uint8_t_a.push_back(static_cast<bool>((i + (i * i) % 31) % 2));
        bits_uint8_t_b.push_back(((i + (i * i) % 31) % 2) ? 0xFF : 0);
        if (bits_uint8_t_a.back()) {
            ranks.back()++;
        }
        ranks.push_back(ranks.back());
    }

    bit_vector_stat a(to_sdsl(bits_bool));
    bit_vector_stat b(to_sdsl(bits_uint8_t_a));
    bit_vector_stat c(to_sdsl(bits_uint8_t_b));

    for (size_t i = 0; i < bits_bool.size(); ++i) {
        EXPECT_EQ(ranks[i], a.rank1(i));
        EXPECT_EQ(ranks[i], b.rank1(i));
        EXPECT_EQ(ranks[i], c.rank1(i));
    }
}

TEST(bit_vector, autocorrelate) {
    std::vector<std::tuple<sdsl::bit_vector, sdsl::bit_vector, uint8_t>> results {
        { sdsl::bit_vector({}), sdsl::bit_vector({}), 5 },
        { sdsl::bit_vector({ 0, 0, 0, 0, 0 }), sdsl::bit_vector({ 0, 0, 0, 0, 0 }), 5 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 5 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 4 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 3 },
        { sdsl::bit_vector({ 1, 1, 1, 1, 1 }), sdsl::bit_vector({ 1, 1, 1, 1, 1 }), 2 },
        { sdsl::bit_vector({ 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0  }),
          sdsl::bit_vector({ 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  }), 3},
    };

    for (const auto &[input, output, offset] : results) {
        sdsl::bit_vector reference(input);

        if (input.size() >= offset) {
            for (size_t i = 0; i < reference.size(); ++i) {
                for (uint8_t j = 1; j < offset && i + j < reference.size(); ++j) {
                    reference[i] = reference[i] && input[i + j];
                }
            }
        }

        ASSERT_EQ(reference, output);
        EXPECT_EQ(reference, autocorrelate(input, offset));
    };
}
