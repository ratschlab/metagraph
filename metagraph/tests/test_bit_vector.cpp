#include "gtest/gtest.h"

#include "bit_vector.hpp"
#include "utils.hpp"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(a, b) (void)0
#endif
#endif


const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bit_vector_dump_test";


void reference_based_test(const bit_vector &vector,
                          const std::vector<bool> &reference) {
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
        EXPECT_EQ(i, vector.rank1(vector.select1(i)));
    }

    EXPECT_EQ(vector[0], vector.rank1(0));
    EXPECT_EQ(vector[0], reference[0]);

    for (size_t i = 1; i < vector.size(); ++i) {
        EXPECT_EQ(vector[i], vector.rank1(i) - vector.rank1(i - 1));
        EXPECT_EQ(vector[i], reference[i]);
    }

    EXPECT_EQ(reference, vector.to_vector());
}


template <class T>
void test_bit_vector_queries() {
    std::unique_ptr<bit_vector> vector { new T() };
    ASSERT_TRUE(vector);

    vector.reset(new T(0, 1));
    ASSERT_TRUE(vector);

    vector.reset(new T(10, 0));
    ASSERT_TRUE(vector);
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
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
        EXPECT_EQ(0u, vector->rank0(i));
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
    vector.reset(new T(init_list));
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);
    reference_based_test(*vector, numbers);
}


TEST(bit_vector_dyn, queries) {
    test_bit_vector_queries<bit_vector_dyn>();
}

TEST(bit_vector_stat, queries) {
    test_bit_vector_queries<bit_vector_stat>();
}

TEST(bit_vector_sd, queries) {
    test_bit_vector_queries<bit_vector_sd>();
}

TEST(bit_vector_rrr, queries) {
    test_bit_vector_queries<bit_vector_rrr<>>();
}

TEST(bit_vector_small, queries) {
    test_bit_vector_queries<bit_vector_small>();
}

TEST(bit_vector_rrr, nonCommonQueries) {
    // Mainly test select0.
    auto vector = std::make_unique<bit_vector_rrr<>>(10, 1);
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
        ASSERT_DEATH(vector->select0(i), "");
    }

    vector.reset(new bit_vector_rrr<>(10, 0));
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(0, (*vector)[i]);
        EXPECT_EQ(i, vector->select0(i + 1));
        EXPECT_EQ(i + 1, vector->rank0(vector->select0(i + 1)));
        EXPECT_EQ(i, vector->select0(vector->rank0(i)));
    }
}

void test_bit_vector_set(bit_vector *vector, std::vector<bool> *numbers) {
    reference_based_test(*vector, *numbers);

    for (size_t i = 0; i < numbers->size(); ++i) {
        bool value = numbers->at(i);

        numbers->at(i) = 1;
        vector->set(i, 1);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = 0;
        vector->set(i, 0);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = 1;
        vector->setBitQuick(i, 1);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = 0;
        vector->setBitQuick(i, 0);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = value;
        vector->set(i, value);
    }
}


TEST(bit_vector_dyn, set) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_dyn(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    test_bit_vector_set(vector, &numbers);

    delete vector;
}

TEST(bit_vector_stat, set) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_stat(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    test_bit_vector_set(vector, &numbers);

    delete vector;
}

TEST(bit_vector_sd, setException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_sd(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}

TEST(bit_vector_rrr, setException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_rrr<>(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}

TEST(bit_vector_small, setException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_small(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}


void test_bit_vector_ins_del(bit_vector *vector, std::vector<bool> *numbers) {
    reference_based_test(*vector, *numbers);

    for (size_t i = 0; i < numbers->size(); ++i) {
        numbers->insert(numbers->begin() + i, i % 16);
        vector->insertBit(i, i % 16);
        reference_based_test(*vector, *numbers);
        numbers->erase(numbers->begin() + i, numbers->begin() + i + 1);
        vector->deleteBit(i);

        numbers->insert(numbers->begin() + i, 0);
        vector->insertBit(i, 0);
        reference_based_test(*vector, *numbers);
        numbers->erase(numbers->begin() + i, numbers->begin() + i + 1);
        vector->deleteBit(i);
    }
}


TEST(bit_vector_dyn, InsertDelete) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_dyn(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    test_bit_vector_ins_del(vector, &numbers);

    delete vector;
}

TEST(bit_vector_stat, InsertDelete) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_stat(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    test_bit_vector_ins_del(vector, &numbers);

    delete vector;
}

TEST(bit_vector_sd, InsertDeleteException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_sd(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_ins_del(vector, &numbers), "");

    delete vector;
}

TEST(bit_vector_rrr, InsertDeleteException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_rrr<>(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_ins_del(vector, &numbers), "");

    delete vector;
}

TEST(bit_vector_small, InsertDeleteException) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    bit_vector *vector = new bit_vector_small(init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_ins_del(vector, &numbers), "");

    delete vector;
}


TEST(bit_vector_dyn, Serialization) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::unique_ptr<bit_vector> vector { new bit_vector_dyn(init_list) };
    std::vector<bool> numbers(init_list);
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

TEST(bit_vector_stat, Serialization) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::unique_ptr<bit_vector> vector { new bit_vector_stat(init_list) };
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector);
    std::ofstream outstream(test_dump_basename, std::ios::binary);
    vector->serialize(outstream);
    outstream.close();

    vector.reset(new bit_vector_stat());
    ASSERT_TRUE(vector);
    std::ifstream instream(test_dump_basename, std::ios::binary);
    ASSERT_TRUE(vector->load(instream));

    reference_based_test(*vector, numbers);
}

TEST(bit_vector_sd, Serialization) {
    std::vector<std::initializer_list<bool>> init_lists = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    for (auto init_list : init_lists) {
        std::vector<bool> numbers(init_list);
        bit_vector *vector = new bit_vector_sd(numbers);
        ASSERT_EQ(16u, numbers.size());
        ASSERT_TRUE(vector);
        std::ofstream outstream(test_dump_basename, std::ios::binary);
        vector->serialize(outstream);
        outstream.close();
        delete vector;

        vector = new bit_vector_sd;
        ASSERT_TRUE(vector);
        std::ifstream instream(test_dump_basename, std::ios::binary);
        ASSERT_TRUE(vector->load(instream));

        reference_based_test(*vector, numbers);

        delete vector;
    }
}

TEST(bit_vector_rrr, Serialization) {
    std::vector<std::initializer_list<bool>> init_lists = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    for (auto init_list : init_lists) {
        sdsl::bit_vector numbers(init_list);
        bit_vector *vector = new bit_vector_rrr<>(numbers);
        ASSERT_EQ(16u, numbers.size());
        ASSERT_TRUE(vector);
        std::ofstream outstream(test_dump_basename, std::ios::binary);
        vector->serialize(outstream);
        outstream.close();
        delete vector;

        vector = new bit_vector_rrr<>;
        ASSERT_TRUE(vector);
        std::ifstream instream(test_dump_basename, std::ios::binary);
        ASSERT_TRUE(vector->load(instream));

        reference_based_test(*vector, std::vector<bool>(numbers.begin(), numbers.end()));

        delete vector;
    }
}

TEST(bit_vector_small, Serialization) {
    std::vector<std::initializer_list<bool>> init_lists = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
        { 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };
    for (auto init_list : init_lists) {
        std::vector<bool> numbers(init_list);
        bit_vector *vector = new bit_vector_small(numbers);
        ASSERT_EQ(16u, numbers.size());
        ASSERT_TRUE(vector);
        std::ofstream outstream(test_dump_basename, std::ios::binary);
        vector->serialize(outstream);
        outstream.close();
        delete vector;

        vector = new bit_vector_small;
        ASSERT_TRUE(vector);
        std::ifstream instream(test_dump_basename, std::ios::binary);
        ASSERT_TRUE(vector->load(instream));

        reference_based_test(*vector, numbers);

        delete vector;
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


TEST(bit_vector_stat, MoveConstructor) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_stat vector(numbers);
    reference_based_test(vector, numbers);
}

TEST(bit_vector_dyn, MoveConstructor) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_dyn vector(numbers);
    reference_based_test(vector, numbers);
}

TEST(bit_vector_sd, MoveConstructor) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_sd vector(numbers);
    reference_based_test(vector, numbers);
}

TEST(bit_vector_rrr, MoveConstructor) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_rrr<> vector(numbers);
    reference_based_test(vector, std::vector<bool>(numbers.begin(), numbers.end()));
}

TEST(bit_vector_small, MoveConstructor) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_small vector(numbers);
    reference_based_test(vector, numbers);
}


TEST(bit_vector_stat, MoveAssignment) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_stat first(numbers);
    bit_vector_stat second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_dyn, MoveAssignment) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_dyn first(numbers);
    bit_vector_dyn second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, MoveAssignment) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_sd first(numbers);
    bit_vector_sd second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_rrr, MoveAssignment) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    sdsl::bit_vector numbers(init_list);
    bit_vector_rrr<> first(numbers);
    bit_vector_rrr<> second;
    second = std::move(first);
    reference_based_test(second, std::vector<bool>(numbers.begin(), numbers.end()));
}

TEST(bit_vector_small, MoveAssignment) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_small first(numbers);
    bit_vector_small second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, MoveAssignmentSparse) {
    std::initializer_list<bool> init_list = { 0, 0, 0, 1, 0, 0, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_sd first(numbers);
    ASSERT_FALSE(first.is_inverted());
    bit_vector_sd second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, MoveAssignmentDense) {
    std::initializer_list<bool> init_list = { 1, 1, 0, 1, 0, 0, 1, 0,
                                              1, 1, 0, 1, 1, 1, 1, 1 };
    std::vector<bool> numbers(init_list);
    bit_vector_sd first(numbers);
    ASSERT_TRUE(first.is_inverted());
    bit_vector_sd second;
    second = std::move(first);
    reference_based_test(second, numbers);
}

TEST(bit_vector_sd, InitializeByBitsSparse) {
    std::vector<uint64_t> set_bits = { 3, 6, 9, 14, 15 };
    std::initializer_list<bool> init_list = { 0, 0, 0, 1, 0, 0, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);
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
    std::vector<bool> numbers(init_list);
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
    utils::ThreadPool thread_pool(3);
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
        vector.insertBit(i, bits[i]);
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

TEST(bit_vector_sd, to_vector_dense) {
    std::vector<bool> vector(10, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_sd bvs(vector);
        ASSERT_TRUE(bvs.is_inverted());
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_sd, to_vector_sparse) {
    std::vector<bool> vector(10);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = (i % 3) == 0;
        bit_vector_sd bvs(vector);
        ASSERT_FALSE(bvs.is_inverted());
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_small, to_vector_dense) {
    std::vector<bool> vector(10, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_small bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_small, to_vector_sparse) {
    std::vector<bool> vector(10);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = (i % 3) == 0;
        bit_vector_small bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}


TEST(bit_vector_rrr, to_vector_dense) {
    sdsl::bit_vector vector(10, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_rrr<> bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, std::vector<bool>(vector.begin(), vector.end()));
    }
}

TEST(bit_vector_rrr, to_vector_sparse) {
    sdsl::bit_vector vector(10);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = (i % 3) == 0;
        bit_vector_rrr<> bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, std::vector<bool>(vector.begin(), vector.end()));
    }
}

TEST(bit_vector_stat, to_vector_dense) {
    std::vector<bool> vector(10, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_stat bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_stat, to_vector_sparse) {
    std::vector<bool> vector(10);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = (i % 3) == 0;
        bit_vector_stat bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_dyn, to_vector_dense) {
    std::vector<bool> vector(10, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_dyn bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_dyn, to_vector_sparse) {
    std::vector<bool> vector(10);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = (i % 3) == 0;
        bit_vector_dyn bvs(vector);
        auto copy = bvs.to_vector();
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}


TEST(bit_vector_small, call_ones_dense) {
    std::vector<bool> vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_small bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_small, call_ones_sparse) {
    std::vector<bool> vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bit_vector_small bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_sd, call_ones_dense) {
    std::vector<bool> vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_sd bvs(vector);
        ASSERT_TRUE(bvs.is_inverted());
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_sd, call_ones_sparse) {
    std::vector<bool> vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bit_vector_sd bvs(vector);
        ASSERT_FALSE(bvs.is_inverted());
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_rrr, call_ones_dense) {
    sdsl::bit_vector vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_rrr<> bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, std::vector<bool>(vector.begin(), vector.end()));
    }
}

TEST(bit_vector_rrr, call_ones_sparse) {
    sdsl::bit_vector vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bit_vector_rrr<> bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, std::vector<bool>(vector.begin(), vector.end()));
    }
}

TEST(bit_vector_dyn, call_ones_dense) {
    std::vector<bool> vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_dyn bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_dyn, call_ones_sparse) {
    std::vector<bool> vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bit_vector_dyn bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_stat, call_ones_dense) {
    std::vector<bool> vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bit_vector_stat bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}

TEST(bit_vector_stat, call_ones_sparse) {
    std::vector<bool> vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bit_vector_stat bvs(vector);
        std::vector<bool> copy(vector.size(), false);
        bvs.call_ones([&copy](auto i) { copy[i] = true; });
        // check if the it inverts when more than half of the bits are set
        ASSERT_EQ(copy, vector);
    }
}
