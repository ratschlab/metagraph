#include "gtest/gtest.h"

#include "wavelet_tree.hpp"


const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bit_vector_dump_test";


void reference_based_test(const wavelet_tree &vector,
                          const std::vector<uint64_t> &reference) {
    for (uint64_t c = 0; c < 16; ++c) {
        uint64_t max_rank = std::count(reference.begin(), reference.end(), c);

        ASSERT_DEATH(vector.select(c, 0), "");
        ASSERT_DEATH(vector.select(c, max_rank + 1), "");

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
}


template <class T>
void test_wavelet_tree_queries() {
    wavelet_tree *vector = new T(4);
    ASSERT_TRUE(vector);
    delete vector;

    vector = new T(4, std::vector<int>(10, 0));
    ASSERT_TRUE(vector);

    EXPECT_EQ(10u, vector->size());
    ASSERT_DEATH((*vector)[vector->size()], "");
    ASSERT_DEATH((*vector)[vector->size() + 1], "");

    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(0u, (*vector)[i]);
        EXPECT_EQ(i + 1, vector->rank(0, i));
        EXPECT_EQ(0u, vector->rank(1, i));
        EXPECT_EQ(0u, vector->rank(2, i));
        EXPECT_EQ(i, vector->select(0, i + 1));
        ASSERT_DEATH(vector->select(1, i + 1), "");
    }
    EXPECT_EQ(0u, vector->rank(2, 0));
    EXPECT_EQ(0u, vector->rank(2, 1'000));
    EXPECT_EQ(10u, vector->rank(0, 1'000));
    ASSERT_DEATH(vector->select(1, 1'000), "");
    ASSERT_DEATH(vector->select(1, 0), "");
    delete vector;

    vector = new T(4, std::vector<int>(10, 2));
    ASSERT_TRUE(vector);

    EXPECT_EQ(10u, vector->size());
    ASSERT_DEATH((*vector)[vector->size()], "");
    ASSERT_DEATH((*vector)[vector->size() + 1], "");

    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(2u, (*vector)[i]);
        EXPECT_EQ(i + 1, vector->rank(2, i));
        EXPECT_EQ(0u, vector->rank(1, i));
        EXPECT_EQ(0u, vector->rank(3, i));
        EXPECT_EQ(i, vector->select(2, i + 1));
        ASSERT_DEATH(vector->select(1, i + 1), "");
    }
    EXPECT_EQ(0u, vector->rank(0, 0));
    EXPECT_EQ(0u, vector->rank(0, 1'000));
    EXPECT_EQ(10u, vector->rank(2, 1'000));
    ASSERT_DEATH(vector->select(3, 1'000), "");
    ASSERT_DEATH(vector->select(1, 0), "");
    delete vector;

    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    vector = new T(4, numbers);
    ASSERT_TRUE(vector);
    reference_based_test(*vector, numbers);

    delete vector;
}


TEST(wavelet_tree_stat, Queries) {
    test_wavelet_tree_queries<wavelet_tree_stat>();
}

TEST(wavelet_tree_small, Queries) {
    test_wavelet_tree_queries<wavelet_tree_small>();
}

TEST(wavelet_tree_dyn, Queries) {
    test_wavelet_tree_queries<wavelet_tree_dyn>();
}


void test_bit_vector_set(wavelet_tree *vector, std::vector<uint64_t> *numbers) {
    reference_based_test(*vector, *numbers);

    for (size_t i = 0; i < numbers->size(); ++i) {
        uint64_t value = numbers->at(i);

        numbers->at(i) = 1;
        vector->set(i, 1);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = 0;
        vector->set(i, 0);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = value;
        vector->set(i, value);
    }
}


TEST(wavelet_tree_stat, Set) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_stat(4, numbers);
    ASSERT_TRUE(vector);

    test_bit_vector_set(vector, &numbers);

    delete vector;
}


TEST(wavelet_tree_small, Set) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_small(4, numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_bit_vector_set(vector, &numbers), "");

    delete vector;
}


TEST(wavelet_tree_dyn, Set) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_dyn(4, numbers);
    ASSERT_TRUE(vector);

    test_bit_vector_set(vector, &numbers);

    delete vector;
}


void test_wavelet_tree_ins_del(wavelet_tree *vector, std::vector<uint64_t> *numbers) {
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


TEST(wavelet_tree_stat, InsertDelete) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_stat(4, numbers);
    ASSERT_TRUE(vector);

    test_wavelet_tree_ins_del(vector, &numbers);

    delete vector;
}


TEST(wavelet_tree_small, InsertDelete) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_small(4, numbers);
    ASSERT_TRUE(vector);

    ASSERT_DEATH(test_wavelet_tree_ins_del(vector, &numbers), "");

    delete vector;
}


TEST(wavelet_tree_dyn, InsertDelete) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_dyn(4, numbers);
    ASSERT_TRUE(vector);

    test_wavelet_tree_ins_del(vector, &numbers);

    delete vector;
}


std::vector<uint64_t> to_std_vector(const sdsl::int_vector<> &int_vector) {
    return std::vector<uint64_t>(int_vector.begin(), int_vector.end());
}

TEST(wavelet_tree_stat, ToVector) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_stat(2, numbers);
    ASSERT_TRUE(vector);
    EXPECT_EQ(numbers, to_std_vector(vector->to_vector()));
    delete vector;
}

TEST(wavelet_tree_small, ToVector) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_small(2, numbers);
    ASSERT_TRUE(vector);
    EXPECT_EQ(numbers, to_std_vector(vector->to_vector()));
    delete vector;
}

TEST(wavelet_tree_dyn, ToVector) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_dyn(2, numbers);
    ASSERT_TRUE(vector);
    EXPECT_EQ(numbers, to_std_vector(vector->to_vector()));
    delete vector;
}


TEST(wavelet_tree_stat, Serialization) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_stat(4, numbers);
    ASSERT_TRUE(vector);
    std::ofstream outstream(test_dump_basename);
    vector->serialise(outstream);
    outstream.close();
    delete vector;

    vector = new wavelet_tree_stat(4);
    ASSERT_TRUE(vector);
    std::ifstream instream(test_dump_basename);
    ASSERT_TRUE(vector->deserialise(instream));

    reference_based_test(*vector, numbers);

    delete vector;
}


TEST(wavelet_tree_small, Serialization) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_small(4, numbers);
    ASSERT_TRUE(vector);
    std::ofstream outstream(test_dump_basename);
    vector->serialise(outstream);
    outstream.close();
    delete vector;

    vector = new wavelet_tree_small(4);
    ASSERT_TRUE(vector);
    std::ifstream instream(test_dump_basename);
    ASSERT_TRUE(vector->deserialise(instream));

    reference_based_test(*vector, numbers);

    delete vector;
}


TEST(wavelet_tree_dyn, Serialization) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree *vector = new wavelet_tree_dyn(4, numbers);
    ASSERT_TRUE(vector);
    std::ofstream outstream(test_dump_basename);
    vector->serialise(outstream);
    outstream.close();
    delete vector;

    vector = new wavelet_tree_dyn(4);
    ASSERT_TRUE(vector);
    std::ifstream instream(test_dump_basename);
    ASSERT_TRUE(vector->deserialise(instream));

    reference_based_test(*vector, numbers);

    delete vector;
}

TEST(wavelet_tree_stat, MoveConstructor) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size());
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    wavelet_tree_stat vector(2, std::move(int_vector));

    reference_based_test(vector, numbers);
}

TEST(wavelet_tree_small, MoveConstructor) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size());
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    wavelet_tree_small vector(2, std::move(int_vector));

    reference_based_test(vector, numbers);
}

TEST(wavelet_tree_stat, MoveAssignment) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size());
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    wavelet_tree_stat vector(1);
    vector = wavelet_tree_stat(2, std::move(int_vector));

    reference_based_test(vector, numbers);
}

TEST(wavelet_tree_small, MoveAssignment) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };

    sdsl::int_vector<> int_vector(numbers.size());
    for (size_t i = 0; i < numbers.size(); ++i) {
        int_vector[i] = numbers[i];
    }

    wavelet_tree_small vector(1);
    vector = wavelet_tree_small(2, std::move(int_vector));

    reference_based_test(vector, numbers);
}

TEST(wavelet_tree_stat, BeyondTheDNA) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree_stat vector(7, numbers);
    EXPECT_EQ(numbers, to_std_vector(vector.to_vector()));
}

TEST(wavelet_tree_small, BeyondTheDNA) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree_small vector(7, numbers);
    EXPECT_EQ(numbers, to_std_vector(vector.to_vector()));
}

TEST(wavelet_tree_dyn, BeyondTheDNA) {
    std::vector<uint64_t> numbers = { 0, 1, 0, 1, 1, 1, 1, 0,
                                      0, 1, 2, 0, 3, 2, 1, 1 };
    wavelet_tree_dyn vector(7, numbers);
    EXPECT_EQ(numbers, to_std_vector(vector.to_vector()));
}
