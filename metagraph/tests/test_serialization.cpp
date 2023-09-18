#include "gtest/gtest.h"

#include <fstream>

#include "common/serialization.hpp"


namespace {

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/vector_dump_test";


TEST(Serialization, SerializationVectorBool) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::vector<bool> numbers(init_list);

    std::ofstream outstream(test_dump_basename, std::ios::binary);
    serialize_number_vector(outstream, numbers);
    outstream.close();

    std::ifstream instream(test_dump_basename, std::ios::binary);
    std::vector<bool> loaded_vector;
    ASSERT_TRUE(load_number_vector(instream, &loaded_vector));

    ASSERT_EQ(numbers, loaded_vector);
}

TEST(Serialization, SerializationVectorUInt8) {
    std::initializer_list<uint8_t> init_list = { 0, 1, 255, 1, 1, 1, 1, 0,
                                                 0, 1, 100, 0, 0, 4, 1, 1 };
    std::vector<uint8_t> numbers(init_list);

    std::ofstream outstream(test_dump_basename, std::ios::binary);
    serialize_number_vector(outstream, numbers);
    outstream.close();

    std::ifstream instream(test_dump_basename, std::ios::binary);
    std::vector<uint8_t> loaded_vector;
    ASSERT_TRUE(load_number_vector(instream, &loaded_vector));

    ASSERT_EQ(numbers, loaded_vector);
}

TEST(Serialization, SerializationVectorUInt64) {
    std::initializer_list<uint64_t> init_list = { 0, 1llu << 63, 255, 1, 1, 1, 1, 0,
                                                  0, 1, 100, 1llu << 63, 0, 4, 1, 1 };
    std::vector<uint64_t> numbers(init_list);

    std::ofstream outstream(test_dump_basename, std::ios::binary);
    serialize_number_vector(outstream, numbers);
    outstream.close();

    std::ifstream instream(test_dump_basename, std::ios::binary);
    std::vector<uint64_t> loaded_vector;
    ASSERT_TRUE(load_number_vector(instream, &loaded_vector));

    ASSERT_EQ(numbers, loaded_vector);
}

TEST(Serialization, SerializationVectorUInt64Effective1Bit) {
    std::initializer_list<uint64_t> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                                  0, 1, 0, 0, 0, 1, 1, 1 };
    std::vector<uint64_t> numbers(init_list);

    std::ofstream outstream(test_dump_basename, std::ios::binary);
    serialize_number_vector(outstream, numbers, 1);
    outstream.close();

    std::ifstream instream(test_dump_basename, std::ios::binary);
    std::vector<uint64_t> loaded_vector;
    ASSERT_TRUE(load_number_vector(instream, &loaded_vector));

    ASSERT_EQ(numbers, loaded_vector);
}


void test_random_vector(size_t length) {
    for (size_t i = 0; i < 10; ++i) {
        std::vector<uint8_t> numbers(length);
        for (size_t j = 0; j < length; ++j) {
            numbers[j] = rand() % 256;
        }

        std::ofstream outstream(test_dump_basename, std::ios::binary);
        serialize_number_vector(outstream, numbers);
        outstream.close();

        {
            std::ifstream instream(test_dump_basename, std::ios::binary);
            ASSERT_EQ(numbers.size(), get_number_vector_size(instream));
            std::vector<uint8_t> loaded_vector;
            ASSERT_TRUE(load_number_vector(instream, &loaded_vector));
            ASSERT_EQ(numbers, loaded_vector);
        }
        {
            std::ifstream instream(test_dump_basename, std::ios::binary);
            std::vector<uint8_t> loaded_vector;
            ASSERT_TRUE(load_number_vector(instream, &loaded_vector));
            ASSERT_EQ(numbers, loaded_vector);
        }
    }
}

TEST(Serialization, SerializationRandomUInt8Vector0) {
    test_random_vector(0);
}

TEST(Serialization, SerializationRandomUInt8Vector10) {
    test_random_vector(10);
}

TEST(Serialization, SerializationRandomUInt8Vector1000) {
    test_random_vector(1000);
}

TEST(Serialization, SerializationRandomUInt8Vector1000000) {
    test_random_vector(1000000);
}

TEST(Serialization, SerializationRandomUInt8Vector10000000) {
    test_random_vector(10000000);
}

} // namespace
