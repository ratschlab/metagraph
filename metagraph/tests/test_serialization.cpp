#include "gtest/gtest.h"

#include <fstream>
#include <brwt/bit_vector.h>
#include <brwt/int_vector.h>

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

TEST(int_vector, SerializeAndLoad) {
    brwt::int_vector seq = {10, 20, 30, 40, 50};

    std::ofstream outstream(test_dump_basename, std::ios::binary);
    seq.serialize(outstream);
    outstream.close();

    std::ifstream instream(test_dump_basename, std::ios::binary);
    brwt::int_vector seq_loaded;
    EXPECT_TRUE(seq_loaded.load(instream));
    instream.close();

    EXPECT_TRUE(seq.size() == seq_loaded.size());
    EXPECT_TRUE(seq.get_bpe() == seq_loaded.get_bpe());
    EXPECT_TRUE(seq == seq_loaded);
}

TEST(bit_vector, SerializeAndLoad) {
    brwt::bit_vector v(100, 0xFF00'4FF4'FF11'33AA);
    std::ofstream outstream(test_dump_basename, std::ios::binary);
    v.serialize(outstream);
    outstream.close();

    std::ifstream instream(test_dump_basename, std::ios::binary);
    brwt::bit_vector v_loaded;
    EXPECT_TRUE(v_loaded.load(instream));
    instream.close();

    EXPECT_EQ(v_loaded.length(), v.length());
    EXPECT_EQ(v_loaded.size(), v.size());
    EXPECT_EQ(v_loaded.num_blocks(), v.num_blocks());
    EXPECT_EQ(v_loaded.allocated_bytes(), v.allocated_bytes());

    for (brwt::bit_vector::size_type i = 0; i < v.num_blocks(); ++i) {
        EXPECT_EQ(v_loaded.get_block(i), v.get_block(i));
    }
}

} // namespace
