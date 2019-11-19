#include "gtest/gtest.h"

#include <fstream>
#include <brwt/bit_vector.h>
#include <brwt/int_vector.h>

#include "serialization.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/vector_dump_test";
const std::string test_dump_basename_bad = test_data_dir + "/vector_dump_test_bad";


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


void test_random_vector_stream(size_t length) {
    for (size_t i = 0; i < 10; ++i) {
        std::vector<uint64_t> numbers(length);
        for (size_t j = 0; j < length; ++j) {
            numbers[j] = rand() % (length * 2);
        }

        VectorFileOutStream outstream(test_dump_basename);
        outstream.write_value(length * 2);
        outstream.write_value(length);

        for (uint64_t num : numbers) {
            outstream.write_value(num);
        }

        outstream.close();

        BitVectorFileInStream instream(test_dump_basename);
        ASSERT_EQ(length * 2, instream.length());
        ASSERT_EQ(length, instream.values_left());

        for (size_t j = 0; j < length; ++j) {
            ASSERT_EQ(length - j, instream.values_left());
            ASSERT_EQ(numbers[j], instream.next_value());
        }

        ASSERT_FALSE(instream.values_left());
    }
}

TEST(Serialization, SerializationRandomUInt64VectorStream0) {
    test_random_vector_stream(0);
}

TEST(Serialization, SerializationRandomUInt64VectorStream10) {
    test_random_vector_stream(10);
}

TEST(Serialization, SerializationRandomUInt64VectorStream256) {
    test_random_vector_stream(256);
}

TEST(Serialization, SerializationRandomUInt64VectorStream1000000) {
    test_random_vector_stream(1000000);
}

TEST(Serialization, SerializationVectorStreamBadRead) {
    ASSERT_THROW(std::make_unique<BitVectorFileInStream>(test_dump_basename_bad),
                 std::ifstream::failure);
}

TEST(Serialization, SerializationVectorStreamBadReadShortFile) {
    VectorFileOutStream outstream(test_dump_basename);
    outstream.write_value(100000);
    outstream.write_value(100000);

    for (size_t i = 0; i < 5; ++i) {
        outstream.write_value(i);
    }

    outstream.close();

    BitVectorFileInStream instream(test_dump_basename);
    for (size_t i = 0; i < 5; ++i) {
        ASSERT_EQ(100000 - i, instream.values_left());
        ASSERT_EQ(i, instream.next_value());
    }

    ASSERT_EQ(100000u - 5, instream.values_left());
    EXPECT_THROW(instream.next_value(), std::ifstream::failure);
}
