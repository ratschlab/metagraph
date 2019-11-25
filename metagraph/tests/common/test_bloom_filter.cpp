#include "gtest/gtest.h"

#include "common/hash/bloom_filter.hpp"
#include "../test_helpers.hpp"

#include <functional>
#include <vector>

#include <sdsl/int_vector.hpp>


const std::hash<size_t> hasher;

TEST(BloomFilter, empty) {
    BloomFilter filter(0, 0, 4);

    // Bloom filter sizes are multiples of 512, empty filters are not allowed
    ASSERT_EQ(512u, filter.size());

    for (size_t i = 0; i < 10000; ++i) {
        EXPECT_FALSE(filter.check(hasher(i)));
    }
}

void insert(sdsl::bit_vector &vector, uint64_t hash, size_t num_hash_functions) {
    if (!vector.size())
        return;

    const auto offset = ((hash % vector.size()) >> 9) << 9;
    uint64_t base = hash, jump = hash >> 32;

    for (size_t i = 0; i < num_hash_functions; ++i) {
        vector[offset + ((base + i * jump) & 511)] = true;
    }
}

bool is_present(const sdsl::bit_vector &vector, uint64_t hash, size_t num_hash_functions) {
    if (!vector.size())
        return true;

    const auto offset = ((hash % vector.size()) >> 9) << 9;
    uint64_t base = hash, jump = hash >> 32;

    for (size_t i = 0; i < num_hash_functions; ++i) {
        if (!vector[offset + ((base + i * jump) & 511)])
            return false;
    }

    return true;
}

double expected_fpr(const BloomFilter &filter, size_t num_elements) {
    return std::pow(double(1.0) - std::pow(double(1.0) - double(1.0) / filter.size(),
                                           num_elements * filter.num_hash_functions()),
                    filter.num_hash_functions());
}

TEST(BloomFilter, check_set_bits) {
    for (double bits_per_element : { 0.1, 1.0, 5.0 }) {
        for (uint64_t num_elements : { 10000, 50000, 100000 }) {
            BloomFilter filter(std::ceil(bits_per_element * num_elements), num_elements, 10);
            sdsl::bit_vector check(filter.size());

            for (size_t i = 0; i < num_elements; ++i) {
                filter.insert(hasher(i));
                insert(check, hasher(i), filter.num_hash_functions());
                ASSERT_EQ(check, filter.data()) << i;
            }

            for (size_t i = 0; i < num_elements; ++i) {
                EXPECT_TRUE(filter.check(hasher(i)));
            }

            uint64_t false_positives = 0;
            for (size_t i = num_elements; i < num_elements + 1000; ++i) {
                auto hash = hasher(i);
                false_positives += is_present(check, hash, filter.num_hash_functions());
                EXPECT_EQ(is_present(check, hash, filter.num_hash_functions()),
                          filter.check(hash));
            }

            TEST_COUT << "Elements: " << num_elements << std::endl
                      << "Bloom filter: " << filter.size() << " bits; "
                      << filter.num_hash_functions() << " hashes" << std::endl
                      << "Expected FPR: " << expected_fpr(filter, num_elements) << std::endl
                      << "False positives: " << double(false_positives) / 1000 << "; "
                      << false_positives << " / 1000" << std::endl;
        }
    }
}

TEST(BloomFilter, check_set_bits_batch) {
    for (double bits_per_element : { 0.1, 1.0, 5.0 }) {
        for (uint64_t num_elements : { 10000, 50000, 100000 }) {
            BloomFilter filter(std::ceil(bits_per_element * num_elements), num_elements, 10);
            sdsl::bit_vector check(filter.size());

            std::vector<uint64_t> hashes(num_elements);
            for (size_t i = 0; i < num_elements; ++i) {
                hashes[i] = hasher(i);
                insert(check, hashes[i], filter.num_hash_functions());
            }

            filter.batch_insert(hashes.data(), num_elements);
            ASSERT_EQ(check, filter.data());

            for (size_t i = 0; i < num_elements; ++i) {
                EXPECT_TRUE(filter.check(hashes[i]));
            }

            uint64_t false_positives = 0;
            for (size_t i = num_elements; i < num_elements + 1000; ++i) {
                auto hash = hasher(i);
                false_positives += is_present(check, hash, filter.num_hash_functions());
                EXPECT_EQ(is_present(check, hash, filter.num_hash_functions()),
                          filter.check(hash));
            }

            TEST_COUT << "Elements: " << num_elements << std::endl
                      << "Bloom filter: " << filter.size() << " bits; "
                      << filter.num_hash_functions() << " hashes" << std::endl
                      << "Expected FPR: " << expected_fpr(filter, num_elements) << std::endl
                      << "False positives: " << double(false_positives) / 1000 << "; "
                      << false_positives << " / 1000" << std::endl;
        }
    }
}
