#include "common/hash/bloom_filter.hpp"
#include "../test_helpers.hpp"

#include <gtest/gtest.h>
#include <sdsl/int_vector.hpp>

#include <random>
#include <vector>


namespace {

using namespace mtg;

// Note: std::hash<uint64_t> is not required by the specification to be a
//       cryographic hash, so use this to generate random numbers as hashes
const struct {
    uint64_t operator()(uint64_t i) const {
        gen.seed(i);
        return gen();
    }

    mutable std::mt19937_64 gen;
} hasher;

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

    const uint64_t offset = ((__uint128_t(hash) * vector.size()) >> (64 + 9)) << 9;
    uint64_t h_hi = hash >> 32;
    uint64_t h_lo = hash & 0xFFFFFFFF;

    for (size_t i = 0; i < num_hash_functions; ++i) {
        vector[offset + ((h_hi + i * h_lo) & 511)] = true;
    }
}

bool is_present(const sdsl::bit_vector &vector, uint64_t hash, size_t num_hash_functions) {
    if (!vector.size())
        return true;

    const uint64_t offset = ((__uint128_t(hash) * vector.size()) >> (64 + 9)) << 9;
    uint64_t h_hi = hash >> 32;
    uint64_t h_lo = hash & 0xFFFFFFFF;

    for (size_t i = 0; i < num_hash_functions; ++i) {
        if (!vector[offset + ((h_hi + i * h_lo) & 511)])
            return false;
    }

    return true;
}

double expected_fpr(const BloomFilter &filter, size_t num_elements) {
    return std::pow(double(1.0) - std::pow(double(1.0) - double(1.0) / filter.size(),
                                           num_elements * filter.num_hash_functions()),
                    filter.num_hash_functions());
}

TEST(BloomFilter, insert_and_check) {
    constexpr uint64_t max_num_hash_functions = 10;

    for (double bits_per_element : { 0.1, 1.0, 5.0 }) {
        for (uint64_t num_elements : { 10000, 50000, 100000 }) {
            BloomFilter filter(std::ceil(bits_per_element * num_elements),
                               num_elements,
                               max_num_hash_functions);

            ASSERT_LE(std::ceil(bits_per_element * num_elements), filter.size());
            ASSERT_EQ(0u, filter.size() % 512);

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

TEST(BloomFilter, batch_insert_and_check) {
    constexpr uint64_t max_num_hash_functions = 10;

    for (double bits_per_element : { 0.1, 1.0, 5.0 }) {
        for (uint64_t num_elements : { 10000, 50000, 100000 }) {
            BloomFilter filter(std::ceil(bits_per_element * num_elements),
                               num_elements,
                               max_num_hash_functions);
            ASSERT_LE(std::ceil(bits_per_element * num_elements), filter.size());
            ASSERT_EQ(0u, filter.size() % 512);

            sdsl::bit_vector check(filter.size());

            std::vector<uint64_t> hashes(num_elements);
            for (size_t i = 0; i < num_elements; ++i) {
                hashes[i] = hasher(i);
                insert(check, hashes[i], filter.num_hash_functions());
            }

            filter.insert(hashes.data(), hashes.data() + num_elements);
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

TEST(BloomFilter, batch_insert_and_batch_check) {
    constexpr uint64_t max_num_hash_functions = 10;

    for (double bits_per_element : { 0.1, 1.0, 5.0 }) {
        for (uint64_t num_elements : { 10000, 50000, 100000 }) {
            BloomFilter filter(std::ceil(bits_per_element * num_elements),
                               num_elements,
                               max_num_hash_functions);
            ASSERT_LE(std::ceil(bits_per_element * num_elements), filter.size());
            ASSERT_EQ(0u, filter.size() % 512);

            sdsl::bit_vector check(filter.size());

            std::vector<uint64_t> hashes(num_elements);
            for (size_t i = 0; i < num_elements; ++i) {
                hashes[i] = hasher(i);
                insert(check, hashes[i], filter.num_hash_functions());
            }

            filter.insert(hashes.data(), hashes.data() + num_elements);
            ASSERT_EQ(check, filter.data());

            EXPECT_EQ(
                num_elements,
                sdsl::util::cnt_one_bits(
                    filter.check(hashes.data(),
                                 hashes.data() + num_elements)
                )
            );

            uint64_t false_positives = 0;
            sdsl::bit_vector checks(1000, false);
            std::vector<uint64_t> next_hashes;
            next_hashes.reserve(1000);
            for (size_t i = num_elements; i < num_elements + 1000; ++i) {
                next_hashes.emplace_back(hasher(i));
                checks[i - num_elements] = is_present(check,
                                                      next_hashes.back(),
                                                      filter.num_hash_functions());
                false_positives += checks[i - num_elements];
                EXPECT_EQ(checks[i - num_elements], filter.check(next_hashes.back()));
            }

            EXPECT_EQ(checks, filter.check(next_hashes.data(), next_hashes.data() + 1000));

            TEST_COUT << "Elements: " << num_elements << std::endl
                      << "Bloom filter: " << filter.size() << " bits; "
                      << filter.num_hash_functions() << " hashes" << std::endl
                      << "Expected FPR: " << expected_fpr(filter, num_elements) << std::endl
                      << "False positives: " << double(false_positives) / 1000 << "; "
                      << false_positives << " / 1000" << std::endl;
        }
    }
}

TEST(BloomFilter, batch_insert_and_batch_check_callback) {
    constexpr uint64_t max_num_hash_functions = 10;

    for (double bits_per_element : { 0.1, 1.0, 5.0 }) {
        for (uint64_t num_elements : { 10000, 50000, 100000 }) {
            BloomFilter filter(std::ceil(bits_per_element * num_elements),
                               num_elements,
                               max_num_hash_functions);
            ASSERT_LE(std::ceil(bits_per_element * num_elements), filter.size());
            ASSERT_EQ(0u, filter.size() % 512);

            sdsl::bit_vector check(filter.size());

            std::vector<uint64_t> hashes(num_elements);
            for (size_t i = 0; i < num_elements; ++i) {
                hashes[i] = hasher(i);
                insert(check, hashes[i], filter.num_hash_functions());
            }

            filter.insert(hashes.data(), hashes.data() + num_elements);
            ASSERT_EQ(check, filter.data());

            auto result = filter.check(hashes.data(), hashes.data() + num_elements);

            EXPECT_EQ(num_elements, sdsl::util::cnt_one_bits(result));

            size_t j = 0;
            filter.check(hashes.data(), hashes.data() + num_elements,
                         [&](size_t i) {
                             ASSERT_GT(result.size(), i);
                             while (j < i) {
                                 EXPECT_FALSE(result[j++]);
                             }

                             EXPECT_TRUE(result[j++]);
                         });

            uint64_t false_positives = 0;
            sdsl::bit_vector checks(1000, false);
            std::vector<uint64_t> next_hashes;
            next_hashes.reserve(1000);
            for (size_t i = num_elements; i < num_elements + 1000; ++i) {
                next_hashes.emplace_back(hasher(i));
                checks[i - num_elements] = is_present(check,
                                                      next_hashes.back(),
                                                      filter.num_hash_functions());
                false_positives += checks[i - num_elements];
                EXPECT_EQ(checks[i - num_elements], filter.check(next_hashes.back()));
            }

            EXPECT_EQ(checks, filter.check(next_hashes.data(), next_hashes.data() + 1000));

            TEST_COUT << "Elements: " << num_elements << std::endl
                      << "Bloom filter: " << filter.size() << " bits; "
                      << filter.num_hash_functions() << " hashes" << std::endl
                      << "Expected FPR: " << expected_fpr(filter, num_elements) << std::endl
                      << "False positives: " << double(false_positives) / 1000 << "; "
                      << false_positives << " / 1000" << std::endl;
        }
    }
}

} // namespace
