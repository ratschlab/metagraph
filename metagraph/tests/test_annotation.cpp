#include <stdio.h>
#include <random>

#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "config.hpp"
#include "annotate.hpp"


std::vector<sdsl::uint256_t> generate_kmers(size_t num) {
    std::vector<sdsl::uint256_t> kmers(num);
    int mod = pow(num, 0.25);
    for (size_t i = 0; i < kmers.size(); ++i) {
        *(reinterpret_cast<uint64_t*>(&kmers[i]))     = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 1) = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 2) = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 3) = rand() % mod;
    }
    return kmers;
}

TEST(Annotate, RandomTestNoFalseNegative) {
    //create annotation
    annotate::BloomFilter<7> bloom(1000);
    annotate::ExactFilter    exact;
    //generate a bunch of kmers
    auto kmers = generate_kmers(1000);
    size_t total = 0, fp = 0;
    for (size_t i = 0; i < kmers.size(); ++i) {
        if (i < kmers.size() / 2) {
            bloom.insert(&kmers[i], &kmers[i] + 1);
            exact.insert(&kmers[i], &kmers[i] + 1);
            ASSERT_TRUE(bloom.find(&kmers[i], &kmers[i] + 1));
            ASSERT_TRUE(exact.find(&kmers[i], &kmers[i] + 1));
        } else {
            if (exact.find(&kmers[i], &kmers[i] + 1)) {
                ASSERT_TRUE(bloom.find(&kmers[i], &kmers[i] + 1));
            }
            if (!exact.find(&kmers[i], &kmers[i] + 1)) {
                total++;
                if (bloom.find(&kmers[i], &kmers[i] + 1)) {
                    fp++;
                }
            }
        }
    }
    //check to make sure the Bloom filter isn't full
    std::cerr << "Total: " << total << " FP: " << fp << std::endl;
    EXPECT_TRUE(total > fp);
}

TEST(Annotate, RandomHashAnnotator) {
    annotate::HashAnnotation<annotate::BloomFilter<7>> bloomhash;
    annotate::HashAnnotation<annotate::ExactFilter>    exacthash;
    size_t num_bits = 5;
    std::vector<size_t> bounds(num_bits);
    std::iota(bounds.begin(), bounds.end(), 0);
    for (size_t i = 0; i < num_bits; ++i) {
        bloomhash.append_bit(1000);
        exacthash.append_bit();
    }
    ASSERT_EQ(bloomhash.size(), num_bits);
    ASSERT_EQ(exacthash.size(), num_bits);
    auto kmers = generate_kmers(1000);
    for (size_t i = 0; i < kmers.size(); ++i) {
        size_t pick_bits = 0;
        if (i < kmers.size()) {
            //insert into random bit positions
            pick_bits = rand() % (1u << num_bits);
            ASSERT_TRUE(pick_bits < (1u << num_bits));
            for (size_t j = 0; j < num_bits; ++j) {
                if (((1lu << j) | pick_bits) == pick_bits) {
                    //insert
                    auto testbloom = bloomhash.insert(&kmers[i], &kmers[i] + 1, j);
                    auto testexact = exacthash.insert(&kmers[i], &kmers[i] + 1, j);
                    //check if it's there
                    testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1, j);
                    testexact = exacthash.find(&kmers[i], &kmers[i] + 1, j);

                    //test OR
                    auto testbloom_merged = testbloom;
                    annotate::HashAnnotation<>::merge_or(testbloom_merged, testexact);
                    ASSERT_EQ(testbloom_merged[0], testbloom[0]);

                    //test bit
                    ASSERT_EQ(testbloom[0] | (1lu << j), testbloom[0]);
                    ASSERT_EQ(testexact[0] | (1lu << j), testexact[0]);

                    //test AND
                    auto testbloom_and = testbloom;
                    annotate::HashAnnotation<>::merge_and(testbloom_merged, testexact);
                    ASSERT_EQ(testbloom_and[0], testexact[0]);
                }
            }
        }
        auto testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1); 
        auto testexact = exacthash.find(&kmers[i], &kmers[i] + 1);
        auto it = testbloom.begin();
        auto jt = testexact.begin();
        for (; it != testbloom.end(); ++it, ++jt) {
            ASSERT_EQ(*it | *jt, *it);
        }
    }
}
