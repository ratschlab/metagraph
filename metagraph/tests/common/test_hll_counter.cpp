#include "common/hash/hll_counter.hpp"
#include "../test_helpers.hpp"

#include <gtest/gtest.h>
#include <sdsl/int_vector.hpp>

#include <random>
#include <vector>
#include <cmath>


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/hll_dump_test";

// Note: std::hash<uint64_t> is not required by the specification to be a
//       cryographic hash, so use this to generate random numbers as hashes
const struct {
    uint64_t operator()(uint64_t i) const {
        gen.seed(i);
        return gen();
    }

    mutable std::mt19937_64 gen;
} hasher;

void test_hll_serialize_load(const HLLCounter &counter) {
    {
        std::ofstream out(test_dump_basename + ".hll");
        counter.serialize(out);
    }

    HLLCounter loaded;
    {
        std::ifstream in(test_dump_basename + ".hll");
        loaded.load(in);
    }

    EXPECT_EQ(counter, loaded);
}

TEST(HLLCounter, empty) {
    HLLCounter counter(0.03);

    EXPECT_EQ(0.0, std::roundl(counter.estimate_cardinality()));

    for (size_t i = 0; i < 10000; ++i) {
        EXPECT_FALSE(counter.check(hasher(i)));
    }

    test_hll_serialize_load(counter);
}

TEST(HLLCounter, one_element) {
    for (size_t i = 0; i < 10000; ++i) {
        HLLCounter counter(0.03);
        counter.insert(hasher(i));
        EXPECT_EQ(1.0, std::roundl(counter.estimate_cardinality()));
        EXPECT_TRUE(counter.check(hasher(i)));

        if (i % 1000)
            continue;

        test_hll_serialize_load(counter);
    }
}

TEST(HLLCounter, many_elements) {
    HLLCounter counter(0.03);
    for (size_t i = 0; i < 10000; ++i) {
        if (i % 4)
            continue;

        counter.insert(hasher(i));
        EXPECT_TRUE(counter.check(hasher(i)));
    }
    EXPECT_GT(0.04, std::abs(2500.0 - std::roundl(counter.estimate_cardinality())) / 2500.0);
    test_hll_serialize_load(counter);
}

TEST(HLLCounter, many_elements_batch_check) {
    HLLCounter counter(0.03);
    std::vector<uint64_t> hashes;
    for (size_t i = 0; i < 10000; ++i) {
        hashes.push_back(hasher(i));
        if (i % 4)
            continue;

        counter.insert(hashes.back());
        EXPECT_TRUE(counter.check(hashes.back()));
    }

    sdsl::bit_vector check = counter.check(hashes.data(), hashes.data() + hashes.size());
    ASSERT_EQ(hashes.size(), check.size());
    for (size_t i = 0; i < hashes.size(); ++i) {
        EXPECT_EQ(counter.check(hashes[i]), check[i]);
    }

    EXPECT_GT(0.04, std::abs(2500.0 - std::roundl(counter.estimate_cardinality())) / 2500.0);
}

TEST(HLLCounter, more_elements) {
    HLLCounter counter(0.03);
    for (size_t i = 0; i < 10000; ++i) {
        counter.insert(hasher(i));
        EXPECT_TRUE(counter.check(hasher(i)));
    }
    EXPECT_GT(0.04, std::abs(10000.0 - std::roundl(counter.estimate_cardinality())) / 10000.0);
    test_hll_serialize_load(counter);
}

TEST(HLLCounter, more_elements_batch) {
    HLLCounter counter_batch(0.03);
    HLLCounter counter(0.03);
    std::vector<uint64_t> hashes;
    hashes.reserve(10000);
    for (size_t i = 0; i < 10000; ++i) {
        hashes.push_back(hasher(i));
        counter.insert(hasher(i));
    }

    counter_batch.insert(hashes.data(), hashes.data() + hashes.size());
    EXPECT_EQ(counter, counter_batch);
}

TEST(HLLCounter, more_elements_intersection) {
    HLLCounter counter_left(0.03);
    HLLCounter counter_right(0.03);
    std::vector<uint64_t> hashes;
    hashes.reserve(10000);
    for (size_t i = 0; i < 10000; ++i) {
        hashes.push_back(hasher(i));
    }

    counter_left.insert(hashes.data(), hashes.data() + 7500);
    counter_right.insert(hashes.data() + 2500, hashes.data() + hashes.size());
    EXPECT_EQ(counter_left.estimate_intersection_cardinality(counter_right),
              counter_right.estimate_intersection_cardinality(counter_left));
    EXPECT_GT(0.04, std::abs(5000.0 - std::roundl(counter_left.estimate_intersection_cardinality(counter_right))) / 5000.0);
}

TEST(HLLCounter, more_elements_union) {
    HLLCounter counter_left(0.03);
    HLLCounter counter_right(0.03);
    std::vector<uint64_t> hashes;
    hashes.reserve(10000);
    for (size_t i = 0; i < 10000; ++i) {
        hashes.push_back(hasher(i));
    }

    counter_left.insert(hashes.data(), hashes.data() + 7500);
    counter_right.insert(hashes.data() + 2500, hashes.data() + hashes.size());
    EXPECT_EQ(counter_left.estimate_union_cardinality(counter_right),
              counter_right.estimate_union_cardinality(counter_left));
    EXPECT_GT(0.04, std::abs(10000.0 - std::roundl(counter_left.estimate_union_cardinality(counter_right))) / 10000.0);
}

TEST(HLLCounter, more_elements_jaccard) {
    HLLCounter counter_left(0.03);
    HLLCounter counter_right(0.03);
    std::vector<uint64_t> hashes;
    hashes.reserve(10000);
    for (size_t i = 0; i < 10000; ++i) {
        hashes.push_back(hasher(i));
    }

    counter_left.insert(hashes.data(), hashes.data() + 7500);
    counter_right.insert(hashes.data() + 2500, hashes.data() + hashes.size());
    EXPECT_EQ(counter_left.estimate_jaccard(counter_right),
              counter_right.estimate_jaccard(counter_left));
    EXPECT_GT(0.04, std::abs(0.5 - counter_left.estimate_jaccard(counter_right)) / 0.5);
}

} // namespace
