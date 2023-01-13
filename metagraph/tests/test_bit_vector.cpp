#include "gtest/gtest.h"

#include <cmath>
#include <cstdlib>

#include "test_helpers.hpp"

#include "common/vectors/sd_vector_disk.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/threads/threading.hpp"
#include "common/data_generation.hpp"


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bit_vector_dump_test";


template <typename Bitmap>
class BitVectorTest : public ::testing::Test { };

typedef ::testing::Types<bit_vector_stat,
                         bit_vector_dyn,
                         bit_vector_sd,
                         bit_vector_rrr<>,
                         bit_vector_il<>,
                         bit_vector_il<4096>,
                         bit_vector_hyb<>,
                         bit_vector_small,
                         bit_vector_smallrank,
                         bit_vector_smart>
        BitVectorTypes;

TYPED_TEST_SUITE(BitVectorTest, BitVectorTypes);


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
    ASSERT_DEBUG_DEATH(vector.next1(vector.size()), "");
    ASSERT_DEBUG_DEATH(vector.next1(vector.size() + 1), "");
    ASSERT_DEBUG_DEATH(vector.next1(vector.size() * 2), "");

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
    ASSERT_DEBUG_DEATH(vector.prev1(vector.size()), "");
    ASSERT_DEBUG_DEATH(vector.prev1(vector.size() + 1), "");
    ASSERT_DEBUG_DEATH(vector.prev1(vector.size() * 2), "");

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

    ASSERT_DEBUG_DEATH(vector.select1(0), "");

    for (size_t i : { 1, 2, 10, 100, 1000 }) {
        ASSERT_DEBUG_DEATH(vector.select1(max_rank + i), "");
        EXPECT_EQ(max_rank, vector.rank1(vector.size() + i - 2))
            << bit_vector_stat(reference);
    }
    ASSERT_DEBUG_DEATH(vector.select1(vector.size() + 1), "");
    ASSERT_DEBUG_DEATH(vector[vector.size()], "");
    ASSERT_DEBUG_DEATH(vector[vector.size() + 1], "");

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
        uint64_t cond_rk = vector.conditional_rank1(i);
        if (cond_rk) {
            EXPECT_TRUE(reference[i]);
            EXPECT_EQ(vector.rank1(i), cond_rk);
        } else {
            EXPECT_FALSE(reference[i]);
        }
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
        ASSERT_DEBUG_DEATH(vector->select1(i), "");
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
    ASSERT_DEBUG_DEATH(vector->select1(1'000), "");
    ASSERT_DEBUG_DEATH(vector->select1(0), "");

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
    for (size_t size : { 1, 2, 3, 4, 5, 50, 51, 52, 54, 100, 200, 300, 1000 }) {
        for (size_t i = 0; i < size; ++i) {
            sdsl::bit_vector bv(size, 0);
            bv[i] = 1;
            TypeParam bit_vector(std::move(bv));
            EXPECT_EQ(i, bit_vector.select1(1));
        }
    }
}

TYPED_TEST(BitVectorTest, select0) {
    // Mainly test select0.
    auto vector = std::make_unique<TypeParam>(10, 1);
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
        ASSERT_DEBUG_DEATH(vector->select0(i), "");
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

void test_bit_vector_set(bit_vector_dyn *vector, sdsl::bit_vector *numbers) {
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
    bit_vector_dyn vector(numbers);

    test_bit_vector_set(&vector, &numbers);
}


void test_bit_vector_ins_del(bit_vector_dyn *vector,
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
    bit_vector_dyn vector(numbers);

    test_bit_vector_ins_del(&vector, numbers);
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

template <class bit_vector_type>
uint64_t space_taken(const bit_vector_type &vec) {
    std::ofstream outstream(test_dump_basename, std::ios::binary);
    vec.serialize(outstream);
    outstream.close();
    std::ifstream instream(test_dump_basename, std::ios::binary | std::ifstream::ate);
    return instream.tellg() * 8;
}

TYPED_TEST(BitVectorTest, PredictedMemoryFootprint) {
    double tolerance = 0.01;

    if constexpr(std::is_same_v<TypeParam, bit_vector_hyb<>>)
        return;

    DataGenerator gen;
    for (uint64_t size : { 1'000'000, 10'000'000 }) {
        for (double density : { .05, .2, .4, .5, .7, .9, .95 }) {
            sdsl::bit_vector bv = gen.generate_random_column(size, density);
            uint64_t footprint = space_taken(TypeParam(bv));
            EXPECT_GE(TypeParam::predict_size(bv.size(), sdsl::util::cnt_one_bits(bv)),
                      footprint * (1 - tolerance)) << "Density: " << density;
            EXPECT_LE(TypeParam::predict_size(bv.size(), sdsl::util::cnt_one_bits(bv)),
                      footprint * (1 + tolerance)) << "Density: " << density;
        }
    }
}

TEST(select_support_mcl, PredictedMemoryFootprint) {
    double tolerance = 0.01;

    DataGenerator gen;
    for (uint64_t size : { 1'000'000, 10'000'000, 100'000'000 }) {
        for (double density : { .05, .2, .4, .5, .7, .9, .98 }) {
            sdsl::bit_vector bv = gen.generate_random_column(size, density);

            uint64_t footprint = space_taken(sdsl::select_support_mcl<1>(&bv));
            EXPECT_GE(footprint_select_support_mcl(bv.size(), sdsl::util::cnt_one_bits(bv)),
                      footprint * (1 - tolerance)) << "Size: " << size << "\tDensity: " << density;
            EXPECT_LE(footprint_select_support_mcl(bv.size(), sdsl::util::cnt_one_bits(bv)),
                      footprint * (1 + tolerance)) << "Size: " << size << "\tDensity: " << density;
        }
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

sdsl::sd_vector_disk<>
construct_sd_vector_disk(sdsl::bit_vector bv, const std::string &fname, size_t offset) {
    sdsl::sd_vector_disk_builder builder(bv.size(), sdsl::util::cnt_one_bits(bv), fname, offset);
    call_ones(bv, [&](uint64_t i) { builder.set(i); });
    return sdsl::sd_vector_disk<>(builder);
}

void compare(const sdsl::sd_vector<> &first, const sdsl::sd_vector_disk<> &second) {
    ASSERT_EQ(first.size(), second.size());

    sdsl::sd_vector<>::select_1_type first_slct(&first);
    sdsl::sd_vector<>::rank_1_type first_rank(&first);

    sdsl::sd_vector_disk<>::select_1_type second_slct(&second);
    sdsl::sd_vector_disk<>::rank_1_type second_rank(&second);

    uint64_t m = first_rank(first.size());

    ASSERT_EQ(first_rank(0), second_rank(0));

    // sequential
    for (size_t i = 0; i < std::min(first.size(), (uint64_t)1000); ++i) {
        ASSERT_EQ(first[i], second[i]);
        ASSERT_EQ(first_rank(i), second_rank(i));
        if (i < m) {
            ASSERT_EQ(first_slct(i + 1), second_slct(i + 1));
        }
    }

    // jumping indexes
    for (size_t i = 0; i < std::min((first.size() + 1) / 2, (uint64_t)1000); ++i) {
        uint64_t j = i;
        ASSERT_EQ(first[j], second[j]);
        ASSERT_EQ(first_rank(j), second_rank(j));
        if (j < m) {
            ASSERT_EQ(first_slct(j + 1), second_slct(j + 1));
        }

        j = (first.size() / 2 + i) % first.size();
        ASSERT_EQ(first[j], second[j]);
        ASSERT_EQ(first_rank(j), second_rank(j));
        if (j < m) {
            ASSERT_EQ(first_slct(j + 1), second_slct(j + 1));
        }
    }
}

TEST(sd_vector_disk, all_tests) {
    std::vector<sdsl::bit_vector> vectors = {
        sdsl::bit_vector({}),
        sdsl::bit_vector({ 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0 }),
        sdsl::bit_vector({ 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 }),
        sdsl::bit_vector({ 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1 }),
        sdsl::bit_vector({ 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1 }),
        sdsl::bit_vector({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 })
    };
    for (size_t len = 1; len < 100; ++len) {
        vectors.emplace_back(len, 0);
        vectors.emplace_back(len, 1);
    }
    vectors.emplace_back(100'000'000, 0);
    vectors.emplace_back(100'000'000, 1);

    for (const sdsl::bit_vector &numbers : vectors) {
        sdsl::sd_vector<> sd_vec(numbers);

        for (size_t offset : { 0, 777 }) {
            {
                auto fname = test_dump_basename + "_sd_disk";
                std::ofstream out(fname, std::ios::binary);
                std::vector<char> buffer(offset, 'X');
                out.write(buffer.data(), offset);
                out.close();

                sdsl::sd_vector_disk<> sd_vec_disk
                        = construct_sd_vector_disk(numbers, fname, offset);
                compare(sd_vec, sd_vec_disk);

                std::ifstream in(fname, std::ios::binary);
                buffer.assign(offset, '-');
                in.read(buffer.data(), offset);
                ASSERT_EQ(std::vector<char>(offset, 'X'), buffer);
            }

            {
                // test load
                auto fname = test_dump_basename + "_sd";
                std::ofstream out(fname, std::ios::binary);
                std::vector<char> buffer(offset, 'X');
                out.write(buffer.data(), offset);
                sd_vec.serialize(out);
                out.close();

                sdsl::sd_vector_disk<> sd_vec_disk;
                sd_vec_disk.load(fname, offset);
                compare(sd_vec, sd_vec_disk);

                std::ifstream in(fname, std::ios::binary);
                buffer.assign(offset, '-');
                in.read(buffer.data(), offset);
                ASSERT_EQ(std::vector<char>(offset, 'X'), buffer);
            }
        }
    }
}


TYPED_TEST(BitVectorTest, add_to_all_zero) {
    for (uint64_t size : { 0, 10, 100, 1000000 }) {
        DataGenerator gen;
        for (double density : { 0.0, 0.5, 1.0 }) {
            sdsl::bit_vector vector(size, false);
            TypeParam bvs(gen.generate_random_column(size, density));
            bvs.add_to(&vector);
            EXPECT_EQ(bvs.to_vector(), vector);
        }
    }
}

TYPED_TEST(BitVectorTest, add_to_all_one) {
    for (uint64_t size : { 0, 10, 100, 1000000 }) {
        DataGenerator gen;
        for (double density : { 0.0, 0.5, 1.0 }) {
            sdsl::bit_vector vector(size, true);
            TypeParam bvs(gen.generate_random_column(size, density));
            bvs.add_to(&vector);
            EXPECT_EQ(sdsl::bit_vector(size, true), vector);
        }
    }
}

TYPED_TEST(BitVectorTest, add_to_same) {
    for (uint64_t size : { 0, 10, 100, 1000000 }) {
        DataGenerator gen;
        for (double density : { 0.0, 0.5, 1.0 }) {
            sdsl::bit_vector vector = gen.generate_random_column(size, density);
            TypeParam bvs(vector);
            bvs.add_to(&vector);
            EXPECT_EQ(bvs.to_vector(), vector);
        }
    }
}

TYPED_TEST(BitVectorTest, add_all_zero) {
    for (uint64_t size : { 0, 10, 100, 1000000 }) {
        DataGenerator gen;
        for (double density : { 0.0, 0.5, 1.0 }) {
            sdsl::bit_vector bv = gen.generate_random_column(size, density);
            sdsl::bit_vector vector = bv;
            TypeParam bvs(sdsl::bit_vector(size, false));
            bvs.add_to(&vector);
            EXPECT_EQ(bv, vector);
        }
    }
}

TYPED_TEST(BitVectorTest, add_all_ones) {
    for (uint64_t size : { 0, 10, 100, 1000000 }) {
        DataGenerator gen;
        for (double density : { 0.0, 0.5, 1.0 }) {
            sdsl::bit_vector vector = gen.generate_random_column(size, density);
            TypeParam bvs(sdsl::bit_vector(size, true));
            bvs.add_to(&vector);
            EXPECT_EQ(sdsl::bit_vector(size, true), vector);
        }
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

} // namespace
