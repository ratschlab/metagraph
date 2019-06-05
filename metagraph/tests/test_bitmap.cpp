#include <random>

#include "gtest/gtest.h"

#define protected public
#define private public

#include "bitmap.hpp"
#include "bit_vector.hpp"
#include "threading.hpp"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(a, b) (void)0
#endif
#endif


const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bitmap_dump_test";


void reference_based_test(const bitmap &vector,
                          const sdsl::bit_vector &reference) {
    ASSERT_DEATH(vector[vector.size()], "");
    ASSERT_DEATH(vector[vector.size() + 1], "");

    EXPECT_EQ(reference.size(), vector.size());
    EXPECT_EQ(sdsl::util::cnt_one_bits(reference), vector.num_set_bits());

    size_t i = 0;
    for (; i + 64 <= vector.size(); i += 64) {
        EXPECT_EQ(reference.get_int(i), vector.get_int(i, 64));
    }
    if (i < vector.size()) {
        EXPECT_EQ(
            reference.get_int(i, vector.size() - i),
            vector.get_int(i, vector.size() - i)
        );
    }

    {
        sdsl::bit_vector copy(vector.size(), 0);
        vector.call_ones([&](auto i) { copy[i] = true; });

        EXPECT_EQ(reference, copy);
    }
    {
        sdsl::bit_vector copy(vector.size(), 0);
        vector.add_to(&copy);

        EXPECT_EQ(reference, copy);
    }
    {
        sdsl::bit_vector copy(vector.size(), 0);
        bitmap_set other(vector.size(), 0);

        other |= vector;
        other.add_to(&copy);

        EXPECT_EQ(reference, copy);
    }
    {
        sdsl::bit_vector copy(vector.size(), 0);
        bitmap_vector other(vector.size(), 0);

        other |= vector;
        other.add_to(&copy);

        EXPECT_EQ(reference, copy);
    }
    {
        sdsl::bit_vector copy(vector.size(), 0);
        bitmap_adaptive other(vector.size(), 0);

        other |= vector;
        other.add_to(&copy);

        EXPECT_EQ(reference, copy);
    }
}

void reference_based_test(const bitmap &vector,
                          const std::vector<bool> &reference) {
    reference_based_test(vector, to_sdsl(reference));
}


template <class T>
void test_bitmap_queries() {
    std::unique_ptr<bitmap> vector { new T() };
    ASSERT_TRUE(vector);

    vector.reset(new T(0, 1));
    ASSERT_TRUE(vector);

    vector.reset(new T(10, 0));
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(0, (*vector)[i]);
    }

    vector.reset(new T(10, 1));
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
    }
}


TEST(bitmap_set, queries) {
    test_bitmap_queries<bitmap_set>();
}

TEST(bitmap_vector, queries) {
    test_bitmap_queries<bitmap_vector>();
}

TEST(bitmap_adaptive, queries) {
    test_bitmap_queries<bitmap_adaptive>();
}


void test_bitmap_set(bitmap *vector, std::vector<bool> *numbers) {
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
        vector->set(i, 1);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = 0;
        vector->set(i, 0);
        reference_based_test(*vector, *numbers);

        numbers->at(i) = value;
        vector->set(i, value);
    }
}


TEST(bitmap_vector, set) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::unique_ptr<bitmap> vector{ new bitmap_vector(init_list) };
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    test_bitmap_set(vector.get(), &numbers);
}

TEST(bitmap_adaptive, set_list) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::unique_ptr<bitmap> vector{ new bitmap_adaptive(init_list) };
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    test_bitmap_set(vector.get(), &numbers);
}

TEST(bitmap_set, set) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::initializer_list<uint64_t> init_bits = { 1, 3, 4, 5, 6, 9, 14, 15 };
    std::unique_ptr<bitmap> vector{ new bitmap_set(init_list.size(), init_bits) };
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    test_bitmap_set(vector.get(), &numbers);
}

TEST(bitmap_adaptive, set_bits) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::initializer_list<uint64_t> init_bits = { 1, 3, 4, 5, 6, 9, 14, 15 };
    std::unique_ptr<bitmap> vector{
        new bitmap_adaptive(init_list.size(), init_bits)
    };
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    test_bitmap_set(vector.get(), &numbers);
}


TEST(bitmap_set, concurrent_reading) {
    ThreadPool thread_pool(3);
    bitmap_set vector(100'000, false);

    std::vector<bool> bits;

    for (size_t i = 0; i < vector.size(); ++i) {
        bits.push_back((i + (i * i) % 31) % 2);
        if (bits.back())
            vector.set(i, true);
    }

    std::vector<uint64_t> indices(vector.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 g(1);
    std::shuffle(indices.begin(), indices.end(), g);

    reference_based_test(vector, bits);

    for (auto i : indices) {
        thread_pool.enqueue([&](auto i) {
            ASSERT_EQ(bits[i], vector[i]);
        }, i);
    }

    thread_pool.join();
}

TEST(bitmap_vector, concurrent_reading) {
    ThreadPool thread_pool(3);
    bitmap_vector vector(100'000, false);

    std::vector<bool> bits;

    for (size_t i = 0; i < vector.size(); ++i) {
        bits.push_back((i + (i * i) % 31) % 2);
        if (bits.back())
            vector.set(i, true);
    }

    std::vector<uint64_t> indices(vector.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 g(1);
    std::shuffle(indices.begin(), indices.end(), g);

    reference_based_test(vector, bits);

    for (auto i : indices) {
        thread_pool.enqueue([&](auto i) {
            ASSERT_EQ(bits[i], vector[i]);
        }, i);
    }

    thread_pool.join();
}

TEST(bitmap_adaptive, concurrent_reading) {
    ThreadPool thread_pool(3);
    bitmap_adaptive vector(100'000, false);

    std::vector<bool> bits;

    for (size_t i = 0; i < vector.size(); ++i) {
        bits.push_back((i + (i * i) % 31) % 2);
        if (bits.back())
            vector.set(i, true);
    }

    std::vector<uint64_t> indices(vector.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 g(1);
    std::shuffle(indices.begin(), indices.end(), g);

    reference_based_test(vector, bits);

    for (auto i : indices) {
        thread_pool.enqueue([&](auto i) {
            ASSERT_EQ(bits[i], vector[i]);
        }, i);
    }

    thread_pool.join();
}

TEST(bitmap_adaptive, concurrent_reading_all_zero) {
    ThreadPool thread_pool(3);
    bitmap_adaptive vector(100'000, false);

    std::vector<bool> bits(vector.size(), false);

    std::vector<uint64_t> indices(vector.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 g(1);
    std::shuffle(indices.begin(), indices.end(), g);

    reference_based_test(vector, bits);

    for (auto i : indices) {
        thread_pool.enqueue([&](auto i) {
            ASSERT_EQ(bits[i], vector[i]);
        }, i);
    }

    thread_pool.join();
}

TEST(bitmap_adaptive, concurrent_reading_all_ones) {
    ThreadPool thread_pool(3);
    bitmap_adaptive vector(100'000, true);

    std::vector<bool> bits(vector.size(), true);

    std::vector<uint64_t> indices(vector.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 g(1);
    std::shuffle(indices.begin(), indices.end(), g);

    reference_based_test(vector, bits);

    for (auto i : indices) {
        thread_pool.enqueue([&](auto i) {
            ASSERT_EQ(bits[i], vector[i]);
        }, i);
    }

    thread_pool.join();
}

template <class Bitmap, class BitVector>
void test_operator_OR() {
    sdsl::bit_vector init(100'000, true);
    init[0] = false;
    init[init.size() / 3] = false;
    init[2 * init.size() / 3] = false;

    Bitmap vector(init.size(), false);

    BitVector other(init);

    vector |= other;

    reference_based_test(vector, other.to_vector());
}

TEST(bitmap_set, or_with_stat) {
    test_operator_OR<bitmap_set, bit_vector_stat>();
}
TEST(bitmap_set, or_with_dyn) {
    test_operator_OR<bitmap_set, bit_vector_dyn>();
}
TEST(bitmap_set, or_with_sd) {
    test_operator_OR<bitmap_set, bit_vector_sd>();
}
TEST(bitmap_set, or_with_smart) {
    test_operator_OR<bitmap_set, bit_vector_smart>();
}
TEST(bitmap_set, or_with_small) {
    test_operator_OR<bitmap_set, bit_vector_small>();
}

TEST(bitmap_vector, or_with_stat) {
    test_operator_OR<bitmap_vector, bit_vector_stat>();
}
TEST(bitmap_vector, or_with_dyn) {
    test_operator_OR<bitmap_vector, bit_vector_dyn>();
}
TEST(bitmap_vector, or_with_sd) {
    test_operator_OR<bitmap_vector, bit_vector_sd>();
}
TEST(bitmap_vector, or_with_smart) {
    test_operator_OR<bitmap_vector, bit_vector_smart>();
}
TEST(bitmap_vector, or_with_small) {
    test_operator_OR<bitmap_vector, bit_vector_small>();
}

TEST(bitmap_adaptive, or_with_stat) {
    test_operator_OR<bitmap_adaptive, bit_vector_stat>();
}
TEST(bitmap_adaptive, or_with_dyn) {
    test_operator_OR<bitmap_adaptive, bit_vector_dyn>();
}
TEST(bitmap_adaptive, or_with_sd) {
    test_operator_OR<bitmap_adaptive, bit_vector_sd>();
}
TEST(bitmap_adaptive, or_with_smart) {
    test_operator_OR<bitmap_adaptive, bit_vector_smart>();
}
TEST(bitmap_adaptive, or_with_small) {
    test_operator_OR<bitmap_adaptive, bit_vector_small>();
}


TEST(bitmap_vector, call_ones_dense) {
    sdsl::bit_vector vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bitmap_vector bmd(vector);
        sdsl::bit_vector copy(vector.size(), false);
        bmd.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_set, call_ones_dense) {
    sdsl::bit_vector vector(50, true);
    std::set<uint64_t> bits;
    for (uint64_t i = 0; i < 50; ++i) {
        bits.emplace_hint(bits.end(), i);
    }

    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        if (!vector[i])
            bits.erase(i);

        bitmap_set bmd(vector.size(), bits);
        sdsl::bit_vector copy(vector.size(), false);
        bmd.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_adaptive, call_ones_dense_vector) {
    sdsl::bit_vector vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        bitmap_adaptive bmd(vector);
        sdsl::bit_vector copy(vector.size(), false);
        bmd.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_set, call_ones_dense_set) {
    sdsl::bit_vector vector(50, true);
    std::set<uint64_t> bits;
    for (uint64_t i = 0; i < 50; ++i) {
        bits.emplace_hint(bits.end(), i);
    }

    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        if (!vector[i])
            bits.erase(i);

        bitmap_adaptive bmd(vector.size(), bits);
        sdsl::bit_vector copy(vector.size(), false);
        bmd.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_vector, call_ones_sparse) {
    sdsl::bit_vector vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bitmap_vector bms(vector);
        sdsl::bit_vector copy(vector.size(), false);
        bms.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_set, call_ones_sparse) {
    sdsl::bit_vector vector(50);
    std::set<uint64_t> bits;
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        if (vector[i])
            bits.emplace_hint(bits.end(), i);

        bitmap_set bms(vector.size(), bits);
        sdsl::bit_vector copy(vector.size(), false);
        bms.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_adaptive, call_ones_sparse_vector) {
    sdsl::bit_vector vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        bitmap_adaptive bms(vector);
        sdsl::bit_vector copy(vector.size(), false);
        bms.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_adaptive, call_ones_sparse_set) {
    sdsl::bit_vector vector(50);
    std::set<uint64_t> bits;
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        if (vector[i])
            bits.emplace_hint(bits.end(), i);

        bitmap_adaptive bms(vector.size(), bits);
        sdsl::bit_vector copy(vector.size(), false);
        bms.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TEST(bitmap_adaptive, switch_type) {
    bitmap_adaptive bm(bitmap_adaptive::kRowCutoff - 2);

    uint64_t i = 0;

    // start as vector
    EXPECT_TRUE(dynamic_cast<const bitmap_vector*>(&bm.data()))
        << bm.size() << " " << i << " " << bm.num_set_bits();

    for (; i + 1 < (bm.size() >> bitmap_adaptive::kMaxNumIndicesLogRatio); ++i) {
        bm.set(i, true);
        EXPECT_TRUE(dynamic_cast<const bitmap_vector*>(&bm.data()))
            << bm.size() << " " << i << " " << bm.num_set_bits();
    }

    // too many bits set to switch back to set after increasing size
    bm.insert_zeros({ bm.size() - 1 });
    EXPECT_TRUE(dynamic_cast<const bitmap_vector*>(&bm.data()))
        << bm.size() << " " << i << " " << bm.num_set_bits();

    bm.insert_zeros({ bm.size() - 1 });
    EXPECT_TRUE(dynamic_cast<const bitmap_vector*>(&bm.data()))
        << bm.size() << " " << i << " " << bm.num_set_bits();

    // start unsetting bits
    ASSERT_LT(0u, i);
    for (--i; bm.num_set_bits()
            > (bm.size() >> (bitmap_adaptive::kMaxNumIndicesLogRatio
                + bitmap_adaptive::kNumIndicesMargin)); --i) {
        ASSERT_LT(0u, i);
        bm.set(i, false);
        EXPECT_TRUE(dynamic_cast<const bitmap_vector*>(&bm.data()))
            << bm.size() << " " << i << " " << bm.num_set_bits();
    }

    // switch to set
    bm.set(i, false);
    EXPECT_TRUE(dynamic_cast<const bitmap_set*>(&bm.data()))
        << bm.size() << " " << i << " " << bm.num_set_bits();

    // start resetting bits back to 1
    for (; bm.num_set_bits() + 1
            < (bm.size() >> bitmap_adaptive::kMaxNumIndicesLogRatio); ++i) {
        bm.set(i, true);
        EXPECT_TRUE(dynamic_cast<const bitmap_set*>(&bm.data()))
            << bm.size() << " " << i << " " << bm.num_set_bits();
    }

    //switch back to vector
    bm.set(i, true);
    EXPECT_TRUE(dynamic_cast<const bitmap_vector*>(&bm.data()))
        << bm.size() << " " << i << " " << bm.num_set_bits();
}
