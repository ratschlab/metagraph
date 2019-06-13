#include <random>

#include "gtest/gtest.h"
#include "test_helpers.hpp"

#define protected public
#define private public

#include "bitmap.hpp"
#include "bit_vector.hpp"
#include "threading.hpp"


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

template <class Bitmap, class Data>
std::unique_ptr<bitmap> build_bitmap(size_t size, const Data &data) {
    return std::unique_ptr<bitmap>(new Bitmap(size, data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_vector, sdsl::bit_vector>(size_t, const sdsl::bit_vector &data) {
    return std::unique_ptr<bitmap>(new bitmap_vector(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_vector,
             std::initializer_list<bool>>(size_t,
                                          const std::initializer_list<bool> &data) {
    return std::unique_ptr<bitmap>(new bitmap_vector(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_vector, std::set<uint64_t>>(size_t size,
                                                const std::set<uint64_t> &data) {
    sdsl::bit_vector vector(size, false);
    for (uint64_t i : data) {
        vector[i] = true;
    }

    return std::unique_ptr<bitmap>(new bitmap_vector(vector));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_vector, std::initializer_list<uint64_t>>(
      size_t size,
      const std::initializer_list<uint64_t> &data) {
    return build_bitmap<bitmap_vector>(size, std::set<uint64_t>(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_set, sdsl::bit_vector>(size_t size, const sdsl::bit_vector &data) {
    assert(size == data.size());

    std::set<uint64_t> bits;
    uint64_t i = 0;
    for (bool b : data) {
        if (b)
            bits.insert(i);

        ++i;
    }

    return std::unique_ptr<bitmap>(new bitmap_set(size, bits));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_set,
             std::initializer_list<bool>>(size_t size,
                                          const std::initializer_list<bool> &data) {
    assert(size == data.size());

    return build_bitmap<bitmap_set>(size, sdsl::bit_vector(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_adaptive, sdsl::bit_vector>(size_t, const sdsl::bit_vector &data) {
    return std::unique_ptr<bitmap>(new bitmap_adaptive(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_adaptive,
             std::initializer_list<bool>>(size_t,
                                          const std::initializer_list<bool> &data) {
    return std::unique_ptr<bitmap>(new bitmap_adaptive(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_lazy, sdsl::bit_vector>(size_t size, const sdsl::bit_vector &data) {
    assert(size == data.size());

    return std::unique_ptr<bitmap>(new bitmap_lazy(
        [data](auto i) { return data[i]; },
        size,
        std::accumulate(data.begin(), data.end(), 0u)
    ));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_lazy,
             std::initializer_list<bool>>(size_t size,
                                          const std::initializer_list<bool> &data) {
    return build_bitmap<bitmap_lazy>(size, sdsl::bit_vector(data));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_lazy,
             std::set<uint64_t>>(size_t size, const std::set<uint64_t> &data) {
    return std::unique_ptr<bitmap>(new bitmap_lazy(
        [data](auto i) { return data.find(i) != data.end(); },
        size,
        data.size()
    ));
}

template <>
std::unique_ptr<bitmap>
build_bitmap<bitmap_lazy, std::initializer_list<uint64_t>>(
      size_t size,
      const std::initializer_list<uint64_t> &data) {
    return build_bitmap<bitmap_lazy>(size, std::set<uint64_t>(data));
}


template <typename Bitmap>
class BitmapTest : public ::testing::Test { };

template <typename Bitmap>
class BitmapDynTest : public BitmapTest<Bitmap> { };

typedef ::testing::Types<bitmap_vector,
                         bitmap_set,
                         bitmap_adaptive,
                         bitmap_lazy> BitmapTypes;

typedef ::testing::Types<bitmap_vector,
                         bitmap_set,
                         bitmap_adaptive> BitmapDynTypes;

TYPED_TEST_CASE(BitmapTest, BitmapTypes);
TYPED_TEST_CASE(BitmapDynTest, BitmapDynTypes);


TYPED_TEST(BitmapTest, queries) {
    test_bitmap_queries<TypeParam>();
}


TYPED_TEST(BitmapTest, list) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    auto vector = build_bitmap<TypeParam>(init_list.size(), init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    reference_based_test(*vector, numbers);
}

TYPED_TEST(BitmapTest, bits) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::initializer_list<uint64_t> init_bits = { 1, 3, 4, 5, 6, 9, 14, 15 };
    auto vector = build_bitmap<TypeParam>(init_list.size(), init_bits);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    reference_based_test(*vector, numbers);
}


TYPED_TEST(BitmapDynTest, set_list) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    auto vector = build_bitmap<TypeParam>(init_list.size(), init_list);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    test_bitmap_set(vector.get(), &numbers);
}

TYPED_TEST(BitmapDynTest, set_bits) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::initializer_list<uint64_t> init_bits = { 1, 3, 4, 5, 6, 9, 14, 15 };
    auto vector = build_bitmap<TypeParam>(init_list.size(), init_bits);
    std::vector<bool> numbers(init_list);
    ASSERT_TRUE(vector.get());

    test_bitmap_set(vector.get(), &numbers);
}


TYPED_TEST(BitmapTest, concurrent_reading) {
    ThreadPool thread_pool(3);
    sdsl::bit_vector bv(100'000, false);

    std::vector<bool> bits;

    for (size_t i = 0; i < bv.size(); ++i) {
        bits.push_back((i + (i * i) % 31) % 2);
        if (bits.back())
            bv[i] = true;
    }

    auto vector = build_bitmap<TypeParam>(bv.size(), bv);

    std::vector<uint64_t> indices(vector->size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 g(1);
    std::shuffle(indices.begin(), indices.end(), g);

    reference_based_test(*vector, bits);

    for (auto i : indices) {
        thread_pool.enqueue([&](auto i) {
            ASSERT_EQ(bits[i], (*vector)[i]);
        }, i);
    }

    thread_pool.join();
}

TYPED_TEST(BitmapDynTest, concurrent_reading) {
    ThreadPool thread_pool(3);
    TypeParam vector(100'000, false);

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

TYPED_TEST(BitmapTest, concurrent_reading_all_zero) {
    ThreadPool thread_pool(3);
    TypeParam vector(100'000, false);

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

TYPED_TEST(BitmapTest, concurrent_reading_all_ones) {
    ThreadPool thread_pool(3);
    TypeParam vector(100'000, true);

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

TYPED_TEST(BitmapDynTest, or_with_stat) {
    test_operator_OR<TypeParam, bit_vector_stat>();
}

TYPED_TEST(BitmapDynTest, or_with_dyn) {
    test_operator_OR<TypeParam, bit_vector_dyn>();
}

TYPED_TEST(BitmapDynTest, or_with_sd) {
    test_operator_OR<TypeParam, bit_vector_sd>();
}

TYPED_TEST(BitmapDynTest, or_with_smart) {
    test_operator_OR<TypeParam, bit_vector_smart>();
}

TYPED_TEST(BitmapDynTest, or_with_small) {
    test_operator_OR<TypeParam, bit_vector_small>();
}


TYPED_TEST(BitmapTest, call_ones_dense_vector) {
    sdsl::bit_vector vector(50, true);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        auto bmd = build_bitmap<TypeParam>(vector.size(), vector);
        sdsl::bit_vector copy(vector.size(), false);
        bmd->call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitmapTest, call_ones_dense_set) {
    sdsl::bit_vector vector(50, true);
    std::set<uint64_t> bits;
    for (uint64_t i = 0; i < 50; ++i) {
        bits.emplace_hint(bits.end(), i);
    }

    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3;
        if (!vector[i])
            bits.erase(i);

        auto bmd = build_bitmap<TypeParam>(vector.size(), bits);
        sdsl::bit_vector copy(vector.size(), false);
        bmd->call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitmapTest, call_ones_sparse_vector) {
    sdsl::bit_vector vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        auto bms = build_bitmap<TypeParam>(vector.size(), vector);
        sdsl::bit_vector copy(vector.size(), false);
        bms->call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitmapTest, call_ones_sparse_set) {
    sdsl::bit_vector vector(50);
    std::set<uint64_t> bits;
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        if (vector[i])
            bits.emplace_hint(bits.end(), i);

        auto bms = build_bitmap<TypeParam>(vector.size(), bits);
        sdsl::bit_vector copy(vector.size(), false);
        bms->call_ones([&copy](auto i) { copy[i] = true; });
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
