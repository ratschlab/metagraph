#include <random>

#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include "common/threads/threading.hpp"

#define private public
#include "common/vectors/bitmap.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vectors/vector_algorithm.hpp"


namespace {

using namespace mtg;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/bitmap_dump_test";


void reference_based_test(const bitmap &vector,
                          const sdsl::bit_vector &reference) {
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


void test_bitmap_set(bitmap_dyn *vector, std::vector<bool> *numbers) {
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


template <class Bitmap>
Bitmap build_bitmap(const sdsl::bit_vector &data) {
    if constexpr(std::is_same_v<Bitmap, bitmap_set>) {
        std::set<uint64_t> bits;
        uint64_t i = 0;
        for (bool b : data) {
            if (b)
                bits.insert(i);

            ++i;
        }

        return bitmap_set(data.size(), bits);

    } else if constexpr(std::is_same_v<Bitmap, bitmap_lazy>) {
        return bitmap_lazy([data](auto i) { return data[i]; },
                           data.size(),
                           std::accumulate(data.begin(), data.end(), 0u));
    } else {
        return Bitmap(data);
    }
}

template <class Bitmap>
Bitmap build_bitmap(const std::initializer_list<bool> &data) {
    return build_bitmap<Bitmap>(sdsl::bit_vector(data));
}

template <class Bitmap>
Bitmap build_bitmap(const std::vector<bool> &data) {
    return build_bitmap<Bitmap>(to_sdsl(data));
}

template <class Bitmap>
Bitmap build_bitmap(size_t size, const std::set<uint64_t> &data) {
    if constexpr(std::is_same_v<Bitmap, bitmap_vector>) {
        sdsl::bit_vector vector(size, false);
        for (uint64_t i : data) {
            vector[i] = true;
        }

        return bitmap_vector(vector);

    } else if constexpr(std::is_same_v<Bitmap, bitmap_lazy>) {
        return bitmap_lazy([data](auto i) { return data.find(i) != data.end(); },
                           size,
                           data.size());
    } else {
        return Bitmap(size, data);
    }
}

template <class Bitmap>
Bitmap build_bitmap(size_t size, const std::initializer_list<uint64_t> &bits) {
    return build_bitmap<Bitmap>(size, std::set<uint64_t>(bits));
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

TYPED_TEST_SUITE(BitmapTest, BitmapTypes);
TYPED_TEST_SUITE(BitmapDynTest, BitmapDynTypes);


TYPED_TEST(BitmapTest, queries) {
    test_bitmap_queries<TypeParam>();
}


TYPED_TEST(BitmapTest, list) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    auto vector = build_bitmap<TypeParam>(init_list);
    std::vector<bool> numbers(init_list);

    reference_based_test(vector, numbers);
}

TYPED_TEST(BitmapTest, bits) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::initializer_list<uint64_t> init_bits = { 1, 3, 4, 5, 6, 9, 14, 15 };
    auto vector = build_bitmap<TypeParam>(init_list.size(), init_bits);
    std::vector<bool> numbers(init_list);

    reference_based_test(vector, numbers);
}


TYPED_TEST(BitmapDynTest, set_list) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    auto vector = build_bitmap<TypeParam>(init_list);
    std::vector<bool> numbers(init_list);

    test_bitmap_set(&vector, &numbers);
}

TYPED_TEST(BitmapDynTest, set_bits) {
    std::initializer_list<bool> init_list = { 0, 1, 0, 1, 1, 1, 1, 0,
                                              0, 1, 0, 0, 0, 0, 1, 1 };
    std::initializer_list<uint64_t> init_bits = { 1, 3, 4, 5, 6, 9, 14, 15 };
    auto vector = build_bitmap<TypeParam>(init_list.size(), init_bits);
    std::vector<bool> numbers(init_list);

    test_bitmap_set(&vector, &numbers);
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

    auto vector = build_bitmap<TypeParam>(bv);

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
        auto bmd = build_bitmap<TypeParam>(vector);
        sdsl::bit_vector copy(vector.size(), false);
        bmd.call_ones([&copy](auto i) { copy[i] = true; });
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
        bmd.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}

TYPED_TEST(BitmapTest, call_ones_sparse_vector) {
    sdsl::bit_vector vector(50);
    for (uint64_t i = 0; i < 10; ++i) {
        vector[i] = i % 3 == 0;
        auto bms = build_bitmap<TypeParam>(vector);
        sdsl::bit_vector copy(vector.size(), false);
        bms.call_ones([&copy](auto i) { copy[i] = true; });
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
        bms.call_ones([&copy](auto i) { copy[i] = true; });
        ASSERT_EQ(copy, vector);
    }
}


std::vector<bool> insert_reference_impl(std::vector<bool> vector,
                                        const std::vector<uint64_t> &new_pos) {
    assert(std::is_sorted(new_pos.begin(), new_pos.end()));
    for (auto i : new_pos) {
        vector.insert(vector.begin() + i, 0);
    }
    return vector;
}

TYPED_TEST(BitmapDynTest, insert_zeros_to_empty) {
    std::vector<bool> vector_init(0);

    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, 4, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, 4, 5, }) }) {

        TypeParam bitmap = build_bitmap<TypeParam>(vector_init);
        bitmap.insert_zeros(new_pos);

        EXPECT_EQ(build_bitmap<TypeParam>(insert_reference_impl(vector_init, new_pos)), bitmap)
            << to_sdsl(insert_reference_impl(vector_init, new_pos));
    }
}

TYPED_TEST(BitmapDynTest, insert_zeros_to_all_zeros) {
    std::vector<bool> vector_init(50, 0);

    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 50, }),
                                 std::vector<uint64_t>({ 50, 51, }),
                                 std::vector<uint64_t>({ 50, 51, 52, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, 51, 52, 53, 54, 55, }) }) {

        TypeParam bitmap = build_bitmap<TypeParam>(vector_init);
        bitmap.insert_zeros(new_pos);

        EXPECT_EQ(build_bitmap<TypeParam>(insert_reference_impl(vector_init, new_pos)), bitmap)
            << to_sdsl(insert_reference_impl(vector_init, new_pos));
    }
}

TYPED_TEST(BitmapDynTest, insert_zeros_to_all_ones) {
    std::vector<bool> vector_init(50, 1);

    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 50, }),
                                 std::vector<uint64_t>({ 50, 51, }),
                                 std::vector<uint64_t>({ 50, 51, 52, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, 51, 52, 53, 54, 55, }) }) {

        TypeParam bitmap = build_bitmap<TypeParam>(vector_init);
        bitmap.insert_zeros(new_pos);

        EXPECT_EQ(build_bitmap<TypeParam>(insert_reference_impl(vector_init, new_pos)), bitmap)
            << to_sdsl(insert_reference_impl(vector_init, new_pos));
    }
}

TYPED_TEST(BitmapTest, operator_eq) {
    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 100000 }) {
        for (bool value : { false, true }) {
            const TypeParam bitmap(size, value);

            EXPECT_EQ(bitmap, bitmap_vector(size, value));
            EXPECT_EQ(bitmap, bitmap_set(size, value));
            EXPECT_EQ(bitmap, bitmap_adaptive(size, value));
            EXPECT_EQ(bitmap, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
            EXPECT_EQ(bitmap, bit_vector_stat(size, value));
            EXPECT_EQ(bitmap, bit_vector_dyn(size, value));
            EXPECT_EQ(bitmap, bit_vector_sd(size, value));
            EXPECT_EQ(bitmap, bit_vector_rrr<>(size, value));
            EXPECT_EQ(bitmap, bit_vector_il<>(size, value));
            EXPECT_EQ(bitmap, bit_vector_hyb<>(size, value));
            EXPECT_EQ(bitmap, bit_vector_small(size, value));
            EXPECT_EQ(bitmap, bit_vector_smart(size, value));
        }
    }
}

TYPED_TEST(BitmapTest, operator_neq) {
    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 100000 }) {
        bool value = false;

        const TypeParam bitmap(size + 1, false);

        EXPECT_NE(bitmap, bitmap_vector(size, value));
        EXPECT_NE(bitmap, bitmap_set(size, value));
        EXPECT_NE(bitmap, bitmap_adaptive(size, value));
        EXPECT_NE(bitmap, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
        EXPECT_NE(bitmap, bit_vector_stat(size, value));
        EXPECT_NE(bitmap, bit_vector_dyn(size, value));
        EXPECT_NE(bitmap, bit_vector_sd(size, value));
        EXPECT_NE(bitmap, bit_vector_rrr<>(size, value));
        EXPECT_NE(bitmap, bit_vector_il<>(size, value));
        EXPECT_NE(bitmap, bit_vector_hyb<>(size, value));
        EXPECT_NE(bitmap, bit_vector_small(size, value));
        EXPECT_NE(bitmap, bit_vector_smart(size, value));
    }

    for (uint64_t size : { 0, 10, 64, 120, 128, 1000, 10000, 100000 }) {
        bool value = false;

        const TypeParam bitmap(size + 1, true);

        EXPECT_NE(bitmap, bitmap_vector(size, value));
        EXPECT_NE(bitmap, bitmap_set(size, value));
        EXPECT_NE(bitmap, bitmap_adaptive(size, value));
        EXPECT_NE(bitmap, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
        EXPECT_NE(bitmap, bit_vector_stat(size, value));
        EXPECT_NE(bitmap, bit_vector_dyn(size, value));
        EXPECT_NE(bitmap, bit_vector_sd(size, value));
        EXPECT_NE(bitmap, bit_vector_rrr<>(size, value));
        EXPECT_NE(bitmap, bit_vector_il<>(size, value));
        EXPECT_NE(bitmap, bit_vector_hyb<>(size, value));
        EXPECT_NE(bitmap, bit_vector_small(size, value));
        EXPECT_NE(bitmap, bit_vector_smart(size, value));
    }

    for (uint64_t size : { 10, 64, 120, 128, 1000, 10000, 100000 }) {
        for (bool value : { false, true }) {
            const TypeParam bitmap(size, !value);

            EXPECT_NE(bitmap, bitmap_vector(size, value));
            EXPECT_NE(bitmap, bitmap_set(size, value));
            EXPECT_NE(bitmap, bitmap_adaptive(size, value));
            EXPECT_NE(bitmap, bitmap_lazy([value](uint64_t) { return value; }, size, size * value));
            EXPECT_NE(bitmap, bit_vector_stat(size, value));
            EXPECT_NE(bitmap, bit_vector_dyn(size, value));
            EXPECT_NE(bitmap, bit_vector_sd(size, value));
            EXPECT_NE(bitmap, bit_vector_rrr<>(size, value));
            EXPECT_NE(bitmap, bit_vector_il<>(size, value));
            EXPECT_NE(bitmap, bit_vector_hyb<>(size, value));
            EXPECT_NE(bitmap, bit_vector_small(size, value));
            EXPECT_NE(bitmap, bit_vector_smart(size, value));
        }
    }
}

TEST(BitmapAdaptiveTest, switch_type) {
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

TEST(sdsl_bit_vector, call_ones_all_zeros) {
    {
        sdsl::bit_vector bv;
        call_ones(bv, [](auto) { ASSERT_FALSE(true); });
    }
    for (size_t size : { 1, 10, 100, 1000, 10000 }) {
        sdsl::bit_vector bv(size, false);
        call_ones(bv, [](auto) { ASSERT_FALSE(true); });
    }
}

TEST(sdsl_bit_vector, call_ones_all_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);
        uint64_t count = 0;
        call_ones(bv, [&](auto i) { ASSERT_EQ(count, i); count++; });
        ASSERT_EQ(size, count);
    }
}

TEST(sdsl_bit_vector, call_ones_all_even_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 2) {
            bv[i] = true;
        }
        uint64_t count = 0;
        call_ones(bv, [&](auto i) { ASSERT_EQ(count * 2, i); count++; });
        ASSERT_EQ((size + 1) / 2, count);
    }
}

TEST(sdsl_bit_vector, call_ones_every_third_one) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = true;
        }
        uint64_t count = 0;
        call_ones(bv, [&](auto i) { ASSERT_EQ(count * 3, i); count++; });
        ASSERT_EQ((size + 2) / 3, count);
    }
}

TEST(sdsl_bit_vector, call_ones_every_third_zero) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = false;
        }
        uint64_t count = 0;
        call_ones(bv, [&](auto i) { ASSERT_TRUE(i % 3); count++; });
        ASSERT_EQ(size - (size + 2) / 3, count);
    }
}

TEST(sdsl_bit_vector, call_ones_in_range_all_zeros) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);

        call_ones(bv, 0, 0, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size / 2, size / 2, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size, size, [](auto) { ASSERT_FALSE(true); });

        call_ones(bv, 0, size, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, 0, size / 2, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size / 2, size, [](auto) { ASSERT_FALSE(true); });
    }
}

TEST(sdsl_bit_vector, call_ones_in_range_all_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);

        call_ones(bv, 0, 0, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size / 2, size / 2, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size, size, [](auto) { ASSERT_FALSE(true); });

        uint64_t count = 0;
        call_ones(bv, 0, size, [&](auto i) { ASSERT_EQ(count, i); count++; });
        ASSERT_EQ(size, count);

        count = 0;
        call_ones(bv, 0, size / 2, [&](auto i) { ASSERT_EQ(count, i); count++; });
        ASSERT_EQ(size / 2, count);

        call_ones(bv, size / 2, size, [&](auto i) { ASSERT_EQ(count, i); count++; });
        ASSERT_EQ(size, count);
    }
}

TEST(sdsl_bit_vector, call_ones_in_range_all_even_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 2) {
            bv[i] = true;
        }

        call_ones(bv, 0, 0, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size / 2, size / 2, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size, size, [](auto) { ASSERT_FALSE(true); });

        uint64_t count = 0;
        call_ones(bv, 0, size, [&](auto i) { ASSERT_EQ(count * 2, i); count++; });
        ASSERT_EQ((size + 1) / 2, count);

        count = 0;
        call_ones(bv, 0, size / 2, [&](auto i) { ASSERT_EQ(count * 2, i); count++; });
        ASSERT_EQ((size / 2 + 1) / 2, count);

        call_ones(bv, size / 2, size, [&](auto i) { ASSERT_EQ(count * 2, i); count++; });
        ASSERT_EQ((size + 1) / 2, count);
    }
}

TEST(sdsl_bit_vector, call_ones_in_range_every_third_one) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = true;
        }

        call_ones(bv, 0, 0, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size / 2, size / 2, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size, size, [](auto) { ASSERT_FALSE(true); });

        uint64_t count = 0;
        call_ones(bv, 0, size, [&](auto i) { ASSERT_EQ(count * 3, i); count++; });
        ASSERT_EQ((size + 2) / 3, count);

        count = 0;
        call_ones(bv, 0, size / 2, [&](auto i) { ASSERT_EQ(count * 3, i); count++; });
        ASSERT_EQ((size / 2 + 2) / 3, count);

        call_ones(bv, size / 2, size, [&](auto i) { ASSERT_EQ(count * 3, i); count++; });
        ASSERT_EQ((size + 2) / 3, count);
    }
}

TEST(sdsl_bit_vector, call_ones_in_range_every_third_zero) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = false;
        }

        call_ones(bv, 0, 0, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size / 2, size / 2, [](auto) { ASSERT_FALSE(true); });
        call_ones(bv, size, size, [](auto) { ASSERT_FALSE(true); });

        uint64_t count = 0;
        call_ones(bv, 0, size, [&](auto i) { ASSERT_TRUE(i % 3); count++; });
        ASSERT_EQ(size - (size + 2) / 3, count);

        count = 0;
        call_ones(bv, 0, size / 2, [&](auto i) { ASSERT_TRUE(i % 3); count++; });
        ASSERT_EQ(size / 2 - (size / 2 + 2) / 3, count);

        call_ones(bv, size / 2, size, [&](auto i) { ASSERT_TRUE(i % 3); count++; });
        ASSERT_EQ(size - (size + 2) / 3, count);
    }
}

void check_call_ones_call_zeros(const sdsl::bit_vector &bv) {
    uint64_t size = bv.size();

    auto flipped = bv;
    flipped.flip();
    ASSERT_EQ(bv.size() - sdsl::util::cnt_one_bits(bv),
              sdsl::util::cnt_one_bits(flipped));

    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse);
    }
    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, 0, 0, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, 0, 0, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse);
    }
    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, size / 2, size / 2, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, size / 2, size / 2, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse);
    }
    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, size, size, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, size, size, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse);
    }
    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, 0, size, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, 0, size, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse);
    }
    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, 0, size / 2, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, 0, size / 2, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse);
    }
    {
        std::vector<uint64_t> called_ones;
        call_ones(bv, size / 2, size, [&](auto i) { called_ones.push_back(i); });
        std::vector<uint64_t> called_zeros_inverse;
        call_zeros(flipped, size / 2, size, [&](auto i) { called_zeros_inverse.push_back(i); });
        ASSERT_EQ(called_ones, called_zeros_inverse) << size;
    }
}

TEST(sdsl_bit_vector, call_ones_call_zeros_all_zeros) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        check_call_ones_call_zeros(bv);
    }
}

TEST(sdsl_bit_vector, call_ones_call_zeros_all_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);
        check_call_ones_call_zeros(bv);
    }
}

TEST(sdsl_bit_vector, call_ones_call_zeros_all_even_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 2) {
            bv[i] = true;
        }
        check_call_ones_call_zeros(bv);
    }
}

TEST(sdsl_bit_vector, call_ones_call_zeros_every_third_one) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = true;
        }
        check_call_ones_call_zeros(bv);
    }
}

TEST(sdsl_bit_vector, call_ones_call_zeros_every_third_zero) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = false;
        }
        check_call_ones_call_zeros(bv);
    }
}

TEST(sdsl_bit_vector, count_ones_all_zeros) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);

        ASSERT_EQ(sdsl::util::cnt_one_bits(bv), count_ones(bv, 0, size));

        ASSERT_EQ(0u, count_ones(bv, 0, 0));
        ASSERT_EQ(0u, count_ones(bv, size / 2, size / 2));
        ASSERT_EQ(0u, count_ones(bv, size, size));

        ASSERT_EQ(0u, count_ones(bv, 0, size));
        ASSERT_EQ(0u, count_ones(bv, 0, size / 2));
        ASSERT_EQ(0u, count_ones(bv, size / 2, size));
    }
}

TEST(sdsl_bit_vector, count_ones_all_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);

        ASSERT_EQ(sdsl::util::cnt_one_bits(bv), count_ones(bv, 0, size));

        ASSERT_EQ(0u, count_ones(bv, 0, 0));
        ASSERT_EQ(0u, count_ones(bv, size / 2, size / 2));
        ASSERT_EQ(0u, count_ones(bv, size, size));

        ASSERT_EQ(size, count_ones(bv, 0, size));
        ASSERT_EQ(size / 2, count_ones(bv, 0, size / 2));
        ASSERT_EQ((size + 1) / 2, count_ones(bv, size / 2, size));
    }
}

TEST(sdsl_bit_vector, count_ones_all_even_ones) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 2) {
            bv[i] = true;
        }

        ASSERT_EQ(sdsl::util::cnt_one_bits(bv), count_ones(bv, 0, size));

        ASSERT_EQ(0u, count_ones(bv, 0, 0));
        ASSERT_EQ(0u, count_ones(bv, size / 2, size / 2));
        ASSERT_EQ(0u, count_ones(bv, size, size));

        ASSERT_EQ((size + 1) / 2, count_ones(bv, 0, size));
        ASSERT_EQ((size / 2 + 1) / 2, count_ones(bv, 0, size / 2));
        ASSERT_EQ((size + 1) / 2 - ((size / 2 + 1) / 2), count_ones(bv, size / 2, size));
    }
}

TEST(sdsl_bit_vector, count_ones_every_third_one) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, false);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = true;
        }

        ASSERT_EQ(sdsl::util::cnt_one_bits(bv), count_ones(bv, 0, size));

        ASSERT_EQ(0u, count_ones(bv, 0, 0));
        ASSERT_EQ(0u, count_ones(bv, size / 2, size / 2));
        ASSERT_EQ(0u, count_ones(bv, size, size));

        ASSERT_EQ((size + 2) / 3, count_ones(bv, 0, size));
        ASSERT_EQ((size / 2 + 2) / 3, count_ones(bv, 0, size / 2));
        ASSERT_EQ((size + 2) / 3 - ((size / 2 + 2) / 3), count_ones(bv, size / 2, size));
    }
}

TEST(sdsl_bit_vector, count_ones_every_third_zero) {
    for (size_t size : { 0, 1, 10, 11, 100, 1000, 10000, 10001 }) {
        sdsl::bit_vector bv(size, true);
        for (size_t i = 0; i < size; i += 3) {
            bv[i] = false;
        }

        ASSERT_EQ(sdsl::util::cnt_one_bits(bv), count_ones(bv, 0, size));

        ASSERT_EQ(0u, count_ones(bv, 0, 0));
        ASSERT_EQ(0u, count_ones(bv, size / 2, size / 2));
        ASSERT_EQ(0u, count_ones(bv, size, size));

        ASSERT_EQ(size - (size + 2) / 3, count_ones(bv, 0, size));
        ASSERT_EQ(size / 2 - (size / 2 + 2) / 3, count_ones(bv, 0, size / 2));
        ASSERT_EQ(size - (size + 2) / 3 - (size / 2 - (size / 2 + 2) / 3), count_ones(bv, size / 2, size));
    }
}

} // namespace
