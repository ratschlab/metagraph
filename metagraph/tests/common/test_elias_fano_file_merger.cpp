#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include <sdsl/uint128_t.hpp>

#include "common/elias_fano_file_merger.hpp"
#include "common/utils/file_utils.hpp"
#include "tests/utils/gtest_patch.hpp"


namespace {

using namespace mtg;

template <typename T>
class EliasFanoFileMergerTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t,
                         sdsl::uint128_t,
                         sdsl::uint256_t,
                         std::pair<uint64_t, uint32_t>>
        ValueTypes;

TYPED_TEST_SUITE(EliasFanoFileMergerTest, ValueTypes);


template <typename T>
void do_encode(const std::vector<T> &values, const std::string &file_name) {
    if (values.empty()) {
        common::EliasFanoEncoder<T> encoder(0, 0, 0, file_name);
        encoder.finish();
        return;
    }
    common::EliasFanoEncoder<T> encoder(values.size(), utils::get_first(values.front()),
                                        utils::get_first(values.back()), file_name);
    std::for_each(values.begin(), values.end(), [&encoder](const T &v) { encoder.add(v); });
    encoder.finish();
}

template <typename T>
std::vector<T>
get_random_values(uint32_t count,
                  std::mt19937 &rng,
                  std::uniform_int_distribution<std::mt19937::result_type> &dist) {
    std::vector<T> values(count);
    for (uint32_t j = 0; j < values.size(); ++j) {
        using F = utils::get_first_type_t<T>;
        F value = j > 0 ? utils::get_first(values[j - 1]) + (F)dist(rng) : (F)dist(rng);
        if constexpr (utils::is_pair_v<T>) {
            uint32_t count = dist(rng);
            values[j] = { value, count };
        } else {
            values[j] = value;
        }
    }
    return values;
}

template <typename T>
void push_back(std::vector<T> &v, const T &value) {
    v.push_back(value);
}

template <typename T, typename C>
void push_back(std::vector<std::pair<T, C>> &v, const T &value) {
    v.emplace_back(value, value);
}

template <typename T, typename C>
void remove_duplicates(std::vector<std::pair<T, C>> *vector) {
    auto first = vector->begin();
    auto last = vector->end();
    auto dest = first;

    while (++first != last) {
        if (first->first == dest->first) {
            dest->second += first->second;
        } else {
            *++dest = std::move(*first);
        }
    }
    vector->erase(++dest, vector->end());
}

TYPED_TEST(EliasFanoFileMergerTest, MergeEmpty) {
    constexpr uint32_t FILE_COUNT = 4;
    std::vector<utils::TempFile> files(FILE_COUNT);
    std::vector<std::string> file_names;
    for (uint32_t i = 0; i < FILE_COUNT; ++i) {
        file_names.push_back(files[i].name());
        do_encode(std::vector<TypeParam>(), file_names.back());
    }
    std::function<void(const TypeParam &v)> on_new_item
            = [](const TypeParam &) { FAIL() << "Should not be called."; };
    common::merge_files(file_names, on_new_item, false);
}

TYPED_TEST(EliasFanoFileMergerTest, MergeIdentical) {
    constexpr uint32_t FILE_COUNT = 4;
    std::vector<utils::TempFile> files(FILE_COUNT);
    std::vector<std::string> file_names;
    for (uint32_t i = 0; i < FILE_COUNT; ++i) {
        file_names.push_back(files[i].name());
        std::vector<TypeParam> values;
        for (uint32_t j = 0; j < 10; j = j + 1) {
            push_back(values, static_cast<utils::get_first_type_t<TypeParam>>(j));
        }
        do_encode(values, file_names.back());
    }
    std::function<void(const TypeParam &v)> on_new_item = [](const TypeParam &v) {
        static uint32_t i = 0;
        EXPECT_EQ(i, utils::get_first(v));
        i++;
    };
    common::merge_files(file_names, on_new_item);
}

TYPED_TEST(EliasFanoFileMergerTest, MergeRandom) {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(4, 10);
    std::uniform_int_distribution<std::mt19937::result_type> dist100(0, 100);

    const uint32_t file_count = dist10(rng);
    std::vector<utils::TempFile> files(file_count);
    std::vector<std::string> file_names;
    std::vector<TypeParam> expected;
    for (uint32_t i = 0; i < file_count; ++i) {
        file_names.push_back(files[i].name());
        std::vector<TypeParam> values
                = get_random_values<TypeParam>(dist100(rng), rng, dist10);
        do_encode(values, files[i].name());
        expected.insert(expected.end(), values.begin(), values.end());
    }
    std::sort(expected.begin(), expected.end(), [](const TypeParam &a, const TypeParam &b) {
        return utils::get_first(a) < utils::get_first(b);
    });
    if constexpr (utils::is_pair_v<TypeParam>) {
        remove_duplicates(&expected);
    } else {
        expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    }
    std::function<void(const TypeParam &v)> on_new_item = [expected](const TypeParam &v) {
        static uint32_t i = 0;
        EXPECT_LT(i, expected.size());
        EXPECT_EQ(expected[i], v) << i;
        i++;
    };
    common::merge_files(file_names, on_new_item);
}

} // namespace
