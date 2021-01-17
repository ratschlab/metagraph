#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "common/file_merger.hpp"
#include "common/utils/file_utils.hpp"


namespace {

using namespace mtg;

template <typename T>
class FileMergerTest : public ::testing::Test {};

typedef ::testing::Types<uint32_t, uint64_t> ValueTypes;

TYPED_TEST_SUITE(FileMergerTest, ValueTypes);


template <typename T>
std::vector<T>
get_random_values(uint32_t count,
                  std::mt19937 &rng,
                  std::uniform_int_distribution<std::mt19937::result_type> &dist) {
    std::vector<T> values(count);
    for (uint32_t j = 0; j < values.size(); ++j) {
        values[j] = j > 0 ? values[j - 1] + (T)dist(rng) : (T)dist(rng);
    }
    return values;
}


TYPED_TEST(FileMergerTest, MergeEmpty) {
    constexpr uint32_t FILE_COUNT = 4;
    std::vector<utils::TempFile> files(FILE_COUNT);
    std::vector<std::string> file_names;
    for (uint32_t i = 0; i < FILE_COUNT; ++i) {
        file_names.push_back(files[i].name());
    }
    std::function<void(const TypeParam &v)> on_new_item
            = [](const TypeParam &) { FAIL() << "Should not be called."; };
    common::merge_files(file_names, on_new_item);
}

TYPED_TEST(FileMergerTest, MergeIdentical) {
    constexpr uint32_t FILE_COUNT = 4;
    std::vector<utils::TempFile> files(FILE_COUNT);
    std::vector<std::string> file_names;
    for (uint32_t i = 0; i < FILE_COUNT; ++i) {
        file_names.push_back(files[i].name());
        std::vector<TypeParam> values;
        for (uint32_t j = 0; j < 10; j = j + 1) {
            values.push_back(j);
        }
        std::ofstream f(files[i].name(), std::ios::binary);
        f.write(reinterpret_cast<char *>(values.data()), values.size() * sizeof(values[0]));
    }
    std::function<void(const TypeParam &v)> on_new_item = [](const TypeParam &v) {
        static uint32_t i = 0;
        ASSERT_EQ(i, v);
        i++;
    };
    common::merge_files(file_names, on_new_item);
}

TYPED_TEST(FileMergerTest, MergeRandom) {
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
        std::ofstream f(files[i].name(), std::ios::binary);
        f.write(reinterpret_cast<char *>(values.data()), values.size() * sizeof(values[0]));
        expected.insert(expected.end(), values.begin(), values.end());
    }
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());

    std::function<void(const TypeParam &v)> on_new_item = [expected](const TypeParam &v) {
        static uint32_t i = 0;
        ASSERT_LT(i, expected.size());
        ASSERT_EQ(expected[i], v) << "Unexpected value for element " << i;
        i++;
    };
    common::merge_files(file_names, on_new_item);
}

} // namespace
