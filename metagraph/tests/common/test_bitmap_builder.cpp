#include "common/vectors/bitmap_builder.hpp"
#include "common/utils/file_utils.hpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <vector>


namespace {

size_t count_entries(const std::filesystem::path &dir) {
    return std::distance(std::filesystem::directory_iterator(dir),
                         std::filesystem::directory_iterator());
}

std::vector<uint64_t> get_ones(bitmap_builder &builder, uint64_t size) {
    auto data = builder.get_initialization_data();
    EXPECT_EQ(size, data.size);
    std::vector<uint64_t> ones;
    data.call_ones([&](uint64_t i) { ones.push_back(i); });
    EXPECT_EQ(ones.size(), data.num_set_bits);
    return ones;
}

class BitmapBuilderSetDisk : public ::testing::Test {
  protected:
    void SetUp() override {
        swap_dir_ = utils::create_temp_dir(std::filesystem::temp_directory_path(),
                                           "test_bitmap_builder");
    }
    void TearDown() override { utils::remove_temp_dir(swap_dir_); }

    std::filesystem::path swap_dir_;
};

TEST_F(BitmapBuilderSetDisk, Empty) {
    bitmap_builder_set_disk builder(100, 1, 8, swap_dir_);
    EXPECT_EQ(std::vector<uint64_t>(), get_ones(builder, 100));
    EXPECT_EQ(0u, count_entries(swap_dir_));
}

TEST_F(BitmapBuilderSetDisk, FitsInBufferWithoutDisk) {
    for (size_t num_threads : { 1, 4 }) {
        bitmap_builder_set_disk builder(1000, num_threads, 100, swap_dir_);
        // unsorted and with duplicates
        std::vector<uint64_t> positions { 999, 7, 42, 7, 0, 42, 500 };
        builder.add_ones(positions.data(), positions.data() + positions.size());
        builder.add_one(3);
        builder.add_one(999);

        EXPECT_EQ(0u, count_entries(swap_dir_));
        EXPECT_EQ(std::vector<uint64_t>({ 0, 3, 7, 42, 500, 999 }),
                  get_ones(builder, 1000));
        EXPECT_EQ(0u, count_entries(swap_dir_));
    }
}

TEST_F(BitmapBuilderSetDisk, NearlyFullBufferStaysInMemory) {
    bitmap_builder_set_disk builder(1000, 1, 100, swap_dir_);
    std::vector<uint64_t> expected;
    // 90 distinct positions inserted one by one into a buffer of 100
    for (uint64_t i = 0; i < 90; ++i) {
        builder.add_one(i * 3);
        expected.push_back(i * 3);
    }
    EXPECT_EQ(0u, count_entries(swap_dir_));
    EXPECT_EQ(expected, get_ones(builder, 1000));
}

TEST_F(BitmapBuilderSetDisk, DuplicatesDoNotForceDisk) {
    bitmap_builder_set_disk builder(1000, 1, 100, swap_dir_);
    // many more insertions than the buffer holds, but only 10 distinct positions
    for (size_t i = 0; i < 5000; ++i) {
        builder.add_one(i % 10);
    }
    EXPECT_EQ(0u, count_entries(swap_dir_));
    EXPECT_EQ(std::vector<uint64_t>({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }),
              get_ones(builder, 1000));
}

TEST_F(BitmapBuilderSetDisk, OverflowGoesToDisk) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t batch_size : { 1, 3, 1000 }) {
            const uint64_t size = 100'000;
            std::vector<uint64_t> expected;
            {
                bitmap_builder_set_disk builder(size, num_threads, 64, swap_dir_);
                std::vector<uint64_t> batch;
                // descending order, every position inserted twice
                for (uint64_t i = 3000; i-- > 0;) {
                    batch.push_back(i * 7);
                    batch.push_back(i * 7);
                    if (batch.size() >= batch_size) {
                        builder.add_ones(batch.data(), batch.data() + batch.size());
                        batch.clear();
                    }
                }
                builder.add_ones(batch.data(), batch.data() + batch.size());
                EXPECT_EQ(1u, count_entries(swap_dir_));

                for (uint64_t i = 0; i < 3000; ++i) {
                    expected.push_back(i * 7);
                }
                EXPECT_EQ(expected, get_ones(builder, size));
            }
            // the temporary directory is removed with the builder
            EXPECT_EQ(0u, count_entries(swap_dir_));
        }
    }
}

} // namespace
