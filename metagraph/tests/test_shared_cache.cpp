#include "gtest/gtest.h"

#include <array>
#include <cstdint>
#include <thread>
#include <vector>

#include "common/shared_cache.hpp"


namespace {

using mtg::common::SharedCache;

TEST(SharedCache, Disabled) {
    SharedCache<uint64_t, uint64_t> cache(0);
    EXPECT_EQ(0u, cache.capacity());
    cache.insert(1, 2);
    EXPECT_FALSE(cache.find(1));
}

TEST(SharedCache, Capacity) {
    EXPECT_EQ(2u, (SharedCache<uint64_t, uint64_t>(1).capacity()));
    EXPECT_EQ(8u, (SharedCache<uint64_t, uint64_t>(5).capacity()));
    EXPECT_EQ(131072u, (SharedCache<uint64_t, uint64_t>(100'000).capacity()));
}

TEST(SharedCache, FindInserted) {
    SharedCache<uint64_t, std::array<uint64_t, 2>> cache(1024);
    EXPECT_FALSE(cache.find(0));
    EXPECT_FALSE(cache.find(5));
    for (uint64_t i = 0; i < 100; ++i) {
        cache.insert(i, { i + 1, i + 2 });
    }
    size_t found = 0;
    for (uint64_t i = 0; i < 100; ++i) {
        if (auto value = cache.find(i)) {
            EXPECT_EQ((std::array<uint64_t, 2>{ i + 1, i + 2 }), *value);
            ++found;
        }
    }
    // a few keys may have replaced each other in a shared slot
    EXPECT_LE(90u, found);
    EXPECT_FALSE(cache.find(100));
}

TEST(SharedCache, MappedTableFindsInserted) {
    // 131,072 slots of 24 B: a table mapped from the kernel (at least 1 MiB), empty until
    // a slot is written
    SharedCache<uint64_t, uint64_t> cache(100'000);
    ASSERT_EQ(131072u, cache.capacity());
    EXPECT_FALSE(cache.find(0));
    EXPECT_FALSE(cache.find(12345));
    for (uint64_t i = 0; i < 1000; ++i) {
        cache.insert(i * 977, i);
    }
    size_t found = 0;
    for (uint64_t i = 0; i < 1000; ++i) {
        if (auto value = cache.find(i * 977)) {
            EXPECT_EQ(i, *value);
            ++found;
        }
    }
    // a few keys may have replaced each other in a shared slot
    EXPECT_LE(990u, found);
    EXPECT_FALSE(cache.find(1000 * 977));
}

TEST(SharedCache, Overwrite) {
    SharedCache<uint64_t, bool> cache(16);
    cache.insert(7, true);
    ASSERT_TRUE(cache.find(7));
    EXPECT_TRUE(*cache.find(7));
    cache.insert(7, false);
    ASSERT_TRUE(cache.find(7));
    EXPECT_FALSE(*cache.find(7));
}

TEST(SharedCache, ReplacesOnCollision) {
    // with two slots, three keys cannot all stay
    SharedCache<uint64_t, uint64_t> cache(2);
    for (uint64_t i = 1; i <= 3; ++i) {
        cache.insert(i, i * 10);
    }
    size_t found = 0;
    for (uint64_t i = 1; i <= 3; ++i) {
        if (auto value = cache.find(i)) {
            EXPECT_EQ(i * 10, *value);
            ++found;
        }
    }
    EXPECT_GE(2u, found);
    EXPECT_LE(1u, found);
}

TEST(SharedCache, ConcurrentReadsSeeWholeEntries) {
    // a small cache, so that threads keep writing the same slots
    SharedCache<uint64_t, std::array<uint64_t, 2>> cache(64);
    std::vector<std::thread> threads;
    for (uint64_t t = 0; t < 4; ++t) {
        threads.emplace_back([&cache, t]() {
            for (uint64_t i = 0; i < 200'000; ++i) {
                uint64_t key = (i * 7 + t) % 1000;
                if (auto value = cache.find(key)) {
                    ASSERT_EQ(key * 3, (*value)[0]);
                    ASSERT_EQ(~key, (*value)[1]);
                } else {
                    cache.insert(key, { key * 3, ~key });
                }
            }
        });
    }
    for (auto &thread : threads) {
        thread.join();
    }
}

} // namespace
