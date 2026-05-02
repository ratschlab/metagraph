#include <atomic>
#include <chrono>
#include <future>
#include <thread>

#include "gtest/gtest.h"

#include "cli/server_cache.hpp"


namespace {

using namespace mtg::cli;
using namespace std::chrono_literals;

// Helper to build a CachedResult with a given size, ignoring response content.
ServerQueryCache::CachedResult make_result(size_t bytes) {
    ServerQueryCache::CachedResult r;
    r.approx_size_bytes = bytes;
    return r;
}


TEST(ServerQueryCache, MissThenHit) {
    ServerQueryCache cache(1ull << 20);
    auto h1 = cache.acquire("k1");
    ASSERT_TRUE(h1.is_miss());
    h1.set_result(make_result(100));
    auto val1 = h1.get();
    ASSERT_NE(val1, nullptr);

    auto h2 = cache.acquire("k1");
    EXPECT_FALSE(h2.is_miss());
    auto val2 = h2.get();
    EXPECT_EQ(val1.get(), val2.get());  // shared CachedResult
}


TEST(ServerQueryCache, ConcurrentDeduplication) {
    ServerQueryCache cache(1ull << 20);

    // Producer thread acquires first, holds the handle, then publishes.
    std::atomic<int> compute_count = 0;
    std::promise<void> waiter_can_start;
    auto producer = std::thread([&]() {
        auto h = cache.acquire("dedup");
        ASSERT_TRUE(h.is_miss());
        compute_count.fetch_add(1);
        // Let the waiter run while we're still computing.
        waiter_can_start.set_value();
        std::this_thread::sleep_for(50ms);
        h.set_result(make_result(42));
        h.get();  // ensure the future is observed before destruction
    });

    waiter_can_start.get_future().wait();
    auto h = cache.acquire("dedup");
    // Second arrival must NOT be a miss — first computation owns the slot.
    EXPECT_FALSE(h.is_miss());
    auto val = h.get();  // blocks until producer publishes
    EXPECT_EQ(val->approx_size_bytes, 42u);
    producer.join();

    // Single computation despite two acquires.
    EXPECT_EQ(compute_count.load(), 1);
}


TEST(ServerQueryCache, EvictionUnderSizePressure) {
    ServerQueryCache cache(/* max */ 100);

    {
        auto h1 = cache.acquire("a");
        h1.set_result(make_result(60));
    }
    {
        auto h2 = cache.acquire("b");
        h2.set_result(make_result(60));
    }
    // Total 120 > 100 → LRU "a" should be evicted.
    EXPECT_FALSE(cache.contains("a"));
    EXPECT_TRUE(cache.contains("b"));
    EXPECT_LE(cache.size_bytes(), 100u);
}


TEST(ServerQueryCache, NeverEvictsEntryWithWaiters) {
    ServerQueryCache cache(/* max */ 100);

    // "a" has a live waiter — should never be evicted even under pressure.
    auto h_a = cache.acquire("a");
    h_a.set_result(make_result(60));

    // Insert "b" and let it become evictable (waiters drop to 0).
    {
        auto h_b = cache.acquire("b");
        h_b.set_result(make_result(60));
    }
    // Total at this point is 120 > 100, but "b" had waiters when its
    // result was published, so it survived the insert-time eviction pass.
    // Adding "c" triggers another eviction sweep: "a" is protected,
    // "b" is now waiterless, so "b" is the LRU victim.
    {
        auto h_c = cache.acquire("c");
        h_c.set_result(make_result(60));
    }
    EXPECT_TRUE(cache.contains("a"));
    EXPECT_FALSE(cache.contains("b"));
}


TEST(ServerQueryCache, FailedTtlExpires) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* failed_ttl */ 50ms);

    {
        auto h = cache.acquire("dies");
        h.set_result(make_result(10));
        h.mark_failed();
    }
    EXPECT_TRUE(cache.contains("dies"));
    std::this_thread::sleep_for(120ms);
    // Touch the cache to trigger opportunistic TTL eviction.
    {
        auto h = cache.acquire("touch");
        h.set_result(make_result(10));
    }
    EXPECT_FALSE(cache.contains("dies"));
}


TEST(ServerQueryCache, FailedWithinTtlIsKept) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* failed_ttl */ 1h);

    {
        auto h = cache.acquire("kept");
        h.set_result(make_result(10));
        h.mark_failed();
    }
    EXPECT_TRUE(cache.contains("kept"));
    auto h = cache.acquire("kept");
    EXPECT_FALSE(h.is_miss());
    EXPECT_EQ(h.get()->approx_size_bytes, 10u);
}


TEST(ServerQueryCache, DeliveredEntryStaysUntilSizePressure) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* failed_ttl */ 1ms);

    {
        auto h = cache.acquire("ok");
        h.set_result(make_result(10));
        h.mark_delivered();
    }
    std::this_thread::sleep_for(20ms);
    // Even though the FAILED-style TTL has elapsed, DELIVERED entries are
    // immune — only LRU pressure can evict them.
    auto h = cache.acquire("ok");
    EXPECT_FALSE(h.is_miss());
}


TEST(ServerQueryCache, DisabledCacheRecomputesEveryTime) {
    ServerQueryCache cache(/* max */ 0);

    auto h1 = cache.acquire("x");
    EXPECT_TRUE(h1.is_miss());
    h1.set_result(make_result(10));

    // Disabled cache: nothing is retained, second acquire is also a miss.
    auto h2 = cache.acquire("x");
    EXPECT_TRUE(h2.is_miss());
    EXPECT_EQ(cache.entry_count(), 0u);
    EXPECT_FALSE(cache.contains("x"));
}


TEST(ServerQueryCache, ProducerExceptionPropagatesToWaiters) {
    ServerQueryCache cache(1ull << 20);

    std::promise<void> waiter_ready;
    auto producer = std::thread([&]() {
        auto h = cache.acquire("boom");
        ASSERT_TRUE(h.is_miss());
        waiter_ready.set_value();
        std::this_thread::sleep_for(50ms);
        h.set_exception(std::make_exception_ptr(std::runtime_error("kaboom")));
    });

    waiter_ready.get_future().wait();
    auto h = cache.acquire("boom");
    EXPECT_FALSE(h.is_miss());
    EXPECT_THROW(h.get(), std::runtime_error);
    producer.join();
}


TEST(ServerQueryCache, AbandonedProducerUnblocksWaiters) {
    ServerQueryCache cache(1ull << 20);

    std::promise<void> waiter_ready;
    auto producer = std::thread([&]() {
        auto h = cache.acquire("orphan");
        ASSERT_TRUE(h.is_miss());
        waiter_ready.set_value();
        std::this_thread::sleep_for(20ms);
        // Drop the handle without ever calling set_result — destructor must
        // satisfy the promise with an exception so waiters don't hang.
    });

    waiter_ready.get_future().wait();
    auto h = cache.acquire("orphan");
    EXPECT_FALSE(h.is_miss());
    EXPECT_THROW(h.get(), std::runtime_error);
    producer.join();
}


} // namespace
