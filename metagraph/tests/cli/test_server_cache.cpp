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


TEST(ServerQueryCache, ProtectedWithinWindowPreservedOverDelivered) {
    // Pass 1 should evict the DELIVERED entry before touching the
    // PROTECTED entry that's still inside its retention window.
    ServerQueryCache cache(/* max */ 100, /* protection_ttl */ 1h);

    {
        auto h_failed = cache.acquire("failed");
        h_failed.set_result(make_result(60));
        h_failed.mark_protected();
    }
    {
        auto h_ok = cache.acquire("ok");
        h_ok.set_result(make_result(40));
        h_ok.mark_delivered();
    }
    // Both waiterless, total=100 (== max). Force a touch under pressure.
    {
        auto h_extra = cache.acquire("extra");
        h_extra.set_result(make_result(40));  // 100 + 40 = 140 > 100
        h_extra.mark_delivered();
    }
    EXPECT_TRUE(cache.contains("failed"));   // protected
    EXPECT_FALSE(cache.contains("ok"));      // sacrificed first
    EXPECT_TRUE(cache.contains("extra"));
}


TEST(ServerQueryCache, ProtectedEntriesEvictedUnderHeavyPressure) {
    // When all waiterless entries are PROTECTED-within-window and the cache is
    // overfull, pass 2 must still evict — otherwise the cache grows without
    // bound.
    ServerQueryCache cache(/* max */ 100, /* protection_ttl */ 1h);

    auto store_failed = [&](const std::string &key, size_t bytes) {
        auto h = cache.acquire(key);
        h.set_result(make_result(bytes));
        h.mark_protected();
    };

    store_failed("f1", 60);
    store_failed("f2", 60);  // 60 + 60 = 120 > 100
    // f1 was protected during f2's set_result (f2 had waiters); f1 is now
    // the only waiterless entry, and pass 2 lets us evict it despite TTL.
    EXPECT_FALSE(cache.contains("f1"));
    EXPECT_TRUE(cache.contains("f2"));
    EXPECT_LE(cache.size_bytes(), 100u);
}


TEST(ServerQueryCache, ProtectedPastWindowGraduatesToMainCache) {
    // After `protection_ttl` of no retry hits, a PROTECTED entry stays
    // in the cache but loses its priority — it's evictable like a
    // DELIVERED entry, only under size pressure (no proactive removal).
    ServerQueryCache cache(/* max */ 100, /* protection_ttl */ 50ms);

    {
        auto h = cache.acquire("aged");
        h.set_result(make_result(60));
        h.mark_protected();
    }
    std::this_thread::sleep_for(120ms);
    // Past the window — entry stays in cache.
    EXPECT_TRUE(cache.contains("aged"));

    // Insert a DELIVERED entry that triggers size pressure. The aged
    // PROTECTED entry is no longer prioritized → evictable in pass 1
    // at normal LRU priority.
    {
        auto h = cache.acquire("new");
        h.set_result(make_result(60));   // 60 + 60 = 120 > 100
        h.mark_delivered();
    }
    EXPECT_FALSE(cache.contains("aged"));
    EXPECT_TRUE(cache.contains("new"));
}


TEST(ServerQueryCache, ProtectedHitRefreshesPriorityWindow) {
    // A cache hit on a PROTECTED entry refreshes its retention TTL —
    // the window is sliding, tracking the upstream's retry pattern.
    ServerQueryCache cache(/* max */ 100, /* protection_ttl */ 100ms);

    {
        auto h = cache.acquire("retried");
        h.set_result(make_result(60));
        h.mark_protected();
    }
    // Sacrificial DELIVERED entry: filling cache to capacity so a later
    // insert triggers pressure with a waiterless eviction candidate
    // available (otherwise pass 1 would have nothing to evict and pass 2
    // would falsely sacrifice the PROTECTED entry we're trying to test).
    {
        auto h = cache.acquire("filler");
        h.set_result(make_result(40));   // total = 100, at cap
        h.mark_delivered();
    }

    std::this_thread::sleep_for(60ms);
    // Mid-window retry: cache hit refreshes ready_at on "retried".
    {
        auto h = cache.acquire("retried");
        ASSERT_FALSE(h.is_miss());
        h.get();
    }
    std::this_thread::sleep_for(60ms);
    // Elapsed = 120ms (past the original 100ms TTL). With the sliding
    // window refresh we're only 60ms into the new window — still protected.

    // Insert an intruder that triggers pressure. Pass 1 should sacrifice
    // "filler" (DELIVERED, waiterless) and leave "retried" alone.
    {
        auto h = cache.acquire("intruder");
        h.set_result(make_result(40));   // total during insert: 140 > 100
        h.mark_delivered();
    }
    EXPECT_TRUE(cache.contains("retried"));    // protected by refreshed window
    EXPECT_FALSE(cache.contains("filler"));    // sacrificed in pass 1
    EXPECT_TRUE(cache.contains("intruder"));
}


TEST(ServerQueryCache, ProtectedWithinWindowIsKept) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* protection_ttl */ 1h);

    {
        auto h = cache.acquire("kept");
        h.set_result(make_result(10));
        h.mark_protected();
    }
    EXPECT_TRUE(cache.contains("kept"));
    auto h = cache.acquire("kept");
    EXPECT_FALSE(h.is_miss());
    EXPECT_EQ(h.get()->approx_size_bytes, 10u);
}


TEST(ServerQueryCache, DeliveredIsSinkStateAgainstLaterFailure) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* protection_ttl */ 50ms);

    // First delivery succeeds — entry becomes DELIVERED.
    {
        auto h = cache.acquire("k");
        h.set_result(make_result(10));
        h.mark_delivered();
    }
    // Second acquirer of the same entry has a delivery failure.
    {
        auto h = cache.acquire("k");
        ASSERT_FALSE(h.is_miss());
        h.mark_protected();  // must NOT downgrade DELIVERED → PROTECTED
    }
    // Wait past the (very short) failed TTL. If the sink semantics held,
    // the entry stays. If they didn't, evict_expired_failed_locked would
    // drop it.
    std::this_thread::sleep_for(120ms);
    {
        // Trigger an eviction sweep via release_waiter.
        auto h = cache.acquire("touch");
        h.set_result(make_result(10));
    }
    EXPECT_TRUE(cache.contains("k"));
}


TEST(ServerQueryCache, DeliveredEntryStaysUntilSizePressure) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* protection_ttl */ 1ms);

    {
        auto h = cache.acquire("ok");
        h.set_result(make_result(10));
        h.mark_delivered();
    }
    std::this_thread::sleep_for(20ms);
    // Even though the protection TTL has elapsed, DELIVERED entries are
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
