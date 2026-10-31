#include <atomic>
#include <chrono>
#include <cmath>
#include <future>
#include <locale>
#include <thread>

#include "gtest/gtest.h"

#include "cli/server_cache.hpp"
#include "cli/server_request.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

// Set up expired windows without sleeps or a clock hook in the production path.
class ServerQueryCacheTestAccess {
  public:
    static void expire(ServerQueryCache &cache, const std::string &key) {
        std::lock_guard<std::mutex> lock(cache.mutex_);
        cache.map_.at(key)->ready_at = std::chrono::steady_clock::now()
                                      - cache.protection_ttl_ - std::chrono::hours(1);
    }
};

} // namespace cli
} // namespace mtg


namespace {

using namespace mtg::cli;
using namespace std::chrono_literals;

// Build a placeholder response string of the requested size in KiB. The
// cache value is the serialized JSON body — for the eviction-policy
// tests we only care about its size, so the contents are arbitrary.
std::string response_of_kib(size_t kib) {
    return std::string(kib * 1024, '\0');
}


TEST(ServerQueryCache, MissThenHit) {
    ServerQueryCache cache(1ull << 20);
    auto h1 = cache.acquire("k1");
    ASSERT_TRUE(h1.is_miss());
    h1.set_result(response_of_kib(100));
    auto val1 = h1.get();
    ASSERT_NE(val1, nullptr);

    auto h2 = cache.acquire("k1");
    EXPECT_FALSE(h2.is_miss());
    auto val2 = h2.get();
    EXPECT_EQ(val1.get(), val2.get());  // both handles share the same response object
}


TEST(ServerQueryCache, ConcurrentDeduplication) {
    ServerQueryCache cache(1ull << 20);

    // Producer thread acquires first, holds the handle, then publishes.
    std::atomic<int> compute_count = 0;
    std::promise<void> waiter_can_start;
    std::promise<void> waiter_attached;
    auto producer = std::thread([&]() {
        auto h = cache.acquire("dedup");
        EXPECT_TRUE(h.is_miss());
        compute_count.fetch_add(1);
        // Let the waiter run while we're still computing.
        waiter_can_start.set_value();
        waiter_attached.get_future().wait();
        h.set_result(response_of_kib(42));
        h.get();  // ensure the future is observed before destruction
    });

    waiter_can_start.get_future().wait();
    auto h = cache.acquire("dedup");
    // Second arrival must NOT be a miss — first computation owns the slot.
    EXPECT_FALSE(h.is_miss());
    waiter_attached.set_value();
    auto val = h.get();  // blocks until producer publishes
    EXPECT_EQ(val->size(), 42 * 1024u);
    producer.join();

    // Single computation despite two acquires.
    EXPECT_EQ(compute_count.load(), 1);
}


TEST(ServerQueryCache, EvictionUnderSizePressure) {
    ServerQueryCache cache(/* max */ 104 * 1024);

    {
        auto h1 = cache.acquire("a");
        h1.set_result(response_of_kib(60));
    }
    {
        auto h2 = cache.acquire("b");
        h2.set_result(response_of_kib(60));
    }
    // The budget fits one 60 KiB result plus metadata, but not two.
    EXPECT_FALSE(cache.contains("a"));
    EXPECT_TRUE(cache.contains("b"));
    EXPECT_LE(cache.size_bytes(), 104 * 1024u);
}


TEST(ServerQueryCache, NeverEvictsEntryWithWaiters) {
    ServerQueryCache cache(/* max */ 104 * 1024);

    // "a" has a live waiter — should never be evicted even under pressure.
    auto h_a = cache.acquire("a");
    h_a.set_result(response_of_kib(60));

    // Insert "b" and let it become evictable (waiters drop to 0).
    {
        auto h_b = cache.acquire("b");
        h_b.set_result(response_of_kib(60));
    }
    // The held entry survives pressure; waiterless entries are evicted.
    {
        auto h_c = cache.acquire("c");
        h_c.set_result(response_of_kib(60));
    }
    EXPECT_TRUE(cache.contains("a"));
    EXPECT_FALSE(cache.contains("b"));
}


TEST(ServerQueryCache, ProtectedWithinWindowPreservedOverDelivered) {
    // Pass 1 should evict the DELIVERED entry before touching the
    // PROTECTED entry that's still inside its retention window.
    ServerQueryCache cache(/* max */ 104 * 1024, /* protection_ttl */ 1h);

    {
        auto h_failed = cache.acquire("failed");
        h_failed.set_result(response_of_kib(60));
        h_failed.mark_protected();
    }
    {
        auto h_ok = cache.acquire("ok");
        h_ok.set_result(response_of_kib(40));
        h_ok.mark_delivered();
    }
    // Both fit with metadata. A third result forces eviction.
    {
        auto h_extra = cache.acquire("extra");
        h_extra.set_result(response_of_kib(40));  // exceeds the budget with three results
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
    ServerQueryCache cache(/* max */ 104 * 1024, /* protection_ttl */ 1h);

    auto store_failed = [&](const std::string &key, size_t bytes) {
        auto h = cache.acquire(key);
        h.set_result(response_of_kib(bytes));
        h.mark_protected();
    };

    store_failed("f1", 60);
    store_failed("f2", 60);  // two results exceed the budget
    // f1 was protected during f2's set_result (f2 had waiters); f1 is now
    // the only waiterless entry, and pass 2 lets us evict it despite TTL.
    EXPECT_FALSE(cache.contains("f1"));
    EXPECT_TRUE(cache.contains("f2"));
    EXPECT_LE(cache.size_bytes(), 104 * 1024u);
}


TEST(ServerQueryCache, ProtectedPastWindowGraduatesToMainCache) {
    ServerQueryCache cache(104 * 1024, 1h);
    {
        auto h = cache.acquire("aged");
        h.set_result(response_of_kib(40));
        h.mark_protected();
    }
    {
        auto h = cache.acquire("delivered");
        h.set_result(response_of_kib(40));
        h.mark_delivered();
    }
    ServerQueryCacheTestAccess::expire(cache, "aged");
    EXPECT_TRUE(cache.contains("aged"));
    {
        auto h = cache.acquire("new");
        h.set_result(response_of_kib(40));
    }
    // Without expiration, the newer delivered entry would be evicted instead.
    EXPECT_FALSE(cache.contains("aged"));
    EXPECT_TRUE(cache.contains("delivered"));
    EXPECT_TRUE(cache.contains("new"));
}

TEST(ServerQueryCache, ProtectedHitRefreshesPriorityWindow) {
    ServerQueryCache cache(104 * 1024, 1h);
    {
        auto h = cache.acquire("retried");
        h.set_result(response_of_kib(40));
        h.mark_protected();
    }
    {
        auto h = cache.acquire("filler");
        h.set_result(response_of_kib(40));
        h.mark_delivered();
    }
    ServerQueryCacheTestAccess::expire(cache, "retried");
    {
        auto h = cache.acquire("retried");
        ASSERT_FALSE(h.is_miss());
        h.get();
    }
    // Make retried the LRU entry again, so recency alone cannot save it.
    { auto h = cache.acquire("filler"); }
    {
        auto h = cache.acquire("intruder");
        h.set_result(response_of_kib(40));
    }
    EXPECT_TRUE(cache.contains("retried"));
    EXPECT_FALSE(cache.contains("filler"));
    EXPECT_TRUE(cache.contains("intruder"));
}


TEST(ServerQueryCache, ProtectedWithinWindowIsKept) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* protection_ttl */ 1h);

    {
        auto h = cache.acquire("kept");
        h.set_result(response_of_kib(10));
        h.mark_protected();
    }
    EXPECT_TRUE(cache.contains("kept"));
    auto h = cache.acquire("kept");
    EXPECT_FALSE(h.is_miss());
    EXPECT_EQ(h.get()->size(), 10 * 1024u);
}


TEST(ServerQueryCache, DeliveredIsSinkStateAgainstLaterFailure) {
    ServerQueryCache cache(/* max */ 104 * 1024, /* protection_ttl */ 1h);
    {
        auto h = cache.acquire("protected");
        h.set_result(response_of_kib(40));
        h.mark_protected();
    }
    {
        auto h = cache.acquire("delivered");
        h.set_result(response_of_kib(40));
        h.mark_delivered();
        h.mark_protected();  // must not downgrade DELIVERED
    }
    {
        auto h = cache.acquire("other");
        h.set_result(response_of_kib(40));
        h.mark_delivered();
    }
    // A downgrade would make pass 2 evict the older protected entry instead.
    EXPECT_TRUE(cache.contains("protected"));
    EXPECT_FALSE(cache.contains("delivered"));
    EXPECT_TRUE(cache.contains("other"));
}

TEST(ServerQueryCache, EvictsMultipleEntriesPastPinnedAndProtectedEntries) {
    ServerQueryCache cache(/* max */ 104 * 1024, /* protection_ttl */ 1h);
    auto pinned = cache.acquire("pinned");
    pinned.set_result(response_of_kib(16));
    {
        auto h = cache.acquire("protected");
        h.set_result(response_of_kib(20));
        h.mark_protected();
    }
    for (const char *key : {"a", "b", "c"}) {
        auto h = cache.acquire(key);
        h.set_result(response_of_kib(20));
        h.mark_delivered();
    }
    {
        auto h = cache.acquire("large");
        h.set_result(response_of_kib(60));
    }
    EXPECT_TRUE(cache.contains("pinned"));
    EXPECT_TRUE(cache.contains("protected"));
    EXPECT_TRUE(cache.contains("large"));
    for (const char *key : {"a", "b", "c"})
        EXPECT_FALSE(cache.contains(key));
    EXPECT_LE(cache.size_bytes(), 104 * 1024u);
}


TEST(ServerQueryCache, OversizedResponsePreservesCachedResultsAndWaiters) {
    ServerQueryCache cache(64 * 1024);
    {
        auto h = cache.acquire("small");
        ASSERT_TRUE(h.set_result(response_of_kib(20)));
    }
    auto producer = cache.acquire("large");
    auto waiter = cache.acquire("large");
    ASSERT_FALSE(waiter.is_miss());
    EXPECT_FALSE(producer.set_result(response_of_kib(128)));
    EXPECT_EQ(producer.get(), waiter.get());
    EXPECT_EQ(waiter.get()->size(), 128 * 1024u);
    EXPECT_FALSE(cache.contains("large"));
    EXPECT_TRUE(cache.contains("small"));
    EXPECT_LE(cache.size_bytes(), 64 * 1024u);
    producer = ServerQueryCache::Handle();
    waiter = ServerQueryCache::Handle();
    EXPECT_TRUE(cache.contains("small"));
}


TEST(ServerQueryCache, DeliveredEntryStaysUntilSizePressure) {
    ServerQueryCache cache(/* max */ 1ull << 20, /* protection_ttl */ 1ms);

    {
        auto h = cache.acquire("ok");
        h.set_result(response_of_kib(10));
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
    h1.set_result(response_of_kib(10));

    // Disabled cache: nothing is retained, second acquire is also a miss.
    auto h2 = cache.acquire("x");
    EXPECT_TRUE(h2.is_miss());
    EXPECT_EQ(cache.entry_count(), 0u);
    EXPECT_FALSE(cache.contains("x"));
}


TEST(ServerQueryCache, ProducerExceptionPropagatesToWaiters) {
    ServerQueryCache cache(1ull << 20);

    std::promise<void> waiter_ready;
    std::promise<void> waiter_attached;
    auto producer = std::thread([&]() {
        auto h = cache.acquire("boom");
        EXPECT_TRUE(h.is_miss());
        waiter_ready.set_value();
        waiter_attached.get_future().wait();
        h.set_exception(std::make_exception_ptr(std::runtime_error("kaboom")));
    });

    waiter_ready.get_future().wait();
    auto h = cache.acquire("boom");
    EXPECT_FALSE(h.is_miss());
    waiter_attached.set_value();
    EXPECT_THROW(h.get(), std::runtime_error);
    producer.join();
}


TEST(ServerQueryCache, AbandonedProducerUnblocksWaiters) {
    ServerQueryCache cache(1ull << 20);

    std::promise<void> waiter_ready;
    std::promise<void> waiter_attached;
    auto producer = std::thread([&]() {
        auto h = cache.acquire("orphan");
        EXPECT_TRUE(h.is_miss());
        waiter_ready.set_value();
        waiter_attached.get_future().wait();
        // Drop the handle without ever calling set_result — destructor must
        // satisfy the promise with an exception so waiters don't hang.
    });

    waiter_ready.get_future().wait();
    auto h = cache.acquire("orphan");
    EXPECT_FALSE(h.is_miss());
    waiter_attached.set_value();
    EXPECT_THROW(h.get(), std::runtime_error);
    producer.join();
}

TEST(ServerQueryCache, FailedComputationCanBeRetriedWhileWaitersRemain) {
    ServerQueryCache cache(1024);
    auto producer = cache.acquire("retry");
    auto waiter = cache.acquire("retry");
    producer.set_exception(std::make_exception_ptr(std::runtime_error("temporary")));

    EXPECT_THROW(producer.get(), std::runtime_error);
    EXPECT_THROW(waiter.get(), std::runtime_error);
    EXPECT_FALSE(cache.contains("retry"));

    auto retry = cache.acquire("retry");
    ASSERT_TRUE(retry.is_miss());
    retry.set_result("success");
    // Late callbacks and destruction from the failed attempt must not affect
    // the new entry with the same key.
    waiter.mark_protected();
    producer.mark_delivered();
    waiter = ServerQueryCache::Handle();
    producer = ServerQueryCache::Handle();
    auto hit = cache.acquire("retry");
    EXPECT_FALSE(hit.is_miss());
    EXPECT_EQ(*hit.get(), "success");
    EXPECT_EQ(cache.entry_count(), 1u);
    EXPECT_GT(cache.size_bytes(), 7u);
}

TEST(ServerQueryCache, AbandonedProducerCanBeRetried) {
    ServerQueryCache cache(1024);
    ServerQueryCache::Handle waiter;
    {
        auto producer = cache.acquire("orphan");
        waiter = cache.acquire("orphan");
    }
    EXPECT_THROW(waiter.get(), std::runtime_error);
    auto retry = cache.acquire("orphan");
    ASSERT_TRUE(retry.is_miss());
    retry.set_result("success");
    EXPECT_EQ(*retry.get(), "success");
}

TEST(ServerQueryCache, MoveAssignmentAbandonsPreviousComputation) {
    ServerQueryCache cache(1024);
    auto producer = cache.acquire("old");
    auto waiter = cache.acquire("old");
    auto replacement = cache.acquire("new");
    replacement.set_result("new result");

    producer = std::move(replacement);
    EXPECT_THROW(waiter.get(), std::runtime_error);
    EXPECT_EQ(*producer.get(), "new result");
    auto retry = cache.acquire("old");
    ASSERT_TRUE(retry.is_miss());
    retry.set_result("retried");
    EXPECT_EQ(cache.entry_count(), 2u);
}

TEST(ServerQueryCache, FailedComputationsDoNotAccumulate) {
    for (size_t budget : {0u, 1024u}) {
        ServerQueryCache cache(budget);
        for (size_t i = 0; i < 100; ++i) {
            auto h = cache.acquire(std::to_string(i));
            h.set_exception(std::make_exception_ptr(std::runtime_error("temporary")));
            EXPECT_EQ(cache.entry_count(), 0u);
            EXPECT_EQ(cache.size_bytes(), 0u);
        }
    }
}

TEST(ServerQueryCache, CacheKeyPreservesFloatingPointParameters) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    Json::Value request;
    request["FASTA"] = ">query\nACGTACGT\n";
    for (const char *parameter : {"discovery_fraction", "min_exact_match",
                                  "max_num_nodes_per_seq_char"}) {
        SCOPED_TRACE(parameter);
        auto first = request, second = request;
        first[parameter] = 0.12345641;
        second[parameter] = std::nextafter(0.12345641, 1.0);
        EXPECT_NE(make_search_cache_key(parse_search_request(first, config), "graph"),
                  make_search_cache_key(parse_search_request(second, config), "graph"));
    }
}

TEST(ServerQueryCache, CacheKeyDoesNotDependOnLocale) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    Json::Value request;
    request["FASTA"] = ">query\nACGTACGT\n";
    request["discovery_fraction"] = 0.25;
    auto key = make_search_cache_key(parse_search_request(request, config), "graph");
    struct CommaDecimal : std::numpunct<char> {
        char do_decimal_point() const override { return ','; }
    };
    struct RestoreLocale {
        std::locale previous = std::locale();
        ~RestoreLocale() { std::locale::global(previous); }
    } restore;
    std::locale::global(std::locale(std::locale::classic(), new CommaDecimal));
    EXPECT_EQ(make_search_cache_key(parse_search_request(request, config), "graph"), key);
}

TEST(ServerSearchRequest, ResolvesDefaultsAndOverridesWithoutChangingServerConfig) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    config.discovery_fraction = 0.25;
    config.alignment_min_exact_match = 0.5;
    config.alignment_max_nodes_per_seq_char = 3.0;
    config.num_top_labels = 7;
    config.query_mode = COORDS;  // HTTP mode defaults to MATCHES, independently of CLI mode.
    config.tmp_dir = "/tmp/request-test";
    Json::Value json;
    json["FASTA"] = ">q\nACGT\n";
    auto defaults = parse_search_request(json, config);
    EXPECT_EQ(defaults.fasta, json["FASTA"].asString());
    EXPECT_DOUBLE_EQ(defaults.config.discovery_fraction, 0.25);
    EXPECT_DOUBLE_EQ(defaults.config.alignment_min_exact_match, 0.5);
    EXPECT_DOUBLE_EQ(defaults.config.alignment_max_nodes_per_seq_char, 3.0);
    EXPECT_EQ(defaults.config.num_top_labels, 7u);
    EXPECT_EQ(defaults.config.query_mode, MATCHES);
    EXPECT_EQ(defaults.config.tmp_dir, config.tmp_dir);
    EXPECT_FALSE(defaults.align);

    json["discovery_fraction"] = 0.75;
    json["min_exact_match"] = 0.125;
    json["max_num_nodes_per_seq_char"] = 9.0;
    json["top_labels"] = 2;
    json["align"] = true;
    auto overrides = parse_search_request(json, config);
    EXPECT_DOUBLE_EQ(overrides.config.discovery_fraction, 0.75);
    EXPECT_DOUBLE_EQ(overrides.config.alignment_min_exact_match, 0.125);
    EXPECT_DOUBLE_EQ(overrides.config.alignment_max_nodes_per_seq_char, 9.0);
    EXPECT_EQ(overrides.config.num_top_labels, 2u);
    EXPECT_TRUE(overrides.align);
    EXPECT_DOUBLE_EQ(config.discovery_fraction, 0.25);
    EXPECT_EQ(config.num_top_labels, 7u);
    EXPECT_EQ(config.query_mode, COORDS);
}

TEST(ServerSearchRequest, ModePrecedenceAndIgnoredFlagsShareCacheIdentity) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    Json::Value json;
    json["FASTA"] = ">q\nACGT\n";
    const char *flags[] = {"query_coords", "query_counts", "with_signature", "abundance_sum"};
    const QueryMode modes[] = {COORDS, COUNTS, SIGNATURE, COUNTS_SUM};
    for (size_t i = 0; i < 4; ++i) {
        auto one = json, several = json;
        one[flags[i]] = true;
        for (size_t j = i; j < 4; ++j)
            several[flags[j]] = true;
        auto request = parse_search_request(several, config);
        EXPECT_EQ(request.config.query_mode, modes[i]);
        EXPECT_EQ(make_search_cache_key(request, "graph"),
                  make_search_cache_key(parse_search_request(one, config), "graph"));
        EXPECT_NE(make_search_cache_key(request, "graph"),
                  make_search_cache_key(parse_search_request(json, config), "graph"));
    }
}

TEST(ServerSearchRequest, EffectiveDefaultsAndTransportFieldsShareCacheIdentity) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    Json::Value json;
    json["FASTA"] = ">q\nACGT\n";
    const auto key = make_search_cache_key(parse_search_request(json, config), "graph");
    json["discovery_fraction"] = config.discovery_fraction;
    json["min_exact_match"] = config.alignment_min_exact_match;
    json["max_num_nodes_per_seq_char"] = config.alignment_max_nodes_per_seq_char;
    json["top_labels"] = Json::UInt64(config.num_top_labels);
    json["align"] = false;
    json["in_ram"] = true;
    EXPECT_EQ(make_search_cache_key(parse_search_request(json, config), "graph"), key);
    json["align"] = true;
    EXPECT_NE(make_search_cache_key(parse_search_request(json, config), "graph"), key);
    json["align"] = false;
    json["top_labels"] = 1;
    EXPECT_NE(make_search_cache_key(parse_search_request(json, config), "graph"), key);
}

TEST(ServerSearchRequest, RejectsInvalidOptionsBeforeCacheLookup) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    Json::Value json;
    EXPECT_THROW(parse_search_request(json, config), std::invalid_argument);
    json["FASTA"] = ">q\nACGT\n";
    for (const char *field : {"discovery_fraction", "min_exact_match"}) {
        for (double value : {-0.1, 1.1}) {
            auto invalid = json;
            invalid[field] = value;
            EXPECT_THROW(parse_search_request(invalid, config), std::invalid_argument);
        }
    }
    for (const char *field : {"top_labels", "max_num_nodes_per_seq_char", "align"}) {
        auto invalid = json;
        invalid[field] = Json::Value(Json::arrayValue);
        EXPECT_THROW(parse_search_request(invalid, config), std::exception);
    }
}


TEST(ServerQueryCache, EmptyResultsConsumeBudget) {
    ServerQueryCache cache(4096);
    for (size_t i = 0; i < 1000; ++i) {
        auto h = cache.acquire(std::to_string(i));
        h.set_result("");
    }
    EXPECT_GT(cache.size_bytes(), 0u);
    EXPECT_LE(cache.size_bytes(), 4096u);
    EXPECT_LT(cache.entry_count(), 128u);
}

TEST(ServerQueryCache, LargeKeysConsumeBudget) {
    ServerQueryCache cache(4096);
    {
        auto h = cache.acquire(std::string(4096, 'A'));
        h.set_result("small");
        EXPECT_EQ(cache.size_bytes(), 0u);  // too large to retain
        EXPECT_EQ(*h.get(), "small");
    }
    EXPECT_EQ(cache.entry_count(), 0u);
    EXPECT_EQ(cache.size_bytes(), 0u);
}

TEST(ServerQueryCache, DistinctFastaHeadersDoNotShareResults) {
    char program[] = "metagraph", command[] = "server_query", input[] = "fixture.csv";
    char *argv[] = {program, command, input};
    Config config(3, argv);
    Json::Value first, second;
    // Equal-length headers that collide under 64-bit libstdc++ std::hash.
    first["FASTA"] = ">aaaaaaabbbbbbbb\nACGTACGT\n";
    second["FASTA"] = ">y23O=^QD5fntS2F\nACGTACGT\n";
    EXPECT_NE(make_search_cache_key(parse_search_request(first, config), "graph"),
              make_search_cache_key(parse_search_request(second, config), "graph"));
}

#ifndef _NO_DEATH_TEST
TEST(ServerCacheConfig, RejectsInvalidSizes) {
    for (const char *value : {"-1", "nan", "inf", "1e309", "1e20", "1GB", "bad"}) {
        SCOPED_TRACE(value);
        char program[] = "metagraph", command[] = "server_query";
        char option[] = "--cache-size", input[] = "fixture.csv";
        char *argv[] = {program, command, option, const_cast<char *>(value), input};
        EXPECT_EXIT({
            Config config(5, argv);
            exit(0);
        }, ::testing::ExitedWithCode(1), "--cache-size must be");
    }
}
#endif

TEST(ServerCacheConfig, AcceptsZeroAndFractionalSizes) {
    for (const char *value : {"0", "0.5", "1"}) {
        char program[] = "metagraph", command[] = "server_query";
        char option[] = "--cache-size", input[] = "fixture.csv";
        char *argv[] = {program, command, option, const_cast<char *>(value), input};
        Config config(5, argv);
        EXPECT_DOUBLE_EQ(config.server_cache_size, std::stod(value));
    }
}


// A multi-graph /search tags its rows with the name a pair was selected under: the cache
// identity holds the names beside the pairs, so a pair listed under two names is two entries
TEST(ServerSearchRequest, PairsIdentityHoldsTheNames) {
    const std::pair<std::string, std::string> pair { "graph.dbg", "anno.annodbg" };
    const std::string left = search_pairs_identity({ { "left", pair } });
    const std::string right = search_pairs_identity({ { "right", pair } });
    EXPECT_NE(left, right);
    EXPECT_EQ(left, search_pairs_identity({ { "left", pair } }));
    // the pairs searched, in their order
    const std::pair<std::string, std::string> other { "other.dbg", "anno.annodbg" };
    EXPECT_NE(search_pairs_identity({ { "a", pair }, { "a", other } }),
              search_pairs_identity({ { "a", other }, { "a", pair } }));
    // every string length-prefixed: no two selections spell one identity
    EXPECT_NE(search_pairs_identity({ { "ab", { "c", "d" } } }),
              search_pairs_identity({ { "a", { "bc", "d" } } }));
    EXPECT_NE(search_pairs_identity({ { "a", { "b", "cd" } } }),
              search_pairs_identity({ { "a", { "bc", "d" } } }));
}

} // namespace
