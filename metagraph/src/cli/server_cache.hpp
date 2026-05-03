#ifndef __SERVER_CACHE_HPP__
#define __SERVER_CACHE_HPP__

#include <atomic>
#include <chrono>
#include <future>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

#include <json/json.h>


namespace mtg {
namespace cli {

class Config;

/**
 * A bounded result cache for the `/search` server endpoint.
 *
 * The cache solves three problems:
 *   1. Concurrent identical requests share one computation
 *      (in-flight dedup via std::shared_future).
 *   2. An upstream retry after a connection drop hits the cached result
 *      instead of recomputing.
 *   3. Failed deliveries (TCP error or no handshake) are retained with
 *      higher priority than normal results for a 2h retry window, so a
 *      flood of small successful requests can't displace responses that
 *      an upstream is still about to retry.
 *
 * Eviction: LRU once the total approximate response size exceeds the
 * configured max. Entries with active waiters (currently being computed
 * or still being read) are never evicted. FAILED-delivery entries have
 * *higher* retention priority than DELIVERED ones for `failed_ttl`
 * (default 2h) — they are sacrificed only after every waiterless
 * DELIVERED (and past-TTL FAILED) entry has been evicted, so a flood of
 * tiny normal requests can't displace responses an upstream is still
 * about to retry. Past `failed_ttl` they are evictable like any other
 * waiterless entry.
 *
 * Cache values are the final response Json::Value. The server's
 * `verbose_output` flag is process-constant (set on startup), so the
 * cached JSON is consistent across all hits within a process; only
 * Accept-Encoding-driven compression is reapplied per request.
 */
class ServerQueryCache {
  public:
    enum class DeliveryState { PENDING, DELIVERED, FAILED };

    struct CachedResult {
        Json::Value response;
        size_t      approx_size_bytes = 0;
    };

    using ResultPtr = std::shared_ptr<const CachedResult>;
    using ResultFuture = std::shared_future<ResultPtr>;

  private:
    // Defined here so `Handle` can refer to `Entry` before its full body.
    struct Entry;

  public:
    /**
     * Handle returned by `acquire()`. RAII-decrements the entry's waiter
     * count on destruction. On a cache miss `is_miss()` is true and the
     * caller must call exactly one of `set_result(...)` / `set_exception(...)`.
     */
    class Handle {
      public:
        Handle() = default;
        Handle(const Handle &) = delete;
        Handle &operator=(const Handle &) = delete;
        Handle(Handle &&) noexcept;
        Handle &operator=(Handle &&) noexcept;
        ~Handle();

        // True if this handle holds an entry (i.e. acquire() ran).
        explicit operator bool() const { return static_cast<bool>(entry_); }

        // True if this caller is responsible for computing the result.
        bool is_miss() const { return static_cast<bool>(producer_); }

        // Wait for the result. Throws whatever the producer set via set_exception.
        ResultPtr get() const;

        // Producer-only: publish the result to all waiters.
        void set_result(CachedResult result);
        void set_exception(std::exception_ptr eptr);

        // Mark the eventual delivery outcome of the response.
        //
        // mark_delivered() is a sink: once an entry has been successfully
        // delivered to any client, it stays DELIVERED — a later mark_failed()
        // on the same entry (e.g. from a duplicate request whose delivery
        // dropped) will NOT resurrect the TTL clock.
        //
        // mark_failed() arms the failed-delivery TTL on the first transition
        // out of PENDING, and (if the producer didn't already publish a
        // result) records the moment for the TTL clock.
        void mark_delivered();
        void mark_failed();

      private:
        friend class ServerQueryCache;

        Handle(ServerQueryCache *cache,
               std::shared_ptr<Entry> entry,
               std::shared_ptr<std::promise<ResultPtr>> producer);

        ServerQueryCache *cache_ = nullptr;
        std::shared_ptr<Entry> entry_;
        std::shared_ptr<std::promise<ResultPtr>> producer_;  // non-null on miss only
    };

    explicit ServerQueryCache(size_t max_size_bytes,
                              std::chrono::nanoseconds failed_ttl = std::chrono::hours(2));

    /**
     * Look up `key` and bump its waiter count. On a miss, the returned
     * handle's `is_miss()` is true and the caller is responsible for
     * computing the result and publishing it via `set_result()`.
     */
    Handle acquire(const std::string &key);

    // Diagnostics / tests.
    size_t size_bytes() const;
    size_t entry_count() const;
    bool contains(const std::string &key) const;
    std::chrono::nanoseconds failed_ttl() const { return failed_ttl_; }

  private:
    struct Entry {
        std::string key;
        ResultFuture future;
        std::atomic<int>      waiters{0};
        std::chrono::steady_clock::time_point ready_at{};
        std::atomic<DeliveryState> delivery{DeliveryState::PENDING};
        size_t                approx_size_bytes = 0;
        // Iterator into ServerQueryCache::lru_, valid while entry is in the cache.
        std::list<std::shared_ptr<Entry>>::iterator lru_pos{};
        bool in_cache = false;
    };

    void release_waiter(const std::shared_ptr<Entry> &entry);
    void on_result_ready(const std::shared_ptr<Entry> &entry, size_t approx_size);
    void on_delivery(const std::shared_ptr<Entry> &entry, DeliveryState state);
    void evict_under_pressure_locked();
    void evict_expired_failed_locked();
    void touch_lru_locked(const std::shared_ptr<Entry> &entry);

    const size_t max_size_bytes_;
    const std::chrono::nanoseconds failed_ttl_;

    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::shared_ptr<Entry>> map_;
    // LRU: front = MRU, back = LRU candidate.
    std::list<std::shared_ptr<Entry>> lru_;
    size_t total_size_bytes_ = 0;
};


/**
 * Build the canonical cache key for a `/search` request given the parsed
 * JSON body, the server's startup Config, and the resolved graph identity.
 *
 * Per-request overrides (discovery_fraction, query_mode flags, top_labels,
 * etc.) are read from `json` directly, mirroring how `process_search_request`
 * resolves them. The key includes only inputs that affect the *semantic*
 * result; it deliberately excludes formatting flags like `verbose_output`
 * and `Accept-Encoding: deflate`.
 */
std::string make_search_cache_key(const Json::Value &json,
                                  const Config &server_config,
                                  const std::string &graph_identity);

} // namespace cli
} // namespace mtg

#endif // __SERVER_CACHE_HPP__
