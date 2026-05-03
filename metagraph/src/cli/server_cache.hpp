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
 * Goals
 * -----
 *   1. Dedup concurrent identical requests so they share one computation.
 *   2. Survive connection drops: an upstream retry after a TCP failure
 *      hits the cached result instead of recomputing a multi-minute query.
 *   3. Bound memory while still favoring entries that look like they're
 *      being retried by a disconnected upstream.
 *
 * Per-request lifecycle
 * ---------------------
 *   1. The handler calls `acquire(key)` and gets a `Handle`. waiters++.
 *   2. On a miss the handler computes the result and calls
 *      `Handle::set_result(...)` (or `set_exception(...)` on a thrown
 *      compute). The `shared_future` underneath unblocks any concurrent
 *      duplicate caller that arrived during the computation.
 *      On a hit the handler simply reads the cached value via
 *      `Handle::get()`.
 *   3. The response is written to the wire. The async-write completion
 *      callback (`on_sent(error_code)`) calls
 *        ec == 0  → `Handle::mark_delivered()`
 *        ec != 0  → `Handle::mark_protected()`
 *   4. `~Handle()` decrements waiters; the last drop allows the entry
 *      to be considered for eviction.
 *
 * Entry states
 * ------------
 *   PENDING    Computation in flight. waiters ≥ 1. Never evicted.
 *   DELIVERED  At least one delivery succeeded. Sink state — a later
 *              `mark_protected()` is a no-op. Evicted only by LRU under
 *              size pressure.
 *   PROTECTED  No delivery has succeeded yet. The entry has elevated
 *              retention priority for `protection_ttl` (default 2 h).
 *              Every cache hit refreshes the window (sliding TTL — the
 *              hit itself is direct evidence the upstream is still
 *              retrying). After 2 h of no hits the entry graduates to
 *              the main cache: state stays PROTECTED but it loses
 *              priority and ages out via normal LRU.
 *
 * Eviction policy (size-pressure sweep, two passes)
 * -------------------------------------------------
 *   Pass 1: walk LRU back→front, evict the oldest waiterless entry that
 *           is NOT PROTECTED-within-window.
 *   Pass 2: still over budget? Walk again, evict any waiterless entry
 *           (PROTECTED-within-window included). Keeps the cache bounded
 *           when PROTECTED entries dominate.
 *   Entries with waiters > 0 are never evicted.
 *
 * Cache value
 * -----------
 *   The final response `Json::Value`. The server's `verbose_output` flag
 *   is process-constant (set on startup), so the cached JSON is
 *   consistent across all hits within a process; only Accept-Encoding-
 *   driven compression is reapplied per request.
 */
class ServerQueryCache {
  public:
    enum class DeliveryState { PENDING, DELIVERED, PROTECTED };

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
     * RAII handle returned by `acquire()`. The handle holds one waiter
     * count on the entry and drops it in the destructor.
     *
     * On a miss (is_miss() == true) this caller is the producer and
     * must call exactly one of `set_result(...)` / `set_exception(...)`
     * to release any concurrent duplicate waiters. On a hit just call
     * `get()` to read the cached value.
     *
     * After the response has been written to the wire, the on_sent
     * callback is expected to call mark_delivered() or mark_protected()
     * exactly once.
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

        // Wait for the result. Throws whatever the producer set via
        // set_exception(). On a miss must be called only after
        // set_result/set_exception, otherwise it deadlocks the producer.
        ResultPtr get() const;

        // Producer-only: publish the result (or thrown exception) to
        // all current and future waiters of this entry's shared_future.
        void set_result(CachedResult result);
        void set_exception(std::exception_ptr eptr);

        // Record the eventual delivery outcome of the response.
        //
        // mark_delivered() is a sink: once any delivery has succeeded
        // the entry stays DELIVERED, and a later mark_protected() on
        // the same entry (e.g. from a duplicate request whose delivery
        // dropped) is a no-op — a successful delivery is a permanent
        // fact.
        //
        // mark_protected() puts the entry under the priority-retention
        // window (or extends it if already there). The window timer is
        // also refreshed on each cache hit, since a hit is itself
        // evidence the upstream is still retrying.
        void mark_delivered();
        void mark_protected();

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
                              std::chrono::nanoseconds protection_ttl = std::chrono::hours(2));

    /**
     * Look up `key` and bump its waiter count.
     *
     * On a miss the returned handle's `is_miss()` is true and the
     * caller must produce the result.
     *
     * On a hit of a PROTECTED entry, the entry's retention window is
     * also refreshed — the hit is direct evidence the upstream is
     * still retrying.
     *
     * If the cache is disabled (max_size_bytes == 0) every call is a
     * miss and the entry is not retained after the handle drops.
     */
    Handle acquire(const std::string &key);

    // Diagnostics / tests.
    size_t size_bytes() const;
    size_t entry_count() const;
    bool contains(const std::string &key) const;
    std::chrono::nanoseconds protection_ttl() const { return protection_ttl_; }

  private:
    struct Entry {
        std::string key;
        // Hinge of the dedup mechanism: producer publishes via promise,
        // every waiter (including late arrivals) reads via .get().
        ResultFuture future;
        // RAII-managed by Handle: producer + each reader contributes 1.
        // An entry with waiters > 0 is never evicted.
        std::atomic<int>      waiters{0};
        // Set when set_result/set_exception/mark_protected first runs;
        // refreshed on every PROTECTED cache hit (sliding window).
        // Used to compute "within protection_ttl" for eviction priority.
        std::chrono::steady_clock::time_point ready_at{};
        std::atomic<DeliveryState> delivery{DeliveryState::PENDING};
        size_t                approx_size_bytes = 0;
        // Iterator into ServerQueryCache::lru_; valid iff in_cache.
        std::list<std::shared_ptr<Entry>>::iterator lru_pos{};
        // false means this entry is detached (evicted, or made by the
        // disabled-cache fast path); accounting/eviction skip it.
        bool in_cache = false;
    };

    void release_waiter(const std::shared_ptr<Entry> &entry);
    void on_result_ready(const std::shared_ptr<Entry> &entry, size_t approx_size);
    void on_delivery(const std::shared_ptr<Entry> &entry, DeliveryState state);
    void evict_under_pressure_locked();
    void touch_lru_locked(const std::shared_ptr<Entry> &entry);

    const size_t max_size_bytes_;
    const std::chrono::nanoseconds protection_ttl_;

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
