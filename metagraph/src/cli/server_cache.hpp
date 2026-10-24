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


namespace mtg {
namespace cli {

struct SearchRequest;

/**
 * Cache serialized /search responses and share in-flight computations.
 * The producer publishes a result or exception; failed entries leave lookup
 * immediately so later requests can retry. Handles pin entries until delivery.
 *
 * Under memory pressure, evict unpinned entries in LRU order. Failed deliveries
 * get preference for protection_ttl after publication or their last cache hit,
 * but are still evicted if necessary. A successful delivery ends that preference
 * permanently, even if another delivery fails later.
 *
 * The budget estimates string allocations and entry overhead. Pinned entries
 * can temporarily exceed it. Oversized results are shared with existing waiters
 * but not retained. A zero budget disables retention and deduplication.
 * Compression is applied per response, after reading the cached JSON body.
 */
class ServerQueryCache {
  public:
    enum class DeliveryState { PENDING, DELIVERED, PROTECTED };

    // The cache value is the serialized JSON response body.
    using ResultPtr = std::shared_ptr<const std::string>;
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
        // Returns whether the result fits the retention budget.
        bool set_result(std::string response);
        void set_exception(std::exception_ptr eptr);

        // Record delivery success or failure. Success is permanent; failures
        // get retention preference without extending the window on their own.
        void mark_delivered();
        void mark_protected();

      private:
        friend class ServerQueryCache;

        Handle(ServerQueryCache *cache,
               std::shared_ptr<Entry> entry,
               std::shared_ptr<std::promise<ResultPtr>> producer);

        void release();

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

  private:
    friend class ServerQueryCacheTestAccess;

    struct Entry {
        std::string key;
        // Hinge of the dedup mechanism: producer publishes via promise,
        // every waiter (including late arrivals) reads via .get().
        ResultFuture future;
        // RAII-managed by Handle: producer + each reader contributes 1.
        // An entry with waiters > 0 is never evicted.
        std::atomic<int>      waiters{0};
        // Set when the producer publishes a result (on_result_ready),
        // and refreshed on every cache hit while
        // the entry is in PROTECTED state (sliding retention window).
        // Used to compute "within protection_ttl" for eviction priority.
        std::chrono::steady_clock::time_point ready_at{};
        std::atomic<DeliveryState> delivery{DeliveryState::PENDING};
        // Charged string allocations and estimated entry overhead.
        // Set once in on_result_ready; used without touching shared_future.
        size_t                size_bytes = 0;
        // Iterator into ServerQueryCache::lru_; valid iff in_cache.
        std::list<std::shared_ptr<Entry>>::iterator lru_pos{};
        // false means this entry is detached (evicted, or made by the
        // disabled-cache fast path); accounting/eviction skip it.
        bool in_cache = false;
    };

    // Allocate a new Entry + paired promise/future. Caller decides
    // whether to attach it to the cache (in_cache, lru_, map_) or use
    // it as a one-shot detached entry (disabled-cache fast path).
    std::pair<std::shared_ptr<std::promise<ResultPtr>>, std::shared_ptr<Entry>>
        make_pending_entry(const std::string &key);

    void release_waiter(const std::shared_ptr<Entry> &entry);
    bool on_result_ready(const std::shared_ptr<Entry> &entry, size_t size);
    void discard(const std::shared_ptr<Entry> &entry);
    void discard_locked(const std::shared_ptr<Entry> &entry);
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


// Key the same resolved options used to execute the search. Transport encoding
// and graph-loading strategy do not change the result and are excluded.
std::string make_search_cache_key(const SearchRequest &request,
                                  const std::string &graph_identity);

} // namespace cli
} // namespace mtg

#endif // __SERVER_CACHE_HPP__
