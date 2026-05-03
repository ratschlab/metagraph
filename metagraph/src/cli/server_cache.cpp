#include "server_cache.hpp"

#include <algorithm>
#include <exception>
#include <sstream>
#include <stdexcept>

#include "config/config.hpp"


namespace mtg {
namespace cli {


// --- Handle -----------------------------------------------------------------

ServerQueryCache::Handle::Handle(ServerQueryCache *cache,
                                 std::shared_ptr<Entry> entry,
                                 std::shared_ptr<std::promise<ResultPtr>> producer)
      : cache_(cache),
        entry_(std::move(entry)),
        producer_(std::move(producer)) {}

ServerQueryCache::Handle::Handle(Handle &&other) noexcept
      : cache_(other.cache_),
        entry_(std::move(other.entry_)),
        producer_(std::move(other.producer_)) {
    other.cache_ = nullptr;
}

ServerQueryCache::Handle &
ServerQueryCache::Handle::operator=(Handle &&other) noexcept {
    if (this != &other) {
        // Release any current entry first.
        if (entry_)
            cache_->release_waiter(entry_);
        cache_ = other.cache_;
        entry_ = std::move(other.entry_);
        producer_ = std::move(other.producer_);
        other.cache_ = nullptr;
    }
    return *this;
}

ServerQueryCache::Handle::~Handle() {
    if (!entry_)
        return;
    // Producer abandoned without publishing (e.g. throw-before-set_result
    // along an unexpected path): satisfy the promise with an exception so
    // duplicate waiters don't block forever on the shared_future.
    if (producer_) {
        try {
            producer_->set_exception(std::make_exception_ptr(
                std::runtime_error("Cache producer abandoned the request")));
        } catch (const std::future_error &) {
            // Already satisfied — fine.
        }
    }
    cache_->release_waiter(entry_);
}

ServerQueryCache::ResultPtr ServerQueryCache::Handle::get() const {
    return entry_->future.get();
}

void ServerQueryCache::Handle::set_result(CachedResult result) {
    auto ptr = std::make_shared<const CachedResult>(std::move(result));
    size_t approx = ptr->approx_size_bytes;
    producer_->set_value(ptr);
    producer_.reset();
    cache_->on_result_ready(entry_, approx);
}

void ServerQueryCache::Handle::set_exception(std::exception_ptr eptr) {
    producer_->set_exception(eptr);
    producer_.reset();
    // Don't remove the entry — concurrent waiters will receive the
    // exception via future.get(). The entry becomes evictable once
    // every waiter drops; subsequent identical requests will hit it,
    // re-throw the same error, and behave consistently.
    cache_->on_result_ready(entry_, /* approx_size = */ 0);
}

void ServerQueryCache::Handle::mark_delivered() {
    cache_->on_delivery(entry_, DeliveryState::DELIVERED);
}

void ServerQueryCache::Handle::mark_protected() {
    cache_->on_delivery(entry_, DeliveryState::PROTECTED);
}


// --- ServerQueryCache -------------------------------------------------------

ServerQueryCache::ServerQueryCache(size_t max_size_bytes,
                                   std::chrono::nanoseconds protection_ttl)
      : max_size_bytes_(max_size_bytes),
        protection_ttl_(protection_ttl) {}

ServerQueryCache::Handle ServerQueryCache::acquire(const std::string &key) {
    if (max_size_bytes_ == 0) {
        // Disabled-cache fast path: every caller gets a fresh, detached
        // entry. No dedup, no retention, no LRU bookkeeping.
        auto producer = std::make_shared<std::promise<ResultPtr>>();
        auto entry = std::make_shared<Entry>();
        entry->key = key;
        entry->future = producer->get_future().share();
        entry->waiters.store(1, std::memory_order_relaxed);
        entry->in_cache = false;
        return Handle(this, std::move(entry), std::move(producer));
    }

    std::lock_guard<std::mutex> lock(mutex_);

    auto it = map_.find(key);
    if (it != map_.end()) {
        // Hit. Bump waiters, move to MRU.
        auto &entry = it->second;
        entry->waiters.fetch_add(1, std::memory_order_relaxed);
        touch_lru_locked(entry);
        // Sliding retention window: a hit on a PROTECTED entry is
        // direct evidence the upstream is still retrying, so we
        // refresh ready_at to extend the priority window.
        if (entry->delivery.load(std::memory_order_acquire) == DeliveryState::PROTECTED) {
            entry->ready_at = std::chrono::steady_clock::now();
        }
        return Handle(this, entry, /* producer = */ nullptr);
    }

    // Miss. Create a fresh entry with a shared_future paired to a
    // promise the caller (the producer) will fulfil via set_result().
    // Eviction-bookkeeping (size, LRU) happens later in on_result_ready.
    auto producer = std::make_shared<std::promise<ResultPtr>>();
    auto entry = std::make_shared<Entry>();
    entry->key = key;
    entry->future = producer->get_future().share();
    entry->waiters.store(1, std::memory_order_relaxed);
    entry->in_cache = true;
    lru_.push_front(entry);
    entry->lru_pos = lru_.begin();
    map_.emplace(key, entry);
    return Handle(this, entry, std::move(producer));
}

void ServerQueryCache::release_waiter(const std::shared_ptr<Entry> &entry) {
    int prev = entry->waiters.fetch_sub(1, std::memory_order_acq_rel);
    if (prev > 1)
        return;  // Still has other waiters.

    // Last waiter dropped: this entry — and possibly other entries that
    // were skipped during a previous insert because they had waiters —
    // are now evictable. Run the sweep so the cache settles back under
    // budget without waiting for the next insert.
    std::lock_guard<std::mutex> lock(mutex_);
    evict_under_pressure_locked();
}

void ServerQueryCache::on_result_ready(const std::shared_ptr<Entry> &entry,
                                       size_t approx_size) {
    // Called by the producer after publishing the value (set_result) or
    // exception (set_exception). Records the moment the result became
    // available, registers its size in the byte budget, and runs the
    // size-pressure sweep — the just-inserted entry itself is protected
    // from eviction in this pass because its producer's Handle is still
    // alive (waiters ≥ 1).
    std::lock_guard<std::mutex> lock(mutex_);
    if (!entry->in_cache)
        return;  // Detached entry from the disabled-cache path.
    entry->ready_at = std::chrono::steady_clock::now();
    entry->approx_size_bytes = approx_size;
    total_size_bytes_ += approx_size;
    evict_under_pressure_locked();
}

void ServerQueryCache::on_delivery(const std::shared_ptr<Entry> &entry,
                                   DeliveryState state) {
    // Called from the on_sent async-write callback exactly once per
    // delivery attempt. The caller passes either DELIVERED (write
    // succeeded) or PROTECTED (write failed → keep the response cached
    // with priority for the configured retry window).
    if (state == DeliveryState::DELIVERED) {
        // Sink: once any delivery has succeeded, the entry is DELIVERED
        // forever. A later mark_protected() (e.g. from a duplicate
        // request whose delivery dropped) is a no-op below.
        entry->delivery.store(DeliveryState::DELIVERED, std::memory_order_release);
        return;
    }
    // state == PROTECTED. CAS-loop the transition so we don't downgrade
    // a DELIVERED entry to PROTECTED.
    DeliveryState current = entry->delivery.load(std::memory_order_acquire);
    bool transitioned_to_protected = false;
    while (current != DeliveryState::DELIVERED) {
        if (entry->delivery.compare_exchange_weak(current, DeliveryState::PROTECTED,
                                                  std::memory_order_acq_rel)) {
            transitioned_to_protected = true;
            break;
        }
    }
    if (transitioned_to_protected) {
        // First entry into PROTECTED: arm ready_at if no result was
        // ever published (set_exception path; on_result_ready would
        // have set it otherwise). Subsequent transitions are no-ops
        // here — the sliding window is refreshed on each cache hit
        // in acquire().
        std::lock_guard<std::mutex> lock(mutex_);
        if (entry->ready_at.time_since_epoch().count() == 0)
            entry->ready_at = std::chrono::steady_clock::now();
    }
}

void ServerQueryCache::touch_lru_locked(const std::shared_ptr<Entry> &entry) {
    if (!entry->in_cache)
        return;
    if (entry->lru_pos != lru_.begin())
        lru_.splice(lru_.begin(), lru_, entry->lru_pos);
    entry->lru_pos = lru_.begin();
}

void ServerQueryCache::evict_under_pressure_locked() {
    auto now = std::chrono::steady_clock::now();
    auto is_within_protection_window = [&](const Entry &e) {
        // PROTECTED entries within the retention window have *higher*
        // priority than DELIVERED ones — they're sacrificed only when
        // there's nothing else to give up.
        return e.delivery.load(std::memory_order_acquire) == DeliveryState::PROTECTED
               && e.ready_at.time_since_epoch().count() != 0
               && now - e.ready_at <= protection_ttl_;
    };

    // Walk LRU back→front (oldest first), evict the first waiterless
    // entry that satisfies `pred`. Returns whether anything was evicted.
    auto evict_one = [&](auto pred) {
        auto rit = std::find_if(lru_.rbegin(), lru_.rend(),
            [&](const std::shared_ptr<Entry> &e) {
                return e->waiters.load(std::memory_order_relaxed) == 0 && pred(*e);
            });
        if (rit == lru_.rend())
            return false;
        auto &entry = *rit;
        total_size_bytes_ -= entry->approx_size_bytes;
        map_.erase(entry->key);
        entry->in_cache = false;
        lru_.erase(std::next(rit).base());
        return true;
    };

    // Pass 1: sacrifice DELIVERED and out-of-window PROTECTED entries
    // first; in-window PROTECTED entries are held back.
    while (total_size_bytes_ > max_size_bytes_
            && evict_one([&](const Entry &e) { return !is_within_protection_window(e); })) {
    }
    // Pass 2: only triggered when in-window PROTECTED entries are the
    // only waiterless candidates left (e.g. cache flooded with retries
    // for many distinct failed requests). Sacrifice them in LRU order
    // so the cache stays bounded.
    while (total_size_bytes_ > max_size_bytes_
            && evict_one([](const Entry &) { return true; })) {
    }
}

size_t ServerQueryCache::size_bytes() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return total_size_bytes_;
}

size_t ServerQueryCache::entry_count() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return map_.size();
}

bool ServerQueryCache::contains(const std::string &key) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return map_.count(key) != 0;
}


// --- key construction -------------------------------------------------------

std::string make_search_cache_key(const Json::Value &json,
                                  const Config &server_config,
                                  const std::string &graph_identity) {
    // Build a stable string key out of every input that affects the
    // semantic result. The FASTA body (the bulk of the request) is
    // collapsed to a hash; flag-sized values are inlined verbatim.
    //
    // Per-request overrides are pulled directly from the JSON the same
    // way process_search_request resolves them, so flags that select
    // different query modes (e.g. abundance_sum vs query_coords) give
    // different keys for the same FASTA.
    //
    // NOTE: keep this resolution logic in sync with the JSON-to-Config
    // mapping in `process_search_request` (server.cpp).
    const std::string &fasta = json["FASTA"].asString();
    size_t fasta_hash = std::hash<std::string>{}(fasta);

    double discovery_fraction
            = json.get("discovery_fraction", server_config.discovery_fraction).asDouble();
    double min_exact_match
            = json.get("min_exact_match", server_config.alignment_min_exact_match).asDouble();
    double max_nodes_per_seq_char
            = json.get("max_num_nodes_per_seq_char",
                       server_config.alignment_max_nodes_per_seq_char).asDouble();
    int top_labels
            = json.get("top_labels", server_config.num_top_labels).asInt();

    int query_mode;
    if (json.get("query_coords", false).asBool()) {
        query_mode = COORDS;
    } else if (json.get("query_counts", false).asBool()) {
        query_mode = COUNTS;
    } else if (json.get("with_signature", false).asBool()) {
        query_mode = SIGNATURE;
    } else if (json.get("abundance_sum", false).asBool()) {
        query_mode = COUNTS_SUM;
    } else {
        query_mode = MATCHES;
    }

    std::ostringstream oss;
    oss << "g:" << graph_identity << '|'
        << "m:" << query_mode << '|'
        << "df:" << discovery_fraction << '|'
        << "ame:" << min_exact_match << '|'
        << "amn:" << max_nodes_per_seq_char << '|'
        << "tl:" << top_labels << '|'
        << "al:" << (json.get("align", false).asBool() ? '1' : '0') << '|'
        << "fl:" << fasta.size() << '|'
        << "fh:" << fasta_hash;
    return oss.str();
}

} // namespace cli
} // namespace mtg
