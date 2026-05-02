#include "server_cache.hpp"

#include <algorithm>
#include <exception>
#include <sstream>
#include <stdexcept>

#include "common/logger.hpp"
#include "config/config.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;


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
    // If we were the producer and never published, fail the promise so
    // any other waiters get an exception rather than hanging.
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
    // Note: we deliberately do NOT remove the entry here. Other waiters
    // (if any) will receive the exception via future.get(); the entry
    // becomes evictable once all waiters drop. For freshly-thrown errors
    // with no other waiters, the entry will be released by ~Handle and
    // pruned by the next eviction sweep.
    cache_->on_result_ready(entry_, /* approx_size = */ 0);
}

void ServerQueryCache::Handle::mark_delivered() {
    cache_->on_delivery(entry_, DeliveryState::DELIVERED);
}

void ServerQueryCache::Handle::mark_failed() {
    cache_->on_delivery(entry_, DeliveryState::FAILED);
}


// --- ServerQueryCache -------------------------------------------------------

ServerQueryCache::ServerQueryCache(size_t max_size_bytes,
                                   std::chrono::nanoseconds failed_ttl)
      : max_size_bytes_(max_size_bytes),
        failed_ttl_(failed_ttl) {}

ServerQueryCache::Handle ServerQueryCache::acquire(const std::string &key) {
    if (max_size_bytes_ == 0) {
        // Cache disabled: produce a one-shot, untracked entry. Every call
        // is a miss; no dedup, no retention.
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
        auto &entry = it->second;
        entry->waiters.fetch_add(1, std::memory_order_relaxed);
        touch_lru_locked(entry);
        return Handle(this, entry, /* producer = */ nullptr);
    }

    // Miss: create a fresh entry with a shared promise/future pair.
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
        return;  // Still has waiters.

    // Last waiter dropped — entry is now eligible for size-pressure / TTL
    // eviction. We don't proactively evict here; the next acquire() that
    // triggers pressure will reclaim space. But we do prune obviously
    // expired FAILED entries opportunistically.
    std::lock_guard<std::mutex> lock(mutex_);
    evict_expired_failed_locked();
}

void ServerQueryCache::on_result_ready(const std::shared_ptr<Entry> &entry,
                                       size_t approx_size) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!entry->in_cache)
        return;
    entry->ready_at = std::chrono::steady_clock::now();
    entry->approx_size_bytes = approx_size;
    total_size_bytes_ += approx_size;
    evict_under_pressure_locked();
}

void ServerQueryCache::on_delivery(const std::shared_ptr<Entry> &entry,
                                   DeliveryState state) {
    entry->delivery.store(state, std::memory_order_release);
    if (state == DeliveryState::FAILED) {
        // Refresh ready_at to start the 2h TTL from the failure observation
        // moment if we don't have one yet (set_result wasn't called).
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
    while (total_size_bytes_ > max_size_bytes_) {
        // Walk LRU from oldest (back) to newest (front). Skip entries with
        // active waiters — they're being computed or read right now.
        auto rit = std::find_if(lru_.rbegin(), lru_.rend(),
            [](const std::shared_ptr<Entry> &e) {
                return e->waiters.load(std::memory_order_relaxed) == 0;
            });
        if (rit == lru_.rend())
            break;  // Every remaining entry has waiters; nothing else to free.

        auto &entry = *rit;
        total_size_bytes_ -= entry->approx_size_bytes;
        map_.erase(entry->key);
        entry->in_cache = false;
        // Convert reverse iterator to forward iterator for erase.
        lru_.erase(std::next(rit).base());
    }
}

void ServerQueryCache::evict_expired_failed_locked() {
    auto now = std::chrono::steady_clock::now();
    for (auto it = lru_.begin(); it != lru_.end(); ) {
        auto &entry = *it;
        bool expired =
            entry->waiters.load(std::memory_order_relaxed) == 0
            && entry->delivery.load(std::memory_order_acquire) != DeliveryState::DELIVERED
            && entry->ready_at.time_since_epoch().count() != 0
            && now - entry->ready_at > failed_ttl_;
        if (!expired) {
            ++it;
            continue;
        }
        total_size_bytes_ -= entry->approx_size_bytes;
        map_.erase(entry->key);
        entry->in_cache = false;
        it = lru_.erase(it);
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
    // We hash the FASTA body (the bulk of the request) and concatenate the
    // remaining flag-sized inputs verbatim. Per-request overrides are read
    // from the JSON the same way `process_search_request` resolves them,
    // so every input that affects the semantic result is in the key.
    //
    // NOTE: keep this in sync with the JSON->Config resolution in
    // `process_search_request` (server.cpp).
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
