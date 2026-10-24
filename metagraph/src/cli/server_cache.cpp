#include "server_cache.hpp"

#include <exception>
#include <iomanip>
#include <limits>
#include <locale>
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
        release();
        cache_ = other.cache_;
        entry_ = std::move(other.entry_);
        producer_ = std::move(other.producer_);
        other.cache_ = nullptr;
    }
    return *this;
}

ServerQueryCache::Handle::~Handle() {
    release();
}

void ServerQueryCache::Handle::release() {
    if (!entry_)
        return;
    // Producer abandoned without publishing (e.g. throw-before-set_result
    // along an unexpected path): satisfy the promise with an exception so
    // duplicate waiters don't block forever on the shared_future.
    if (producer_) {
        try {
            set_exception(std::make_exception_ptr(
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

bool ServerQueryCache::Handle::set_result(std::string response) {
    auto ptr = std::make_shared<const std::string>(std::move(response));
    size_t bytes = ptr->capacity() + 1;
    producer_->set_value(ptr);
    producer_.reset();
    return cache_->on_result_ready(entry_, bytes);
}

void ServerQueryCache::Handle::set_exception(std::exception_ptr eptr) {
    // Existing waiters own this entry. Remove it from lookup before waking
    // them so retries can start a new computation, even while they unwind.
    cache_->discard(entry_);
    producer_->set_exception(eptr);
    producer_.reset();
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

std::pair<std::shared_ptr<std::promise<ServerQueryCache::ResultPtr>>,
          std::shared_ptr<ServerQueryCache::Entry>>
ServerQueryCache::make_pending_entry(const std::string &key) {
    auto producer = std::make_shared<std::promise<ResultPtr>>();
    auto entry = std::make_shared<Entry>();
    entry->key = key;
    entry->future = producer->get_future().share();
    entry->waiters.store(1, std::memory_order_relaxed);
    return { std::move(producer), std::move(entry) };
}

ServerQueryCache::Handle ServerQueryCache::acquire(const std::string &key) {
    if (max_size_bytes_ == 0) {
        // Disabled-cache fast path: every caller gets a fresh, detached
        // entry. No dedup, no retention, no LRU bookkeeping.
        auto [producer, entry] = make_pending_entry(key);
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

    // Miss. Create a fresh entry, attach it to the cache, and hand the
    // promise to the caller (the producer) to fulfil via set_result().
    // Eviction-bookkeeping (size accounting) happens later in
    // on_result_ready, when the producer publishes.
    auto [producer, entry] = make_pending_entry(key);
    entry->in_cache = true;
    lru_.push_front(entry);
    entry->lru_pos = lru_.begin();
    map_.emplace(key, entry);
    return Handle(this, std::move(entry), std::move(producer));
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

bool ServerQueryCache::on_result_ready(const std::shared_ptr<Entry> &entry,
                                       size_t size) {
    // Called by the producer after publishing the value. Records when it became
    // available, registers its size in the byte budget, and runs the
    // size-pressure sweep — the just-inserted entry itself is protected
    // from eviction in this pass because its producer's Handle is still
    // alive (waiters ≥ 1).
    std::lock_guard<std::mutex> lock(mutex_);
    if (!entry->in_cache)
        return false;  // Detached entry from the disabled-cache path.
    entry->ready_at = std::chrono::steady_clock::now();
    // Charge both stored keys, the response allocation and an estimate of
    // entry/container bookkeeping, so even empty responses consume budget.
    size_t charged_size = size + 2 * (entry->key.capacity() + 1)
                          + sizeof(Entry) + sizeof(decltype(map_)::value_type)
                          + sizeof(std::string) + 8 * sizeof(void *);
    // An entry that cannot fit would evict useful results, then itself.
    if (charged_size > max_size_bytes_) {
        discard_locked(entry);
        return false;
    }
    entry->size_bytes = charged_size;
    total_size_bytes_ += entry->size_bytes;
    evict_under_pressure_locked();
    return true;
}

void ServerQueryCache::discard(const std::shared_ptr<Entry> &entry) {
    std::lock_guard<std::mutex> lock(mutex_);
    discard_locked(entry);
}

void ServerQueryCache::discard_locked(const std::shared_ptr<Entry> &entry) {
    if (!entry->in_cache)
        return;
    total_size_bytes_ -= entry->size_bytes;
    map_.erase(entry->key);
    lru_.erase(entry->lru_pos);
    entry->in_cache = false;
}

void ServerQueryCache::on_delivery(const std::shared_ptr<Entry> &entry,
                                   DeliveryState state) {
    if (state == DeliveryState::DELIVERED) {
        entry->delivery.store(DeliveryState::DELIVERED, std::memory_order_release);
    } else {
        // A failed delivery cannot undo an earlier success, even concurrently.
        auto expected = DeliveryState::PENDING;
        entry->delivery.compare_exchange_strong(expected, DeliveryState::PROTECTED,
                                                 std::memory_order_acq_rel);
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

    // Each pass walks oldest first, skipping entries still in use.
    auto evict = [&](bool keep_protected) {
        for (auto it = lru_.end(); total_size_bytes_ > max_size_bytes_ && it != lru_.begin();) {
            --it;
            auto &entry = *it;
            if (entry->waiters.load(std::memory_order_relaxed) != 0
                    || (keep_protected && is_within_protection_window(*entry)))
                continue;
            total_size_bytes_ -= entry->size_bytes;
            map_.erase(entry->key);
            entry->in_cache = false;
            it = lru_.erase(it);
        }
    };
    // Prefer keeping failed deliveries, but evict them too if needed to fit.
    evict(true);
    evict(false);
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
    // semantic result. Keep the complete FASTA so a hash collision in
    // the lookup table cannot return a different request's response.
    //
    // Per-request overrides are pulled directly from the JSON the same
    // way process_search_request resolves them, so flags that select
    // different query modes (e.g. abundance_sum vs query_coords) give
    // different keys for the same FASTA.
    //
    // NOTE: keep this resolution logic in sync with the JSON-to-Config
    // mapping in `process_search_request` (server.cpp).
    const std::string &fasta = json["FASTA"].asString();

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
    oss.imbue(std::locale::classic());
    oss << std::setprecision(std::numeric_limits<double>::max_digits10);
    oss << "g:" << graph_identity.size() << ':' << graph_identity << '|'
        << "m:" << query_mode << '|'
        << "df:" << discovery_fraction << '|'
        << "ame:" << min_exact_match << '|'
        << "amn:" << max_nodes_per_seq_char << '|'
        << "tl:" << top_labels << '|'
        << "al:" << (json.get("align", false).asBool() ? '1' : '0') << '|'
        << "f:" << fasta.size() << ':' << fasta;
    return oss.str();
}

} // namespace cli
} // namespace mtg
