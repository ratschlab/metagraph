#ifndef __MTG_CACHES__
#define __MTG_CACHES__

#include <list>
#include <optional>
#include <mutex>
#include <shared_mutex>
#include <utility>

#include <tsl/hopscotch_map.h>

namespace mtg {
namespace common {

/**
 * A generic locking LRU cache.
 * The underlying data structures are a list and a map.
 * The list should store (Key, Value) pairs and support splice without invalidating iterators.
 * The map type should map from each key to an iterator in this list.
 * Each Put and TryGet operation marks the associated (key, value) pair as being
 * recently used.
 */
template <typename Key,
          typename Value,
          typename Storage = std::list<std::pair<Key, Value>>,
          typename Map = tsl::hopscotch_map<Key, typename Storage::iterator>>
class LRUCache {
  public:
    LRUCache(size_t max_size) : max_size_(max_size) {}

    // If the key is present in the cache, return its associated value, otherwise,
    // return std::nullopt. This also marks the (key, value) pair as being recently used.
    std::optional<Value> TryGet(const Key &key) {
        typename Map::iterator find;
        std::optional<Value> ret_val;
        {
            std::shared_lock<std::shared_mutex> read_lock(mu_);
            find = finder_.find(key);
            if (find == finder_.end())
                return std::nullopt;

            // no need to promote the key, so we can avoid a unique_lock
            if (find->second == values_.begin())
                return find->second->second;

            ret_val = find->second->second;
        }

        std::unique_lock<std::shared_mutex> promote_lock(mu_);
        find = finder_.find(key);

        // If a Put was called between releasing the shared_lock and acquiring
        // the unique_lock, then we need to add the (key, value) pair back in
        if (find == finder_.end()) {
            Put_nolock(key, *ret_val);
            find = finder_.find(key);
        }

        values_.splice(values_.begin(), values_, find->second);
        return ret_val;
    }

    // Add a (key, value) pair to the front of the cache.
    void Put(const Key &key, Value value) {
        std::unique_lock<std::shared_mutex> put_lock(mu_);
        Put_nolock(key, std::move(value));
    }

    void Clear() {
        std::unique_lock<std::shared_mutex> clear_lock(mu_);
        values_.clear();
        finder_.clear();
    }

  private:
    size_t max_size_;

    Storage values_;
    Map finder_;

    std::shared_mutex mu_;

    void Put_nolock(const Key &key, Value value) {
        auto [it, inserted] = finder_.emplace(key, values_.end());
        if (inserted) {
            values_.insert(values_.begin(), std::make_pair(key, std::move(value)));
            it.value() = values_.begin();
        } else {
            it.value()->second = std::move(value);
            values_.splice(values_.begin(), values_, it->second);
        }

        if (values_.size() > max_size_) {
            finder_.erase(values_.back().first);
            values_.pop_back();
        }
    }
};

} // namespace common
} // namespace mtg

#endif // __MTG_CACHES__
