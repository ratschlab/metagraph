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

// A fast LRU cache
template <typename Key, typename Value>
class LRUCache {
  public:
    LRUCache(size_t max_size) : max_size_(max_size) {}

    std::optional<Value> TryGet(const Key &key) {
        typename Map::iterator find;
        {
            std::shared_lock<std::shared_mutex> read_lock(mu_);
            find = finder_.find(key);
            if (find == finder_.end())
                return std::nullopt;

            // no need to promote the key, so we can avoid a unique_lock
            if (find->second == values_.begin())
                return find->second->second;
        }

        std::unique_lock<std::shared_mutex> promote_lock(mu_);
        find = finder_.find(key);
        values_.splice(values_.begin(), values_, find->second);
        return find->second->second;
    }

    void Put(const Key &key, Value value) {
        std::unique_lock<std::shared_mutex> put_lock(mu_);

        typename Map::iterator find = finder_.find(key);
        if (find != finder_.end()) {
            find->second->second = std::move(value);
            values_.splice(values_.begin(), values_, find->second);
        } else {
            values_.insert(values_.begin(), std::make_pair(key, std::move(value)));
            finder_[key] = values_.begin();
        }

        if (values_.size() > max_size_) {
            finder_.erase(values_.back().first);
            values_.pop_back();
        }
    }

    void Clear() {
        std::unique_lock<std::shared_mutex> clear_lock(mu_);
        values_.clear();
        finder_.clear();
    }

  private:
    size_t max_size_;

    // TODO: switch to a more efficient storage vector
    typedef std::list<std::pair<Key, Value>> Storage;
    typedef tsl::hopscotch_map<Key, typename Storage::iterator> Map;

    Storage values_;
    Map finder_;

    std::shared_mutex mu_;
};

} // namespace common
} // namespace mtg

#endif // __MTG_CACHES__
