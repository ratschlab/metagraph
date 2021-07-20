#ifndef __MTG_CACHES__
#define __MTG_CACHES__

#include <list>
#include <optional>
#include <utility>

#include <tsl/hopscotch_map.h>

namespace mtg {
namespace common {

// A fast LRU cache
template <typename Key, typename Value>
class ThreadUnsafeLRUCache {
  public:
    ThreadUnsafeLRUCache(size_t max_size) : max_size_(max_size) {}

    std::optional<Value> TryGet(const Key &key) {
        auto find = finder_.find(key);
        if (find == finder_.end())
            return std::nullopt;

        values_.splice(values_.begin(), values_, find->second);
        return find->second->second;
    }

    void Put(const Key &key, Value value) {
        auto find = finder_.find(key);
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
        values_.clear();
        finder_.clear();
    }

  private:
    size_t max_size_;

    // TODO: switch to a more efficient storage vector
    std::list<std::pair<Key, Value>> values_;

    tsl::hopscotch_map<Key, typename std::list<std::pair<Key, Value>>::iterator> finder_;
};

} // namespace common
} // namespace mtg

#endif // __MTG_CACHES__
