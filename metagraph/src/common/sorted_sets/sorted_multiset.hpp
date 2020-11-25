#ifndef __SORTED_MULTISET_HPP__
#define __SORTED_MULTISET_HPP__

#include <cassert>
#include <mutex>
#include <shared_mutex>

#include "common/logger.hpp"
#include "common/vector.hpp"


namespace mtg {
namespace common {

// Thread safe data storage for counting
template <typename T,
          typename C = uint8_t,
          class Container = Vector<std::pair<T, C>>>
class SortedMultiset {
  public:
    static_assert(std::is_same_v<std::pair<T, C>, typename Container::value_type>);

    typedef T key_type;
    typedef C count_type;
    typedef std::pair<T, C> value_type;
    typedef Container storage_type;
    typedef Container result_type;

    SortedMultiset(size_t num_threads = 1, size_t max_num_elements = 0)
          : num_threads_(num_threads) {
        reserve(max_num_elements);
    }

    static constexpr uint64_t max_count() { return std::numeric_limits<C>::max(); }

    template <class Iterator>
    inline void insert(Iterator begin, Iterator end);

    void reserve(size_t size);

    size_t buffer_size() const { return data_.capacity(); }

    result_type& data();

    void clear();

  private:
    void shrink_data();

    void sort_and_merge_counts();

    void try_reserve(size_t size, size_t min_size = 0);

    storage_type data_;
    size_t num_threads_;

    // indicate the end of the preprocessed distinct and sorted values
    uint64_t sorted_end_ = 0;

    mutable std::mutex mutex_resize_;
    mutable std::shared_timed_mutex mutex_copy_;
};

template <typename T, typename C, class Container>
template <class Iterator>
void SortedMultiset<T, C, Container>::insert(Iterator begin, Iterator end) {
    assert(begin <= end);

    uint64_t batch_size = end - begin;

    if (!batch_size)
        return;

    // acquire the mutex to restrict the number of writing threads
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);

    if (data_.size() + batch_size > data_.capacity()) {
        std::unique_lock<std::shared_timed_mutex> reallocate_lock(mutex_copy_);

        if (data_.size() != sorted_end_)
            shrink_data();

        try {
            try_reserve(data_.size() + data_.size() / 2,
                        data_.size() + batch_size);
        } catch (const std::bad_alloc &exception) {
            logger->error("Can't reallocate. Not enough memory");
            exit(1);
        }
    }

    size_t offset = data_.size();
    data_.resize(data_.size() + batch_size);

    std::shared_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    resize_lock.unlock();

    if constexpr(std::is_same_v<T, std::decay_t<decltype(*begin)>>) {
        std::transform(begin, end, data_.begin() + offset,
            [](const T &value) { return std::make_pair(value, C(1)); });
    } else {
        std::copy(begin, end, data_.begin() + offset);
    }
}

} // namespace common
} // namespace mtg

#endif // __SORTED_MULTISET_HPP__
