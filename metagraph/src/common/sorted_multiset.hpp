#ifndef __SORTED_MULTISET_HPP__
#define __SORTED_MULTISET_HPP__

#include <mutex>
#include <shared_mutex>
#include <iostream>
#include <vector>
#include <cassert>

#include <ips4o.hpp>

#include "logger.hpp"

namespace mg {
namespace common {

// Thread safe data storage for counting
template <typename T,
          typename C = uint8_t,
          class Container = std::vector<std::pair<T, C>>>
class SortedMultiset {
  public:
    static_assert(std::is_same_v<std::pair<T, C>, typename Container::value_type>);

    typedef T key_type;
    typedef C count_type;
    typedef std::pair<T, C> value_type;
    typedef Container storage_type;
    typedef Container result_type;

    SortedMultiset(std::function<void(storage_type*)> cleanup = [](storage_type*) {},
                   size_t num_threads = 1, size_t max_num_elements = 0)
      : num_threads_(num_threads), cleanup_(cleanup) {
        reserve(max_num_elements);
    }

    ~SortedMultiset() {}

    static constexpr uint64_t max_count() { return std::numeric_limits<C>::max(); }

    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        assert(begin <= end);

        uint64_t batch_size = end - begin;

        if (!batch_size)
            return;

        // acquire the mutex to restrict the number of writing threads
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);

        if (data_.size() + batch_size > data_.capacity()) {
            std::unique_lock<std::shared_timed_mutex> reallocate_lock(mutex_copy_);

            shrink_data();

            try {
                try_reserve(data_.size() + data_.size() / 2,
                            data_.size() + batch_size);
            } catch (const std::bad_alloc &exception) {
                std::cerr << "ERROR: Can't reallocate. Not enough memory" << std::endl;
                exit(1);
            }
        }

        size_t offset = data_.size();
        data_.resize(data_.size() + batch_size);

        std::shared_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        resize_lock.unlock();

        if constexpr(std::is_same<T, std::remove_cv_t<
                                     std::remove_reference_t<
                                        decltype(*begin)>>>::value) {
            std::transform(begin, end, data_.begin() + offset,
                [](const T &value) { return std::make_pair(value, C(1)); });
        } else {
            std::copy(begin, end, data_.begin() + offset);
        }
    }

    void reserve(size_t size) {
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);
        std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        try_reserve(size);
    }

    size_t buffer_size() {
        return data_.capacity();
    }

    result_type& data() {
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);
        std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        if (sorted_end_ != data_.size()) {
            sort_and_merge_duplicates();
            sorted_end_ = data_.size();
        }

        return data_;
    }

    void clear() {
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);
        std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        data_ = decltype(data_)();
        sorted_end_ = 0;
    }

  private:
    void shrink_data() {
        logger->trace("Allocated capacity exceeded, erasing duplicate values...");

        size_t old_size = data_.size();
        sort_and_merge_duplicates();
        sorted_end_ = data_.size();

        logger->trace("...done. Size reduced from {} to {}, {}MiB", old_size,
                      data_.size(), (data_.size() * sizeof(value_type) >> 20));
    }

    void sort_and_merge_duplicates() {
        if (!data_.size())
            return;

        ips4o::parallel::sort(data_.begin(), data_.end(),
            [](const value_type &first, const value_type &second) {
                return first.first < second.first;
            },
            num_threads_
        );

        auto first = data_.begin();
        auto last = data_.end();

        auto dest = first;

        while (++first != last) {
            if (first->first == dest->first) {
                if (first->second < max_count() - dest->second) {
                    dest->second += first->second;
                } else {
                    dest->second = max_count();
                }
            } else {
                *++dest = std::move(*first);;
            }
        }

        data_.erase(++dest, data_.end());

        cleanup_(&data_);
    }

    void try_reserve(size_t size, size_t min_size = 0) {
        size = std::max(size, min_size);

        while (size > min_size) {
            try {
                data_.reserve(size);
                return;
            } catch (const std::bad_alloc &exception) {
                size = min_size + (size - min_size) * 2 / 3;
            }
        }
        data_.reserve(min_size);
    }

    storage_type data_;
    size_t num_threads_;

    std::function<void(storage_type*)> cleanup_;

    // indicate the end of the preprocessed distinct and sorted values
    uint64_t sorted_end_ = 0;

    mutable std::mutex mutex_resize_;
    mutable std::shared_timed_mutex mutex_copy_;
};

} // namespace common
} // namespace mg

#endif // __SORTED_MULTISET_HPP__
