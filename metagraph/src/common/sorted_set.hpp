#ifndef __SORTED_SET_HPP__
#define __SORTED_SET_HPP__

#include <cassert>
#include <iostream>
#include <mutex>
#include <shared_mutex>
#include <vector>

#include <ips4o.hpp>

#include "common/logger.hpp"
#include "common/vector.hpp"

namespace mg {
namespace common {

// Thread safe data storage to extract distinct elements
template <typename T, class Container = Vector<T>>
class SortedSet {
  public:
    static_assert(std::is_same_v<T, typename Container::value_type>);

    typedef T key_type;
    typedef T value_type;
    typedef Container storage_type;
    typedef Container result_type;

    SortedSet(std::function<void(storage_type *)> cleanup = [](storage_type *) {},
              size_t num_threads = 1,
              size_t reserved_num_elements = 0)
          : num_threads_(num_threads), cleanup_(cleanup) {
        reserve(reserved_num_elements);
    }

    ~SortedSet() {}

    template <class Iterator>
    void insert(Iterator begin, Iterator end) {
        assert(begin <= end);

        uint64_t batch_size = end - begin;

        // acquire the mutex to restrict the number of writing threads
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);

        if (data_.size() + batch_size > data_.capacity()) {
            std::unique_lock<std::shared_timed_mutex> reallocate_lock(mutex_copy_);

            shrink_data();

            try {
                try_reserve(data_.size() + data_.size() / 2,
                            data_.size() + batch_size);
            } catch (const std::bad_alloc &exception) {
                logger->error("ERROR: Can't reallocate. Not enough memory");
                exit(1);
            }
        }

        size_t offset = data_.size();
        data_.resize(data_.size() + batch_size);

        std::shared_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        resize_lock.unlock();

        std::copy(begin, end, data_.begin() + offset);
    }

    void reserve(size_t size) {
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);
        std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        try_reserve(size);
    }

    size_t buffer_size() const { return data_.capacity(); }

    result_type& data() {
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);
        std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        if (sorted_end_ != data_.size()) {
            sort_and_remove_duplicates(&data_, num_threads_);
            sorted_end_ = data_.size();
        }

        return data_;
    }

    void clear() {
        std::unique_lock<std::mutex> resize_lock(mutex_resize_);
        std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

        data_ = storage_type();
        sorted_end_ = 0;
    }

    template <class Array>
    void sort_and_remove_duplicates(Array *vector, size_t num_threads) const {
        assert(vector);

        ips4o::parallel::sort(vector->begin(), vector->end(),
                              std::less<typename Array::value_type>(),
                              num_threads);
        // remove duplicates
        auto unique_end = std::unique(vector->begin(), vector->end());
        vector->erase(unique_end, vector->end());

        cleanup_(vector);
    }

  private:
    void shrink_data() {
        logger->trace("Allocated capacity exceeded, erase duplicate values...");

        size_t old_size = data_.size();
        sort_and_remove_duplicates(&data_, num_threads_);
        sorted_end_ = data_.size();

        logger->trace("Erasing duplicate values done. Size reduced from {} to {}, {}MiB",
                      old_size, data_.size(), (data_.size() * sizeof(T) >> 20));
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

    std::function<void(storage_type *)> cleanup_;

    // indicate the end of the preprocessed distinct and sorted values
    uint64_t sorted_end_ = 0;

    mutable std::mutex mutex_resize_;
    mutable std::shared_timed_mutex mutex_copy_;
};

} // namespace common
} // namespace mg

#endif // __SORTED_SET_HPP__
