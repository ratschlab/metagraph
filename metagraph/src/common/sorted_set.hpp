#ifndef __SORTED_SET_HPP__
#define __SORTED_SET_HPP__

#include <cassert>
#include <mutex>
#include <shared_mutex>

#include "common/logger.hpp"
#include "common/vector.hpp"


namespace mtg {
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

    SortedSet(size_t num_threads = 1, size_t reserved_num_elements = 0)
          : num_threads_(num_threads) {
        reserve(reserved_num_elements);
    }

    template <class Iterator>
    inline void insert(Iterator begin, Iterator end);

    void reserve(size_t size);

    size_t buffer_size() const { return data_.capacity(); }

    result_type& data();

    void clear();

  private:
    void shrink_data();

    void sort_and_remove_duplicates();

    void try_reserve(size_t size, size_t min_size = 0);

    storage_type data_;
    size_t num_threads_;

    // indicate the end of the preprocessed distinct and sorted values
    uint64_t sorted_end_ = 0;

    mutable std::mutex mutex_resize_;
    mutable std::shared_timed_mutex mutex_copy_;
};

template <typename T, class Container>
template <class Iterator>
void SortedSet<T, Container>::insert(Iterator begin, Iterator end) {
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
            logger->error("Can't reallocate. Not enough memory");
            exit(1);
        }
    }

    size_t offset = data_.size();
    data_.resize(data_.size() + batch_size);

    std::shared_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    resize_lock.unlock();

    std::copy(begin, end, data_.begin() + offset);
}

} // namespace common
} // namespace mtg

#endif // __SORTED_SET_HPP__
