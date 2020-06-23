#include "sorted_set.hpp"

#include <ips4o.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace common {

template <typename T, class Container>
void SortedSet<T, Container>::reserve(size_t size) {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    try_reserve(size);
}

template <typename T, class Container>
typename SortedSet<T, Container>::result_type& SortedSet<T, Container>::data() {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    if (sorted_end_ != data_.size()) {
        sort_and_remove_duplicates();
        sorted_end_ = data_.size();
    }

    return data_;
}

template <typename T, class Container>
void SortedSet<T, Container>::clear() {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    data_ = storage_type();
    sorted_end_ = 0;
}

template <typename T, class Container>
void SortedSet<T, Container>::sort_and_remove_duplicates() {
    ips4o::parallel::sort(data_.begin(), data_.end(),
                          std::less<value_type>(), num_threads_);
    // remove duplicates
    auto unique_end = std::unique(data_.begin(), data_.end());
    data_.erase(unique_end, data_.end());
}

template <typename T, class Container>
void SortedSet<T, Container>::shrink_data() {
    logger->trace("Allocated capacity exceeded, erase duplicate values...");

    size_t old_size = data_.size();
    sort_and_remove_duplicates();
    sorted_end_ = data_.size();

    logger->trace("Erasing duplicate values done. Size reduced from {} to {}, {} MiB",
                  old_size, data_.size(), (data_.size() * sizeof(T) >> 20));
}

template <typename T, class Container>
void SortedSet<T, Container>::try_reserve(size_t size, size_t min_size) {
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

template class SortedSet<uint64_t>;
template class SortedSet<sdsl::uint128_t>;
template class SortedSet<sdsl::uint256_t>;

} // namespace common
} // namespace mtg
