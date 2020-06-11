#include "sorted_multiset.hpp"

#include <ips4o.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace common {

template <typename T, typename C, class Container>
void SortedMultiset<T, C, Container>::reserve(size_t size) {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    try_reserve(size);
}

template <typename T, typename C, class Container>
typename SortedMultiset<T, C, Container>::result_type&
SortedMultiset<T, C, Container>::data() {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    if (sorted_end_ != data_.size()) {
        sort_and_merge_counts();
        sorted_end_ = data_.size();
    }

    return data_;
}

template <typename T, typename C, class Container>
void SortedMultiset<T, C, Container>::clear() {
    std::unique_lock<std::mutex> resize_lock(mutex_resize_);
    std::unique_lock<std::shared_timed_mutex> copy_lock(mutex_copy_);

    data_ = decltype(data_)();
    sorted_end_ = 0;
}

template <typename T, typename C, class Container>
void SortedMultiset<T, C, Container>::shrink_data() {
    logger->trace("Allocated capacity exceeded, erasing duplicate values...");

    size_t old_size = data_.size();
    sort_and_merge_counts();
    sorted_end_ = data_.size();

    logger->trace("...done. Size reduced from {} to {}, {} MiB", old_size,
                  data_.size(), (data_.size() * sizeof(value_type) >> 20));
}

template <class Container>
void sort_and_merge(size_t num_threads, Container *data) {
    if (data->empty())
        return;

    ips4o::parallel::sort(data->begin(), data->end(),
        [](const auto &first, const auto &second) {
            return first.first < second.first;
        },
        num_threads);

    auto first = data->begin();
    auto last = data->end();

    auto dest = first;

    using T = typename Container::value_type;
    constexpr auto max = std::numeric_limits<typename T::second_type>::max();
    while (++first != last) {
        if (first->first == dest->first) {
            if (first->second < max - dest->second) {
                dest->second += first->second;
            } else {
                dest->second = max;
            }
        } else {
            *++dest = std::move(*first);
        }
    }

    data->erase(++dest, data->end());
}

template <typename T, typename C, class Container>
void SortedMultiset<T, C, Container>::sort_and_merge_counts() {
  sort_and_merge(num_threads_, &data_);
}

template <typename T, typename C, class Container>
void SortedMultiset<T, C, Container>::try_reserve(size_t size, size_t min_size) {
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

template class SortedMultiset<uint64_t, uint8_t>;
template class SortedMultiset<sdsl::uint128_t, uint8_t>;
template class SortedMultiset<sdsl::uint256_t, uint8_t>;
template class SortedMultiset<uint64_t, uint16_t>;
template class SortedMultiset<sdsl::uint128_t, uint16_t>;
template class SortedMultiset<sdsl::uint256_t, uint16_t>;
template class SortedMultiset<uint64_t, uint32_t>;
template class SortedMultiset<sdsl::uint128_t, uint32_t>;
template class SortedMultiset<sdsl::uint256_t, uint32_t>;

template void sort_and_merge(size_t, Vector<std::pair<uint64_t, uint8_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<sdsl::uint128_t, uint8_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<sdsl::uint256_t, uint8_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<uint64_t, uint16_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<sdsl::uint128_t, uint16_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<sdsl::uint256_t, uint16_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<uint64_t, uint32_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<sdsl::uint128_t, uint32_t>> *);
template void sort_and_merge(size_t, Vector<std::pair<sdsl::uint256_t, uint32_t>> *);
} // namespace common
} // namespace mtg
