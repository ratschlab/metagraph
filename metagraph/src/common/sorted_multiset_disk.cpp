#include "sorted_multiset_disk.hpp"

#include <cassert>

#include <ips4o.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mg {
namespace common {

template <typename T, typename C>
void SortedMultisetDisk<T, C>::sort_and_remove_duplicates(storage_type *vector,
                                                    size_t num_threads) const {
    assert(vector);
    ips4o::parallel::sort(
            vector->begin(), vector->end(),
            [](const value_type &first, const value_type &second) {
                return first.first < second.first;
            },
            num_threads);

    auto first = vector->begin();
    auto last = vector->end();

    auto dest = first;

    while (++first != last) {
        if (first->first == dest->first) {
            if (first->second < max_count() - dest->second) {
                dest->second += first->second;
            } else {
                dest->second = max_count();
            }
        } else {
            *++dest = std::move(*first);
        }
    }

    vector->erase(++dest, this->data_.end());

    this->cleanup_(vector);
}

template class SortedMultisetDisk<uint64_t, uint8_t>;
template class SortedMultisetDisk<sdsl::uint128_t, uint8_t>;
template class SortedMultisetDisk<sdsl::uint256_t, uint8_t>;
template class SortedMultisetDisk<uint64_t, uint16_t>;
template class SortedMultisetDisk<sdsl::uint128_t, uint16_t>;
template class SortedMultisetDisk<sdsl::uint256_t, uint16_t>;
template class SortedMultisetDisk<uint64_t, uint32_t>;
template class SortedMultisetDisk<sdsl::uint128_t, uint32_t>;
template class SortedMultisetDisk<sdsl::uint256_t, uint32_t>;

} // namespace common
} // namespace mg
