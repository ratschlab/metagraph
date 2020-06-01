#include "sorted_multiset_disk.hpp"

#include <cassert>

#include <ips4o.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace common {

template <typename T, typename C>
void SortedMultisetDisk<T, C>::sort_and_dedupe() {
    ips4o::parallel::sort(
            this->data_.begin(), this->data_.end(),
            [](const value_type &first, const value_type &second) {
                return first.first < second.first;
            },
            this->num_threads_);

    auto first = this->data_.begin();
    auto last = this->data_.end();

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

    this->data_.erase(++dest, this->data_.end());
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
} // namespace mtg
