#include "sorted_set_disk.hpp"

#include <cassert>

#include <ips4o.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace common {

template <typename T>
void SortedSetDisk<T>::sort_and_dedupe() {
    ips4o::parallel::sort(this->data_.begin(), this->data_.end(),
                          std::less<value_type>(), this->num_threads_);
    // remove duplicates
    this->data_.erase(std::unique(this->data_.begin(), this->data_.end()),
                      this->data_.end());
}

template class SortedSetDisk<uint64_t>;
template class SortedSetDisk<sdsl::uint128_t>;
template class SortedSetDisk<sdsl::uint256_t>;

} // namespace common
} // namespace mtg
