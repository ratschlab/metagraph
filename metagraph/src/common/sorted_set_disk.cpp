#include "sorted_set_disk.hpp"

#include <cassert>

#include <ips4o.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mg {
namespace common {

template <typename T>
void SortedSetDisk<T>::sort_and_remove_duplicates(storage_type *vector,
                                                  size_t num_threads) const {
    assert(vector);

    ips4o::parallel::sort(vector->begin(), vector->end(), std::less<value_type>(),
                          num_threads);
    // remove duplicates
    auto unique_end = std::unique(vector->begin(), vector->end());
    vector->erase(unique_end, vector->end());

    this->cleanup_(vector); // typically removes source dummy k-mers
}

template class SortedSetDisk<uint64_t>;
template class SortedSetDisk<sdsl::uint128_t>;
template class SortedSetDisk<sdsl::uint256_t>;

} // namespace common
} // namespace mg
