#include "test_kmer_helpers.hpp"

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


template <typename TypeParam>
void left_shift(TypeParam *num, size_t bits) {
    (*num) = num->operator<<(bits);
}

template void left_shift<sdsl::uint128_t>(sdsl::uint128_t*, size_t);
template void left_shift<sdsl::uint256_t>(sdsl::uint256_t*, size_t);

template<>
void left_shift(uint64_t *num, size_t bits) {
    (*num) = (*num) << bits;
}
