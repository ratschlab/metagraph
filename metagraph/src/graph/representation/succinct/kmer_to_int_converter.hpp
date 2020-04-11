#pragma once
#include <cstddef>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "kmer/kmer_extractor.hpp"

// Utility methods for converting between a KMERBoss and its integer representation
// None of the methods below should generate any code, so they shouldn't add any overhead

// Convert from a 64/128/256-bit integer, to the corresponding KMERBoss type
template <typename T, typename = void>
struct get_kmer {};
template <>
struct get_kmer<uint64_t> { using type = KmerExtractorBOSS::Kmer64;};
template <>
struct get_kmer<sdsl::uint128_t> { using type = KmerExtractorBOSS::Kmer128;};
template <>
struct get_kmer<sdsl::uint256_t> { using type = KmerExtractorBOSS::Kmer256;};
template <typename T>
using get_kmer_t = typename get_kmer<T>::type;
template <typename T>
struct get_kmer<T, void_t<typename T::first_type>> {
using type = std::pair<get_kmer_t<typename T::first_type>, typename T::second_type>;
};

/** Converts an integer value to its corresponding KMERBoss value */
template <typename KMER_INT>
inline const get_kmer_t<KMER_INT>& to_kmer(const KMER_INT &kmer_int) {
    return *reinterpret_cast<const get_kmer_t<KMER_INT> *>(&kmer_int);
}

template <typename KMER_INT>
inline get_kmer_t<KMER_INT>& to_kmer(KMER_INT &kmer_int) {
    return *reinterpret_cast<get_kmer_t<KMER_INT> *>(&kmer_int);
}
