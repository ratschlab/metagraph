#pragma once

#include <cstddef>

/** Compile-time utils for conversion between k-mers and integers. */

namespace mtg {
namespace kmer {

template <typename T, typename = void>
struct get_int { using type = typename T::WordType; };

template <typename T>
using get_int_t = typename get_int<T>::type;

template <typename T>
struct get_int<T, std::void_t<typename T::first_type>> {
    using type = std::pair<get_int_t<typename T::first_type>, typename T::second_type>;
};


template <typename KMER, typename T, typename = void>
struct get_kmer {
    static_assert(std::is_same_v<T, typename KMER::WordType>);
    using type = KMER;
};

template <typename KMER, typename T>
using get_kmer_t = typename get_kmer<KMER, T>::type;

template <typename KMER, typename T>
struct get_kmer<KMER, T, std::void_t<typename T::second_type>> {
    static_assert(std::is_same_v<typename T::first_type, typename KMER::WordType>);
    using type = std::pair<KMER, typename T::second_type>;
};

} // namespace kmer
} // namespace mtg
