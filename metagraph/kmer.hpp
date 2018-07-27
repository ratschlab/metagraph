#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

typedef uint64_t KMerCharType;

#if _PROTEIN_GRAPH
const int kBitsPerChar = 5;
#elif _DNA_CASE_SENSITIVE_GRAPH
const int kBitsPerChar = 4;
#elif _DNA_GRAPH
const int kBitsPerChar = 3;
#else
#endif


template <typename G>
class KMer {
    template <typename U>
    friend std::ostream& operator<<(std::ostream &os, const KMer<U> &kmer);

  public:
    typedef G KMerWordType;

    KMer() {}
    template <typename T>
    KMer(const T &arr, size_t k) : seq_(pack_kmer(arr, k)) {}

    template <class Map, class String>
    KMer(const String &seq, Map &&to_alphabet);

    explicit KMer(KMerWordType &&seq) : seq_(seq) {}
    explicit KMer(const KMerWordType &seq) : seq_(seq) {}

    bool operator<(const KMer &other) const { return seq_ < other.seq_; }
    bool operator==(const KMer &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMer &other) const { return seq_ != other.seq_; }

    inline KMerCharType operator[](size_t i) const;

    static bool compare_suffix(const KMer &k1,
                               const KMer &k2, size_t minus = 0);

    std::string to_string(const std::string &alphabet) const;

    template<typename T>
    static inline KMerWordType pack_kmer(const T &arr, size_t k);

    /**
     * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
     * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
     *      = s[7] << k + (kmer & mask) >> 1 + s[8].
     */
    static inline  void update_kmer(size_t k,
                                    KMerCharType edge_label,
                                    KMerCharType last,
                                    KMerWordType *kmer);
  private:
    inline KMerCharType get_digit(size_t i) const;

    static const KMerCharType kFirstCharMask;

    KMerWordType seq_; // kmer sequence
};


template <typename G>
KMerCharType KMer<G>::operator[](size_t i) const {
    assert(get_digit(i) > 0);
    return get_digit(i) - 1;
}

template <typename G>
template <typename T>
G KMer<G>::pack_kmer(const T &arr, size_t k) {
    if (k * kBitsPerChar > sizeof(KMerWordType) * 8 || k < 2) {
        std::cerr << "ERROR: Too large k-mer size: must be between 2 and "
                  << sizeof(KMerWordType) * 8 / kBitsPerChar << std::endl;
        exit(1);
    }

    KMerWordType result(0);

    for (int i = k - 2; i >= 0; --i) {

        assert(static_cast<uint64_t>(arr[i] + 1) < (1llu << kBitsPerChar)
                 && "Alphabet size too big for the given number of bits");

        result = result << kBitsPerChar;
        result += arr[i] + 1;
    }
    result = result << kBitsPerChar;
    result += arr[k - 1] + 1;
    return result;
}

template <typename G>
template <class Map, class String>
KMer<G>::KMer(const String &seq, Map &&to_alphabet) {
    std::vector<uint8_t> arr(seq.size());
    std::transform(seq.begin(), seq.end(), arr.begin(), to_alphabet);
    *this = KMer(arr.data(), arr.size());
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G>
void KMer<G>::update_kmer(size_t k,
                          KMerCharType edge_label,
                          KMerCharType last,
                          KMerWordType *kmer) {
    *kmer = *kmer >> kBitsPerChar;
    *kmer += KMerWordType(last + 1) << static_cast<int>(kBitsPerChar * k);
    *kmer |= kFirstCharMask;
    *kmer -= kFirstCharMask;
    *kmer += edge_label + 1;
}

template <typename G>
KMerCharType KMer<G>::get_digit(size_t i) const {
    static_assert(kBitsPerChar <= 64, "too large digit");
    assert(kBitsPerChar * (i + 1) <= sizeof(KMerWordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}


#endif // __KMER_HPP__
