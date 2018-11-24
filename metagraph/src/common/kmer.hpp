#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>


template <typename G, int L>
class KMer {
    template <typename U, int P>
    friend std::ostream& operator<<(std::ostream &os, const KMer<U, P> &kmer);

  public:
    typedef G KMerWordType;
    typedef uint64_t KMerCharType;
    static const int kBitsPerChar = L;

    KMer() {}
    template <typename V>
    KMer(const V &arr, size_t k) : seq_(pack_kmer(arr, k)) {}
    template <typename T>
    KMer(const std::vector<T> &arr) : KMer(arr, arr.size()) {}

    explicit KMer(KMerWordType &&seq) : seq_(seq) {}
    explicit KMer(const KMerWordType &seq) : seq_(seq) {}

    bool operator<(const KMer &other) const { return seq_ < other.seq_; }
    bool operator==(const KMer &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMer &other) const { return seq_ != other.seq_; }

    inline KMerCharType operator[](size_t i) const;

    static inline bool compare_suffix(const KMer &k1,
                                      const KMer &k2, size_t minus = 0);

    std::string to_string(const std::string &alphabet) const;

    template <typename T>
    static inline KMerWordType pack_kmer(const T &arr, size_t k);

    /**
     * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
     * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
     *      = s[7] << k + (kmer & mask) >> 1 + s[8].
     */
    static inline void update_kmer(size_t k,
                                   KMerCharType edge_label,
                                   KMerCharType last,
                                   KMerWordType *kmer);

    inline KMer<G, L> prev_kmer(size_t k, KMerCharType first_char) const;
    inline void next_kmer(size_t k, KMerCharType edge_label);

  private:
    inline KMerCharType get_digit(size_t i) const;

    static const KMerCharType kFirstCharMask;

    KMerWordType seq_; // kmer sequence
};


template <typename G, int L>
typename KMer<G, L>::KMerCharType KMer<G, L>::operator[](size_t i) const {
    assert(get_digit(i) > 0);
    return get_digit(i) - 1;
}

template <typename G, int L>
template <typename T>
G KMer<G, L>::pack_kmer(const T &arr, size_t k) {
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

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMer<G, L>::update_kmer(size_t k,
                             KMerCharType edge_label,
                             KMerCharType last,
                             KMerWordType *kmer) {
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    *kmer = *kmer >> kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    *kmer += KMerWordType(last + 1) << static_cast<int>(kBitsPerChar * (k - 1));
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    *kmer |= kFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]1111
    *kmer -= kFirstCharMask - (edge_label + 1);
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMer<G, L>::next_kmer(size_t k, KMerCharType edge_label) {
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    KMerWordType last = seq_ & KMerWordType((1llu << kBitsPerChar) - 1);
    seq_ = seq_ >> kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ += last << static_cast<int>(kBitsPerChar * k);
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ |= kFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]1111
    seq_ -= kFirstCharMask - (edge_label + 1);
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
}

template <typename G, int L>
KMer<G, L> KMer<G, L>::prev_kmer(size_t k, KMerCharType first_char) const {
    KMerWordType kmer = seq_;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
    kmer |= kFirstCharMask;
    kmer -= kFirstCharMask - (first_char + 1);
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    int shift = kBitsPerChar * (k - 1);
    KMerWordType last_char = seq_ >> shift;
    kmer -= last_char << shift;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    kmer = kmer << kBitsPerChar;
    // s[6]s[5]s[4]s[3]s[2]s[1]0000
    kmer += last_char;
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    return KMer<G, L>(kmer);
}

template <typename G, int L>
typename KMer<G, L>::KMerCharType KMer<G, L>::get_digit(size_t i) const {
    static_assert(kBitsPerChar <= 64, "too large digit");
    assert(kBitsPerChar * (i + 1) <= sizeof(KMerWordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}

template <typename G, int L>
bool KMer<G, L>::compare_suffix(const KMer &k1, const KMer &k2, size_t minus) {
    return k1.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar)
             == k2.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar);
}

#endif // __KMER_HPP__
