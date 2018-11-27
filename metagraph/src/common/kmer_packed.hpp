#ifndef __KMER_PACKED_HPP__
#define __KMER_PACKED_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>


template <typename G, int L>
class KMerPacked {
    template <typename U, int P>
    friend std::ostream& operator<<(std::ostream &os, const KMerPacked<U, P> &kmer);

  public:
    typedef G KMerWordType;
    typedef uint64_t KMerCharType;
    static const int kBitsPerChar = L;

    KMerPacked() {}
    template <typename V>
    KMerPacked(const V &arr, size_t k) : seq_(pack_kmer(arr, k)) {}
    template <typename T>
    KMerPacked(const std::vector<T> &arr) : KMerPacked(arr, arr.size()) {}

    explicit KMerPacked(KMerWordType &&seq) : seq_(seq) {}
    explicit KMerPacked(const KMerWordType &seq) : seq_(seq) {}

    bool operator<(const KMerPacked &other) const { return seq_ < other.seq_; }
    bool operator==(const KMerPacked &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMerPacked &other) const { return seq_ != other.seq_; }

    inline KMerCharType operator[](size_t i) const;

    static inline bool compare_suffix(const KMerPacked &k1,
                                      const KMerPacked &k2, size_t minus = 0);

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

    inline KMerPacked<G, L> prev_kmer(size_t k, KMerCharType first_char) const;
    inline void next_kmer(size_t k, KMerCharType edge_label);

    inline KMerWordType& data() { return seq_; }
    inline KMerWordType data() const { return seq_; }

    inline uint32_t get_k() const { return get_k(seq_); }

  private:
    static void set_length_bit(KMerWordType *kmer, uint64_t k);
    static void unset_length_bit(KMerWordType *kmer, uint64_t k);
    static uint32_t get_k(const KMerWordType &kmer);
    inline KMerCharType get_digit(size_t i) const;

    static const KMerCharType kFirstCharMask;

    KMerWordType seq_; // kmer sequence
};


template <typename G, int L>
typename KMerPacked<G, L>::KMerCharType KMerPacked<G, L>::operator[](size_t i) const {
    return get_digit(i);
}

template <typename G, int L>
template <typename T>
G KMerPacked<G, L>::pack_kmer(const T &arr, size_t k) {
    if (k + 1 > sizeof(KMerWordType) * 8 / kBitsPerChar || k < 2) {
        std::cerr << "ERROR: Too large k-mer size: must be between 2 and "
                  << (sizeof(KMerWordType) * 8 / kBitsPerChar - 1) << std::endl;
        exit(1);
    }

    KMerWordType result(1);

    for (int i = k - 2; i >= 0; --i) {

        assert(static_cast<uint64_t>(arr[i]) < (1llu << kBitsPerChar)
                 && "Alphabet size too big for the given number of bits");

        result = result << kBitsPerChar;
        result |= arr[i];
    }
    result = result << kBitsPerChar;
    result |= arr[k - 1];
    assert(get_k(result) == k);
    return result;
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMerPacked<G, L>::update_kmer(size_t k,
                             KMerCharType edge_label,
                             KMerCharType last,
                             KMerWordType *kmer) {
    assert(k == get_k(*kmer));
    assert(KMerWordType(last) == (*kmer & KMerWordType(kFirstCharMask)));
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    unset_length_bit(kmer, k);
    *kmer = *kmer >> kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    *kmer += KMerWordType(last) << static_cast<int>(kBitsPerChar * (k - 1));
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    *kmer |= kFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]1111
    *kmer -= kFirstCharMask - edge_label;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
    set_length_bit(kmer, k);
    assert(k == get_k(*kmer));
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMerPacked<G, L>::next_kmer(size_t k, KMerCharType edge_label) {
    assert(k == get_k(seq_));
    unset_length_bit(&seq_, k);
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    KMerWordType last = seq_ & KMerWordType((1llu << kBitsPerChar) - 1);
    seq_ = seq_ >> kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ += last << static_cast<int>(kBitsPerChar * k);
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ |= kFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]1111
    seq_ -= kFirstCharMask - edge_label;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
    set_length_bit(&seq_, k);
    assert(k == get_k(seq_));
}

template <typename G, int L>
KMerPacked<G, L> KMerPacked<G, L>::prev_kmer(size_t k, KMerCharType first_char) const {
    assert(k == get_k(seq_));
    KMerWordType kmer = seq_;
    unset_length_bit(&kmer, k);
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
    kmer |= kFirstCharMask;
    kmer -= kFirstCharMask - first_char;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    KMerWordType last_char = seq_ >> static_cast<int>(kBitsPerChar * k);
    kmer -= last_char << static_cast<int>(kBitsPerChar * k);
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    kmer = kmer << kBitsPerChar;
    // s[6]s[5]s[4]s[3]s[2]s[1]0000
    kmer += last_char;
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    set_length_bit(&kmer, k);
    assert(k == get_k(kmer));
    return KMerPacked<G, L>(kmer);
}

template <typename G, int L>
typename KMerPacked<G, L>::KMerCharType KMerPacked<G, L>::get_digit(size_t i) const {
    static_assert(kBitsPerChar <= 64, "too large digit");
    assert(kBitsPerChar * (i + 1) <= sizeof(KMerWordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}

template <typename G, int L>
bool KMerPacked<G, L>::compare_suffix(const KMerPacked &k1, const KMerPacked &k2, size_t minus) {
    return k1.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar)
             == k2.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar);
}

#endif // __KMER_PACKED_HPP__
