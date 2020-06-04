#ifndef __KMER_HPP__
#define __KMER_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace kmer {

/**
 * Models a kmer (https://en.wikipedia.org/wiki/K-mer). Each character in the k-mer
 * uses L bits of the internal representation type G (typically a 64, 128 or 256 bit
 * integer). The last character of the k-mer uses the most significant bits of
 * G, while the first will use the least significant ones. In other words, the memory
 * layout of the k-mer "ACGT" is "TGCA". This way, ordering k-mers in co-lex order can be
 * achieved by simply comparing the internal representations in G.
 *
 * @tparam G the type storing a kmer, typically a 64/128/256 bit integer
 * @tparam L the number of bits storing a k-mer character (2 bits for DNA)
 */
template <typename G, int L>
class KMer {
  public:
    typedef G WordType;
    typedef uint64_t CharType;
    static constexpr int kBitsPerChar = L;

    /** Construct a default (uninitialized) k-mer. */
    KMer() {}

    /**
     * Construct a k-mer with the given size from the given array
     * @tparam V an indexed data structure (e.g. std::vector or L[])
     * @param arr unpacked k-mer
     * @param k k-mer length
     */
    template <typename V>
    KMer(const V &arr, size_t k);
    /**
     * Construct a k-mer from the given vector
     * @tparam T k-mer character type
     * @param arr unpacked k-mer
     */
    template <typename T>
    KMer(const std::vector<T> &arr) : KMer(arr, arr.size()) {}

    KMer(WordType&& seq) noexcept : seq_(seq) {}
    explicit KMer(const WordType &seq) noexcept : seq_(seq) {}

    // corresponds to the co-lex order of k-mers
    bool operator<(const KMer &other) const { return seq_ < other.seq_; }
    bool operator<=(const KMer &other) const { return seq_ <= other.seq_; }
    bool operator>(const KMer &other) const { return seq_ > other.seq_; }
    bool operator>=(const KMer &other) const { return seq_ >= other.seq_; }
    bool operator==(const KMer &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMer &other) const { return seq_ != other.seq_; }

    /** Return the character at position #i in the kmer. Undefined behavior if #i is
     * out of range. */
    inline CharType operator[](size_t i) const;

    /**
     * Return the human-readable representation of the kmer using the given alphabet.
     */
    std::string to_string(size_t k, const std::string &alphabet) const;

    /**
     * Construct the next k-mer for s[7]s[6]s[5]s[4]s[3]s[2]s[1].
     * next = s[8]s[7]s[6]s[5]s[4]s[3]s[2]
     *      = ( s[8] << k ) + ( kmer >> 1 ).
     * @param k the k-mer size
     * @param edge_label the new last character
     */
    inline void to_next(size_t k, WordType edge_label);
    /**
     * Replaces the current k-mer with its predecessor. The last k-1 characters of the
     * predecessors are identical to the first k-1 characters of the successors.
     * @param k the k-mer size
     * @param first_char the first character in the predecessor
     */
    inline void to_prev(size_t k, CharType first_char);

    inline const WordType& data() const { return seq_; }

    // TODO: remove this from here
    template <typename T>
    inline static bool match_suffix(const T *kmer, size_t k, const std::vector<T> &suffix) {
        assert(k > 0);
        assert(k >= suffix.size());
        return suffix.empty()
                || std::equal(suffix.begin(), suffix.end(), kmer + k - suffix.size());
    }

    void print_hex(std::ostream &os) const;

  private:
    /** Bit mask for extracting the first character in packed kmers. */
    static constexpr CharType kFirstCharMask = (1ull << kBitsPerChar) - 1;
    static inline const WordType kAllSetMask = ~(WordType(0ull));
    /** Packed k-mer representation */
    WordType seq_;
};


template <typename G, int L>
template <typename V>
KMer<G, L>::KMer(const V &arr, size_t k) : seq_(0) {
    assert(k * kBitsPerChar <= sizeof(WordType) * 8 && k >= 1);

    for (int i = k - 1; i > 0; --i) {
        assert(static_cast<CharType>(arr[i]) <= kFirstCharMask
                 && "Too small Digit size for representing the character");

        seq_ |= arr[i];
        seq_ <<= kBitsPerChar;
    }

    assert(static_cast<CharType>(arr[0]) <= kFirstCharMask
            && "Too small Digit size for representing the character");

    seq_ |= arr[0];
}

template <typename G, int L>
void KMer<G, L>::to_next(size_t k, WordType edge_label) {
    assert(edge_label <= static_cast<WordType>(kFirstCharMask));
    assert(k * kBitsPerChar <= sizeof(WordType) * 8);
    seq_ >>= kBitsPerChar;
    seq_ |= edge_label << static_cast<int>(kBitsPerChar * (k - 1));
}

template <typename G, int L>
void KMer<G, L>::to_prev(size_t k, CharType first_char) {
    assert(k * kBitsPerChar <= sizeof(WordType) * 8);
    seq_ <<= kBitsPerChar;
    seq_ &= kAllSetMask >> (sizeof(WordType) * 8 - kBitsPerChar * k);
    seq_ |= first_char;
}

template <typename G, int L>
typename KMer<G, L>::CharType KMer<G, L>::operator[](size_t i) const {
    static_assert(kBitsPerChar <= 64, "Too large digit!");
    assert(kBitsPerChar * (i + 1) <= sizeof(WordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}


template <typename G, int L>
std::ostream& operator<<(std::ostream &os, const KMer<G, L> &kmer) {
    kmer.print_hex(os);
    return os;
}

} // namespace kmer
} // namespace mtg

#endif // __KMER_HPP__
