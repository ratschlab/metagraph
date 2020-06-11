#ifndef __KMER_BOSS_HPP__
#define __KMER_BOSS_HPP__

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
 * Models a kmer (https://en.wikipedia.org/wiki/K-mer) that is stored in a BOSS table.
 * Just like #KMer, each character in the k-mer uses L bits of the internal representation
 * type G (typically a 64, 128 or 256 bit integer) and characters are stored in reverse
 * order.
 * The  only  difference from #KMer is that #KmerBoss leaves the last character in the
 * least significant bits of G. For example, the k-mer "ACGT" is stored in memory as
 * "GCAT", while #KMer stores it as "TGCA".
 *
 * @tparam G the type storing a kmer, typically a 64/128/256 bit integer
 * @tparam L the number of bits storing a k-mer character (2 bits for DNA)
 */
template <typename G, int L>
class KMerBOSS {
  public:
    typedef G WordType;
    typedef uint64_t CharType;
    static constexpr int kBitsPerChar = L;
    static constexpr CharType kFirstCharMask = (1ull << kBitsPerChar) - 1;

    /** Construct a default, uninitialized, BOSS k-mer. */
    KMerBOSS() {}
    /**
     * Construct a BOSS k-mer with the given size from the given array
     * @tparam V an indexed data structure (e.g. std::vector or L[])
     * @param arr the k-mer characters
     * @param k k-mer length
     */
    template <typename V>
    KMerBOSS(const V &arr, size_t k);
    /**
     * Construct a BOSS k-mer from the given vector
     * @tparam T k-mer character type
     * @param arr unpacked k-mer (e.g. 'ACATGGAA')
     */
    template <typename T>
    KMerBOSS(const std::vector<T> &arr) : KMerBOSS(arr, arr.size()) {}

    KMerBOSS(WordType&& seq) noexcept : seq_(seq) {}
    explicit KMerBOSS(const WordType &seq) noexcept : seq_(seq) {}

    // corresponds to the BOSS (co-lex, one-swapped) order of k-mers
    bool operator<(const KMerBOSS &other) const { return seq_ < other.seq_; }
    bool operator<=(const KMerBOSS &other) const { return seq_ <= other.seq_; }
    bool operator>(const KMerBOSS &other) const { return seq_ > other.seq_; }
    bool operator>=(const KMerBOSS &other) const { return seq_ >= other.seq_; }
    bool operator==(const KMerBOSS &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMerBOSS &other) const { return seq_ != other.seq_; }

    /**
     * Returns the character at position i. For the k-mer a_1|a_2|...|a_(k-1)|a_k, which is laid out
     * in memory as a_(k-1)|a_(k-2)|...|a_1|a_k, kmer[0]=a_k, kmer[1]=a_1, ..., kmer[k-1]=a_(k-1)
     */
    inline CharType operator[](size_t i) const;

    /**
     * Compares k-mers without one last and |minus| first characters.
     * Examples: For s[6]s[5]s[4]s[3]s[2]s[1]s[7],
     *                  compares s[6]s[5]s[4]s[3]s[2]s[1] if minus = 0.
     * In general, checks if s[minus+1]...s[k-1] are the same for both kmers.
     */
    static inline bool compare_suffix(const KMerBOSS &k1,
                                      const KMerBOSS &k2, size_t minus = 0);

    static inline bool less(const KMerBOSS &k1, const KMerBOSS &k2, size_t minus = 0);

    std::string to_string(size_t k, const std::string &alphabet) const;

    /**
     * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
     * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
     *      = s[7] << k + (kmer & mask) >> 1 + s[8].
     * @param k k-mer size
     * @param new_last the character to append to the current k-mer
     * @param old_last the last character in the current k-mer (provided, if available,
     * for efficiency reasons)
     */
    inline void to_next(size_t k, CharType new_last, WordType old_last);
    inline void to_next(size_t k, CharType new_last);
    /**
     * Replaces the current k-mer with its predecessor. The last k-1 characters of the
     * predecessors are identical to the first k-1 characters of the successors.
     * @param k the k-mer size
     * @param first_char the first character in the predecessor
     */
    inline void to_prev(size_t k, CharType new_first);

    inline const WordType& data() const { return seq_; }

    template <typename T>
    inline static bool match_suffix(const T *kmer, size_t k, const std::vector<T> &suffix) {
        assert(k > 1);
        assert(k > suffix.size());
        return suffix.empty()
                || std::equal(suffix.begin(), suffix.end(), kmer + k - suffix.size() - 1);
    }

    void print_hex(std::ostream &os) const;

  private:
    static inline const WordType kAllButFirstCharMask = ~(WordType(kFirstCharMask));
    static inline const WordType kAllSetMask = ~(WordType(0ull));
    WordType seq_; // kmer sequence
};


template <typename G, int L>
template <typename V>
KMerBOSS<G, L>::KMerBOSS(const V &arr, size_t k) : seq_(0) {
    assert(k * kBitsPerChar <= sizeof(WordType) * 8 && k >= 2);

    for (int i = k - 2; i >= 0; --i) {
        assert(static_cast<uint64_t>(arr[i]) <= kFirstCharMask
                 && "Too small Digit size for representing the character");

        seq_ |= arr[i];
        seq_ <<= kBitsPerChar;
    }

    assert(static_cast<uint64_t>(arr[k - 1]) <= kFirstCharMask
            && "Too small Digit size for representing the character");

    seq_ |= arr[k - 1];
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMerBOSS<G, L>::to_next(size_t k, CharType new_last, WordType old_last) {
    assert(old_last == (seq_ & kFirstCharMask));
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    seq_ >>= kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ |= old_last << static_cast<int>(kBitsPerChar * (k - 1));
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ &= kAllButFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]0000
    seq_ |= new_last;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMerBOSS<G, L>::to_next(size_t k, CharType new_last) {
    to_next(k, new_last, seq_ & static_cast<WordType>(kFirstCharMask));
}

template <typename G, int L>
void KMerBOSS<G, L>::to_prev(size_t k, CharType new_first) {
    const int shift = kBitsPerChar * (k - 1);
    WordType last_char = seq_ >> shift;
    //     s[7]s[6]s[5]s[4]s[3]s[2]s[8]
    seq_ &= kAllButFirstCharMask;
    //     s[7]s[6]s[5]s[4]s[3]s[2]0000
    seq_ |= new_first;
    //     s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ <<= kBitsPerChar;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]0000
    seq_ &= kAllSetMask >> (sizeof(WordType) * 8 - kBitsPerChar * k);
    //     s[6]s[5]s[4]s[3]s[2]s[1]0000
    seq_ |= last_char;
    //     s[6]s[5]s[4]s[3]s[2]s[1]s[7]
}

template <typename G, int L>
typename KMerBOSS<G, L>::CharType KMerBOSS<G, L>::operator[](size_t i) const {
    static_assert(kBitsPerChar <= 64, "Too large digit!");
    assert(kBitsPerChar * (i + 1) <= sizeof(WordType) * 8);
    return static_cast<CharType>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}

/**
 * Compares k-mers without one last and |minus| first characters.
 * Examples: For s[6]s[5]s[4]s[3]s[2]s[1]s[7],
 *                  compares s[6]s[5]s[4]s[3]s[2]s[1] if minus = 0.
 * In general, checks if s[minus+1]...s[k-1] are the same for both kmers.
 */
template <typename G, int L>
bool KMerBOSS<G, L>::compare_suffix(const KMerBOSS &k1, const KMerBOSS &k2, size_t minus) {
    return k1.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar)
             == k2.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar);
}

template <typename G, int L>
bool KMerBOSS<G, L>::less(const KMerBOSS &k1, const KMerBOSS &k2, size_t minus) {
    return k1.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar)
            < k2.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar);
}


template <typename G, int L>
std::ostream& operator<<(std::ostream &os, const KMerBOSS<G, L> &kmer) {
    kmer.print_hex(os);
    return os;
}

} // namespace kmer
} // namespace mtg

#endif // __KMER_BOSS_HPP__
