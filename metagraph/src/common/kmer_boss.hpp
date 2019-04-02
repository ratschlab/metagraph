#ifndef __KMER_BOSS_HPP__
#define __KMER_BOSS_HPP__

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>


template <typename G, int L>
class KMerBOSS {
    template <typename U, int P>
    friend std::ostream& operator<<(std::ostream &os, const KMerBOSS<U, P> &kmer);

  public:
    typedef G WordType;
    typedef uint64_t CharType;
    static constexpr int kBitsPerChar = L;

    KMerBOSS() {}
    template <typename V>
    KMerBOSS(const V &arr, size_t k);
    template <typename T>
    KMerBOSS(const std::vector<T> &arr) : KMerBOSS(arr, arr.size()) {}

    explicit KMerBOSS(WordType &&seq) : seq_(seq) {}
    explicit KMerBOSS(const WordType &seq) : seq_(seq) {}

    // corresponds to the BOSS (co-lex, one-swapped) order of k-mers
    bool operator<(const KMerBOSS &other) const { return seq_ < other.seq_; }
    bool operator<=(const KMerBOSS &other) const { return seq_ <= other.seq_; }
    bool operator>(const KMerBOSS &other) const { return seq_ > other.seq_; }
    bool operator>=(const KMerBOSS &other) const { return seq_ >= other.seq_; }
    bool operator==(const KMerBOSS &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMerBOSS &other) const { return seq_ != other.seq_; }

    inline CharType operator[](size_t i) const;

    /**
     * Compares k-mers without one last and |minus| first characters.
     * Examples: For s[6]s[5]s[4]s[3]s[2]s[1]s[7],
     *                  compares s[6]s[5]s[4]s[3]s[2]s[1] if minus = 0.
     * In general, checks if s[minus+1]...s[k-1] are the same for both kmers.
     */
    static inline bool compare_suffix(const KMerBOSS &k1,
                                      const KMerBOSS &k2, size_t minus = 0);

    std::string to_string(size_t k, const std::string &alphabet) const;

    /**
     * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
     * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
     *      = s[7] << k + (kmer & mask) >> 1 + s[8].
     */
    inline void to_next(size_t k, CharType new_last, CharType old_last);
    inline void to_next(size_t k, CharType new_last);
    inline void to_prev(size_t k, CharType new_first);

    inline const WordType& data() const { return seq_; }

  private:
    static const CharType kFirstCharMask;
    WordType seq_; // kmer sequence
};


template <typename G, int L>
template <typename V>
KMerBOSS<G, L>::KMerBOSS(const V &arr, size_t k) : seq_(0) {
    if (k * kBitsPerChar > sizeof(WordType) * 8 || k < 2) {
        std::cerr << "ERROR: Invalid k-mer size: passed "
                  << k << " but must be between 2 and "
                  << sizeof(WordType) * 8 / kBitsPerChar << std::endl;
        exit(1);
    }

    for (int i = k - 2; i >= 0; --i) {
        assert(static_cast<uint64_t>(arr[i]) <= kFirstCharMask
                 && "Too small Digit size for representing the character");

        seq_ |= arr[i];
        seq_ = seq_ << kBitsPerChar;
    }
    seq_ |= arr[k - 1];
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMerBOSS<G, L>::to_next(size_t k, CharType new_last, CharType old_last) {
    assert(old_last == (seq_ & WordType(kFirstCharMask)));
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    seq_ = seq_ >> kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ += WordType(old_last) << static_cast<int>(kBitsPerChar * (k - 1));
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ |= kFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]1111
    seq_ -= kFirstCharMask - new_last;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
}

/**
 * Construct the next k-mer for s[6]s[5]s[4]s[3]s[2]s[1]s[7].
 * next = s[7]s[6]s[5]s[4]s[3]s[2]s[8]
 *      = s[7] << k + (kmer & mask) >> 1 + s[8].
 */
template <typename G, int L>
void KMerBOSS<G, L>::to_next(size_t k, CharType new_last) {
    WordType old_last = seq_ & WordType(kFirstCharMask);
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
    seq_ = seq_ >> kBitsPerChar;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ += old_last << static_cast<int>(kBitsPerChar * (k - 1));
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ |= kFirstCharMask;
    // s[7]s[6]s[5]s[4]s[3]s[2]1111
    seq_ -= kFirstCharMask - new_last;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
}

template <typename G, int L>
void KMerBOSS<G, L>::to_prev(size_t k, CharType new_first) {
    const int shift = kBitsPerChar * (k - 1);
    WordType last_char = seq_ >> shift;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[8]
    seq_ |= kFirstCharMask;
    seq_ -= kFirstCharMask - new_first;
    // s[7]s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ -= last_char << shift;
    // 0000s[6]s[5]s[4]s[3]s[2]s[1]
    seq_ = seq_ << kBitsPerChar;
    // s[6]s[5]s[4]s[3]s[2]s[1]0000
    seq_ += last_char;
    // s[6]s[5]s[4]s[3]s[2]s[1]s[7]
}

template <typename G, int L>
typename KMerBOSS<G, L>::CharType KMerBOSS<G, L>::operator[](size_t i) const {
    static_assert(kBitsPerChar <= 64, "Too large digit!");
    assert(kBitsPerChar * (i + 1) <= sizeof(WordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
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

#endif // __KMER_BOSS_HPP__
