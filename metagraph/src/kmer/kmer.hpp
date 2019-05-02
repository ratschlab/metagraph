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
    typedef G WordType;
    typedef uint64_t CharType;
    static constexpr int kBitsPerChar = L;

    KMer() {}
    template <typename V>
    KMer(const V &arr, size_t k);
    template <typename T>
    KMer(const std::vector<T> &arr) : KMer(arr, arr.size()) {}

    explicit KMer(WordType &&seq) : seq_(seq) {}
    explicit KMer(const WordType &seq) : seq_(seq) {}

    // corresponds to the co-lex order of k-mers
    bool operator<(const KMer &other) const { return seq_ < other.seq_; }
    bool operator<=(const KMer &other) const { return seq_ <= other.seq_; }
    bool operator>(const KMer &other) const { return seq_ > other.seq_; }
    bool operator>=(const KMer &other) const { return seq_ >= other.seq_; }
    bool operator==(const KMer &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMer &other) const { return seq_ != other.seq_; }

    inline CharType operator[](size_t i) const;

    std::string to_string(size_t k, const std::string &alphabet) const;

    /**
     * Construct the next k-mer for s[7]s[6]s[5]s[4]s[3]s[2]s[1].
     * next = s[8]s[7]s[6]s[5]s[4]s[3]s[2]
     *      = ( s[8] << k ) + ( kmer >> 1 ).
     */
    inline void to_next(size_t k, CharType edge_label);
    inline void to_prev(size_t k, CharType first_char);

    inline const WordType& data() const { return seq_; }

    template <typename T>
    inline static bool match_suffix(const T *kmer, size_t k, const std::vector<T> &suffix) {
        assert(k > 0);
        assert(k >= suffix.size());
        return suffix.empty()
                || std::equal(suffix.begin(), suffix.end(), kmer + k - suffix.size());
    }

  private:
    static const CharType kFirstCharMask;
    WordType seq_; // kmer sequence
};


template <typename G, int L>
template <typename V>
KMer<G, L>::KMer(const V &arr, size_t k) : seq_(0) {
    if (k * kBitsPerChar > sizeof(WordType) * 8 || k < 1) {
        std::cerr << "ERROR: Invalid k-mer size "
                  << k << ": must be between 1 and "
                  << sizeof(WordType) * 8 / kBitsPerChar << std::endl;
        exit(1);
    }

    for (int i = k - 1; i >= 0; --i) {
        assert(static_cast<uint64_t>(arr[i]) <= kFirstCharMask
                 && "Too small Digit size for representing the character");

        seq_ = seq_ << kBitsPerChar;
        seq_ |= arr[i];
    }
}

template <typename G, int L>
void KMer<G, L>::to_next(size_t k, CharType edge_label) {
    seq_ = seq_ >> kBitsPerChar;
    seq_ |= WordType(edge_label) << static_cast<int>(kBitsPerChar * (k - 1));
}

template <typename G, int L>
void KMer<G, L>::to_prev(size_t k, CharType first_char) {
    seq_ = seq_ << kBitsPerChar;
    seq_ = seq_ & ((WordType(1llu) << static_cast<int>(kBitsPerChar * k)) - WordType(1));
    seq_ |= first_char;
}

template <typename G, int L>
typename KMer<G, L>::CharType KMer<G, L>::operator[](size_t i) const {
    static_assert(kBitsPerChar <= 64, "Too large digit!");
    assert(kBitsPerChar * (i + 1) <= sizeof(WordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}

#endif // __KMER_HPP__
