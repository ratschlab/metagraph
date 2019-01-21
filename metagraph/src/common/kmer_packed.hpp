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
    static constexpr int kBitsPerChar = L;

    KMerPacked() {}
    template <typename V>
    KMerPacked(const V &arr, size_t k);
    template <typename T>
    KMerPacked(const std::vector<T> &arr) : KMerPacked(arr, arr.size()) {}

    explicit KMerPacked(KMerWordType &&seq) : seq_(seq) {}
    explicit KMerPacked(const KMerWordType &seq) : seq_(seq) {}

    // corresponds to the co-lex order of k-mers
    bool operator<(const KMerPacked &other) const { return seq_ < other.seq_; }
    bool operator<=(const KMerPacked &other) const { return seq_ <= other.seq_; }
    bool operator>(const KMerPacked &other) const { return seq_ > other.seq_; }
    bool operator>=(const KMerPacked &other) const { return seq_ >= other.seq_; }
    bool operator==(const KMerPacked &other) const { return seq_ == other.seq_; }
    bool operator!=(const KMerPacked &other) const { return seq_ != other.seq_; }

    inline KMerCharType operator[](size_t i) const;

    std::string to_string(size_t k, const std::string &alphabet) const;

    /**
     * Construct the next k-mer for s[7]s[6]s[5]s[4]s[3]s[2]s[1].
     * next = s[8]s[7]s[6]s[5]s[4]s[3]s[2]
     *      = ( s[8] << k ) + ( kmer >> 1 ).
     */
    inline void to_next(size_t k, KMerCharType edge_label);
    inline void to_prev(size_t k, KMerCharType first_char);

    inline const KMerWordType& data() const { return seq_; }

  private:
    static const KMerCharType kFirstCharMask;
    KMerWordType seq_; // kmer sequence
};


template <typename G, int L>
template <typename V>
KMerPacked<G, L>::KMerPacked(const V &arr, size_t k) : seq_(0) {
    if (k * kBitsPerChar > sizeof(KMerWordType) * 8 || k < 1) {
        std::cerr << "ERROR: Invalid k-mer size "
                  << k << ": must be between 1 and "
                  << sizeof(KMerWordType) * 8 / kBitsPerChar << std::endl;
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
void KMerPacked<G, L>::to_next(size_t k, KMerCharType edge_label) {
    seq_ = seq_ >> kBitsPerChar;
    seq_ |= KMerWordType(edge_label) << static_cast<int>(kBitsPerChar * (k - 1));
}

template <typename G, int L>
void KMerPacked<G, L>::to_prev(size_t k, KMerCharType first_char) {
    seq_ = seq_ << kBitsPerChar;
    seq_ = seq_ & ((KMerWordType(1llu) << static_cast<int>(kBitsPerChar * k)) - KMerWordType(1));
    seq_ |= first_char;
}

template <typename G, int L>
typename KMerPacked<G, L>::KMerCharType KMerPacked<G, L>::operator[](size_t i) const {
    static_assert(kBitsPerChar <= 64, "Too large digit!");
    assert(kBitsPerChar * (i + 1) <= sizeof(KMerWordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}

#endif // __KMER_PACKED_HPP__
