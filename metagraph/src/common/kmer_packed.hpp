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

    std::string to_string(size_t k, const std::string &alphabet) const;

    KMerPacked<G, L>
    reverse_complement(size_t k, const std::vector<KMerCharType> &complement_map) const;

    template <typename T>
    static inline KMerWordType pack_kmer(const T &arr, size_t k);

    static inline void update_kmer(size_t k,
                                   KMerCharType edge_label,
                                   KMerCharType last,
                                   KMerWordType *kmer);

    inline KMerPacked<G, L> prev_kmer(size_t k, KMerCharType first_char) const;
    inline void next_kmer(size_t k, KMerCharType edge_label);

    inline KMerWordType& data() { return seq_; }
    inline KMerWordType data() const { return seq_; }

  private:
    inline KMerCharType get_digit(size_t i) const;

    static const KMerCharType kFirstCharMask;

    KMerWordType seq_; // kmer sequence
};


template <typename G, int L>
typename KMerPacked<G, L>::KMerCharType KMerPacked<G, L>::operator[](size_t i) const {
    return get_digit(i);
}

template <typename G, int L>
KMerPacked<G, L>
KMerPacked<G, L>
::reverse_complement(size_t k,
                     const std::vector<KMerCharType> &complement_map) const {
    KMerWordType pack(0);
    KMerWordType storage = seq_;
    KMerWordType mask = kFirstCharMask;
    for (size_t i = 0; i < k; ++i) {
        pack = pack << kBitsPerChar;
        pack |= complement_map.at(storage & mask);
        storage = storage >> kBitsPerChar;
    }
    assert(storage == static_cast<KMerWordType>(0));
    return KMerPacked<G, L>(pack);
}

template <typename G, int L>
template <typename T>
G KMerPacked<G, L>::pack_kmer(const T &arr, size_t k) {
    if (k > sizeof(KMerWordType) * 8 / kBitsPerChar || k < 2) {
        std::cerr << "Invalid k-mer size "
                  << k
                  << ": must be between 2 and "
                  << (sizeof(KMerWordType) * 8 / kBitsPerChar) << std::endl;
        exit(1);
    }

    KMerWordType result(0);

    for (int i = k - 1; i >= 0; --i) {

        if (static_cast<uint64_t>(arr[i]) >= (1llu << kBitsPerChar)) {
            std::cerr << "Alphabet size too big for the given number of bits" << std::endl;
            exit(1);
        }

        result = result << kBitsPerChar;
        result |= arr[i];
    }
    return result;
}

template <typename G, int L>
void KMerPacked<G, L>::update_kmer(size_t k,
                                   KMerCharType edge_label,
                                   KMerCharType,
                                   KMerWordType *kmer) {
    *kmer = *kmer >> kBitsPerChar;
    *kmer |= KMerWordType(edge_label) << static_cast<int>(kBitsPerChar * (k - 1));
}

template <typename G, int L>
void KMerPacked<G, L>::next_kmer(size_t k, KMerCharType edge_label) {
    update_kmer(k, edge_label, 0, &seq_);
}

template <typename G, int L>
KMerPacked<G, L> KMerPacked<G, L>::prev_kmer(size_t k, KMerCharType first_char) const {
    KMerWordType kmer = seq_;
    kmer = kmer << kBitsPerChar;
    kmer = kmer
        & ((KMerWordType(1llu) << static_cast<int>(kBitsPerChar * k)) - KMerWordType(1));
    kmer |= first_char;
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
    return k1.seq_ >> static_cast<int>(minus * kBitsPerChar)
             == k2.seq_ >> static_cast<int>(minus * kBitsPerChar);
}

#endif // __KMER_PACKED_HPP__
