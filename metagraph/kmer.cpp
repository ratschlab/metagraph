#include "kmer.hpp"

#include <sdsl/uint256_t.hpp>

template <typename G>
const KMerCharType KMer<G>::kFirstCharMask = (1llu << kBitsPerChar) - 1;


template <typename G>
KMerCharType KMer<G>::get_digit(size_t i) const {
    static_assert(kBitsPerChar <= 64, "too large digit");
    assert(kBitsPerChar * (i + 1) <= sizeof(KMerWordType) * 8);
    return static_cast<uint64_t>(seq_ >> static_cast<int>(kBitsPerChar * i))
             & kFirstCharMask;
}

template <typename G>
bool KMer<G>::compare_suffix(const KMer &k1, const KMer &k2, size_t minus) {
    return k1.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar)
             == k2.seq_ >> static_cast<int>((minus + 1) * kBitsPerChar);
}

template <typename G>
KMerCharType KMer<G>::operator[](size_t i) const {
    assert(get_digit(i) > 0);
    return get_digit(i) - 1;
}

template <typename G>
std::string KMer<G>::to_string(const std::string &alphabet) const {
    std::string seq;

    const size_t max_len = sizeof(KMerWordType) * 8 / kBitsPerChar;
    seq.reserve(max_len);

    KMerCharType cur;
    for (size_t i = 0; i < max_len && (cur = get_digit(i)); ++i) {
        seq.push_back(alphabet.at(cur - 1));
    }
    return seq;
}

template <typename G>
std::ostream& operator<<(std::ostream &os, const KMer<G> &kmer) {
    return os << sdsl::uint256_t(kmer.seq_);
}

template std::ostream&
operator<<<uint64_t>(std::ostream &os, const KMer<uint64_t> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t>(std::ostream &os, const KMer<sdsl::uint128_t> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t>(std::ostream &os, const KMer<sdsl::uint256_t> &kmer);


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

template class KMer<uint64_t>;
template class KMer<sdsl::uint128_t>;
template class KMer<sdsl::uint256_t>;
