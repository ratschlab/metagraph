#include "kmer_packed.hpp"

#include <sdsl/uint256_t.hpp>
#include <sdsl/bits.hpp>


template <typename G, int L>
uint32_t
KMerPacked<G, L>::get_k(const KMerWordType &kmer) {
    auto limb = reinterpret_cast<uint64_t const*>(&kmer);
    if (sizeof(kmer) == 32) {
        if (limb[3]) {
            return (sdsl::bits::hi(limb[3]) + 192) / kBitsPerChar;;
        } else if (limb[2]) {
            return (sdsl::bits::hi(limb[2]) + 128) / kBitsPerChar;
        } else if (limb[1]) {
            return (sdsl::bits::hi(limb[1]) + 64) / kBitsPerChar;
        } else {
            return sdsl::bits::hi(limb[0]) / kBitsPerChar;
        }
    } else if (sizeof(kmer) == 16) {
        return limb[1]
            ? (sdsl::bits::hi(limb[1]) + 64) / kBitsPerChar
            : sdsl::bits::hi(limb[0]) / kBitsPerChar;
    } else {
        return sdsl::bits::hi(limb[0]) / kBitsPerChar;
    }
}

template <typename G, int L>
void
KMerPacked<G, L>::set_length_bit(KMerWordType *kmer, uint64_t k) {
    auto limb = reinterpret_cast<uint64_t*>(kmer);
    size_t length = k * kBitsPerChar;
    limb[length >> 6] |= 1llu << (length & 63);
}

template <typename G, int L>
void
KMerPacked<G, L>::unset_length_bit(KMerWordType *kmer, uint64_t k) {
    auto limb = reinterpret_cast<uint64_t*>(kmer);
    size_t length = k * kBitsPerChar;
    limb[length >> 6] &= (1llu << (length & 63)) - 1;
}

template <typename G, int L>
const typename KMerPacked<G, L>::KMerCharType
KMerPacked<G, L>::kFirstCharMask = (1llu << kBitsPerChar) - 1;

template <typename G, int L>
std::string
KMerPacked<G, L>::to_string(const std::string &alphabet) const {
    std::string seq(get_k(), '\0');
    for (size_t i = 1; i < seq.length(); ++i) {
        seq[i - 1] = alphabet.at(get_digit(i));
    }
    seq[seq.length() - 1] = alphabet.at(get_digit(0));

    return seq;
}

template <typename G, int L>
std::ostream& operator<<(std::ostream &os, const KMerPacked<G, L> &kmer) {
    return os << sdsl::uint256_t(kmer.seq_);
}

template std::ostream&
operator<<<uint64_t, 2>(std::ostream &os, const KMerPacked<uint64_t, 2> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t, 2>(std::ostream &os, const KMerPacked<sdsl::uint128_t, 2> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t, 2>(std::ostream &os, const KMerPacked<sdsl::uint256_t, 2> &kmer);


template class KMerPacked<uint64_t, 2>;
template class KMerPacked<sdsl::uint128_t, 2>;
template class KMerPacked<sdsl::uint256_t, 2>;


template std::ostream&
operator<<<uint64_t, 3>(std::ostream &os, const KMerPacked<uint64_t, 3> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t, 3>(std::ostream &os, const KMerPacked<sdsl::uint128_t, 3> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t, 3>(std::ostream &os, const KMerPacked<sdsl::uint256_t, 3> &kmer);

template class KMerPacked<uint64_t, 3>;
template class KMerPacked<sdsl::uint128_t, 3>;
template class KMerPacked<sdsl::uint256_t, 3>;


template std::ostream&
operator<<<uint64_t, 5>(std::ostream &os, const KMerPacked<uint64_t, 5> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t, 5>(std::ostream &os, const KMerPacked<sdsl::uint128_t, 5> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t, 5>(std::ostream &os, const KMerPacked<sdsl::uint256_t, 5> &kmer);

template class KMerPacked<uint64_t, 5>;
template class KMerPacked<sdsl::uint128_t, 5>;
template class KMerPacked<sdsl::uint256_t, 5>;
