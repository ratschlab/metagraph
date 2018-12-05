#include "kmer_packed.hpp"

#include <sdsl/uint256_t.hpp>
#include <sdsl/bits.hpp>


template <typename G, int L>
const typename KMerPacked<G, L>::KMerCharType
KMerPacked<G, L>::kFirstCharMask = (1llu << kBitsPerChar) - 1;

template <typename G, int L>
std::string
KMerPacked<G, L>::to_string(size_t k, const std::string &alphabet) const {
    std::string seq(k, '\0');
    for (size_t i = 0; i < seq.length(); ++i) {
        seq[i] = alphabet.at(get_digit(i));
    }

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
