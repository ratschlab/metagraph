#include "kmer_boss.hpp"

#include <sdsl/uint256_t.hpp>


template <typename G, int L>
const typename KMerBOSS<G, L>::CharType
KMerBOSS<G, L>::kFirstCharMask = (1llu << kBitsPerChar) - 1;

template <typename G, int L>
std::string KMerBOSS<G, L>::to_string(size_t k, const std::string &alphabet) const {
    std::string seq(k, '\0');

    for (size_t i = 0; i + 1 < k; ++i) {
        assert(operator[](i + 1) < alphabet.size());
        seq[i] = alphabet[operator[](i + 1)];
    }
    assert(operator[](0) < alphabet.size());
    seq[k - 1] = alphabet[operator[](0)];

    return seq;
}

template <typename G, int L>
std::ostream& operator<<(std::ostream &os, const KMerBOSS<G, L> &kmer) {
    return os << sdsl::uint256_t(kmer.seq_);
}

template std::ostream&
operator<<<uint64_t, 3>(std::ostream &os, const KMerBOSS<uint64_t, 3> &kmer);
template std::ostream&
operator<<<sdsl::uint128_t, 3>(std::ostream &os, const KMerBOSS<sdsl::uint128_t, 3> &kmer);
template std::ostream&
operator<<<sdsl::uint256_t, 3>(std::ostream &os, const KMerBOSS<sdsl::uint256_t, 3> &kmer);

template class KMerBOSS<uint64_t, 3>;
template class KMerBOSS<sdsl::uint128_t, 3>;
template class KMerBOSS<sdsl::uint256_t, 3>;

template std::ostream&
operator<<<uint64_t, 4>(std::ostream &os, const KMerBOSS<uint64_t, 4> &kmer);
template std::ostream&
operator<<<sdsl::uint128_t, 4>(std::ostream &os, const KMerBOSS<sdsl::uint128_t, 4> &kmer);
template std::ostream&
operator<<<sdsl::uint256_t, 4>(std::ostream &os, const KMerBOSS<sdsl::uint256_t, 4> &kmer);

template class KMerBOSS<uint64_t, 4>;
template class KMerBOSS<sdsl::uint128_t, 4>;
template class KMerBOSS<sdsl::uint256_t, 4>;

template std::ostream&
operator<<<uint64_t, 5>(std::ostream &os, const KMerBOSS<uint64_t, 5> &kmer);
template std::ostream&
operator<<<sdsl::uint128_t, 5>(std::ostream &os, const KMerBOSS<sdsl::uint128_t, 5> &kmer);
template std::ostream&
operator<<<sdsl::uint256_t, 5>(std::ostream &os, const KMerBOSS<sdsl::uint256_t, 5> &kmer);

template class KMerBOSS<uint64_t, 5>;
template class KMerBOSS<sdsl::uint128_t, 5>;
template class KMerBOSS<sdsl::uint256_t, 5>;
