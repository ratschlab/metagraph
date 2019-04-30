#include "kmer.hpp"

#include <sdsl/uint256_t.hpp>
#include <sdsl/bits.hpp>


template <typename G, int L>
const typename KMer<G, L>::CharType
KMer<G, L>::kFirstCharMask = (1llu << kBitsPerChar) - 1;

template <typename G, int L>
std::string
KMer<G, L>::to_string(size_t k, const std::string &alphabet) const {
    std::string seq(k, '\0');

    for (size_t i = 0; i < k; ++i) {
        assert(operator[](i) < alphabet.size());
        seq[i] = alphabet[operator[](i)];
    }

    return seq;
}

template <typename G, int L>
std::ostream& operator<<(std::ostream &os, const KMer<G, L> &kmer) {
    return os << sdsl::uint256_t(kmer.seq_);
}

template std::ostream&
operator<<<uint64_t, 2>(std::ostream &os, const KMer<uint64_t, 2> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t, 2>(std::ostream &os, const KMer<sdsl::uint128_t, 2> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t, 2>(std::ostream &os, const KMer<sdsl::uint256_t, 2> &kmer);


template class KMer<uint64_t, 2>;
template class KMer<sdsl::uint128_t, 2>;
template class KMer<sdsl::uint256_t, 2>;


template std::ostream&
operator<<<uint64_t, 3>(std::ostream &os, const KMer<uint64_t, 3> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t, 3>(std::ostream &os, const KMer<sdsl::uint128_t, 3> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t, 3>(std::ostream &os, const KMer<sdsl::uint256_t, 3> &kmer);

template class KMer<uint64_t, 3>;
template class KMer<sdsl::uint128_t, 3>;
template class KMer<sdsl::uint256_t, 3>;


template std::ostream&
operator<<<uint64_t, 5>(std::ostream &os, const KMer<uint64_t, 5> &kmer);

template std::ostream&
operator<<<sdsl::uint128_t, 5>(std::ostream &os, const KMer<sdsl::uint128_t, 5> &kmer);

template std::ostream&
operator<<<sdsl::uint256_t, 5>(std::ostream &os, const KMer<sdsl::uint256_t, 5> &kmer);

template class KMer<uint64_t, 5>;
template class KMer<sdsl::uint128_t, 5>;
template class KMer<sdsl::uint256_t, 5>;
