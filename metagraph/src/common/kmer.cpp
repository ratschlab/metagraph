#include "kmer.hpp"

#include <sdsl/uint256_t.hpp>

template <typename G>
const KMerCharType KMer<G>::kFirstCharMask = (1llu << kBitsPerChar) - 1;

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


template class KMer<uint64_t>;
template class KMer<sdsl::uint128_t>;
template class KMer<sdsl::uint256_t>;
