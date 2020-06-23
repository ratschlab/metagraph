#include "kmer_boss.hpp"

#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace kmer {

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
void KMerBOSS<G, L>::print_hex(std::ostream &os) const {
    os << sdsl::uint256_t(seq_);
}

template class KMerBOSS<uint64_t, 3>;
template class KMerBOSS<sdsl::uint128_t, 3>;
template class KMerBOSS<sdsl::uint256_t, 3>;

template class KMerBOSS<uint64_t, 4>;
template class KMerBOSS<sdsl::uint128_t, 4>;
template class KMerBOSS<sdsl::uint256_t, 4>;

template class KMerBOSS<uint64_t, 5>;
template class KMerBOSS<sdsl::uint128_t, 5>;
template class KMerBOSS<sdsl::uint256_t, 5>;

} // namespace kmer
} // namespace mtg
