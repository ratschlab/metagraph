#include "kmer.hpp"

#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace kmer {

template <typename G, int L>
std::string KMer<G, L>::to_string(size_t k, const std::string &alphabet) const {
    std::string seq(k, '\0');

    for (size_t i = 0; i < k; ++i) {
        assert(operator[](i) < alphabet.size());
        seq[i] = alphabet[operator[](i)];
    }

    return seq;
}

template <typename G, int L>
void KMer<G, L>::print_hex(std::ostream &os) const {
    os << sdsl::uint256_t(seq_);
}

template class KMer<uint64_t, 2>;
template class KMer<sdsl::uint128_t, 2>;
template class KMer<sdsl::uint256_t, 2>;

template class KMer<uint64_t, 3>;
template class KMer<sdsl::uint128_t, 3>;
template class KMer<sdsl::uint256_t, 3>;

template class KMer<uint64_t, 4>;
template class KMer<sdsl::uint128_t, 4>;
template class KMer<sdsl::uint256_t, 4>;

template class KMer<uint64_t, 5>;
template class KMer<sdsl::uint128_t, 5>;
template class KMer<sdsl::uint256_t, 5>;

} // namespace kmer
} // namespace mtg
