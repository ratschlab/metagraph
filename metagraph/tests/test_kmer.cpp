#include "gtest/gtest.h"

#include <vector>
#include <functional>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "kmer.hpp"


template <typename T>
T encode(char c, const T *char_map) {
    assert(static_cast<size_t>(c) < 128);
    return char_map[static_cast<size_t>(c)];
}

template <typename T>
char decode(T a, const std::string &alphabet) {
    assert(a < alphabet.size());
    return alphabet[a];
}

template <typename T>
std::vector<T> encode(const std::string &sequence, const T *char_map) {
    std::vector<T> encoded;
    std::transform(sequence.begin(), sequence.end(),
                   std::back_inserter(encoded),
                   [&](char c) { return encode(c, char_map); });
    assert(encoded.size() == sequence.size());

    return encoded;
}

template <typename T>
std::string decode(const std::vector<T> &encoded, const std::string &alphabet) {
    std::string decoded;
    std::transform(encoded.begin(), encoded.end(),
                   std::back_inserter(decoded),
                   [&](T a) { return decode(a, alphabet); });
    assert(decoded.size() == encoded.size());

    return decoded;
}

typedef uint8_t TAlphabet;


 // Nucleotide
std::vector<TAlphabet> encode_nucleotide(const std::string &sequence) {
    const TAlphabet kCharToNucleotide[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    return encode(sequence, kCharToNucleotide);
}

template <typename G, int L>
std::string decode_nucleotide(const KMer<G, L> &kmer) {
    return kmer.to_string("ACGTN");
}



template <typename G, int L>
void test_kmer_codec(const std::string &sequence,
                     const std::function<std::vector<TAlphabet>(const std::string&)> &encode,
                     const std::function<std::string(const KMer<G, L>&)> &decode) {
    const auto encoded = encode(sequence);

    for (uint64_t k = 2; k < sizeof(G) * 8 / L; ++k) {
        if (k > encoded.size())
            break;

        ASSERT_LE(k, encoded.size());
        std::vector<KMer<G, L>> kmers;
        auto packed = KMer<G, L>::pack_kmer(encoded.data(), k);
        for (uint64_t i = 0; i + k <= encoded.size(); ++i) {
            kmers.emplace_back(std::vector<TAlphabet>(encoded.begin() + i,
                                                      encoded.begin() + i + k));
            ASSERT_EQ(sequence.substr(i, k), decode(kmers.back()))
                 << k << " " << i;

            KMer<G, L> kmer_alt(encoded.data() + i, k);
            ASSERT_EQ(kmers.back(), kmer_alt) << k << " " << i;

            KMer<G, L> kmer_packed(packed);
            ASSERT_EQ(kmers.back(), kmer_packed) << k << " " << i;
            if (i + k < encoded.size())
                KMer<G, L>::update_kmer(
                    k,
                    encoded[i + k],
                    encoded[i + k - 1],
                    &packed);
        }
    }
}

TEST(Kmer, nucleotide_alphabet_pack_64) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    test_kmer_codec<uint64_t, 3>(sequence, encode_nucleotide, decode_nucleotide<uint64_t, 3>);
}

TEST(Kmer, nucleotide_alphabet_pack_128) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    test_kmer_codec<sdsl::uint128_t, 3>(sequence, encode_nucleotide, decode_nucleotide<sdsl::uint128_t, 3>);
}

TEST(Kmer, nucleotide_alphabet_pack_256) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    test_kmer_codec<sdsl::uint256_t, 3>(sequence, encode_nucleotide, decode_nucleotide<sdsl::uint256_t, 3>);
}
