#include "gtest/gtest.h"

#include <vector>
#include <functional>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "kmer_packed.hpp"


template <typename T>
T encode_c(char c, const T *char_map);

template <typename T>
char decode_c(T a, const std::string &alphabet);

template <typename T>
std::vector<T> encode_c(const std::string &sequence, const T *char_map);

template <typename T>
std::string decode_c(const std::vector<T> &encoded, const std::string &alphabet);

typedef uint8_t TAlphabet;

// Nucleotide 2 bit
std::vector<TAlphabet> encode_nucleotide_2bit(const std::string &sequence) {
    const TAlphabet kCharToNucleotide[128] = {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
    };

    return encode_c(sequence, kCharToNucleotide);
}

template <typename G, int L>
std::string decode_nucleotide_2bit(const KMerPacked<G, L> &kmer, size_t k) {
    return kmer.to_string(k, "ACGT");
}


template <typename G, int L>
void test_kmer_codec(const std::string &sequence,
                     const std::function<std::vector<TAlphabet>(const std::string&)> &encoder,
                     const std::function<std::string(const KMerPacked<G, L>&, size_t)> &decoder) {
    const auto encoded = encoder(sequence);

    for (uint64_t k = 2; k < sizeof(G) * 8 / L; ++k) {
        if (k > encoded.size())
            break;

        ASSERT_LE(k, encoded.size());
        std::vector<KMerPacked<G, L>> kmers;
        auto packed = KMerPacked<G, L>::pack_kmer(encoded.data(), k);
        for (uint64_t i = 0; i + k <= encoded.size(); ++i) {
            kmers.emplace_back(std::vector<TAlphabet>(encoded.begin() + i,
                                                      encoded.begin() + i + k));
            assert(i + k <= sequence.size());
            //ASSERT_EQ(k, kmers.back().get_k()) << sequence.substr(i, k) << " " << k << " " << i;
            ASSERT_EQ(sequence.substr(i, k), decoder(kmers.back(), k))
                 << sequence.substr(i, k) << " " << k << " " << i;

            KMerPacked<G, L> kmer_alt(encoded.data() + i, k);
            //ASSERT_EQ(k, kmer_alt.get_k()) << sequence.substr(i, k) << " " << k << " " << i;
            ASSERT_EQ(kmers.back(), kmer_alt) << sequence.substr(i, k) << " " << k << " " << i;

            KMerPacked<G, L> kmer_packed(packed);
            //ASSERT_EQ(k, kmer_packed.get_k()) << sequence.substr(i, k) << " " << k << " " << i;
            ASSERT_EQ(kmers.back(), kmer_packed) << sequence.substr(i, k) << " " << k << " " << i;
            if (i + k < encoded.size())
                KMerPacked<G, L>::update_kmer(
                    k,
                    encoded[i + k],
                    encoded[i + k - 1],
                    &packed);
        }
    }
}

TEST(KmerPacked, nucleotide_alphabet_pack_6_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<uint64_t, 2>(sequence, encode_nucleotide_2bit, decode_nucleotide_2bit<uint64_t, 2>);
}

TEST(KmerPacked, nucleotide_alphabet_pack_128_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint128_t, 2>(sequence, encode_nucleotide_2bit, decode_nucleotide_2bit<sdsl::uint128_t, 2>);
}

TEST(KmerPacked, nucleotide_alphabet_pack_256_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint256_t, 2>(sequence, encode_nucleotide_2bit, decode_nucleotide_2bit<sdsl::uint256_t, 2>);
}

TEST(KmerPacked, nucleotide_alphabet_pack_6) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<uint64_t, 3>(sequence, encode_nucleotide_2bit, decode_nucleotide_2bit<uint64_t, 3>);
}

TEST(KmerPacked, nucleotide_alphabet_pack_128) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint128_t, 3>(sequence, encode_nucleotide_2bit, decode_nucleotide_2bit<sdsl::uint128_t, 3>);
}

TEST(KmerPacked, nucleotide_alphabet_pack_256) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint256_t, 3>(sequence, encode_nucleotide_2bit, decode_nucleotide_2bit<sdsl::uint256_t, 3>);
}
