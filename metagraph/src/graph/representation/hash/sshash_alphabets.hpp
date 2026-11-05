#pragma once
#include <array>
#include <bit>
#include <kmer.hpp>
#include "kmer/alphabets.hpp"

namespace mtg::graph::detail {

// Use MetaGraph alphabet tables with SSHash's packed storage.
template <typename Word, auto &Alphabet, auto &Encode, auto &Complement>
struct AlphabetKmer : sshash::alpha_kmer_t<Word, std::bit_width(sizeof(Alphabet) - 2u), Alphabet> {
    using base = sshash::alpha_kmer_t<Word, std::bit_width(sizeof(Alphabet) - 2u), Alphabet>;
    using base::base;
    static constexpr bool has_reverse_complement = true;

    static uint64_t char_to_uint(unsigned char c) {
        if (c >= sizeof(Encode)) return uint64_t(-1);
        auto code = Encode[c];
        if (code >= base::alphabet_size) return uint64_t(-1);
        char letter = Alphabet[code];
        // Preserve strict SSHash validation: MetaGraph maps other symbols to N/X.
        bool same_letter = c == letter || (c >= 'a' && c <= 'z' && c - 'a' + 'A' == letter);
        bool rna_alias = ((c == 'U' && letter == 'T') || (c == 'u' && (letter == 'T' || letter == 't')));
        return same_letter || rna_alias ? code : uint64_t(-1);
    }
    static bool is_valid(char c) { return char_to_uint(static_cast<unsigned char>(c)) != uint64_t(-1); }

    inline static const auto canonicalize_basepair_reverse_map = [] {
        std::array<char, 256> map{};
        for (unsigned c = 0; c < map.size(); ++c) {
            auto code = char_to_uint(c);
            if (code != uint64_t(-1)) map[c] = Alphabet[Complement[code]];
        }
        return map;
    }();

    void reverse_complement_inplace(uint64_t k) override {
        assert(k <= base::max_k);
        AlphabetKmer result;
        for (uint64_t i = 0; i < k; ++i) result.set(k - 1 - i, Complement[base::at(i)]);
        *this = result;
    }
    static uint64_t reverse_complement_mmer(uint64_t value, uint64_t m) {
        assert(m <= base::max_m);
        AlphabetKmer<uint64_t, Alphabet, Encode, Complement> packed(value);
        packed.reverse_complement_inplace(m);
        return uint64_t(packed);
    }
    static void compute_reverse_complement(char const* input, char* output, uint64_t n) {
        for (uint64_t i = 0; i < n; ++i)
            output[n - 1 - i] = canonicalize_basepair_reverse_map[static_cast<unsigned char>(input[i])];
    }
};

template <typename W>
using DNA5Kmer = AlphabetKmer<W, kmer::alphabets::kAlphabetDNA5,
                            kmer::alphabets::kCharToDNA, kmer::alphabets::kComplementMapDNA>;

template <typename W>
using CaseSensitiveDNAKmer = AlphabetKmer<W, kmer::alphabets::kAlphabetDNACaseSent,
                                        kmer::alphabets::kCharToDNACaseSent,
                                        kmer::alphabets::kComplementMapDNACaseSent>;

} // namespace mtg::graph::detail
