#ifndef __SSHASH_ALPHABETS_HPP__
#define __SSHASH_ALPHABETS_HPP__

#include <array>
#include <kmer.hpp>

namespace mtg::graph::detail {

// Extended DNA alphabets. N is a literal symbol and complements itself.
inline constexpr char nucleotides5[] = "ACGTN";
inline constexpr char nucleotides_case[] = "ACGTNacgt";
template <typename Kmer, bool CaseSensitive>
struct ExtendedDNAKmer : sshash::alpha_kmer_t<Kmer, CaseSensitive ? 4 : 3,
                                            CaseSensitive ? nucleotides_case : nucleotides5> {
    using base = sshash::alpha_kmer_t<Kmer, CaseSensitive ? 4 : 3,
                                     CaseSensitive ? nucleotides_case : nucleotides5>;
    using base::base;
    static constexpr bool has_reverse_complement = true;

    static constexpr uint64_t char_to_uint(unsigned char c) {
        if constexpr (CaseSensitive) {
            switch (c) {
                case 'a': return 5;
                case 'c': return 6;
                case 'g': return 7;
                case 't': case 'u': return 8;
            }
        }
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': case 'U': case 'u': return 3;
            case 'N': case 'n': return 4;
            default: return uint64_t(-1);
        }
    }
    static constexpr bool is_valid(char c) {
        return char_to_uint(static_cast<unsigned char>(c)) != uint64_t(-1);
    }
    static constexpr uint64_t complement(uint64_t c) {
        return c < 4 ? 3 - c : (CaseSensitive && c > 4 ? 13 - c : 4);
    }

    // Streaming lookup uses this map when rolling the reverse complement.
    inline static constexpr auto canonicalize_basepair_reverse_map = [] {
        std::array<char, 256> map{};
        for (unsigned c = 0; c < 256; ++c) {
            auto code = char_to_uint(c);
            if (code != uint64_t(-1)) map[c] = base::alphabet[complement(code)];
        }
        return map;
    }();

    void reverse_complement_inplace(uint64_t k) override {
        assert(k <= base::max_k);
        ExtendedDNAKmer result;
        for (uint64_t i = 0; i < k; ++i) result.set(k - i - 1, complement(base::at(i)));
        *this = result;
    }
    static uint64_t reverse_complement_mmer(uint64_t value, uint64_t m) {
        assert(m <= base::max_m);
        uint64_t result = 0;
        for (uint64_t i = 0; i < m; ++i) {
            result = (result << base::bits_per_char)
                   | complement(value & ((uint64_t(1) << base::bits_per_char) - 1));
            value >>= base::bits_per_char;
        }
        return result;
    }
    static void compute_reverse_complement(char const* input, char* output, uint64_t size) {
        for (uint64_t i = 0; i < size; ++i)
            output[size - i - 1] = canonicalize_basepair_reverse_map[
                static_cast<unsigned char>(input[i])];
    }
};

template <typename Kmer>
using DNA5Kmer = ExtendedDNAKmer<Kmer, false>;

template <typename Kmer>
using CaseSensitiveDNAKmer = ExtendedDNAKmer<Kmer, true>;

} // namespace mtg::graph::detail

#endif // __SSHASH_ALPHABETS_HPP__
