#ifndef __SSHASH_ALPHABETS_HPP__
#define __SSHASH_ALPHABETS_HPP__

#include <array>
#include <kmer.hpp>

namespace mtg::graph::detail {

// Five-letter DNA. N is a literal symbol and complements itself.
inline constexpr char nucleotides5[] = "ACGTN";
template <typename Kmer>
struct DNA5Kmer : sshash::alpha_kmer_t<Kmer, 3, nucleotides5> {
    using base = sshash::alpha_kmer_t<Kmer, 3, nucleotides5>;
    using base::base;
    static constexpr bool has_reverse_complement = true;

    static constexpr uint64_t char_to_uint(unsigned char c) {
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
    static constexpr uint64_t complement(uint64_t c) { return c < 4 ? 3 - c : 4; }

    // Streaming lookup uses this map when rolling the reverse complement.
    inline static constexpr auto canonicalize_basepair_reverse_map = [] {
        std::array<char, 256> map{};
        for (unsigned c = 0; c < 256; ++c) {
            auto code = char_to_uint(c);
            if (code != uint64_t(-1)) map[c] = nucleotides5[complement(code)];
        }
        return map;
    }();

    void reverse_complement_inplace(uint64_t k) override {
        assert(k <= base::max_k);
        DNA5Kmer result;
        for (uint64_t i = 0; i < k; ++i) result.set(k - i - 1, complement(base::at(i)));
        *this = result;
    }
    static uint64_t reverse_complement_mmer(uint64_t value, uint64_t m) {
        assert(m <= base::max_m);
        uint64_t result = 0;
        for (uint64_t i = 0; i < m; ++i) {
            result = (result << 3) | complement(value & 7);
            value >>= 3;
        }
        return result;
    }
    static void compute_reverse_complement(char const* input, char* output, uint64_t size) {
        for (uint64_t i = 0; i < size; ++i)
            output[size - i - 1] = canonicalize_basepair_reverse_map[
                static_cast<unsigned char>(input[i])];
    }
};

} // namespace mtg::graph::detail

#endif // __SSHASH_ALPHABETS_HPP__
