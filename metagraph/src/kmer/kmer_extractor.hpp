#ifndef __KMER_EXTRACTOR_HPP__
#define __KMER_EXTRACTOR_HPP__

#include <cstdint>
#include <string>
#include <vector>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "common/vector.hpp"
#include "common/utils/string_utils.hpp"
#include "kmer.hpp"
#include "kmer_boss.hpp"
#include "alphabets.hpp"


namespace mtg {
namespace kmer {

/**
 * Extracts k-mers to be placed in a BOSS table from sequences. Each sequence is
 * prepended with (k-1) $ characters (aka 'dummy prefix') and terminated with a $ to
 * facilitate graph traversal.
 * A sequence may sometimes contain multiple segments separated by a character not in
 * the alphabet (DNA sequences are typically separated by an N, e.g. AAATNGGGC
 * resulting in the segments AAAT and GGGC).
 */
class KmerExtractorBOSS {
  public:

    #if _PROTEIN_GRAPH
    static constexpr size_t bits_per_char = alphabets::kBOSSBitsPerCharProtein;
    #elif _DNA_CASE_SENSITIVE_GRAPH
    static constexpr size_t bits_per_char = alphabets::kBOSSBitsPerCharDNACaseSent;
    #elif _DNA5_GRAPH
    static constexpr size_t bits_per_char = alphabets::kBOSSBitsPerCharDNA5;
    #elif _DNA_GRAPH
    static constexpr size_t bits_per_char = alphabets::kBOSSBitsPerCharDNA;
    #else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
    #endif

    typedef KMerBOSS<uint64_t, bits_per_char> Kmer64;
    typedef KMerBOSS<sdsl::uint128_t, bits_per_char> Kmer128;
    typedef KMerBOSS<sdsl::uint256_t, bits_per_char> Kmer256;

    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    static const std::vector<TAlphabet> kComplementCode;

    KmerExtractorBOSS();

    /**
     * Break the sequence into kmers and add them to the kmer storage.
     * Adds only valid k-mers.
     */
    template <class KMER>
    static void sequence_to_kmers(std::string_view sequence,
                                  size_t k,
                                  const std::vector<TAlphabet> &suffix,
                                  Vector<KMER> *kmers,
                                  bool canonical_mode = false);

    /**
     * Basic conversion of a sequence to a k-mer. Sentinel characters in the
     * input sequence are also encoded.
     */
    template <class KMER>
    static KMER sequence_to_kmer(std::string_view sequence, bool canonical_mode = false);

    template <class KMER>
    static std::string kmer_to_sequence(const KMER &kmer, size_t k) {
        return kmer.to_string(k, alphabet);
    }

    // map input character to k-mer character
    static TAlphabet encode(char s);
    static std::vector<TAlphabet> encode(std::string_view sequence);
    // map k-mer character to input character
    static char decode(TAlphabet c);
    static std::string decode(const std::vector<TAlphabet> &sequence);

    static void reverse_complement(std::vector<TAlphabet> *sequence);
    static TAlphabet complement(TAlphabet c);

    /**
     * Generate all valid suffixes of the given length
     * ACGT, ACG$, ..., $$$$ -- valid
     * AC$T, A$$T, ..., $AAA -- invalid
     */
    static std::vector<std::string> generate_suffixes(size_t len);

    static const std::string alphabet;

  private:
    static const TAlphabet *kCharToNucleotide;
};


template <const uint8_t BitsPerChar>
class KmerExtractorT {
  public:
    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    static constexpr uint8_t bits_per_char = BitsPerChar;
    const std::string alphabet;

    // k-mers
    typedef KMer<uint64_t, bits_per_char> Kmer64;
    typedef KMer<sdsl::uint128_t, bits_per_char> Kmer128;
    typedef KMer<sdsl::uint256_t, bits_per_char> Kmer256;
    // k-mers with the BOSS layout
    typedef KMerBOSS<uint64_t, bits_per_char> KmerBOSS64;
    typedef KMerBOSS<sdsl::uint128_t, bits_per_char> KmerBOSS128;
    typedef KMerBOSS<sdsl::uint256_t, bits_per_char> KmerBOSS256;

    #if _PROTEIN_GRAPH
    KmerExtractorT(const char Alphabet[] = alphabets::kAlphabetProtein,
                   const uint8_t CharToCode[128] = alphabets::kCharToProtein,
                   const std::vector<uint8_t> &complement_code = alphabets::kComplementMapProtein);
    #elif _DNA_CASE_SENSITIVE_GRAPH
    KmerExtractorT(const char Alphabet[] = alphabets::kAlphabetDNACaseSent,
                   const uint8_t CharToCode[128] = alphabets::kCharToDNACaseSent,
                   const std::vector<uint8_t> &complement_code = alphabets::kComplementMapDNACaseSent);
    #elif _DNA5_GRAPH
    KmerExtractorT(const char Alphabet[] = alphabets::kAlphabetDNA5,
                   const uint8_t CharToCode[128] = alphabets::kCharToDNA,
                   const std::vector<uint8_t> &complement_code = alphabets::kComplementMapDNA);
    #elif _DNA_GRAPH
    KmerExtractorT(const char Alphabet[] = alphabets::kAlphabetDNA,
                   const uint8_t CharToCode[128] = alphabets::kCharToDNA,
                   const std::vector<uint8_t> &complement_code = alphabets::kComplementMapDNA);
    #else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
    #endif

    /**
     * Break the sequence into kmers and add them to the kmer collector. If suffix is
     * not empty, only kmers with the given suffix are added.
     */
    template <class KMER>
    void sequence_to_kmers(std::string_view sequence,
                           size_t k,
                           const std::vector<TAlphabet> &suffix,
                           Vector<KMER> *kmers,
                           bool canonical_mode = false) const;

    /**
     * Extract all k-mers from sequence.
     * Returned pairs are k-mers and flags: `true` for the valid k-mers
     * and `false` for invalid ones.
     */
    template <class KMER>
    Vector<std::pair<KMER, bool>>
    sequence_to_kmers(std::string_view sequence,
                      size_t k,
                      bool canonical_mode = false,
                      const std::vector<TAlphabet> &suffix = {}) const;

    template <class KMER>
    std::string kmer_to_sequence(const KMER &kmer, size_t k) const {
        return kmer.to_string(k, alphabet);
    }

    // map input character to k-mer character
    TAlphabet encode(char s) const;
    std::vector<TAlphabet> encode(std::string_view sequence) const;
    // map k-mer character to input character
    char decode(TAlphabet c) const;
    std::string decode(const std::vector<TAlphabet> &sequence) const;

    std::vector<std::string> generate_suffixes(size_t len) const;

    const std::vector<TAlphabet>& complement_code() {
        return complement_code_;
    }

  private:
    const TAlphabet *char_to_code_;
    const std::vector<TAlphabet> complement_code_;
};

#if _PROTEIN_GRAPH
typedef KmerExtractorT<alphabets::kBitsPerCharProtein> KmerExtractor2Bit;
#elif _DNA_CASE_SENSITIVE_GRAPH
typedef KmerExtractorT<alphabets::kBitsPerCharDNACaseSent> KmerExtractor2Bit;
#elif _DNA5_GRAPH
typedef KmerExtractorT<alphabets::kBitsPerCharDNA5> KmerExtractor2Bit;
#elif _DNA_GRAPH
typedef KmerExtractorT<alphabets::kBitsPerCharDNA> KmerExtractor2Bit;
#else
static_assert(false,
    "Define an alphabet: either "
    "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
);
#endif

} // namespace kmer
} // namespace mtg

#endif // __KMER_EXTRACTOR_HPP__
