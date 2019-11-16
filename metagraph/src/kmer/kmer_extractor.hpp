#ifndef __KMER_EXTRACTOR_HPP__
#define __KMER_EXTRACTOR_HPP__

#include <cstdint>
#include <string>
#include <vector>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>
#include <sdsl/int_vector.hpp>

#include "utils/string_utils.hpp"
#include "utils/vectors.hpp"
#include "kmer.hpp"
#include "kmer_boss.hpp"
#include "alphabets.hpp"


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

    KmerExtractorBOSS();

    /**
     * Break the sequence into kmers and add them to the kmer storage.
     */
    template <class KMER>
    static void sequence_to_kmers(const std::string &sequence,
                                  size_t k,
                                  const std::vector<TAlphabet> &suffix,
                                  Vector<KMER> *kmers,
                                  bool canonical_mode = false);

    template <class KMER>
    static Vector<KMER> sequence_to_kmers(const std::string &sequence,
                                          size_t k,
                                          bool canonical_mode = false,
                                          const std::vector<TAlphabet> &suffix = {});

    static sdsl::bit_vector valid_kmers(const std::string &sequence, size_t k);

    template <class KMER>
    static std::string kmer_to_sequence(const KMER &kmer, size_t k) {
        return kmer.to_string(k, alphabet);
    }

    // map input character to k-mer character
    static TAlphabet encode(char s);
    static std::vector<TAlphabet> encode(const std::string &sequence);
    // map k-mer character to input character
    static char decode(TAlphabet c);
    static std::string decode(const std::vector<TAlphabet> &sequence);

    static std::vector<TAlphabet>
    reverse_complement(const std::vector<TAlphabet> &sequence);
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
    static const std::vector<TAlphabet> kComplementCode;
};


template <const uint8_t BitsPerChar>
class KmerExtractor2BitT {
  public:
    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    static constexpr uint8_t bits_per_char = BitsPerChar;
    const std::string alphabet;

    // k-mer
    template <class T>
    using Kmer = KMer<T, bits_per_char>;

    typedef Kmer<uint64_t> Kmer64;
    typedef Kmer<sdsl::uint128_t> Kmer128;
    typedef Kmer<sdsl::uint256_t> Kmer256;

    KmerExtractor2BitT(const char Alphabet[] = alphabets::kAlphabetDNA,
                       const uint8_t CharToCode[128] = alphabets::kCharToDNA,
                       const std::vector<uint8_t> &complement_code = alphabets::kComplementMapDNA);

    /**
     * Break the sequence into kmers and add them to the kmer storage.
     */
    template <class T>
    void sequence_to_kmers(const std::string &sequence,
                           size_t k,
                           const std::vector<TAlphabet> &suffix,
                           Vector<Kmer<T>> *kmers,
                           bool canonical_mode = false) const;
    template <class KMER>
    Vector<KMER> sequence_to_kmers(const std::string &sequence,
                                   size_t k,
                                   bool canonical_mode = false,
                                   const std::vector<TAlphabet> &suffix = {}) const;

    sdsl::bit_vector valid_kmers(const std::string &sequence, size_t k) const;

    template <class T>
    std::string kmer_to_sequence(const Kmer<T> &kmer, size_t k) const {
        return kmer.to_string(k, alphabet);
    }

    std::string reverse_complement(const std::string &sequence) const;

    // map input character to k-mer character
    TAlphabet encode(char s) const;
    std::vector<TAlphabet> encode(const std::string &sequence) const;
    // map k-mer character to input character
    char decode(TAlphabet c) const;
    std::string decode(const std::vector<TAlphabet> &sequence) const;

    std::vector<std::string> generate_suffixes(size_t len) const;

  private:
    const TAlphabet *char_to_code_;
    const std::vector<TAlphabet> complement_code_;
};

typedef KmerExtractor2BitT<alphabets::kBitsPerCharDNA> KmerExtractor2Bit;

#endif // __KMER_EXTRACTOR_HPP__
