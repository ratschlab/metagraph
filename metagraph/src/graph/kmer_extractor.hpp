#ifndef __KMER_EXTRACTOR_HPP__
#define __KMER_EXTRACTOR_HPP__

#include <cstdint>
#include <string>
#include <vector>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "utils.hpp"
#include "kmer.hpp"
#include "kmer_packed.hpp"


class KmerExtractor {
  public:

    #if _PROTEIN_GRAPH
    static constexpr size_t kLogSigma = 5;
    #elif _DNA_CASE_SENSITIVE_GRAPH
    static constexpr size_t kLogSigma = 4;
    #elif _DNA_GRAPH
    static constexpr size_t kLogSigma = 3;
    #else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
    #endif

    typedef KMer<uint64_t, kLogSigma> Kmer64;
    typedef KMer<sdsl::uint128_t, kLogSigma> Kmer128;
    typedef KMer<sdsl::uint256_t, kLogSigma> Kmer256;

    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    KmerExtractor();

    /**
     * Break the sequence to kmers and extend the temporary kmers storage.
     */
    template <class KMER>
    static void sequence_to_kmers(const std::string &sequence,
                                  size_t k,
                                  const std::vector<TAlphabet> &suffix,
                                  Vector<KMER> *kmers,
                                  bool canonical_mode = false);

    template <class KMER>
    static std::string kmer_to_sequence(const KMER &kmer, size_t k);

    // map input character to k-mer character
    static TAlphabet encode(char s);
    static std::vector<TAlphabet> encode(const std::string &sequence);
    // map k-mer character to input character
    static char decode(TAlphabet c);
    static std::string decode(const std::vector<TAlphabet> &sequence);

    static const std::string alphabet;

  private:
    static const TAlphabet kCharToNucleotide[128];
};


class KmerExtractor2Bit {
  public:
    static constexpr size_t kLogSigma = 2;

    typedef KMerPacked<uint64_t, kLogSigma> Kmer64;
    typedef KMerPacked<sdsl::uint128_t, kLogSigma> Kmer128;
    typedef KMerPacked<sdsl::uint256_t, kLogSigma> Kmer256;

    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    KmerExtractor2Bit();

    template <class KMER>
    static void sequence_to_kmers(const std::string &sequence,
                                  size_t k,
                                  const std::vector<TAlphabet> &suffix,
                                  Vector<KMER> *kmers,
                                  bool canonical_mode = false);

    template <class KMER>
    static std::string kmer_to_sequence(const KMER &kmer, size_t k);

    // map input character to k-mer character
    static TAlphabet encode(char s);
    static std::vector<TAlphabet> encode(const std::string &sequence);
    // map k-mer character to input character
    static char decode(TAlphabet c);
    static std::string decode(const std::vector<TAlphabet> &sequence);

    static const std::string alphabet;

  private:
    static const TAlphabet kCharToNucleotide[128];
};

#endif // __KMER_EXTRACTOR_HPP__
