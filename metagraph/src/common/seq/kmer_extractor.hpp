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
#include "alphabets.hpp"


class KmerExtractor {
  public:

    #if _PROTEIN_GRAPH
    static constexpr size_t kLogSigma = 5;
    #elif _DNA_CASE_SENSITIVE_GRAPH
    static constexpr size_t kLogSigma = 4;
    #elif _DNA_GRAPH
    static constexpr size_t kLogSigma = 3;
    #elif _DNA4_GRAPH
    static constexpr size_t kLogSigma = 3;
    #else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA4_GRAPH, _DNA_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
    #endif

    typedef KMer<uint64_t, kLogSigma> Kmer64;
    typedef KMer<sdsl::uint128_t, kLogSigma> Kmer128;
    typedef KMer<sdsl::uint256_t, kLogSigma> Kmer256;

    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    KmerExtractor();

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

    static const std::string alphabet;

  private:
    static const TAlphabet *kCharToNucleotide;
};


template <const uint8_t LogSigma>
class KmerExtractor2BitT {
  public:
    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;

    static constexpr uint8_t kLogSigma = LogSigma;
    const std::string alphabet;

    // k-mer
    template <class T>
    using Kmer = KMerPacked<T, kLogSigma>;

    typedef Kmer<uint64_t> Kmer64;
    typedef Kmer<sdsl::uint128_t> Kmer128;
    typedef Kmer<sdsl::uint256_t> Kmer256;

    KmerExtractor2BitT(const char Alphabet[] = alphabets::kAlphabetDNA4,
                       const uint8_t CharToCode[128] = alphabets::kCharToDNA4,
                       const std::vector<uint8_t> &complement_code = alphabets::kCanonicalMapDNA4);

    /**
     * Break the sequence into kmers and add them to the kmer storage.
     */
    template <class T>
    void sequence_to_kmers(const std::string &sequence,
                           size_t k,
                           const std::vector<TAlphabet> &suffix,
                           Vector<Kmer<T>> *kmers,
                           bool canonical_mode = false) const;
    template <class T>
    Vector<Kmer<T>> sequence_to_kmers(const std::string &sequence,
                                      size_t k,
                                      bool canonical_mode = false,
                                      const std::vector<TAlphabet> &suffix = {}) const;

    template <class T>
    std::string kmer_to_sequence(const Kmer<T> &kmer, size_t k) const {
        return kmer.to_string(k, alphabet);
    }

    // map input character to k-mer character
    TAlphabet encode(char s) const;
    std::vector<TAlphabet> encode(const std::string &sequence) const;
    // map k-mer character to input character
    char decode(TAlphabet c) const;
    std::string decode(const std::vector<TAlphabet> &sequence) const;

  private:
    const TAlphabet *char_to_code_;
    const std::vector<TAlphabet> complement_code_;
};


typedef KmerExtractor2BitT<alphabets::kLogSigmaDNA4> KmerExtractor2Bit;


#endif // __KMER_EXTRACTOR_HPP__
