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

class DBG_succ;


class KmerExtractor {
    friend class ::DBG_succ;

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
                                  Vector<KMER> *kmers);

    // extract k-mers from sequence
    template <class KMER>
    static void sequence_to_kmers(const std::vector<TAlphabet> &sequence,
                                  size_t k,
                                  const std::vector<TAlphabet> &suffix,
                                  Vector<KMER> *kmers);

    template <class KMER>
    static std::string kmer_to_sequence(const KMER &kmer);

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
    static constexpr size_t kLogSigma = 2;

  public:
    // alphabet for k-mer representation
    typedef uint8_t TAlphabet;
    // kmer type
    typedef KMerPacked<uint64_t, kLogSigma> Kmer;

    KmerExtractor2Bit();

    // extract k-mers from sequence
    Vector<Kmer> sequence_to_kmers(const std::string &sequence,
                                   size_t k) const;

    std::string kmer_to_sequence(const Kmer &kmer) const;

    // map input character to k-mer character
    TAlphabet encode(char s) const;
    // map k-mer character to input character
    char decode(TAlphabet c) const;

  private:
    std::vector<TAlphabet> encode(const std::string &sequence) const;
    std::string decode(const std::vector<TAlphabet> &sequence) const;

    Vector<Kmer>
    sequence_to_kmers(std::vector<TAlphabet>&& seq, size_t k) const;

    static const std::string alphabet;
    static const TAlphabet kCharToNucleotide[128];
};


#endif // __KMER_EXTRACTOR_HPP__
