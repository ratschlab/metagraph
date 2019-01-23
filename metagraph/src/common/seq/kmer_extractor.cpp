#include "kmer_extractor.hpp"

#include <algorithm>
#include <cstdlib>


namespace extractor {
/*
 * Helper functions
 */

template <typename TAlphabet>
inline char decode(TAlphabet c, const std::string &alphabet) {
    assert(c < alphabet.size());
    return alphabet[c];
}

template <typename TAlphabet>
inline std::string decode(const std::vector<TAlphabet> &sequence,
                          const std::string &alphabet) {
    std::string str(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(), str.begin(),
        [&alphabet](TAlphabet x) { return decode(x, alphabet); }
    );
    return str;
}

template <typename TAlphabet>
inline TAlphabet encode(char s, const TAlphabet kCharToNucleotide[]) {
    return s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                  : kCharToNucleotide[0];
}

template <typename TAlphabet>
inline std::vector<TAlphabet> encode(const std::string &sequence,
                                     const TAlphabet kCharToNucleotide[]) {
    std::vector<TAlphabet> seq_encoded(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(), seq_encoded.begin(),
        [&](char c) { return encode(c, kCharToNucleotide); }
    );
    return seq_encoded;
}


/**
 * Break the sequence into kmers and add them to the kmer storage.
 */
template <class KMER, typename TAlphabet>
inline void sequence_to_kmers(size_t k,
                              const std::vector<TAlphabet> &seq,
                              const std::vector<TAlphabet> &rev_comp,
                              const std::vector<TAlphabet> &suffix,
                              Vector<KMER> *kmers) {
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        if (!rev_comp.size()) {
            if (std::equal(suffix.begin(), suffix.end(), &seq[i + k - suffix.size()]))
                kmers->emplace_back(&seq[i], k);

            continue;
        }

        bool suffix_matched_forward = std::equal(suffix.begin(), suffix.end(),
                                                 &seq[i + k - suffix.size()]);
        bool suffix_matched_reverse = std::equal(suffix.begin(), suffix.end(),
                                                 &rev_comp[i + k - suffix.size()]);
        if (!suffix_matched_forward && !suffix_matched_reverse)
            continue;

        auto forward = KMER(&seq[i], k);
        auto reverse = KMER(&rev_comp[i], k);

        if (forward < reverse) {
            if (suffix_matched_forward)
                kmers->push_back(forward);
        } else {
            if (suffix_matched_reverse)
                kmers->push_back(reverse);
        }
    }
}

/**
 * Break the sequence into kmers and add them to the kmer storage.
 */
template <class KMER, typename TAlphabet>
inline void sequence_to_kmers_slide(size_t k,
                                    const std::vector<TAlphabet> &seq,
                                    const std::vector<TAlphabet> &rev_comp,
                                    const std::vector<TAlphabet> &suffix,
                                    Vector<KMER> *kmers) {
    // initialize and add the first kmer from sequence
    auto kmer = KMER(seq, k);
    KMER rev;
    if (rev_comp.size())
        rev = KMER(&rev_comp[seq.size() - k], k);

    if (std::equal(suffix.begin(), suffix.end(),
                   &seq[k - suffix.size()])) {
        if (rev_comp.empty() || rev >= kmer) {
            kmers->emplace_back(kmer);
        } else if (std::equal(suffix.begin(), suffix.end(),
                              &rev_comp[seq.size() - suffix.size()])) {
            kmers->emplace_back(rev);
        }
    } else if (rev_comp.size()
                    && rev < kmer
                    && std::equal(suffix.begin(), suffix.end(),
                                  &rev_comp[seq.size() - suffix.size()])) {
        kmers->emplace_back(rev);
    }

    // add all other kmers
    for (size_t i = 1; i + k <= seq.size(); ++i) {
        kmer.to_next(k, seq[i + k - 1]);//, seq[i + k - 2]);
        if (rev_comp.size())
            rev.to_prev(k, rev_comp[seq.size() - i - k]);

        if (std::equal(suffix.begin(), suffix.end(),
                       &seq[i + k - suffix.size()])) {
            if (rev_comp.empty() || rev >= kmer) {
                kmers->emplace_back(kmer);
            } else if (std::equal(suffix.begin(), suffix.end(),
                               &rev_comp[seq.size() - i - suffix.size()])) {
                kmers->emplace_back(rev);
            }
        } else if (rev_comp.size()
                        && rev < kmer
                        && std::equal(suffix.begin(), suffix.end(),
                                      &rev_comp[seq.size() - i - suffix.size()])) {
            kmers->emplace_back(rev);
        }
    }
}


/**
 * Break the sequence into kmers and add them to the kmer storage.
 */
template <class KMER, typename TAlphabet>
inline void sequence_to_kmers(const std::vector<TAlphabet> &seq,
                              size_t k,
                              const std::vector<TAlphabet> &suffix,
                              Vector<KMER> *kmers,
                              const std::vector<uint8_t> &canonical_map) {
    assert(k);
    assert(suffix.size() < k);
    assert(kmers);

    if (seq.size() < k)
        return;

    std::vector<TAlphabet> rev_comp;
    if (canonical_map.size()) {
        rev_comp.reserve(seq.size());
        std::transform(seq.rbegin(), seq.rend(), std::back_inserter(rev_comp),
                       [&](const auto c) -> TAlphabet { return canonical_map.at(c); });
    }

    // based on performance comparison
    // for KMer::pack_kmer and KMer::update_kmer
    if (suffix.size() > 1) {
        sequence_to_kmers(k, seq, rev_comp, suffix, kmers);
    } else {
        sequence_to_kmers_slide(k, seq, rev_comp, suffix, kmers);
    }
}

} // namespace extractor


/**
 * KmerExtractor
 */
#if _PROTEIN_GRAPH
    const std::string KmerExtractor::alphabet = alphabets::kBOSSAlphabetProtein;
    const KmerExtractor::TAlphabet *KmerExtractor::kCharToNucleotide = alphabets::kBOSSCharToProtein;
    const std::vector<uint8_t> canonical_map = alphabets::kBOSSCanonicalMapProtein;
#elif _DNA_CASE_SENSITIVE_GRAPH
    const std::string KmerExtractor::alphabet = alphabets::kBOSSAlphabetDNACaseSent;
    const KmerExtractor::TAlphabet *KmerExtractor::kCharToNucleotide = alphabets::kBOSSCharToDNACaseSent;
    const std::vector<uint8_t> canonical_map = alphabets::kBOSSCanonicalMapDNACaseSent;
#elif _DNA_GRAPH
    const std::string KmerExtractor::alphabet = alphabets::kBOSSAlphabetDNA;
    const KmerExtractor::TAlphabet *KmerExtractor::kCharToNucleotide = alphabets::kBOSSCharToDNA;
    const std::vector<uint8_t> canonical_map = alphabets::kBOSSCanonicalMapDNA;
#elif _DNA4_GRAPH
    const std::string KmerExtractor::alphabet = alphabets::kBOSSAlphabetDNA4;
    const KmerExtractor::TAlphabet *KmerExtractor::kCharToNucleotide = alphabets::kBOSSCharToDNA4;
    const std::vector<uint8_t> canonical_map = alphabets::kBOSSCanonicalMapDNA4;
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA4_GRAPH, _DNA_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

static_assert(KmerExtractor::kLogSigma <= sizeof(KmerExtractor::TAlphabet) * 8,
              "Choose type for TAlphabet properly");


KmerExtractor::KmerExtractor() {
    assert(alphabet.size() <= (1llu << kLogSigma));
}

KmerExtractor::TAlphabet KmerExtractor::encode(char s) {
    assert(extractor::encode(s, kCharToNucleotide) < alphabet.size());
    return extractor::encode(s, kCharToNucleotide);
}

char KmerExtractor::decode(TAlphabet c) {
    return extractor::decode(c, alphabet);
}

std::vector<KmerExtractor::TAlphabet>
KmerExtractor::encode(const std::string &sequence) {
    #ifndef NDEBUG
    auto encoded = extractor::encode(sequence, kCharToNucleotide);
    assert(std::all_of(encoded.begin(), encoded.end(),
        [&](TAlphabet c) { return c < alphabet.size(); }
    ));
    #endif
    return extractor::encode(sequence, kCharToNucleotide);
}

std::string KmerExtractor::decode(const std::vector<TAlphabet> &sequence) {
    return extractor::decode(sequence, alphabet);
}

/**
 * Break the sequence into kmers and add them to the kmer storage.
 */
template <typename KMER>
void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers,
                                      bool canonical_mode) {
    assert(k);
    // suffix does not include the last character
    assert(suffix.size() < k);

    if (sequence.size() < k)
        return;

    // encode sequence
    size_t dummy_prefix_size = suffix.size() > 0 ? k - 1 : 1;

    std::vector<TAlphabet> seq(sequence.size() + dummy_prefix_size + 1, 0);

    std::transform(sequence.begin(), sequence.end(), &seq[dummy_prefix_size],
        [](char c) { return encode(c); }
    );

    extractor::sequence_to_kmers(seq, k, suffix, kmers,
        canonical_mode ? canonical_map : std::vector<uint8_t>()
    );
}

template
void KmerExtractor::sequence_to_kmers(const std::string&,
                                      size_t,
                                      const std::vector<TAlphabet>&,
                                      Vector<Kmer64>*,
                                      bool);
template
void KmerExtractor::sequence_to_kmers(const std::string&,
                                      size_t,
                                      const std::vector<TAlphabet>&,
                                      Vector<Kmer128>*,
                                      bool);
template
void KmerExtractor::sequence_to_kmers(const std::string&,
                                      size_t,
                                      const std::vector<TAlphabet>&,
                                      Vector<Kmer256>*,
                                      bool);


/**
 * KmerExtractor2Bit
 */

#define KmerExtractor2BitTDecl(...) \
template <const uint8_t LogSigma> \
__VA_ARGS__ KmerExtractor2BitT<LogSigma>

KmerExtractor2BitTDecl()
::KmerExtractor2BitT(const char alph[],
                     const uint8_t char_to_code[128],
                     const std::vector<uint8_t> &complement_code)
      : alphabet(alph),
        char_to_code_(char_to_code),
        complement_code_(complement_code) {
    static_assert(kLogSigma <= sizeof(TAlphabet) * 8,
                  "Choose type for TAlphabet properly");

    assert(alphabet.size() <= (1llu << kLogSigma));
}

KmerExtractor2BitTDecl(KmerExtractor::TAlphabet)
::encode(char s) const {
    assert(extractor::encode(s, char_to_code_) < alphabet.size());
    return extractor::encode(s, char_to_code_);
}

KmerExtractor2BitTDecl(char)
::decode(TAlphabet c) const {
    return extractor::decode(c, alphabet);
}

KmerExtractor2BitTDecl(std::vector<KmerExtractor2Bit::TAlphabet>)
::encode(const std::string &sequence) const {
    #ifndef NDEBUG
    auto encoded = extractor::encode(sequence, char_to_code_);
    assert(std::all_of(encoded.begin(), encoded.end(),
        [&](TAlphabet c) { return c < alphabet.size(); }
    ));
    #endif
    return extractor::encode(sequence, char_to_code_);
}

KmerExtractor2BitTDecl(std::string)
::decode(const std::vector<TAlphabet> &sequence) const {
    return extractor::decode(sequence, alphabet);
}

/**
 * Break the sequence into kmers and add them to the kmer storage.
 */
KmerExtractor2BitTDecl(template <typename T> void)
::sequence_to_kmers(const std::string &sequence,
                    size_t k,
                    const std::vector<TAlphabet> &suffix,
                    Vector<Kmer<T>> *kmers,
                    bool canonical_mode) const {
    assert(k);
    // suffix does not include the last character
    assert(suffix.size() < k);

    if (sequence.size() < k)
        return;

    extractor::sequence_to_kmers(
        encode(sequence), k, suffix, kmers,
        canonical_mode ? complement_code_
                       : std::vector<TAlphabet>()
    );
}

template class KmerExtractor2BitT<alphabets::kLogSigmaDNA4>;

#define ExplicitInstantiation_sequence_to_kmers(T) \
template \
void KmerExtractor2Bit \
::sequence_to_kmers<T>(const std::string&, \
                              size_t, \
                              const std::vector<TAlphabet>&, \
                              Vector<Kmer<T>>*, \
                              bool) const;

ExplicitInstantiation_sequence_to_kmers(uint64_t)
ExplicitInstantiation_sequence_to_kmers(sdsl::uint128_t)
ExplicitInstantiation_sequence_to_kmers(sdsl::uint256_t)
