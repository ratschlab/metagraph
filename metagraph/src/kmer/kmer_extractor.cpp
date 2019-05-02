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
inline std::vector<TAlphabet> encode(std::string::const_iterator begin,
                                     std::string::const_iterator end,
                                     const TAlphabet kCharToNucleotide[]) {
    std::vector<TAlphabet> seq_encoded(end - begin, 0);
    std::transform(begin, end, seq_encoded.begin(),
        [&](char c) { return encode(c, kCharToNucleotide); }
    );
    return seq_encoded;
}

template <typename TAlphabet>
std::vector<TAlphabet>
inline reverse_complement(const std::vector<TAlphabet> &sequence,
                          const std::vector<uint8_t> &canonical_map) {
    std::vector<TAlphabet> rev_comp(sequence.size());
    std::transform(sequence.rbegin(), sequence.rend(), rev_comp.begin(),
                   [&](const auto c) -> TAlphabet { return canonical_map.at(c); });
    return rev_comp;
}


/*
 * k-mer extractors
 */

template <class KMER, typename TAlphabet>
inline void sequence_to_kmers(size_t k,
                              const std::vector<TAlphabet> &seq,
                              const std::vector<TAlphabet> &suffix,
                              Vector<KMER> *kmers) {
    assert(kmers);
    assert(seq.size() >= k);

    for (size_t i = 0; i + k <= seq.size(); ++i) {
        if (KMER::match_suffix(&seq[i], k, suffix))
            kmers->emplace_back(&seq[i], k);
    }
}

template <class KMER, typename TAlphabet>
inline void sequence_to_kmers_slide(size_t k,
                                    const std::vector<TAlphabet> &seq,
                                    const std::vector<TAlphabet> &suffix,
                                    Vector<KMER> *kmers) {
    // initialize and add the first kmer from sequence
    assert(kmers);
    assert(seq.size() >= k);

    KMER kmer(seq.data(), k);
    if (KMER::match_suffix(seq.data(), k, suffix))
        kmers->push_back(kmer);

    // add all other kmers
    for (size_t i = k; i < seq.size(); ++i) {
        kmer.to_next(k, seq[i]);
        if (KMER::match_suffix(&seq[i - k + 1], k, suffix))
            kmers->push_back(kmer);
    }
}

// extract canonical (lexicographically smallest) k-mers
template <class KMER, typename TAlphabet>
inline void sequence_to_kmers_canonical(size_t k,
                                        const std::vector<TAlphabet> &seq,
                                        const std::vector<TAlphabet> &rev_comp,
                                        const std::vector<TAlphabet> &suffix,
                                        Vector<KMER> *kmers) {
    assert(kmers);
    assert(seq.size() >= k);
    assert(seq.size() == rev_comp.size());

    for (size_t i = 0; i + k <= seq.size(); ++i) {
        bool suffix_matched_forward = KMER::match_suffix(&seq[i], k, suffix);
        bool suffix_matched_reverse = KMER::match_suffix(&rev_comp[rev_comp.size() - i - k], k, suffix);

        if (!suffix_matched_forward && !suffix_matched_reverse)
            continue;

        KMER forward(&seq[i], k);
        KMER reverse(&rev_comp[rev_comp.size() - i - k], k);

        if (forward <= reverse) {
            if (suffix_matched_forward)
                kmers->push_back(forward);
        } else {
            if (suffix_matched_reverse)
                kmers->push_back(reverse);
        }
    }
}

template <class KMER, typename TAlphabet>
inline void push_back_smallest(const KMER &kmer_first, const TAlphabet *first,
                               const KMER &kmer_second, const TAlphabet *second,
                               size_t k,
                               const std::vector<TAlphabet> &suffix,
                               Vector<KMER> *kmers) {
    if (kmer_first <= kmer_second) {
        if (KMER::match_suffix(first, k, suffix))
            kmers->push_back(kmer_first);
    } else {
        if (KMER::match_suffix(second, k, suffix))
            kmers->push_back(kmer_second);
    }
}

// extract canonical (lexicographically smallest) k-mers
template <class KMER, typename TAlphabet>
inline void sequence_to_kmers_canonical_slide(size_t k,
                                              const std::vector<TAlphabet> &seq,
                                              const std::vector<TAlphabet> &rev_comp,
                                              const std::vector<TAlphabet> &suffix,
                                              Vector<KMER> *kmers) {
    assert(kmers);
    assert(seq.size() >= k);
    assert(seq.size() == rev_comp.size());

    // initialize and add the first kmer from sequence
    KMER kmer(seq.data(), k);
    KMER rev(&rev_comp[seq.size() - k], k);

    push_back_smallest(kmer, &seq[0],
                       rev, &rev_comp[rev_comp.size() - k],
                       k,
                       suffix, kmers);

    // add all other kmers
    for (size_t forward_last = k, reverse_first = seq.size() - k - 1;
                                        forward_last < seq.size();
                                            ++forward_last, --reverse_first) {
        kmer.to_next(k, seq[forward_last]);
        rev.to_prev(k, rev_comp[reverse_first]);

        push_back_smallest(kmer, &seq[forward_last - k + 1],
                           rev, &rev_comp[reverse_first],
                           k,
                           suffix, kmers);
    }
}


/**
 * Break the sequence into k-mers and add them to the kmer storage.
 */

template <class KMER, typename TAlphabet>
inline void sequence_to_kmers(const std::vector<TAlphabet> &seq,
                              size_t k,
                              const std::vector<TAlphabet> &suffix,
                              Vector<KMER> *kmers,
                              const std::vector<uint8_t> &canonical_map) {
    assert(kmers);
    assert(suffix.size() <= k);

    if (seq.size() < k)
        return;

    if (canonical_map.empty()) {
        // based on performance comparison
        // for KMer::pack_kmer and KMer::update_kmer
        if (suffix.size() > 1) {
            sequence_to_kmers(k, seq, suffix, kmers);
        } else {
            sequence_to_kmers_slide(k, seq, suffix, kmers);
        }
    } else {
        auto rev_comp = reverse_complement(seq, canonical_map);
        if (suffix.size() > 1) {
            sequence_to_kmers_canonical(k, seq, rev_comp, suffix, kmers);
        } else {
            sequence_to_kmers_canonical_slide(k, seq, rev_comp, suffix, kmers);
        }
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
    auto encoded = extractor::encode(sequence.begin(), sequence.end(), kCharToNucleotide);
    assert(std::all_of(encoded.begin(), encoded.end(),
        [&](TAlphabet c) { return c < alphabet.size(); }
    ));
    #endif
    return extractor::encode(sequence.begin(), sequence.end(), kCharToNucleotide);
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
    assert(kmers);
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

template <typename KMER>
Vector<KMER> KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                              size_t k,
                                              bool canonical_mode,
                                              const std::vector<TAlphabet> &suffix) {
    Vector<KMER> kmers;

    if (sequence.length() < k)
        return kmers;

    kmers.reserve(sequence.length() + 1 - k);
    sequence_to_kmers(sequence, k, suffix, &kmers, canonical_mode);
    return kmers;
}

template
Vector<KmerExtractor::Kmer64>
KmerExtractor::sequence_to_kmers(const std::string&,
                                 size_t,
                                 bool,
                                 const std::vector<TAlphabet>&);
template
Vector<KmerExtractor::Kmer128>
KmerExtractor::sequence_to_kmers(const std::string&,
                                 size_t,
                                 bool,
                                 const std::vector<TAlphabet>&);
template
Vector<KmerExtractor::Kmer256>
KmerExtractor::sequence_to_kmers(const std::string&,
                                 size_t,
                                 bool,
                                 const std::vector<TAlphabet>&);


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
    auto encoded = extractor::encode(sequence.begin(), sequence.end(), char_to_code_);
    assert(std::all_of(encoded.begin(), encoded.end(),
        [&](TAlphabet c) { return c < alphabet.size(); }
    ));
    #endif
    return extractor::encode(sequence.begin(), sequence.end(), char_to_code_);
}

KmerExtractor2BitTDecl(std::string)
::decode(const std::vector<TAlphabet> &sequence) const {
    return extractor::decode(sequence, alphabet);
}

KmerExtractor2BitTDecl(std::string)
::reverse_complement(const std::string &sequence) const {
    return decode(extractor::reverse_complement(encode(sequence), complement_code_));
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
    assert(kmers);
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

KmerExtractor2BitTDecl(template <typename KMER> Vector<KMER>)
::sequence_to_kmers(const std::string &sequence,
                    size_t k,
                    bool canonical_mode,
                    const std::vector<TAlphabet> &suffix) const {
    Vector<KMER> kmers;

    if (sequence.length() < k)
        return kmers;

    kmers.reserve(sequence.length() + 1 - k);
    sequence_to_kmers(sequence, k, suffix, &kmers, canonical_mode);
    return kmers;
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

#define ExplicitInstantiation_sequence_to_kmers_vector(T) \
template \
Vector<KmerExtractor2Bit::Kmer<T>> KmerExtractor2Bit \
::sequence_to_kmers<KmerExtractor2Bit::Kmer<T>>(const std::string&, \
                                                size_t, \
                                                bool, \
                                                const std::vector<TAlphabet>&) const;

ExplicitInstantiation_sequence_to_kmers_vector(uint64_t)
ExplicitInstantiation_sequence_to_kmers_vector(sdsl::uint128_t)
ExplicitInstantiation_sequence_to_kmers_vector(sdsl::uint256_t)
