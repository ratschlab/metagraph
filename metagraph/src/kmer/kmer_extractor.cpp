#include "kmer_extractor.hpp"

#include <algorithm>
#include <cstdlib>


namespace extractor {

/*
 * Helper functions
 */

template <typename TAlphabet>
inline char decode(TAlphabet c, const std::string &alphabet) {
    return alphabet.at(c);
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
inline std::vector<TAlphabet> encode(std::string::const_iterator __restrict begin,
                                     std::string::const_iterator __restrict end,
                                     const TAlphabet kCharToNucleotide[]) {
    std::vector<TAlphabet> seq_encoded(end - begin, 0);
    std::transform(begin, end, seq_encoded.begin(),
        [&](char c) { return encode(c, kCharToNucleotide); }
    );
    return seq_encoded;
}

template <typename TAlphabet>
std::vector<TAlphabet>
inline reverse_complement(const TAlphabet * __restrict begin,
                          const TAlphabet * __restrict end,
                          const std::vector<uint8_t> &canonical_map) {
    assert(end >= begin);
    std::vector<TAlphabet> rev_comp(end - begin);
    std::transform(begin, end, rev_comp.rbegin(),
                   [&](const auto c) -> TAlphabet { return canonical_map.at(c); });
    return rev_comp;
}


/*
 * k-mer extractors
 */

template <class KMER, typename TAlphabet>
inline void __sequence_to_kmers(const TAlphabet * __restrict begin,
                                const TAlphabet * __restrict end,
                                size_t k,
                                const std::vector<TAlphabet> &suffix,
                                Vector<KMER> *kmers) {
    assert(kmers);
    assert(end >= begin);
    // done to ensure that begin + k doesn't overflow
    assert(k <= static_cast<size_t>(end - begin));

    for (auto it = begin; it + k <= end; ++it) {
        if (KMER::match_suffix(it, k, suffix))
            kmers->emplace_back(it, k);
    }
}

template <class KMER, typename TAlphabet>
inline void __sequence_to_kmers_slide(const TAlphabet * __restrict begin,
                                      const TAlphabet * __restrict end,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers) {
    // initialize and add the first kmer from sequence
    assert(kmers);
    assert(end >= begin);
    // done to ensure that begin + k doesn't overflow
    assert(k <= static_cast<size_t>(end - begin));

    KMER kmer(begin, k);
    if (KMER::match_suffix(begin, k, suffix))
        kmers->push_back(kmer);

    // add all other kmers
    for (auto it = begin + 1; it + k <= end; ++it) {
        kmer.to_next(k, *(it + (k - 1)));
        if (KMER::match_suffix(it, k, suffix))
            kmers->push_back(kmer);
    }
}

// extract canonical (lexicographically smallest) k-mers
template <class KMER, typename TAlphabet>
inline void __sequence_to_kmers_canonical(const TAlphabet * __restrict begin,
                                          const TAlphabet * __restrict end,
                                          const std::vector<TAlphabet> &rev_comp,
                                          size_t k,
                                          const std::vector<TAlphabet> &suffix,
                                          Vector<KMER> *kmers) {
    assert(kmers);
    assert(end - begin);
    assert(rev_comp.size() == static_cast<size_t>(end - begin));
    assert(k <= rev_comp.size());
    std::ignore = end;

    for (size_t i = 0; i + k <= rev_comp.size(); ++i) {
        bool suffix_matched_forward = KMER::match_suffix(begin + i, k, suffix);
        bool suffix_matched_reverse = KMER::match_suffix(&rev_comp[rev_comp.size() - i - k], k, suffix);

        if (!suffix_matched_forward && !suffix_matched_reverse)
            continue;

        KMER forward(begin + i, k);
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
inline void __push_back_smallest(const KMER &kmer_first,
                                 const TAlphabet * __restrict first,
                                 const KMER &kmer_second,
                                 const TAlphabet * __restrict second,
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
inline void __sequence_to_kmers_canonical_slide(const TAlphabet * __restrict begin,
                                                const TAlphabet * __restrict end,
                                                const std::vector<TAlphabet> &rev_comp,
                                                size_t k,
                                                const std::vector<TAlphabet> &suffix,
                                                Vector<KMER> *kmers) {
    assert(kmers);
    assert(end - begin);
    assert(rev_comp.size() == static_cast<size_t>(end - begin));
    assert(k <= rev_comp.size());
    std::ignore = end;

    // initialize and add the first kmer from sequence
    KMER kmer(begin, k);
    KMER rev(&rev_comp[rev_comp.size() - k], k);

    __push_back_smallest(kmer, begin,
                         rev, &rev_comp[rev_comp.size() - k],
                         k,
                         suffix, kmers);

    // add all other kmers
    for (size_t forward_last = k,
                reverse_first = rev_comp.size() - k - 1;
                                        forward_last < rev_comp.size();
                                            ++forward_last, --reverse_first) {
        kmer.to_next(k, begin[forward_last]);
        rev.to_prev(k, rev_comp[reverse_first]);

        __push_back_smallest(kmer, begin + forward_last - (k - 1),
                             rev, &rev_comp[reverse_first],
                             k,
                             suffix, kmers);
    }
}


/**
 * Break the sequence into k-mers and add them to the kmer storage.
 */

template <class KMER, typename TAlphabet>
inline void sequence_to_kmers(const TAlphabet * __restrict begin,
                              const TAlphabet * __restrict end,
                              size_t k,
                              const std::vector<TAlphabet> &suffix,
                              Vector<KMER> *kmers,
                              const std::vector<uint8_t> &canonical_map) {
    assert(k);
    assert(kmers);
    assert(suffix.size() <= k);
    assert(end >= begin);

    // done to avoid begin + k overflowing
    if (k > static_cast<size_t>(end - begin))
        return;

    if (canonical_map.empty()) {
        // based on performance comparison
        // for KMer::pack_kmer and KMer::update_kmer
        if (suffix.size() > 1) {
            __sequence_to_kmers(begin, end, k, suffix, kmers);
        } else {
            __sequence_to_kmers_slide(begin, end, k, suffix, kmers);
        }
    } else {
        auto rev_comp = reverse_complement(begin, end, canonical_map);
        if (suffix.size() > 1) {
            __sequence_to_kmers_canonical(begin, end, rev_comp, k, suffix, kmers);
        } else {
            __sequence_to_kmers_canonical_slide(begin, end, rev_comp, k, suffix, kmers);
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
#elif _DNA5_GRAPH
    const std::string KmerExtractor::alphabet = alphabets::kBOSSAlphabetDNA5;
    const KmerExtractor::TAlphabet *KmerExtractor::kCharToNucleotide = alphabets::kBOSSCharToDNA;
    const std::vector<uint8_t> canonical_map = alphabets::kBOSSCanonicalMapDNA;
#elif _DNA_GRAPH
    const std::string KmerExtractor::alphabet = alphabets::kBOSSAlphabetDNA;
    const KmerExtractor::TAlphabet *KmerExtractor::kCharToNucleotide = alphabets::kBOSSCharToDNA;
    const std::vector<uint8_t> canonical_map = alphabets::kBOSSCanonicalMapDNA;
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

static_assert(KmerExtractor::bits_per_char <= sizeof(KmerExtractor::TAlphabet) * 8,
              "Choose type for TAlphabet properly");


KmerExtractor::KmerExtractor() {
    assert(alphabet.size() <= (1llu << bits_per_char));
}

KmerExtractor::TAlphabet KmerExtractor::encode(char s) {
    return extractor::encode(s, kCharToNucleotide);
}

char KmerExtractor::decode(TAlphabet c) {
    return extractor::decode(c, alphabet);
}

std::vector<KmerExtractor::TAlphabet>
KmerExtractor::encode(const std::string &sequence) {
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
    const size_t dummy_prefix_size = suffix.size() > 0 ? k - 1 : 1;

    std::vector<TAlphabet> seq(sequence.size() + dummy_prefix_size + 1, 0);

    std::transform(sequence.begin(), sequence.end(), &seq[dummy_prefix_size],
        [](char c) { return encode(c); }
    );

    TAlphabet * __restrict begin_segm = seq.data();
    TAlphabet * __restrict end_segm = seq.data() + dummy_prefix_size;
    TAlphabet * __restrict end = seq.data() + seq.size();

    while (begin_segm + dummy_prefix_size + k + 1 <= end) {
        assert(end >= end_segm);
        end_segm = std::find_if(end_segm, end,
            [&](auto c) { return c >= alphabet.size(); }
        );

        // ****AAANAAAA* ->
        // ****AAA*AAAA*
        if (end_segm < end) {
            assert(*end_segm >= alphabet.size());
            ++end_segm;
        }

        if (begin_segm + dummy_prefix_size + k + 1 <= end_segm) {
            // set dummy prefix for the next segment
            // ****AAA*AAAA*
            // ********AAAA*
            std::fill(begin_segm, begin_segm + dummy_prefix_size, 0);
            *(end_segm - 1) = 0;
            extractor::sequence_to_kmers(begin_segm, end_segm, k, suffix, kmers,
                canonical_mode ? canonical_map : std::vector<uint8_t>()
            );
        }
        begin_segm = end_segm - dummy_prefix_size;
    }
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

std::vector<std::string> KmerExtractor::generate_suffixes(size_t len) {
    std::vector<std::string> valid_suffixes;

    const char sentinel_char = alphabet[0];

    for (auto&& suffix : utils::generate_strings(alphabet, len)) {
        size_t last = suffix.rfind(sentinel_char);
        if (last == std::string::npos
                || suffix.substr(0, last + 1)
                    == std::string(last + 1, sentinel_char))
            valid_suffixes.push_back(std::move(suffix));
    }

    return valid_suffixes;
}


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
    static_assert(bits_per_char <= sizeof(TAlphabet) * 8,
                  "Choose type for TAlphabet properly");

    assert(alphabet.size() <= (1llu << bits_per_char));
}

KmerExtractor2BitTDecl(KmerExtractor::TAlphabet)
::encode(char s) const {
    return extractor::encode(s, char_to_code_);
}

KmerExtractor2BitTDecl(char)
::decode(TAlphabet c) const {
    return extractor::decode(c, alphabet);
}

KmerExtractor2BitTDecl(std::vector<KmerExtractor2Bit::TAlphabet>)
::encode(const std::string &sequence) const {
    return extractor::encode(sequence.begin(), sequence.end(), char_to_code_);
}

KmerExtractor2BitTDecl(std::string)
::decode(const std::vector<TAlphabet> &sequence) const {
    return extractor::decode(sequence, alphabet);
}

KmerExtractor2BitTDecl(std::string)
::reverse_complement(const std::string &sequence) const {
    std::string rev(sequence.size(), 0);
    std::transform(sequence.rbegin(), sequence.rend(), rev.begin(),
        [&](char c) {
            auto code = encode(c);
            return code >= alphabet.size() ? code : decode(complement_code_[code]);
        }
    );
    return rev;
}

KmerExtractor2BitTDecl(std::vector<std::string>)
::generate_suffixes(size_t len) const {
    std::vector<std::string> result;
    for (auto&& suffix : utils::generate_strings(alphabet, len)) {
        result.push_back(std::move(suffix));
    }
    return result;
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
    assert(suffix.size() <= k);

    if (sequence.size() < k)
        return;

    auto seq = encode(sequence);

    TAlphabet * __restrict begin_segm = seq.data();
    TAlphabet * __restrict end_segm = seq.data();
    TAlphabet * __restrict end = seq.data() + seq.size();

    do {
        assert(end >= end_segm);
        end_segm = std::find_if(end_segm, end,
            [&](auto c) { return c >= alphabet.size(); }
        );
        if (begin_segm + k <= end_segm) {
            extractor::sequence_to_kmers(
                begin_segm, end_segm, k, suffix, kmers,
                canonical_mode ? complement_code_
                               : std::vector<TAlphabet>()
            );
        }
        begin_segm = ++end_segm;
    } while (begin_segm + k <= end);
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

template class KmerExtractor2BitT<alphabets::kBitsPerCharDNA>;


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
