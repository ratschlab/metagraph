#include "kmer_extractor.hpp"

#include <algorithm>
#include <cstdlib>

#include "utils/algorithms.hpp"


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
inline std::vector<TAlphabet> encode(std::string::const_iterator begin,
                                     std::string::const_iterator end,
                                     const TAlphabet kCharToNucleotide[]) {
    assert(begin <= end);
    std::vector<TAlphabet> seq_encoded(end - begin, 0);
    std::transform(begin, end, seq_encoded.begin(),
        [&](char c) { return encode(c, kCharToNucleotide); }
    );
    return seq_encoded;
}

template <typename Iterator>
inline void reverse_complement(Iterator begin,
                               Iterator end,
                               const std::vector<uint8_t> &complement_code) {
    assert(end >= begin);
    assert(std::all_of(begin, end, [&](auto c) { return c < complement_code.size(); }));

    while (begin < --end) {
        auto temp_value = complement_code[static_cast<int>(*begin)];
        *begin = complement_code[static_cast<int>(*end)];
        *end = temp_value;
        ++begin;
    }

    if (begin == end)
        *begin = complement_code[static_cast<int>(*begin)];
}


/**
 * k-mer extractors
 */

template <class KMER, typename TAlphabet, typename Callback, typename Call>
inline void __sequence_to_kmers(const TAlphabet *begin,
                                const TAlphabet *end,
                                size_t k,
                                const std::vector<TAlphabet> &suffix,
                                Callback callback,
                                Call skip) {
    assert(end >= begin);
    // done to ensure that begin + k doesn't overflow
    assert(k <= static_cast<size_t>(end - begin));

    for (auto kmer_begin = begin; kmer_begin + k <= end; ++kmer_begin) {
        if (!skip() && KMER::match_suffix(kmer_begin, k, suffix))
            callback(KMER(kmer_begin, k));
    }
}

template <class KMER, typename TAlphabet, typename Callback, typename Call>
inline void __sequence_to_kmers_slide(const TAlphabet *begin,
                                      const TAlphabet *end,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Callback callback,
                                      Call skip) {
    assert(end >= begin);
    // done to ensure that begin + k doesn't overflow
    assert(k <= static_cast<size_t>(end - begin));

    // initialize and call the first kmer from sequence
    KMER kmer(begin, k);
    if (!skip() && KMER::match_suffix(begin, k, suffix))
        callback(kmer);

    // call all other kmers
    for (auto kmer_begin = begin + 1; kmer_begin + k <= end; ++kmer_begin) {
        kmer.to_next(k, kmer_begin[k - 1]);
        if (!skip() && KMER::match_suffix(kmer_begin, k, suffix))
            callback(kmer);
    }
}

// extract canonical (lexicographically smallest) k-mers
template <class KMER, typename TAlphabet, typename Callback, typename Call>
inline void __sequence_to_kmers_canonical(const TAlphabet *seq,
                                          const std::vector<TAlphabet> &rev_comp,
                                          size_t k,
                                          const std::vector<TAlphabet> &suffix,
                                          Callback callback,
                                          Call skip) {
    assert(rev_comp.size() >= k);

    for (size_t i = 0; i + k <= rev_comp.size(); ++i) {
        if (skip())
            continue;

        bool suffix_matched_forward = KMER::match_suffix(seq + i, k, suffix);
        bool suffix_matched_reverse = KMER::match_suffix(&rev_comp[rev_comp.size() - i - k], k, suffix);

        if (!suffix_matched_forward && !suffix_matched_reverse)
            continue;

        KMER forward(seq + i, k);
        KMER reverse(&rev_comp[rev_comp.size() - i - k], k);

        if (forward <= reverse) {
            if (suffix_matched_forward)
                callback(forward);
        } else {
            if (suffix_matched_reverse)
                callback(reverse);
        }
    }
}

template <class KMER, typename TAlphabet, typename Callback, typename Call>
inline void __call_smallest(const KMER &kmer_first,
                            const TAlphabet *first,
                            const KMER &kmer_second,
                            const TAlphabet *second,
                            size_t k,
                            const std::vector<TAlphabet> &suffix,
                            Callback callback,
                            Call skip) {
    if (skip())
        return;

    if (kmer_first <= kmer_second) {
        if (KMER::match_suffix(first, k, suffix))
            callback(kmer_first);
    } else {
        if (KMER::match_suffix(second, k, suffix))
            callback(kmer_second);
    }
}

// extract canonical (lexicographically smallest) k-mers
template <class KMER, typename TAlphabet, typename Callback, typename Call>
inline void __sequence_to_kmers_canonical_slide(const TAlphabet *seq,
                                                const std::vector<TAlphabet> &rev_comp,
                                                size_t k,
                                                const std::vector<TAlphabet> &suffix,
                                                Callback callback,
                                                Call skip) {
    assert(rev_comp.size() >= k);

    // initialize and call the first kmer from sequence
    KMER kmer(seq, k);
    KMER rev(&rev_comp[rev_comp.size() - k], k);

    __call_smallest<KMER>(kmer, seq,
                          rev, &rev_comp[rev_comp.size() - k],
                          k,
                          suffix, callback, skip);

    // call all other kmers
    for (size_t forward_last = k,
                reverse_first = rev_comp.size() - k - 1;
                                        forward_last < rev_comp.size();
                                            ++forward_last, --reverse_first) {
        kmer.to_next(k, seq[forward_last]);
        rev.to_prev(k, rev_comp[reverse_first]);

        __call_smallest<KMER>(kmer, seq + forward_last - (k - 1),
                              rev, &rev_comp[reverse_first],
                              k,
                              suffix, callback, skip);
    }
}


/**
 * Break the sequence into k-mers and call them.
 */
template <class KMER, typename TAlphabet, typename Callback, typename Call>
inline void sequence_to_kmers(const TAlphabet *begin,
                              const TAlphabet *end,
                              size_t k,
                              const std::vector<TAlphabet> &suffix,
                              Callback callback,
                              const std::vector<uint8_t> &complement_code,
                              Call skip) {
    assert(k);
    assert(suffix.size() <= k);
    assert(end >= begin);

    // done to avoid begin + k overflowing
    if (k > static_cast<size_t>(end - begin))
        return;

    if (complement_code.empty()) {
        // based on performance comparison
        // for KMer::pack_kmer and KMer::update_kmer
        if (suffix.size() > 1) {
            __sequence_to_kmers<KMER>(begin, end, k, suffix, callback, skip);
        } else {
            __sequence_to_kmers_slide<KMER>(begin, end, k, suffix, callback, skip);
        }
    } else {
        std::vector<TAlphabet> rev_comp(begin, end);
        reverse_complement(rev_comp.begin(), rev_comp.end(), complement_code);
        assert(rev_comp.size() == static_cast<size_t>(end - begin));

        if (suffix.size() > 1) {
            __sequence_to_kmers_canonical<KMER>(begin, rev_comp, k, suffix, callback, skip);
        } else {
            __sequence_to_kmers_canonical_slide<KMER>(begin, rev_comp, k, suffix, callback, skip);
        }
    }
}

} // namespace extractor


/**
 * KmerExtractorBOSS
 */

#if _PROTEIN_GRAPH
    const std::string KmerExtractorBOSS::alphabet = alphabets::kBOSSAlphabetProtein;
    const KmerExtractorBOSS::TAlphabet *KmerExtractorBOSS::kCharToNucleotide = alphabets::kBOSSCharToProtein;
    const std::vector<uint8_t> KmerExtractorBOSS::kComplementCode = alphabets::kBOSSComplementMapProtein;
#elif _DNA_CASE_SENSITIVE_GRAPH
    const std::string KmerExtractorBOSS::alphabet = alphabets::kBOSSAlphabetDNACaseSent;
    const KmerExtractorBOSS::TAlphabet *KmerExtractorBOSS::kCharToNucleotide = alphabets::kBOSSCharToDNACaseSent;
    const std::vector<uint8_t> KmerExtractorBOSS::kComplementCode = alphabets::kBOSSComplementMapDNACaseSent;
#elif _DNA5_GRAPH
    const std::string KmerExtractorBOSS::alphabet = alphabets::kBOSSAlphabetDNA5;
    const KmerExtractorBOSS::TAlphabet *KmerExtractorBOSS::kCharToNucleotide = alphabets::kBOSSCharToDNA;
    const std::vector<uint8_t> KmerExtractorBOSS::kComplementCode = alphabets::kBOSSComplementMapDNA;
#elif _DNA_GRAPH
    const std::string KmerExtractorBOSS::alphabet = alphabets::kBOSSAlphabetDNA;
    const KmerExtractorBOSS::TAlphabet *KmerExtractorBOSS::kCharToNucleotide = alphabets::kBOSSCharToDNA;
    const std::vector<uint8_t> KmerExtractorBOSS::kComplementCode = alphabets::kBOSSComplementMapDNA;
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

static_assert(KmerExtractorBOSS::bits_per_char <= sizeof(KmerExtractorBOSS::TAlphabet) * 8,
              "Choose type for TAlphabet properly");


KmerExtractorBOSS::KmerExtractorBOSS() {
    assert(alphabet.size() <= (1llu << bits_per_char));
}

KmerExtractorBOSS::TAlphabet KmerExtractorBOSS::encode(char s) {
    return extractor::encode(s, kCharToNucleotide);
}

char KmerExtractorBOSS::decode(TAlphabet c) {
    return extractor::decode(c, alphabet);
}

std::vector<KmerExtractorBOSS::TAlphabet>
KmerExtractorBOSS::encode(const std::string &sequence) {
    return extractor::encode(sequence.begin(), sequence.end(), kCharToNucleotide);
}

std::string KmerExtractorBOSS::decode(const std::vector<TAlphabet> &sequence) {
    return extractor::decode(sequence, alphabet);
}

void KmerExtractorBOSS::reverse_complement(std::vector<TAlphabet> *sequence) {
    assert(sequence);
    return extractor::reverse_complement(sequence->begin(), sequence->end(), kComplementCode);
}

KmerExtractorBOSS::TAlphabet KmerExtractorBOSS::complement(TAlphabet c) {
    assert(c < kComplementCode.size());
    return kComplementCode[c];
}

/**
 * Break the sequence into k-mers and add them to the k-mer storage.
 * @param sequence sequence to be broken into k-mers
 * @param k the k-mer length
 * @param suffix if not empty, only k-mers that match this suffix are kept
 * @param[out] kmers output parameter for the resulting k-mers
 * @param canonical_mode if true, extracts canonical (lexicographically smaller vs. the
 * reverse complement) k-mers
 */
template <typename KMER>
void KmerExtractorBOSS::sequence_to_kmers(const std::string &sequence,
                                          size_t k,
                                          const std::vector<TAlphabet> &suffix,
                                          Vector<KMER> *kmers,
                                          bool canonical_mode) {
    assert(kmers);
    assert(k);
    assert(suffix.size() < k && "suffix does not include the last character");

    if (sequence.size() < k)
        return;

    // encode sequence
    const size_t dummy_prefix_size = suffix.size() > 0 ? k - 1 : 1;

    std::vector<TAlphabet> seq(sequence.size() + dummy_prefix_size + 1, 0);

    std::transform(sequence.begin(), sequence.end(), &seq[dummy_prefix_size],
        [](char c) { return encode(c); }
    );

    TAlphabet *begin_segm = seq.data();
    TAlphabet *end_segm = seq.data() + dummy_prefix_size;
    TAlphabet *end = seq.data() + seq.size();

    while (begin_segm + dummy_prefix_size + k + 1 <= end) {
        assert(end >= end_segm);
        // DNA segments may be stored continuously separated by a character outside of
        // the valid alphabet (usually 'N'). Each segment is treated as a distinct DNA
        // read and must be prepended with a dummy prefix
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
            extractor::sequence_to_kmers<KMER>(begin_segm, end_segm, k, suffix,
                [&kmers](auto kmer) { kmers->push_back(kmer); },
                canonical_mode ? kComplementCode : std::vector<uint8_t>(),
                []() { return false; }
            );
        }
        begin_segm = end_segm - dummy_prefix_size;
    }
}

template
void KmerExtractorBOSS::sequence_to_kmers(const std::string&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer64>*,
                                          bool);
template
void KmerExtractorBOSS::sequence_to_kmers(const std::string&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer128>*,
                                          bool);
template
void KmerExtractorBOSS::sequence_to_kmers(const std::string&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer256>*,
                                          bool);

std::vector<std::string> KmerExtractorBOSS::generate_suffixes(size_t len) {
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

KmerExtractor2BitTDecl(typename KmerExtractor2BitT<LogSigma>::TAlphabet)
::encode(char s) const {
    return extractor::encode(s, char_to_code_);
}

KmerExtractor2BitTDecl(char)
::decode(TAlphabet c) const {
    return extractor::decode(c, alphabet);
}

KmerExtractor2BitTDecl(std::vector<typename KmerExtractor2BitT<LogSigma>::TAlphabet>)
::encode(const std::string &sequence) const {
    return extractor::encode(sequence.begin(), sequence.end(), char_to_code_);
}

KmerExtractor2BitTDecl(std::string)
::decode(const std::vector<TAlphabet> &sequence) const {
    return extractor::decode(sequence, alphabet);
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
    assert(suffix.size() <= k);

    if (sequence.size() < k)
        return;

    auto seq = encode(sequence);

    assert(std::all_of(seq.begin(), seq.end(),
                       [&](auto c) { return c <= alphabet.size(); }));

    // Mark where (k+1)-mers with invalid characters end
    // Example for (k+1)=3: [X]***[X]****[X]***
    //              ---->   [111]0[111]00[111]0
    auto invalid = utils::drag_and_mark_segments(seq, alphabet.size(), k);
    // Set invalid characters to zero so that k-mers don't overflow.
    // k-mers containing these invalid characters are invalid and will
    // be skipped anyway.
    std::replace(seq.begin(), seq.end(), static_cast<TAlphabet>(alphabet.size()),
                                         static_cast<TAlphabet>(0));
    size_t i = k - 1;

    extractor::sequence_to_kmers<Kmer<T>>(
        seq.data(), seq.data() + seq.size(), k, suffix,
        [&kmers](auto kmer) { kmers->push_back(kmer); },
        canonical_mode ? complement_code_ : std::vector<uint8_t>(),
        [&]() { return invalid[i++]; }
    );
}

KmerExtractor2BitTDecl(template <typename KMER> Vector<std::pair<KMER, bool>>)
::sequence_to_kmers(const std::string &sequence,
                    size_t k,
                    bool canonical_mode,
                    const std::vector<TAlphabet> &suffix) const {
    assert(k);
    assert(suffix.size() <= k);

    if (sequence.size() < k)
        return {};

    Vector<std::pair<KMER, bool>> kmers;
    kmers.reserve(sequence.length() + 1 - k);

    auto seq = encode(sequence);

    assert(std::all_of(seq.begin(), seq.end(),
                       [&](auto c) { return c <= alphabet.size(); }));

    // Mark where (k+1)-mers with invalid characters end
    // Example for (k+1)=3: [X]***[X]****[X]***
    //              ---->   [111]0[111]00[111]0
    auto invalid = utils::drag_and_mark_segments(seq, alphabet.size(), k);
    // Set invalid characters to zero so that k-mers don't overflow.
    // k-mers containing these invalid characters are invalid and will
    // be replaced anyway.
    std::replace(seq.begin(), seq.end(), static_cast<TAlphabet>(alphabet.size()),
                                         static_cast<TAlphabet>(0));
    size_t i = k - 1;

    extractor::sequence_to_kmers<KMER>(
        seq.data(), seq.data() + seq.size(), k, suffix,
        [&kmers](auto kmer) { kmers.emplace_back(kmer, true); },
        canonical_mode ? complement_code_ : std::vector<uint8_t>(),
        [&]() {
            if (!invalid[i++])
                return false;

            kmers.emplace_back(KMER(), false);
            return true;
        }
    );

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
Vector<std::pair<KmerExtractor2Bit::Kmer<T>, bool>> KmerExtractor2Bit \
::sequence_to_kmers<KmerExtractor2Bit::Kmer<T>>(const std::string&, \
                                                size_t, \
                                                bool, \
                                                const std::vector<TAlphabet>&) const;

ExplicitInstantiation_sequence_to_kmers_vector(uint64_t)
ExplicitInstantiation_sequence_to_kmers_vector(sdsl::uint128_t)
ExplicitInstantiation_sequence_to_kmers_vector(sdsl::uint256_t)
