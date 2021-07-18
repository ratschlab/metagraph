#include "kmer_extractor.hpp"

#include <algorithm>
#include <cstdlib>

#include "common/algorithms.hpp"


namespace {

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
inline std::vector<TAlphabet> encode(std::string_view sequence,
                                     const TAlphabet kCharToNucleotide[]) {
    std::vector<TAlphabet> seq_encoded(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(), seq_encoded.begin(),
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
                                const Callback &callback,
                                const Call &skip) {
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
                                      const Callback &callback,
                                      const Call &skip) {
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
                                          const Callback &callback,
                                          const Call &skip) {
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
                            const Callback &callback,
                            const Call &skip) {
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
                                                const Callback &callback,
                                                const Call &skip) {
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
                              const Callback &callback,
                              const std::vector<uint8_t> &complement_code,
                              const Call &skip) {
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

} // namespace


namespace mtg {
namespace kmer {

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
    return ::encode(s, kCharToNucleotide);
}

char KmerExtractorBOSS::decode(TAlphabet c) {
    return ::decode(c, alphabet);
}

std::vector<KmerExtractorBOSS::TAlphabet>
KmerExtractorBOSS::encode(std::string_view sequence) {
    return ::encode(sequence, kCharToNucleotide);
}

std::string KmerExtractorBOSS::decode(const std::vector<TAlphabet> &sequence) {
    return ::decode(sequence, alphabet);
}

void KmerExtractorBOSS::reverse_complement(std::vector<TAlphabet> *sequence) {
    assert(sequence);
    return ::reverse_complement(sequence->begin(), sequence->end(), kComplementCode);
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
void KmerExtractorBOSS::sequence_to_kmers(std::string_view sequence,
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
    const size_t dummy_prefix_size = suffix.empty() ? 0 : k - 1;
    const size_t dummy_suffix_size = suffix.empty() ? 0 : 1;

    std::vector<TAlphabet> seq(dummy_prefix_size
                                    + sequence.size() + 1, alphabet.size());

    std::transform(sequence.begin(), sequence.end(), &seq[dummy_prefix_size],
        [](char c) { return encode(c); }
    );

    TAlphabet *end_segm = seq.data() + dummy_prefix_size;
    TAlphabet *end = seq.data() + seq.size();

    do {
        assert(end >= end_segm);

        TAlphabet *begin_segm = end_segm - dummy_prefix_size;

        // DNA segments may be stored continuously separated by a character outside of
        // the valid alphabet (usually 'N'). Each segment is treated as a distinct DNA
        // read and will be prepended with a dummy prefix if #suffix is not empty
        while (*end_segm < alphabet.size()) {
            end_segm++;
        }
        assert(end_segm < end);

        if (begin_segm + dummy_prefix_size + k <= end_segm) {
            if (dummy_suffix_size) {
                assert(*end_segm >= alphabet.size());
                // set the dummy suffix
                // ***NAAANAAAA* ->
                // ***NAAA*AAAA*
                //     ---^
                *end_segm = 0;

                assert(dummy_prefix_size);
                // set the dummy prefix for the next segment
                // ***NAAANAAAA* ->
                // ****AAA*AAAA*
                //   ^^---^
                std::fill(begin_segm, begin_segm + dummy_prefix_size, 0);
            }

            ::sequence_to_kmers<KMER>(begin_segm, end_segm + dummy_suffix_size, k, suffix,
                [&kmers](auto kmer) { kmers->push_back(kmer); },
                canonical_mode ? kComplementCode : std::vector<uint8_t>(),
                []() { return false; }
            );
        }

    } while (++end_segm + k + dummy_suffix_size <= end);
}

template
void KmerExtractorBOSS::sequence_to_kmers(std::string_view,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer64>*,
                                          bool);
template
void KmerExtractorBOSS::sequence_to_kmers(std::string_view,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer128>*,
                                          bool);
template
void KmerExtractorBOSS::sequence_to_kmers(std::string_view,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer256>*,
                                          bool);

template <typename KMER>
KMER KmerExtractorBOSS::sequence_to_kmer(std::string_view sequence,
                                         bool canonical_mode) {
    // encode sequence
    std::vector<TAlphabet> seq(sequence.size(), alphabet.size());

    std::transform(sequence.begin(), sequence.end(), seq.begin(), [](char c) {
        return c != alphabet[0] ? encode(c) : 0;
    });

    if (std::find(seq.begin(), seq.end(), alphabet.size()) != seq.end())
        throw std::runtime_error("Invalid characters encoded");

    KMER ret_val{0};
#ifndef NDEBUG
    bool encoded = false;
#endif
    ::sequence_to_kmers<KMER>(seq.data(), seq.data() + seq.size(), seq.size(), {},
        [&](auto kmer) {
#ifndef NDEBUG
            encoded = true;
#endif
            ret_val = kmer;
        },
        canonical_mode ? kComplementCode : std::vector<uint8_t>(),
        []() { return false; }
    );

    assert(encoded);

    return ret_val;
}

template KmerExtractorBOSS::Kmer64 KmerExtractorBOSS::sequence_to_kmer<KmerExtractorBOSS::Kmer64>(std::string_view, bool);
template KmerExtractorBOSS::Kmer128 KmerExtractorBOSS::sequence_to_kmer<KmerExtractorBOSS::Kmer128>(std::string_view, bool);
template KmerExtractorBOSS::Kmer256 KmerExtractorBOSS::sequence_to_kmer<KmerExtractorBOSS::Kmer256>(std::string_view, bool);


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

#define KmerExtractorTDecl(...) \
template <const uint8_t LogSigma> \
__VA_ARGS__ KmerExtractorT<LogSigma>

KmerExtractorTDecl()
::KmerExtractorT(const char alph[],
                 const uint8_t char_to_code[128],
                 const std::vector<uint8_t> &complement_code)
      : alphabet(alph),
        char_to_code_(char_to_code),
        complement_code_(complement_code) {
    static_assert(bits_per_char <= sizeof(TAlphabet) * 8,
                  "Choose type for TAlphabet properly");

    assert(alphabet.size() <= (1llu << bits_per_char));
}

KmerExtractorTDecl(typename KmerExtractorT<LogSigma>::TAlphabet)
::encode(char s) const {
    return ::encode(s, char_to_code_);
}

KmerExtractorTDecl(char)
::decode(TAlphabet c) const {
    return ::decode(c, alphabet);
}

KmerExtractorTDecl(std::vector<typename KmerExtractorT<LogSigma>::TAlphabet>)
::encode(std::string_view sequence) const {
    return ::encode(sequence, char_to_code_);
}

KmerExtractorTDecl(std::string)
::decode(const std::vector<TAlphabet> &sequence) const {
    return ::decode(sequence, alphabet);
}

KmerExtractorTDecl(std::vector<std::string>)
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
KmerExtractorTDecl(template <typename KMER> void)
::sequence_to_kmers(std::string_view sequence,
                    size_t k,
                    const std::vector<TAlphabet> &suffix,
                    Vector<KMER> *kmers,
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

    ::sequence_to_kmers<KMER>(
        seq.data(), seq.data() + seq.size(), k, suffix,
        [&kmers](auto kmer) { kmers->push_back(kmer); },
        canonical_mode ? complement_code_ : std::vector<uint8_t>(),
        [&]() { return invalid[i++]; }
    );
}

KmerExtractorTDecl(template <typename KMER> Vector<std::pair<KMER, bool>>)
::sequence_to_kmers(std::string_view sequence,
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

    ::sequence_to_kmers<KMER>(
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

#if _PROTEIN_GRAPH
template class KmerExtractorT<alphabets::kBitsPerCharProtein>;
#elif _DNA_CASE_SENSITIVE_GRAPH
template class KmerExtractorT<alphabets::kBitsPerCharDNACaseSent>;
#elif _DNA5_GRAPH
template class KmerExtractorT<alphabets::kBitsPerCharDNA5>;
#elif _DNA_GRAPH
template class KmerExtractorT<alphabets::kBitsPerCharDNA>;
#else
static_assert(false,
    "Define an alphabet: either "
    "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
);
#endif


#define ExplicitInstantiation_sequence_to_kmers(KMER) \
template \
void KmerExtractor2Bit \
::sequence_to_kmers<KMER>(std::string_view, size_t, \
                          const std::vector<TAlphabet>&, Vector<KMER>*, bool) const; \
template \
Vector<std::pair<KMER, bool>> KmerExtractor2Bit \
::sequence_to_kmers<KMER>(std::string_view, size_t, \
                          bool, const std::vector<TAlphabet>&) const;

ExplicitInstantiation_sequence_to_kmers(KmerExtractor2Bit::Kmer64)
ExplicitInstantiation_sequence_to_kmers(KmerExtractor2Bit::Kmer128)
ExplicitInstantiation_sequence_to_kmers(KmerExtractor2Bit::Kmer256)

ExplicitInstantiation_sequence_to_kmers(KmerExtractor2Bit::KmerBOSS64)
ExplicitInstantiation_sequence_to_kmers(KmerExtractor2Bit::KmerBOSS128)
ExplicitInstantiation_sequence_to_kmers(KmerExtractor2Bit::KmerBOSS256)

} // namespace kmer
} // namespace mtg
