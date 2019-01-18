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
    std::transform(sequence.begin(), sequence.end(),
                   str.begin(),
        [&alphabet](TAlphabet x) {
            return decode(x, alphabet);
        }
    );
    return str;
}

template <typename TAlphabet>
inline TAlphabet encode(char s,
                        const TAlphabet *kCharToNucleotide) {
    return s >= 0
        ? kCharToNucleotide[static_cast<size_t>(s)]
        : kCharToNucleotide[0];
}

template <typename TAlphabet>
inline std::vector<TAlphabet> encode(const std::string &sequence,
                                     const TAlphabet *kCharToNucleotide) {
    std::vector<TAlphabet> seq_encoded(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(),
        [&](char c) {
            return encode(c, kCharToNucleotide);
        }
    );
    return seq_encoded;
}

template <class KMER, typename TAlphabet>
inline void sequence_to_kmers(const std::vector<TAlphabet> &seq,
                              size_t k,
                              const std::vector<TAlphabet> &suffix,
                              Vector<KMER> *kmers,
                              const std::vector<uint64_t> &canonical_map) {
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
        typename KMER::KMerWordType forward, rev;
        for (size_t i = 0; i + k <= seq.size(); ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k - suffix.size()])) {
                if (canonical_map.empty()) {
                    kmers->emplace_back(&seq[i], k);
                } else {
                    forward = KMER::pack_kmer(&seq[i], k);
                    rev = KMER::pack_kmer(&rev_comp[seq.size() - i - k], k);
                    if (rev >= forward) {
                        kmers->emplace_back(forward);
                    } else if (std::equal(suffix.begin(), suffix.end(),
                                          &rev_comp[seq.size() - i - suffix.size()])) {
                        kmers->emplace_back(rev);
                    }
                }
            } else if (canonical_map.size()
                && std::equal(suffix.begin(), suffix.end(),
                              &rev_comp[seq.size() - i - suffix.size()])) {
                forward = KMER::pack_kmer(&seq[i], k);
                rev = KMER::pack_kmer(&rev_comp[seq.size() - i - k], k);
                if (rev < forward)
                    kmers->emplace_back(rev);
            }
        }
    } else {
        // initialize and add the first kmer from sequence
        auto kmer = KMER::pack_kmer(seq, k);
        KMER rev;
        if (canonical_map.size())
            rev = KMER(KMER::pack_kmer(&rev_comp[seq.size() - k], k));

        if (std::equal(suffix.begin(), suffix.end(),
                       &seq[k - suffix.size()])) {
            if (canonical_map.empty() || rev.data() >= kmer) {
                kmers->emplace_back(kmer);
            } else if (std::equal(suffix.begin(), suffix.end(),
                                  &rev_comp[seq.size() - suffix.size()])) {
                kmers->emplace_back(rev);
            }
        } else if (canonical_map.size()
            && rev.data() < kmer
            && std::equal(suffix.begin(), suffix.end(), &rev_comp[seq.size() - suffix.size()])) {
            kmers->emplace_back(rev);
        }

        // add all other kmers
        for (size_t i = 1; i + k <= seq.size(); ++i) {
            KMER::update_kmer(k, seq[i + k - 1], seq[i + k - 2], &kmer);
            if (canonical_map.size())
                rev = rev.prev_kmer(k, rev_comp[seq.size() - i - k]);

            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k - suffix.size()])) {
                if (canonical_map.empty() || rev.data() >= kmer) {
                    kmers->emplace_back(kmer);
                } else if (std::equal(suffix.begin(), suffix.end(),
                                   &rev_comp[seq.size() - i - suffix.size()])) {
                    kmers->emplace_back(rev);
                }
            } else if (canonical_map.size()
                && rev.data() < kmer
                && std::equal(suffix.begin(), suffix.end(), &rev_comp[seq.size() - i - suffix.size()])) {
                kmers->emplace_back(rev);
            }
        }
    }
}

template <class KMER>
inline std::string kmer_to_sequence(const KMER &kmer,
                                    const std::string &alphabet,
                                    size_t k) {
    return kmer.to_string(k, alphabet);
}

} // namespace extractor


/*
 * KmerExtractor
 */
#if _PROTEIN_GRAPH
    const std::string KmerExtractor::alphabet = "$ABCDEFGHIJKLMNOPQRSTUVWYZX";
    const KmerExtractor::TAlphabet KmerExtractor::kCharToNucleotide[128] = {
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
        26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
        16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26,
        26,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15,
        16, 17, 18, 19,  20, 21, 22, 23,  26, 24, 25, 26,  26, 26, 26, 26
    };
    const std::vector<uint64_t> canonical_map;
#elif _DNA_CASE_SENSITIVE_GRAPH
    //for case-specific DNA and RNA (U <-> T) data
    const std::string KmerExtractor::alphabet = "$ACGTNacgt";
    const KmerExtractor::TAlphabet KmerExtractor::kCharToNucleotide[128] = {
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 6, 5, 7,  5, 5, 5, 8,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  9, 9, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
    };
    const std::vector<uint64_t> canonical_map = {0, 4, 3, 2, 1, 0, 8, 7, 6, 5};
#elif _DNA_GRAPH
    //for DNA and RNA (U <-> T) alphabets
    const std::string KmerExtractor::alphabet = "$ACGTN";
    const KmerExtractor::TAlphabet KmerExtractor::kCharToNucleotide[128] = {
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 4, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
    };
    const std::vector<uint64_t> canonical_map = {0, 4, 3, 2, 1, 0};
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
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
        [&](TAlphabet c) {
            return c < alphabet.size();
        }
    ));
    #endif
    return extractor::encode(sequence, kCharToNucleotide);
}

std::string KmerExtractor::decode(const std::vector<TAlphabet> &sequence) {
    return extractor::decode(sequence, alphabet);
}

/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 * suffix does not take into account the last character.
 */
template <class KMER>
void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers,
                                      const std::vector<uint64_t> &canonical_map) {
    extractor::sequence_to_kmers(seq, k, suffix, kmers, canonical_map);
}

template
void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet>&,
                                      size_t,
                                      const std::vector<TAlphabet>&,
                                      Vector<Kmer64>*,
                                      const std::vector<uint64_t>&);
template
void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet>&,
                                      size_t,
                                      const std::vector<TAlphabet>&,
                                      Vector<Kmer128>*,
                                      const std::vector<uint64_t>&);
template
void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet>&,
                                      size_t,
                                      const std::vector<TAlphabet>&,
                                      Vector<Kmer256>*,
                                      const std::vector<uint64_t>&);

template <class KMER>
std::string KmerExtractor::kmer_to_sequence(const KMER &kmer, size_t k) {
    return extractor::kmer_to_sequence(kmer, alphabet, k);
}

template std::string KmerExtractor::kmer_to_sequence(const Kmer64&, size_t);
template std::string KmerExtractor::kmer_to_sequence(const Kmer128&, size_t);
template std::string KmerExtractor::kmer_to_sequence(const Kmer256&, size_t);


/**
 * Break the sequence to kmers and extend the temporary kmers storage.
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

    std::transform(sequence.begin(), sequence.end(),
                   &seq[dummy_prefix_size],
                   [](char c) { return encode(c); });

    sequence_to_kmers(seq, k, suffix, kmers,
        canonical_mode ? canonical_map : std::vector<uint64_t>{});
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

/*
 * KmerExtractor2Bit
 */
const std::string KmerExtractor2Bit::alphabet = "ACGT";
const KmerExtractor2Bit::TAlphabet KmerExtractor2Bit::kCharToNucleotide[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 0, 1,  0, 4, 4, 2,  4, 4, 4, 0,  4, 0, 0, 4,
    4, 4, 0, 0,  3, 3, 0, 0,  4, 0, 4, 0,  0, 0, 0, 0,
    4, 0, 0, 1,  0, 4, 4, 2,  4, 4, 4, 0,  4, 0, 0, 4,
    4, 4, 0, 0,  3, 3, 0, 0,  4, 0, 4, 0,  0, 0, 0, 0
};
const std::vector<uint64_t> canonical_map_2bit = {3,2,1,0};

KmerExtractor2Bit::KmerExtractor2Bit() {
    assert(alphabet.size() <= (1llu << kLogSigma));
}

KmerExtractor::TAlphabet KmerExtractor2Bit::encode(char s) {
    return extractor::encode(s, kCharToNucleotide);
}

char KmerExtractor2Bit::decode(TAlphabet c) {
    return extractor::decode(c, alphabet);
}

std::vector<KmerExtractor2Bit::TAlphabet>
KmerExtractor2Bit::encode(const std::string &sequence) {
    auto seq_encoded = extractor::encode(sequence, kCharToNucleotide);
    return std::any_of(seq_encoded.begin(), seq_encoded.end(),
                       [](TAlphabet c) { return c >= alphabet.size(); })
        ? std::vector<TAlphabet>{}
        : seq_encoded;
}

std::string KmerExtractor2Bit::decode(const std::vector<TAlphabet> &sequence) {
    return extractor::decode(sequence, alphabet);
}


/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 * suffix does not take into account the last character.
 */
template <class KMER>
void KmerExtractor2Bit::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                          size_t k,
                                          const std::vector<TAlphabet> &suffix,
                                          Vector<KMER> *kmers,
                                          const std::vector<uint64_t> &canonical_map) {
    extractor::sequence_to_kmers(seq, k, suffix, kmers, canonical_map);
}

template
void KmerExtractor2Bit::sequence_to_kmers(const std::vector<TAlphabet>&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer64>*,
                                          const std::vector<uint64_t>&);
template
void KmerExtractor2Bit::sequence_to_kmers(const std::vector<TAlphabet>&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer128>*,
                                          const std::vector<uint64_t>&);
template
void KmerExtractor2Bit::sequence_to_kmers(const std::vector<TAlphabet>&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer256>*,
                                          const std::vector<uint64_t>&);


template <class KMER>
std::string KmerExtractor2Bit::kmer_to_sequence(const KMER &kmer, size_t k) {
    return extractor::kmer_to_sequence(kmer, alphabet, k);
}

template std::string KmerExtractor2Bit::kmer_to_sequence(const Kmer64&, size_t);
template std::string KmerExtractor2Bit::kmer_to_sequence(const Kmer128&, size_t);
template std::string KmerExtractor2Bit::kmer_to_sequence(const Kmer256&, size_t);


/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
template <typename KMER>
void KmerExtractor2Bit::sequence_to_kmers(const std::string &sequence,
                                          size_t k,
                                          const std::vector<TAlphabet> &suffix,
                                          Vector<KMER> *kmers,
                                          bool canonical_mode) {
    assert(k);
    // suffix does not include the last character
    assert(suffix.size() < k);

    if (sequence.size() < k)
        return;

    sequence_to_kmers(encode(sequence), k, suffix, kmers,
        canonical_mode ? canonical_map_2bit : std::vector<uint64_t>{});
}

template
void KmerExtractor2Bit::sequence_to_kmers(const std::string&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer64>*,
                                          bool);
template
void KmerExtractor2Bit::sequence_to_kmers(const std::string&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer128>*,
                                          bool);
template
void KmerExtractor2Bit::sequence_to_kmers(const std::string&,
                                          size_t,
                                          const std::vector<TAlphabet>&,
                                          Vector<Kmer256>*,
                                          bool);
