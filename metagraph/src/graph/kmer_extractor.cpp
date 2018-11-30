#include "kmer_extractor.hpp"

#include <algorithm>
#include <cstdlib>


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
    assert((s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                   : kCharToNucleotide[0]) < alphabet.size());

    return s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                  : kCharToNucleotide[0];
}

char KmerExtractor::decode(TAlphabet c) {
    assert(c < alphabet.size());
    return alphabet[c];
}

std::vector<KmerExtractor::TAlphabet>
KmerExtractor::encode(const std::string &sequence) {
    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), [](char c) { return encode(c); });
    return seq_encoded;
}

std::string KmerExtractor::decode(const std::vector<TAlphabet> &sequence) {
    std::string str(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(),
                   str.begin(), [](TAlphabet x) { return decode(x); });
    return str;
}

/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 * suffix does not take into account the last character.
 */
template <class KMER>
void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers) {
    assert(k);
    assert(suffix.size() < k);
    assert(kmers);

    if (seq.size() < k)
        return;

    // based on performance comparison
    // for KMer::pack_kmer and KMer::update_kmer
    if (suffix.size() > 1) {
        for (size_t i = 0; i < seq.size() - k + 1; ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k - 1] - suffix.size())) {
                kmers->emplace_back(&seq[i], k);
            }
        }
    } else {
        // initialize and add the first kmer from sequence
        auto kmer = KMER::pack_kmer(seq.data(), k);

        if (std::equal(suffix.begin(), suffix.end(),
                       &seq[k - 1] - suffix.size())) {
            kmers->emplace_back(kmer);
        }

        // add all other kmers
        for (size_t i = 1; i < seq.size() - k + 1; ++i) {
            KMER::update_kmer(k, seq[i + k - 1], seq[i + k - 2], &kmer);

            if (std::equal(suffix.begin(), suffix.end(),
                           &seq[i + k - 1] - suffix.size())) {
                kmers->emplace_back(kmer);
            }
        }
    }
}

template void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer64> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer128> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::vector<TAlphabet> &seq,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer256> *kmers);

template <class KMER>
std::string KmerExtractor::kmer_to_sequence(const KMER &kmer) {
    return kmer.to_string(alphabet);
}

template std::string KmerExtractor::kmer_to_sequence(const Kmer64 &kmer);
template std::string KmerExtractor::kmer_to_sequence(const Kmer128 &kmer);
template std::string KmerExtractor::kmer_to_sequence(const Kmer256 &kmer);

/**
 * Break the sequence to kmers and extend the temporary kmers storage.
 */
template <typename KMER>
void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                      size_t k,
                                      const std::vector<TAlphabet> &suffix,
                                      Vector<KMER> *kmers) {
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

    KmerExtractor::sequence_to_kmers(seq, k, suffix, kmers);
}

template void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer64> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer128> *kmers);
template void KmerExtractor::sequence_to_kmers(const std::string &sequence,
                                               size_t k,
                                               const std::vector<TAlphabet> &suffix,
                                               Vector<Kmer256> *kmers);


const std::string KmerExtractor2Bit::alphabet = "ACGT";
const KmerExtractor2Bit::TAlphabet KmerExtractor2Bit::kCharToNucleotide[128] = {
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

KmerExtractor2Bit::KmerExtractor2Bit() {
    assert(alphabet.size() <= (1llu << kLogSigma));
}

KmerExtractor2Bit::TAlphabet KmerExtractor2Bit::encode(char s) const {
    assert((s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                   : kCharToNucleotide[0]) < alphabet.size());

    return s >= 0 ? kCharToNucleotide[static_cast<size_t>(s)]
                  : kCharToNucleotide[0];
}

char KmerExtractor2Bit::decode(TAlphabet c) const {
    assert(c < alphabet.size());
    return alphabet[c];
}

std::vector<KmerExtractor2Bit::TAlphabet>
KmerExtractor2Bit::encode(const std::string &sequence) const {
    std::vector<TAlphabet> seq_encoded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   seq_encoded.begin(), [this](char c) { return encode(c); });
    return seq_encoded;
}

std::string KmerExtractor2Bit::decode(const std::vector<TAlphabet> &sequence) const {
    std::string str(sequence.size(), 0);
    std::transform(sequence.begin(), sequence.end(),
                   str.begin(), [this](TAlphabet x) { return decode(x); });
    return str;
}

Vector<KmerExtractor2Bit::Kmer>
KmerExtractor2Bit
::sequence_to_kmers(std::vector<TAlphabet>&& seq, size_t k) const {
    assert(k);

    if (seq.size() < k)
        return {};

    Vector<Kmer> kmers;
    kmers.reserve(seq.size() - k + 1);

    // initialize and add the first kmer from sequence
    auto kmer = Kmer::pack_kmer(seq.data(), k);
    kmers.emplace_back(kmer);

    // add all other kmers
    for (size_t i = 1; i < seq.size() - k + 1; ++i) {
        Kmer::update_kmer(k, seq[i + k - 1], seq[i + k - 2], &kmer);
        kmers.emplace_back(kmer);
    }
    return kmers;
}

Vector<KmerExtractor2Bit::Kmer>
KmerExtractor2Bit::sequence_to_kmers(const std::string &sequence, size_t k) const {
    assert(k);
    return sequence_to_kmers(encode(sequence), k);
}

std::string KmerExtractor2Bit::kmer_to_sequence(const Kmer &kmer) const {
    return kmer.to_string(alphabet);
}
