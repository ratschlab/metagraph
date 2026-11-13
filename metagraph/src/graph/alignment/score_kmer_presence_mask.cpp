#include "score_kmer_presence_mask.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>

#include "common/aligned_vector.hpp"
#include "common/vectors/vector_algorithm.hpp"


namespace mtg {
namespace graph {
namespace align {

namespace {

// The lengths of the blocks of zeros (table[0]) and of ones (table[1]) of a mask of
// |size| positions after its ones are kept only where the next |offset| - 1 positions
// within the mask are ones too. Every block but the last one counts |correction| more.
// |call_runs| gives the runs of ones of the mask (ascending, maximal). A run ending
// before the end of the mask loses its last |offset| - 1 positions (all of them if it
// is not longer), and a mask shorter than |offset| stays as it is.
std::array<AlignedVector<size_t>, 2>
tabulate_runs(size_t size, const CallPresentRuns &call_runs, uint8_t offset, size_t correction) {
    std::array<AlignedVector<size_t>, 2> table;
    if (!size)
        return table;

    const size_t shrink = offset - 1;
    bool last_block = false;
    size_t last_size = 0;
    bool started = false;
    // start the next block of zeros or ones, ending the previous one
    auto block = [&](bool bit, size_t length) {
        if (started)
            table[last_block].push_back(last_size + correction);
        started = true;
        last_block = bit;
        last_size = length;
    };
    size_t pos = 0;
    call_runs([&](size_t begin, size_t end) {
        assert(begin >= pos && begin < end && end <= size);
        size_t kept_end = end;
        if (size >= offset && end < size)
            kept_end = end - begin > shrink ? end - shrink : begin;
        if (kept_end == begin)
            return;
        if (begin > pos)
            block(false, begin - pos);
        block(true, kept_end - begin);
        pos = kept_end;
    });
    if (pos < size)
        block(false, size - pos);
    table[last_block].push_back(last_size);
    return table;
}

// the score of a mask of |size| k-mers from the table of its blocks
int32_t score_table(size_t k,
                    size_t size,
                    const std::array<AlignedVector<size_t>, 2> &score_counter,
                    int32_t kmer_adjust,
                    int32_t match_score,
                    int32_t mismatch_score) {
    const size_t sequence_length = size + k - 1;
    const int32_t SNP_t = k + kmer_adjust;

    double score = std::accumulate(score_counter[1].begin(),
                                   score_counter[1].end(), 0) * match_score;
    if (score == 0)
        return 0;

    if (score_counter[0].empty())
        return score * sequence_length / size;

    for (double count : score_counter[0]) {
        // A penalty function used in BIGSI
        double min_N_snps = count / SNP_t;
        double max_N_snps = std::max(count - SNP_t + 1, min_N_snps);
        double mean_N_snps = max_N_snps * 0.05 + min_N_snps;

        assert(count >= mean_N_snps);

        double mean_penalty = mean_N_snps * mismatch_score;
        score += (count - mean_penalty) * match_score - mean_penalty;
    }

    return std::max(score * sequence_length / size, 0.);
}

const int32_t kmer_adjust = 3;

} // namespace

int32_t score_kmer_presence_mask(size_t k,
                                 size_t size,
                                 const CallPresentRuns &call_present_runs,
                                 int32_t match_score,
                                 int32_t mismatch_score) {
    if (!size)
        return 0;

    return score_table(k, size, tabulate_runs(size, call_present_runs, kmer_adjust, 1),
                       kmer_adjust, match_score, mismatch_score);
}

int32_t score_kmer_presence_mask(size_t k,
                                 const sdsl::bit_vector &kmer_presence_mask,
                                 int32_t match_score,
                                 int32_t mismatch_score) {
    return score_kmer_presence_mask(k, kmer_presence_mask.size(), [&](const auto &callback) {
        call_runs_of_ones(kmer_presence_mask, callback);
    }, match_score, mismatch_score);
}

} // namespace align
} // namespace graph
} // namespace mtg
