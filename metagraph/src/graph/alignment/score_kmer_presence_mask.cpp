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

std::array<AlignedVector<size_t>, 2>
tabulate_score(const sdsl::bit_vector &presence, size_t correction = 0) {
    std::array<AlignedVector<size_t>, 2> table;
    table[0].reserve(presence.size());
    table[1].reserve(presence.size());

    if (!presence.size())
        return table;

    bool last_block = presence[0];
    size_t last_size = 1;
    size_t i = 1;

    if (presence.size() >= 64) {
        const auto word = *presence.data();
        if (!word || word == 0xFFFFFFFFFFFFFFFF) {
            last_size = 64;
            i = 64;
        }
    }

    for ( ; i < presence.size(); ++i) {
        if (!(i & 0x3F) && i + 64 <= presence.size()) {
            // if at a word boundary and the next word is either all zeros or
            // all ones
            const uint64_t word = presence.get_int(i);
            if (!word || word == 0xFFFFFFFFFFFFFFFF) {
                if ((!word && last_block) || (word == 0xFFFFFFFFFFFFFFFF && !last_block)) {
                    table[last_block].push_back(last_size + correction);
                    last_block = !last_block;
                    last_size = 0;
                }

                last_size += 64;
                i += 63;
                continue;
            }
        }

        if (last_block == presence[i]) {
            ++last_size;
        } else {
            table[last_block].push_back(last_size + correction);
            last_block = !last_block;
            last_size = 1;
        }
    }

    table[last_block].push_back(last_size);

    assert(std::accumulate(table[0].begin(), table[0].end(), size_t(0))
         + std::accumulate(table[1].begin(), table[1].end(), size_t(0))
         - correction * (table[0].size() + table[1].size() - 1)
        == presence.size());

    assert(std::accumulate(table[1].begin(), table[1].end(), size_t(0))
         - correction * (table[1].size() - last_block)
        == sdsl::util::cnt_one_bits(presence));

#ifndef NDEBUG
    if (correction == 1) {
        std::array<AlignedVector<size_t>, 2> check;

        size_t cnt = 1;
        for (size_t i = 0; i < presence.size(); ++i) {
            if (i < presence.size() - 1) {
                ++cnt;
                if (presence[i] != presence[i + 1]) {
                    check[presence[i]].push_back(cnt);
                    cnt = 1;
                }
            } else {
                check[presence[i]].push_back(cnt);
            }
        }

        assert(check[0] == table[0]);
        assert(check[1] == table[1]);
    }
#endif // NDEBUG

    return table;
}

// tabulate_score(autocorrelate(mask, offset), correction) of the mask of |size| positions whose
// ones are the runs [begin, end) |call_runs| gives (ascending, none empty and no two adjacent),
// built from the runs: autocorrelate keeps a one if the next |offset| - 1 positions within
// the mask are ones too, so a run that ends before the mask's end loses its last |offset| - 1
// positions (and vanishes if it has no more), and keeps a mask shorter than |offset| as it is;
// tabulate_score lists the lengths of the alternating runs of zeros and ones, each but the last
// plus |correction|.
std::array<AlignedVector<size_t>, 2>
tabulate_runs(size_t size,
              const CallPresentRuns &call_runs,
              uint8_t offset,
              size_t correction) {
    std::array<AlignedVector<size_t>, 2> table;
    if (!size)
        return table;

    const size_t shrink = offset - 1;
    bool last_block = false;
    size_t last_size = 0;
    bool started = false;
    // the next block, which ends the one before it
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

// The score of a mask of |size| k-mers from the table tabulate_score gives of its
// autocorrelation (by |kmer_adjust|, correction 1)
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
                                 const sdsl::bit_vector &kmer_presence_mask,
                                 int32_t match_score,
                                 int32_t mismatch_score) {
    if (!kmer_presence_mask.size())
        return 0;

    return score_table(k, kmer_presence_mask.size(),
                       tabulate_score(autocorrelate(kmer_presence_mask, kmer_adjust), 1),
                       kmer_adjust, match_score, mismatch_score);
}

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

} // namespace align
} // namespace graph
} // namespace mtg
