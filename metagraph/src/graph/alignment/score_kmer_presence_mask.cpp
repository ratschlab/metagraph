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

} // namespace

int32_t score_kmer_presence_mask(size_t k,
                                 const sdsl::bit_vector &kmer_presence_mask,
                                 int32_t match_score,
                                 int32_t mismatch_score) {
    if (!kmer_presence_mask.size())
        return 0;

    const int32_t kmer_adjust = 3;

    const size_t sequence_length = kmer_presence_mask.size() + k - 1;
    const int32_t SNP_t = k + kmer_adjust;

    auto score_counter = tabulate_score(autocorrelate(kmer_presence_mask, kmer_adjust), 1);

    double score = std::accumulate(score_counter[1].begin(),
                                   score_counter[1].end(), 0) * match_score;
    if (score == 0)
        return 0;

    if (score_counter[0].empty())
        return score * sequence_length / kmer_presence_mask.size();

    for (double count : score_counter[0]) {
        // A penalty function used in BIGSI
        double min_N_snps = count / SNP_t;
        double max_N_snps = std::max(count - SNP_t + 1, min_N_snps);
        double mean_N_snps = max_N_snps * 0.05 + min_N_snps;

        assert(count >= mean_N_snps);

        double mean_penalty = mean_N_snps * mismatch_score;
        score += (count - mean_penalty) * match_score - mean_penalty;
    }

    return std::max(score * sequence_length / kmer_presence_mask.size(), 0.);
}

} // namespace align
} // namespace graph
} // namespace mtg
