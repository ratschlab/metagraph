#include "aligner_chainer.hpp"

#include <unordered_set>

#include "dbg_aligner.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"

#include "common/aligned_vector.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/algorithms.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/rc_dbg.hpp"
#include "graph/graph_extensions/node_rc.hpp"
#include "graph/graph_extensions/hll_wrapper.hpp"

namespace mtg {
namespace graph {
namespace align {

using common::logger;

typedef DeBruijnGraph::node_index node_index;

constexpr uint32_t nid = std::numeric_limits<uint32_t>::max();

struct TableElem {
    Alignment::Column label;
    long long coordinate;
    int seed_clipping;
    int seed_end;
    score_t chain_score;
    uint32_t current_seed_index;

    TableElem(Alignment::Column c, long long coordinate, int seed_clipping,
              int seed_end, score_t chain_score, uint32 current_seed_index)
          : label(c), coordinate(coordinate), seed_clipping(seed_clipping),
            seed_end(seed_end), chain_score(chain_score), current_seed_index(current_seed_index) {}
} __attribute__((aligned(32)));
static_assert(sizeof(TableElem) == 32);

inline constexpr bool operator>(const TableElem &a, const TableElem &b) {
    return std::tie(a.label, a.coordinate, a.seed_clipping, a.seed_end)
        > std::tie(b.label, b.coordinate, b.seed_clipping, b.seed_end);
}

typedef AlignedVector<TableElem> ChainDPTable;

std::tuple<ChainDPTable, AlignedVector<int32_t>, size_t, size_t>
chain_seeds(const IDBGAligner &aligner,
            const DBGAlignerConfig &config,
            std::string_view query,
            std::vector<Seed> &seeds);

struct ChainHash {
    inline std::size_t operator()(const Chain &chain) const {
        uint64_t hash = 0;
        for (const auto &[aln, dist] : chain) {
            for (node_index node : aln.get_nodes()) {
                hash ^= node + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            hash ^= dist + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

std::pair<size_t, size_t>
call_seed_chains_both_strands(const IDBGAligner &aligner,
                              std::string_view forward,
                              std::string_view reverse,
                              const DBGAlignerConfig &config,
                              std::vector<Seed>&& fwd_seeds,
                              std::vector<Seed>&& bwd_seeds,
                              const std::function<void(Chain&&, score_t)> &callback,
                              const std::function<bool(Alignment::Column)> &skip_column) {
    fwd_seeds.erase(std::remove_if(fwd_seeds.begin(), fwd_seeds.end(),
                                   [](const auto &a) { return a.empty() || !a.label_columns; }),
                    fwd_seeds.end());
    bwd_seeds.erase(std::remove_if(bwd_seeds.begin(), bwd_seeds.end(),
                                   [](const auto &a) { return a.empty() || !a.label_columns; }),
                    bwd_seeds.end());

    if (fwd_seeds.empty() && bwd_seeds.empty())
        return { 0, 0 };


    bool has_labels = dynamic_cast<const ILabeledAligner*>(&aligner);

    // filter out empty seeds
    std::vector<Seed> both_seeds[2];
    both_seeds[0].reserve(fwd_seeds.size());
    both_seeds[1].reserve(bwd_seeds.size());
    for (auto&& fwd_seed : fwd_seeds) {
        if (!fwd_seed.empty())
            both_seeds[0].emplace_back(std::move(fwd_seed));
    }
    for (auto&& bwd_seed : bwd_seeds) {
        if (!bwd_seed.empty())
            both_seeds[1].emplace_back(std::move(bwd_seed));
    }

    fwd_seeds = std::vector<Seed>();
    bwd_seeds = std::vector<Seed>();

    // perform chaining on the forward, and the reverse-complement seeds
    ChainDPTable dp_tables[2];
    AlignedVector<int32_t> seed_backtraces[2];

    logger->trace("Chaining forward seeds");
    size_t num_seeds;
    size_t num_nodes;
    std::tie(dp_tables[0], seed_backtraces[0], num_seeds, num_nodes)
        = chain_seeds(aligner, config, forward, both_seeds[0]);

    logger->trace("Chaining reverse complement seeds");
    size_t num_seeds_bwd;
    size_t num_nodes_bwd;
    std::tie(dp_tables[1], seed_backtraces[1], num_seeds_bwd, num_nodes_bwd)
        = chain_seeds(aligner, config, reverse, both_seeds[1]);

    num_seeds += num_seeds_bwd;
    num_nodes += num_nodes_bwd;

    // construct chains by backtracking
    std::vector<std::tuple<score_t, uint32_t, ssize_t>> starts;
    starts.reserve(dp_tables[0].size() + dp_tables[1].size());
    sdsl::bit_vector both_used[2] {
        sdsl::bit_vector(dp_tables[0].size(), false),
        sdsl::bit_vector(dp_tables[1].size(), false)
    };

    for (size_t i = 0; i < dp_tables[0].size(); ++i) {
        starts.emplace_back(dp_tables[0][i].chain_score, 0, -static_cast<ssize_t>(i));
    }
    for (size_t i = 0; i < dp_tables[1].size(); ++i) {
        starts.emplace_back(dp_tables[1][i].chain_score, 1, -static_cast<ssize_t>(i));
    }

    if (starts.empty()) {
        logger->trace("No chains found");
        return std::make_pair(num_seeds, num_nodes);
    }

    std::sort(starts.begin(), starts.end(), std::greater<decltype(starts)::value_type>());

    score_t last_chain_score = std::numeric_limits<score_t>::min();
    std::unordered_multiset<Chain, ChainHash> chains;

    auto flush_chains = [&]() {
        assert(chains.size());
        auto it = chains.begin();
        Chain last_chain = *it;
        for (++it; it != chains.end(); ++it) {
            const Chain &chain = *it;
            if (chain != last_chain) {
                callback(std::move(last_chain), last_chain_score);
                last_chain = *it;
                continue;
            }

            // if this chain has the same seeds as the last one, merge their coordinate sets
            for (size_t i = 0; i < chain.size(); ++i) {
                Vector<Alignment::Column> columns;
                const Vector<Alignment::Column> &last_columns = last_chain[i].first.get_columns();
                const Vector<Alignment::Column> &cur_columns = chain[i].first.get_columns();
                if (chain[i].first.label_coordinates.size()) {
                    Alignment::CoordinateSet coord_union;
                    auto add_col_coords = [&](auto col, auto &coords) {
                        columns.push_back(col);
                        coord_union.emplace_back(std::move(coords));
                    };
                    utils::match_indexed_values(
                        last_columns.begin(), last_columns.end(),
                        last_chain[i].first.label_coordinates.begin(),
                        cur_columns.begin(), cur_columns.end(),
                        chain[i].first.label_coordinates.begin(),
                        [&](auto col, const auto &coords, const auto &other_coords) {
                            columns.push_back(col);
                            coord_union.emplace_back();
                            std::set_union(coords.begin(), coords.end(),
                                           other_coords.begin(), other_coords.end(),
                                           std::back_inserter(coord_union.back()));
                        },
                        add_col_coords, add_col_coords
                    );
                    std::swap(last_chain[i].first.label_coordinates, coord_union);
                } else {
                    std::set_union(last_columns.begin(), last_columns.end(),
                                   cur_columns.begin(), cur_columns.end(),
                                   std::back_inserter(columns));
                }
                last_chain[i].first.set_columns(std::move(columns));
            }
        }

        callback(std::move(last_chain), last_chain_score);

        chains.clear();
    };

    for (const auto &[chain_score, j, neg_i] : starts) {
        auto &used = both_used[j];
        uint32_t i = -neg_i;
        if (used[i])
            continue;

        const auto &dp_table = dp_tables[j];
        const auto &seeds = both_seeds[j];
        const auto &seed_backtrace = seed_backtraces[j];

        // iterate through the DP table, adding seeds to the chain
        std::vector<std::pair<Seed, int64_t>> chain_seeds;

        while (i != nid) {
            const auto &[label, coord, clipping, end, score, seed_i] = dp_table[i];
            if (skip_column(label))
                break;

            used[i] = true;
            chain_seeds.emplace_back(seeds[seed_i], coord);
            if (has_labels) {
                chain_seeds.back().first.set_columns(Vector<Alignment::Column>{ label });
                if (aligner.has_coordinates()) {
                    chain_seeds.back().first.label_coordinates.resize(1);
                    chain_seeds.back().first.label_coordinates[0].assign(1, coord);
                }
            } else {
                chain_seeds.back().first.label_encoder = nullptr;
            }
            i = seed_backtrace[i];
        }

        if (chain_seeds.empty())
            continue;

        // clean chain by merging overlapping seeds
        if (aligner.has_coordinates()) {
            for (size_t i = chain_seeds.size() - 1; i > 0; --i) {
                auto &cur_seed = chain_seeds[i].first;
                auto &prev_seed = chain_seeds[i - 1].first;

                assert(cur_seed.size());
                assert(prev_seed.size());
                assert(prev_seed.get_clipping() < cur_seed.get_clipping());
                assert(prev_seed.get_end_clipping() > cur_seed.get_end_clipping());

                size_t prev_end = prev_seed.get_clipping()
                                    + prev_seed.get_query_view().size();
                if (prev_end > cur_seed.get_clipping()) {
                    // they overlap
                    size_t coord_dist = cur_seed.label_coordinates[0][0]
                                            + cur_seed.get_query_view().size()
                                            - prev_seed.label_coordinates[0][0]
                                            - prev_seed.get_query_view().size();
                    size_t dist = cur_seed.get_clipping()
                                    + cur_seed.get_query_view().size() - prev_end;

                    if (dist == coord_dist && cur_seed.get_nodes().size() >= dist) {
                        prev_seed.expand({ cur_seed.get_nodes().end() - dist,
                                           cur_seed.get_nodes().end() });
                        cur_seed = Seed();
                    }
                }
            }
        }

        chain_seeds.erase(std::remove_if(chain_seeds.begin(), chain_seeds.end(),
                                         [](const auto &a) { return a.first.empty(); }),
                          chain_seeds.end());
        assert(chain_seeds.size());

        for (size_t i = chain_seeds.size() - 1; i > 0; --i) {
            assert(chain_seeds[i].first.get_clipping()
                        > chain_seeds[i - 1].first.get_clipping());
            assert(chain_seeds[i].first.get_end_clipping()
                        < chain_seeds[i - 1].first.get_end_clipping());
            chain_seeds[i].second -= chain_seeds[i - 1].second;
            assert(chain_seeds[i].second > 0);
        }

        chain_seeds[0].second = 0;
        if (!chain_seeds[0].first.label_columns)
            continue;

        Chain chain;
        chain.reserve(chain_seeds.size());
        std::transform(chain_seeds.begin(), chain_seeds.end(), std::back_inserter(chain),
                       [&](const auto &c) {
                           return std::make_pair(Alignment(c.first, config), c.second);
                       });

        if (chains.empty()) {
            chains.emplace(std::move(chain));
            last_chain_score = chain_score;
            continue;
        }

        if (chain_score == last_chain_score) {
            chains.emplace(std::move(chain));
            continue;
        }

        flush_chains();
        chains.emplace(std::move(chain));
        last_chain_score = chain_score;
    }

    flush_chains();

    return std::make_pair(num_seeds, num_nodes);
}

std::tuple<ChainDPTable, AlignedVector<int32_t>, size_t, size_t>
chain_seeds(const IDBGAligner &aligner,
            const DBGAlignerConfig &config,
            std::string_view query,
            std::vector<Seed> &seeds) {
    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    if (seeds.empty() || !labeled_aligner || !aligner.has_coordinates())
        return {};

    size_t num_nodes = 0;

    ssize_t query_size = query.size();

    ChainDPTable dp_table;
    dp_table.reserve(seeds.size());
    std::reverse(seeds.begin(), seeds.end());

    tsl::hopscotch_map<Alignment::Column, size_t> label_sizes;

    for (size_t i = 0; i < seeds.size(); ++i) {
        const Vector<Alignment::Column> &columns = seeds[i].get_columns();
        for (size_t j = 0; j < seeds[i].label_coordinates.size(); ++j) {
            Alignment::Column c = columns[j];
            auto rbegin = seeds[i].label_coordinates[j].rbegin();
            auto rend = rbegin + std::min(seeds[i].label_coordinates[j].size(),
                                          config.max_num_seeds_per_locus);
            std::for_each(rbegin, rend, [&](ssize_t coord) {
                ++label_sizes[c];
                dp_table.emplace_back(c, coord, seeds[i].get_clipping(),
                                      seeds[i].get_clipping() + seeds[i].get_query_view().size(),
                                      seeds[i].get_query_view().size(), i);
            });
        }
        seeds[i].label_columns = 0;
        seeds[i].label_coordinates = Alignment::CoordinateSet{};
    }

    dp_table.reserve(dp_table.size() + 9);

    size_t num_seeds = dp_table.size();
    AlignedVector<int32_t> backtrace(dp_table.size(), nid);
    if (dp_table.empty())
        return std::make_tuple(std::move(dp_table), std::move(backtrace), num_seeds, num_nodes);

    logger->trace("Sorting {} anchors", dp_table.size());
    // sort seeds by label, then by decreasing reference coordinate
    std::sort(dp_table.begin(), dp_table.end(), std::greater<TableElem>());
    logger->trace("Chaining anchors");

    size_t bandwidth = 65;

    // scoring function derived from minimap2
    // https://academic.oup.com/bioinformatics/article/34/18/3094/4994778
    float sl = static_cast<float>(config.min_seed_length) * 0.01;

    size_t cur_label_end = 0;
    size_t i = 0;
    while (cur_label_end < dp_table.size()) {
        cur_label_end += label_sizes[dp_table[i].label];
        for ( ; i < cur_label_end; ++i) {
            const auto &[prev_label, prev_coord, prev_clipping, prev_end,
                         prev_score, prev_seed_i] = dp_table[i];

            if (!prev_clipping)
                continue;

            size_t it_end = std::min(bandwidth, cur_label_end - i) + i;
            ssize_t coord_cutoff = prev_coord - query_size;
#if __AVX2__
            const __m256i coord_cutoff_v = _mm256_set1_epi64x(coord_cutoff);
            const __m256i prev_coord_v = _mm256_set1_epi64x(prev_coord);
            const __m256i prev_clipping_v = _mm256_set1_epi32(prev_clipping);
            const __m256i query_size_v = _mm256_set1_epi32(query_size);
            const __m256i prev_score_v = _mm256_set1_epi32(prev_score);
            const __m256i it_end_v = _mm256_set1_epi32(it_end - 1);
            const __m256i i_v = _mm256_set1_epi32(i);
            auto epi64_to_epi32 = [](__m256i v) {
                return _mm256_castsi256_si128(_mm256_permute4x64_epi64(_mm256_shuffle_epi32(v, 8), 8));
            };

            __m256i j_v = _mm256_add_epi32(_mm256_set1_epi32(i + 1), _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0));
            for (size_t j = i + 1; true; j += 8) {
                // if (coord_cutoff > coord || j > it_end - 1)
                //     break;
                __m256i coord_1_v = _mm256_i32gather_epi64(&dp_table[j].coordinate, _mm_set_epi32(12, 8, 4, 0), 8);
                __m256i coord_2_v = _mm256_i32gather_epi64(&dp_table[j + 4].coordinate, _mm_set_epi32(12, 8, 4, 0), 8);
                __m256i coord_1_mask = _mm256_cmpgt_epi64(coord_cutoff_v, coord_1_v);
                __m256i coord_2_mask = _mm256_cmpgt_epi64(coord_cutoff_v, coord_2_v);
                __m256i coord_neg_mask = _mm256_blend_epi32(coord_1_mask, coord_2_mask, 0b10101010);
                __m256i j_neg_mask = _mm256_cmpgt_epi32(j_v, it_end_v);
                coord_neg_mask = _mm256_or_si256(j_neg_mask, coord_neg_mask);

                // int32_t dist = prev_clipping - clipping;
                __m256i clipping_v = _mm256_i32gather_epi32(&dp_table[j].seed_clipping, _mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0), 4);
                __m256i dist_v = _mm256_sub_epi32(prev_clipping_v, clipping_v);

                // int32_t coord_dist = prev_coord - coord;
                // a[0:32],b[0:32],a[64:96],b[64:96],a[128:160],b[128:160],a[192:224],b[192:224]
                __m128i coord_dist_1_v = epi64_to_epi32(_mm256_sub_epi64(prev_coord_v, coord_1_v));
                __m128i coord_dist_2_v = epi64_to_epi32(_mm256_sub_epi64(prev_coord_v, coord_2_v));
                __m256i coord_dist_v = _mm256_set_m128i(coord_dist_2_v, coord_dist_1_v);

                __m256i dist_mask = _mm256_cmpgt_epi32(dist_v, _mm256_setzero_si256());
                __m256i dmax = _mm256_max_epi32(dist_v, coord_dist_v);
                __m256i dmax_mask = _mm256_cmpgt_epi32(query_size_v, dmax);
                dist_mask = _mm256_and_si256(dist_mask, dmax_mask);
                dist_mask = _mm256_andnot_si256(coord_neg_mask, dist_mask);

                // if (dist > 0 && std::max(dist, coord_dist) < query_size) {
                if (_mm256_movemask_epi8(dist_mask)) {
                    // score_t match = std::min({ dist, coord_dist, end - clipping });
                    __m256i match_v = _mm256_min_epi32(dist_v, coord_dist_v);
                    __m256i end_v = _mm256_i32gather_epi32(&dp_table[j].seed_end, _mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0), 4);
                    __m256i length_v = _mm256_sub_epi32(end_v, clipping_v);
                    match_v = _mm256_min_epi32(match_v, length_v);

                    // score_t cur_score = prev_score + match;
                    __m256i cur_score_v = _mm256_add_epi32(prev_score_v, match_v);

                    // float coord_diff = std::abs(coord_dist - dist);
                    __m256i coord_diff = _mm256_sub_epi32(coord_dist_v, dist_v);
                    coord_diff = _mm256_abs_epi32(coord_diff);
                    __m256 coord_diff_f = _mm256_cvtepi32_ps(coord_diff);

                    // float linear_penalty = coord_diff * sl;
                    __m256 linear_penalty_v = _mm256_mul_ps(coord_diff_f, _mm256_set1_ps(sl));

                    // from:
                    // https://github.com/IntelLabs/Trans-Omics-Acceleration-Library/blob/master/src/dynamic-programming/parallel_chaining_v2_22.h
                    auto mg_log2_avx2 = [](__m256i dd_v) -> __m256 // NB: this doesn't work when x<2
                    {
                        // -------- constant vectors --------------------
                        __m256i v255 = _mm256_set1_epi32(255);
                        __m256i v128 = _mm256_set1_epi32(128);
                        __m256i shift1 =  _mm256_set1_epi32(~(255 << 23));
                        __m256i shift2 =  _mm256_set1_epi32(127 << 23);
                            __m256 fc1 = _mm256_set1_ps(-0.34484843f);
                            __m256 fc2 = _mm256_set1_ps(2.02466578f);
                            __m256 fc3 = _mm256_set1_ps(0.67487759f);
                        // ---------------------------------------------

                        __m256 dd_v_f = _mm256_cvtepi32_ps(dd_v);
                            __m256i dd_v_i = _mm256_castps_si256(dd_v_f);

                        __m256i log2_v_i = _mm256_sub_epi32 (_mm256_and_si256( _mm256_srli_epi32(dd_v_i, 23), v255) , v128);

                        dd_v_i = _mm256_and_si256(dd_v_i, shift1);
                        dd_v_i = _mm256_add_epi32(dd_v_i, shift2);

                        dd_v_f = _mm256_castsi256_ps(dd_v_i);

                        __m256 t1 =_mm256_add_ps (_mm256_mul_ps(fc1, dd_v_f), fc2);
                        __m256 t2 = _mm256_sub_ps(_mm256_mul_ps(t1, dd_v_f), fc3);

                        __m256 log2_v_f = _mm256_add_ps(_mm256_cvtepi32_ps(log2_v_i), t2);

                        return log2_v_f;
                    };

                    // float log_penalty = log2(coord_diff + 1) * 0.5;
                    __m256 log_penalty_v = mg_log2_avx2(_mm256_add_epi32(coord_diff, _mm256_set1_epi32(1)));
                    log_penalty_v = _mm256_mul_ps(log_penalty_v, _mm256_set1_ps(0.5));

                    // cur_score -= linear_penalty + log_penalty;
                    __m256 gap_penalty_f = _mm256_add_ps(linear_penalty_v, log_penalty_v);
                    __m256i gap_penalty_v = _mm256_cvtps_epi32(gap_penalty_f);
                    __m256i dist_cutoff_mask = _mm256_cmpgt_epi32(coord_diff, _mm256_setzero_si256());
                    gap_penalty_v = _mm256_blendv_epi8(_mm256_setzero_si256(), gap_penalty_v, dist_cutoff_mask);
                    cur_score_v = _mm256_sub_epi32(cur_score_v, gap_penalty_v);

                    // if (cur_score >= score) {
                    //     score = cur_score;
                    //     backtrace[j] = i;
                    // }
                    __m256i old_scores_v = _mm256_i32gather_epi32(&dp_table[j].chain_score, _mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0), 4);
                    __m256i score_neg_mask = _mm256_cmpgt_epi32(old_scores_v, cur_score_v);
                    __m256i mask = _mm256_andnot_si256(score_neg_mask, dist_mask);

                    cur_score_v = _mm256_blendv_epi8(old_scores_v, cur_score_v, mask);
                    _mm256_maskstore_epi32(&backtrace[j], mask, i_v);

                    // note: _mm256_i32scatter_epi32 not supported in AVX2
                    score_t cur_scores[8] __attribute__((aligned(32)));
                    _mm256_store_si256((__m256i*)cur_scores, cur_score_v);
                    auto *dp_table_o = &dp_table[j];
                    dp_table_o[0].chain_score = cur_scores[0];
                    dp_table_o[1].chain_score = cur_scores[1];
                    dp_table_o[2].chain_score = cur_scores[2];
                    dp_table_o[3].chain_score = cur_scores[3];
                    dp_table_o[4].chain_score = cur_scores[4];
                    dp_table_o[5].chain_score = cur_scores[5];
                    dp_table_o[6].chain_score = cur_scores[6];
                    dp_table_o[7].chain_score = cur_scores[7];
                }

                if (_mm256_movemask_epi8(coord_neg_mask))
                    break;

                j_v = _mm256_add_epi32(j_v, _mm256_set1_epi32(8));
            }
#else
            for (size_t j = i + 1; j < it_end; ++j) {
                auto &[label, coord, clipping, end, score, seed_i] = dp_table[j];
                assert(label == prev_label);
                if (coord_cutoff > coord)
                    break;

                int32_t dist = prev_clipping - clipping;
                int32_t coord_dist = prev_coord - coord;
                if (dist > 0 && std::max(dist, coord_dist) < query_size) {
                    score_t match = std::min({ dist, coord_dist, end - clipping });
                    score_t cur_score = prev_score + match;
                    if (coord_dist != dist) {
                        float coord_diff = std::abs(coord_dist - dist);
                        float linear_penalty = coord_diff * sl;
                        float log_penalty = log2(coord_diff + 1) * 0.5;
                        cur_score -= linear_penalty + log_penalty;
                    }
                    if (cur_score >= score) {
                        score = cur_score;
                        backtrace[j] = i;
                    }
                }
            }
#endif
        }
    }

    return std::make_tuple(std::move(dp_table), std::move(backtrace), num_seeds, num_nodes);
}

Alignment::LabelChangeScorer make_label_change_scorer(const IDBGAligner &aligner);

template <class AlignmentCompare, class LabelChangeScorer>
void construct_alignment_chain(const IDBGAligner &aligner,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> *best_score,
                               const std::function<void(Alignment&&)> &callback,
                               LabelChangeScorer &get_label_change_score);

template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(const IDBGAligner &aligner,
                                        std::vector<Alignment>&& alignments,
                                        std::string_view query,
                                        std::string_view rc_query) {
    const DBGAlignerConfig &config = aligner.get_config();

    if (alignments.size() < 2 || !config.post_chain_alignments) {
        DEBUG_LOG("Too few alignments found, nothing to chain.");
        return std::move(alignments);
    }

    for (const auto &a : alignments) {
        if (a.label_coordinates.size())
            throw std::runtime_error("Post-chaining alignments with coordinates not supported");
    }

    DBGAlignerConfig no_chain_config { config };
    no_chain_config.post_chain_alignments = false;
    AlignmentAggregator<AlignmentCompare> aggregator(no_chain_config);

    alignments.erase(std::remove_if(alignments.begin(), alignments.end(), [&](Alignment &a) {
        // TODO: handle offset case later
        if ((!a.get_clipping() && !a.get_end_clipping()) || a.get_offset()) {
            aggregator.add_alignment(std::move(a));
            return true;
        }

        return false;
    }), alignments.end());

    if (alignments.empty())
        return aggregator.get_alignments();

    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_orientation(),
                               a.get_clipping() + a.get_query_view().size(),
                               a.get_clipping(),
                               b.get_score(),
                               a.get_sequence().size())
            < std::make_tuple(b.get_orientation(),
                              b.get_clipping() + b.get_query_view().size(),
                              b.get_clipping(),
                              a.get_score(),
                              b.get_sequence().size());
    });

    logger->trace("Chaining alignments:\n{}", fmt::join(alignments, "\t\n"));

    auto get_label_change_score = make_label_change_scorer(aligner);

    auto run = [&](std::string_view this_query, auto begin, auto end) {
        std::vector<score_t> best_score(this_query.size() + 1, 0);
        for (auto it = begin; it != end; ++it) {
            size_t end_pos = it->get_query_view().data() + it->get_query_view().size()
                                - this_query.data();
            if (it->get_score() > best_score[end_pos]) {
                best_score[end_pos] = it->get_score();
                construct_alignment_chain<AlignmentCompare>(
                    aligner, this_query, Alignment(*it), it + 1, end, &best_score,
                    [&](Alignment&& chain) { aggregator.add_alignment(std::move(chain)); },
                    get_label_change_score
                );
            }
        }
    };

    // recursively construct chains
    auto split_it = std::find_if(alignments.begin(), alignments.end(),
                                 [](const auto &a) { return a.get_orientation(); });
    run(query, alignments.begin(), split_it);
    run(rc_query, split_it, alignments.end());

    return aggregator.get_alignments();
}

// TODO: rewrite this to not use recursion
template <class AlignmentCompare>
void construct_alignment_chain(const IDBGAligner &aligner,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> *best_score,
                               const std::function<void(Alignment&&)> &callback,
                               Alignment::LabelChangeScorer &get_label_change_score) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.get_query_view().data();
    const char *chain_end = chain.get_query_view().data() + chain.get_query_view().size();
    if (begin == end || chain_end == query.data() + query.size()) {
        callback(std::move(chain));
        return;
    }

    const DBGAlignerConfig &config = aligner.get_config();
    const DeBruijnGraph &graph = aligner.get_graph();
    size_t node_overlap = graph.get_k() - 1;

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        const char *next_begin = it->get_query_view().data();
        const char *next_end = it->get_query_view().data() + it->get_query_view().size();

        assert(chain_begin - chain.get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.get_orientation());

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        Alignment next_chain = chain;
        Alignment aln = *it;

        ssize_t query_overlap = chain_end - next_begin;

        if (query_overlap <= 0) {
            // no overlap
            aln.insert_gap_prefix(-query_overlap, node_overlap, config);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t num_trimmed_from_last_match
                = aln.trim_query_prefix(query_overlap, node_overlap, config);

            if (aln.empty())
                continue;

            // try to connect the two alignments
            size_t gap_fill_length = graph.get_k() - aln.get_offset();
            std::vector<node_index> nodes;
            std::string_view spelling(aln.get_sequence().data(), gap_fill_length);
            if (aln.get_sequence().size() >= gap_fill_length) {
                graph.traverse(next_chain.get_nodes().back(),
                               spelling.data(), spelling.data() + gap_fill_length,
                               [&](node_index next) { nodes.emplace_back(next); });
            }

            bool found = false;
            if (nodes.size() == gap_fill_length && nodes.back() == aln.get_nodes()[0]) {
                found = true;

                // if the gap is too long, extend aln's offset
                if (gap_fill_length > 1) {
                    nodes.pop_back();
                    std::vector<size_t> node_columns;
                    std::vector<score_t> node_scores;
                    if (const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner)) {
                        AnnotationBuffer &buffer = labeled_aligner->get_annotation_buffer();
                        buffer.queue_path(std::vector<node_index>(nodes));
                        buffer.fetch_queued_annotations();

                        const auto *cur_columns = &aln.get_columns(0);

                        node_columns.reserve(nodes.size());
                        node_scores.reserve(nodes.size());
                        bool had_diff = false;
                        auto jt = spelling.rbegin();
                        for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
                            assert(jt != spelling.rend());
                            const auto *next_columns = buffer.get_labels(*it);
                            assert(next_columns);
                            Vector<Alignment::Column> inter;
                            Vector<Alignment::Column> diff;
                            utils::set_intersection_difference(next_columns->begin(),
                                                               next_columns->end(),
                                                               cur_columns->begin(),
                                                               cur_columns->end(),
                                                               std::back_inserter(inter),
                                                               std::back_inserter(diff));
                            if (inter.size()) {
                                node_columns.emplace_back(buffer.cache_column_set(
                                    std::move(inter)
                                ));
                                node_scores.emplace_back(0);
                            } else {
                                assert(diff.size());
                                score_t label_change_score
                                    = get_label_change_score(*it, std::string_view((jt + graph.get_k() + 1).base(), graph.get_k()), *jt, *cur_columns, diff);
                                if (label_change_score == config.ninf) {
                                    found = false;
                                    break;
                                }

                                node_columns.emplace_back(buffer.cache_column_set(
                                    std::move(diff)
                                ));
                                node_scores.emplace_back(label_change_score);
                                had_diff = true;
                            }
                            cur_columns = &buffer.get_cached_column_set(node_columns.back());
                            ++jt;
                        }
                        if (found) {
                            std::reverse(node_columns.begin(), node_columns.end());
                            if (had_diff) {
                                std::reverse(node_scores.begin(), node_scores.end());
                            } else {
                                node_scores.clear();
                            }
                        }
                    }

                    if (found) {
                        aln.extend_offset(std::move(nodes), std::move(node_columns),
                                          std::move(node_scores));
                    }
                }
            }

            // if connection failed, then insert dummy graph nodes to connect the two
            if (!found) {
                size_t overlap = std::min(
                    static_cast<size_t>((chain.get_cigar().data().end() - 2)->second),
                    num_trimmed_from_last_match
                );

                if (aln.get_sequence().size() <= node_overlap
                        || (aln.get_cigar().data().begin()
                            + static_cast<bool>(aln.get_clipping()))->first != Cigar::MATCH) {
                    continue;
                }

                assert(aln.get_query_view().data()
                    == chain.get_query_view().data() + chain.get_query_view().size());

                assert(overlap < node_overlap);
                aln.insert_gap_prefix(-overlap, node_overlap, config);
            } else {
                aln.trim_clipping();
            }
        }

        assert(!aln.empty());

        // use append instead of splice because any clipping in aln represents
        // internally clipped characters
        next_chain.trim_end_clipping();
        bool changed = next_chain.append(
                std::move(aln), get_label_change_score, config.label_change_union);
        assert(next_chain.is_valid(graph, &config));
        score_t next_score = next_chain.get_score();
        if (next_score <= (*best_score)[next_end - query.data()])
            continue;

        (*best_score)[next_end - query.data()] = next_score;
        if (next_chain.size()) {
            construct_alignment_chain<AlignmentCompare>(
                    aligner, query, std::move(next_chain),
                    it + 1, end, best_score, callback, get_label_change_score);
            called |= changed;
        }
    }

    if (!called)
        callback(std::move(chain));
}

typedef boss::BOSS::edge_index edge_index;
typedef std::pair<edge_index, edge_index> BOSSRange;
typedef boss::BOSS::TAlphabet TAlphabet;

void call_boss_edges(const DBGSuccinct &dbg_succ,
                     const BOSSRange &node_range,
                     const std::function<void(node_index, const SmallVector<node_index>&)> &callback) {
    const boss::BOSS &boss = dbg_succ.get_boss();
    std::vector<BOSSRange> next_ranges(boss.alph_size);
    std::vector<BOSSRange> next_is(boss.alph_size);
    boss.call_tightened_ranges(node_range.first, node_range.second, [&](edge_index first, edge_index second, TAlphabet s) {
        next_ranges[s] = std::make_pair(first, second);
        next_is[s].second = boss.succ_last(first);
        next_is[s].first = boss.pred_last(next_is[s].second - 1) + 1;
    });

    TAlphabet s = boss.get_W(node_range.first);
    edge_index succ_W = 0;
    edge_index succ_Wp = 0;
    if (node_range.first + 1 < boss.get_W().size()) {
        succ_W = s < boss.alph_size ? node_range.first + 1 : boss.succ_W(node_range.first + 1, s + boss.alph_size);
        succ_Wp = s >= boss.alph_size ? node_range.first + 1 : boss.succ_W(node_range.first + 1, s);
    }
    for (edge_index node_w = node_range.first; node_w <= node_range.second; ++node_w) {
        if (!s) {
            if (node_w + 1 <= node_range.second) {
                s = boss.get_W(node_w + 1);
                if (s < boss.alph_size) {
                    succ_W = boss.succ_W(node_w + 1, s);
                } else {
                    succ_Wp = boss.succ_W(node_w + 1, s);
                }
            }

            continue;
        }

        TAlphabet last_c = s % boss.alph_size;

        auto &[i_start, i] = next_is[last_c];
        assert(boss.fwd(node_w, last_c) == i);
        assert(i <= next_ranges[last_c].second);

        if (node_index n = dbg_succ.boss_to_kmer_index(node_w)) {
            SmallVector<node_index> edges(boss.alph_size);
            for (edge_index next_edge = i_start; next_edge <= i; ++next_edge) {
                if (node_index m = dbg_succ.boss_to_kmer_index(next_edge))
                    edges[boss.get_W(next_edge) % boss.alph_size] = m;
            }
            callback(n, edges);
        }

        assert(node_w + 1 < boss.get_W().size() || node_w == node_range.second);

        if (node_w + 1 < boss.get_W().size()) {
            if (succ_Wp >= succ_W) {
                i_start = i + 1;
                i = boss.succ_last(i_start);
            }
            s = boss.get_W(node_w + 1);
            if (s < boss.alph_size) {
                succ_W = boss.succ_W(node_w + 1, s);
            } else {
                succ_Wp = boss.succ_W(node_w + 1, s);
            }
        }
    }
}

typedef SmallVector<std::pair<node_index, TAlphabet>> Edges;
typedef SmallVector<std::pair<node_index, Edges>> ParallelEdges;

std::pair<ParallelEdges, size_t> get_parallel_edges(const DeBruijnGraph &graph,
                                                    const DBGSuccinct &dbg_succ,
                                                    const CanonicalDBG *canonical,
                                                    std::string_view seq,
                                                    size_t max_edits) {
    const auto &boss = dbg_succ.get_boss();
    auto encoded = boss.encode(seq);
    std::vector<std::tuple<BOSSRange, size_t, size_t, size_t>> traversal_stack;
    ParallelEdges result;

    std::vector<TAlphabet> encoded_rc;
    if (canonical) {
        std::string seq_rc(seq);
        ::reverse_complement(seq_rc.begin(), seq_rc.end());
        encoded_rc = boss.encode(seq_rc);
    }

    for (TAlphabet s = 1; s < boss.alph_size; ++s) {
        edge_index first = boss.get_F(s) + 1 < boss.get_W().size()
             ? boss.get_F(s) + 1
             : boss.get_W().size(); // lower bound
        edge_index last = s + 1 < boss.alph_size
             ? boss.get_F(s + 1)
             : boss.get_W().size() - 1;
        if (first <= last) {
            traversal_stack.emplace_back(std::make_pair(first, last), 1,
                                         s != encoded[0],
                                         encoded_rc.size() ? s != encoded_rc[0] : max_edits + 1);
        }
    }

    size_t traversal_steps = 2;
    while (traversal_stack.size()) {
        auto [range, i, edits, edits_rc] = traversal_stack.back();
        traversal_stack.pop_back();

        if (i + 1 == encoded.size()) {
            auto add_parallel_node = [&](node_index parallel_node) {
                bool found = false;
                result.emplace_back(parallel_node, Edges{});
                graph.call_outgoing_kmers(parallel_node, [&](node_index next, char c) {
                    TAlphabet s = boss.encode(c);
                    found = true;
                    result.back().second.emplace_back(next, s);
                });

                if (!found)
                    result.back().second.emplace_back(DeBruijnGraph::npos, 0);
            };

            if (edits == max_edits) {
                if (encoded.back() != boss.alph_size) {
                    if (auto edge = boss.pick_edge(range.second, encoded.back())) {
                        if (auto parallel_node = dbg_succ.boss_to_kmer_index(edge))
                            add_parallel_node(parallel_node);
                    }
                }
            } else if (edits < max_edits) {
                boss.call_outgoing(range.second, [&](edge_index pedge) {
                    if (auto parallel_node = dbg_succ.boss_to_kmer_index(pedge))
                        add_parallel_node(parallel_node);
                });
            }

            if (edits_rc == max_edits) {
                if (encoded_rc.back() != boss.alph_size) {
                    if (auto edge = boss.pick_edge(range.second, encoded_rc.back())) {
                        if (auto parallel_node = dbg_succ.boss_to_kmer_index(edge))
                            add_parallel_node(canonical->reverse_complement(parallel_node));
                    }
                }
            } else if (edits_rc < max_edits) {
                boss.call_outgoing(range.second, [&](edge_index pedge) {
                    if (auto parallel_node = dbg_succ.boss_to_kmer_index(pedge))
                        add_parallel_node(canonical->reverse_complement(parallel_node));
                });
            }

            continue;
        }

        if (range.first == range.second) {
            ++traversal_steps;
            if (TAlphabet s = boss.get_W(range.second) % boss.alph_size) {
                size_t next_edits = edits <= max_edits ? edits + (s != encoded[i]) : max_edits + 1;
                size_t next_edits_rc = edits_rc <= max_edits ? edits_rc + (s != encoded_rc[i]) : max_edits + 1;
                if (next_edits <= max_edits || next_edits_rc <= max_edits) {
                    range.second = boss.fwd(range.second, s);
                    range.first = boss.pred_last(range.second - 1) + 1;
                    traversal_stack.emplace_back(std::move(range), i + 1, next_edits, next_edits_rc);
                }
            }
        } else if (edits < max_edits || edits_rc < max_edits) {
            traversal_steps += 2 * (boss.alph_size - 1);
            boss.call_tightened_ranges(range.first, range.second, [&](edge_index first, edge_index second, TAlphabet s) {
                size_t next_edits = edits <= max_edits ? edits + (s != encoded[i]) : max_edits + 1;
                size_t next_edits_rc = edits_rc <= max_edits ? edits_rc + (s != encoded_rc[i]) : max_edits + 1;
                if (next_edits <= max_edits || next_edits_rc <= max_edits)
                    traversal_stack.emplace_back(std::make_pair(first, second), i + 1, next_edits, next_edits_rc);
            });
        } else {
            if (edits <= max_edits) {
                auto next_range = range;
                traversal_steps += 2;
                if (boss.tighten_range(&next_range.first, &next_range.second, encoded[i])) {
                    traversal_stack.emplace_back(std::move(next_range), i + 1, edits,
                                                 edits_rc <= max_edits ? edits_rc + (encoded[i] != encoded_rc[i]) : max_edits + 1);
                }
            }

            if (edits_rc <= max_edits && encoded[i] != encoded_rc[i]) {
                auto next_range = range;
                traversal_steps += 2;
                if (boss.tighten_range(&next_range.first, &next_range.second, encoded_rc[i])) {
                    traversal_stack.emplace_back(std::move(next_range), i + 1,
                                                 edits <= max_edits ? edits + (encoded[i] != encoded_rc[i]) : max_edits + 1,
                                                 edits_rc);
                }
            }
        }
    }

    return { std::move(result), traversal_steps };
}

Alignment::LabelChangeScorer
make_label_change_scorer(const IDBGAligner &aligner) {
    const DBGAlignerConfig &config = aligner.get_config();

    const DeBruijnGraph *graph = &aligner.get_graph();
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(graph);
    const DeBruijnGraph *base_graph = canonical
        ? &canonical->get_graph() : graph;

    const DBGSuccinct *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);
    const HLLWrapper<> *hll_wrapper = dbg_succ ? dbg_succ->get_extension_threadsafe<HLLWrapper<>>() : nullptr;
    if (!dbg_succ || !hll_wrapper || config.label_change_score != config.ninf) {
        return [ninf=config.ninf,score=config.label_change_score](auto,
                                                                  auto,
                                                                  char c,
                                                                  const auto &ref_columns,
                                                                  const auto &diff_columns) {
            return c != boss::BOSS::kSentinel && ref_columns.size() && diff_columns.size()
                ? score : ninf;
        };
    }

    if (!config.label_change_edit_distance) {
        return [&config,ninf=config.ninf,hll_wrapper](auto, auto, char c,
                                                      const auto &ref_columns,
                                                      const auto &diff_columns) -> score_t {
            if (c == boss::BOSS::kSentinel || ref_columns.empty() || diff_columns.empty())
                return ninf;

            const auto &hll = hll_wrapper->data();
            auto [union_est, inter_est]
                = hll.estimate_column_union_intersection_cardinality(ref_columns,
                                                                     diff_columns);
            assert(union_est > 0.0);
            return inter_est > 0.0
                ? score_t(log2(inter_est) - log2(union_est)) * config.score_matrix[c][c]
                : ninf;
        };
    }

    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    if (!labeled_aligner)
        return [ninf=config.ninf](auto&&...) { return ninf; };

    const auto *node_rc = dbg_succ->get_extension_threadsafe<NodeRC>();

    constexpr double dninf = std::numeric_limits<double>::min();

    return [&config,max_edits=config.label_change_edit_distance,ninf=config.ninf,
            labeled_aligner,graph,dbg_succ,canonical,node_rc,hll_wrapper,
            cache=tsl::hopscotch_map<node_index,ParallelEdges>()](node_index node,
                                                                  std::string_view seq,
                                                                  char c,
                                                                  const auto &ref_columns,
                                                                  const auto &diff_columns) mutable -> score_t {
        if (c == boss::BOSS::kSentinel || ref_columns.empty() || diff_columns.empty())
            return ninf;

        const boss::BOSS &boss = dbg_succ->get_boss();
        TAlphabet s = boss.encode(c);

        if (s == boss.alph_size)
            return ninf;

        size_t num_nodes_covered = 0;
        size_t traversal_steps = 0;
        auto [parallel_edges, inserted] = cache.try_emplace(node, ParallelEdges{});
        AnnotationBuffer &annotation_buffer = labeled_aligner->get_annotation_buffer();
        if (inserted) {
            std::tie(parallel_edges.value(), traversal_steps)
                = get_parallel_edges(*graph, *dbg_succ, canonical, seq, max_edits);
            std::vector<node_index> nodes_to_queue;
            for (const auto &[parallel_node, edges] : parallel_edges->second) {
                nodes_to_queue.push_back(parallel_node);
                for (const auto &[next, s] : edges) {
                    if (next)
                        nodes_to_queue.push_back(next);
                }
            }

            num_nodes_covered = nodes_to_queue.size();
            annotation_buffer.queue_path(std::move(nodes_to_queue));
            annotation_buffer.fetch_queued_annotations();
        }

        assert(parallel_edges->size());

        size_t total = 0;
        size_t matches = 0;
        for (const auto &[n, edges] : parallel_edges->second) {
            auto [n_labels, n_coords] = annotation_buffer.get_labels_and_coords(n);
            assert(n_labels);

            bool found = false;
            for (const auto &[m, s_cur] : edges) {
                if (s != s_cur)
                    continue;

                found = true;

                if (!m) {
                    if (n_coords) {
                        for (const auto &coords : *n_coords) {
                            total += coords.size();
                        }
                    } else {
                        total += n_labels->size();
                    }
                    continue;
                }

                auto [m_labels, m_coords] = annotation_buffer.get_labels_and_coords(m);
                assert(m_labels);
                if (m_coords) {
                    utils::match_indexed_values(
                        n_labels->begin(), n_labels->end(), n_coords->begin(),
                        m_labels->begin(), m_labels->end(), m_coords->begin(),
                        [&](auto, const auto &coords, const auto &other_coords) {
                            total += coords.size();
                            matches += std::min(coords.size(), other_coords.size());
                        },
                        [&](auto, const auto &coords) { total += coords.size(); },
                        [](auto&&...) {}
                    );
                } else {
                    total += n_labels->size();
                    matches += utils::count_intersection(n_labels->begin(), n_labels->end(),
                                                         m_labels->begin(), m_labels->end());
                }
            }

            if (!found) {
                if (n_coords) {
                    for (const auto &coords : *n_coords) {
                        total += coords.size();
                    }
                } else {
                    total += n_labels->size();
                }
            }
        }

        double pseudocount = 0.5;
        double outscore = log2(static_cast<double>(matches + pseudocount))
                            - log2(static_cast<double>(total + pseudocount));

        const auto &hll = hll_wrapper->data();
        double best_jaccard = 0.0;
        for (auto c : ref_columns) {
            for (auto d : diff_columns) {
                auto [union_est, inter_est]
                    = hll.estimate_column_union_intersection_cardinality(c, d);
                double jaccard = inter_est > 0.0 ? inter_est / union_est : 0.0;
                best_jaccard = std::max(best_jaccard, jaccard);
                DEBUG_LOG("Label Jaccard {} -> {}: {}",
                    annotation_buffer.get_annotator().get_label_encoder().decode(c)
                    annotation_buffer.get_annotator().get_label_encoder().decode(d)
                    jaccard);
            }
        }

        if (best_jaccard > 0.0) {
            DEBUG_LOG("Best Jaccard: {}", log2(best_jaccard));
            outscore += log2(best_jaccard);
        } else {
            outscore = dninf;
        }

        score_t ret_val = outscore != dninf
            ? std::floor(outscore) * config.score_matrix[c][c]
            : ninf;

        logger->trace("Label change score: {}, {} / {}; "
                      "Covered {} / {} nodes after {} traversal steps.",
                      ret_val, matches, total, num_nodes_covered,
                      dbg_succ->num_nodes(), traversal_steps);
        return ret_val;
    };
}

template
std::vector<Alignment> chain_alignments<LocalAlignmentLess>(const IDBGAligner&,
                                                            std::vector<Alignment>&&,
                                                            std::string_view,
                                                            std::string_view);

} // namespace align
} // namespace graph
} // namespace mtg
