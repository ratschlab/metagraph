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

template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(const IDBGAligner &aligner,
                                        std::vector<Alignment>&& alignments) {
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
        // TODO: handle offset and mixed label cases later
        if ((!a.get_clipping() && !a.get_end_clipping()) || a.get_offset() || a.label_column_diffs.size()) {
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

    const HLLWrapper<> *hll_wrapper = aligner.get_graph().get_extension_threadsafe<HLLWrapper<>>();
    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    auto get_label_change_scores = [&](Alignment::Columns a_col, Alignment::Columns b_col)
            -> std::vector<std::pair<Alignment::Columns, score_t>> {
        if (a_col == b_col)
            return { std::make_pair(a_col, 0) };

        const auto &a_cols = labeled_aligner->get_annotation_buffer().get_cached_column_set(a_col);
        const auto &b_cols = labeled_aligner->get_annotation_buffer().get_cached_column_set(b_col);
        Vector<Alignment::Column> inter;
        Vector<Alignment::Column> diff;
        utils::set_intersection_difference(b_cols.begin(), b_cols.end(),
                                           a_cols.begin(), a_cols.end(),
                                           std::back_inserter(inter),
                                           std::back_inserter(diff));
        if (inter.size()) {
            return { std::make_pair(
                labeled_aligner->get_annotation_buffer().cache_column_set(std::move(inter)), 0
            ) };
        }

        if (!hll_wrapper || config.label_change_score != config.ninf) {
            return { std::make_pair(labeled_aligner->get_annotation_buffer().cache_column_set(std::move(diff)),
                                    config.label_change_score) };
        }

        tsl::hopscotch_map<score_t, VectorSet<Alignment::Column>> scores;
        for (Alignment::Column d : diff) {
            for (Alignment::Column c : a_cols) {
                auto [a_size, b_size, union_size]
                    = hll_wrapper->data().estimate_column_sizes_union_cardinality(c, d);

                uint64_t size_sum = a_size + b_size;
                if (union_size < size_sum)
                    continue;

                double dbsize = b_size;
                score_t label_change_score = (log2(std::min(dbsize, size_sum - union_size)) - log2(dbsize));
                scores[label_change_score].emplace(d);
            }
        }

        std::vector<std::pair<size_t, score_t>> results;

        // TODO: a structured for-loop here crashes g++-8
        for (const auto &score_diff : scores) {
            const auto &diff = score_diff.second;
            assert(std::is_sorted(diff.begin(), diff.end()));
            results.emplace_back(labeled_aligner->get_annotation_buffer().cache_column_set(diff.begin(), diff.end()),
                                 score_diff.first);
        }

        return results;
    };

    typedef std::tuple<score_t /* best score */,
                       score_t /* extra score */,
                       ssize_t /* gap */,
                       size_t /* prefix trim */,
                       size_t /* last table i */,
                       Alignment::Columns /* last labels */,
                       size_t /* prev suffix trim */> TableVal;
    typedef tsl::hopscotch_map<Alignment::Columns, TableVal> Table;
    std::vector<VectorMap<size_t, Table>> chain_table(alignments.size());
    const DeBruijnGraph &graph = aligner.get_graph();
    ssize_t node_overlap = graph.get_k() - 1;

    auto last_is_indel = [](const Cigar &a) {
        auto it = a.data().rbegin();
        if (it->first == Cigar::CLIPPED)
            ++it;

        return it != a.data().rend() && (it->first == Cigar::INSERTION || it->first == Cigar::DELETION);
    };

    auto first_is_indel = [](const Cigar &a) {
        auto it = a.data().begin();
        if (it->first == Cigar::CLIPPED)
            ++it;

        return it != a.data().end() && (it->first == Cigar::INSERTION || it->first == Cigar::DELETION);
    };


    for (size_t i = 0; i < alignments.size(); ++i) {
        Alignment a = alignments[i];
        std::vector<std::pair<size_t, Table>> vals;
        ssize_t max_suffix_trim = 0;
        const char *a_end = a.get_query_view().data() + a.get_query_view().size();
        for (size_t j = i + 1; j < alignments.size(); ++j) {
            max_suffix_trim = std::max(a_end - alignments[j].get_query_view().data(),
                                       max_suffix_trim);
        }
        for (ssize_t trim = 0; a.size() && trim <= max_suffix_trim; a.trim_query_suffix(1, config), ++trim) {
            if (last_is_indel(a.get_cigar()))
                continue;

            vals.emplace_back(trim, Table{});
            vals.back().second.try_emplace(a.label_columns, TableVal {
                a.get_score(), 0, 0, 0, std::numeric_limits<size_t>::max(), 0, 0
            });
        }
        chain_table[i].insert(vals.begin(), vals.end());
        assert(std::is_sorted(chain_table[i].begin(), chain_table[i].end(),
                              [&](const auto &a, const auto &b) {
                                  return a.first < b.first;
                              }));
        assert(chain_table[i].size() == vals.size());
    }

    for (size_t i = 0; i < alignments.size() - 1; ++i) {
        if (alignments[i].get_orientation() != alignments[i + 1].get_orientation())
            continue;

        for (const auto &[a_suffix_trim, tab] : chain_table[i]) {
            Alignment a = alignments[i];
            a.trim_query_suffix(a_suffix_trim, config);
            assert(a.size());
            assert(a.is_valid(graph, &config));
            assert(!last_is_indel(a.get_cigar()));

            for (const auto &[cur_columns, vals] : tab) {
                const auto &[a_score, a_extra_score, a_gap, a_prefix_trim,
                             a_last_table_i, a_last_labels, a_last_suffix_trim] = vals;
                const char *chain_begin = a.get_query_view().data() + a_prefix_trim;
                const char *chain_end = a.get_query_view().data() + a.get_query_view().size();

                for (size_t j = i + 1; j < alignments.size(); ++j) {
                    if (a.get_orientation() != alignments[j].get_orientation())
                        break;

                    assert(std::is_sorted(chain_table[j].begin(), chain_table[j].end(),
                                          [&](const auto &a, const auto &b) {
                                              return a.first < b.first;
                                          }));

                    size_t last_b_suffix_trim = 0;
                    Alignment b_base = alignments[j];
                    for (auto it = chain_table[j].begin(); it != chain_table[j].end(); ++it) {
                        size_t b_suffix_trim = it->first;
                        const char *next_begin = b_base.get_query_view().data();
                        const char *next_end = next_begin + b_base.get_query_view().size() - b_suffix_trim + last_b_suffix_trim;
                        Alignment::Columns b_col = b_base.label_columns;
                        if (next_begin < chain_begin && cur_columns == b_col)
                            continue;

                        if (next_end <= chain_end)
                            continue;

                        size_t b_prefix_trim = std::max((ptrdiff_t)0, chain_begin - next_begin);
                        next_begin += b_prefix_trim;

                        ssize_t gap = next_begin - chain_end;

                        assert(b_suffix_trim >= last_b_suffix_trim);
                        if (b_suffix_trim > last_b_suffix_trim) {
                            b_base.trim_query_suffix(b_suffix_trim - last_b_suffix_trim, config);
                            last_b_suffix_trim = b_suffix_trim;
                            assert(!last_is_indel(b_base.get_cigar()));
                        }

                        assert(b_base.size());
                        assert(a.get_query_view().data() - a.get_clipping()
                                == next_begin - b_base.get_clipping());
                        assert(next_begin + b_base.get_query_view().size() == next_end);
                        assert(b_base.is_valid(graph, &config));

                        auto label_change_scores = get_label_change_scores(cur_columns, b_col);

                        score_t gap_score = 0;
                        score_t b_score = b_base.get_score();
                        char c = b_base.get_sequence()[0];

                        if (gap >= 0) {
                            // no overlap
                            gap_score = config.gap_opening_penalty;
                            if (gap > 0)
                                gap_score += config.gap_opening_penalty + (gap - 1) * config.gap_extension_penalty;

                        } else {
                            // overlap
                            Alignment b = b_base;

                            b_prefix_trim -= gap;

                            // TODO: this is expensive
                            b.trim_query_prefix(b_prefix_trim, node_overlap, config);

                            if (b.empty() || first_is_indel(b.get_cigar()))
                                continue;

                            assert(b_base.get_sequence().size() >= b.get_sequence().size());
                            std::string_view b_prefix(
                                b_base.get_sequence().data(),
                                b_base.get_sequence().size() - b.get_sequence().size()
                            );
                            assert(b_prefix.data() + b_prefix.size()
                                    <= b_base.get_sequence().data() + b_base.get_sequence().size());

                            if (b_prefix.size() > static_cast<size_t>(node_overlap))
                                b_prefix.remove_prefix(b_prefix.size() - node_overlap);

                            auto [a_mm, b_mm] = std::mismatch(a.get_sequence().rbegin(),
                                                              a.get_sequence().rend(),
                                                              b_prefix.rbegin(),
                                                              b_prefix.rend());
                            gap = -(b_mm - b_prefix.rbegin());
                            assert(gap <= 0);

                            if (gap == 0) {
                                b.trim_offset();
                                if (b.get_offset())
                                    continue;

                                gap_score = config.gap_opening_penalty;
                            } else if (node_overlap + gap > 0) {
                                ssize_t gapoffset = gap + static_cast<ssize_t>(b.get_offset());
                                if (gapoffset < 0 || static_cast<ssize_t>(b.size()) <= gapoffset)
                                    continue;
                            }

                            b_score = b.get_score();
                            c = b.get_sequence()[0];

#ifndef NDEBUG
                            Alignment prev = a;
                            Alignment cur = b;
                            prev.trim_query_prefix(a_prefix_trim, node_overlap, config);

                            if (gap == 0) {
                                cur.trim_offset();
                                assert(!cur.get_offset());
                                cur.trim_clipping();
                            }

                            if (gap == 0 || node_overlap + gap > 0) {
                                cur.insert_gap_prefix(gap, node_overlap, config);
                            } else {
                                cur.trim_clipping();
                            }

                            prev.trim_end_clipping();
                            prev.append(std::move(cur));
                            assert(prev.size());
                            assert(prev.is_valid(graph, &config));
#endif
                            // std::cerr << "sss\t" << a_score + b_score + gap_score << "\t" << a << "\t" << b << std::endl;
                        }

                        Table &b_tab = it.value();
                        score_t next_base_score = a_score + b_score + gap_score;

                        score_t lambda = config.score_matrix[c][c];
                        for (auto &[cols, lc_score] : label_change_scores) {
                            lc_score *= lambda;
                            score_t next_score = next_base_score + lc_score;
                            lc_score += gap_score;
                            auto find = b_tab.find(cols);
                            if (find == b_tab.end()) {
                                b_tab.try_emplace(cols, TableVal{
                                    next_score, lc_score, gap, b_prefix_trim, i, cur_columns, a_suffix_trim
                                });
                            } else if (next_score > std::get<0>(find->second)) {
                                find.value() = std::tie(
                                    next_score, lc_score, gap, b_prefix_trim, i, cur_columns, a_suffix_trim
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    DEBUG_LOG("Done chaining");

    std::vector<std::tuple<score_t, size_t, Alignment::Columns, size_t>> indices;
    sdsl::bit_vector used(chain_table.size(), false);
    for (size_t i = 0; i < chain_table.size(); ++i) {
        for (const auto &[suffix_trim, tab] : chain_table[i]) {
            for (const auto &[cols, vals] : tab) {
                indices.emplace_back(std::get<0>(vals), suffix_trim, cols, i);
            }
        }
    }
    std::sort(indices.begin(), indices.end());
    for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
        auto [start_score, suffix_trim, cols, i] = *it;
        if (used[i])
            continue;

        used[i] = true;
        Alignment cur = alignments[i];
        auto [score, extra_score, gap, prefix_trim,
              last_i, last_cols, last_suffix_trim] = chain_table[i][suffix_trim][cols];

        cur.trim_query_suffix(suffix_trim, config);
        assert(cur.size());

        cur.trim_query_prefix(prefix_trim, node_overlap, config);
        assert(cur.size());

        assert(cur.is_valid(graph, &config));

        while (last_i != std::numeric_limits<size_t>::max()) {
            i = last_i;
            suffix_trim = last_suffix_trim;
            cols = last_cols;
            ssize_t last_gap = gap;
            std::tie(score, extra_score, gap, prefix_trim,
                     last_i, last_cols, last_suffix_trim) = chain_table[i][suffix_trim][cols];
            used[i] = true;

            Alignment prev = alignments[i];
            prev.trim_query_suffix(suffix_trim, config);
            assert(prev.size());

            prev.trim_query_prefix(prefix_trim, node_overlap, config);
            assert(prev.size());

            if (last_gap >= 0 || node_overlap + last_gap > 0) {
                if (last_gap == 0 && cur.get_offset()) {
                    cur.trim_offset();
                    assert(!cur.get_offset());
                    cur.trim_clipping();
                }
                cur.insert_gap_prefix(last_gap, node_overlap, config);
            } else {
                cur.trim_clipping();
            }

            prev.trim_end_clipping();
            prev.append(std::move(cur), extra_score);
            assert(prev.size());
            std::swap(cur, prev);
            assert(cur.is_valid(graph, &config));
        }

        assert(cur.get_score() == start_score);

        aggregator.add_alignment(std::move(cur));
    }

    return aggregator.get_alignments();
}

template
std::vector<Alignment> chain_alignments<LocalAlignmentLess>(const IDBGAligner&,
                                                            std::vector<Alignment>&&);

} // namespace align
} // namespace graph
} // namespace mtg
