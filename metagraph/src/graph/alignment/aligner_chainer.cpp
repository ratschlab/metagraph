#include "aligner_chainer.hpp"

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

bool last_is_indel(const Cigar &a) {
    auto it = a.data().rbegin();
    if (it->first == Cigar::CLIPPED)
        ++it;

    return it != a.data().rend() && (it->first == Cigar::INSERTION || it->first == Cigar::DELETION);
}

bool first_is_indel(const Cigar &a) {
    auto it = a.data().begin();
    if (it->first == Cigar::CLIPPED)
        ++it;

    return it != a.data().end() && (it->first == Cigar::INSERTION || it->first == Cigar::DELETION);
}

std::pair<size_t, size_t> call_alignment_chains(const IDBGAligner &aligner,
                                                const std::vector<Alignment> &alignments,
                                                const std::function<void(Alignment&&)> &callback,
                                                std::string_view query,
                                                std::string_view query_rc);

template <class AlignmentCompare>
std::tuple<std::vector<Alignment>, size_t, size_t>
chain_alignments(const IDBGAligner &aligner,
                 std::vector<Alignment>&& alignments,
                 std::string_view query,
                 std::string_view query_rc) {
    const DBGAlignerConfig &config = aligner.get_config();

    if (alignments.size() < 2) {
        DEBUG_LOG("Too few alignments found, nothing to chain.");
        return std::make_tuple(std::move(alignments), 0, 0);
    }

    for (const auto &a : alignments) {
        if (a.label_coordinates.size())
            throw std::runtime_error("Post-chaining alignments with coordinates not supported");
    }

    DBGAlignerConfig no_chain_config { config };
    no_chain_config.post_chain_alignments = false;
    AlignmentAggregator<AlignmentCompare> aggregator(no_chain_config);

    std::vector<Alignment> split_alignments[2];

    DEBUG_LOG("Chaining alignments:");
    for (size_t i = 0; i < alignments.size(); ++i) {
        auto &a = alignments[i];
        assert(a.size());

        // TODO: handle offset and mixed label cases later
        if ((!a.get_clipping() && !a.get_end_clipping())
                || a.size() == 1 || a.label_column_diffs.size()) {
            aggregator.add_alignment(std::move(a));
        } else {
            // if (evalues.size()) {
                // logger->trace("\t{}\t{}", evalues[i], a);
            // } else {
                DEBUG_LOG("\t{}", a);
            // }
            size_t orientation = a.get_orientation();
            split_alignments[orientation].emplace_back(std::move(a));
        }
    }

    if (split_alignments[0].empty() && split_alignments[1].empty())
        return std::make_tuple(aggregator.get_alignments(), 0, 0);

    auto acomp = [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_clipping() + a.get_query_view().size(),
                               a.get_clipping(),
                               b.get_score(),
                               a.get_sequence().size())
            < std::make_tuple(b.get_clipping() + b.get_query_view().size(),
                              b.get_clipping(),
                              a.get_score(),
                              b.get_sequence().size());
    };

    std::sort(split_alignments[0].begin(), split_alignments[0].end(), acomp);
    std::sort(split_alignments[1].begin(), split_alignments[1].end(), acomp);

    auto callback = [&](Alignment&& a) {
        aggregator.add_alignment(std::move(a));
    };

    auto [num_extensions, num_explored_nodes]
        = call_alignment_chains(aligner, split_alignments[0], callback, query, query_rc);
    auto [n_ext, n_exp] = call_alignment_chains(aligner, split_alignments[1], callback, query_rc, query);
    num_extensions += n_ext;
    num_explored_nodes += n_exp;
    return std::make_tuple(aggregator.get_alignments(), num_extensions, num_explored_nodes);
}

std::pair<size_t, size_t> call_alignment_chains(const IDBGAligner &aligner,
                                                const std::vector<Alignment> &alignments,
                                                const std::function<void(Alignment&&)> &callback,
                                                std::string_view query,
                                                std::string_view query_rc) {
    if (alignments.empty())
        return {};

    assert(std::all_of(alignments.begin(), alignments.end(), [](const auto &a) { return a.size(); }));
    assert(std::all_of(alignments.begin(), alignments.end(), [](const auto &a) { return !a.get_offset(); }));

    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    const DBGAlignerConfig &config = aligner.get_config();
    const DeBruijnGraph &graph = aligner.get_graph();
    ssize_t node_overlap = graph.get_k() - 1;

    DEBUG_LOG("Preprocessing alignments");
    std::vector<const char *> min_next_char;
    min_next_char.reserve(alignments.size());
    for (auto jt = alignments.begin() + 1; jt != alignments.end(); ++jt) {
        min_next_char.push_back(jt->get_query_view().data());
    }
    min_next_char.push_back(std::numeric_limits<const char *>::max());
    for (auto it = min_next_char.rbegin() + 1; it != min_next_char.rend(); ++it) {
        *it = std::min(*it, *(it - 1));
    }

    std::vector<const char *> max_prev_char;
    max_prev_char.reserve(alignments.size());
    max_prev_char.push_back(std::numeric_limits<const char *>::min());
    for (auto jt = alignments.begin(); jt != alignments.end() - 1; ++jt) {
        max_prev_char.push_back(jt->get_query_view().data() + jt->get_query_view().size());
    }
    for (auto it = max_prev_char.begin() + 1; it != max_prev_char.end(); ++it) {
        *it = std::max(*it, *(it - 1));
    }

    std::vector<std::vector<std::pair<score_t, size_t>>> prefix_score_size_diffs(alignments.size());
    std::vector<std::vector<std::pair<score_t, size_t>>> suffix_score_size_diffs(alignments.size());
    for (size_t i = 0; i < alignments.size(); ++i) {
        if (i) {
            for (Alignment a = alignments[i]; a.size(); a.trim_query_prefix(1, node_overlap, config)) {
                assert(alignments[i].get_sequence().size() >= a.get_sequence().size());
                score_t score = first_is_indel(a.get_cigar()) ? DBGAlignerConfig::ninf : alignments[i].get_score() - a.get_score();
                prefix_score_size_diffs[i].emplace_back(score, alignments[i].get_sequence().size() - a.get_sequence().size());
                if (a.get_query_view().data() > max_prev_char[i] && !first_is_indel(a.get_cigar()))
                    break;
            }
        } else {
            prefix_score_size_diffs[i].emplace_back(0, 0);
        }

        if (i + 1 < alignments.size()) {
            for (Alignment a = alignments[i]; a.size(); a.trim_query_suffix(1, config)) {
                assert(alignments[i].get_sequence().size() >= a.get_sequence().size());
                score_t score = last_is_indel(a.get_cigar()) ? DBGAlignerConfig::ninf : alignments[i].get_score() - a.get_score();
                suffix_score_size_diffs[i].emplace_back(score, alignments[i].get_sequence().size() - a.get_sequence().size());
                if (a.get_query_view().data() + a.get_query_view().size() < min_next_char[i] && last_is_indel(a.get_cigar()))
                    break;
            }
        } else {
            suffix_score_size_diffs[i].emplace_back(0, 0);
        }
    }

    DEBUG_LOG("Done preprocessing alignments, starting chaining");

    typedef std::tuple<score_t /* best score */,
                       score_t /* extra score */,
                       ssize_t /* gap */,
                       size_t /* prefix trim */,
                       size_t /* last table i */,
                       Alignment::Columns /* last labels */,
                       size_t /* last suffix trim */> TableVal;
    typedef tsl::hopscotch_map<Alignment::Columns, TableVal> Table;
    std::vector<Table> chain_table(alignments.size());

    for (size_t i = 0; i < alignments.size() - 1; ++i) {
        auto find = chain_table[i].find(alignments[i].label_columns);
        if (find == chain_table[i].end()) {
            chain_table[i][alignments[i].label_columns] = TableVal{
                alignments[i].get_score(), 0, 0, 0,
                std::numeric_limits<size_t>::max(),
                0, 0
            };
        } else if (alignments[i].get_score() > std::get<0>(find->second)) {
            find.value() = TableVal{
                alignments[i].get_score(), 0, 0, 0,
                std::numeric_limits<size_t>::max(),
                0, 0
            };
        }
        for (const auto &[cur_columns, val] : chain_table[i]) {
            score_t a_score = std::get<0>(val);
            size_t a_prefix_trim = std::get<3>(val);
            for (size_t j = i + 1; j < alignments.size(); ++j) {
                const char *a_begin = alignments[i].get_query_view().data() + a_prefix_trim;
                const char *a_end = alignments[i].get_query_view().data() + alignments[i].get_query_view().size();
                const char *b_begin = alignments[j].get_query_view().data();
                const char *b_end = b_begin + alignments[j].get_query_view().size();

                if (b_end <= a_end)
                    continue;

                size_t b_prefix_trim = 0;
                if (b_begin < a_begin) {
                    b_prefix_trim = a_begin - b_begin;
                    b_begin += b_prefix_trim;
                }

                if (b_prefix_trim >= prefix_score_size_diffs[j].size())
                    continue;

                auto [prefix_score_diff, prefix_size_diff] = prefix_score_size_diffs[j][b_prefix_trim];
                while (b_prefix_trim < prefix_score_size_diffs[j].size() && prefix_score_diff == DBGAlignerConfig::ninf) {
                    std::tie(prefix_score_diff, prefix_size_diff) = prefix_score_size_diffs[j][++b_prefix_trim];
                    ++b_begin;
                }

                if (b_prefix_trim == prefix_score_size_diffs[j].size() || b_begin == b_end)
                    continue;

                size_t b_offset = std::min(prefix_size_diff, graph.get_k() - 1);
                size_t b_size = alignments[j].get_sequence().size() - prefix_size_diff - graph.get_k() + 1 + b_offset;

                ssize_t gap = b_begin - a_end;

                score_t b_score = DBGAlignerConfig::ninf;
                char c = alignments[j].get_sequence()[prefix_size_diff];
                size_t a_suffix_trim = 0;
                if (gap >= 0) {
                    // no overlap
                    if (b_offset >= b_size)
                        continue;

                    score_t gap_score = config.gap_opening_penalty;
                    if (gap > 0)
                        gap_score += config.gap_opening_penalty + (gap - 1) * config.gap_extension_penalty;

                    b_score = a_score + alignments[j].get_score() - prefix_score_diff + gap_score;
                } else {
                    // overlap
                    gap = std::max({
                        gap,
                        static_cast<ssize_t>(b_prefix_trim)
                            - static_cast<ssize_t>(prefix_score_size_diffs[j].size()) + 1
                    });

                    b_score = DBGAlignerConfig::ninf;
                    ssize_t best_mismatch = 0;
                    ssize_t best_gap = gap;
                    for (ssize_t cur_gap = 0; cur_gap >= gap; --cur_gap) {
                        size_t a_suffix_trim = cur_gap - gap;

                        if (a_suffix_trim >= suffix_score_size_diffs[i].size())
                            continue;

                        if (a_suffix_trim + a_prefix_trim >= alignments[i].get_query_view().size())
                            continue;

                        assert(alignments[i].get_query_view().data() + alignments[i].get_query_view().size() - (cur_gap - gap)
                            == alignments[j].get_query_view().data() + (b_prefix_trim - cur_gap));
                        auto [b_prefix_score_diff, b_seq_prefix_trim] = prefix_score_size_diffs[j][b_prefix_trim - cur_gap];
                        auto [a_suffix_score_diff, a_seq_suffix_trim] = suffix_score_size_diffs[i][a_suffix_trim];
                        if (std::min(a_suffix_score_diff, b_prefix_score_diff) == DBGAlignerConfig::ninf)
                            continue;

                        assert(a_prefix_trim < prefix_score_size_diffs[i].size());
                        size_t a_seq_prefix_trim = prefix_score_size_diffs[i][a_prefix_trim].second;

                        if (alignments[i].get_sequence().size() - a_seq_suffix_trim - a_seq_prefix_trim < graph.get_k()
                                || alignments[j].get_sequence().size() - b_seq_prefix_trim < graph.get_k()) {
                            continue;
                        }

                        score_t cur_score = a_score - a_suffix_score_diff + alignments[j].get_score() - b_prefix_score_diff;

                        size_t a_seq_extra_trim = std::min(b_seq_prefix_trim,
                                                           alignments[i].get_sequence().size() - a_seq_suffix_trim);

                        size_t cur_b_offset = std::min(b_seq_prefix_trim, graph.get_k() - 1);
                        size_t cur_b_size = alignments[j].get_sequence().size() - b_seq_prefix_trim - graph.get_k() + 1 + cur_b_offset;
                        ssize_t cur_mismatch = 0;
                        if (!b_seq_prefix_trim || !a_seq_extra_trim) {
                            if (cur_b_offset >= cur_b_size)
                                continue;

                            cur_score += config.gap_opening_penalty;
                        } else {
                            std::string_view b_prefix(alignments[j].get_sequence().begin(),
                                                      b_seq_prefix_trim);

                            std::string_view a_suffix(alignments[i].get_sequence().end() - a_seq_suffix_trim - a_seq_extra_trim,
                                                      a_seq_extra_trim);

                            if (b_prefix.size() >= graph.get_k()) {
                                b_prefix.remove_prefix(b_prefix.size() - node_overlap);
                                a_suffix.remove_prefix(a_suffix.size() - node_overlap);
                            }

                            // TODO: optimize this later if needed
                            auto [a_mm, b_mm] = std::mismatch(a_suffix.rbegin(), a_suffix.rend(),
                                                              b_prefix.rbegin(), b_prefix.rend());

                            cur_mismatch = -(b_mm - b_prefix.rbegin());

                            if (cur_mismatch == 0) {
                                if (cur_b_offset >= cur_b_size)
                                    continue;

                                cur_score += config.gap_opening_penalty;
                            } else if (cur_mismatch > -node_overlap && cur_b_size <= cur_b_offset + cur_mismatch) {
                                continue;
                            }

                        }

                        if (cur_score > b_score) {
                            b_score = cur_score;
                            best_gap = cur_gap;
                            best_mismatch = cur_mismatch;
                        }
                    }

                    if (b_score == DBGAlignerConfig::ninf)
                        continue;

                    b_prefix_trim -= best_gap;
                    b_begin -= best_gap;

                    auto [b_prefix_score_diff, b_seq_prefix_trim] = prefix_score_size_diffs[j][b_prefix_trim];
                    b_offset = std::min(b_seq_prefix_trim, graph.get_k() - 1);
                    b_size = alignments[j].get_sequence().size() - b_seq_prefix_trim - graph.get_k() + 1 + b_offset;

                    c = alignments[j].get_sequence()[b_seq_prefix_trim];

                    a_suffix_trim = best_gap - gap;
                    gap = best_mismatch;
                    assert(gap <= 0);
                    assert(gap == 0 || gap <= -node_overlap || b_size > b_offset + gap);
                }

#ifndef NDEBUG
                Alignment a = alignments[i];
                if (a_suffix_trim) {
                    a.trim_query_suffix(a_suffix_trim, config);
                    assert(a.size());
                    assert(a.is_valid(graph, &config));
                }

                if (a_prefix_trim) {
                    a.trim_query_prefix(a_prefix_trim, node_overlap, config);
                    assert(a.size());
                    assert(a.is_valid(graph, &config));
                }

                Alignment b = alignments[j];
                if (b_prefix_trim) {
                    b.trim_query_prefix(b_prefix_trim, node_overlap, config);
                    assert(b.size());
                    assert(b.is_valid(graph, &config));
                }

                assert(b.get_offset() == b_offset);
                assert(b.size() == b_size);

                if (gap > -node_overlap) {
                    b.insert_gap_prefix(gap, node_overlap, config);
                    assert(b.size());
                } else {
                    b.trim_clipping();
                    assert(b.size());
                }

                assert(a.size());
                a.trim_end_clipping();
                a.append(std::move(b), a.label_columns != alignments[j].label_columns);
                assert(a.size());
                assert(a.is_valid(graph, &config));
#endif

                ILabeledAligner::LabelChangeScores label_change_scores;
                if (!labeled_aligner) {
                    label_change_scores.emplace_back(0, 0);
                } else {
                    label_change_scores = labeled_aligner->get_label_change_scores(
                        cur_columns, alignments[j].label_columns,
                        labeled_aligner->get_annotation_buffer().get_hll_wrapper()
                    );
                }

                auto &b_tab = chain_table[j];
                score_t lambda = config.score_matrix[c][c];
                for (auto &[cols, lc_score] : label_change_scores) {
                    lc_score *= lambda;
                    assert(lc_score <= 0);
                    score_t next_score = b_score + lc_score;
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

    DEBUG_LOG("Done chaining, backtracking");

    std::vector<std::tuple<score_t, size_t, Alignment::Columns>> indices;
    sdsl::bit_vector used(chain_table.size(), false);
    for (size_t i = 0; i < chain_table.size(); ++i) {
        for (const auto &[cols, vals] : chain_table[i]) {
            indices.emplace_back(std::get<0>(vals), cols, i);
        }
    }
    std::sort(indices.begin(), indices.end());

    tsl::hopscotch_map<size_t, Alignment> front_ext_map;
    size_t num_extensions = 0;
    size_t num_explored_nodes = 0;
    for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
        auto [start_score, cols, i] = *it;
        if (used[i])
            continue;

        used[i] = true;
        Alignment cur = alignments[i];

        auto [score, extra_score, gap, prefix_trim,
              last_i, last_cols, last_suffix_trim] = chain_table[i][cols];

        if (prefix_trim) {
            cur.trim_query_prefix(prefix_trim, node_overlap, config);
            assert(cur.size());
            assert(cur.is_valid(graph, &config));
        }

        if (cur.get_end_clipping()) {
            DEBUG_LOG("Extending chain end");
            auto extender = aligner.make_extender(query);
            auto [left, next] = split_seed(graph, config, cur);
            auto extensions = extender->get_extensions(next,
                config.ninf,                         // min_path_score
                true,                                // force_fixed_seed
                0,                                   // target_length
                DeBruijnGraph::npos,                 // target_node
                false
            );

            if (extensions.size() && extensions[0].get_end_clipping() < cur.get_end_clipping()) {
                assert(extensions[0].get_nodes().front() == next.get_nodes().front());
                DEBUG_LOG("left:\n\t{}\next\n\t{}", left, extensions[0]);
                left.splice(std::move(extensions[0]));
                assert(cur.get_offset() == left.get_offset());
                std::swap(left, cur);
                assert(cur.is_valid(graph, &config));
                DEBUG_LOG("Extended successfully:\n\t{}", cur);
            }
            ++num_extensions;
            num_explored_nodes += extender->num_explored_nodes();
        }

        while (last_i != std::numeric_limits<size_t>::max()) {
            Alignment prev = alignments[last_i];
            used[last_i] = true;

            if (std::get<4>(chain_table[last_i][last_cols]) == std::numeric_limits<size_t>::max()) {
                assert(!std::get<3>(chain_table[last_i][last_cols]));
                if (prev.get_clipping()) {
                    DEBUG_LOG("Extending chain front");
                    auto find = front_ext_map.find(i);
                    if (find != front_ext_map.end()) {
                        if (find->second.size()) {
                            // start_score += find->second.get_score() - prev.get_score();
                            prev = find->second;
                        }
                    } else {
                        // extend backwards
                        RCDBG rc_dbg(std::shared_ptr<const DeBruijnGraph>(
                                        std::shared_ptr<const DeBruijnGraph>(), &graph));
                        const DeBruijnGraph &rc_graph = graph.get_mode() != DeBruijnGraph::CANONICAL
                            ? rc_dbg : graph;

                        auto rev = prev;
                        rev.reverse_complement(rc_graph, query_rc);
                        if (rev.size() && rev.get_nodes().back()) {
                            assert(rev.get_end_clipping());
                            auto extender_rc = aligner.make_extender(query_rc);
                            extender_rc->set_graph(rc_graph);
                            auto extensions = extender_rc->get_extensions(rev, config.ninf, true);
                            ++num_extensions;
                            num_explored_nodes += extender_rc->num_explored_nodes();
                            if (extensions.size() && extensions[0].get_score() <= prev.get_score())
                                extensions.clear();

                            if (extensions.size()) {
                                extensions[0].reverse_complement(rc_graph, query);
                                if (extensions[0].empty())
                                    extensions.clear();
                            }

                            if (extensions.size()) {
                                front_ext_map[i] = extensions[0];
                                // start_score += extensions[0].get_score() - prev.get_score();
                                assert(extensions[0].get_end_clipping() == prev.get_end_clipping());
                                DEBUG_LOG("\tsuccess");
                                std::swap(prev, extensions[0]);
                            } else {
                                front_ext_map[i] = Alignment();
                            }
                        }
                    }
                }
            }
            if (last_suffix_trim) {
                // score_t cur_prev_score = prev.get_score();
                prev.trim_query_suffix(last_suffix_trim, config);
                if (last_is_indel(prev.get_cigar())) {
                    DEBUG_LOG("Invalid front extension");
                    // start_score -= cur_prev_score - alignments[last_i].get_score();
                    prev = alignments[last_i];
                    prev.trim_query_suffix(last_suffix_trim, config);
                }
                assert(prev.size());
                assert(prev.is_valid(graph, &config));
            }

            if (gap > -node_overlap) {
                cur.insert_gap_prefix(gap, node_overlap, config);
            } else {
                cur.trim_clipping();
            }

            i = last_i;
            cols = last_cols;
            std::tie(score, extra_score, gap, prefix_trim,
                     last_i, last_cols, last_suffix_trim) = chain_table[i][cols];

            if (prefix_trim) {
                prev.trim_query_prefix(prefix_trim, node_overlap, config);
                assert(prev.size());
                assert(prev.is_valid(graph, &config));
            }

            prev.trim_end_clipping();
            prev.append(std::move(cur), extra_score);
            assert(prev.size());
            std::swap(cur, prev);
            assert(cur.is_valid(graph, &config));
            DEBUG_LOG("Partial alignment: {}\n\t{}", extra_score, cur);
        }

        // assert(cur.get_score() == start_score);

        assert(cur.get_clipping() + cur.get_query_view().size() + cur.get_end_clipping()
            == query.size());
        assert(cur.is_valid(graph, &config));
        DEBUG_LOG("Full alignment:\n\t{}", cur);
        callback(std::move(cur));
    }

    return std::make_pair(num_extensions, num_explored_nodes);
}

template
std::tuple<std::vector<Alignment>, size_t, size_t>
chain_alignments<LocalAlignmentLess>(const IDBGAligner&,
                                     std::vector<Alignment>&&,
                                     std::string_view,
                                     std::string_view);

} // namespace align
} // namespace graph
} // namespace mtg
