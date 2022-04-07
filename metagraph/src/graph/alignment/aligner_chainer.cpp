#include "aligner_chainer.hpp"

#include "aligner_seeder_methods.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"

#include "common/aligned_vector.hpp"
#include "common/utils/simd_utils.hpp"

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
                                   [](const auto &a) {
                                       return a.empty() || a.label_columns.empty();
                                   }),
                    fwd_seeds.end());
    bwd_seeds.erase(std::remove_if(bwd_seeds.begin(), bwd_seeds.end(),
                                   [](const auto &a) {
                                       return a.empty() || a.label_columns.empty();
                                   }),
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

    std::sort(starts.begin(), starts.end(), std::greater<decltype(starts)::value_type>());

    Chain last_chain;
    score_t last_chain_score = std::numeric_limits<score_t>::min();

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
            chain_seeds.back().first.label_columns.clear();
            chain_seeds.back().first.label_coordinates.clear();
            if (has_labels) {
                chain_seeds.back().first.label_columns.assign(1, label);
                if (aligner.has_coordinates()) {
                    chain_seeds.back().first.label_coordinates.resize(1);
                    chain_seeds.back().first.label_coordinates[0].assign(1, coord);
                }
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

        Chain chain;
        chain.reserve(chain_seeds.size());
        std::transform(chain_seeds.begin(), chain_seeds.end(), std::back_inserter(chain),
                       [&](const auto &c) {
                           return std::make_pair(Alignment(c.first, config), c.second);
                       });

        if (last_chain.empty()) {
            std::swap(last_chain, chain);
            last_chain_score = chain_score;
            continue;
        }

        if (last_chain != chain) {
            callback(std::move(last_chain), last_chain_score);
            std::swap(last_chain, chain);
            last_chain_score = chain_score;
            continue;
        }

        if (chain[0].first.label_columns.empty())
            continue;

        // if this chain has the same seeds as the last one, merge their coordinate sets
        for (size_t i = 0; i < chain.size(); ++i) {
            Alignment::Columns columns;
            if (chain[i].first.label_coordinates.size()) {
                assert(last_chain[i].first.label_columns.size()
                        == last_chain[i].first.label_coordinates.size());
                assert(chain[i].first.label_columns.size()
                        == chain[i].first.label_coordinates.size());
                Alignment::CoordinateSet coord_union;
                auto add_col_coords = [&](auto col, const auto &coords) {
                    columns.push_back(col);
                    coord_union.push_back(coords);
                };
                utils::match_indexed_values(
                    last_chain[i].first.label_columns.begin(),
                    last_chain[i].first.label_columns.end(),
                    last_chain[i].first.label_coordinates.begin(),
                    chain[i].first.label_columns.begin(),
                    chain[i].first.label_columns.end(),
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
                assert(chain[i].first.label_columns.size());
                std::set_union(last_chain[i].first.label_columns.begin(),
                               last_chain[i].first.label_columns.end(),
                               chain[i].first.label_columns.begin(),
                               chain[i].first.label_columns.end(),
                               std::back_inserter(columns));
            }
            std::swap(last_chain[i].first.label_columns, columns);
        }
    }

    if (last_chain.size())
        callback(std::move(last_chain), last_chain_score);

    return std::make_pair(num_seeds, num_nodes);
}

inline __m256i avx2_log2_epi32(__m256i v) {
    //avx2_lzcnt_epi32
    // prevent value from being rounded up to the next power of two
    v = _mm256_andnot_si256(_mm256_srli_epi32(v, 8), v); // keep 8 MSB

    v = _mm256_castps_si256(_mm256_cvtepi32_ps(v)); // convert an integer to float
    v = _mm256_srli_epi32(v, 23); // shift down the exponent
    v = _mm256_subs_epu16(_mm256_set1_epi32(158), v); // undo bias
    v = _mm256_min_epi16(v, _mm256_set1_epi32(32)); // clamp at 32

    // log2
    v = _mm256_sub_epi32(_mm256_set1_epi32(31), v);

    return v;
}

std::tuple<ChainDPTable, AlignedVector<int32_t>, size_t, size_t>
chain_seeds(const IDBGAligner &aligner,
            const DBGAlignerConfig &config,
            std::string_view query,
            std::vector<Seed> &seeds) {
    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    if (seeds.empty() || !labeled_aligner || !aligner.has_coordinates())
        return {};

    AnnotationBuffer &buffer = labeled_aligner->get_annotation_buffer();

    size_t num_nodes = 0;

    ssize_t query_size = query.size();

    ChainDPTable dp_table;
    dp_table.reserve(seeds.size());

    logger->trace("Sorting seeds");
    std::sort(seeds.begin(), seeds.end(), [](const auto &a, const auto &b) {
        return a.get_clipping() > b.get_clipping();
    });

    tsl::hopscotch_map<Alignment::Column, size_t> label_sizes;

    for (size_t i = 0; i < seeds.size(); ++i) {
        if (!(i % config.row_batch_size)) {
            size_t end = std::min(i + config.row_batch_size, seeds.size());
            std::vector<node_index> nodes;
            nodes.reserve(end - i);
            for (size_t j = i; j < end; ++j) {
                nodes.emplace_back(seeds[j].get_nodes()[0]);
            }
            buffer.prefetch_coords(nodes);
        }

        auto fetch_coords = buffer.get_coords(seeds[i].get_nodes()[0], false,
                                              &seeds[i].label_columns);
        assert(fetch_coords);

        for (size_t j = 0; j < fetch_coords->size(); ++j) {
            Alignment::Column c = seeds[i].label_columns[j];
            size_t num_seeds = std::min((*fetch_coords)[j].size(),
                                        config.max_num_seeds_per_locus);
            auto rbegin = (*fetch_coords)[j].rbegin();
            auto rend = rbegin + num_seeds;
            std::for_each(rbegin, rend, [&](ssize_t coord) {
                ++label_sizes[c];
                dp_table.emplace_back(c, coord + seeds[i].get_offset(), seeds[i].get_clipping(),
                                      seeds[i].get_clipping() + seeds[i].get_query_view().size(),
                                      seeds[i].get_query_view().size(), i);
            });
        }
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

    size_t bandwidth = 201;

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
                    auto mg_log2_avx2 = [](__m256 dd_v_f) -> __m256 // NB: this doesn't work when x<2
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

                    // float log_penalty = log2(coord_diff) * 0.5;
                    __m256 log_penalty_v = mg_log2_avx2(coord_diff_f);
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
                        float log_penalty = log2(coord_diff) * 0.5;
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
void construct_alignment_chain(size_t node_overlap,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> *best_score,
                               const std::function<void(Alignment&&)> &callback);

template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(std::vector<Alignment>&& alignments,
                                        std::string_view query,
                                        std::string_view rc_query,
                                        const DBGAlignerConfig &config,
                                        size_t node_overlap) {
    if (alignments.size() < 2 || !config.post_chain_alignments)
        return std::move(alignments);

    for (const auto &a : alignments) {
        if (a.label_coordinates.size())
            throw std::runtime_error("Post-chaining alignments with coordinates not supported");
    }

    DBGAlignerConfig no_chain_config { config };
    no_chain_config.post_chain_alignments = false;
    AlignmentAggregator<AlignmentCompare> aggregator(no_chain_config);

    alignments.erase(std::remove_if(alignments.begin(), alignments.end(), [&](Alignment &a) {
        if (!a.get_clipping() && !a.get_end_clipping()) {
            aggregator.add_alignment(std::move(a));
            return true;
        }

        return false;
    }), alignments.end());

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

    DEBUG_LOG("Chaining alignments:\n{}", fmt::join(alignments, "\t\n"));


    auto run = [&](std::string_view this_query, auto begin, auto end) {
        std::vector<score_t> best_score(this_query.size() + 1, 0);
        for (auto it = begin; it != end; ++it) {
            size_t end_pos = it->get_query_view().data() + it->get_query_view().size()
                                - this_query.data();
            if (it->get_score() > best_score[end_pos]) {
                best_score[end_pos] = it->get_score();
                construct_alignment_chain<AlignmentCompare>(
                    node_overlap, config, this_query, Alignment(*it), it + 1, end, &best_score,
                    [&](Alignment&& chain) { aggregator.add_alignment(std::move(chain)); }
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
void construct_alignment_chain(size_t node_overlap,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> *best_score,
                               const std::function<void(Alignment&&)> &callback) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.get_query_view().data();
    const char *chain_end = chain.get_query_view().data() + chain.get_query_view().size();
    if (begin == end || chain_end == query.data() + query.size()) {
        callback(std::move(chain));
        return;
    }

    score_t score = chain.get_score();

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        // TODO: handle this case later
        if (it->get_offset())
            continue;

        const char *next_begin = it->get_query_view().data();
        const char *next_end = it->get_query_view().data() + it->get_query_view().size();

        assert(chain_begin - chain.get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.get_orientation());

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        if (chain.label_columns.size()
                && !utils::share_element(it->label_columns.begin(),
                                         it->label_columns.end(),
                                         chain.label_columns.begin(),
                                         chain.label_columns.end())) {
            continue;
        }

        Alignment aln = *it;

        if (next_begin >= chain_end) {
            // no overlap
            aln.insert_gap_prefix(next_begin - chain_end, node_overlap, config);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t overlap = std::min(
                static_cast<size_t>((chain.get_cigar().data().end() - 2)->second),
                aln.trim_query_prefix(chain_end - it->get_query_view().data(),
                                      node_overlap, config)
            );

            if (aln.empty() || aln.get_sequence().size() <= node_overlap
                    || (aln.get_cigar().data().begin()
                            + static_cast<bool>(aln.get_clipping()))->first != Cigar::MATCH) {
                continue;
            }

            assert(aln.get_query_view().data()
                    == chain.get_query_view().data() + chain.get_query_view().size());

            if (overlap < node_overlap) {
                aln.insert_gap_prefix(-overlap, node_overlap, config);
            } else {
                aln.trim_clipping();
            }
        }

        assert(!aln.empty());

        score_t next_score = score + aln.get_score();
        if (next_score <= (*best_score)[next_end - query.data()])
            continue;

        (*best_score)[next_end - query.data()] = next_score;
        // use append instead of splice because any clipping in aln represents
        // internally clipped characters
        Alignment next_chain = chain;
        next_chain.trim_end_clipping();
        bool changed = next_chain.append(std::move(aln));
        if (next_chain.size()) {
            assert(next_chain.get_score() == next_score);
            construct_alignment_chain<AlignmentCompare>(
                    node_overlap, config, query, std::move(next_chain),
                    it + 1, end, best_score, callback);
            called |= changed;
        }
    }

    if (!called)
        callback(std::move(chain));
}

template
std::vector<Alignment> chain_alignments<LocalAlignmentLess>(std::vector<Alignment>&&,
                                                            std::string_view,
                                                            std::string_view,
                                                            const DBGAlignerConfig&,
                                                            size_t);

} // namespace align
} // namespace graph
} // namespace mtg
