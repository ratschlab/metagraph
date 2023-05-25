#include "aligner_chainer.hpp"

#include <x86/svml.h>

#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"
#include "chainer.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"
#include "graph/graph_extensions/path_index.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg {
namespace graph {
namespace align {

using common::logger;

typedef DeBruijnGraph::node_index node_index;

constexpr uint32_t nid = std::numeric_limits<uint32_t>::max();

struct TableElem {
    Alignment::Column label;
    int64_t coordinate;
    int32_t seed_clipping;
    int32_t seed_end;
    score_t chain_score;
    uint32_t current_seed_index;

    TableElem(Alignment::Column c, int64_t coordinate, int32_t seed_clipping,
              int32_t seed_end, score_t chain_score, uint32_t current_seed_index)
          : label(c), coordinate(coordinate), seed_clipping(seed_clipping),
            seed_end(seed_end), chain_score(chain_score), current_seed_index(current_seed_index) {}
} SIMDE_ALIGN_TO_32;
static_assert(sizeof(TableElem) == 32);

inline constexpr bool operator>(const TableElem &a, const TableElem &b) {
    return std::tie(a.label, a.coordinate, a.seed_clipping, a.seed_end)
        > std::tie(b.label, b.coordinate, b.seed_clipping, b.seed_end);
}

typedef AlignedVector<TableElem> ChainDPTable;

std::tuple<ChainDPTable, AlignedVector<int32_t>, size_t, size_t>
chain_seeds(const DBGAlignerConfig &config,
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
                              const std::function<bool(Alignment::Column)> &skip_column,
                              const std::function<bool()> &terminate) {
    fwd_seeds.erase(std::remove_if(fwd_seeds.begin(), fwd_seeds.end(),
                                   [](const auto &a) { return a.empty() || !a.label_columns; }),
                    fwd_seeds.end());
    bwd_seeds.erase(std::remove_if(bwd_seeds.begin(), bwd_seeds.end(),
                                   [](const auto &a) { return a.empty() || !a.label_columns; }),
                    bwd_seeds.end());

    if (terminate() || (fwd_seeds.empty() && bwd_seeds.empty()))
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
        = chain_seeds(config, forward, both_seeds[0]);

    logger->trace("Chaining reverse complement seeds");
    size_t num_seeds_bwd;
    size_t num_nodes_bwd;
    std::tie(dp_tables[1], seed_backtraces[1], num_seeds_bwd, num_nodes_bwd)
        = chain_seeds(config, reverse, both_seeds[1]);

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

    bool coverage_too_low = false;

    auto flush_chains = [&]() {
        if (coverage_too_low)
            return;

        assert(chains.size());
        auto it = chains.begin();
        Chain last_chain = *it;
        for (++it; it != chains.end(); ++it) {
            const Chain &chain = *it;
            if (chain != last_chain) {
                double exact_match_fraction
                    = static_cast<double>(get_num_char_matches_in_seeds(last_chain.begin(),
                                                                        last_chain.end()))
                        / forward.size();

                if (exact_match_fraction < config.min_exact_match) {
                    coverage_too_low = true;
                    return;
                }

                callback(std::move(last_chain), last_chain_score);
                last_chain = *it;
                continue;
            }

            // if this chain has the same seeds as the last one, merge their coordinate sets
            for (size_t i = 0; i < chain.size(); ++i) {
                Vector<Alignment::Column> columns;
                const auto &last_columns = last_chain[i].first.get_columns();
                const auto &cur_columns = chain[i].first.get_columns();
                if (chain[i].first.label_coordinates.size()) {
                    assert(last_chain[i].first.label_columns
                            == last_chain[i].first.label_coordinates.size());
                    assert(chain[i].first.label_columns
                            == chain[i].first.label_coordinates.size());
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
                    assert(chain[i].first.label_columns);
                    std::set_union(last_columns.begin(), last_columns.end(),
                                   cur_columns.begin(), cur_columns.end(),
                                   std::back_inserter(columns));
                }
                last_chain[i].first.set_columns(std::move(columns));
            }
        }

        double exact_match_fraction
            = static_cast<double>(get_num_char_matches_in_seeds(last_chain.begin(),
                                                                last_chain.end()))
                / forward.size();

        if (exact_match_fraction < config.min_exact_match) {
            coverage_too_low = true;
            return;
        }

        callback(std::move(last_chain), last_chain_score);

        chains.clear();
    };

    for (const auto &[chain_score, j, neg_i] : starts) {
        if (coverage_too_low || terminate())
            break;

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
            chain_seeds.back().first.label_encoder = seeds[seed_i].label_encoder;
            if (has_labels) {
                chain_seeds.back().first.set_columns(Vector<Alignment::Column>{ label });
                chain_seeds.back().first.label_coordinates.resize(1);
                chain_seeds.back().first.label_coordinates[0].assign(1, coord);
            }
            i = seed_backtrace[i];
        }

        if (chain_seeds.empty())
            continue;

        // clean chain by merging overlapping seeds
        for (size_t i = chain_seeds.size() - 1; i > 0; --i) {
            auto &cur_seed = chain_seeds[i].first;
            auto &prev_seed = chain_seeds[i - 1].first;

            assert(cur_seed.size());
            assert(prev_seed.size());
            assert(prev_seed.get_clipping() <= cur_seed.get_clipping());
            assert(prev_seed.get_end_clipping() >= cur_seed.get_end_clipping());

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
chain_seeds(const DBGAlignerConfig &config,
            std::string_view query,
            std::vector<Seed> &seeds) {
    if (seeds.empty())
        return {};

    if (std::any_of(seeds.begin(), seeds.end(),
                    [](const auto &a) { return a.label_coordinates.empty(); })) {
        throw std::runtime_error("Chaining only supported for seeds with coordinates");
    }

    size_t num_nodes = 0;

    ssize_t query_size = query.size();

    ChainDPTable dp_table;
    dp_table.reserve(seeds.size());
    std::reverse(seeds.begin(), seeds.end());

    tsl::hopscotch_map<Alignment::Column, size_t> label_sizes;

    for (size_t i = 0; i < seeds.size(); ++i) {
        const auto &columns = seeds[i].get_columns();
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

            const simde__m256i coord_cutoff_v = simde_mm256_set1_epi64x(coord_cutoff);
            const simde__m256i prev_coord_v = simde_mm256_set1_epi64x(prev_coord);
            const simde__m256i prev_clipping_v = simde_mm256_set1_epi32(prev_clipping);
            const simde__m256i query_size_v = simde_mm256_set1_epi32(query_size);
            const simde__m256i prev_score_v = simde_mm256_set1_epi32(prev_score);
            const simde__m256i it_end_v = simde_mm256_set1_epi32(it_end - 1);
            const simde__m256i i_v = simde_mm256_set1_epi32(i);
            auto epi64_to_epi32 = [](simde__m256i v) {
                return simde_mm256_castsi256_si128(simde_mm256_permute4x64_epi64(simde_mm256_shuffle_epi32(v, 8), 8));
            };

            simde__m256i j_v = simde_mm256_add_epi32(simde_mm256_set1_epi32(i + 1), simde_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0));
            for (size_t j = i + 1; true; j += 8) {
                // if (coord_cutoff > coord || j > it_end - 1)
                //     break;
                simde__m256i coord_1_v = simde_mm256_i32gather_epi64(&dp_table[j].coordinate, simde_mm_set_epi32(12, 8, 4, 0), 8);
                simde__m256i coord_2_v = simde_mm256_i32gather_epi64(&dp_table[j + 4].coordinate, simde_mm_set_epi32(12, 8, 4, 0), 8);
                simde__m256i coord_1_mask = simde_mm256_cmpgt_epi64(coord_cutoff_v, coord_1_v);
                simde__m256i coord_2_mask = simde_mm256_cmpgt_epi64(coord_cutoff_v, coord_2_v);
                simde__m256i coord_neg_mask = simde_mm256_blend_epi32(coord_1_mask, coord_2_mask, 0b10101010);
                simde__m256i j_neg_mask = simde_mm256_cmpgt_epi32(j_v, it_end_v);
                coord_neg_mask = simde_mm256_or_si256(j_neg_mask, coord_neg_mask);

                // int32_t dist = prev_clipping - clipping;
                simde__m256i clipping_v = simde_mm256_i32gather_epi32(&dp_table[j].seed_clipping, simde_mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0), 4);
                simde__m256i dist_v = simde_mm256_sub_epi32(prev_clipping_v, clipping_v);

                // int32_t coord_dist = prev_coord - coord;
                // a[0:32],b[0:32],a[64:96],b[64:96],a[128:160],b[128:160],a[192:224],b[192:224]
                simde__m128i coord_dist_1_v = epi64_to_epi32(simde_mm256_sub_epi64(prev_coord_v, coord_1_v));
                simde__m128i coord_dist_2_v = epi64_to_epi32(simde_mm256_sub_epi64(prev_coord_v, coord_2_v));
                simde__m256i coord_dist_v = simde_mm256_set_m128i(coord_dist_2_v, coord_dist_1_v);

                simde__m256i dist_mask = simde_mm256_cmpgt_epi32(dist_v, simde_mm256_setzero_si256());
                simde__m256i dmax = simde_mm256_max_epi32(dist_v, coord_dist_v);
                simde__m256i dmax_mask = simde_mm256_cmpgt_epi32(query_size_v, dmax);
                dist_mask = simde_mm256_and_si256(dist_mask, dmax_mask);
                dist_mask = simde_mm256_andnot_si256(coord_neg_mask, dist_mask);

                // if (dist > 0 && std::max(dist, coord_dist) < query_size) {
                if (simde_mm256_movemask_epi8(dist_mask)) {
                    // score_t match = std::min({ dist, coord_dist, end - clipping });
                    simde__m256i match_v = simde_mm256_min_epi32(dist_v, coord_dist_v);
                    simde__m256i end_v = simde_mm256_i32gather_epi32(&dp_table[j].seed_end, simde_mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0), 4);
                    simde__m256i length_v = simde_mm256_sub_epi32(end_v, clipping_v);
                    match_v = simde_mm256_min_epi32(match_v, length_v);

                    // score_t cur_score = prev_score + match;
                    simde__m256i cur_score_v = simde_mm256_add_epi32(prev_score_v, match_v);

                    // float coord_diff = std::abs(coord_dist - dist);
                    simde__m256i coord_diff = simde_mm256_sub_epi32(coord_dist_v, dist_v);
                    coord_diff = simde_mm256_abs_epi32(coord_diff);
                    simde__m256 coord_diff_f = simde_mm256_cvtepi32_ps(coord_diff);

                    // float linear_penalty = coord_diff * sl;
                    simde__m256 linear_penalty_v = simde_mm256_mul_ps(coord_diff_f, simde_mm256_set1_ps(sl));

                    // float log_penalty = log2(coord_diff + 1) * 0.5;
                    simde__m256 log_penalty_v = simde_mm256_log2_ps(simde_mm256_cvtepi32_ps(simde_mm256_add_epi32(coord_diff, simde_mm256_set1_epi32(1))));
                    log_penalty_v = simde_mm256_mul_ps(log_penalty_v, simde_mm256_set1_ps(0.5));

                    // cur_score -= linear_penalty + log_penalty;
                    simde__m256 gap_penalty_f = simde_mm256_add_ps(linear_penalty_v, log_penalty_v);
                    simde__m256i gap_penalty_v = simde_mm256_cvtps_epi32(gap_penalty_f);
                    simde__m256i dist_cutoff_mask = simde_mm256_cmpgt_epi32(coord_diff, simde_mm256_setzero_si256());
                    gap_penalty_v = simde_mm256_blendv_epi8(simde_mm256_setzero_si256(), gap_penalty_v, dist_cutoff_mask);
                    cur_score_v = simde_mm256_sub_epi32(cur_score_v, gap_penalty_v);

                    // if (cur_score >= score) {
                    //     score = cur_score;
                    //     backtrace[j] = i;
                    // }
                    simde__m256i old_scores_v = simde_mm256_i32gather_epi32(&dp_table[j].chain_score, simde_mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0), 4);
                    simde__m256i score_neg_mask = simde_mm256_cmpgt_epi32(old_scores_v, cur_score_v);
                    simde__m256i mask = simde_mm256_andnot_si256(score_neg_mask, dist_mask);

                    cur_score_v = simde_mm256_blendv_epi8(old_scores_v, cur_score_v, mask);
                    simde_mm256_maskstore_epi32(&backtrace[j], mask, i_v);

                    // TODO: simde_mm256_i32scatter_epi32 not implemented yet
                    score_t cur_scores[8] SIMDE_ALIGN_TO_32;
                    simde_mm256_store_si256((simde__m256i*)cur_scores, cur_score_v);
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

                if (simde_mm256_movemask_epi8(coord_neg_mask))
                    break;

                j_v = simde_mm256_add_epi32(j_v, simde_mm256_set1_epi32(8));
            }
#if 0
            // reference implementation
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

std::tuple<size_t, size_t, size_t>
chain_and_filter_seeds(const IDBGAligner &aligner,
                       std::shared_ptr<ISeeder> &seeder,
                       SeedFilteringExtender&& extender,
                       SeedFilteringExtender&& bwd_extender,
                       bool allow_label_change) {
    allow_label_change = false;
    size_t query_size = extender.get_query().size();
    const DeBruijnGraph &graph_ = aligner.get_graph();
    const DBGAlignerConfig &config_ = aligner.get_config();
    auto column_path_index = graph_.get_extension_threadsafe<ColumnPathIndex>();
    if (!column_path_index)
        return {};

    const auto &label_encoder = column_path_index->get_anno_graph().get_annotator().get_label_encoder();
    const auto &in_anchors = seeder->get_seeds();
    if (in_anchors.empty())
        return {};

    assert(is_sorted(in_anchors.begin(), in_anchors.end(), [](const auto &a, const auto &b) {
        return a.get_orientation() < b.get_orientation();
    }));

    size_t num_seeds = 0;
    size_t num_extensions = 0;
    size_t num_explored_nodes = 0;

    num_seeds += in_anchors.size();
    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    tsl::hopscotch_map<std::string_view::const_iterator, size_t> end_counter;

    std::vector<Alignment> alignments;
    std::vector<Seed> seeds;
    seeds.reserve(in_anchors.size());
    size_t orientation_change = std::numeric_limits<size_t>::max();
    if (!labeled_aligner) {
        for (const auto &in_anchor : in_anchors) {
            if (!in_anchor.get_clipping() && !in_anchor.get_end_clipping()) {
                alignments.emplace_back(in_anchor, config_);
            } else {
                if (seeds.size() && in_anchor.get_orientation() != seeds.back().get_orientation())
                    orientation_change = seeds.size();

                seeds.emplace_back(in_anchor);
           }
        }
    } else {
        for (const auto &in_anchor : in_anchors) {
            if (!in_anchor.get_clipping() && !in_anchor.get_end_clipping()) {
                alignments.emplace_back(in_anchor, config_);
                continue;

            }
            ++end_counter[in_anchor.get_query_view().end()];
            const auto &cols = in_anchor.get_columns();
            for (auto col : cols) {
                if (seeds.size() && in_anchor.get_orientation() != seeds.back().get_orientation())
                    orientation_change = seeds.size();

                auto &anchor = seeds.emplace_back(in_anchor);
                anchor.set_columns({ col });
            }
        }
    }

    if (seeds.size() <= 1) {
        if (seeds.size() == 1 && seeds[0].get_query_view().size() > config_.min_seed_length)
            alignments.emplace_back(seeds[0], config_);

        seeder = std::make_unique<ManualSeeder>(std::move(alignments), query_size);
        return {};
    }

    orientation_change = std::min(orientation_change, seeds.size());

    auto preprocess_anchors = [&](auto begin, auto end) {
        if (begin == end)
            return;

        std::sort(begin, end, [](const auto &a, const auto &b) {
            return std::pair(b.label_columns, a.get_query_view().end())
                    > std::pair(a.label_columns, b.get_query_view().end());
        });

        if (end_counter.empty())
            return;

        // merge into MUMs
        for (auto i = end - 1; i != begin; --i) {
            if (i->label_columns == (i - 1)->label_columns
                    && i->get_query_view().end() + 1 == (i - 1)->get_query_view().end()
                    && end_counter[i->get_query_view().end()] == 1
                    && end_counter[i->get_query_view().end() + 1] == 1
                    && graph_.traverse(i->get_nodes().back(),
                                       (i - 1)->get_query_view().front())
                        == (i - 1)->get_nodes()[0]) {
                // we have a MUM
                assert(i->get_query_view().end() < (i - 1)->get_query_view().end());
                i->expand((i - 1)->get_nodes());
                std::swap(*i, *(i - 1));
                *i = Seed();
            }
        }
    };

    preprocess_anchors(seeds.begin(), seeds.begin() + orientation_change);
    preprocess_anchors(seeds.begin() + orientation_change, seeds.end());

    seeds.erase(std::remove_if(seeds.begin(), seeds.end(),
                               [](const auto &a) { return a.empty(); }),
                seeds.end());

    tsl::hopscotch_map<node_index, std::vector<node_index>> out_nodes;
    std::vector<node_index> nodes;
    std::vector<std::pair<size_t, size_t>> anchor_ends;
    anchor_ends.reserve(seeds.size());
    nodes.reserve(seeds.size());
    out_nodes.reserve(seeds.size());
    for (const auto &seed : seeds) {
        graph_.call_outgoing_kmers(seed.get_nodes().back(), [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel)
                out_nodes[seed.get_nodes().back()].emplace_back(next);
        });

        DEBUG_LOG("Anchor: {}", Alignment(seed, config_));

        auto &[anchor_front_idx, anchor_back_idx] = anchor_ends.emplace_back(
            nodes.size(), nodes.size()
        );
        nodes.emplace_back(seed.get_nodes().front());
        if (seed.get_nodes().size() > 1) {
            ++anchor_back_idx;
            nodes.emplace_back(seed.get_nodes().back());
        }
    }

    tsl::hopscotch_map<size_t, ColumnPathIndex::NodesInfo> node_col_coords;
    for (auto &[label, nodes_info] : column_path_index->get_chain_info(nodes)) {
        auto encode = label.size() ? label_encoder.encode(label)
                                   : std::numeric_limits<Seed::Column>::max();
        node_col_coords.emplace(encode, std::move(nodes_info));
    }

    float sl = -static_cast<float>(config_.min_seed_length) * 0.01;

    sdsl::bit_vector matching_pos[2] {
        sdsl::bit_vector(query_size, false),
        sdsl::bit_vector(query_size, false)
    };
    bool terminate = false;
    score_t best_score = std::numeric_limits<score_t>::min();

    std::unordered_multiset<Chain, ChainHash> chains;
    score_t last_chain_score = std::numeric_limits<score_t>::min();
    auto flush_chains = [&]() {
        if (chains.empty())
            return;

        auto callback = [&](Chain&& chain, score_t chain_score, auto &extender, auto &bwd_extender) {
            bool orientation = chain[0].first.get_orientation();
            best_score = std::max(best_score, chain_score);
            bool score_too_low = (chain_score < best_score * config_.rel_score_cutoff);
            if (score_too_low) {
                terminate = true;
                return;
            }

            logger->trace("Chain: {}", chain_score);
            for (const auto &[aln, dist] : chain) {
                logger->trace("\t{} (dist: {})", aln, dist);
            }

            bool added = false;

            aligner.extend_chain(std::move(chain), extender, [&](Alignment&& aln) {
                std::vector<Alignment> alns;
                if (!aln.get_end_clipping()) {
                    DEBUG_LOG("\t\t{}", aln);
                    added |= aln.get_cigar().mark_exact_matches(matching_pos[orientation]);
                    alns.emplace_back(std::move(aln));
                } else {
                    alns = extender.get_extensions(aln, 0, true);
                }

                for (auto&& ext : alns) {
                    if (!ext.get_clipping()) {
                        added |= ext.get_cigar().mark_exact_matches(matching_pos[orientation]);
                        DEBUG_LOG("\t\t{}", ext);
                        alignments.emplace_back(std::move(ext));
                    } else {
                        bwd_extender.rc_extend_rc(ext, [&](Alignment&& aln) {
                            assert(aln.is_valid(graph_, &config_));
                            for (node_index node : aln.get_nodes()) {
                                extender.filter_nodes(node, aln.get_clipping(),
                                                      query_size - aln.get_end_clipping());
                            }

                            added |= aln.get_cigar().mark_exact_matches(matching_pos[orientation]);
                            DEBUG_LOG("\t\t{}", aln);
                            alignments.emplace_back(std::move(aln));
                        }, true, 0);
                    }
                }
            }, true);

            terminate |= score_too_low && !added;
        };

        auto it = chains.begin();
        Chain last_chain = *it;
        for (++it; it != chains.end(); ++it) {
            const Chain &chain = *it;
            if (chain != last_chain) {
                bool orientation = last_chain[0].first.get_orientation();
                callback(std::move(last_chain), last_chain_score,
                         orientation ? bwd_extender : extender,
                         orientation ? extender : bwd_extender);
                last_chain = *it;
                continue;
            }

            // if this chain has the same anchors as the last one, merge their labels
            for (size_t i = 0; i < chain.size(); ++i) {
                assert(chain[i].first.label_columns);
                Vector<Alignment::Column> columns;
                const auto &last_columns = last_chain[i].first.get_columns();
                const auto &cur_columns = chain[i].first.get_columns();
                std::set_union(last_columns.begin(), last_columns.end(),
                               cur_columns.begin(), cur_columns.end(),
                               std::back_inserter(columns));
                last_chain[i].first.set_columns(std::move(columns));
            }
        }

        bool orientation = last_chain[0].first.get_orientation();
        callback(std::move(last_chain), last_chain_score,
                 orientation ? bwd_extender : extender,
                 orientation ? extender : bwd_extender);

        chains.clear();
    };

    logger->trace("Chaining {} anchors for a query of length {}", seeds.size(), query_size);
    chain_anchors<Seed>(config_, seeds.data(), seeds.data() + seeds.size(),
        [&](const Seed &a_i,
            ssize_t max_coord_dist,
            const Seed *begin,
            const Seed *end,
            auto chain_scores,
            const auto &update_score) {
            size_t node_i = anchor_ends[&a_i - seeds.data()].second;
            const auto &score_i = std::get<0>(*(chain_scores - (begin - seeds.data()) + (&a_i - seeds.data())));
            std::string_view query_i = a_i.get_query_view();
            auto col_i = a_i.get_columns()[0];

            assert(node_col_coords.count(col_i));
            assert(node_col_coords[col_i].size() == nodes.size());
            const auto &coords_i_back = node_col_coords[col_i][node_i];
            assert(std::get<0>(coords_i_back.first));

            --chain_scores;
            std::for_each(begin, end, [&](const Seed &a_j) {
                // try to connect a_i to a_j
                assert(&a_i != &a_j);
                assert(a_i.get_orientation() == a_j.get_orientation());
                ++chain_scores;

                std::string_view query_j = a_j.get_query_view();

                // want a_i.begin < a_j.begin && a_i.end < a_j.end
                if (query_i.begin() >= query_j.begin() || query_i.end() >= query_j.end())
                    return;

                ssize_t dist = query_j.end() - query_i.end();
                score_t label_change_score = 0;
                auto col_j = a_j.get_columns()[0];
                if (labeled_aligner) {
                    if (!allow_label_change && col_i != col_j)
                        return;

                    label_change_score = labeled_aligner->get_label_change_score(col_i, col_j);
                    if (label_change_score == DBGAlignerConfig::ninf)
                        return;
                }

                ssize_t num_added = query_j.end() - std::max(query_j.begin(), query_i.end());
                assert(num_added >= 0);

                score_t score_j = std::get<0>(*chain_scores);
                score_t base_added_score = score_j + num_added + label_change_score;
                if (base_added_score <= score_i)
                    return;

                bool overlap = (query_i.end() > query_j.begin());

                if (num_added == 1) {
                    auto find = out_nodes.find(a_i.get_nodes().back());
                    if (find != out_nodes.end()) {
                        int64_t coord_dist = 1;
                        float gap = dist - coord_dist;
                        bool found = false;
                        for (node_index next : find->second) {
                            if (next == a_j.get_nodes().front()) {
                                score_t gap_cost = ceil(sl * gap - log2(gap + 1) * 0.5);
                                assert(gap > 0 || gap_cost == 0);
                                found |= update_score(base_added_score + gap_cost,
                                                      &a_j, coord_dist);
                                break;
                            }
                        }

                        if (found && gap == 0)
                            return;
                    }
                }

                size_t node_j = anchor_ends[&a_j - seeds.data()].first;

                assert(node_col_coords.count(col_j));
                assert(node_col_coords[col_j].size() == nodes.size());
                const auto &coords_j_front = node_col_coords[col_j][node_j];
                assert(std::get<0>(coords_j_front.first));

                const std::string &label = col_i == std::numeric_limits<Seed::Column>::max()
                    ? "" : label_encoder.decode(col_i);

                column_path_index->call_distances(label,
                                                  coords_i_back, coords_j_front,
                                                  [&](int64_t coord_dist) {
                    if (num_added > coord_dist || (!num_added && coord_dist == 0))
                        return;

                    float gap = std::abs(coord_dist - dist);
                    if (gap != 0 && overlap)
                        return;

                    score_t gap_cost = ceil(sl * gap - log2(gap + 1) * 0.5);
                    assert(gap > 0 || gap_cost == 0);

                    update_score(base_added_score + gap_cost, &a_j, coord_dist);
                }, max_coord_dist, 2);
            });
        },
        [&](const auto &chain, score_t score) {
            if (chain.empty() || (chain.size() == 1 && chain[0].first->get_query_view().size() <= config_.min_seed_length))
                return false;

            if (last_chain_score != score) {
                flush_chains();
                last_chain_score = score;
            }

            std::vector<std::pair<Seed, size_t>> merged_chain;
            merged_chain.reserve(chain.size());
            for (const auto &[seed_ptr, dist] : chain) {
                merged_chain.emplace_back(*seed_ptr, dist);
            }

            for (auto jt = merged_chain.rbegin(); jt + 1 != merged_chain.rend(); ++jt) {
                auto &prev_seed = (jt + 1)->first;
                auto &cur_seed = jt->first;
                auto prev_end = prev_seed.get_query_view().end();
                if (prev_seed.get_clipping() + 1 == cur_seed.get_clipping()
                        && prev_end + 1 == cur_seed.get_query_view().end() && jt->second == 1
                        && cur_seed.label_columns == prev_seed.label_columns) {
                    assert(graph_.traverse(prev_seed.get_nodes().back(), *prev_end) == cur_seed.get_nodes()[0]);
                    prev_seed.expand(cur_seed.get_nodes());
                    cur_seed = Seed();
                }
            }

            Chain cur_chain;
            for (const auto &[seed, dist] : merged_chain) {
                if (!seed.empty())
                    cur_chain.emplace_back(Alignment(seed, config_), dist);
            }

            assert(std::all_of(cur_chain.begin() + 1, cur_chain.end(),
                               [&](const auto &a) { return a.second > 0; }));

            chains.emplace(std::move(cur_chain));
            return true;
        },
        false,
        [](const auto*, auto&&, size_t, const auto&) {},
        [](auto&&) {},
        [&terminate]() { return terminate; },
        false,
        config_.max_dist_between_seeds,
        config_.max_gap_shrinking_factor
    );

    if (!terminate)
        flush_chains();

    num_extensions += extender.num_extensions() + bwd_extender.num_extensions();
    num_explored_nodes += extender.num_explored_nodes() + bwd_extender.num_explored_nodes();
    seeder = std::make_unique<ManualSeeder>(
        std::move(alignments),
        std::max(sdsl::util::cnt_one_bits(matching_pos[0]),
                 sdsl::util::cnt_one_bits(matching_pos[1]))
    );

    return std::make_tuple(num_seeds, num_extensions, num_explored_nodes);
}

void chain_alignments(const IDBGAligner &aligner,
                      const std::vector<Alignment> &alignments,
                      const std::function<void(Alignment&&)> &callback) {
    assert(std::is_sorted(alignments.begin(), alignments.end(),
                          [](const auto &a, const auto &b) {
                              return a.get_orientation() < b.get_orientation();
                          }));

    const auto &config = aligner.get_config();
    if (!config.allow_jump && !config.allow_label_change)
        return;

    if (alignments.size() <= 1
            || (alignments.size() == 2
                && alignments[1].get_orientation() != alignments[0].get_orientation())) {
        return;
    }

    if (std::any_of(alignments.begin(), alignments.end(),
                    [](const auto &a) { return !a.get_clipping() && !a.get_end_clipping(); })) {
        return;
    }

    const DeBruijnGraph &graph = aligner.get_graph();
    std::vector<std::vector<score_t>> per_char_scores_prefix;
    per_char_scores_prefix.reserve(alignments.size());

    tsl::hopscotch_map<std::string_view::const_iterator, size_t> end_counter;

    // preprocess alignments
    for (size_t i = 0; i < alignments.size(); ++i) {
        const auto &alignment = alignments[i];
        DEBUG_LOG("Alignment {}:\t{}", i, alignment);
        std::string_view query = alignment.get_query_view();
        auto &prefix_scores_with_deletions
            = per_char_scores_prefix.emplace_back(std::vector<score_t>(query.size() + 1, 0));

        auto cur = alignment;
        auto it = prefix_scores_with_deletions.begin();
        while (cur.size()) {
            cur.trim_query_prefix(1, graph.get_k() - 1, config);
            ++it;
            assert(it != prefix_scores_with_deletions.end());
            *it = alignment.get_score() - cur.get_score();
        }
        assert(prefix_scores_with_deletions.back() == alignment.get_score());
    }

    size_t seed_size = std::min(config.min_seed_length, graph.get_k());

    struct Anchor {
        std::string_view::const_iterator end;
        std::string_view::const_iterator begin;
        uint64_t index;
        int64_t aln_index_back;
        int64_t aln_index_front;
        std::string_view::const_iterator aln_begin;
        std::string_view::const_iterator aln_end;

        uint32_t last;
        uint64_t mem_length;
    };

    std::vector<Anchor> anchors;
    size_t orientation_change = std::numeric_limits<size_t>::max();

    for (size_t i = 0; i < alignments.size(); ++i) {
        const auto &alignment = alignments[i];
        if (i && alignments[i - 1].get_orientation() != alignment.get_orientation())
            orientation_change = anchors.size();

        auto add_anchor = [&](auto begin, auto end, ssize_t node_i) {
            ++end_counter[end];
            anchors.emplace_back(Anchor{
                .end = end,
                .begin = begin,
                .index = i,
                .aln_index_back = node_i,
                .aln_index_front = node_i,
                .aln_begin = alignment.get_query_view().begin(),
                .aln_end = alignment.get_query_view().end(),
                .last = std::numeric_limits<uint32_t>::max(),
                .mem_length = static_cast<uint64_t>(end - begin),
            });
        };

        auto cur = alignment;
        for ( ; cur.get_nodes().size() > 1; cur.trim_query_suffix(1, config)) {
            auto it = cur.get_cigar().data().rbegin();
            if (it->first == Cigar::CLIPPED)
                ++it;

            assert(it != cur.get_cigar().data().rend());
            if (it->first == Cigar::MATCH && it->second >= seed_size) {
                auto end = cur.get_query_view().end();
                auto begin = end - seed_size;
                ssize_t node_i = cur.get_nodes().size() - 1;
                add_anchor(begin, end, node_i);
            }
        }

        if (cur.get_nodes().size() != 1)
            continue;

        auto it = cur.get_cigar().data().rbegin();
        if (it->first == Cigar::CLIPPED)
            ++it;

        assert(it != cur.get_cigar().data().rend());
        if (it->first == Cigar::INSERTION)
            continue;

        if (it->first == Cigar::MATCH && it->second >= seed_size) {
            auto end = cur.get_query_view().end();
            auto begin = end - seed_size;
            ssize_t node_i = 0;
            add_anchor(begin, end, node_i);
        }

        for ( ; cur.get_query_view().size() > seed_size; cur.trim_query_prefix(1, graph.get_k() - 1, config)) {
            auto jt = cur.get_cigar().data().begin();
            if (jt->first == Cigar::CLIPPED)
                ++jt;

            if (jt->first == Cigar::MATCH && jt->second >= seed_size) {
                auto begin = cur.get_query_view().begin();
                auto end = begin + seed_size;
                ssize_t node_i = -static_cast<ssize_t>(cur.get_sequence().size()) + seed_size;
                add_anchor(begin, end, node_i);
            }
        }
    }

    orientation_change = std::min(orientation_change, anchors.size());

    if (orientation_change <= 1 && anchors.size() - orientation_change <= 1)
        return;

    auto preprocess_anchors = [&](auto begin, auto end) {
        if (begin == end)
            return;

        std::sort(begin, end, [](const auto &a, const auto &b) {
            return std::tie(a.end, a.aln_begin) > std::tie(b.end, b.aln_begin);
        });
        auto rbegin = std::make_reverse_iterator(end);
        auto rend = std::make_reverse_iterator(begin);
        for (auto it = rbegin; it + 1 != rend; ++it) {
            assert(alignments[it->index].get_orientation()
                    == alignments[(it + 1)->index].get_orientation());
            if ((it + 1)->index == it->index
                    && it->aln_index_back + 1 == (it + 1)->aln_index_front
                    && it->end + 1 == (it + 1)->end
                    && end_counter[it->end] == 1
                    && end_counter[it->end + 1] == 1) {
                // we have a MUM
                (it + 1)->aln_index_front = it->aln_index_front;
                (it + 1)->begin = it->begin;
                (it + 1)->mem_length = (it + 1)->end - (it + 1)->begin;

                // clear out this anchor
                it->index = std::numeric_limits<uint64_t>::max();
            }
        }
    };
    preprocess_anchors(anchors.begin(), anchors.begin() + orientation_change);
    preprocess_anchors(anchors.begin() + orientation_change, anchors.end());

    anchors.erase(std::remove_if(anchors.begin(), anchors.end(),
                                 [&](const auto &a) {
                                     return a.index == std::numeric_limits<uint64_t>::max();
                                 }),
                  anchors.end());

    struct AnchorExtraInfo {
        uint64_t index;
        int64_t aln_index_back;
        int64_t aln_index_front;

        int64_t last_dist;
        uint64_t mem_length;
        score_t label_change_score;
    };
    std::vector<Alignment> anchor_alns;
    std::vector<AnchorExtraInfo> anchor_extra_info;
    anchor_alns.reserve(anchors.size());
    anchor_extra_info.reserve(anchors.size());

    for (const auto &anchor : anchors) {
        auto &aln = anchor_alns.emplace_back(alignments[anchor.index]);
        if (aln.get_offset() != graph.get_k() - 1) {
            aln.extend_offset(std::vector<node_index>(graph.get_k() - 1 - aln.get_offset(),
                                                      DeBruijnGraph::npos));
        }

        aln.trim_query_suffix(aln.get_query_view().end() - anchor.end, config);
        aln.trim_query_prefix(anchor.begin - aln.get_query_view().begin(), graph.get_k() - 1, config);

        DEBUG_LOG("Seq: {}\tAnchor: {}", anchor.index, aln);
        anchor_extra_info.emplace_back(AnchorExtraInfo{
            .index = anchor.index,
            .aln_index_back = anchor.aln_index_back,
            .aln_index_front = anchor.aln_index_front,
            .last_dist = 0,
            .mem_length = anchor.mem_length,
            .label_change_score = DBGAlignerConfig::ninf,
        });
    }

    size_t num_found = 0;
    score_t node_insert = config.node_insertion_penalty;
    score_t gap_open = config.gap_opening_penalty;
    score_t gap_ext = config.gap_extension_penalty;
    assert(gap_open < 0);
    assert(gap_ext < 0);
    assert(gap_ext >= gap_open);
    assert(node_insert < 0);

    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);

    size_t last_index;
    size_t last_anchor;
    score_t chain_score;
    Alignment start_back_aln;
    chain_anchors<Alignment>(config, anchor_alns.data(), anchor_alns.data() + anchor_alns.size(),
        [&](const Alignment &a_i,
            ssize_t,
            const Alignment *begin,
            const Alignment *end,
            auto chain_scores,
            const auto &update_score) {
            auto &info_i = anchor_extra_info[&a_i - anchor_alns.data()];
            score_t &score_i = std::get<0>(*(
                chain_scores - (begin - anchor_alns.data()) + (&a_i - anchor_alns.data())
            ));
            std::string_view full_query_i = alignments[info_i.index].get_query_view();
            std::string_view query_i = a_i.get_query_view();
            const auto &prefix_scores_with_deletions_i = per_char_scores_prefix[info_i.index];

            --chain_scores;
            std::for_each(begin, end, [&](const Alignment &a_j) {
                // try to connect a_i to a_j
                ++chain_scores;
                assert(a_i.get_orientation() == a_j.get_orientation());
                if (&a_i == &a_j)
                    return;

                const auto &info_j = anchor_extra_info[&a_j - anchor_alns.data()];

                const auto &prefix_scores_with_deletions_j = per_char_scores_prefix[info_j.index];
                std::string_view query_j = a_j.get_query_view();
                std::string_view full_query_j = alignments[info_j.index].get_query_view();

                auto [score_j, last, last_dist] = *chain_scores;
                bool is_start = (last == anchor_alns.data() + anchor_alns.size());

                if (is_start) {
                    score_j = alignments[info_j.index].get_score()
                                - prefix_scores_with_deletions_j[query_j.begin() - full_query_j.begin()];
                }

                if (info_i.index == info_j.index) {
                    assert(info_j.aln_index_back >= info_i.aln_index_back);
                    score_t updated_score = is_start
                        ? alignments[info_i.index].get_score()
                                - prefix_scores_with_deletions_i[query_i.begin() - full_query_i.begin()]
                        : score_j + prefix_scores_with_deletions_i[query_j.begin() - full_query_j.begin()]
                            - prefix_scores_with_deletions_i[query_i.begin() - full_query_i.begin()];

                    if (update_score(updated_score, &a_j, 0)) {
                        size_t num_added = info_j.aln_index_front - info_i.aln_index_front;
                        info_i.mem_length = info_j.mem_length + num_added;
                        info_i.label_change_score = 0;
                    }

                    return;
                }

                auto get_label_change_score = [&](auto a_i_col, auto a_j_col,
                                                  std::string_view c) {
                    if (!labeled_aligner || a_i_col == a_j_col)
                        return 0;

                    score_t label_change_score = DBGAlignerConfig::ninf;

                    if (!a_i_col || !a_j_col
                            || (!config.allow_label_change && a_i_col != a_j_col)) {
                        return label_change_score;
                    }

                    score_t match_score = config.match_score(c);

                    for (auto&& [labels, lc_score, is_subset] : labeled_aligner->get_label_change_scores(a_i_col, a_j_col)) {
                        label_change_score = std::max(label_change_score,
                                                      lc_score * match_score);
                    }

                    return label_change_score;
                };

                if (config.allow_jump && full_query_i.end() <= full_query_j.begin()) {
                    // completely disjoint
                    score_t gap = full_query_j.begin() - full_query_i.end();
                    if (info_j.mem_length >= graph.get_k()) {
                        score_t gap_cost = node_insert + gap_open;
                        if (gap > 0)
                            gap_cost += gap_open + (gap - 1) * gap_ext;

                        assert(gap_cost < 0);

                        score_t base_updated_score = score_j + gap_cost
                            + alignments[info_i.index].get_score()
                            - prefix_scores_with_deletions_i[query_i.begin() - full_query_i.begin()];

                        if (base_updated_score <= score_i)
                            return;

                        score_t label_change_score = get_label_change_score(
                            a_i.label_column_diffs.size() ? a_i.label_column_diffs.back()
                                                          : a_i.label_columns,
                            a_j.label_columns,
                            std::string_view(full_query_j.begin(), 1)
                        );

                        score_t updated_score = base_updated_score + label_change_score;

                        if (update_score(updated_score, &a_j, 0)) {
                            info_i.mem_length = query_i.size();
                            info_i.label_change_score = label_change_score;
                        }
                    }

                    return;
                }

                score_t gap = query_j.begin() - query_i.end();
                if (config.allow_jump && gap >= 0) {
                    // alignments overlap, but there's no overlapping k-mer
                    return;
                }

                if (query_j.end() != query_i.end())
                    return;

                if (info_i.aln_index_front < static_cast<ssize_t>(alignments[info_i.index].get_offset()) + 1)
                    return;

                score_t base_updated_score = score_j + a_i.get_score() - a_j.get_score();

                auto update_score_with_labels = [&]() {
                    if (base_updated_score <= score_i)
                        return;

                    score_t label_change_score = get_label_change_score(
                        a_i.label_column_diffs.size() ? a_i.label_column_diffs.back()
                                                      : a_i.label_columns,
                        a_j.label_columns,
                        std::string_view(query_j.begin(), 1)
                    );

                    score_t updated_score = base_updated_score + label_change_score;
                    if (update_score(updated_score, &a_j, 0)) {
                        info_i.mem_length = query_i.size();
                        info_i.label_change_score = label_change_score;
                    }
                };

                if (a_i.get_nodes().back() == a_j.get_nodes().back()
                        && info_j.mem_length > query_j.size()) {
                    // perfect overlap, easy to connect
                    assert(query_i.size() == query_j.size());
                    update_score_with_labels();
                    return;
                }

                if (config.allow_jump && info_j.mem_length >= graph.get_k()) {
                    assert(query_i.end() > query_j.begin());
                    base_updated_score += node_insert;
                    update_score_with_labels();
                }
            });
        },
        [&](const auto &chain, score_t score) {
            if (chain.size() <= 1)
                return false;

            chain_score = score;
            DEBUG_LOG("Chain: {}", score);

            bool all_equal = true;
            DEBUG_LOG("\t{} (aln: {}; length: {})",
                      *chain[0].first,
                      anchor_extra_info[chain[0].first - anchor_alns.data()].index,
                      anchor_extra_info[chain[0].first - anchor_alns.data()].mem_length);
            for (size_t i = 1; i < chain.size(); ++i) {
                const auto &info = anchor_extra_info[chain[i].first - anchor_alns.data()];
                DEBUG_LOG("\t{} (aln: {}; dist: {}; length: {})",
                          *chain[i].first, info.index, chain[i].second, info.mem_length);
                all_equal &= (info.index
                                == anchor_extra_info[chain[i - 1].first - anchor_alns.data()].index);
            }

            if (all_equal) {
                DEBUG_LOG("\tSkipping: all from same alignment");
                return false;
            }

            last_anchor = chain.back().first - anchor_alns.data();
            last_index = anchor_extra_info[last_anchor].index;
            const Alignment *start = chain[0].first;
            const auto &start_extra_info = anchor_extra_info[start - anchor_alns.data()];
            if (start_extra_info.mem_length < graph.get_k()) {
                DEBUG_LOG("\tSkipping: first alignment fragment too short ({} < {})",
                          start_extra_info.mem_length, graph.get_k());
                return false;
            }

            start_back_aln = alignments[anchor_extra_info[chain.back().first - anchor_alns.data()].index];

            return true;
        },
        true,
        [&](const Alignment *first, Alignment&& cur, size_t, const auto &callback) {
            if (start_back_aln.size()) {
                std::swap(cur, start_back_aln);
                start_back_aln = Alignment();
            }

            ssize_t overlap = first->get_query_view().end() - anchor_alns[last_anchor].get_query_view().begin();
            last_anchor = first - anchor_alns.data();
            const auto &first_extra_info = anchor_extra_info[last_anchor];

            if (last_index == first_extra_info.index) {
                DEBUG_LOG("\tCurrent: {}", cur);
                callback(std::move(cur));
                return;
            }

            last_index = first_extra_info.index;

            Alignment alignment = alignments[last_index];
            DEBUG_LOG("\tMerging in: {}", alignment);
            assert(alignment.get_query_view().begin() <= first->get_query_view().begin());
            assert(alignment.get_query_view().end() >= first->get_query_view().end());
            if (overlap <= 0) {
                assert(alignment.get_query_view().end() <= cur.get_query_view().begin() && "Not implemented");
                cur.insert_gap_prefix(cur.get_query_view().begin() - alignment.get_query_view().end(), graph.get_k() - 1, config);
                assert(cur.size());
                // assert(cur.is_valid(graph, &config));
            } else {
                cur.trim_query_prefix(anchor_alns[last_anchor].get_query_view().begin() - cur.get_query_view().begin(),
                                      graph.get_k() - 1, config);
                assert(cur.get_query_view().begin() == anchor_alns[last_anchor].get_query_view().begin());

                assert(first->get_query_view().begin() == cur.get_query_view().begin());
                cur.extend_offset(std::vector<node_index>(graph.get_k() - 1 - cur.get_offset(),
                                                          DeBruijnGraph::npos));
                bool insert_gap_prefix = (cur.get_nodes()[overlap - 1] != first->get_nodes().back());

                cur.trim_query_prefix(overlap, graph.get_k() - 1, config, false);
                assert(cur.size());
                assert(cur.is_valid(graph, &config));

                if (insert_gap_prefix) {
                    cur.insert_gap_prefix(-overlap, graph.get_k() - 1, config);
                    assert(cur.size());
                    assert(cur.is_valid(graph, &config));
                }

                alignment.extend_offset(std::vector<node_index>(graph.get_k() - 1 - alignment.get_offset(),
                                                                DeBruijnGraph::npos));
                alignment.trim_query_suffix(alignment.get_query_view().end() - cur.get_query_view().begin(),
                                            config);
                assert(alignment.size());
            }

            alignment.splice(std::move(cur), first_extra_info.label_change_score);
            DEBUG_LOG("\tCurrent: {}", alignment);
            assert(alignment.size());
            assert(alignment.is_valid(graph, &config));
            callback(std::move(alignment));
        },
        [&](Alignment&& aln) {
            ++num_found;
            aln.trim_offset();
            DEBUG_LOG("\tFinal: {}\t{}", chain_score, aln);
            assert(anchor_alns[last_anchor].get_query_view().begin() >= aln.get_query_view().begin());
            assert(aln.get_score()
                - per_char_scores_prefix[last_index][anchor_alns[last_anchor].get_query_view().begin() - aln.get_query_view().begin()] == chain_score);
            callback(std::move(aln));
        },
        [&]() { return num_found >= config.num_alternative_paths; },
        true,
        config.max_dist_between_seeds,
        config.max_gap_shrinking_factor
    );
}

} // namespace align
} // namespace graph
} // namespace mtg
