#include "aligner_chainer.hpp"

#include <unordered_set>

#include <x86/svml.h>

#include "aligner_seeder_methods.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"
#include "chainer.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"
#include "graph/graph_extensions/graph_topology.hpp"

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

typedef std::vector<TableElem> ChainDPTable;

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
                if (chain[i].first.label_coordinates.size()) {
                    assert(last_chain[i].first.get_columns().size()
                            == last_chain[i].first.label_coordinates.size());
                    assert(chain[i].first.get_columns().size()
                            == chain[i].first.label_coordinates.size());
                    Alignment::CoordinateSet coord_union;
                    auto add_col_coords = [&](auto col, auto &coords) {
                        columns.push_back(col);
                        coord_union.emplace_back(std::move(coords));
                    };
                    utils::match_indexed_values(
                        last_chain[i].first.get_columns().begin(),
                        last_chain[i].first.get_columns().end(),
                        last_chain[i].first.label_coordinates.begin(),
                        chain[i].first.get_columns().begin(),
                        chain[i].first.get_columns().end(),
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
                    std::set_union(last_chain[i].first.get_columns().begin(),
                                   last_chain[i].first.get_columns().end(),
                                   chain[i].first.get_columns().begin(),
                                   chain[i].first.get_columns().end(),
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
                chain_seeds.back().first.set_columns(Vector<Alignment::Column>(1, label));
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
        const auto &seed_columns = seeds[i].get_columns();
        for (size_t j = 0; j < seeds[i].label_coordinates.size(); ++j) {
            Alignment::Column c = seed_columns[j];
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

    size_t bandwidth = config.chaining_bandwidth;

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

std::pair<size_t, size_t>
cluster_seeds(const IDBGAligner &aligner,
              std::string_view forward,
              std::string_view reverse,
              const DBGAlignerConfig &config,
              std::vector<Seed>&& fwd_seeds,
              std::vector<Seed>&& bwd_seeds,
              const std::function<void(Alignment&&)> &callback,
              const std::function<void(Seed&&)> &discarded_seed_callback,
              const std::function<bool(Alignment::Column)> &skip_column) {
    const auto *seq_annotator = aligner.get_seq_annotator();
    if (!seq_annotator) {
        for (const auto &seed : fwd_seeds) {
            callback(Alignment(seed, config));
        }

        for (const auto &seed : bwd_seeds) {
            callback(Alignment(seed, config));
        }

        return std::make_pair(0, 0);
    }

    const auto &graph = aligner.get_graph();

    std::vector<node_index> nodes;
    nodes.reserve(fwd_seeds.size() + bwd_seeds.size());
    for (auto &fwd_seed : fwd_seeds) {
        assert(fwd_seed.get_nodes().back());
        nodes.emplace_back(fwd_seed.get_nodes().back());
    }
    for (auto &bwd_seed : bwd_seeds) {
        assert(bwd_seed.get_nodes().back());
        nodes.emplace_back(bwd_seed.get_nodes().back());
    }
    assert(nodes.size() == fwd_seeds.size() + bwd_seeds.size());

    size_t seed_count = 0;

    struct Anchor {
        Anchor(const Seed &seed, size_t tuple_id, const DBGAlignerConfig &config)
            : seed(seed), tuple_id(tuple_id), used(false), score_(seed.get_score(config)) {}

        score_t get_score(const DBGAlignerConfig &) const { return score_; }

        size_t get_clipping() const { return seed.get_clipping(); }
        size_t get_end_clipping() const { return seed.get_end_clipping(); }
        std::string_view get_query_view() const { return seed.get_query_view(); }
        bool get_orientation() const { return seed.get_orientation(); }

        Seed seed;
        size_t tuple_id;
        std::vector<std::pair<size_t, int64_t>> coords;
        mutable bool used;
        score_t score_;
    };

    tsl::hopscotch_map<Alignment::Column, std::vector<Anchor>> clustered_seed_map;

    logger->trace("Fetching node coordinates");
    auto [seq_ids, coords] = seq_annotator->get_seq_ids(nodes);
    assert(seq_ids.size() == nodes.size());
    assert(coords.size() == nodes.size());

    if (nodes.size()) {
        logger->trace("Clustering anchors by label");
        auto it = seq_ids.begin();
        auto jt = coords.begin();

        auto parse_set = [&](const auto &seeds) {
            for (const auto &seed : seeds) {
                assert(it != seq_ids.end());
                assert(jt != coords.end());
                assert(it->size() == jt->size());
                size_t tuple_id = jt - coords.begin();
                for (size_t i = 0; i < it->size(); ++i) {
                    const auto &[col, seq_id] = (*it)[i];
                    assert(col == (*jt)[i].first);
                    assert(seq_id.size() == 1);
                    const auto &seqs = seq_id[0];

                    const auto &coords = (*jt)[i].second;
                    assert(seqs.size() == coords.size());

                    auto &cur = clustered_seed_map[col].emplace_back(seed, tuple_id, config).coords;

                    cur.reserve(coords.size());
                    for (size_t j = 0; j < coords.size(); ++j) {
                        cur.emplace_back(seqs[j], coords[j]);
                    }
                }
                ++it;
                ++jt;
            }
        };

        parse_set(fwd_seeds);
        parse_set(bwd_seeds);

        for (auto it = clustered_seed_map.begin(); it != clustered_seed_map.end(); ++it) {
            auto &cluster = it.value();
            std::sort(cluster.begin(), cluster.end(), [&](const auto &a, const auto &b) {
                return std::make_pair(a.get_orientation(), a.seed.get_query_view().end())
                        < std::make_pair(b.get_orientation(), b.seed.get_query_view().end());
            });
        }

        assert(it == seq_ids.end());
        assert(jt == coords.end());
    }

    if (clustered_seed_map.empty()) {
        for (const auto &seed : fwd_seeds) {
            callback(Alignment(seed, config));
        }

        for (const auto &seed : bwd_seeds) {
            callback(Alignment(seed, config));
        }

        return std::make_pair(0, 0);
    }

    logger->trace("Sorting anchor clusters");

    std::vector<std::pair<Alignment::Column, score_t>> max_scores;
    max_scores.reserve(clustered_seed_map.size());
    for (const auto &[col, anchors] : clustered_seed_map) {
        sdsl::bit_vector covered_fwd(forward.size(), false);
        sdsl::bit_vector covered_bwd(reverse.size(), false);

        for (const Anchor &anchor : anchors) {
            if (!anchor.get_orientation()) {
                std::fill(covered_fwd.begin() + anchor.get_clipping(),
                          covered_fwd.end() - anchor.get_end_clipping(),
                          true);
            } else {
                std::fill(covered_bwd.begin() + anchor.get_clipping(),
                          covered_bwd.end() - anchor.get_end_clipping(),
                          true);
            }
        }

        score_t max_score_fwd = 0;
        score_t max_score_bwd = 0;
        auto it = covered_fwd.begin();
        while (it != covered_fwd.end()) {
            it = std::find(it, covered_fwd.end(), true);
            auto next_it = std::find(it, covered_fwd.end(), false);
            if (next_it > it) {
                std::string_view window(forward.data() + (it - covered_fwd.begin()),
                                        next_it - it);
                max_score_fwd += config.match_score(window);
            }

            it = next_it;
        }
        it = covered_bwd.begin();
        while (it != covered_bwd.end()) {
            it = std::find(it, covered_bwd.end(), true);
            auto next_it = std::find(it, covered_bwd.end(), false);
            if (next_it > it) {
                std::string_view window(reverse.data() + (it - covered_bwd.begin()),
                                        next_it - it);
                max_score_bwd += config.match_score(window);
            }

            it = next_it;
        }

        max_scores.emplace_back(col, std::max(max_score_fwd, max_score_bwd));
    }

    std::sort(max_scores.begin(), max_scores.end(), utils::GreaterSecond());

    std::vector<std::pair<Alignment::Column, std::vector<Anchor>>> clustered_seeds;
    clustered_seeds.reserve(clustered_seed_map.size());

    for (const auto &[col, max_score] : max_scores) {
        auto it = clustered_seed_map.find(col);
        assert(it != clustered_seed_map.end());
        clustered_seeds.emplace_back(it->first, std::move(it.value()));
    }

    size_t max_num_seeds = config.num_alternative_paths;

    logger->trace("Chaining anchors from {} labels", clustered_seed_map.size());

    std::vector<std::tuple<Alignment::Column,
                        AnchorChain<Anchor>,
                        std::vector<score_t>>> best_chains;

    score_t best_last_score = 0;
    size_t num_skipped = 0;
    auto max_score_it = max_scores.begin();
    for (auto &[col, anchors] : clustered_seeds) {
        assert(max_score_it != max_scores.end());
        assert(col == max_score_it->first);
        if (skip_column(col) || max_score_it->second < best_last_score) {
            ++num_skipped;
            ++max_score_it;
            continue;
        }

        size_t count = 0;
        score_t last_score = std::numeric_limits<score_t>::max();

        chain_anchors<Anchor>(config, anchors.data(), anchors.data() + anchors.size(),
            [&config](const Anchor &a_j,
                            ssize_t,
                            const Anchor *begin,
                            const Anchor *end,
                            auto *chain_scores,
                            const auto &update_score) {
                std::string_view query_j = a_j.get_query_view();
                const score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));

                std::for_each(begin, end, [&](const Anchor &a_i) {
                    std::string_view query_i = a_i.get_query_view();

                    score_t dist = query_j.end() - query_i.end();

                    if (dist <= 0 || query_i.begin() >= query_j.begin()) {
                        ++chain_scores;
                        return;
                    }

                    score_t base_score = std::get<0>(*chain_scores);

                    if (query_i.end() <= query_j.begin()) {
                        base_score += a_j.get_score(config);
                    } else {
                        std::string_view ext(query_i.data() + query_i.size(),
                                             query_j.end() - query_i.end());
                        base_score += config.match_score(ext);

                        if (!a_j.get_end_clipping())
                            base_score += config.right_end_bonus;
                    }


                    if (base_score <= score_j) {
                        ++chain_scores;
                        return;
                    }

                    auto it = a_i.coords.begin();
                    auto jt = a_j.coords.begin();

                    score_t min_diff = std::numeric_limits<score_t>::max();
                    while (it != a_i.coords.end() && jt != a_j.coords.end()) {
                        const auto &[s_i, c_i] = *it;
                        const auto &[s_j, c_j] = *jt;
                        if (s_i < s_j) {
                            ++it;
                        } else if (s_i > s_j || c_i >= c_j) {
                            ++jt;
                        } else {
                            // this coordinate pair works
                            score_t score = base_score;
                            score_t coord_dist = c_j - c_i;

                            if (coord_dist != dist) {
                                score_t diff = std::abs(coord_dist - dist);
                                min_diff = std::min(min_diff, diff);
                                score += config.gap_opening_penalty
                                        + (diff - 1) * config.gap_extension_penalty;
                            } else {
                                min_diff = 0;
                            }

                            update_score(score, &a_i, coord_dist);

                            if (min_diff == 0)
                                break;

                            ++it;
                        }
                    }

                    ++chain_scores;
                });
            },
            [&](const AnchorChain<Anchor> &chain, const std::vector<score_t> &score_traceback) {
                if (score_traceback.back() < last_score) {
                    if (++count > max_num_seeds) {
                        assert(last_score != std::numeric_limits<score_t>::max());
                        best_last_score = std::max(best_last_score, last_score);
                        last_score = score_traceback.back();
                        return false;
                    } else {
                        last_score = score_traceback.back();
                    }
                }

                best_chains.emplace_back(col, chain, score_traceback);
                return true;
            },
            false,
            [](const auto*, const auto*, auto&&, auto, auto, const auto&) {},
            [](auto&&) {},
            [&]() { return count > max_num_seeds; }
        );

        ++max_score_it;
    }

    assert(max_score_it == max_scores.end());

    if (num_skipped)
        logger->trace("Skipping {} columns", num_skipped);

    if (best_chains.empty()) {
        logger->trace("Reduced {} seeds down to {}", fwd_seeds.size() + bwd_seeds.size(), seed_count);
        return std::make_pair(0, 0);
    }

    struct AnchorChainHash {
        inline std::size_t operator()(const AnchorChain<Anchor> &chain) const {
            uint64_t hash = 0;
            for (const auto &[aln, dist] : chain) {
                for (node_index node : aln->seed.get_nodes()) {
                    hash ^= node + 0x9e3779b9 + (hash << 6) + (hash >> 2);
                }
                hash ^= dist + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };

    struct AnchorChainEqual {
        inline std::size_t operator()(const AnchorChain<Anchor> &a,
                                    const AnchorChain<Anchor> &b) const {
            if (a.size() != b.size())
                return false;

            if (a.get_score() != b.get_score())
                return false;

            for (size_t i = 0; i < a.size(); ++i) {
                if (a[i].second != b[i].second || !(a[i].first->seed == b[i].first->seed))
                    return false;
            }

            return true;
        }
    };

    assert(reverse.empty() || reverse.size() == forward.size());
    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    if (labeled_aligner) {
        auto *anno_buffer = &labeled_aligner->get_annotation_buffer();
        tsl::hopscotch_map<AnchorChain<Anchor>, Vector<Alignment::Column>,
                            AnchorChainHash, AnchorChainEqual> chain_map;
        tsl::hopscotch_map<Alignment::Column, sdsl::bit_vector> col_coverage;

        for (const auto &[col, chain, score_traceback] : best_chains) {
            if (skip_column(col))
                continue;

            auto &coverage = col_coverage[col];
            if (coverage.empty())
                coverage = sdsl::bit_vector(forward.size(), false);

            auto cov_begin = coverage.begin() + chain.get_clipping();
            auto cov_end = cov_begin + chain.get_query_view().size();
            if (std::all_of(cov_begin, cov_end, std::logical_not())) {
                std::fill(cov_begin, cov_end, true);
                chain_map[chain].emplace_back(col);
            }
        }

        std::vector<std::tuple<Alignment::Column, AnchorChain<Anchor>, std::vector<score_t>>> merged_best_chains;
        for (auto it = chain_map.begin(); it != chain_map.end(); ++it) {
            auto &chain = it->first;
            auto &cols = it.value();
            std::sort(cols.begin(), cols.end());

            size_t size = chain.size();
            merged_best_chains.emplace_back(
                anno_buffer->cache_column_set(cols.begin(), cols.end()),
                std::move(chain),
                std::vector<score_t>(size, 0)
            );
        }

        std::swap(merged_best_chains, best_chains);
    }

    std::sort(best_chains.begin(), best_chains.end(), [&](const auto &a, const auto &b) {
        const auto &chain_a = std::get<1>(a);
        const auto &chain_b = std::get<1>(b);
        return std::make_tuple(chain_a.get_score(), chain_b.get_clipping(), chain_b.get_orientation())
                > std::make_tuple(chain_b.get_score(), chain_a.get_clipping(), chain_a.get_orientation());
    });

    if (best_chains.size() > max_num_seeds) {
        logger->trace("Found {} initial chains with scores from {} to {}",
                        best_chains.size(),
                        std::get<1>(best_chains.front()).get_score(),
                        std::get<1>(best_chains.back()).get_score());
        size_t count = 0;
        for (auto it = best_chains.begin() + 1; it != best_chains.end(); ++it) {
            score_t score_a = std::get<1>(*(it - 1)).get_score();
            score_t score_b = std::get<1>(*it).get_score();
            if (score_a > score_b) {
                ++count;
                if (count == max_num_seeds) {
                    best_chains.erase(it, best_chains.end());
                    break;
                }
            }
        }
    }

    logger->trace("Found {} chains with scores from {} to {}",
                  best_chains.size(),
                  std::get<1>(best_chains.front()).get_score(),
                  std::get<1>(best_chains.back()).get_score());

    size_t num_extensions = 0;
    size_t num_explored_nodes = 0;

    for (const auto &[col, chain, score_traceback] : best_chains) {
        extend_chain<Anchor>(chain, score_traceback,
            [&](const Anchor *next,
                const Anchor *last_anchor,
                Alignment&& cur,
                size_t coord_dist,
                score_t,
                const auto &continue_callback) {

                if (cur.empty()) {
                    assert(!coord_dist);
                    assert(!last_anchor);

                    next->used = true;
                    cur = Alignment(next->seed, config);
                    cur.label_columns = col;
                    assert(cur.is_valid(graph, &config));
                    continue_callback(std::move(cur));
                    return;
                }

                assert(coord_dist);
                std::string_view query = !next->get_orientation() ? forward : reverse;

                Vector<Alignment::Column> inter_columns;
                Alignment::CoordinateSet inter_coords;
                const auto &last_coords = coords[last_anchor->tuple_id];
                const auto &next_coords = coords[next->tuple_id];
                {
                    auto it = last_coords.begin();
                    auto jt = next_coords.begin();
                    while (it != last_coords.end() && jt != next_coords.end()) {
                        if (it->first < jt->first) {
                            ++it;
                        } else if (it->first > jt->first) {
                            ++jt;
                        } else {
                            inter_columns.emplace_back(it->first);
                            auto &coords = inter_coords.emplace_back();
                            auto iit = it->second.begin();
                            auto jjt = jt->second.begin();
                            while (iit != it->second.end() && jjt != jt->second.end()) {
                                if (*iit + coord_dist < *jjt) {
                                    ++iit;
                                } else if (*iit + coord_dist > *jjt) {
                                    ++jjt;
                                } else {
                                    coords.emplace_back(*iit);
                                    ++iit;
                                    ++jjt;
                                }
                            }

                            if (coords.empty()) {
                                inter_columns.pop_back();
                                inter_coords.pop_back();
                            }

                            ++it;
                            ++jt;
                        }
                    }
                }
                assert(inter_columns.size());

                Vector<Alignment::Column> start_columns;
                start_columns.reserve(inter_columns.size());
                {
                    const auto &cur_columns = cur.get_columns();
                    ssize_t coord_offset = graph.get_k() - cur.get_sequence().size();
                    cur.label_coordinates.clear();
                    auto it = cur_columns.begin();
                    auto jt = inter_columns.begin();
                    auto jjt = inter_coords.begin();
                    while (it != cur_columns.end() && jt != inter_columns.end()) {
                        if (*it < *jt) {
                            ++it;
                        } else if (*it > *jt) {
                            ++jt;
                            ++jjt;
                        } else {
                            start_columns.emplace_back(*it);
                            auto &coords = cur.label_coordinates.emplace_back();
                            coords.reserve(jjt->size());
                            for (int64_t coord : *jjt) {
                                coords.emplace_back(coord + coord_offset);
                            }
                            ++it;
                            ++jt;
                            ++jjt;
                        }
                    }
                }

                if (start_columns.empty()) {
                    cur.trim_offset();
                    ++seed_count;
                    assert(cur.is_valid(graph, &config));
                    callback(std::move(cur));
                    return;
                }

                cur.set_columns(std::move(start_columns));
                auto [first, seed] = split_seed(graph, config, cur);

                auto extender = aligner.build_extender(query);

                auto extensions = extender->get_extensions(
                    seed,
                    config.ninf,
                    true,
                    coord_dist + seed.get_sequence().size(),
                    next->seed.get_nodes().back(),
                    false,
                    next->get_end_clipping(),
                    std::max(config.xdrop, 100) - config.xdrop
                );

                ++num_extensions;
                num_explored_nodes += extender->num_explored_nodes();

                if (extensions.empty()) {
                    cur.trim_offset();
                    ++seed_count;
                    assert(cur.is_valid(graph, &config));
                    callback(std::move(cur));
                } else {
                    for (auto&& ext : extensions) {
                        if (ext.get_end_clipping() == next->get_end_clipping()
                                && ext.get_clipping() == seed.get_clipping()
                                && ext.get_nodes().back() == next->seed.get_nodes().back()) {
                            next->used = true;
                            if (first.size()) {
                                Alignment cur_first = first;
                                cur_first.label_columns = ext.label_columns;
                                assert(cur_first.is_valid(graph, &config));
                                cur_first.splice(std::move(ext));
                                assert(cur_first.size());
                                assert(cur_first.is_valid(graph, &config));
                                std::swap(cur_first, ext);
                            }

                            ext.label_coordinates.clear();
                            continue_callback(std::move(ext));
                        }
                    }
                }
            },
            [&](Alignment&& aln) {
                aln.trim_offset();
                ++seed_count;
                assert(aln.is_valid(graph, &config));
                callback(std::move(aln));
            }
        );
    }

    logger->trace("Reduced {} seeds down to {}",
                  fwd_seeds.size() + bwd_seeds.size(),
                  seed_count);

    if (!config.ignore_discarded_seeds) {
        for (auto &[col, anchors] : clustered_seeds) {
            for (auto &anchor : anchors) {
                if (!anchor.used)
                    discarded_seed_callback(std::move(anchor.seed));
            }
        }
    }

    return std::make_pair(num_extensions, num_explored_nodes);
}


Alignment::score_t get_label_change_score(const AnnotationBuffer *anno_buffer,
                                          Alignment::Column col_a,
                                          Alignment::Column col_b) {
    if (col_a == col_b)
        return 0;

    assert(anno_buffer);
    return anno_buffer->get_label_change_score(col_a, col_b);
};

void chain_alignments(const IDBGAligner &aligner,
                      std::vector<Alignment>&& alignments,
                      const std::function<bool(Alignment::Column, size_t, score_t)> &start_backtrack,
                      const std::function<void(Alignment&&)> &callback,
                      const std::function<bool()> &terminate) {
    if (terminate())
        return;

    const auto &config = aligner.get_config();

    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        return a.get_orientation() < b.get_orientation();
    });

    if (alignments.size() <= 1
            || (alignments.size() == 2
                && alignments[1].get_orientation() != alignments[0].get_orientation())) {
        return;
    }

    const DeBruijnGraph &graph = aligner.get_graph();
    std::string_view query = alignments[0].get_full_query_view();

    struct Anchor {
        std::string_view::const_iterator end;
        std::string_view::const_iterator begin;
        size_t index;
        size_t spelling_length;
        bool orientation;
        size_t clipping;
        size_t end_clipping;
        ssize_t node_idx;
        score_t score;
        Alignment::Column col;

        std::string_view get_query_view() const {
            return std::string_view(begin, end - begin);
        }

        bool get_orientation() const { return orientation; }

        size_t get_clipping() const { return clipping; }
        size_t get_end_clipping() const { return end_clipping; }

        score_t get_score(const DBGAlignerConfig&) const { return score; }
    };

    size_t seed_size = std::min(config.min_seed_length, graph.get_k());

    // preprocess alignments
    size_t orientation_change = 0;
    std::vector<Anchor> anchors;

    for (size_t i = 0; i < alignments.size(); ++i) {
        const auto &alignment = alignments[i];
        bool is_fwd_orientation = !alignment.get_orientation();
        DEBUG_LOG("Alignment {}: {}\t{}\t{}",
                  i, alignment.get_query_view(), alignment.get_nodes().size(), alignment);
        auto full_query_begin = alignment.get_full_query_view().begin();

        auto cur = alignment;
        while (cur.size() > 1 && std::min(cur.get_query_view().size(), cur.get_sequence().size()) >= seed_size) {
            auto it = cur.get_cigar().data().rbegin();
            if (it->first == Cigar::CLIPPED)
                ++it;

            switch (it->first) {
                case Cigar::INSERTION:
                case Cigar::MISMATCH: {
                    if (it->first == Cigar::MISMATCH && it->second >= cur.size()) {
                        cur.trim_query_suffix(cur.size() - 1, config, false);
                    } else {
                        cur.trim_query_suffix(it->second, config, false);
                    }
                    assert(cur.size());
                } break;
                case Cigar::DELETION:
                case Cigar::MATCH: {
                    size_t num_to_trim = it->second;
                    if (it->first == Cigar::MATCH && it->second >= seed_size) {
                        for (size_t num_matches = it->second; num_matches >= seed_size && cur.size() > 1; --num_matches, --num_to_trim, cur.trim_query_suffix(1, config, false)) {
                            DEBUG_LOG("Anchor from: {}\t{}\t{}", i, cur.get_nodes().size() - 1, cur);
                            orientation_change += is_fwd_orientation;

                            auto end = cur.get_query_view().end();
                            assert(anchors.empty() || anchors.back().index != i || anchors.back().end > end);

                            size_t spelling_length = cur.get_sequence().size();
                            assert(anchors.empty() || anchors.back().index != i || anchors.back().spelling_length > spelling_length);

                            const auto &a_i = anchors.emplace_back(Anchor{
                                .end = end,
                                .begin = end - seed_size,
                                .index = i,
                                .spelling_length = spelling_length,
                                .orientation = alignment.get_orientation(),
                                .clipping = static_cast<size_t>(end - full_query_begin) - seed_size,
                                .end_clipping = cur.get_end_clipping(),
                                .node_idx = static_cast<ssize_t>(cur.get_nodes().size()) - 1,
                                .score = cur.get_score(),
                                .col = std::numeric_limits<Alignment::Column>::max()
                            });
                            std::ignore = a_i;

                            assert(cur.get_nodes().back() == alignment.get_nodes()[a_i.node_idx]);
                            assert(graph.get_node_sequence(cur.get_nodes().back()).substr(graph.get_k() - seed_size)
                                    == std::string_view(a_i.begin, a_i.end - a_i.begin));
                        }
                    }

                    if (num_to_trim >= cur.size()) {
                        cur.trim_reference_suffix(cur.size() - 1, config, false);
                    } else {
                        cur.trim_reference_suffix(num_to_trim, config, false);
                    }

                    assert(cur.size());
                } break;
                case Cigar::NODE_INSERTION:
                case Cigar::CLIPPED: {
                    assert(false);
                }
            }
        }

        assert(cur.size() == 1);
        ssize_t spelling_size = cur.get_sequence().size();
        cur.extend_offset(std::vector<node_index>(spelling_size - cur.size(),
                                                  DeBruijnGraph::npos));
        assert(cur.get_offset() == graph.get_k() - 1);

        for ( ; std::min(cur.get_query_view().size(), cur.get_sequence().size()) >= seed_size; cur.trim_query_suffix(1, config)) {
            auto it = cur.get_cigar().data().rbegin();
            if (it->first == Cigar::CLIPPED)
                ++it;

            if (it->first == Cigar::MATCH && it->second >= seed_size) {
                ssize_t node_idx = static_cast<ssize_t>(cur.get_sequence().size()) - spelling_size;
                DEBUG_LOG("Anchor from: {}\t{}\t{}", i, node_idx, cur);
                orientation_change += is_fwd_orientation;

                auto end = cur.get_query_view().end();
                assert(anchors.empty() || anchors.back().index != i || anchors.back().end > end);

                size_t spelling_length = cur.get_sequence().size();
                assert(anchors.empty() || anchors.back().index != i || anchors.back().spelling_length > spelling_length);

                anchors.emplace_back(Anchor{
                    .end = end,
                    .begin = end - seed_size,
                    .index = i,
                    .spelling_length = spelling_length,
                    .orientation = alignment.get_orientation(),
                    .clipping = static_cast<size_t>(end - full_query_begin) - seed_size,
                    .end_clipping = cur.get_end_clipping(),
                    .node_idx = node_idx,
                    .score = cur.get_score(),
                    .col = std::numeric_limits<Alignment::Column>::max()
                });
            }
        }
    }

    auto preprocess_range = [&](auto begin, auto end) {
        if (begin == end)
            return;

        std::sort(begin, end, [](const Anchor &a, const Anchor &b) { return a.col < b.col; });

        std::vector<tsl::hopscotch_set<size_t>> end_counters(query.size() + 1);

        auto last_it = begin;
        while (last_it != end) {
            auto it = std::find_if(last_it, end, [&](const Anchor &a) {
                return a.col != last_it->col;
            });

            for (auto &c : end_counters) {
                c.clear();
            }

            std::for_each(last_it, it, [&](const Anchor &a) {
                end_counters[a.end_clipping + 1].emplace(a.index);
            });

            std::for_each(last_it, it, [&](Anchor &a) {
                if (end_counters[a.end_clipping + 1].size() == 1
                        && end_counters[a.end_clipping].count(a.index)) {
                    a.index = std::numeric_limits<size_t>::max();
                }
            });

            last_it = it;
        }
    };

    preprocess_range(anchors.begin(), anchors.begin() + orientation_change);
    preprocess_range(anchors.begin() + orientation_change, anchors.end());

    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    AnnotationBuffer *anno_buffer = nullptr;
    bool allow_label_change = false;

    if (labeled_aligner) {
        anno_buffer = &labeled_aligner->get_annotation_buffer();
        allow_label_change = anno_buffer->allow_label_change();

        std::vector<Anchor> split_anchors;
        for (auto &a : anchors) {
            if (a.index != std::numeric_limits<size_t>::max()) {
                assert(alignments[a.index].label_columns);
                assert(alignments[a.index].label_column_diffs.empty());
                for (auto c : alignments[a.index].get_columns()) {
                    assert(c != std::numeric_limits<Alignment::Column>::max());
                    split_anchors.emplace_back(a);
                    split_anchors.back().col = c;
                }
            }
        }
        std::swap(split_anchors, anchors);

    } else {
        anchors.erase(std::remove_if(anchors.begin(), anchors.end(),
                                     [](const auto &a) {
                                         return a.index == std::numeric_limits<uint64_t>::max();
                                     }),
                      anchors.end());
    }

    if (!allow_label_change) {
        std::sort(anchors.begin(), anchors.end(), [](const auto &a, const auto &b) {
            return std::tie(a.col, a.orientation, a.end) < std::tie(b.col, b.orientation, b.end);
        });
    } else {
        std::sort(anchors.begin(), anchors.end(), [](const auto &a, const auto &b) {
            return std::tie(a.orientation, a.end) < std::tie(b.orientation, b.end);
        });
    }

    // construct index to seed map
    using EndMap = tsl::hopscotch_map<std::string_view::const_iterator, const Anchor*>;
    std::vector<tsl::hopscotch_map<Alignment::Column, EndMap>> coord_map(alignments.size());
    for (const auto &a : anchors) {
        auto &bucket = coord_map[a.index][a.col];
        assert(bucket.find(a.end) == bucket.end());
        bucket[a.end] = &a;
    }

    logger->trace("Chaining alignments using {} anchors for a query of length {}",
                  anchors.size(), query.size());
    DEBUG_LOG("\tAllowing label changes: {}", allow_label_change);

    score_t node_insert = config.node_insertion_penalty;
    score_t gap_open = config.gap_opening_penalty;
    score_t gap_ext = config.gap_extension_penalty;
    assert(gap_open < 0);
    assert(gap_ext < 0);
    assert(gap_ext >= gap_open);
    assert(node_insert < 0);

    auto get_overlapping_prev_bucket = [&](const EndMap &bucket, const Anchor &a_i) -> const Anchor* {
        // when we want to connect a_i -> a_j
        // check if there is a previous seed a_last from in a_j.index s.t. a_last.end == a_i.end
        auto find_last_col = bucket.find(a_i.end);
        return find_last_col != bucket.end() ? find_last_col->second : nullptr;
    };

    auto get_overlapping_prev = [&](const Anchor &a_j, const Anchor &a_i) -> const Anchor* {
        // when we want to connect a_i -> a_j
        // check if there is a previous seed a_last from in a_j.index s.t. a_last.end == a_i.end
        assert(coord_map[a_j.index].count(a_j.col));
        return get_overlapping_prev_bucket(coord_map[a_j.index].find(a_j.col)->second, a_i);
    };

    auto last_anchor_it = anchors.data();
    while (!terminate() && last_anchor_it != anchors.data() + anchors.size()) {
        auto anchor_it = anchors.data() + anchors.size();

        if (!allow_label_change) {
            anchor_it = std::find_if(last_anchor_it, anchor_it, [&](const Anchor &a) {
                return a.col != last_anchor_it->col;
            });
        }

        score_t chain_score = 0;
        AnchorChain<Anchor> last_chain(0);
        Alignment::Columns col_idx = 0;
        score_t full_score = 0;

        chain_anchors<Anchor>(config, last_anchor_it, anchor_it,
            [&](const Anchor &a_j,
                ssize_t,
                const Anchor *begin,
                const Anchor *end,
                auto *chain_scores,
                const auto &update_score) {

                const auto &[score_j, last_j, last_dist_j] = *(chain_scores + (end - begin));
                const Alignment &full_j = alignments[a_j.index];
                std::string_view full_query_j = full_j.get_query_view();
                std::string_view query_j(a_j.begin, a_j.end - a_j.begin);
                const EndMap *bucket = nullptr;

                --chain_scores;

                std::for_each(begin, end, [&](const Anchor &a_i) {
                    // try to connect a_i -> a_j
                    ++chain_scores;
                    assert(a_i.end <= a_j.end);
                    if (a_i.end == a_j.end)
                        return;

                    auto [score_i, last_i, last_dist_i] = *chain_scores;
                    if (last_i == anchor_it) {
                        assert(!last_dist_i);
                        last_dist_i = a_i.spelling_length;
                    }

                    if (a_i.index == a_j.index) {
                        assert(a_j.spelling_length > a_i.spelling_length);
                        size_t added_length = a_j.spelling_length - a_i.spelling_length;
                        update_score(score_i + a_j.score - a_i.score, &a_i,
                                     last_dist_i + added_length);
                        return;
                    }

                    if (last_dist_i < graph.get_k())
                        return;

                    assert(a_i.node_idx >= 0);

                    const Alignment &full_i = alignments[a_i.index];
                    std::string_view full_query_i = full_i.get_query_view();
                    std::string_view query_i(a_i.begin, a_i.end - a_i.begin);
                    if (!bucket)
                        bucket = &coord_map[a_j.index].find(a_j.col)->second;

                    if (const Anchor *a_last = get_overlapping_prev_bucket(*bucket, a_i)) {
                        // we want to connect a_i -> a_last(on j) -> a_j
                        assert(a_last != &a_j);
                        assert(a_last->index == a_j.index);
                        assert(a_last->spelling_length < a_j.spelling_length);
                        assert(a_last->end == a_i.end);

                        // previous score, plus score of a_last -> a_j
                        score_t updated_score = score_i + a_j.score - a_last->score;

                        // calculate score of jumping from a_i -> a_last
                        // if a_i and a_last are not the same node, add a node insertion
                        size_t overlap = seed_size;

                        if (a_last->node_idx < 0
                                || full_j.get_nodes()[a_last->node_idx] != full_i.get_nodes()[a_i.node_idx]) {
                            auto rbegin_i = full_i.get_sequence().rbegin() + (full_i.get_sequence().size() - a_i.spelling_length);
                            auto rend_i = rbegin_i + std::min(graph.get_k() - 1, static_cast<size_t>(full_i.get_sequence().rend() - rbegin_i));

                            auto rbegin_last = full_j.get_sequence().rbegin() + (full_j.get_sequence().size() - a_last->spelling_length);
                            auto rend_last = rbegin_last + std::min(graph.get_k() - 1, static_cast<size_t>(full_j.get_sequence().rend() - rbegin_last));

                            overlap = std::mismatch(rbegin_i, rend_i, rbegin_last, rend_last).first - rbegin_i;

                            // if there is not enough sequence left in the alignment, we can't connect them
                            if (full_j.get_sequence().size() - a_last->spelling_length + overlap < graph.get_k())
                                return;

                            updated_score += node_insert;
                        } else {
                            overlap = graph.get_k();
                        }

                        if (updated_score > score_j) {
                            updated_score += get_label_change_score(anno_buffer, a_i.col, a_last->col);
                            update_score(updated_score, &a_i,
                                        overlap + (a_j.spelling_length - a_last->spelling_length));
                        }
                    } else if (full_query_i.end() <= full_query_j.begin()) {
                        // completely disjoint and a_i is at the end of full_i
                        score_t gap = full_query_j.begin() - full_query_i.end();
                        score_t gap_cost = node_insert + gap_open;
                        if (gap > 0)
                            gap_cost += gap_open + (gap - 1) * gap_ext;

                        assert(gap_cost < 0);

                        score_t updated_score = score_i + full_i.get_score() - a_i.score + gap_cost + a_j.score;

                        if (updated_score > score_j) {
                            updated_score += get_label_change_score(anno_buffer, a_i.col, a_j.col);
                            update_score(updated_score, &a_i, a_j.spelling_length);
                        }
                    }
                });
            },
            [&](const AnchorChain<Anchor> &chain, const std::vector<score_t> &score_traceback) {
                assert(chain.size());
                assert(score_traceback.size());

#ifndef NDEBUG
                DEBUG_LOG("Chain: {}\t{}", chain.size(), chain.get_score());
                for (const auto &[anchor, dist] : chain) {
                    DEBUG_LOG("\t{}\t{}\t{}S{}={}S\t{}", dist,
                        anchor->index,
                        anchor->clipping,anchor->end - anchor->begin,anchor->end_clipping,
                        std::string_view(anchor->begin, anchor->end-anchor->begin)
                    );
                }
#endif

                const auto &last_anchor = *chain.back().first;
                const auto &last_aln = alignments[last_anchor.index];

                if (chain.back().second + last_aln.get_sequence().size() - last_anchor.spelling_length < graph.get_k())
                    return false;

                score_t score = score_traceback.back();
                if (chain_score == score && std::equal(chain.begin(), chain.end(),
                                                       last_chain.begin(), last_chain.end(),
                                                       [](const auto &a, const auto &b) {
                                                           return a.first->index == b.first->index
                                                                    && a.first->col == b.first->col;
                                                       })) {
                    return false;
                }

                if (chain.size() > 1 && std::all_of(chain.begin() + 1, chain.end(),
                        [&](const auto &a) {
                            return a.first->index == chain.front().first->index;
                        })) {
                    return false;
                }

                full_score = score + last_aln.get_score() - last_anchor.get_score(config);

                size_t aln_size = 0;
                for (const auto &[ptr, d] : chain) {
                    aln_size += d;
                }

                if (start_backtrack(chain[0].first->col, aln_size, full_score)) {
                    last_chain = chain;
                    chain_score = score;
                    if (labeled_aligner) {
                        assert(anno_buffer);
                        col_idx = anno_buffer->cache_column_set(1, chain[0].first->col);
                    }

                    return true;
                } else {
                    return false;
                }
            },
            true /* extend_anchors */,
            [&](const Anchor *next,
                const Anchor *last_anchor,
                Alignment&& cur,
                size_t dist,
                score_t score_up_to_now,
                const auto &callback) {

                if (last_anchor && next->col != last_anchor->col) {
                    assert(anno_buffer);
                    DEBUG_LOG("\tSwitching {} -> {}",
                              anno_buffer->get_annotator().get_label_encoder().decode(last_anchor->col),
                              anno_buffer->get_annotator().get_label_encoder().decode(next->col));
                    col_idx = anno_buffer->cache_column_set(1, next->col);
                }

#ifndef NDEBUG
                auto check_aln = [&](Alignment aln, score_t label_change_score = DBGAlignerConfig::ninf) {
                    assert(aln.size());
                    assert(next->end <= aln.get_query_view().end());
                    aln.trim_query_suffix(aln.get_query_view().end() - next->end, config);
                    DEBUG_LOG("\tScore to now: {}\tChain: {}\tNode insertion penalty: {}\tLabel change score: {}",
                              score_up_to_now, aln, node_insert, label_change_score);
                    return aln.get_score() == score_up_to_now;
                };
#endif

                if (cur.empty()) {
                    assert(!last_anchor);
                    Alignment alignment = alignments[next->index];
                    alignment.label_columns = col_idx;

                    size_t added = alignment.get_sequence().size() - alignment.size();
                    alignment.extend_offset(std::vector<node_index>(added, DeBruijnGraph::npos));

                    DEBUG_LOG("\tStarting: {}\tfrom {}", alignment, alignments[next->index]);
                    assert(check_aln(alignment));
                    callback(std::move(alignment));
                    return;
                }

                assert(dist);
                assert(last_anchor);
                if (next->index == last_anchor->index) {
                    assert(next->spelling_length > last_anchor->spelling_length);
                    assert(next->col == last_anchor->col);
                    assert(check_aln(cur));
                    callback(std::move(cur));
                    return;
                }

                std::ignore = dist;

                Alignment alignment = alignments[next->index];
                alignment.label_columns = col_idx;

                DEBUG_LOG("\t\tcur: {}", cur);
                DEBUG_LOG("\t\tnxt: {}", alignment);

                if (const Anchor *a_o = get_overlapping_prev(*next, *last_anchor)) {
                    assert(dist > seed_size);
                    assert(last_anchor->end == a_o->end);
                    assert(next->col == a_o->col);

                    cur.trim_query_suffix(cur.get_query_view().end() - a_o->end, config);
                    assert(cur.size());
                    assert(cur.is_valid(graph, &config));

                    alignment.trim_query_prefix(a_o->end - alignment.get_query_view().begin(),
                                                graph.get_k() - 1, config, false);
                    assert(alignment.size());
                    alignment.extend_offset(std::vector<node_index>(
                        alignment.get_sequence().size() - alignment.size(),
                        DeBruijnGraph::npos
                    ));
                    assert(alignment.is_valid(graph, &config));

                    if (a_o->node_idx < 0 || cur.get_nodes().back() != alignments[a_o->index].get_nodes()[a_o->node_idx]) {
#ifndef NDEBUG
                        auto rbegin_i = alignments[last_anchor->index].get_sequence().rbegin() + (alignments[last_anchor->index].get_sequence().size() - last_anchor->spelling_length);
                        auto rend_i = rbegin_i + std::min(graph.get_k() - 1, static_cast<size_t>(alignments[last_anchor->index].get_sequence().rend() - rbegin_i));

                        auto rbegin_last = alignments[next->index].get_sequence().rbegin() + (alignments[next->index].get_sequence().size() - a_o->spelling_length);
                        auto rend_last = rbegin_last + std::min(graph.get_k() - 1, static_cast<size_t>(alignments[next->index].get_sequence().rend() - rbegin_last));
                        ssize_t test_overlap = std::mismatch(rbegin_i, rend_i, rbegin_last, rend_last).first - rbegin_i;
                        assert(dist == test_overlap + (next->spelling_length - a_o->spelling_length));
#endif
                        // dist == overlap + (next->spelling_length - a_o->spelling_length)
                        alignment.insert_gap_prefix(
                            static_cast<ssize_t>(next->spelling_length - a_o->spelling_length) - dist,
                            graph.get_k() - 1, config
                        );
                        assert(alignment.size());
                    }
                } else {
                    // no overlap
                    assert(cur.get_query_view().end() <= alignment.get_query_view().begin());
                    assert(dist == next->spelling_length);
                    assert(cur.get_query_view().end() == alignments[last_anchor->index].get_query_view().end());

                    alignment.insert_gap_prefix(
                        alignment.get_query_view().begin() - cur.get_query_view().end(),
                        graph.get_k() - 1, config
                    );
                    assert(alignment.size());
                }

                DEBUG_LOG("\t\tA: {}", cur);
                DEBUG_LOG("\t\tB: {}", alignment);

                score_t label_change_score = DBGAlignerConfig::ninf;
                if (next->col != last_anchor->col) {
                    assert(anno_buffer);
                    assert(allow_label_change);
                    label_change_score = get_label_change_score(anno_buffer, last_anchor->col, next->col);
                    DEBUG_LOG("\t\t\tLabel change: {} ({}) -> {} ({})\t{}",
                        last_anchor->col, anno_buffer->get_annotator().get_label_encoder().decode(last_anchor->col),
                        next->col, anno_buffer->get_annotator().get_label_encoder().decode(next->col),
                        label_change_score);
                }

                assert(next->col == last_anchor->col ||
                        label_change_score != DBGAlignerConfig::ninf);

                cur.splice(std::move(alignment), label_change_score);
                DEBUG_LOG("\tCurrent: {}", cur);
                assert(cur.size());
                assert(cur.is_valid(graph, &config));
                assert(cur.get_end_clipping() == alignments[next->index].get_end_clipping());
                assert(check_aln(cur, label_change_score));
                callback(std::move(cur));
            },
            [&](Alignment&& aln) {
                aln.trim_offset();
                DEBUG_LOG("\tFinal: {}\t{}", chain_score, aln);
                assert(aln.size());
                assert(aln.get_score() == full_score);
                callback(std::move(aln));
            },
            terminate
        );

        last_anchor_it = anchor_it;
    }
}

} // namespace align
} // namespace graph
} // namespace mtg
