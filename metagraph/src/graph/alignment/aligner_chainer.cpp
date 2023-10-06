#include "aligner_chainer.hpp"

#include <unordered_set>

#include <x86/svml.h>

#include "aligner_seeder_methods.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"
#include "chainer.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"

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
        uint64_t index;
        size_t spelling_length;
        bool orientation;
        uint64_t clipping;
        uint64_t end_clipping;
        int64_t node_idx;
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
    std::vector<std::vector<score_t>> per_char_scores_prefix;
    std::vector<std::vector<score_t>> per_char_scores_prefix_del;
    per_char_scores_prefix.reserve(alignments.size());
    per_char_scores_prefix_del.reserve(alignments.size());

    for (size_t i = 0; i < alignments.size(); ++i) {
        const auto &alignment = alignments[i];
        bool is_fwd_orientation = !alignment.get_orientation();
        DEBUG_LOG("Alignment {}: {}\t{}\t{}",
                  i, alignment.get_query_view(), alignment.get_nodes().size(), alignment);
        std::string_view query = alignment.get_query_view();

        auto &prefix_scores_with_deletions = per_char_scores_prefix.emplace_back();
        prefix_scores_with_deletions.reserve(query.size() + 1);
        auto &prefix_scores_without_deletions = per_char_scores_prefix_del.emplace_back();
        prefix_scores_without_deletions.reserve(query.size() + 1);

        ssize_t start_node_idx = static_cast<ssize_t>(alignment.get_offset())
                                    - graph.get_k() + seed_size;

        for (auto cur = alignment; cur.size(); cur.trim_query_prefix(1, graph.get_k() - 1, config, false)) {
            prefix_scores_without_deletions.emplace_back(cur.get_score());
            auto it = cur.get_cigar().data().begin();
            assert(it != cur.get_cigar().data().end());
            if (it->first == Cigar::CLIPPED) {
                ++it;
                assert(it != cur.get_cigar().data().end());
            }

            if (it->first == Cigar::DELETION) {
                cur.trim_reference_prefix(it->second, graph.get_k() - 1, config, false);
                it = cur.get_cigar().data().begin();
                assert(it != cur.get_cigar().data().end());
                if (it->first == Cigar::CLIPPED) {
                    ++it;
                    assert(it != cur.get_cigar().data().end());
                }
            }

            ssize_t node_idx = start_node_idx + alignment.get_sequence().size()
                                - cur.get_sequence().size();
            prefix_scores_with_deletions.emplace_back(cur.get_score());
            if (it->first == Cigar::MATCH && it->second >= seed_size) {
                orientation_change += is_fwd_orientation;
                DEBUG_LOG("Anchor from: {}\t{}", i, cur);
                anchors.emplace_back(Anchor{
                    .end = cur.get_query_view().begin() + seed_size,
                    .begin = cur.get_query_view().begin(),
                    .index = i,
                    .spelling_length = cur.get_sequence().size(),
                    .orientation = alignment.get_orientation(),
                    .clipping = cur.get_clipping(),
                    .end_clipping = alignment.get_full_query_view().end()
                                        - cur.get_query_view().begin() - seed_size,
                    .node_idx = node_idx,
                    .score = cur.get_score(),
                    .col = std::numeric_limits<Alignment::Column>::max()
                });

#ifndef NDEBUG
                const auto &a_i = anchors.back();
                if (a_i.node_idx >= 0) {
                    assert(static_cast<size_t>(a_i.node_idx) < alignment.size());
                    assert(graph.get_node_sequence(alignment.get_nodes()[a_i.node_idx]).substr(graph.get_k() - seed_size)
                            == std::string_view(a_i.begin, a_i.end - a_i.begin));
                }
#endif
            }
        }

        prefix_scores_with_deletions.emplace_back(0);
        prefix_scores_without_deletions.emplace_back(0);
    }

    auto preprocess_range = [&](auto begin, auto end) {
        if (begin == end)
            return;

        std::sort(begin, end, [](const Anchor &a, const Anchor &b) {
            return std::tie(b.col, a.end, a.begin) > std::tie(a.col, b.end, b.begin);
        });

        auto last_it = begin;
        std::vector<tsl::hopscotch_set<size_t>> end_counters(query.size() + 1);

        while (last_it != end) {
            auto it = last_it + 1;
            while (it != end && it->col == last_it->col) {
                ++it;
            }

            for (auto &c : end_counters) {
                c.clear();
            }

            std::for_each(last_it, it, [&](const Anchor &a) {
                end_counters[a.end_clipping].emplace(a.index);
            });

            std::for_each(last_it, it, [&](Anchor &a) {
                if (end_counters[a.end_clipping].size() == 1
                        && end_counters[a.end_clipping + 1].count(a.index)) {
                    a.index = std::numeric_limits<uint64_t>::max();
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
        if (allow_label_change)
            logger->trace("\tAllowing label changes");

        std::vector<Anchor> split_anchors;
        for (auto &a : anchors) {
            if (a.index != std::numeric_limits<uint64_t>::max()) {
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
            return std::tie(b.col, b.orientation, a.end) > std::tie(a.col, a.orientation, b.end);
        });
    } else {
        std::sort(anchors.begin(), anchors.end(), [](const auto &a, const auto &b) {
            return std::tie(b.orientation, a.end) > std::tie(a.orientation, b.end);
        });
    }

    // construct index to seed map
    std::vector<tsl::hopscotch_map<std::string_view::const_iterator,
                                   tsl::hopscotch_map<Alignment::Column, size_t>>> coord_map(alignments.size());
    for (size_t i = 0; i < anchors.size(); ++i) {
        const auto &a = anchors[i];
        auto &bucket = coord_map[a.index][a.end];
        assert(bucket.find(a.col) == bucket.end());
        bucket[a.col] = i;
    }

    logger->trace("Chaining alignments using {} anchors for a query of length {}",
                  anchors.size(), query.size());

    score_t node_insert = config.node_insertion_penalty;
    score_t gap_open = config.gap_opening_penalty;
    score_t gap_ext = config.gap_extension_penalty;
    assert(gap_open < 0);
    assert(gap_ext < 0);
    assert(gap_ext >= gap_open);
    assert(node_insert < 0);

    auto get_overlapping_prev = [&](const Anchor &a_i, const Anchor &a_j) -> size_t {
        // check if there is a previous seed a_last from in a_i.index s.t. a_last.end == a_j.end
        auto find_last = coord_map[a_i.index].find(a_j.end);
        if (find_last == coord_map[a_i.index].end())
            return std::numeric_limits<size_t>::max();

        auto find_last_col = find_last->second.find(a_i.col);

        if (find_last_col == find_last->second.end())
            return std::numeric_limits<size_t>::max();

        size_t last_overlap_i = find_last_col->second;
        const Anchor &a_last = anchors[last_overlap_i];

        assert(a_last.index == a_i.index);
        assert(a_last.end == a_j.end);

        if (a_last.node_idx < 0 || a_i.spelling_length <= a_last.spelling_length)
            return std::numeric_limits<size_t>::max();

        return last_overlap_i;
    };

    auto last_anchor_it = anchors.data();
    while (!terminate() && last_anchor_it != anchors.data() + anchors.size()) {
        auto anchor_it = last_anchor_it + 1;
        while (anchor_it != anchors.data() + anchors.size()
                && (allow_label_change || anchor_it->col == last_anchor_it->col)
                && anchor_it->orientation == last_anchor_it->orientation) {
            ++anchor_it;
        }

        const Anchor *last_anchor;
        score_t chain_score = 0;
        AnchorChain<Anchor> last_chain;
        Alignment::Columns col_idx = 0;
        score_t full_score = 0;

        chain_anchors<Anchor>(config, last_anchor_it, anchor_it,
            [&](const Anchor &a_i, ssize_t, const Anchor *begin, const Anchor *end, auto *chain_scores, const auto &update_score) {
                const auto *chain_scores_i = chain_scores - (begin - last_anchor_it) + (&a_i - last_anchor_it);

                const Alignment &full_i = alignments[a_i.index];
                std::string_view full_query_i = full_i.get_query_view();
                std::string_view query_i(a_i.begin, a_i.end - a_i.begin);

                --chain_scores;
                std::for_each(begin, end, [&](const Anchor &a_j) {
                    // try to connect a_i -> a_j
                    assert(&a_i != &a_j);
                    ++chain_scores;

                    const Alignment &full_j = alignments[a_j.index];
                    std::string_view full_query_j = full_j.get_query_view();
                    std::string_view query_j(a_j.begin, a_j.end - a_j.begin);

                    if (query_j.end() == query_i.end())
                        return;

                    auto [score_j, last, last_dist] = *chain_scores;
                    if (last == anchor_it) {
                        assert(last_dist == std::numeric_limits<size_t>::max());
                        last_dist = -a_j.spelling_length;
                    }

                    if (a_i.col != a_j.col && (!allow_label_change || a_i.index == a_j.index))
                        return;

                    if (a_i.index == a_j.index) {
                        if (a_i.spelling_length <= a_j.spelling_length) {
                            // TODO: is this test still valid?
                            // assert(a_i.col != a_j.col);
                            return;
                        }

                        size_t added_length = a_i.spelling_length - a_j.spelling_length;
                        update_score(score_j + a_i.score - a_j.score,
                                     &a_j, last_dist - added_length);
                        return;
                    }

                    if (-last_dist < graph.get_k())
                        return;

                    const auto &[score_i, last_j, last_dist_i] = *chain_scores_i;
                    size_t last_overlap_i = get_overlapping_prev(a_i, a_j);
                    if (last_overlap_i != std::numeric_limits<size_t>::max()) {
                        const Anchor &a_last = anchors[last_overlap_i];

                        // a_i -> a_last score
                        score_t a_i_updated_score = a_i.score - a_last.score;

                        // compute rest of the score
                        score_t score_seed_j = a_j.score
                            - per_char_scores_prefix_del[a_j.index][a_j.end - full_j.get_query_view().begin()];

                        const Alignment &full_last = alignments[a_last.index];
                        score_t score_seed_last = a_last.score
                            - per_char_scores_prefix_del[a_last.index][a_last.end - full_last.get_query_view().begin()];

                        score_t updated_score = score_j - score_seed_j + score_seed_last + a_i_updated_score;

                        if (a_j.node_idx < 0 || full_last.get_nodes()[a_last.node_idx] != full_j.get_nodes()[a_j.node_idx]) {
                            // nodes are different, add a node insertion penalty
                            updated_score += node_insert;
                        }

                        if (updated_score < score_i)
                            return;

                        updated_score += get_label_change_score(anno_buffer, a_last.col, a_j.col);

                        update_score(updated_score, &a_j,
                                     -seed_size - (a_i.spelling_length - a_last.spelling_length));
                    } else if (full_query_i.end() <= full_query_j.begin()
                            && a_j.clipping == full_j.get_clipping()
                            && a_i.spelling_length >= graph.get_k()) {
                        // completely disjoint
                        score_t gap = full_query_j.begin() - full_query_i.end();
                        score_t gap_cost = node_insert + gap_open;
                        if (gap > 0)
                            gap_cost += gap_open + (gap - 1) * gap_ext;

                        assert(gap_cost < 0);

                        score_t base_updated_score = score_j + gap_cost + a_i.score;

                        if (base_updated_score < score_i)
                            return;

                        base_updated_score += get_label_change_score(anno_buffer, a_i.col, a_j.col);
                        update_score(base_updated_score, &a_j, -a_i.spelling_length);
                    }
                });
            },
            [&](const AnchorChain<Anchor> &chain, score_t score) {
                assert(chain.size());
                if (chain_score == score && std::equal(chain.begin(), chain.end(),
                                                       last_chain.begin(), last_chain.end(),
                                                       [](const auto &a, const auto &b) {
                                                           return a.first->index == b.first->index
                                                                    && a.first->col == b.first->col;
                                                       })) {
                    return false;
                }

                if (chain.size() > 1) {
                    if (-chain[1].second < graph.get_k())
                        return false;

                    if (std::all_of(chain.begin() + 1, chain.end(),
                                    [&](const auto &a) {
                                        return a.first->index == chain.front().first->index;
                                    })) {
                        return false;
                    }
                }

                const auto &first_anchor = *chain.front().first;
                const auto &first_aln = alignments[first_anchor.index];
                full_score = score
                    + first_aln.get_score()
                    - per_char_scores_prefix[first_anchor.index][first_anchor.begin - first_aln.get_query_view().begin()];

                size_t aln_size = 0;
                for (const auto &[ptr, d] : chain) {
                    aln_size += -d;
                }

                if (start_backtrack(chain[0].first->col, aln_size, full_score)) {
                    last_chain = chain;
                    chain_score = score;
                    DEBUG_LOG("Chain: {}", score);
                    last_anchor = chain.back().first;
                    if (labeled_aligner) {
                        assert(anno_buffer);
                        col_idx = anno_buffer->cache_column_set(1, last_anchor->col);
                    }

                    return true;
                } else {
                    return false;
                }
            },
            true /* extend_anchors */,
            [&](const Anchor *first,
                Alignment&& cur,
                size_t dist,
                score_t score_up_to_now,
                const auto &callback) {

                Alignment alignment = alignments[first->index];
                if (first->col != last_anchor->col) {
                    assert(anno_buffer);
                    DEBUG_LOG("\tSwitching {} -> {}",
                              anno_buffer->get_annotator().get_label_encoder().decode(last_anchor->col),
                              anno_buffer->get_annotator().get_label_encoder().decode(first->col));
                    col_idx = anno_buffer->cache_column_set(1, first->col);
                }

                alignment.label_columns = col_idx;

                auto check_aln = [&](Alignment aln, score_t label_change_score = 0) {
#ifndef NDEBUG
                    assert(first->begin >= aln.get_query_view().begin());
                    aln.trim_query_prefix(first->begin - aln.get_query_view().begin(),
                                          graph.get_k() - 1,
                                          config);
                    DEBUG_LOG("\tScore to now: {}\tScore of chain: {}\tNode insertion penalty: {}\tLabel change score: {}",
                              score_up_to_now, aln.get_score(), node_insert, label_change_score);
                    assert(aln.get_score() == score_up_to_now);
#else
                    std::ignore = aln;
                    std::ignore = score_up_to_now;
                    std::ignore = label_change_score;
#endif
                };

                if (cur.empty()) {
                    assert(first == last_anchor);
                    DEBUG_LOG("\tStarting: {}\tfrom {}", alignment, alignments[first->index]);
                    check_aln(alignment);
                    callback(std::move(alignment));
                    return;
                }

                if (first->index == last_anchor->index) {
                    assert(first->col == last_anchor->col);
                    last_anchor = first;
                    check_aln(cur);
                    callback(std::move(cur));
                    return;
                }

                std::ignore = dist;
                assert(-dist >= graph.get_k());

                DEBUG_LOG("\t\taln: {}", alignment);
                DEBUG_LOG("\t\tcur: {}", cur);

                size_t last_overlap_i = get_overlapping_prev(*first, *last_anchor);
                if (last_overlap_i == std::numeric_limits<size_t>::max()) {
                    // no overlap
                    assert(alignment.get_query_view().end() <= cur.get_query_view().begin());
                    assert(dist == -first->spelling_length);
                    assert(last_anchor->begin == cur.get_query_view().begin());
                    cur.insert_gap_prefix(
                        cur.get_query_view().begin() - alignment.get_query_view().end(),
                        graph.get_k() - 1, config
                    );
                    assert(cur.size());
                } else {
                    assert(-dist > seed_size);
                    const Anchor &a_o = anchors[last_overlap_i];

                    assert(last_anchor->end == a_o.end);
                    assert(first->col == a_o.col);
                    alignment.extend_offset(std::vector<node_index>(graph.get_k() - 1 - alignment.get_offset(),
                                                                    DeBruijnGraph::npos),
                                            std::vector<size_t>(graph.get_k() - 1 - alignment.get_offset(), 0));
                    alignment.trim_query_suffix(alignment.get_query_view().end() - a_o.end,
                                                config);
                    assert(alignment.size());
                    assert(a_o.node_idx >= 0);
                    assert(alignment.get_nodes().back()
                        == alignments[a_o.index].get_nodes()[a_o.node_idx]);
                    // assert(alignment.is_valid(graph, &config));

                    cur.extend_offset(std::vector<node_index>(graph.get_k() - 1 - cur.get_offset(),
                                                              DeBruijnGraph::npos),
                                      std::vector<size_t>(graph.get_k() - 1 - alignment.get_offset(), 0));
                    cur.trim_query_prefix(a_o.end - cur.get_query_view().begin(),
                                          graph.get_k() - 1,
                                          config,
                                          false);
                    assert(cur.size());
                    assert(cur.is_valid(graph, &config));
                    node_index cur_front = last_anchor->node_idx >= 0
                        ? alignments[last_anchor->index].get_nodes()[last_anchor->node_idx]
                        : DeBruijnGraph::npos;

                    if (alignment.get_nodes().back() != cur_front) {
                        cur.insert_gap_prefix(-seed_size, graph.get_k() - 1, config);
                        assert(cur.size());
                    } else {
                        assert(cur_front);
                    }
                }

                DEBUG_LOG("\t\tA: {}", alignment);
                DEBUG_LOG("\t\tB: {}", cur);

                score_t label_change_score = DBGAlignerConfig::ninf;
                if (first->col != last_anchor->col) {
                    assert(anno_buffer);
                    assert(allow_label_change);
                    label_change_score = get_label_change_score(anno_buffer, first->col, last_anchor->col);
                    DEBUG_LOG("\t\t\tLabel change: {} ({}) -> {} ({})\t{}",
                        first->col, anno_buffer->get_annotator().get_label_encoder().decode(first->col),
                        last_anchor->col, anno_buffer->get_annotator().get_label_encoder().decode(last_anchor->col),
                        label_change_score);
                }

                assert(first->col == last_anchor->col ||
                        label_change_score != DBGAlignerConfig::ninf);

                last_anchor = first;

                alignment.splice(std::move(cur), label_change_score);
                DEBUG_LOG("\tCurrent: {}", alignment);
                assert(alignment.size());
                assert(alignment.is_valid(graph, &config));
                assert(alignment.get_clipping() == alignments[first->index].get_clipping());
                check_aln(alignment, label_change_score);
                callback(std::move(alignment));
            },
            [&](Alignment&& aln) {
                aln.trim_offset();

                DEBUG_LOG("\tFinal: {}\tfull_score: {}\t{}", chain_score, full_score, aln);
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
