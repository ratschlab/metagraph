#include "aligner_seeder_methods.hpp"

#include "sdust.h"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef Alignment::score_t score_t;
typedef Alignment::node_index node_index;
typedef boss::BOSS::edge_index edge_index;


#if ! _PROTEIN_GRAPH
inline bool is_low_complexity(std::string_view s, int T = 20, int W = 64) {
    int n;
    std::unique_ptr<uint64_t> r { sdust(0, (const uint8_t*)s.data(), s.size(), T, W, &n) };
    return n > 0;
}
#else
inline bool is_low_complexity(std::string_view, int = 20, int = 64) {
    // TODO: implement a checker here
    return false;
}
#endif

ExactSeeder::ExactSeeder(const DeBruijnGraph &graph,
                         std::string_view query,
                         bool orientation,
                         std::vector<node_index>&& nodes,
                         const DBGAlignerConfig &config)
      : graph_(graph),
        query_(query),
        orientation_(orientation),
        query_nodes_(std::move(nodes)),
        config_(config),
        num_matching_(num_exact_matching()) { assert(config_.check_config_scores()); }

size_t ExactSeeder::num_exact_matching() const {
    size_t num_matching = 0;
    size_t last_match_count = 0;
    for (auto it = query_nodes_.begin(); it != query_nodes_.end(); ++it) {
        if (*it) {
            auto jt = std::find(it + 1, query_nodes_.end(), node_index());
            num_matching += graph_.get_k() + std::distance(it, jt) - 1 - last_match_count;
            last_match_count = graph_.get_k();
            it = jt - 1;
        } else if (last_match_count) {
            --last_match_count;
        }
    }
    assert(num_matching <= query_.size());

    return num_matching;
}

auto ExactSeeder::get_seeds() const -> std::vector<Seed> {
    size_t k = graph_.get_k();
    assert(k >= config_.min_seed_length);

    if (num_matching_ < config_.min_exact_match * query_.size())
        return {};

    std::vector<Seed> seeds;

    if (config_.max_seed_length < k)
        return seeds;

    size_t end_clipping = query_.size() - k;
    for (size_t i = 0; i < query_nodes_.size(); ++i, --end_clipping) {
        if (query_nodes_[i] != DeBruijnGraph::npos) {
            assert(i + k <= query_.size());
            std::string_view query_window = query_.substr(i, k);
            if (config_.no_seed_complexity_filter || !is_low_complexity(query_window)) {
                seeds.emplace_back(query_window,
                                   std::vector<node_index>{ query_nodes_[i] },
                                   orientation_, 0, i, end_clipping);
            }
        }
    }

    return seeds;
}

void suffix_to_prefix(const DBGSuccinct &dbg_succ,
                      const std::pair<edge_index, edge_index> &index_range,
                      size_t length,
                      const std::function<void(node_index)> &callback) {
    const auto &boss = dbg_succ.get_boss();
    assert(length);
    assert(length < dbg_succ.get_k());

    auto call_range = [&](const std::pair<edge_index, edge_index> &range) {
#ifndef NDEBUG
        std::string suffix = boss.get_node_str(index_range.second).substr(boss.get_k() - length);
#endif
        for (edge_index i = range.first; i <= range.second; ++i) {
            assert(boss.get_node_str(i).substr(0, length) == suffix);
            if (node_index node = dbg_succ.boss_to_kmer_index(i)) {
                assert(dbg_succ.get_node_sequence(node).substr(0, length) == suffix);
                callback(node);
            }
        }
    };

    if (length == boss.get_k()) {
        call_range(index_range);
        return;
    }

    std::vector<std::pair<std::pair<edge_index, edge_index>, size_t>> range_stack;
    range_stack.emplace_back(index_range, length);
    while (range_stack.size()) {
        auto [cur_range, cur_length] = std::move(range_stack.back());
        range_stack.pop_back();
        if (cur_length == boss.get_k()) {
            call_range(cur_range);
            continue;
        }

        ++cur_length;

        boss.call_tightened_ranges(cur_range.first, cur_range.second, [&](auto first, auto second, auto) {
            range_stack.emplace_back(std::make_pair(first, second), cur_length);
        });
    }
}

const DBGSuccinct& get_base_dbg_succ(const DeBruijnGraph *graph) {
    if (const auto *wrapper = dynamic_cast<const DBGWrapper<>*>(graph))
        graph = &wrapper->get_graph();

    try {
        return dynamic_cast<const DBGSuccinct&>(*graph);
    } catch (const std::bad_cast &e) {
        logger->error("SuffixSeeder can be used only with succinct graph representation");
        throw e;
    }
}

template <class BOSSRange>
void append_range_nodes(const DBGSuccinct &dbg_succ,
                        size_t min_seed_length,
                        size_t max_seed_length,
                        size_t max_num_seeds_per_locus,
                        bool no_seed_complexity_filter,
                        std::string_view query,
                        const std::function<void(size_t, size_t, node_index)> &callback,
                        const CanonicalDBG *canonical,
                        const std::unique_ptr<boss::BOSS> &sample_boss,
                        const std::vector<BOSSRange> &matching_ranges) {
    const boss::BOSS &boss = dbg_succ.get_boss();
    auto encoded = boss.encode(query);
    size_t seed_length_cutoff = std::min(boss.get_k(), max_seed_length);
    std::pair<edge_index, edge_index> range;
    size_t length = 0;
    size_t num_positions = query.size() - min_seed_length + 1;
    std::vector<std::pair<edge_index, edge_index>> edge_matched_ranges;

    if (sample_boss) {
        std::vector<edge_index> edges = sample_boss->map_to_edges(encoded);
        assert(edges.size() == num_positions);
        edge_matched_ranges.reserve(num_positions);
        std::transform(edges.begin(), edges.end(), std::back_inserter(edge_matched_ranges),
                       [&](edge_index e) { return matching_ranges[e]; });
    }

    size_t i = 0;
    if (edge_matched_ranges[0].first) {
        range = edge_matched_ranges[0];
        length = min_seed_length;
        i = 1;
    }

    for ( ; i <= num_positions; ++i) {
        auto begin = encoded.begin() + i;
        if (length) {
            assert(length <= boss.get_k());
            assert(length >= min_seed_length);
            size_t j = i - (length - min_seed_length) - 1;
            assert(j < query.size());
            assert(j + length == i + min_seed_length - 1);
            assert(j + length <= query.size());

            assert(std::string_view(query.data() + j, length)
                == boss.get_node_str(range.second).substr(boss.get_k() - length));

            auto end = encoded.begin() + j + length;
            if (length < seed_length_cutoff && end < encoded.end()) {
                if (boss.tighten_range(&range.first, &range.second, *end)) {
                    ++length;
                    continue;
                }
            }

            if (range.second - range.first + 1 <= max_num_seeds_per_locus) {
                std::string_view query_window(query.data() + j, length);
                if (no_seed_complexity_filter || !is_low_complexity(query_window)) {
                    if (!canonical) {
                        call_ones(boss.get_last(), range.first, range.second + 1,
                                  [&](edge_index edge) {
                            assert(boss.get_node_str(edge).substr(boss.get_k() - length)
                                == query_window);
                            boss.call_incoming_to_target(boss.bwd(edge),
                                                         boss.get_node_last_value(edge),
                                [&](edge_index e) {
                                    if (node_index n = dbg_succ.boss_to_kmer_index(e))
                                        callback(j, length, n);
                                }
                            );
                        });
                    } else {
                        // matching is done query prefix -> node suffix
                        // e.g.,
                        // k = 6;
                        // rev: rev_end_pos = 8
                        //     j
                        //     ****--      <-- start position in forward depends on match length
                        // GCTAGCATCTGAGAGGGGA fwd
                        // TCCCCTCTCAGATGCTAGC rc
                        //          --****
                        //            i    <-- match returned from call
                        suffix_to_prefix(dbg_succ, range, length,
                                         [&](node_index n) {
                            assert(std::string_view(query.data() + j, length)
                                == dbg_succ.get_node_sequence(n).substr(0, length));
                            callback(query.size() - length - j, length,
                                     canonical->reverse_complement(n));
                        });
                    }
                }
            }
        }

        if (i < num_positions) {
            if (edge_matched_ranges.size()) {
                if (edge_matched_ranges[i].first) {
                    range = edge_matched_ranges[i];
                    length = min_seed_length;
                } else {
                    length = 0;
                }
            } else {
                auto end = begin + min_seed_length;
                auto [first, last, match_end] = boss.index_range(begin, end);
                if (match_end == end) {
                    range = std::make_pair(first, last);
                    length = min_seed_length;
                } else {
                    length = 0;
                }
            }

            assert(!length || std::string_view(query.data() + i, length)
                == boss.get_node_str(range.second).substr(boss.get_k() - length));
        }
    }
}

template <class BaseSeeder>
void SuffixSeeder<BaseSeeder>
::generate_seeds(const std::unique_ptr<boss::BOSS> &sample_boss,
                 const std::vector<BOSSRange> &matching_ranges) {
    static_assert(std::is_base_of_v<ExactSeeder, BaseSeeder>);

    if (this->query_.size() < this->config_.min_seed_length)
        return;

    if (this->config_.min_seed_length >= this->graph_.get_k()) {
        seeds_ = this->BaseSeeder::get_seeds();
        return;
    }

    const DBGSuccinct &dbg_succ = get_base_dbg_succ(&this->graph_);
    auto add_seed = [&](size_t clipping, size_t length, node_index n) {
        seeds_.emplace_back(std::string_view(this->query_.data() + clipping, length),
                            std::vector<node_index>{ n }, this->orientation_,
                            dbg_succ.get_k() - length, clipping,
                            this->query_.size() - clipping - length);
    };

    append_range_nodes(dbg_succ, this->config_.min_seed_length,
                       this->config_.max_seed_length,
                       this->config_.max_num_seeds_per_locus,
                       this->config_.no_seed_complexity_filter, this->query_, add_seed,
                       nullptr, sample_boss, matching_ranges);

    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_)) {
        // find sub-k matches in the reverse complement
        std::string query_rc(this->query_);
        reverse_complement(query_rc.begin(), query_rc.end());
        append_range_nodes(dbg_succ, this->config_.min_seed_length,
                           this->config_.max_seed_length,
                           this->config_.max_num_seeds_per_locus,
                           this->config_.no_seed_complexity_filter, query_rc,
                           add_seed, canonical, sample_boss, matching_ranges);
    }

    std::sort(seeds_.begin(), seeds_.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_clipping(), b.get_end_clipping(), a.get_nodes()[0])
            < std::make_tuple(b.get_clipping(), a.get_end_clipping(), b.get_nodes()[0]);
    });
}

auto MEMSeeder::get_seeds() const -> std::vector<Seed> {
    size_t k = graph_.get_k();

    if (k >= config_.max_seed_length)
        return ExactSeeder::get_seeds();

    if (num_matching_ < config_.min_exact_match * query_.size())
        return {};

    std::vector<uint8_t> query_node_flags(query_nodes_.size(), 0);
    for (size_t i = 0; i < query_node_flags.size(); ++i) {
        if (query_nodes_[i] != DeBruijnGraph::npos) {
            // the second bit indicates that a node has been found, while the
            // first bit indicates if the node is a maximal exact match terminus
            query_node_flags[i] = 2 |
                (i + 1 == query_nodes_.size()
                    || query_nodes_[i + 1] == DeBruijnGraph::npos
                    || get_mem_terminator()[query_nodes_[i]]);

        }
    }

    std::vector<Seed> seeds;

    // find start of MEM
    auto it = query_node_flags.begin();
    while ((it = std::find_if(it, query_node_flags.end(),
                              [](uint8_t flags) { return flags & 2; }))
            != query_node_flags.end()) {
        // find end of MEM
        auto next = std::find_if(
            it, query_node_flags.end(),
            [](uint8_t flags) { return (flags & 1) == 1 || (flags & 2) == 0; }
        );

        if (next != query_node_flags.end() && ((*next) & 2))
            ++next;

        assert(next > it);
        assert(next <= query_node_flags.end());

        size_t i = it - query_node_flags.begin();
        assert(it == query_node_flags.end()
                || query_nodes_[i] != DeBruijnGraph::npos);

        size_t mem_length = (next - it) + k - 1;
        assert(i + mem_length <= query_.size());

        if (mem_length >= config_.min_seed_length) {
            const char *begin_it = query_.data() + i;
            auto node_begin_it = query_nodes_.begin() + i;
            auto node_end_it = node_begin_it + (next - it);
            assert(std::find(node_begin_it, node_end_it, DeBruijnGraph::npos)
                    == node_end_it);

            seeds.emplace_back(std::string_view(begin_it, mem_length),
                               std::vector<node_index>{ node_begin_it, node_end_it },
                               orientation_, 0, i, query_.size() - i - mem_length);
        }

        it = next;
    }

    return seeds;
}

template class SuffixSeeder<ExactSeeder>;
template class SuffixSeeder<UniMEMSeeder>;

} // namespace align
} // namespace graph
} // namespace mtg
