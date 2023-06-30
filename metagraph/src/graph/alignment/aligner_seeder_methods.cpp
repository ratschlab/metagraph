#include "aligner_seeder_methods.hpp"

#include <sdust.h>
#include <tsl/hopscotch_set.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef Alignment::score_t score_t;


#if ! _PROTEIN_GRAPH
bool is_low_complexity(std::string_view s, int T, int W) {
    int n;
    std::unique_ptr<uint64_t, decltype(std::free)*> r {
        sdust(0, (const uint8_t*)s.data(), s.size(), T, W, &n),
        std::free
    };
    return n > 0;
}
#else
bool is_low_complexity(std::string_view, int, int) {
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

    std::vector<Seed> seeds;

    if (config_.max_seed_length < k)
        return seeds;

    size_t end_clipping = query_.size() - k;
    for (size_t i = 0; i < query_nodes_.size(); ++i, --end_clipping) {
        if (query_nodes_[i] != DeBruijnGraph::npos) {
            assert(i + k <= query_.size());
            std::string_view query_window = query_.substr(i, k);
            seeds.emplace_back(query_window,
                               std::vector<node_index>{ query_nodes_[i] },
                               orientation_, 0, i, end_clipping);
        }
    }

    return seeds;
}

template <class BOSSEdgeRange>
void suffix_to_prefix(const DBGSuccinct &dbg_succ,
                      const BOSSEdgeRange &index_range,
                      const std::function<void(DBGSuccinct::node_index)> &callback) {
    const auto &boss = dbg_succ.get_boss();
    assert(std::get<2>(index_range));
    assert(std::get<2>(index_range) < dbg_succ.get_k());

    auto call_nodes_in_range = [&](const BOSSEdgeRange &final_range) {
        const auto &[first, last, seed_length] = final_range;
        assert(seed_length == boss.get_k());
        for (boss::BOSS::edge_index i = first; i <= last; ++i) {
            DBGSuccinct::node_index node = dbg_succ.boss_to_kmer_index(i);
            if (node)
                callback(node);
        }
    };

    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(index_range);
        return;
    }

    std::vector<BOSSEdgeRange> range_stack { index_range };

    while (range_stack.size()) {
        BOSSEdgeRange cur_range = std::move(range_stack.back());
        range_stack.pop_back();
        assert(std::get<2>(cur_range) < boss.get_k());
        ++std::get<2>(cur_range);

        for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
            auto next_range = cur_range;
            auto &[first, last, seed_length] = next_range;

            if (boss.tighten_range(&first, &last, s)) {
                if (seed_length == boss.get_k()) {
                    call_nodes_in_range(next_range);
                } else {
                    range_stack.emplace_back(std::move(next_range));
                }
            }
        }
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

template <class BaseSeeder>
void SuffixSeeder<BaseSeeder>::generate_seeds() {
    typedef typename BaseSeeder::node_index node_index;

    // this method assumes that seeds from the BaseSeeder are exact match only
    static_assert(std::is_base_of_v<ExactSeeder, BaseSeeder>);

    if (this->query_.size() < this->config_.min_seed_length)
        return;

    if (this->config_.min_seed_length >= this->graph_.get_k()) {
        seeds_ = this->BaseSeeder::get_seeds();
        return;
    }

    const DBGSuccinct &dbg_succ = get_base_dbg_succ(&this->graph_);
    if (dbg_succ.get_mask())
        logger->warn("Graph has a dummy k-mer mask. Seeds containing dummy k-mers will be missed.");

    bool found_first = false;
    std::vector<std::vector<Seed>> found_seeds(this->query_.size() - this->config_.min_seed_length + 1);
    size_t total_seed_count = 0;
    if (this->query_.size() >= this->graph_.get_k()) {
        if (this->config_.max_seed_length >= this->graph_.get_k()) {
            assert(this->query_nodes_.size()
                == this->query_.size() - this->graph_.get_k() + 1);
            for (auto &seed : this->BaseSeeder::get_seeds()) {
                found_first |= !seed.get_clipping();
                auto &bucket = found_seeds[seed.get_end_clipping()];
                bucket.emplace_back(std::move(seed));
                ++total_seed_count;
            }
        } else {
            std::string_view window(this->query_.data(), this->graph_.get_k());
            auto first_path = map_to_nodes_sequentially(this->graph_, window);
            assert(first_path.size() == 1);
            if (first_path[0]) {
                found_first = true;
                size_t end_clipping = this->query_.size() - window.size();
                found_seeds[end_clipping].emplace_back(
                    window, std::move(first_path), this->orientation_,
                    this->graph_.get_k() - window.size(),
                    0, end_clipping
                );
                ++total_seed_count;
            }
        }
    }

    auto add_seeds = [&](size_t i, size_t max_seed_length) {
        std::string_view max_window(this->query_.data() + i, max_seed_length);
        dbg_succ.call_nodes_with_suffix_matching_longest_prefix(max_window,
            [&](node_index alt_node, size_t seed_len) {
                std::string_view window(this->query_.data() + i, seed_len);
                size_t end_clipping = this->query_.size() - i - window.size();
                auto &bucket = found_seeds[end_clipping];
                if (bucket.size()) {
                    if (seed_len < bucket[0].get_query_view().size())
                        return;

                    if (seed_len > bucket[0].get_query_view().size()) {
                        total_seed_count -= bucket.size();
                        bucket.clear();
                    }
                }

                bucket.emplace_back(window, std::vector<node_index>{ alt_node },
                                    this->orientation_,
                                    this->graph_.get_k() - window.size(),
                                    i, end_clipping);
                found_first |= !i;
                ++total_seed_count;
            },
            this->config_.min_seed_length
        );
    };

    size_t max_seed_length = std::min(this->graph_.get_k() - 1,
                                      this->config_.max_seed_length);
    size_t i = found_first ? this->graph_.get_k() - max_seed_length : 0;
    for ( ; i + max_seed_length <= this->query_.size(); ++i) {
        add_seeds(i, max_seed_length);
    }

    assert(i == this->query_.size() - max_seed_length + 1);
    if (this->config_.min_seed_length < max_seed_length) {
        size_t cur_length = max_seed_length;
        for ( ; i + this->config_.min_seed_length <= this->query_.size(); ++i) {
            add_seeds(i, --cur_length);
        }
    }

    if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
        const auto *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_);
        assert(canonical);
        std::string query_rc(this->query_);
        ::reverse_complement(query_rc.begin(), query_rc.end());
        auto add_seeds = [&](size_t i, size_t max_seed_length) {
            std::string_view max_window_rc(query_rc.data() + i, max_seed_length);
            tsl::hopscotch_set<node_index> found_nodes;

            const auto &boss = dbg_succ.get_boss();
            auto encoded = boss.encode(max_window_rc);
            auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
            size_t seed_len = end - encoded.begin();
            if (seed_len < this->config_.min_seed_length)
                return;

            size_t clipping = this->query_.size() - i - seed_len;
            if (found_first && !clipping)
                return;

            std::string_view window(this->query_.data() + clipping, seed_len);
            auto &bucket = found_seeds[i];
            if (bucket.size()) {
                if (seed_len < bucket[0].get_query_view().size())
                    return;

                if (seed_len > bucket[0].get_query_view().size()) {
                    total_seed_count -= bucket.size();
                    bucket.clear();
                }
            }

            suffix_to_prefix(dbg_succ,
                std::make_tuple(boss.pred_last(first - 1) + 1, last, seed_len),
                [&](node_index alt_node) {
                    found_nodes.emplace(canonical->reverse_complement(alt_node));
                }
            );

            for (node_index alt_node : found_nodes) {
                assert(this->graph_.get_node_sequence(alt_node).substr(
                    this->graph_.get_k() - window.size()) == window);
                bucket.emplace_back(window, std::vector<node_index>{ alt_node },
                                    this->orientation_,
                                    this->graph_.get_k() - window.size(),
                                    clipping, i);
                found_first |= !clipping;
                ++total_seed_count;
            }
        };

        size_t i = 0;
        for ( ; i + max_seed_length <= query_rc.size(); ++i) {
            add_seeds(i, max_seed_length);
        }

        assert(i == this->query_.size() - max_seed_length + 1);
        if (this->config_.min_seed_length < max_seed_length) {
            size_t cur_length = max_seed_length;
            for ( ; i + this->config_.min_seed_length <= this->query_.size(); ++i) {
                add_seeds(i, --cur_length);
            }
        }
    }

    auto first_front_match = std::find_if(found_seeds.begin(), found_seeds.end(),
        [](const auto &bucket) {
            return bucket.size() && !bucket[0].get_clipping();
        }
    );

    if (first_front_match != found_seeds.end()) {
        std::for_each(first_front_match + 1, found_seeds.end(), [](auto &bucket) {
            bucket.clear();
        });
    }

    seeds_.clear();
    seeds_.reserve(total_seed_count);
    for (auto &bucket : found_seeds) {
        seeds_.insert(seeds_.end(),
                      std::make_move_iterator(bucket.begin()),
                      std::make_move_iterator(bucket.end()));
    }

    this->num_matching_ = get_num_char_matches_in_seeds(seeds_.begin(), seeds_.end());
}

auto MEMSeeder::get_seeds() const -> std::vector<Seed> {
    size_t k = graph_.get_k();

    if (k >= config_.max_seed_length)
        return ExactSeeder::get_seeds();

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

template <typename It>
It merge_into_unitig_mums(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          It begin,
                          It end,
                          ssize_t min_seed_size,
                          size_t max_seed_size) {
    if (begin == end)
        return end;

    ssize_t graph_k = graph.get_k();
    std::sort(begin, end, [](const auto &a, const auto &b) {
        return std::pair(a.get_query_view().end(), a.get_query_view().begin())
            > std::pair(b.get_query_view().end(), b.get_query_view().begin());
    });

    // first, discard redundant seeds
    for (auto i = begin; i + 1 != end; ++i) {
        Seed &a_i = *(i + 1);
        Seed &a_j = *i;

        if (a_i.get_end_clipping() != a_j.get_end_clipping())
            continue;

        const auto &nodes_i = a_i.get_nodes();
        const auto &nodes_j = a_j.get_nodes();
        if (a_i.get_clipping() == a_j.get_clipping() && a_i.get_offset() == a_j.get_offset()
                && nodes_i == nodes_j) {
            // these are the same alignment, merge their annotations
            if (a_i.label_columns.empty() || a_j.label_columns.empty()) {
                if (a_i.label_columns.empty())
                    std::swap(a_i, a_j);

                a_j = Seed();
                continue;
            }

            assert(a_i.label_coordinates.empty() == a_j.label_coordinates.empty());

            Alignment::Columns merged_columns;
            if (a_i.label_coordinates.empty()) {
                std::set_union(a_i.label_columns.begin(), a_i.label_columns.end(),
                               a_j.label_columns.begin(), a_j.label_columns.end(),
                               std::back_inserter(merged_columns));
            } else {
                Alignment::CoordinateSet merged_coords;
                auto add_diff = [&](auto label, const auto &c) {
                    merged_columns.emplace_back(label);
                    merged_coords.emplace_back(c);
                };
                utils::match_indexed_values(a_i.label_columns.begin(), a_i.label_columns.end(),
                                            a_i.label_coordinates.begin(),
                                            a_j.label_columns.begin(), a_j.label_columns.end(),
                                            a_j.label_coordinates.begin(),
                    [&](auto label, const auto &c1, const auto &c2) {
                        merged_columns.emplace_back(label);
                        auto &c = merged_coords.emplace_back();
                        std::set_union(c1.begin(), c1.end(), c2.begin(), c2.end(),
                                       std::back_inserter(c));
                    },
                    add_diff,
                    add_diff
                );
                std::swap(a_i.label_coordinates, merged_coords);
            }

            std::swap(a_i.label_columns, merged_columns);
            a_j = Seed();
            continue;
        }

        if (a_i.label_columns != a_j.label_columns)
            continue;

        std::string_view query_i = a_i.get_query_view();
        std::string_view query_j = a_j.get_query_view();

        if (nodes_j.back() == nodes_i.back()) {
            if (query_j.size() > query_i.size())
                std::swap(a_i, a_j);

            a_j = Seed();
        }
    }

    end = std::remove_if(begin, end, [](const auto &a) { return a.empty(); });

    size_t query_size = begin->get_clipping() + begin->get_end_clipping()
                        + begin->get_query_view().size();
    sdsl::int_vector<2> end_counter(query_size, 0);
    std::for_each(begin, end, [&](const auto &a) {
        size_t i = a.get_end_clipping();
        if (end_counter[i] < 2)
            ++end_counter[i];
    });
    for (auto i = begin; i + 1 != end; ++i) {
        // try to merge a_i to a_j
        Seed &a_i = *(i + 1);
        if (a_i.get_query_view().size() >= max_seed_size)
            continue;

        Seed &a_j = *i;

        if (a_i.label_columns != a_j.label_columns)
            continue;

        bool coordinates_consistent = true;
        assert(a_i.label_coordinates.size() == a_j.label_coordinates.size());
        auto jt = a_j.label_coordinates.begin();
        for (auto &tuple : a_i.label_coordinates) {
            assert(jt != a_j.label_coordinates.end());
            if (tuple.size() != jt->size()) {
                coordinates_consistent = false;
                break;
            }
            ++jt;
        }

        if (!coordinates_consistent)
            continue;

        assert(jt == a_j.label_coordinates.end());

        const auto &nodes_i = a_i.get_nodes();
        const auto &nodes_j = a_j.get_nodes();
        std::string_view query_i = a_i.get_query_view();
        std::string_view query_j = a_j.get_query_view();

        // alignments are disjoint
        if (query_i.end() <= query_j.begin())
            continue;

        ssize_t num_added = query_j.end() - std::max(query_j.begin(), query_i.end());
        ssize_t overlap = query_i.end() - query_j.begin();
        if (num_added < 0 || overlap < min_seed_size - 1)
            continue;

        if (num_added == 0) {
            if (nodes_i.back() == nodes_j.back()) {
                if (query_j.size() > query_i.size())
                    std::swap(a_i, a_j);

                a_j = Seed();
            }
            continue;
        }

        // we want query_j.begin() + graph_k - a_j.get_offset() + x == query_i.end() + 1
        // ->      graph_k - a_j.get_offset() + x == overlap + 1
        // -> x == overlap + 1 + a_j.get_offset() - graph_k
        ssize_t a_j_node_idx = overlap + 1 + static_cast<ssize_t>(a_j.get_offset()) - graph_k;
        assert(a_j_node_idx < static_cast<ssize_t>(nodes_j.size()));

        if (a_j_node_idx < 0)
            continue;

        int64_t coord_dist = nodes_j.size() - a_j_node_idx;
        int64_t dist = query_j.end() - query_i.end();
        if (coord_dist != dist)
            continue;

        bool unique = true;
        for (size_t i = a_j.get_end_clipping(); i < a_i.get_end_clipping(); ++i) {
            if (end_counter[i] == 2) {
                unique = false;
                break;
            }
        }

        if (!unique)
            continue;

        if (graph.has_multiple_outgoing(nodes_i.back())
                || !graph.has_single_incoming(nodes_i.back()))
            continue;

        assert(overlap < graph_k - 1
                || graph.traverse(nodes_i.back(), *query_i.end()) == nodes_j[a_j_node_idx]);

        if (overlap < graph_k - 1 && graph.traverse(nodes_i.back(), *query_i.end())
                                        != nodes_j[a_j_node_idx])
            continue;

        jt = a_j.label_coordinates.begin();
        if (!coordinates_consistent)
            continue;

        for (auto &tuple : a_i.label_coordinates) {
            assert(jt != a_j.label_coordinates.end());
            assert(tuple.size() == jt->size());

            auto jt_c = jt->begin();
            for (ssize_t c : tuple) {
                assert(jt_c != jt->end());

                if (c + static_cast<ssize_t>(nodes_i.size()) != *jt_c + a_j_node_idx) {
                    coordinates_consistent = false;
                    break;
                }

                ++jt_c;
            }

            if (!coordinates_consistent)
                break;

            assert(jt_c == jt->end());
            ++jt;
        }

        assert(jt == a_j.label_coordinates.end());

        // we have a MUM
        a_i.expand(std::vector<Alignment::node_index>(nodes_j.begin() + a_j_node_idx,
                                                      nodes_j.end()));
        assert(Alignment(a_i, config).is_valid(graph, &config));
        a_j = Seed();
    }

    return std::remove_if(begin, end, [](const auto &a) { return a.empty(); });
}

template Seed* merge_into_unitig_mums(const DeBruijnGraph &,
                                      const DBGAlignerConfig &,
                                      Seed*,
                                      Seed*,
                                      ssize_t,
                                      size_t);
template std::vector<Seed>::iterator merge_into_unitig_mums(const DeBruijnGraph &,
                                                            const DBGAlignerConfig &,
                                                            std::vector<Seed>::iterator,
                                                            std::vector<Seed>::iterator,
                                                            ssize_t,
                                                            size_t);

} // namespace align
} // namespace graph
} // namespace mtg
