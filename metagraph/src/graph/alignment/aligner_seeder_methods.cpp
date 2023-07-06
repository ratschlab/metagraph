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
typedef boss::BOSS::edge_index edge_index;
typedef boss::BOSS::TAlphabet TAlphabet;


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
    assert(std::get<0>(index_range));
    assert(std::get<1>(index_range));
    assert(std::get<2>(index_range));
    assert(std::get<2>(index_range) < dbg_succ.get_k());

#ifndef NDEBUG
    size_t offset = boss.get_k() - std::get<2>(index_range);
    std::string check_str = boss.get_node_str(std::get<0>(index_range)).substr(offset);
    assert(std::get<0>(index_range) == 1
        || boss.get_node_str(std::get<0>(index_range) - 1).substr(offset) != check_str);

    assert(boss.get_node_str(std::get<1>(index_range)).substr(offset) == check_str);
    assert(std::get<1>(index_range) == boss.get_W().size() - 1
        || boss.get_node_str(std::get<1>(index_range) + 1).substr(offset) != check_str);
#endif

    auto call_nodes_in_range = [&](const BOSSEdgeRange &final_range) {
        const auto &[first, last, seed_length] = final_range;
        assert(seed_length == boss.get_k());
        for (boss::BOSS::edge_index i = first; i <= last; ++i) {
            assert(boss.get_node_str(i).substr(0, std::get<2>(index_range)) == check_str);
            if (auto node = dbg_succ.boss_to_kmer_index(i)) {
                assert(dbg_succ.get_node_sequence(node).substr(0, std::get<2>(index_range))
                    == check_str);
                callback(node);
            }
        }
    };

    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(index_range);
        return;
    }

    std::vector<BOSSEdgeRange> range_stack;
    range_stack.emplace_back(index_range);

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
    assert(this->config_.min_seed_length);
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

    seeds_.clear();
    sdsl::bit_vector matching(this->query_.size(), false);

    std::vector<std::tuple<edge_index, edge_index, TAlphabet, edge_index, edge_index>> ranges(
        this->query_.size() - this->config_.min_seed_length + 1
    );

    const auto &boss = dbg_succ.get_boss();
    for (size_t i = 0; i + this->config_.min_seed_length <= this->query_.size(); ++i) {
        auto &[first, last, s, first_rc, last_rc] = ranges[i];

        std::string_view window(this->query_.data() + i, this->config_.min_seed_length);
        s = boss.encode(window.back()) % boss.alph_size;
        if (!s)
            continue;

        bool low_complexity = this->config_.seed_complexity_filter
            ? is_low_complexity(window)
            : false;
        bool found = false;

        std::string_view window_prefix(window.data(), window.size() - 1);
        auto encoded = boss.encode(window_prefix);
        auto end = encoded.begin();

        std::tie(first, last, end) = boss.index_range(encoded.begin(), encoded.end());

        if (end == encoded.end()) {
            found = true;
            low_complexity |= (last - first > this->config_.max_num_seeds_per_locus);
            first = !low_complexity ? boss.pred_last(first - 1) + 1 : 0;
            if (first && i + this->config_.min_seed_length < this->query_.size()) {
                auto first_test = first;
                auto last_test = last;
                if (boss.tighten_range(&first_test, &last_test, s)) {
                    auto next_s = boss.encode(this->query_[i + this->config_.min_seed_length]) % boss.alph_size;
                    if (next_s && first_test == last_test && boss.tighten_range(&first_test, &last_test, next_s)) {
                        first = 0;
                    }

                } else {
                    first = 0;
                }
            }
        } else {
            first = 0;
        }

        if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY && (!low_complexity || !found)) {
            assert(dynamic_cast<const CanonicalDBG*>(&this->graph_));
            std::string window_rc(window);
            ::reverse_complement(window_rc.begin(), window_rc.end());
            auto encoded = boss.encode(window_rc);
            auto end = encoded.begin();
            std::tie(first_rc, last_rc, end) = boss.index_range(encoded.begin(), encoded.end());
            if (end == encoded.end()) {
                found = true;
                low_complexity |= (last_rc - first_rc > this->config_.max_num_seeds_per_locus);
                first_rc = !low_complexity ? boss.pred_last(first_rc - 1) + 1 : 0;
                if (first_rc && i && first_rc == last_rc) {
                    auto prev_s = complement(boss.encode(this->query_[i - 1]) % boss.alph_size);
                    if (prev_s && boss.get_minus_k_value(first_rc, this->config_.min_seed_length).first == prev_s) {
                        first_rc = 0;
                    }
                }
            } else {
                first_rc = 0;
            }
        }

        if (found) {
            for (size_t j = i; j < i + this->config_.min_seed_length; ++j) {
                matching[j] = true;
            }
        }
    }

    this->num_matching_ = sdsl::util::cnt_one_bits(matching);
    if (this->num_matching_ < this->query_.size() * this->config_.min_exact_match) {
        this->num_matching_ = 0;
        return;
    }

    std::vector<tsl::hopscotch_set<node_index>> found_nodes(ranges.size());
    for (size_t i = 0; i < ranges.size(); ++i) {
#ifndef NDEBUG
        std::string_view window(this->query_.data() + i, this->config_.min_seed_length);
#endif
        auto [first, last, s, first_rc, last_rc] = ranges[i];

        std::vector<node_index> nodes;
        if (first) {
            for (auto e = boss.succ_W(first, s); e <= last; e = boss.succ_W(e + 1, s)) {
                if (auto node = dbg_succ.boss_to_kmer_index(e)) {
                    assert(dbg_succ.get_node_sequence(node).substr(dbg_succ.get_k() - window.size())
                            == window);
                    nodes.emplace_back(node);
                }

                if (e + 1 == boss.get_W().size())
                    break;
            }

            s += boss.alph_size;
            for (auto e = boss.succ_W(first, s); e <= last; e = boss.succ_W(e + 1, s)) {
                if (auto node = dbg_succ.boss_to_kmer_index(e)) {
                    assert(dbg_succ.get_node_sequence(node).substr(dbg_succ.get_k() - window.size())
                            == window);
                    nodes.emplace_back(node);
                }

                if (e + 1 == boss.get_W().size())
                    break;
            }
        }

        if (first_rc) {
            const auto *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_);
            assert(canonical);
            suffix_to_prefix(dbg_succ,
                std::make_tuple(first_rc, last_rc, this->config_.min_seed_length),
                [&](node_index node) {
                    node = canonical->reverse_complement(node);
                    assert(canonical->get_node_sequence(node).substr(dbg_succ.get_k() - window.size())
                            == window);
                    nodes.emplace_back(node);
                }
            );
        }

        for (node_index node : nodes) {
            assert(node);
            if (!found_nodes[i].emplace(node).second)
                continue;

            std::vector<node_index> path;
            path.emplace_back(node);
            size_t end_i = i + this->config_.min_seed_length;
            if (this->config_.max_seed_length > this->config_.min_seed_length
                    && end_i < this->query_.size()) {
                std::string_view rest(this->query_.data() + end_i,
                                      this->query_.size() - end_i);
                this->graph_.traverse(node, rest.begin(), rest.end(),
                    [&](node_index next) {
                        found_nodes[i + path.size()].emplace(next);
                        path.emplace_back(next);
                    },
                    [&]() {
                        return this->config_.min_seed_length + path.size() - 1
                                    >= this->config_.max_seed_length
                            || this->graph_.has_multiple_outgoing(path.back())
                            || !this->graph_.has_single_incoming(path.back());
                    }
                );
            }

            std::string_view seed_window(this->query_.data() + i,
                                         this->config_.min_seed_length + path.size() - 1);

            size_t offset = this->graph_.get_k() - this->config_.min_seed_length;
            if (path.size() > 1 && offset) {
                if (path.size() - 1 <= offset) {
                    offset -= path.size() - 1;
                    path.assign(1, path.back());
                } else {
                    path.erase(path.begin(), path.begin() + offset);
                    offset = 0;
                }
            }

            seeds_.emplace_back(
                seed_window,
                std::move(path),
                this->orientation_,
                offset,
                i,
                this->query_.size() - i - seed_window.size()
            );
        }
    }

    if (seeds_.empty())
        this->num_matching_ = 0;
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

    assert(std::all_of(begin, end, [](const auto &a) { return a.get_nodes().size(); }));

    using seed_t = std::remove_reference_t<decltype(*begin)>;

    if constexpr(std::is_same_v<seed_t, Alignment>) {
        // first, move all inexact matches to the front and ignore them
        begin = std::partition(begin, end, [](const auto &a) {
            const auto &cigar = a.get_cigar().data();
            auto c_begin = cigar.begin();
            auto c_end = cigar.end();
            assert(c_begin != c_end);

            if (c_begin->first == Cigar::CLIPPED)
                ++c_begin;

            assert(c_begin != c_end);

            if ((c_end - 1)->first == Cigar::CLIPPED)
                --c_end;

            return c_end != c_begin + 1 || c_begin->first != Cigar::MATCH;
        });

        if (begin == end)
            return end;
    }

    ssize_t graph_k = graph.get_k();
    std::sort(begin, end, [](const auto &a, const auto &b) {
        return std::pair(a.get_query_view().end(), a.get_query_view().begin())
            > std::pair(b.get_query_view().end(), b.get_query_view().begin());
    });

    static_assert((std::is_same_v<seed_t, Seed> || std::is_same_v<seed_t, Alignment>)
        && "Only implemented for Seed and Alignment"
    );

    auto clear_seed = [](auto &seed) { seed = seed_t(); };

    // first, discard redundant seeds
    for (auto i = begin; i + 1 != end; ++i) {
        auto &a_i = *(i + 1);
        auto &a_j = *i;

        const auto &nodes_i = a_i.get_nodes();
        const auto &nodes_j = a_j.get_nodes();
        assert(nodes_i.size());
        assert(nodes_j.size());

        if (a_i.get_end_clipping() != a_j.get_end_clipping())
            continue;

        if (a_i.get_clipping() == a_j.get_clipping() && a_i.get_offset() == a_j.get_offset()
                && nodes_i == nodes_j) {
            // these are the same alignment, merge their annotations
            if (a_i.label_columns.empty() || a_j.label_columns.empty()) {
                if (a_i.label_columns.empty())
                    std::swap(a_i, a_j);

                clear_seed(a_j);
                assert(a_i.get_nodes().size());
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
            clear_seed(a_j);
            assert(a_i.get_nodes().size());
            continue;
        }

        if (a_i.label_columns != a_j.label_columns)
            continue;

        std::string_view query_i = a_i.get_query_view();
        std::string_view query_j = a_j.get_query_view();

        assert(nodes_i.size());
        assert(nodes_j.size());
        if (nodes_j.back() == nodes_i.back()) {
            if (query_j.size() > query_i.size())
                std::swap(a_i, a_j);

            clear_seed(a_j);
            assert(a_i.get_nodes().size());
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
        auto &a_i = *(i + 1);
        if (a_i.get_query_view().size() >= max_seed_size)
            continue;

        auto &a_j = *i;

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

                clear_seed(a_j);
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

        char next_c = *(query_i.data() + query_i.size());
        assert(overlap < graph_k - 1
                || graph.traverse(nodes_i.back(), next_c) == nodes_j[a_j_node_idx]);

        if (overlap < graph_k - 1 && graph.traverse(nodes_i.back(), next_c)
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
        std::vector<Alignment::node_index> added_nodes(nodes_j.begin() + a_j_node_idx, nodes_j.end());
        if constexpr(std::is_same_v<seed_t, Seed>) {
            a_i.expand(std::move(added_nodes));
            assert(Alignment(a_i, config).is_valid(graph, &config));
            clear_seed(a_j);
        }

        if constexpr(std::is_same_v<seed_t, Alignment>) {
            std::string_view added_query(query_j.data() + query_j.size() - added_nodes.size(), added_nodes.size());
            Alignment inserted_seed(
                Seed(added_query,
                     std::move(added_nodes),
                     a_j.get_orientation(),
                     graph.get_k() - 1,
                     a_j.get_clipping() + query_j.size() - added_query.size(),
                     a_j.get_end_clipping()),
                config
            );
            inserted_seed.label_columns = a_j.label_columns;
            inserted_seed.label_coordinates = a_j.label_coordinates;
            size_t coord_diff = inserted_seed.get_clipping() - a_j.get_clipping();
            for (auto &tuple : inserted_seed.label_coordinates) {
                for (auto &c : tuple) {
                    c += coord_diff;
                }
            }
            assert(inserted_seed.is_valid(graph, &config));
            a_i.splice(std::move(inserted_seed));
            assert(a_i.is_valid(graph, &config));
            clear_seed(a_j);
        }
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

template std::vector<Alignment>::iterator merge_into_unitig_mums(const DeBruijnGraph &,
                                                                 const DBGAlignerConfig &,
                                                                 std::vector<Alignment>::iterator,
                                                                 std::vector<Alignment>::iterator,
                                                                 ssize_t,
                                                                 size_t);

} // namespace align
} // namespace graph
} // namespace mtg
