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
                      std::string_view rest,
                      const BOSSEdgeRange &index_range,
                      const std::function<void(DBGSuccinct::node_index, size_t)> &callback) {
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

    auto call_nodes_in_range = [&](size_t num_exact_match, const BOSSEdgeRange &final_range) {
        const auto &[first, last, seed_length] = final_range;
        assert(seed_length == boss.get_k());
        assert(num_exact_match <= seed_length);
        for (boss::BOSS::edge_index i = first; i <= last; ++i) {
            assert(boss.get_node_str(i).substr(0, std::get<2>(index_range)) == check_str);
            if (auto node = dbg_succ.boss_to_kmer_index(i)) {
                assert(dbg_succ.get_node_sequence(node).substr(0, std::get<2>(index_range))
                    == check_str);
                callback(node, num_exact_match);
            }
        }
    };

    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(boss.get_k(), index_range);
        return;
    }

    auto encoded = boss.encode(rest);
    std::vector<std::tuple<size_t, bool, BOSSEdgeRange>> range_stack;
    range_stack.emplace_back(0, true, index_range);

    while (range_stack.size()) {
        auto [num_extra_match, is_exact_match, cur_range] = std::move(range_stack.back());
        range_stack.pop_back();
        assert(std::get<2>(cur_range) < boss.get_k());
        ++std::get<2>(cur_range);

        for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
            auto next_range = cur_range;
            auto &[first, last, seed_length] = next_range;

            if (boss.tighten_range(&first, &last, s)) {
                if (seed_length == boss.get_k()) {
                    call_nodes_in_range(std::get<2>(index_range) + num_extra_match, next_range);
                } else {
                    bool next_exact_match = is_exact_match && (s == encoded[num_extra_match]);
                    range_stack.emplace_back(
                        num_extra_match + next_exact_match,
                        next_exact_match,
                        std::move(next_range)
                    );
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
    const auto &boss = dbg_succ.get_boss();

    sdsl::bit_vector matched(this->query_.size(), false);

    auto generate_from_query = [&](std::string_view query, auto find_nodes, bool is_rc) {
        std::vector<std::vector<std::pair<edge_index, edge_index>>> ranges(
            query.size() - this->config_.min_seed_length + 1
        );

        auto encoded = boss.encode(query);
        for (size_t i = 0; i + this->config_.min_seed_length <= query.size(); ++i) {
            auto begin = encoded.begin() + i;
            auto end = begin + this->config_.min_seed_length - 1;

            if (!((*end) % boss.alph_size))
                continue;

            auto [first, last, it] = boss.index_range(begin, end);

            if (it != end)
                continue;

            first = boss.pred_last(first - 1) + 1;
            auto last_it = std::min({ begin + dbg_succ.get_k(),
                                      encoded.end(),
                                      begin + this->config_.max_seed_length });
            for (size_t j = i; it != last_it; ++j, ++it) {
                assert(it <= begin + boss.get_k());
                edge_index first_next = first;
                edge_index last_next = last;
                if (boss.tighten_range(&first_next, &last_next, *it)) {
                    if (ranges[j].size() <= j - i)
                        ranges[j].resize(j - i + 1);

                    ranges[j][j - i] = std::make_pair(first, last);

                    // TODO: how do we deal with this for the rc strand?
                    if (is_rc)
                        break;

                    first = first_next;
                    last = last_next;
                } else {
                    break;
                }
            }

            if (ranges[i].size()) {
                if (is_rc) {
                    std::fill(matched.end() - i - this->config_.min_seed_length,
                              matched.end() - i,
                              true);
                } else {
                    std::fill(matched.begin() + i,
                              matched.begin() + i + this->config_.min_seed_length,
                              true);
                }
            }
        }

        for (size_t i = 0; i < ranges.size(); ++i) {
            if (ranges[i].empty())
                continue;

            assert(!is_rc || ranges[i].size() == 1);

            size_t added_length = 0;

            auto s = encoded[i + this->config_.min_seed_length - 1];
            for (auto begin = ranges[i].begin(); begin + 1 != ranges[i].end(); ++begin, ++added_length) {
                auto [first, last] = *begin;
                assert(first);
                assert(last);

                auto [first_next, last_next] = *(begin + 1);
                assert(first <= first_next);
                assert(last >= last_next);

                std::string_view seed_window(query.data() + i - added_length,
                                             this->config_.min_seed_length + added_length);

                if (this->config_.seed_complexity_filter && is_low_complexity(seed_window))
                    continue;

                if (first != first_next) {
                    find_nodes(query, i, seed_window, first, first_next - 1, s);
                    find_nodes(query, i, seed_window, first, first_next - 1, s + boss.alph_size);
                }

                if (last_next != last) {
                    find_nodes(query, i, seed_window, last_next + 1, last, s);
                    find_nodes(query, i, seed_window, last_next + 1, last, s + boss.alph_size);
                }
            }

            std::string_view seed_window(query.data() + i - added_length,
                                         this->config_.min_seed_length + added_length);

            if (this->config_.seed_complexity_filter && is_low_complexity(seed_window))
                return;

            auto [first, last] = ranges[i].back();
            assert(first);
            assert(last);
            find_nodes(query, i, seed_window, first, last, s);
            find_nodes(query, i, seed_window, first, last, s + boss.alph_size);
        }
    };

    auto add_seed = [&](std::string_view query, size_t i, std::string_view seed_window, node_index node) {
        assert(node);
        size_t added_length = seed_window.size() - this->config_.min_seed_length;
        assert(i >= added_length);
        std::vector<node_index> path;
        path.emplace_back(node);
        assert(this->config_.min_seed_length + added_length <= this->graph_.get_k());
        size_t offset = this->graph_.get_k() - this->config_.min_seed_length - added_length;
        assert(i - added_length < query.size());
        assert(query.size() - (i - added_length) - seed_window.size() < query.size());
        seeds_.emplace_back(seed_window,
                            std::move(path),
                            this->orientation_,
                            offset,
                            i - added_length,
                            query.size() - (i - added_length) - seed_window.size());
        assert(Alignment(seeds_.back(), this->config_).is_valid(this->graph_, &this->config_));
    };

    auto find_nodes_fwd = [&](std::string_view query, size_t i, std::string_view seed_window, auto first, auto last, auto s) {
        for (auto e = boss.succ_W(first, s); e <= last; e = boss.succ_W(e + 1, s)) {
            if (auto node = dbg_succ.boss_to_kmer_index(e))
                add_seed(query, i, seed_window, node);

            if (e + 1 == boss.get_W().size())
                break;
        }
    };

    generate_from_query(this->query_, find_nodes_fwd, false);

    if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
        const auto &canonical = static_cast<const CanonicalDBG&>(this->graph_);
        std::string query_rc(this->query_);
        ::reverse_complement(query_rc.begin(), query_rc.end());
        std::vector<tsl::hopscotch_map<node_index, size_t>> nodes(
            this->query_.size() - this->config_.min_seed_length + 1
        );
        auto find_nodes_bwd = [&](std::string_view, size_t i, std::string_view rc_seed_window, auto first, auto last, auto s) {
            assert(rc_seed_window.size() == this->config_.min_seed_length);
            if (s >= boss.alph_size)
                return;

            bool check = boss.tighten_range(&first, &last, s);
            std::ignore = check;
            assert(check);
            assert(boss.get_node_str(first).substr(boss.get_k() - rc_seed_window.size())
                == rc_seed_window);

            std::string_view rest(rc_seed_window.data() + rc_seed_window.size(),
                                  boss.get_k() - rc_seed_window.size());
            i = this->query_.size() - (i + rc_seed_window.size());

            suffix_to_prefix(dbg_succ,
                rest,
                std::make_tuple(first, last, rc_seed_window.size()),
                [&](node_index node, size_t num_matches) {
                    assert(num_matches >= this->config_.min_seed_length);
                    assert(num_matches <= boss.get_k());
                    node = canonical.reverse_complement(node);
                    size_t added_length = num_matches - this->config_.min_seed_length;
                    std::string_view seed_window(this->query_.data() + i - added_length,
                                                 num_matches);
                    assert(canonical.get_node_sequence(node).substr(dbg_succ.get_k() - num_matches)
                        == seed_window);
                    size_t end_clipping = this->query_.size() - (i - added_length) - seed_window.size();
                    auto it = nodes[end_clipping].try_emplace(node, num_matches).first;
                    it.value() = std::max(it.value(), num_matches);
                }
            );
        };
        generate_from_query(query_rc, find_nodes_bwd, true);
        for (size_t end_clipping = 0; end_clipping < nodes.size(); ++end_clipping) {
            for (const auto &[node, seed_length] : nodes[end_clipping]) {
                size_t clipping = this->query_.size() - end_clipping - seed_length;
                std::string_view seed_window(this->query_.data() + clipping,
                                             seed_length);
                size_t num_added = seed_length - this->config_.min_seed_length;
                add_seed(this->query_, clipping + num_added, seed_window, node);
            }
        }
    }

    this->num_matching_ = seeds_.empty() ? 0 : sdsl::util::cnt_one_bits(matched);

    if (this->num_matching_ < this->query_.size() * this->config_.min_exact_match) {
        this->num_matching_ = 0;
        seeds_.clear();
    }
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
