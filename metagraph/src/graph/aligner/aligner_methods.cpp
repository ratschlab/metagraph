#include "aligner_methods.hpp"

#include <unordered_set>
#include <string_view>

#include "dbg_succinct.hpp"
#include "bounded_priority_queue.hpp"
#include "utils.hpp"


template <typename NodeType>
typename Alignment<NodeType>::DPTable
initialize_dp_table(NodeType start_node,
                    char start_char,
                    typename Alignment<NodeType>::score_t initial_score,
                    typename Alignment<NodeType>::score_t min_score,
                    size_t size,
                    int8_t gap_opening_penalty,
                    int8_t gap_extension_penalty) {
    typename Alignment<NodeType>::DPTable dp_table;

    // Initialize first column
    auto& table_init = dp_table[start_node];
    std::get<0>(table_init).resize(size, min_score);
    std::get<1>(table_init).resize(size);
    std::get<2>(table_init) = start_char;
    std::get<3>(table_init) = 0;

    std::get<0>(table_init).front() = initial_score;
    std::get<1>(table_init).front() = std::make_pair(Cigar::Operator::MATCH,
                                                     DeBruijnGraph::npos);

    if (size > 1 && std::get<0>(table_init).front()
                        + gap_opening_penalty > min_score) {
        std::get<0>(table_init)[1] = std::get<0>(table_init).front()
            + gap_opening_penalty;

        std::get<1>(table_init)[1] = std::make_pair(Cigar::Operator::INSERTION,
                                                    start_node);
    }

    for (size_t i = 2; i < size; ++i) {
        if (std::get<0>(table_init)[i - 1] + gap_extension_penalty <= min_score)
            break;

        std::get<0>(table_init)[i] = std::get<0>(table_init)[i - 1]
            + gap_extension_penalty;

        std::get<1>(table_init)[i] = std::make_pair(Cigar::Operator::INSERTION,
                                                    start_node);
    }

    return dp_table;
}

bool early_cutoff(const DeBruijnGraph &graph,
                  const DBGAlignerConfig &config,
                  const char *seed_begin,
                  const char *seed_end) {
    if (seed_begin == seed_end)
        return true;

    if (config.min_seed_length > std::min(graph.get_k(),
                                          static_cast<size_t>(seed_end - seed_begin)))
        return true;

    return false;
}


template <typename NodeType>
std::vector<Alignment<NodeType>> exact_seeder(const DeBruijnGraph &graph,
                                              const DBGAlignerConfig &config,
                                              const char *seed_begin,
                                              const char *seed_end,
                                              size_t clipping,
                                              bool orientation) {
    assert(seed_end >= seed_begin);

    if (early_cutoff(graph, config, seed_begin, seed_end)
            || config.max_seed_length < graph.get_k()
            || static_cast<size_t>(seed_end - seed_begin) < graph.get_k())
        return {};

    auto exact = graph.kmer_to_node(std::string(seed_begin,
                                                seed_begin + graph.get_k()));

    if (exact == DeBruijnGraph::npos)
        return {};

    auto match_score = config.match_score(seed_begin,
                                          seed_begin + graph.get_k());

    if (match_score <= config.min_cell_score)
        return {};

    return { Alignment<NodeType>(seed_begin,
                                 seed_begin + graph.get_k(),
                                 { exact },
                                 match_score,
                                 clipping,
                                 orientation) };
}

template <typename NodeType>
std::vector<Alignment<NodeType>> suffix_seeder(const DeBruijnGraph &graph,
                                               const DBGAlignerConfig &config,
                                               const char *seed_begin,
                                               const char *seed_end,
                                               size_t clipping,
                                               bool orientation) {
    assert(seed_end >= seed_begin);

    if (early_cutoff(graph, config, seed_begin, seed_end))
        return {};

    auto seeds = exact_seeder<NodeType>(graph,
                                        config,
                                        seed_begin,
                                        seed_end,
                                        clipping,
                                        orientation);

    if (seeds.size())
        return seeds;

    const auto* dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
    if (!dbg_succ)
        return {};

    auto max_seed_length = std::min(static_cast<size_t>(seed_end - seed_begin),
                                    std::min(config.max_seed_length, graph.get_k()));
    dbg_succ->call_nodes_with_suffix(
        seed_begin,
        seed_begin + max_seed_length,
        [&](auto node, uint64_t seed_length) {
            assert(node != DeBruijnGraph::npos);

            auto match_score = config.match_score(seed_begin,
                                                  seed_begin + seed_length);

            if (match_score > config.min_cell_score)
                seeds.emplace_back(seed_begin,
                                   seed_begin + seed_length,
                                   std::vector<NodeType>{ node },
                                   match_score,
                                   clipping,
                                   orientation,
                                   graph.get_k() - seed_length);
        },
        config.min_seed_length,
        config.max_num_seeds_per_locus
    );

    return seeds;
}


template <typename NodeType, class Compare>
void default_extender(const DeBruijnGraph &graph,
                      const Alignment<NodeType> &path,
                      std::vector<Alignment<NodeType>>* next_paths,
                      const char* sequence_end,
                      const DBGAlignerConfig &config,
                      bool orientation,
                      typename Alignment<NodeType>::score_t min_path_score) {
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;
    typedef typename DBGAlignment::DPTable DPTable;
    typedef typename DBGAlignment::Step Step;
    typedef typename DBGAlignment::Column Column;

    // this extender only works if at least one character has been matched
    assert(path.get_query_end() > path.get_query_begin());

    const auto align_start = &*path.get_query_end();
    size_t size = sequence_end - align_start + 1;

    // branch and bound
    if (path.get_score() + config.match_score(align_start,
                                              sequence_end) < min_path_score)
        return;

    // used for branch and bound checks below
    std::vector<score_t> partial_sum(size - 1);
    std::string_view extension(align_start, sequence_end - align_start);
    auto jt = extension.rbegin();
    partial_sum.back() = config.get_row(*jt)[*jt++];
    for (auto it = partial_sum.rbegin() + 1; it != partial_sum.rend(); ++it) {
        *it = *(it - 1) + config.get_row(*jt)[*jt++];
    }

    assert(partial_sum.front() == config.match_score(align_start, sequence_end));
    assert(partial_sum.back() == config.get_row(extension.back())[extension.back()]);

    // keep track of which columns to use next
    struct PriorityFunction {
        bool operator()(typename DPTable::value_type* a,
                        typename DPTable::value_type* b) const {
            return compare_(a->second, b->second);
        }

        const Compare compare_ = Compare();
    };

    BoundedPriorityQueue<typename DPTable::value_type*,
                         std::vector<typename DPTable::value_type*>,
                         PriorityFunction> columns_to_update(
        config.queue_size
    );

    auto dp_table = initialize_dp_table(path.back(),
                                        *(align_start - 1),
                                        path.get_score(),
                                        config.min_cell_score,
                                        size,
                                        config.gap_opening_penalty,
                                        config.gap_extension_penalty);
    assert(dp_table.size());

    // keep track of node and position in column to start backtracking
    // store raw pointer since they are not invalidated by emplace
    auto start_node = &*dp_table.begin();

    // for storage of intermediate values
    std::vector<typename DPTable::value_type*> out_columns;
    std::vector<node_index> in_nodes;
    std::vector<std::pair<Step, score_t>> updates(size);
    std::vector<int8_t> char_scores(size);
    std::vector<Cigar::Operator> match_ops(size);
    std::vector<int8_t> gap_scores(size);

    // dynamic programming
    assert(dp_table.find(path.back()) != dp_table.end());
    //columns_to_update.emplace(std::get<0>(dp_table[path.back()]).front(),
    //                          path.back());
    columns_to_update.emplace(&*dp_table.begin());
    while (columns_to_update.size()) {
        //auto cur_node = std::get<1>(columns_to_update.pop_top());
        auto cur_node = columns_to_update.pop_top()->first;
        // get next columns
        out_columns.clear();
        graph.call_outgoing_kmers(
            cur_node,
            [&](auto next_node, char c) {
                auto emplace = dp_table.emplace(
                    next_node,
                    Column { std::vector<score_t>(size, config.min_cell_score),
                             std::vector<Step>(size),
                             c,
                             0 }
                );

                // if the emplace didn't happen, correct the stored character
                if (!emplace.second)
                    std::get<2>(emplace.first->second) = c;

                out_columns.emplace_back(&*emplace.first);
            }
        );

        // update columns
        for (auto &iter : out_columns) {
            auto next_node = iter->first;
            auto& next_column = iter->second;
            auto& next = std::get<0>(next_column);
            auto& next_step = std::get<1>(next_column);
            auto next_char = std::get<2>(next_column);

            // for now, vertical banding is not done
            assert(next.size() == size);
            assert(updates.size() == size);
            std::transform(next_step.begin(),
                           next_step.end(),
                           next.begin(),
                           updates.begin(),
                           [](const Step &step, score_t score) {
                               return std::make_pair(step, score);
                           });

            // speed hack for linear portions of the graph
            if (LIKELY(graph.indegree(next_node) == 1)) {
                in_nodes.assign(1, cur_node);
            } else {
                in_nodes.clear();
                graph.adjacent_incoming_nodes(next_node, [&](auto i) { in_nodes.push_back(i); });
            }

            // match and deletion scores
            for (const auto &prev_node : in_nodes) {
                // the value of the node last character stored here is a
                // placeholder which is later corrected by call_outgoing_kmers above
                auto incoming = dp_table.emplace(
                    prev_node,
                    Column { std::vector<score_t>(size, config.min_cell_score),
                             std::vector<Step>(size),
                             '\0',
                             0 }
                ).first;

                const auto& current = std::get<0>(incoming->second);
                const auto& current_step = std::get<1>(incoming->second);

                // match
                std::transform(align_start,
                               align_start + size - 1,
                               std::next(char_scores.begin()),
                               [row = config.get_row(next_char)](char c) {
                                   return row[c];
                               });

                std::transform(align_start,
                               align_start + size - 1,
                               std::next(match_ops.begin()),
                               [op_row = Cigar::get_op_row(next_char)](char c) {
                                   return op_row[c];
                               });

                for (size_t i = 1; i < size; ++i) {
                    // prevent underflow if min_cell_score == MIN_INT
                    if (char_scores[i] < 0
                            && current[i - 1] == config.min_cell_score)
                        continue;

                    if (current[i - 1] + char_scores[i] > std::get<1>(updates[i]))
                        updates[i] = std::make_pair(Step(match_ops[i], prev_node),
                                                    current[i - 1] + char_scores[i]);
                }


                // delete
                std::transform(
                    current_step.begin(),
                    current_step.end(),
                    gap_scores.begin(),
                    [open = config.gap_opening_penalty,
                     ext = config.gap_extension_penalty](const auto &cigar_tuple) {
                        return std::get<0>(cigar_tuple) == Cigar::Operator::DELETION
                            ? ext : open;
                    }
                );

                for (size_t i = 0; i < size; ++i) {
                    // disabled check for deletion after insertion
                    // if (std::get<0>(current_step[i]) == Cigar::Operator::INSERTION)
                    //     continue;
                    if (current[i] == config.min_cell_score)
                        continue;

                    if (current[i] + gap_scores[i] > std::get<1>(updates[i]))
                        updates[i] = std::make_pair(
                            Step(Cigar::Operator::DELETION, prev_node),
                            current[i] + gap_scores[i]
                        );
                }
            }

            // compute insert scores
            for (size_t i = 1; i < size; ++i) {
                auto prev_op = std::get<0>(std::get<0>(updates[i - 1]));
                // disabled check for insertion after deletion
                // if (prev_op == Cigar::Operator::DELETION)
                //     continue;

                score_t insert_score = std::get<1>(updates[i - 1])
                    + (prev_op == Cigar::Operator::INSERTION
                        ? config.gap_extension_penalty
                        : config.gap_opening_penalty);

                if (insert_score > std::get<1>(updates[i]))
                    updates[i] = std::make_pair(
                        Step(Cigar::Operator::INSERTION, next_node),
                        insert_score
                    );
            }

            // update scores
            bool updated = false;
            auto max_pos = next.begin() + std::get<3>(next_column);
            for (size_t i = 0; i < size; ++i) {
                if (std::get<1>(updates[i]) <= next[i])
                    continue;

                std::tie(next_step[i], next[i]) = updates[i];
                updated = true;
                if (next[i] > *max_pos)
                    max_pos = next.begin() + i;
            }

            if (updated) {
                // store max pos
                auto& cur_best_pos = std::get<3>(next_column);
                cur_best_pos = max_pos - next.begin();

                auto best_score = std::get<0>(start_node->second).at(
                    std::get<3>(start_node->second)
                );

                if (*max_pos > best_score) {
                    start_node = &*iter;
                    best_score = *max_pos;
                }

                // branch and bound

                // starting from the position of the best score, check if each
                // position can be extended to a better alignment than the best
                // alignment
                auto is_not_extendable = [best_score](auto a, auto b) {
                    return a + b < best_score;
                };

                auto check_first_half = std::equal(
                    std::make_reverse_iterator(partial_sum.begin() + cur_best_pos + 1),
                    partial_sum.rbegin(),
                    std::make_reverse_iterator(std::get<0>(next_column).begin()
                                                   + cur_best_pos + 1),
                    is_not_extendable
                );

                if (!check_first_half
                        || !std::equal(partial_sum.begin() + cur_best_pos + 1,
                                       partial_sum.end(),
                                       std::get<0>(next_column).begin() + cur_best_pos + 1,
                                       is_not_extendable))
                    columns_to_update.emplace(&*iter);
            }
        }
    }

    auto best_pos = std::get<3>(start_node->second);
    auto best_score = std::get<0>(start_node->second).at(best_pos);
    assert(best_score > config.min_cell_score);

    //no good path found
    if (UNLIKELY(start_node->first == DeBruijnGraph::npos
            || (start_node->first == path.back() && !best_pos)
            || best_score < min_path_score))
        return;

#ifndef NDEBUG
    // make sure that start_node points to the best column
    auto max_element = std::max_element(
        dp_table.begin(), dp_table.end(),
        [](const auto &a, const auto &b) {
            return std::get<0>(a.second).at(std::get<3>(a.second))
                 < std::get<0>(b.second).at(std::get<3>(b.second));
        }
    );
    assert(best_score
        == std::get<0>(max_element->second).at(std::get<3>(max_element->second)));
#endif

    assert(std::get<0>(std::get<1>(start_node->second).at(best_pos))
        == Cigar::Operator::MATCH);

    // avoid sorting column iterators if we're only interested in the top path
    if (config.num_alternative_paths == 1) {
        next_paths->emplace_back(dp_table,
                                 start_node,
                                 best_pos,
                                 best_score - path.get_score(),
                                 align_start,
                                 orientation);

        if (next_paths->back().empty() && !next_paths->back().get_query_begin()) {
            next_paths->pop_back();
        } else {
            assert(next_paths->back().get_score() + path.get_score() == best_score);
        }

        return;
    }

    // get alternative alignments

    // store visited nodes in paths to avoid returning subalignments
    std::unordered_set<node_index> visited_nodes;

    // store the best score with each iterator in a pair to avoid std::get<3>
    // computations while sorting
    std::vector<std::pair<score_t, const typename DPTable::value_type*>> starts;
    starts.reserve(dp_table.size());
    for (auto it = dp_table.cbegin(); it != dp_table.cend(); ++it) {
        auto max_pos = std::get<3>(it->second);
        auto score = std::get<0>(it->second).at(max_pos);
        auto last_op = std::get<0>(std::get<1>(it->second).at(max_pos));

        if (score > min_path_score && last_op == Cigar::Operator::MATCH)
            starts.emplace_back(score, &*it);
    }

    if (starts.empty())
        return;

    std::sort(starts.begin(), starts.end(),
              [](const auto &a, const auto &b) { return a.first > b.first; });

    assert(best_score == starts.front().first);
    assert(start_node == starts.front().second);
    for (const auto &pair : starts) {
        if (next_paths->size() >= config.num_alternative_paths)
            break;

        auto column_it = pair.second;

        // ignore if the current point is a subalignment of one already visited
        if (visited_nodes.find(column_it->first) != visited_nodes.end())
            continue;

        auto best_pos = std::get<3>(column_it->second);
        auto best_score = std::get<0>(column_it->second).at(best_pos);
        assert(best_score > config.min_cell_score);
        assert(best_score > min_path_score);

        next_paths->emplace_back(dp_table,
                                 column_it,
                                 best_pos,
                                 best_score - path.get_score(),
                                 align_start,
                                 orientation);

        if (next_paths->back().empty() && !next_paths->back().get_query_begin()) {
            next_paths->pop_back();
        } else {
            visited_nodes.insert(next_paths->back().begin(), next_paths->back().end());
        }
    }
}

template <typename NodeType>
Seeder<NodeType>
build_mem_seeder(const std::vector<NodeType> &nodes,
                 std::function<bool(NodeType,
                                    const DeBruijnGraph &)> stop_matching,
                 const DeBruijnGraph &graph) {
    typedef Alignment<NodeType> DBGAlignment;

    auto stop_matching_mapped = [stop_matching, &graph](auto node) {
        return node == DeBruijnGraph::npos || stop_matching(node, graph);
    };

    return [nodes,
            stop_matching_mapped](const DeBruijnGraph &graph,
                                  const DBGAlignerConfig &config,
                                  const char* seed_begin,
                                  const char* seed_end,
                                  size_t clipping,
                                  bool orientation) -> std::vector<DBGAlignment> {
        assert(seed_end >= seed_begin);

        if (seed_begin == seed_end)
            return {};

        assert(nodes.begin() + clipping < nodes.end());
        assert(static_cast<size_t>(seed_end - seed_begin) + clipping
                   == nodes.size() + graph.get_k() - 1);

        auto start = std::find_if(
            nodes.begin() + clipping,
            nodes.end(),
            [](auto node) { return node != DeBruijnGraph::npos; }
        );

        if (start != nodes.begin() + clipping)
            return {};

        assert(*start != DeBruijnGraph::npos);
        auto next = std::find_if(start, nodes.end(), stop_matching_mapped);

        if (UNLIKELY(next != nodes.end() && *next != DeBruijnGraph::npos))
            next++;

        assert(next != start);

        auto end_it = next == nodes.end()
            ? seed_end
            : seed_begin + (next - nodes.begin() - clipping) + graph.get_k() - 1;

        assert(end_it >= seed_begin);
        assert(end_it <= seed_end);
        assert(static_cast<size_t>(end_it - seed_begin) >= graph.get_k());

        auto match_score = config.match_score(seed_begin, end_it);

        if (match_score <= config.min_cell_score)
            return {};

        return { DBGAlignment(seed_begin,
                              end_it,
                              std::vector<NodeType>(start, next),
                              config.match_score(seed_begin, end_it),
                              clipping,
                              orientation) };
    };
}


template <typename NodeType>
Seeder<NodeType> build_unimem_seeder(const std::vector<NodeType> &nodes,
                                     const DeBruijnGraph &graph) {
    return build_mem_seeder<NodeType>(
        nodes,
        [](auto node, const DeBruijnGraph &graph) {
            return graph.outdegree(node) > 1 || graph.indegree(node) > 1;
        },
        graph
    );
}


template std::vector<Alignment<DeBruijnGraph::node_index>>
exact_seeder<DeBruijnGraph::node_index>(
    const DeBruijnGraph&,
    const DBGAlignerConfig&,
    const char*,
    const char*,
    size_t,
    bool
);

template std::vector<Alignment<DeBruijnGraph::node_index>>
suffix_seeder<DeBruijnGraph::node_index>(
    const DeBruijnGraph&,
    const DBGAlignerConfig&,
    const char*,
    const char*,
    size_t,
    bool
);

template void
default_extender<DeBruijnGraph::node_index>(
    const DeBruijnGraph&,
    const Alignment<DeBruijnGraph::node_index>&,
    std::vector<Alignment<DeBruijnGraph::node_index>>*,
    const char*,
    const DBGAlignerConfig&,
    bool,
    typename Alignment<DeBruijnGraph::node_index>::score_t
);

template Seeder<DeBruijnGraph::node_index>
build_mem_seeder(
    const std::vector<DeBruijnGraph::node_index> &,
    std::function<bool(DeBruijnGraph::node_index,
                       const DeBruijnGraph &)>,
    const DeBruijnGraph &
);

template Seeder<DeBruijnGraph::node_index>
build_unimem_seeder(
    const std::vector<DeBruijnGraph::node_index> &,
    const DeBruijnGraph &
);
