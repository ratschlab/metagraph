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
    auto& table_init = dp_table.emplace(
        start_node,
        typename Alignment<NodeType>::Column {
            std::vector<typename Alignment<NodeType>::score_t>(size, min_score),
            std::vector<typename Alignment<NodeType>::Step>(size),
            start_char,
            0
        }
    ).first->second;

    table_init.scores.front() = initial_score;
    table_init.steps.front().cigar_op = Cigar::Operator::MATCH;
    table_init.steps.front().prev_node = DeBruijnGraph::npos;

    if (size > 1 && table_init.scores.front() + gap_opening_penalty > min_score) {
        table_init.scores[1] = table_init.scores.front() + gap_opening_penalty;
        table_init.steps[1].cigar_op = Cigar::Operator::INSERTION;
        table_init.steps[1].prev_node = start_node;
    }

    for (size_t i = 2; i < size; ++i) {
        if (table_init.scores[i - 1] + gap_extension_penalty <= min_score)
            break;

        table_init.scores[i] = table_init.scores[i - 1] + gap_extension_penalty;
        table_init.steps[i].cigar_op = Cigar::Operator::INSERTION;
        table_init.steps[i].prev_node = start_node;
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
std::vector<Alignment<NodeType>> ExactSeeder<NodeType>
::operator()(const char *seed_begin,
             const char *seed_end,
             size_t clipping,
             bool orientation) const {
    assert(seed_end >= seed_begin);

    if (early_cutoff(graph_, config_, seed_begin, seed_end)
            || config_.max_seed_length < graph_.get_k()
            || static_cast<size_t>(seed_end - seed_begin) < graph_.get_k())
        return {};

    auto exact = graph_.kmer_to_node(std::string(seed_begin,
                                                 seed_begin + graph_.get_k()));

    if (exact == DeBruijnGraph::npos)
        return {};

    auto match_score = config_.match_score(seed_begin,
                                           seed_begin + graph_.get_k());

    if (match_score <= config_.min_cell_score)
        return {};

    return { Alignment<NodeType>(seed_begin,
                                 seed_begin + graph_.get_k(),
                                 { exact },
                                 match_score,
                                 clipping,
                                 orientation) };
}

template <typename NodeType>
std::vector<Alignment<NodeType>> SuffixSeeder<NodeType>
::operator()(const char *seed_begin,
             const char *seed_end,
             size_t clipping,
             bool orientation) const {
    if (!dynamic_cast<const DBGSuccinct*>(&get_graph()))
        throw std::runtime_error("Only implemented for DBGSuccinct");

    const auto &graph = get_graph();
    const auto &config = get_config();

    assert(seed_end >= seed_begin);

    if (early_cutoff(graph, config, seed_begin, seed_end))
        return {};

    auto seeds = exact_seeder_(seed_begin, seed_end, clipping, orientation);

    if (seeds.size())
        return seeds;

    auto max_seed_length = std::min(static_cast<size_t>(seed_end - seed_begin),
                                    std::min(config.max_seed_length, graph.get_k()));

    const auto &dbg_succ = dynamic_cast<const DBGSuccinct&>(graph);
    dbg_succ.call_nodes_with_suffix(
        seed_begin,
        seed_begin + max_seed_length,
        [&](auto node, uint64_t seed_length) {
            assert(node != DeBruijnGraph::npos);

            auto match_score = config.match_score(seed_begin,
                                                  seed_begin + seed_length);

            if (match_score > get_config().min_cell_score)
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

    assert(seeds.size() <= config.max_num_seeds_per_locus);

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

    // stop path early if it can't be better than the min_path_score
    if (path.get_score() + config.match_score(align_start,
                                              sequence_end) < min_path_score)
        return;

    // used for branch and bound checks below
    std::vector<score_t> partial_sum(size);
    std::string_view extension(align_start - 1, sequence_end - align_start + 1);
    auto jt = extension.rbegin();
    partial_sum.back() = config.get_row(*jt)[*jt++];
    for (auto it = partial_sum.rbegin() + 1; it != partial_sum.rend(); ++it) {
        assert(jt != extension.rend());
        *it = *(it - 1) + config.get_row(*jt)[*jt++];
    }

    assert(partial_sum.front() == config.match_score(align_start - 1, sequence_end));
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
    columns_to_update.emplace(&*dp_table.begin());
    while (columns_to_update.size()) {
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
                    emplace.first->second.last_char = c;

                out_columns.emplace_back(&*emplace.first);
            }
        );

        // update columns
        for (auto &iter : out_columns) {
            auto next_node = iter->first;
            auto& next_column = iter->second;

            // for now, vertical banding is not done
            assert(next_column.scores.size() == size);
            assert(updates.size() == size);
            std::transform(next_column.steps.begin(),
                           next_column.steps.end(),
                           next_column.scores.begin(),
                           updates.begin(),
                           [](const Step &step, score_t score) {
                               return std::make_pair(step, score);
                           });

            // speed hack for linear portions of the graph
            if (LIKELY(graph.indegree(next_node) == 1)) {
                in_nodes.assign(1, cur_node);
            } else {
                in_nodes.clear();
                graph.adjacent_incoming_nodes(next_node,
                                              [&](auto i) { in_nodes.push_back(i); });
            }

            // match and deletion scores
            for (const auto &prev_node : in_nodes) {
                // the value of the node last character stored here is a
                // placeholder which is later corrected by call_outgoing_kmers above
                auto& incoming = dp_table.emplace(
                    prev_node,
                    Column { std::vector<score_t>(size, config.min_cell_score),
                             std::vector<Step>(size),
                             '\0',
                             0 }
                ).first->second;

                // match
                std::transform(align_start,
                               align_start + size - 1,
                               std::next(char_scores.begin()),
                               [row = config.get_row(next_column.last_char)](char c) {
                                   return row[c];
                               });

                std::transform(align_start,
                               align_start + size - 1,
                               std::next(match_ops.begin()),
                               [op_row = Cigar::get_op_row(next_column.last_char)](char c) {
                                   return op_row[c];
                               });

                for (size_t i = 1; i < size; ++i) {
                    // prevent underflow if min_cell_score == MIN_INT
                    if (char_scores[i] < 0
                            && incoming.scores[i - 1] == config.min_cell_score)
                        continue;

                    if (incoming.scores[i - 1] + char_scores[i] > std::get<1>(updates[i]))
                        updates[i] = std::make_pair(
                            Step { match_ops[i], prev_node },
                            incoming.scores[i - 1] + char_scores[i]
                        );
                }


                // delete
                std::transform(
                    incoming.steps.begin(),
                    incoming.steps.end(),
                    gap_scores.begin(),
                    [open = config.gap_opening_penalty,
                     ext = config.gap_extension_penalty](const auto &cigar_tuple) {
                        return cigar_tuple.cigar_op == Cigar::Operator::DELETION
                            ? ext : open;
                    }
                );

                for (size_t i = 0; i < size; ++i) {
                    // disabled check for deletion after insertion
                    // if (incoming.steps[i].cigar_op == Cigar::Operator::INSERTION)
                    //     continue;
                    if (incoming.scores[i] == config.min_cell_score)
                        continue;

                    if (incoming.scores[i] + gap_scores[i] > std::get<1>(updates[i]))
                        updates[i] = std::make_pair(
                            Step { Cigar::Operator::DELETION, prev_node },
                            incoming.scores[i] + gap_scores[i]
                        );
                }
            }

            // compute insert scores
            for (size_t i = 1; i < size; ++i) {
                auto prev_op = std::get<0>(updates[i - 1]).cigar_op;
                // disabled check for insertion after deletion
                // if (prev_op == Cigar::Operator::DELETION)
                //     continue;

                score_t insert_score = std::get<1>(updates[i - 1])
                    + (prev_op == Cigar::Operator::INSERTION
                        ? config.gap_extension_penalty
                        : config.gap_opening_penalty);

                if (insert_score > std::get<1>(updates[i]))
                    updates[i] = std::make_pair(
                        Step { Cigar::Operator::INSERTION, next_node },
                        insert_score
                    );
            }

            // update scores
            bool updated = false;
            auto max_pos = next_column.scores.begin() + next_column.best_pos;
            for (size_t i = 0; i < size; ++i) {
                if (std::get<1>(updates[i]) <= next_column.scores[i])
                    continue;

                std::tie(next_column.steps[i], next_column.scores[i]) = updates[i];
                updated = true;
                if (next_column.scores[i] > *max_pos)
                    max_pos = next_column.scores.begin() + i;
            }

            if (updated) {
                // store max pos
                next_column.best_pos = max_pos - next_column.scores.begin();

                if (*max_pos > start_node->second.best_score())
                    start_node = &*iter;

                // branch and bound

                // starting from the position of the best score, check if each
                // position can be extended to a better alignment than the best
                // alignment
                auto is_not_extendable =
                    [best_score = std::max(start_node->second.best_score(),
                                           min_path_score)](auto a, auto b) {
                        return a + b < best_score;
                    };

                auto start_it = partial_sum.begin() + next_column.best_pos + 1;
                auto next_column_start_it = next_column.scores.begin()
                    + next_column.best_pos + 1;

                // the first value of partial_sum includes the last character
                // of the seed, so it should be excluded from this check
                auto first_half_end_it = next_column.best_pos >= config.bandwidth
                    ? std::make_reverse_iterator(start_it - config.bandwidth)
                    : partial_sum.rend() - 1;

                assert(first_half_end_it != partial_sum.rend());

                auto second_half_end_it = size - next_column.best_pos > config.bandwidth
                    ? start_it + config.bandwidth
                    : partial_sum.end();


                if (!std::equal(std::make_reverse_iterator(start_it),
                                first_half_end_it,
                                std::make_reverse_iterator(next_column_start_it),
                                is_not_extendable)
                        || !std::equal(start_it,
                                       second_half_end_it,
                                       next_column_start_it,
                                       is_not_extendable))
                    columns_to_update.emplace(&*iter);
            }
        }
    }

    assert(start_node->second.best_score() > config.min_cell_score);

    //no good path found
    if (UNLIKELY(start_node->first == DeBruijnGraph::npos
            || (start_node->first == path.back() && !start_node->second.best_pos)
            || start_node->second.best_score() < min_path_score))
        return;

    // check to make sure that start_node stores the best starting point
    assert(start_node->second.best_score()
        == std::max_element(dp_table.begin(),
                            dp_table.end(),
                            [](const auto &a, const auto &b) {
                                return a.second < b.second;
                            })->second.best_score());

    assert(start_node->second.best_step().cigar_op == Cigar::Operator::MATCH);

    // avoid sorting column iterators if we're only interested in the top path
    if (config.num_alternative_paths == 1) {
        next_paths->emplace_back(dp_table,
                                 start_node,
                                 start_node->second.best_pos,
                                 start_node->second.best_score() - path.get_score(),
                                 align_start,
                                 orientation);

        if (next_paths->back().empty() && !next_paths->back().get_query_begin()) {
            next_paths->pop_back();
        } else {
            assert(next_paths->back().get_score() + path.get_score()
                == start_node->second.best_score());
        }

        return;
    }

    // get alternative alignments

    // store visited nodes in paths to avoid returning subalignments
    std::unordered_set<node_index> visited_nodes;

    std::vector<const typename DPTable::value_type*> starts;
    starts.reserve(dp_table.size());
    for (auto it = dp_table.cbegin(); it != dp_table.cend(); ++it) {
        if (it->second.best_score() > min_path_score
                && it->second.best_step().cigar_op == Cigar::Operator::MATCH)
            starts.emplace_back(&*it);
    }

    if (starts.empty())
        return;

    std::sort(starts.begin(), starts.end(),
              [](const auto &a, const auto &b) {
                  return a->second.best_score() > b->second.best_score();
              });

    assert(start_node == starts.front());
    for (const auto &column_it : starts) {
        if (next_paths->size() >= config.num_alternative_paths)
            break;

        // ignore if the current point is a subalignment of one already visited
        if (visited_nodes.find(column_it->first) != visited_nodes.end())
            continue;

        next_paths->emplace_back(dp_table,
                                 column_it,
                                 column_it->second.best_pos,
                                 column_it->second.best_score() - path.get_score(),
                                 align_start,
                                 orientation);

        if (next_paths->back().empty() && !next_paths->back().get_query_begin()) {
            next_paths->pop_back();
        } else {
            visited_nodes.insert(next_paths->back().begin(), next_paths->back().end());
        }
    }
}

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

template <typename NodeType>
std::vector<Alignment<NodeType>> MEMSeeder<NodeType>
::operator()(const char *seed_begin,
             const char *seed_end,
             size_t clipping,
             bool orientation) const {
    if (query_nodes_.empty())
        throw std::runtime_error("MEMSeeder uninitialized");

    if (orientation != orientation_)
        throw std::runtime_error("wrong orientation passed");

    assert(seed_end >= seed_begin);

    if (seed_begin == seed_end)
        return {};

    assert(static_cast<size_t>(seed_end - seed_begin) + clipping
               == query_nodes_.size() + graph_.get_k() - 1);

    auto start = std::find_if(
        query_nodes_.begin() + clipping,
        query_nodes_.end(),
        [](auto node) { return node != DeBruijnGraph::npos; }
    );

    if (start != query_nodes_.begin() + clipping)
        return {};

    assert(*start != DeBruijnGraph::npos);
    auto next = std::find_if(start,
                             query_nodes_.end(),
                             [&](auto i) { return i == DeBruijnGraph::npos
                                               || (*is_mem_terminus_)[i]; });

    if (UNLIKELY(next != query_nodes_.end() && *next != DeBruijnGraph::npos))
        next++;

    assert(next != start);

    auto end_it = next == query_nodes_.end()
        ? seed_end
        : seed_begin + (next - query_nodes_.begin() - clipping) + graph_.get_k() - 1;

    assert(end_it >= seed_begin);
    assert(end_it <= seed_end);
    assert(static_cast<size_t>(end_it - seed_begin) >= graph_.get_k());

    auto match_score = config_.match_score(seed_begin, end_it);

    if (match_score <= config_.min_cell_score)
        return {};

    return { Alignment<NodeType>(seed_begin,
                                 end_it,
                                 std::vector<NodeType>(start, next),
                                 config_.match_score(seed_begin, end_it),
                                 clipping,
                                 orientation) };
}

template class ExactSeeder<DeBruijnGraph::node_index>;
template class SuffixSeeder<DeBruijnGraph::node_index>;
template class MEMSeeder<DeBruijnGraph::node_index>;
template class UniMEMSeeder<DeBruijnGraph::node_index>;
