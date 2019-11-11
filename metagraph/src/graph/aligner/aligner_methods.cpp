#include "aligner_methods.hpp"

#include <unordered_set>
#include <string_view>

#include "dbg_succinct.hpp"
#include "bounded_priority_queue.hpp"


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
std::vector<typename DefaultColumnExtender<NodeType, Compare>::DBGAlignment>
DefaultColumnExtender<NodeType, Compare>
::operator()(const DBGAlignment &path,
             const char* sequence_end,
             bool orientation,
             score_t min_path_score) const {
    // this extender only works if at least one character has been matched
    assert(path.get_query_end() > path.get_query_begin());
    assert(sequence_end >= path.get_query_end());

    const auto *align_start = path.get_query_end();
    size_t size = sequence_end - align_start + 1;

    // compute perfect match scores for all suffixes
    // used for branch and bound checks below
    std::vector<score_t> partial_sum(size);
    std::transform(align_start - 1, sequence_end,
                   partial_sum.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sum.rbegin(), partial_sum.rend(), partial_sum.rbegin());
    assert(config_.match_score(align_start - 1, sequence_end) == partial_sum.front());
    assert(config_.get_row(*(sequence_end - 1))[*(sequence_end - 1)] == partial_sum.back());

    // stop path early if it can't be better than the min_path_score
    if (path.get_score() + partial_sum.at(1) < min_path_score)
        return {};

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
                         PriorityFunction> columns_to_update(config_.queue_size);

    auto dp_table = initialize_dp_table(path.back(),
                                        *(align_start - 1),
                                        path.get_score(),
                                        config_.min_cell_score,
                                        size,
                                        config_.gap_opening_penalty,
                                        config_.gap_extension_penalty);
    assert(dp_table.size() == 1);
    assert(dp_table.find(path.back()) != dp_table.end());

    // for storage of intermediate values
    std::vector<typename DPTable::value_type*> out_columns;
    std::vector<node_index> in_nodes;
    std::vector<std::pair<Step, score_t>> updates;
    std::vector<int8_t> char_scores;
    std::vector<Cigar::Operator> match_ops;
    std::vector<int8_t> gap_scores;

    // dynamic programming
    // keep track of node and position in column to start backtracking
    // store raw pointer since they are not invalidated by emplace
    auto *start_node = &*dp_table.begin();
    columns_to_update.emplace(start_node);
    while (columns_to_update.size()) {
        const auto cur_col = columns_to_update.pop_top();
        auto cur_node = cur_col->first;
        // get next columns
        out_columns.clear();
        graph_.call_outgoing_kmers(
            cur_node,
            [&](auto next_node, char c) {
                auto emplace = dp_table.emplace(
                    next_node,
                    Column { std::vector<score_t>(size, config_.min_cell_score),
                             std::vector<Step>(size),
                             c,
                             cur_col->second.best_pos + 1 != size
                                ? cur_col->second.best_pos + 1
                                : cur_col->second.best_pos }
                );

                // if the emplace didn't happen, correct the stored character
                if (!emplace.second)
                    emplace.first->second.last_char = c;

                assert(emplace.first != dp_table.end());
                out_columns.emplace_back(&*emplace.first);
            }
        );

        // update columns
        for (auto iter : out_columns) {
            auto next_node = iter->first;
            auto& next_column = iter->second;

            // speed hack for linear portions of the graph
            if (LIKELY(graph_.indegree(next_node) == 1)) {
                in_nodes.assign(1, cur_node);
            } else {
                in_nodes.clear();
                graph_.adjacent_incoming_nodes(next_node,
                                               [&](auto i) { in_nodes.push_back(i); });
            }

            size_t overall_begin = size;
            size_t overall_end = 0;

            // find incoming nodes to check for alignment extension
            // set boundaries for vertical band
            for (const auto &prev_node : in_nodes) {
                const auto& best_pos = dp_table.emplace(
                    prev_node,
                    Column { std::vector<score_t>(size, config_.min_cell_score),
                             std::vector<Step>(size),
                             '\0',
                             next_column.best_pos ? next_column.best_pos - 1 : 0 }
                ).first->second.best_pos;

                if (overall_begin)
                    overall_begin = std::min(
                        overall_begin,
                        best_pos >= config_.bandwidth ? best_pos - config_.bandwidth : 0
                    );

                if (overall_end < size)
                    overall_end = std::max(
                        overall_end,
                        config_.bandwidth <= size - best_pos
                            ? best_pos + config_.bandwidth
                            : size
                    );
            }

            if (overall_begin)
                overall_begin = std::min(
                    overall_begin,
                    next_column.best_pos >= config_.bandwidth
                        ? next_column.best_pos - config_.bandwidth
                        : 0
                );

            if (overall_end < size)
                overall_end = std::max(
                    overall_end,
                    config_.bandwidth <= size - next_column.best_pos
                        ? next_column.best_pos + config_.bandwidth
                        : size
                );

            if (overall_begin >= overall_end)
                return {};

            assert(overall_begin <= next_column.best_pos);
            assert(overall_end > next_column.best_pos);

            updates.resize(overall_end - overall_begin);
            for (auto &u : updates) {
                u.first.prev_node = DeBruijnGraph::npos;
                u.second = config_.min_cell_score;
            }

            // match and deletion scores
            for (const auto &prev_node : in_nodes) {
                // the value of the node last character stored here is a
                // placeholder which is later corrected by call_outgoing_kmers above
                const auto& incoming = dp_table[prev_node];

                size_t begin = incoming.best_pos >= config_.bandwidth
                    ? incoming.best_pos - config_.bandwidth : 0;
                size_t end = config_.bandwidth <= size - incoming.best_pos
                    ? incoming.best_pos + config_.bandwidth : size;

                // match
                assert(end);

                char_scores.resize(end - begin);
                const auto &row = config_.get_row(next_column.last_char);
                std::transform(align_start + begin,
                               align_start + end - 1,
                               char_scores.begin() + 1,
                               [&row](char c) { return row[c]; });

                match_ops.resize(end - begin);
                const auto &op_row = Cigar::get_op_row(next_column.last_char);
                std::transform(align_start + begin,
                               align_start + end - 1,
                               match_ops.begin() + 1,
                               [&op_row](char c) { return op_row[c]; });

                for (size_t i = begin + 1; i < end; ++i) {
                    // prevent underflow if min_cell_score == MIN_INT
                    if (char_scores[i - begin] < 0
                            && incoming.scores[i - 1] == config_.min_cell_score)
                        continue;

                    if (incoming.scores[i - 1] + char_scores[i - begin]
                            > std::get<1>(updates[i - overall_begin]))
                        updates[i - overall_begin] = std::make_pair(
                            Step { match_ops[i - begin], prev_node },
                            incoming.scores[i - 1] + char_scores[i - begin]
                        );
                }


                // delete
                gap_scores.resize(end - begin);
                std::transform(
                    incoming.steps.begin() + begin,
                    incoming.steps.begin() + end,
                    gap_scores.begin(),
                    [this](const auto &cigar_tuple) {
                        return cigar_tuple.cigar_op == Cigar::Operator::DELETION
                            ? config_.gap_extension_penalty
                            : config_.gap_opening_penalty;
                    }
                );

                for (size_t i = begin; i < end; ++i) {
                    // disabled check for deletion after insertion
                    // if (incoming.steps[i].cigar_op == Cigar::Operator::INSERTION)
                    //     continue;
                    if (incoming.scores[i] == config_.min_cell_score)
                        continue;

                    if (incoming.scores[i] + gap_scores[i - begin]
                            > std::get<1>(updates[i - overall_begin])) {
                        updates[i - overall_begin].second
                            = incoming.scores[i] + gap_scores[i - begin];
                        updates[i - overall_begin].first.cigar_op
                            = Cigar::Operator::DELETION;
                        updates[i - overall_begin].first.prev_node = prev_node;
                    }
                }
            }

            // compute insert scores
            for (size_t i = 1; i < overall_end - overall_begin; ++i) {
                auto prev_op = std::get<0>(updates[i - 1]).cigar_op;
                // disabled check for insertion after deletion
                // if (prev_op == Cigar::Operator::DELETION)
                //     continue;

                score_t insert_score = std::get<1>(updates[i - 1])
                    + (prev_op == Cigar::Operator::INSERTION
                        ? config_.gap_extension_penalty
                        : config_.gap_opening_penalty);

                if (insert_score > std::get<1>(updates[i])) {
                    updates[i].second = insert_score;
                    updates[i].first.cigar_op = Cigar::Operator::INSERTION;
                    updates[i].first.prev_node = next_node;
                }
            }

            // update scores
            bool updated = false;
            auto max_pos = next_column.scores.begin() + next_column.best_pos;
            for (size_t i = overall_begin; i < overall_end; ++i) {
                if (std::get<1>(updates[i - overall_begin]) <= next_column.scores[i])
                    continue;

                std::tie(next_column.steps[i],
                         next_column.scores[i]) = updates[i - overall_begin];
                updated = true;

                if (next_column.scores[i] > *max_pos)
                    max_pos = next_column.scores.begin() + i;
            }

            if (updated) {
                // store max pos
                next_column.best_pos = max_pos - next_column.scores.begin();

                if (*max_pos > start_node->second.best_score())
                    start_node = iter;

                // branch and bound
                // TODO: this cuts off too early (before the scores have converged)
                //       so the code below has to be used to compute correct scores
                auto best_score = std::max(start_node->second.best_score(),
                                           min_path_score);

                if (!std::equal(partial_sum.begin() + overall_begin,
                                partial_sum.begin() + overall_end,
                                next_column.scores.begin() + overall_begin,
                                [&](auto a, auto b) { return a + b < best_score; }))
                    columns_to_update.emplace(iter);
            }
        }
    }

    assert(start_node->second.best_score() > config_.min_cell_score);

    //no good path found
    if (UNLIKELY(start_node->first == DeBruijnGraph::npos
            || (start_node->first == path.back() && !start_node->second.best_pos)
            || start_node->second.best_score() < min_path_score))
        return {};

    // check to make sure that start_node stores the best starting point
    assert(start_node->second.best_score()
        == std::max_element(dp_table.begin(),
                            dp_table.end(),
                            [](const auto &a, const auto &b) {
                                return a.second < b.second;
                            })->second.best_score());

    assert(start_node->second.best_step().cigar_op == Cigar::Operator::MATCH);

    // avoid sorting column iterators if we're only interested in the top path
    std::vector<DBGAlignment> next_paths;
    if (config_.num_alternative_paths == 1) {
        next_paths.emplace_back(dp_table,
                                start_node,
                                start_node->second.best_pos,
                                start_node->second.best_score() - path.get_score(),
                                align_start,
                                orientation,
                                graph_.get_k() - 1);

        if (next_paths.back().empty() && !next_paths.back().get_query_begin()) {
            next_paths.pop_back();
        } else {
            assert(next_paths.back().get_score() + path.get_score()
                == start_node->second.best_score());

            // TODO: remove this when the branch and bound is set to only consider
            //       converged scores
            next_paths.back().recompute_score(config_);

            assert(next_paths.back().is_valid(graph_, &config_));
        }

        return next_paths;
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
        return {};

    std::sort(starts.begin(), starts.end(),
              [](const auto &a, const auto &b) {
                  return a->second.best_score() > b->second.best_score();
              });

    assert(start_node == starts.front());
    for (const auto &column_it : starts) {
        if (next_paths.size() >= config_.num_alternative_paths)
            break;

        // ignore if the current point is a subalignment of one already visited
        if (visited_nodes.find(column_it->first) != visited_nodes.end())
            continue;

        next_paths.emplace_back(dp_table,
                                column_it,
                                column_it->second.best_pos,
                                column_it->second.best_score() - path.get_score(),
                                align_start,
                                orientation,
                                graph_.get_k() - 1);

        if (next_paths.back().empty() && !next_paths.back().get_query_begin()) {
            next_paths.pop_back();
        } else {
            visited_nodes.insert(next_paths.back().begin(), next_paths.back().end());

            // TODO: remove this when the branch and bound is set to only consider
            //       converged scores
            next_paths.back().recompute_score(config_);

            assert(next_paths.back().is_valid(graph_, &config_));
        }
    }

    return next_paths;
}

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

template class ExactSeeder<>;
template class SuffixSeeder<>;
template class MEMSeeder<>;
template class UniMEMSeeder<>;
template class DefaultColumnExtender<>;
