#include "aligner_methods.hpp"

#include "common/bounded_priority_queue.hpp"
#include "utils/algorithms.hpp"


/*
 * Helpers for DefaultColumnExtender::operator()
 */

// since insertion invalidates references, return a vector of pointers
template <typename NodeType,
          typename Column = typename DPTable<NodeType>::Column,
          typename score_t = typename DPTable<NodeType>::score_t,
          typename Step = typename DPTable<NodeType>::Step>
std::vector<typename DPTable<NodeType>::value_type*>
get_outgoing_columns(const DeBruijnGraph &graph,
                     DPTable<NodeType> &dp_table,
                     NodeType cur_node,
                     size_t size,
                     size_t best_pos,
                     score_t min_cell_score) {
    std::vector<typename DPTable<NodeType>::value_type*> out_columns;

    graph.call_outgoing_kmers(
        cur_node,
        [&](auto next_node, char c) {
            auto find = dp_table.find(next_node);
            if (find != dp_table.end()) {
                // if the entry was found, correct the stored character
                find->second.last_char = c;
            } else {
                find = dp_table.emplace(
                    next_node,
                    Column { std::vector<score_t>(size, min_cell_score),
                             std::vector<Step>(size),
                             c,
                             best_pos + 1 != size ? best_pos + 1 : best_pos }
                ).first;
            }

            assert(find != dp_table.end());
            out_columns.emplace_back(&*find);
        }
    );

    return out_columns;
}

template <typename NodeType>
std::pair<size_t, size_t>
get_column_boundaries(const DeBruijnGraph &graph,
                      const NodeType &node,
                      DPTable<NodeType> &dp_table,
                      std::vector<NodeType> &in_nodes,
                      size_t size,
                      size_t prev_best_pos,
                      size_t bandwidth) {
    in_nodes.clear();
    graph.adjacent_incoming_nodes(node, [&](auto i) {
        if (dp_table.find(i) != dp_table.end())
            in_nodes.push_back(i);
    });

    size_t overall_begin = size;
    size_t overall_end = 0;

    // find incoming nodes to check for alignment extension
    // set boundaries for vertical band
    for (const auto &prev_node : in_nodes) {
        assert(dp_table.find(prev_node) != dp_table.end());
        size_t best_pos = dp_table.find(prev_node)->second.best_pos;

        if (overall_begin) {
            overall_begin = std::min(
                overall_begin,
                best_pos >= bandwidth ? best_pos - bandwidth : 0
            );
        }

        if (overall_end < size) {
            overall_end = std::max(
                overall_end,
                bandwidth <= size - best_pos ? best_pos + bandwidth : size
            );
        }
    }

    if (overall_begin) {
        overall_begin = std::min(
            overall_begin,
            prev_best_pos >= bandwidth ? prev_best_pos - bandwidth : 0
        );
    }

    if (overall_end < size) {
        overall_end = std::max(
            overall_end,
            bandwidth <= size - prev_best_pos ? prev_best_pos + bandwidth : size
        );
    }

    assert(overall_begin <= prev_best_pos);
    assert(overall_end > prev_best_pos);

    return std::make_pair(overall_begin, overall_end);
}

template <typename ScoreType,
          typename OpType,
          typename ScoreRowType,
          typename OpRowType>
void compute_match_scores(const char *align_begin,
                          const char *align_end,
                          std::vector<ScoreType> &char_scores,
                          std::vector<OpType> &match_ops,
                          const ScoreRowType &row,
                          const OpRowType &op_row) {
    static_assert(std::is_same<ScoreType, typename ScoreRowType::value_type>::value);
    static_assert(std::is_same<OpType, typename OpRowType::value_type>::value);
    static_assert(sizeof(ScoreType) == 1);
    static_assert(sizeof(OpType) == 1);
    assert(align_end >= align_begin);

    // compute match scores
    char_scores.resize(align_end - align_begin + 1);
    match_ops.resize(align_end - align_begin + 1);
    std::transform(align_begin, align_end,
                   char_scores.begin() + 1,
                   [&row](char c) { return row[c]; });

    std::transform(align_begin, align_end,
                   match_ops.begin() + 1,
                   [&op_row](char c) { return op_row[c]; });
}

template <typename UpdateVector,
          typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t,
          typename Step = typename Alignment<NodeType>::Step>
void compute_match_delete_updates(const DBGAlignerConfig &config,
                                  UpdateVector *updates,
                                  const NodeType &prev_node,
                                  const score_t *incoming_scores,
                                  const Step *incoming_steps,
                                  const int8_t *char_scores,
                                  const Cigar::Operator *match_ops,
                                  std::vector<int8_t> &gap_scores,
                                  const score_t min_cell_score,
                                  size_t length) {
    // compute delete scores
    gap_scores.resize(length);
    std::transform(incoming_steps,
                   incoming_steps + length,
                   gap_scores.begin(),
                   [&config](const auto &cigar_tuple) {
                       return cigar_tuple.cigar_op == Cigar::Operator::DELETION
                           ? config.gap_extension_penalty
                           : config.gap_opening_penalty;
                   });

    auto gap_it = gap_scores.begin();

    if (*incoming_scores != min_cell_score
            && *incoming_scores + *gap_it > std::get<1>(*updates)) {
        updates->first.cigar_op = Cigar::Operator::DELETION;
        updates->first.prev_node = prev_node;
        updates->second = *incoming_scores + *gap_it;
    }

    ++incoming_scores;
    ++char_scores;
    ++match_ops;
    ++gap_it;
    ++updates;

    for (size_t i = 1; i < length; ++i) {
        // prevent underflow if min_cell_score == MIN_INT
        if ((*char_scores >= 0 || *(incoming_scores - 1) != min_cell_score)
                && (*(incoming_scores - 1) + *char_scores > std::get<1>(*updates))) {
            updates->first.cigar_op = *match_ops;
            updates->first.prev_node = prev_node;
            updates->second = *(incoming_scores - 1) + *char_scores;
        }

        // TODO: enable check for deletion after insertion?
        if ((*incoming_scores != min_cell_score)
                && (*incoming_scores + *gap_it > std::get<1>(*updates))) {
            updates->first.cigar_op = Cigar::Operator::DELETION;
            updates->first.prev_node = prev_node;
            updates->second = *incoming_scores + *gap_it;
        }

        ++incoming_scores;
        ++char_scores;
        ++match_ops;
        ++gap_it;
        ++updates;
    }
}


/*
 * DefaultColumnExtender::operator()
 */

template <typename NodeType, class Compare>
std::vector<typename DefaultColumnExtender<NodeType, Compare>::DBGAlignment>
DefaultColumnExtender<NodeType, Compare>
::operator()(const DBGAlignment &path,
             const char *sequence_end,
             const score_t *match_score_begin,
             bool orientation,
             score_t min_path_score) const {
    // this extender only works if at least one character has been matched
    assert(path.get_query_end() > path.get_query_begin());
    assert(sequence_end >= path.get_query_end());

    const auto *align_start = path.get_query_end();
    size_t size = sequence_end - align_start + 1;

    assert(config_.match_score(align_start - 1, sequence_end) == *match_score_begin);
    assert(config_.get_row(*(sequence_end - 1))[*(sequence_end - 1)]
        == *(match_score_begin + size - 1));

    // stop path early if it can't be better than the min_path_score
    if (path.get_score() + *(match_score_begin + 1) < min_path_score)
        return {};

    // keep track of which columns to use next
    BoundedPriorityQueue<typename DPTable::value_type*,
                         std::vector<typename DPTable::value_type*>,
                         ColumnPriorityFunction> columns_to_update(config_.queue_size);

    DPTable dp_table(path.back(),
                     *(align_start - 1),
                     path.get_score(),
                     config_.min_cell_score,
                     size,
                     config_.gap_opening_penalty,
                     config_.gap_extension_penalty);

    assert(dp_table.size() == 1);

    // for storage of intermediate values
    std::vector<int8_t> char_scores;
    std::vector<Cigar::Operator> match_ops;
    std::vector<int8_t> gap_scores;
    std::vector<node_index> in_nodes;

    std::vector<std::pair<Step, score_t>> updates;

    // dynamic programming
    // keep track of node and position in column to start backtracking
    // store raw pointer since they are not invalidated by emplace
    auto *start_node = &*dp_table.begin();
    columns_to_update.emplace(start_node);
    while (columns_to_update.size()) {
        const auto* cur_col = columns_to_update.pop_top();
        auto cur_node = cur_col->first;

        // get next columns
        auto out_columns = get_outgoing_columns(graph_,
                                                dp_table,
                                                cur_node,
                                                size,
                                                cur_col->second.best_pos,
                                                config_.min_cell_score);

        // update columns
        for (const auto &iter : out_columns) {
            auto next_node = iter->first;
            auto &next_column = iter->second;

            auto [overall_begin, overall_end] = get_column_boundaries(
                graph_,
                next_node,
                dp_table,
                in_nodes,
                size,
                next_column.best_pos,
                config_.bandwidth
            );

            updates.assign(overall_end - overall_begin, { {}, config_.min_cell_score });

            compute_match_scores(align_start + overall_begin,
                                 align_start + overall_end - 1,
                                 char_scores, match_ops,
                                 config_.get_row(next_column.last_char),
                                 Cigar::get_op_row(next_column.last_char));

            // match and deletion scores
            for (const auto &prev_node : in_nodes) {
                assert(dp_table.find(prev_node) != dp_table.end());
                const auto &incoming = dp_table.find(prev_node)->second;

                size_t begin = incoming.best_pos >= config_.bandwidth
                    ? incoming.best_pos - config_.bandwidth : 0;
                size_t end = config_.bandwidth <= size - incoming.best_pos
                    ? incoming.best_pos + config_.bandwidth : size;

                assert(end > begin);

                compute_match_delete_updates(
                    config_,
                    updates.data() + (begin - overall_begin),
                    prev_node,
                    incoming.scores.data() + begin,
                    incoming.steps.data() + begin,
                    char_scores.data() + (begin - overall_begin),
                    match_ops.data() + (begin - overall_begin),
                    gap_scores,
                    config_.min_cell_score,
                    end - begin
                );
            }

            // compute insert score and update overall scores
            auto max_pos = next_column.scores.begin() + next_column.best_pos;
            bool updated = false;

            if (std::get<1>(updates[0]) > next_column.scores[overall_begin]) {
                std::tie(next_column.steps[overall_begin],
                         next_column.scores[overall_begin]) = updates[0];
                updated = true;

                if (next_column.scores[overall_begin] > *max_pos)
                    max_pos = next_column.scores.begin() + overall_begin;
            }

            for (size_t i = 1; i < updates.size(); ++i) {
                // TODO: check for insertion after deletion?
                score_t insert_score = std::get<1>(updates[i - 1])
                    + (std::get<0>(updates[i - 1]).cigar_op == Cigar::Operator::INSERTION
                        ? config_.gap_extension_penalty
                        : config_.gap_opening_penalty);

                if (insert_score > std::get<1>(updates[i])) {
                    updates[i].first.cigar_op = Cigar::Operator::INSERTION;
                    updates[i].first.prev_node = next_node;
                    updates[i].second = insert_score;
                }

                if (std::get<1>(updates[i]) <= next_column.scores[overall_begin + i])
                    continue;

                std::tie(next_column.steps[overall_begin + i],
                         next_column.scores[overall_begin + i]) = updates[i];

                updated = true;

                if (next_column.scores[overall_begin + i] > *max_pos)
                    max_pos = next_column.scores.begin() + overall_begin + i;
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

                if (!std::equal(match_score_begin + overall_begin,
                                match_score_begin + overall_end,
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
        == std::max_element(dp_table.begin(), dp_table.end(),
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

        if (UNLIKELY(next_paths.back().empty() && !next_paths.back().get_query_begin())) {
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
    return dp_table.extract_alignments(graph_,
                                       config_,
                                       path.get_score(),
                                       align_start,
                                       orientation,
                                       min_path_score);
}


template class DefaultColumnExtender<>;
