#include "aligner_methods.hpp"

#ifdef __AVX2__
#include <emmintrin.h>
#include <immintrin.h>
#endif

#include <Eigen/StdVector>

#include "common/bounded_priority_queue.hpp"
#include "utils/algorithms.hpp"


template <typename T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;


/*
 * Helpers for DefaultColumnExtender::operator()
 */

// since insertion invalidates references, return a vector of pointers
template <typename NodeType,
          typename Column = typename DPTable<NodeType>::Column,
          typename score_t = typename DPTable<NodeType>::score_t>
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
                             std::vector<Cigar::Operator>(size),
                             std::vector<NodeType>(size, DeBruijnGraph::npos),
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

void compute_delete_scores(const Cigar::Operator *incoming_ops_begin,
                           const Cigar::Operator *incoming_ops_end,
                           std::vector<int8_t> &gap_scores,
                           int8_t gap_opening_penalty,
                           int8_t gap_extension_penalty) {
    assert(incoming_ops_end >= incoming_ops_begin);

    // compute delete scores
    gap_scores.resize(incoming_ops_end - incoming_ops_begin);
    std::transform(incoming_ops_begin,
                   incoming_ops_end,
                   gap_scores.begin(),
                   [&](const auto &op) {
                       return op == Cigar::Operator::DELETION
                           ? gap_extension_penalty
                           : gap_opening_penalty;
                   });
}

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t>
void compute_match_delete_updates(score_t *update_scores,
                                  NodeType *update_prevs,
                                  Cigar::Operator *update_ops,
                                  const NodeType &prev_node,
                                  const score_t *incoming_scores,
                                  const int8_t *char_scores,
                                  const Cigar::Operator *match_ops,
                                  const int8_t *gap_it,
                                  size_t length) {
    // handle first element (i.e., no match update possible)
    if (*incoming_scores + *gap_it > *update_scores) {
        *update_ops = Cigar::Operator::DELETION;
        *update_prevs = prev_node;
        *update_scores = *incoming_scores + *gap_it;
    }

    ++incoming_scores;
    ++char_scores;
    ++match_ops;
    ++gap_it;
    ++update_scores;
    ++update_prevs;
    ++update_ops;

    size_t i = 1;

#ifdef __AVX2__
    static_assert(sizeof(Cigar::Operator) == sizeof(int32_t));
    static_assert(sizeof(score_t) == sizeof(int32_t));
    static_assert(sizeof(NodeType) == sizeof(int64_t));

    __m256i prev_packed = _mm256_set1_epi64x(prev_node);
    __m256i del_packed = _mm256_set1_epi32(Cigar::Operator::DELETION);
    for (; i + 8 <= length; i += 8) {
        // compute match and delete scores
        __m256i match_scores_packed = _mm256_add_epi32(
            _mm256_cvtepi8_epi32(_mm_set1_epi64(*(const __m64*)char_scores)),
            _mm256_loadu_si256((__m256i*)(incoming_scores - 1))
        );

        __m256i delete_scores_packed = _mm256_add_epi32(
            _mm256_cvtepi8_epi32(_mm_set1_epi64(*(const __m64*)gap_it)),
            _mm256_loadu_si256((__m256i*)incoming_scores) // TODO: avoid a reload
        );

        // update scores
        __m256i max_scores = _mm256_max_epi32(match_scores_packed, delete_scores_packed);
        __m256i both_cmp = _mm256_cmpgt_epi32(
            max_scores,
            _mm256_loadu_si256((__m256i*)update_scores)
        );
        _mm256_maskstore_epi32(update_scores, both_cmp, max_scores);

        // update ops
        // pick match operation if match score >= gap score
        // pick delete operation if gap score > match score
        _mm256_maskstore_epi32(
            (int32_t*)update_ops,
            both_cmp,
            _mm256_blendv_epi8(_mm256_loadu_si256((__m256i*)match_ops),
                               del_packed,
                               _mm256_cmpgt_epi32(delete_scores_packed,
                                                  match_scores_packed))
        );

        // update prev nodes
        // TODO: this can be done with one AVX512 instruction
        __m128i *cmp_array = (__m128i*)&both_cmp;
        _mm256_maskstore_epi64((long long int*)update_prevs,
                               _mm256_cvtepi32_epi64(cmp_array[0]),
                               prev_packed);
        _mm256_maskstore_epi64((long long int*)update_prevs + 4,
                               _mm256_cvtepi32_epi64(cmp_array[1]),
                               prev_packed);

        update_scores += 8;
        update_prevs += 8;
        update_ops += 8;
        incoming_scores += 8;
        char_scores += 8;
        match_ops += 8;
        gap_it += 8;
    }

    // clean up after AVX2 instructions
    _mm256_zeroupper();
#endif

    // handle rest without AVX2
    for (; i < length; ++i) {
        if (*(incoming_scores - 1) + *char_scores > *update_scores) {
            *update_ops = *match_ops;
            *update_prevs = prev_node;
            *update_scores = *(incoming_scores - 1) + *char_scores;
        }

        // TODO: enable check for deletion after insertion?
        if (*incoming_scores + *gap_it > *update_scores) {
            *update_ops = Cigar::Operator::DELETION;
            *update_prevs = prev_node;
            *update_scores = *incoming_scores + *gap_it;
        }

        ++incoming_scores;
        ++char_scores;
        ++match_ops;
        ++gap_it;
        ++update_scores;
        ++update_prevs;
        ++update_ops;
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

    AlignedVector<score_t> update_scores;
    AlignedVector<Cigar::Operator> update_ops;
    AlignedVector<node_index> update_prevs;

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

            update_scores.assign(overall_end - overall_begin, config_.min_cell_score);
            update_prevs.assign(overall_end - overall_begin, DeBruijnGraph::npos);
            update_ops.assign(overall_end - overall_begin, Cigar::Operator::CLIPPED);

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

                compute_delete_scores(incoming.ops.data() + begin,
                                      incoming.ops.data() + end,
                                      gap_scores,
                                      config_.gap_opening_penalty,
                                      config_.gap_extension_penalty);

                compute_match_delete_updates(
                    update_scores.data() + (begin - overall_begin),
                    update_prevs.data() + (begin - overall_begin),
                    update_ops.data() + (begin - overall_begin),
                    prev_node,
                    incoming.scores.data() + begin,
                    char_scores.data() + (begin - overall_begin),
                    match_ops.data() + (begin - overall_begin),
                    gap_scores.data(),
                    end - begin
                );
            }

            // compute insert score and update overall scores
            auto max_pos = next_column.scores.begin() + next_column.best_pos;
            bool updated = false;

            score_t insert_score;
            for (size_t i = 1; i < update_scores.size(); ++i) {
                // TODO: check for insertion after deletion?
                insert_score = update_scores[i - 1]
                    + (update_ops[i - 1] == Cigar::Operator::INSERTION
                        ? config_.gap_extension_penalty
                        : config_.gap_opening_penalty);

                if (insert_score > update_scores[i]) {
                    update_ops[i] = Cigar::Operator::INSERTION;
                    update_prevs[i] = next_node;
                    update_scores[i] = insert_score;
                }
            }

            size_t i = 0;
#ifdef __AVX2__
            for (; i + 8 <= update_scores.size(); i += 8) {
                // load update scores
                __m256i updates = _mm256_load_si256((__m256i*)&update_scores[i]);

                // compare updates to column
                __m256i cmp = _mm256_cmpgt_epi32(
                    updates,
                    _mm256_loadu_si256((__m256i*)&next_column.scores[overall_begin + i])
                );

                updated |= bool(_mm256_movemask_epi8(cmp));

                // store updates in column
                _mm256_maskstore_epi32(
                    &next_column.scores[overall_begin + i],
                    cmp,
                    updates
                );
                _mm256_maskstore_epi32(
                    (int32_t*)&next_column.ops[overall_begin + i],
                    cmp,
                    _mm256_load_si256((__m256i*)&update_ops[i])
                );

                // TODO: this can be done with one AVX512 instruction
                __m128i *cmp_array = (__m128i*)&cmp;
                _mm256_maskstore_epi64(
                    (long long int*)&next_column.prev_nodes[overall_begin + i],
                    _mm256_cvtepi32_epi64(cmp_array[0]),
                    _mm256_load_si256((__m256i*)&update_prevs[i])
                );
                _mm256_maskstore_epi64(
                    (long long int*)&next_column.prev_nodes[overall_begin + i + 4],
                    _mm256_cvtepi32_epi64(cmp_array[1]),
                    _mm256_load_si256((__m256i*)&update_prevs[i + 4])
                );
            }

            // clean up after AVX2 instructions
            _mm256_zeroupper();
#endif
            for (; i < update_scores.size(); ++i) {
                if (update_scores[i] > next_column.scores[overall_begin + i]) {
                    next_column.ops[overall_begin + i] = update_ops[i];
                    next_column.prev_nodes[overall_begin + i] = update_prevs[i];
                    next_column.scores[overall_begin + i] = update_scores[i];

                    updated = true;
                }
            }

            if (updated) {
                // store max pos
                next_column.best_pos = std::max_element(
                    next_column.scores.begin() + overall_begin,
                    next_column.scores.begin() + overall_end
                ) - next_column.scores.begin();

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
