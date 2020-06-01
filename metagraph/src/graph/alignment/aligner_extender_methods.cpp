#include "aligner_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "common/logger.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/vectors/aligned_vector.hpp"

using mtg::common::logger;


template <typename NodeType>
void DefaultColumnExtender<NodeType>
::initialize_query(const std::string_view query) {
    this->query = query;
    partial_sums_.resize(query.size());
    std::transform(query.begin(), query.end(),
                   partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query) == partial_sums_.front());
    assert(config_.get_row(query.back())[query.back()] == partial_sums_.back());

    profile_score.clear();
    profile_op.clear();
    for (char c : graph_->alphabet()) {
        auto &profile_score_row = profile_score.emplace(c, query.size() + 8).first.value();
        auto &profile_op_row = profile_op.emplace(c, query.size() + 8).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = Cigar::get_op_row(c);

#ifdef __AVX2__

        __m256i *score_cast = (__m256i*)profile_score_row.data();
        __m256i *op_cast = (__m256i*)profile_op_row.data();
        for (size_t i = 0; i < query.size(); i += 8) {
            __m256i str_p = expandepu8_epi32(*(uint64_t*)&query[i]);
            _mm256_store_si256(score_cast++, _mm256_i32gather_epi32(row.data(), str_p, 4));
            _mm256_store_si256(op_cast++, _mm256_i32gather_epi32((int*)op_row.data(), str_p, 4));
        }

#else

        std::transform(query.begin(), query.end(), profile_score_row.begin(),
                       [&row](char q) { return row[q]; });
        std::transform(query.begin(), query.end(), profile_op_row.begin(),
                       [&op_row](char q) { return op_row[q]; });

#endif

    }
}

std::pair<size_t, size_t> get_band(size_t size, size_t best_pos, size_t bandwidth) {
    assert(best_pos < size);
    auto begin = best_pos >= bandwidth ? best_pos - bandwidth : 0;

    // align begin + 1 to 32-byte boundary
    if (begin > 7)
        begin = (begin & 0xFFFFFFFFFFFFFFF8) - 1;

    auto end = bandwidth <= size - best_pos ? best_pos + bandwidth : size;
    end = std::min(begin + ((end - begin + 7) / 8) * 8, size);

    assert(begin <= best_pos);
    assert(end > best_pos);
    assert(begin < end);

    return std::make_pair(begin, end);
}

template <typename NodeType>
std::pair<typename DPTable<NodeType>::iterator, bool> DefaultColumnExtender<NodeType>
::emplace_node(NodeType node, NodeType, char c, size_t size,
               size_t best_pos, size_t last_priority_pos) {
    auto find = dp_table.find(node);
    if (find == dp_table.end()) {
        return dp_table.emplace(node,
                                typename DPTable<NodeType>::Column(
                                    size, config_.min_cell_score, c,
                                    best_pos + 1 != size ? best_pos + 1 : best_pos,
                                    last_priority_pos
                                ));
    } else {
        auto [node_begin, node_end] = get_band(size,
                                               find->second.best_pos,
                                               config_.bandwidth);

        overlapping_range_ |= ((begin <= node_begin && end > node_begin)
            || (begin < node_end && end >= node_end));
    }

    return std::make_pair(find, false);
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>
::add_seed(size_t clipping) {
    assert(path_->get_cigar().back().first == Cigar::Operator::MATCH
        || path_->get_cigar().back().first == Cigar::Operator::MISMATCH);

    return dp_table.add_seed(start_node, *(align_start - 1),
                             config_.get_row(*(align_start - 1))[path_->get_sequence().back()],
                             path_->get_score(), config_.min_cell_score, size, 0,
                             config_.gap_opening_penalty, config_.gap_extension_penalty,
                             clipping);
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    // this extender only works if at least one character has been matched
    assert(path.get_query_end() > path.get_query().data());
    assert(path.get_query_end() > query.data());
    assert(query.data() + query.size() > path.get_query_end());

    align_start = path.get_query_end();
    size = query.data() + query.size() - align_start + 1;
    match_score_begin = partial_sums_.data() + (align_start - 1 - query.data());

    assert(config_.match_score(std::string_view(align_start - 1, size))
        == *match_score_begin);
    assert(config_.get_row(query.back())[query.back()] == match_score_begin[size - 1]);

    start_node = path.back();
    this->path_ = &path;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::check_and_push(ColumnRef&& next_column) {
    const auto &[next_node, best_score_update, converged] = next_column;

    // always push the next column if it hasn't converged
    if (!converged) {
        columns_to_update.emplace(std::move(next_column));
        return;
    }

    assert(dp_table.find(next_node) != dp_table.end());

    // TODO: does the procedure below provably ensure that non-converged columns
    // are not dropped?
    const auto &column = dp_table.find(next_node)->second;

    // Ignore if there is no way it can be extended to an optimal alignment.
    if (xdrop_cutoff - column.best_score() > config_.xdrop
            || std::equal(match_score_begin + begin,
                          match_score_begin + end,
                          column.scores.data() + begin,
                          [&](auto a, auto b) { return a + b < score_cutoff; })) {
        return;
    }

    // if the queue has space, push the next column
    if (columns_to_update.size() < config_.queue_size) {
        columns_to_update.emplace(std::move(next_column));
        return;
    }

    const ColumnRef &bottom = columns_to_update.minimum();

    if (!utils::LessSecond()(bottom, next_column))
        return;

    if (std::get<2>(bottom) || std::get<1>(bottom)
            != dp_table.find(std::get<0>(bottom))->second.last_priority_value()) {
        // if the bottom has converged, or it is an invalidated reference
        // (it's last priority value has changed), then replace the bottom element
        columns_to_update.update(columns_to_update.begin(), std::move(next_column));
        return;
    } else {
        // otherwise, push
        columns_to_update.emplace(std::move(next_column));
        return;
    }
}


/*
 * Helpers for column score updating
 */

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t>
inline void update_del_scores(const DBGAlignerConfig &config,
                              score_t *update_scores,
                              NodeType *update_prevs,
                              Cigar::Operator *update_ops,
                              const NodeType &node,
                              int32_t *updated_mask,
                              size_t length) {
    for (size_t i = 1; i < length; ++i) {
        score_t del_score = std::max(config.min_cell_score,
            update_scores[i - 1] + (update_ops[i - 1] == Cigar::Operator::DELETION
                ? config.gap_extension_penalty
                : config.gap_opening_penalty
        ));

        if (del_score > update_scores[i]) {
            while (i < length && del_score > update_scores[i]) {
                update_scores[i] = del_score;
                update_ops[i] = Cigar::Operator::DELETION;
                update_prevs[i] = node;

                if (updated_mask)
                    updated_mask[i] = 0xFFFFFFFF;

                del_score += config.gap_extension_penalty;
                ++i;
            }
            --i;
        }
    }
}

#ifdef __AVX2__

inline void compute_HE_avx2(size_t length,
                            __m256i prev_node,
                            __m256i gap_opening_penalty,
                            __m256i gap_extension_penalty,
                            int32_t *update_scores,
                            int32_t *update_gap_scores,
                            long long int *update_prevs,
                            int32_t *update_ops,
                            const int32_t *incoming_scores,
                            const int32_t *incoming_gap_scores,
                            const int32_t *incoming_ops,
                            const int32_t *profile_scores,
                            const int32_t *profile_ops,
                            int32_t *updated_mask,
                            __m256i min_cell_score) {
    assert(update_scores != incoming_scores);
    assert(update_gap_scores != incoming_gap_scores);
    assert(update_ops != incoming_ops);

    __m256i insert_p = _mm256_set1_epi32(Cigar::Operator::INSERTION);
    for (size_t i = 1; i < length; i += 8) {
        __m256i H_orig = _mm256_loadu_si256((__m256i*)&update_scores[i]);

        // check match
        __m256i incoming_p = _mm256_loadu_si256((__m256i*)&incoming_scores[i - 1]);
        __m256i match_score = _mm256_add_epi32(
            incoming_p,
            _mm256_loadu_si256((__m256i*)&profile_scores[i])
        );

        __m256i H = _mm256_max_epi32(H_orig, match_score);

        // check insert
        __m256i update_score = _mm256_max_epi32(min_cell_score, _mm256_max_epi32(
            _mm256_add_epi32(rshiftpushback_epi32(incoming_p, incoming_scores[i + 7]),
                             gap_opening_penalty),
            _mm256_blendv_epi8(
                min_cell_score,
                _mm256_add_epi32(_mm256_loadu_si256((__m256i*)&incoming_gap_scores[i]),
                                 gap_extension_penalty),
                _mm256_cmpeq_epi32(_mm256_loadu_si256((__m256i*)&incoming_ops[i]),
                                   insert_p)
            )
        ));


        __m256i update_cmp = _mm256_cmpgt_epi32(update_score, H);
        H = _mm256_max_epi32(H, update_score);

        // update scores
        _mm256_storeu_si256((__m256i*)&update_scores[i], H);
        _mm256_storeu_si256((__m256i*)&update_gap_scores[i], update_score);

        __m256i both_cmp = _mm256_cmpgt_epi32(H, H_orig);
        _mm256_maskstore_epi32(&update_ops[i], both_cmp,
            _mm256_blendv_epi8(_mm256_loadu_si256((__m256i*)&profile_ops[i]),
                               insert_p,
                               update_cmp)
        );

        _mm256_maskstore_epi32(&updated_mask[i], both_cmp, both_cmp);

        __m128i *cmp = (__m128i*)&both_cmp;
        _mm256_maskstore_epi64(&update_prevs[i], _mm256_cvtepi32_epi64(cmp[0]), prev_node);
        _mm256_maskstore_epi64(&update_prevs[i + 4], _mm256_cvtepi32_epi64(cmp[1]), prev_node);
    }
}

#endif

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t>
inline void compute_updates(const DBGAlignerConfig &config,
                            score_t *update_scores,
                            score_t *update_gap_scores,
                            NodeType *update_prevs,
                            Cigar::Operator *update_ops,
                            const NodeType &prev_node,
                            const NodeType &node,
                            const score_t *incoming_scores,
                            const score_t *incoming_gap_scores,
                            const Cigar::Operator *incoming_ops,
                            const score_t *profile_scores,
                            const Cigar::Operator *profile_ops,
                            AlignedVector<int32_t> &updated_mask,
                            size_t length) {
    assert(length);
    score_t min_cell_score = config.min_cell_score;

    // handle first element (i.e., no match update possible)
    score_t update_score = std::max(min_cell_score, std::max(
        incoming_scores[0] + config.gap_opening_penalty,
        incoming_ops[0] == Cigar::Operator::INSERTION
            ? incoming_gap_scores[0] + config.gap_extension_penalty
            : min_cell_score
    ));

    if (update_score > update_scores[0]) {
        update_scores[0] = update_score;
        update_ops[0] = Cigar::Operator::INSERTION;
        update_prevs[0] = prev_node;
        updated_mask[0] = 0xFFFFFFFF;
    }

    update_gap_scores[0] = update_score;

    size_t i = 1;

#ifdef __AVX2__

    // ensure sizes are the same before casting for AVX2 instructions
    static_assert(sizeof(NodeType) == sizeof(long long int));
    static_assert(sizeof(Cigar::Operator) == sizeof(int32_t));

    if (prev_node != node) {
        // update 8 scores at a time
        compute_HE_avx2(length,
                        _mm256_set1_epi64x(prev_node),
                        _mm256_set1_epi32(config.gap_opening_penalty),
                        _mm256_set1_epi32(config.gap_extension_penalty),
                        update_scores,
                        update_gap_scores,
                        reinterpret_cast<long long int*>(update_prevs),
                        reinterpret_cast<int32_t*>(update_ops),
                        incoming_scores,
                        incoming_gap_scores,
                        reinterpret_cast<const int32_t*>(incoming_ops),
                        reinterpret_cast<const int32_t*>(profile_scores),
                        reinterpret_cast<const int32_t*>(profile_ops),
                        updated_mask.data(),
                        _mm256_set1_epi32(min_cell_score));
        i = length;
    }

#endif

    for (; i < length; ++i) {
        // store score updates
        score_t H = update_scores[i];
        Cigar::Operator H_op;
        NodeType H_prev = DeBruijnGraph::npos;

        // check match
        score_t match_score = incoming_scores[i - 1] + profile_scores[i];
        if (match_score > H) {
            H = match_score;
            H_op = profile_ops[i];
            H_prev = prev_node;
        }

        // check insert
        score_t update_score = std::max(min_cell_score, std::max(
            incoming_scores[i] + config.gap_opening_penalty,
            incoming_ops[i] == Cigar::Operator::INSERTION
                ? incoming_gap_scores[i] + config.gap_extension_penalty
                : min_cell_score
        ));

        if (update_score > H) {
            H = update_score;
            H_op = Cigar::Operator::INSERTION;
            H_prev = prev_node;
        }

        update_gap_scores[i] = update_score;

        // update scores
        if (H_prev != DeBruijnGraph::npos) {
            update_scores[i] = H;
            update_ops[i] = H_op;
            update_prevs[i] = H_prev;
            updated_mask[i] = 0xFFFFFFFF;
        }
    }

    update_del_scores(config,
                      update_scores,
                      update_prevs,
                      update_ops,
                      node,
                      updated_mask.data(),
                      length);
}


template <typename NodeType>
void DefaultColumnExtender<NodeType>
::operator()(std::function<void(DBGAlignment&&, NodeType)> callback,
             score_t min_path_score) {
    assert(graph_);
    assert(columns_to_update.empty());

    const auto &path = get_seed();

    // stop path early if it can't be better than the min_path_score
    if (path.get_score() + match_score_begin[1] < min_path_score)
        return;

    size_t query_clipping = path.get_clipping() + path.get_query().size() - 1;

    if (dp_table.size()
            && query_clipping >= dp_table.get_query_offset()
            && std::all_of(path.begin(), path.end(),
                           [&](auto i) { return dp_table.find(i) != dp_table.end(); })) {
        auto find = dp_table.find(path.back());
        if (find->second.scores.at(query_clipping - dp_table.get_query_offset())
                   >= path.get_score()) {
            return;
        }
    }

    reset();
    if (!add_seed(query_clipping))
        return;

    start_score = dp_table.find(start_node)->second.best_score();
    score_cutoff = std::max(start_score, min_path_score);
    begin = 0;
    end = size;
    xdrop_cutoff = start_score;

    check_and_push(ColumnRef(start_node, start_score, false));

    if (columns_to_update.empty())
        return;

    extend_main(callback, min_path_score);
}

template <typename NodeType>
std::deque<std::pair<NodeType, char>> DefaultColumnExtender<NodeType>
::fork_extension(NodeType node,
                 std::function<void(DBGAlignment&&, NodeType)>,
                 score_t) {
    overlapping_range_ = false;
    std::deque<std::pair<DeBruijnGraph::node_index, char>> out_columns;
    graph_->call_outgoing_kmers(node, [&](auto next_node, char c) {
        if (c != '$') {
            if (next_node == node) {
                out_columns.emplace_front(next_node, c);
            } else {
                out_columns.emplace_back(next_node, c);
            }
        }
    });

    return out_columns;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::extend_main(std::function<void(DBGAlignment&&, NodeType)> callback,
              score_t min_path_score) {
    assert(start_score == dp_table.best_score().second);

    while (columns_to_update.size()) {
        ColumnRef top = columns_to_update.top();
        NodeType node = std::get<0>(top);
        score_t best_score_update = std::get<1>(top);
        const auto &cur_col = dp_table.find(node)->second;

        // if this happens, then it means that the column was in the priority
        // queue multiple times, so we don't need to consider it again
        if (best_score_update != cur_col.last_priority_value()) {
            columns_to_update.pop();
            continue;
        }

        auto out_columns = fork_extension(node, callback, min_path_score);

        assert(std::all_of(out_columns.begin(), out_columns.end(), [&](const auto &pair) {
            return graph_->traverse(node, pair.second) == pair.first;
        }));

        columns_to_update.pop();

        update_columns(node, out_columns, min_path_score);
    }

    assert(start_score > config_.min_cell_score);

    // no good path found
    if (start_node == SequenceGraph::npos
            || start_score == get_seed().get_score()
            || score_cutoff > start_score)
        return;

    // check to make sure that start_node stores the best starting point
    assert(start_score == dp_table.best_score().second);

    if (dp_table.find(start_node)->second.best_op() != Cigar::Operator::MATCH)
        logger->trace("best alignment does not end with a MATCH");

    // get all alignments
    dp_table.extract_alignments(*graph_,
                                config_,
                                std::string_view(align_start, size - 1),
                                callback,
                                min_path_score,
                                get_seed(),
                                &start_node);
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::update_columns(NodeType incoming_node,
                 const std::deque<std::pair<NodeType, char>> &out_columns,
                 score_t min_path_score) {
    // set boundaries for vertical band
    auto *incoming = &dp_table.find(incoming_node).value();

    if (dp_table.size() == 1 && out_columns.size() && out_columns.front().first != incoming_node) {
        update_del_scores(config_,
                          incoming->scores.data(),
                          incoming->prev_nodes.data(),
                          incoming->ops.data(),
                          incoming_node,
                          nullptr,
                          size);
    }

    std::tie(begin, end) = get_band(size, incoming->best_pos, config_.bandwidth);

    for (const auto &[next_node, c] : out_columns) {
        auto emplace = emplace_node(next_node, incoming_node, c, size);

        // emplace_node may have invalidated incoming, so update the pointer
        if (emplace.second)
            incoming = &dp_table.find(incoming_node).value();

        auto &next_column = emplace.first.value();

        // store the mask indicating which cells were updated
        // this is padded to ensure that the vectorized code doesn't access
        // out of bounds
        AlignedVector<int32_t> updated_mask(end - begin + 8, false);

        assert(next_column.scores.size() == next_column.gap_scores.size());
        assert(incoming->scores.size() == incoming->gap_scores.size());
        compute_updates(
            config_,
            next_column.scores.data() + begin,
            next_column.gap_scores.data() + begin,
            next_column.prev_nodes.data() + begin,
            next_column.ops.data() + begin,
            incoming_node,
            next_node,
            incoming->scores.data() + begin,
            incoming->gap_scores.data() + begin,
            incoming->ops.data() + begin,
            profile_score[next_column.last_char].data() + query.size() - size + begin,
            profile_op[next_column.last_char].data() + query.size() - size + begin,
            updated_mask,
            end - begin
        );

        // Find the maximum changed value
        const score_t *best_update = nullptr;
        for (size_t i = begin, j = 0; i < end; ++i, ++j) {
            if (updated_mask[j] && (!best_update || next_column.scores[i] > *best_update))
                best_update = &next_column.scores[i];
        }

        // put column back in the priority queue if it's updated
        if (best_update) {
            assert(updated_mask[best_update - &next_column.scores[begin]]);

            next_column.last_priority_pos = best_update - next_column.scores.data();

            if (*best_update > next_column.best_score()) {
                next_column.best_pos = best_update - next_column.scores.data();

                assert(next_column.best_pos >= begin);
                assert(next_column.best_pos < end);
                assert(next_column.best_pos < size);
            }

            assert(*best_update == next_column.best_score()
                || next_column.best_pos < begin
                || next_column.best_pos >= end
                || !updated_mask[next_column.best_pos - begin]);

            // update global max score
            if (*best_update > start_score) {
                start_node = next_node;
                start_score = *best_update;
                xdrop_cutoff = std::max(start_score, xdrop_cutoff);
                assert(start_score == dp_table.best_score().second);
                score_cutoff = std::max(start_score, min_path_score);
            }

            check_and_push(ColumnRef(next_node, *best_update, !overlapping_range_));
        }
    }
}


template class DefaultColumnExtender<>;
