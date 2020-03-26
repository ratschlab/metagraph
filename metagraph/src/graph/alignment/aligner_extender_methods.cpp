#include "aligner_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "common/logger.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/vectors/aligned_vector.hpp"

using mg::common::logger;


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
    }

    return std::make_pair(find, false);
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>
::add_seed(size_t clipping) {
    return dp_table.add_seed(start_node, *(align_start - 1), path_->get_score(),
                             config_.min_cell_score, size, 0,
                             config_.gap_opening_penalty, config_.gap_extension_penalty,
                             clipping);
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    // this extender only works if at least one character has been matched
    assert(path.get_query_end() > path.get_query().data());
    assert(path.get_query_end() > query.data());
    assert(query.data() + query.size() >= path.get_query_end());

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
auto DefaultColumnExtender<NodeType>
::extendable(const ColumnRef &next_column) const -> typename ColumnQueue::Decision {
    assert(dp_table.find(next_column.first) != dp_table.end());

    const auto &column = dp_table.find(next_column.first)->second;

    // Ignore if there is no way it can be extended to an optimal alignment.
    // TODO: this cuts off too early (before the scores have converged)
    //       so path scores have to be recomputed after alignment
    if (xdrop_cutoff - column.best_score() > config_.xdrop
            || std::equal(match_score_begin + begin,
                          match_score_begin + end,
                          column.scores.data() + begin,
                          [&](auto a, auto b) { return a + b < score_cutoff; })) {
        return ColumnQueue::Decision::IGNORE;
    }

    if (columns_to_update.size() < config_.queue_size)
        return ColumnQueue::Decision::ADD;

    return columns_to_update.compare(columns_to_update.bottom(), next_column)
        ? ColumnQueue::Decision::REPLACE_BOTTOM
        : ColumnQueue::Decision::IGNORE;
}


template <typename NodeType>
std::vector<std::pair<NodeType, char>> DefaultColumnExtender<NodeType>
::fork_extension(NodeType node,
                 std::function<void(DBGAlignment&&, NodeType)>,
                 score_t) {
    std::vector<std::pair<DeBruijnGraph::node_index, char>> out_columns;
    graph_->call_outgoing_kmers(node, [&](auto next_node, char c) {
        out_columns.emplace_back(next_node, c);
    });

    return out_columns;
}


/*
 * Helpers for DefaultColumnExtender::operator()
 */

#ifdef __AVX2__

template <class ScoreRow, class OpRow>
inline bool compute_updates_avx2(size_t length,
                                 int64_t prev_node,
                                 int64_t node,
                                 int32_t gap_opening_penalty,
                                 int32_t gap_extension_penalty,
                                 int32_t *update_scores,
                                 long long int *update_prevs,
                                 int32_t *update_ops,
                                 const int32_t *incoming_scores,
                                 const int32_t *incoming_ops,
                                 const ScoreRow &row,
                                 const OpRow &op_row,
                                 const char *align_start,
                                 int32_t *updated_mask) {
    // Note: the char_scores and match_ops arrays are missing the first cell
    // (so the indices are shifted down one) to ensure byte alignment
    bool updated = false;

    for (size_t i = 1; i < length; i += 8) {
        // load characters of input string into 32-bit cells
        __m256i str_packed = expandepu8_epi32(*(long long int*)&align_start[i - 1]);

        // load scores from incoming column
        __m256i incoming_packed = _mm256_loadu_si256((__m256i*)&incoming_scores[i - 1]);

        // TODO: use encoded sequences so the gather can be replaced with a shuffle
        // compute match scores
        __m256i match_scores_packed = _mm256_add_epi32(
            _mm256_i32gather_epi32(row.data(), str_packed, 4), // gather match scores
            incoming_packed
        );

        // compute insertion scores
        __m256i ins_packed = _mm256_set1_epi32(Cigar::Operator::INSERTION);
        __m256i incoming_ops_packed = _mm256_loadu_si256((__m256i*)&incoming_ops[i]);
        __m256i incoming_is_insert = _mm256_cmpeq_epi32(incoming_ops_packed, ins_packed);
        __m256i incoming_is_delete = _mm256_cmpgt_epi32(incoming_ops_packed, ins_packed);

        __m256i incoming_shift = rshiftpushback_epi32(incoming_packed, incoming_scores[i + 7]);
        __m256i insert_scores_packed = _mm256_add_epi32(
            incoming_shift,
            _mm256_blendv_epi8(_mm256_set1_epi32(gap_opening_penalty),
                               _mm256_set1_epi32(gap_extension_penalty),
                               incoming_is_insert)
        );

        // disallow insertion after deletion
        static_assert(Cigar::Operator::DELETION > Cigar::Operator::INSERTION);
        insert_scores_packed = _mm256_blendv_epi8(
            insert_scores_packed,
            _mm256_set1_epi32(std::numeric_limits<int32_t>::min()),
            incoming_is_delete
        );

        // update scores
        __m256i max_scores = _mm256_max_epi32(match_scores_packed, insert_scores_packed);


        // mark which indices are valid
        __m256i idx = _mm256_cmpgt_epi32(
            _mm256_set1_epi32(length),
            _mm256_add_epi32(_mm256_set1_epi32(i),
                             _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7))
        );
        // mark which indices should be updated
        __m256i both_cmp = _mm256_and_si256(
            _mm256_cmpgt_epi32(max_scores, _mm256_loadu_si256((__m256i*)&update_scores[i])),
            idx
        );
        _mm256_maskstore_epi32(&update_scores[i], both_cmp, max_scores);

        // set mask indicating which scores were updated
        _mm256_maskstore_epi32(&updated_mask[i], both_cmp, both_cmp);
        updated |= static_cast<bool>(_mm256_movemask_epi8(both_cmp));

        // TODO: use encoded sequences so the gather can be replaced with a shuffle
        // update ops
        // pick match operation if match score >= gap score
        // pick delete operation if gap score > match score
        _mm256_maskstore_epi32(
            &update_ops[i],
            both_cmp,
            _mm256_blendv_epi8(_mm256_i32gather_epi32((int*)op_row.data(), str_packed, 4),
                               ins_packed,
                               _mm256_cmpgt_epi32(insert_scores_packed, match_scores_packed))
        );

        // update prev nodes
        // TODO: this can be done with one AVX512 instruction
        __m256i prev_packed = _mm256_set1_epi64x(prev_node);
        __m128i *cmp_array = (__m128i*)&both_cmp;
        _mm256_maskstore_epi64(&update_prevs[i],
                               _mm256_cvtepi32_epi64(cmp_array[0]),
                               prev_packed);
        _mm256_maskstore_epi64(&update_prevs[i + 4],
                               _mm256_cvtepi32_epi64(cmp_array[1]),
                               prev_packed);
    }

    for (size_t i = 1; i < length; ++i) {
        // disallow DELETION after INSERTION
        if (update_ops[i - 1] != Cigar::Operator::INSERTION) {
            score_t gap_score = update_scores[i - 1]
                + (update_ops[i - 1] == Cigar::Operator::DELETION
                    ? gap_extension_penalty
                    : gap_opening_penalty);

            if (gap_score > update_scores[i]) {
                while (i < length && gap_score > update_scores[i]) {
                    updated_mask[i] = 0xFFFFFFFF;
                    updated = true;
                    update_ops[i] = Cigar::Operator::DELETION;
                    update_prevs[i] = node;
                    update_scores[i] = gap_score;
                    ++i;
                    gap_score += gap_extension_penalty;
                }
                --i;
            }
        }
    }

    return updated;
}

#endif

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t,
          class ScoreRow,
          class OpRow>
inline bool compute_updates(const DBGAlignerConfig &config,
                            score_t *update_scores,
                            NodeType *update_prevs,
                            Cigar::Operator *update_ops,
                            const NodeType &prev_node,
                            const NodeType &node,
                            const score_t *incoming_scores,
                            const Cigar::Operator *incoming_ops,
                            const ScoreRow &row,
                            const OpRow &op_row,
                            const char *align_start,
                            AlignedVector<int32_t> &updated_mask,
                            size_t length) {
    // Note: the char_scores and match_ops arrays are missing the first cell
    // (so the indices are shifted down one) to ensure byte alignment
    bool updated = false;

    // handle first element (i.e., no match update possible)
    score_t gap_score = incoming_scores[0] + (incoming_ops[0] == Cigar::Operator::INSERTION
        ? config.gap_extension_penalty : config.gap_opening_penalty);

    if (gap_score > update_scores[0]) {
        updated_mask[0] = 0xFFFFFFFF;
        updated = true;
        update_ops[0] = Cigar::Operator::INSERTION;
        update_prevs[0] = prev_node;
        update_scores[0] = gap_score;
    }

#ifdef __AVX2__

    // ensure sizes are the same before casting for AVX2 instructions
    static_assert(sizeof(NodeType) == sizeof(long long int));
    static_assert(sizeof(Cigar::Operator) == sizeof(int32_t));

    // update 8 scores at a time
    updated |= compute_updates_avx2(length,
                                    prev_node,
                                    node,
                                    config.gap_opening_penalty,
                                    config.gap_extension_penalty,
                                    update_scores,
                                    reinterpret_cast<long long int*>(update_prevs),
                                    reinterpret_cast<int32_t*>(update_ops),
                                    incoming_scores,
                                    reinterpret_cast<const int32_t*>(incoming_ops),
                                    row,
                                    op_row,
                                    align_start,
                                    updated_mask.data());

#else

    for (size_t i = 1; i < length; ++i) {
        if (incoming_scores[i - 1] + row[align_start[i - 1]] > update_scores[i]) {
            updated_mask[i] = 0xFFFFFFFF;
            updated = true;
            update_ops[i] = op_row[align_start[i - 1]];
            update_prevs[i] = prev_node;
            update_scores[i] = incoming_scores[i - 1] + row[align_start[i - 1]];
        }

        // disallow INSERTION after DELETION?
        if (incoming_ops[i] != Cigar::Operator::DELETION) {
            gap_score = incoming_scores[i] + (incoming_ops[i] == Cigar::Operator::INSERTION
                ? config.gap_extension_penalty
                : config.gap_opening_penalty);

            if (gap_score > update_scores[i]) {
                updated_mask[i] = 0xFFFFFFFF;
                updated = true;
                update_ops[i] = Cigar::Operator::INSERTION;
                update_prevs[i] = prev_node;
                update_scores[i] = gap_score;
            }
        }

        // disallow DELETION after INSERTION
        if (update_ops[i - 1] != Cigar::Operator::INSERTION) {
            gap_score = update_scores[i - 1]
                + (update_ops[i - 1] == Cigar::Operator::DELETION
                    ? config_.gap_extension_penalty
                    : config_.gap_opening_penalty);

            if (gap_score > update_scores[i]) {
                updated_mask[i] = 0xFFFFFFFF;
                updated = true;
                update_ops[i] = Cigar::Operator::DELETION;
                update_prevs[i] = node;
                update_scores[i] = gap_score;
            }
        }
    }

#endif

    return updated;
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
    end = partial_sums_.size();
    xdrop_cutoff = start_score;

    columns_to_update.emplace(start_node, start_score);

    if (columns_to_update.empty())
        return;

    extend_main(callback, min_path_score);
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::extend_main(std::function<void(DBGAlignment&&, NodeType)> callback,
              score_t min_path_score) {
    assert(start_score == dp_table.best_score().second);

    while (columns_to_update.size()) {
        const auto next_pair = columns_to_update.top();
        const auto &cur_col = dp_table.find(next_pair.first)->second;

        // if this happens, then it means that the column was in the priority
        // queue multiple times, so we don't need to consider it again
        if (next_pair.second != cur_col.last_priority_value()) {
            columns_to_update.pop();
            continue;
        }

        auto out_columns = fork_extension(next_pair.first, callback, min_path_score);
        columns_to_update.pop();

        update_columns(next_pair.first, out_columns, min_path_score);
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
                                get_seed().get_orientation(),
                                min_path_score,
                                &start_node);
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::update_columns(NodeType incoming_node,
                 const std::vector<std::pair<NodeType, char>> &out_columns,
                 score_t min_path_score) {
    // set boundaries for vertical band
    auto *incoming = &dp_table.find(incoming_node).value();
    size_t best_pos = incoming->best_pos;
    assert(best_pos < size);
    begin = best_pos >= config_.bandwidth ? best_pos - config_.bandwidth : 0;
    end = config_.bandwidth <= size - best_pos ? best_pos + config_.bandwidth : size;

    assert(begin <= best_pos);
    assert(end > best_pos);
    assert(begin < end);

    for (const auto &[next_node, c] : out_columns) {
        auto emplace = emplace_node(next_node, incoming_node, c, size);

        // emplace_node may have invalidated incoming, so update the pointer
        if (emplace.second)
            incoming = &dp_table.find(incoming_node).value();

        auto &next_column = emplace.first.value();

        // store the mask indicating which cells were updated
        AlignedVector<int32_t> updated_mask(end - begin + 8, false);

        bool updated = compute_updates(config_,
                                       next_column.scores.data() + begin,
                                       next_column.prev_nodes.data() + begin,
                                       next_column.ops.data() + begin,
                                       incoming_node,
                                       next_node,
                                       incoming->scores.data() + begin,
                                       incoming->ops.data() + begin,
                                       config_.get_row(next_column.last_char),
                                       Cigar::get_op_row(next_column.last_char),
                                       align_start + begin,
                                       updated_mask,
                                       end - begin);

        // Find the maximum changed value
        const score_t *best_update = nullptr;
        for (size_t i = begin, j = 0; i < end; ++i, ++j) {
            if (updated_mask[j] && (!best_update || next_column.scores[i] > *best_update))
                best_update = &next_column.scores[i];
        }

        assert(updated == static_cast<bool>(best_update));

        // put column back in the priority queue if it's updated
        if (updated) {
            assert(updated_mask[best_update - &next_column.scores[begin]]);

            next_column.last_priority_pos = best_update - next_column.scores.data();

            if (*best_update > next_column.best_score()) {
                next_column.best_pos = best_update - next_column.scores.data();

                assert(next_column.best_pos >= begin);
                assert(next_column.best_pos < end);
                assert(next_column.best_pos < next_column.size());
                assert(next_column.best_score()
                    == *std::max_element(next_column.scores.begin(), next_column.scores.end()));
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

            columns_to_update.emplace(next_node, *best_update);
        }
    }
}


template class DefaultColumnExtender<>;
