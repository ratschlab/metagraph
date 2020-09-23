#include "aligner_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "common/logger.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;


template <typename NodeType>
void DefaultColumnExtender<NodeType>
::initialize_query(const std::string_view query) {
    this->query = query;
    max_num_nodes = config_.max_nodes_per_seq_char < std::numeric_limits<double>::max()
                        ? std::ceil(config_.max_nodes_per_seq_char * query.size())
                        : std::numeric_limits<size_t>::max();

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

template <class Vector>
std::pair<size_t, size_t> get_column_boundaries(const Vector &scores,
                                                size_t size,
                                                size_t best_pos,
                                                DBGAlignerConfig::score_t xdrop_cutoff,
                                                size_t bandwidth) {
    if (xdrop_cutoff > scores[best_pos])
        return std::make_pair(best_pos, best_pos);

    assert(best_pos < size);
    size_t begin = best_pos >= bandwidth ? best_pos - bandwidth : 0;
    size_t end = bandwidth <= size - best_pos ? best_pos + bandwidth : size;

    while (begin < end && scores[begin] < xdrop_cutoff) {
        ++begin;
    }

    if (begin == size)
        return std::make_pair(begin, begin);

    // ensure that the next position is included in the range [begin, end)
    size_t cur_end = best_pos + 2;
    while (cur_end < end && scores[cur_end] >= xdrop_cutoff) {
        ++cur_end;
    }
    end = cur_end;

    assert(end > best_pos);

    if (begin >= end)
        return std::make_pair(begin, begin);

    // align begin + 1 to 32-byte boundary
    if (begin > 7)
        begin = (begin & 0xFFFFFFFFFFFFFFF8) - 1;

    // round up to nearest multiple of 8
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
        auto [node_begin, node_end] = get_column_boundaries(find->second.scores,
                                                            size,
                                                            find->second.best_pos,
                                                            xdrop_cutoff,
                                                            config_.bandwidth);

        if (node_begin != node_end)
            overlapping_range_ |= (begin < node_end && end > node_begin);
    }

    return std::make_pair(find, false);
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>
::add_seed(size_t clipping) {
    assert(path_->get_cigar().back().first == Cigar::MATCH
        || path_->get_cigar().back().first == Cigar::MISMATCH);

    return dp_table.add_seed(*path_, config_, size, 0, clipping);
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
    if (xdrop_cutoff > column.best_score()
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

    assert(!columns_to_update.empty());

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
                              size_t length,
                              score_t xdrop_cutoff) {
    for (size_t i = 1; i < length; ++i) {
        score_t del_score = std::max(config.min_cell_score,
            update_scores[i - 1] + (update_ops[i - 1] == Cigar::DELETION
                ? config.gap_extension_penalty
                : config.gap_opening_penalty
        ));

        if (del_score >= xdrop_cutoff && del_score > update_scores[i]) {
            while (i < length && del_score > update_scores[i]) {
                update_scores[i] = del_score;
                update_ops[i] = Cigar::DELETION;
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
                            int32_t *update_gap_count,
                            long long int *update_gap_prevs,
                            const int32_t *incoming_scores,
                            const int32_t *incoming_gap_scores,
                            const int32_t *profile_scores,
                            const int32_t *profile_ops,
                            const int32_t *incoming_gap_count,
                            int32_t *updated_mask,
                            __m256i xdrop_cutoff) {
    assert(update_scores != incoming_scores);
    assert(update_gap_scores != incoming_gap_scores);

    __m256i insert_p = _mm256_set1_epi32(Cigar::INSERTION);
    for (size_t i = 1; i < length; i += 8) {
        // load previous values for cells to update
        __m256i H_orig = _mm256_loadu_si256((__m256i*)&update_scores[i]);
        __m256i gap_orig = _mm256_loadu_si256((__m256i*)&update_gap_scores[i]);

        // compute match score
        __m256i incoming_p = _mm256_loadu_si256((__m256i*)&incoming_scores[i - 1]);
        __m256i match_score = _mm256_add_epi32(
            incoming_p,
            _mm256_loadu_si256((__m256i*)&profile_scores[i])
        );

        // compute score for cell update
        __m256i H = _mm256_max_epi32(H_orig, match_score);

        // compute insert score
        __m256i update_score_open = _mm256_add_epi32(
            rshiftpushback_epi32(incoming_p, incoming_scores[i + 7]),
            gap_opening_penalty
        );
        __m256i update_score_extend = _mm256_add_epi32(
            _mm256_loadu_si256((__m256i*)&incoming_gap_scores[i]),
            gap_extension_penalty
        );
        __m256i update_score = _mm256_max_epi32(update_score_open, update_score_extend);

        // compute updated gap size count
        __m256i ones = _mm256_set1_epi32(1);
        __m256i incoming_count = _mm256_add_epi32(
            _mm256_loadu_si256((__m256i*)&incoming_gap_count[i]),
            ones
        );

        // compute score for cell update. check if inserting a gap improves the update
        __m256i update_cmp = _mm256_cmpgt_epi32(update_score, H);
        H = _mm256_max_epi32(H, update_score);

        // determine which indices satisfy the x-drop criteria
        __m256i xdrop_cmp = _mm256_cmpgt_epi32(H, xdrop_cutoff);

        // mask indices which are out of bounds
        xdrop_cmp = _mm256_and_si256(_mm256_cmpgt_epi32(_mm256_set1_epi32(length),
                                     _mm256_add_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7),
                                                      _mm256_set1_epi32(i))), xdrop_cmp);

        // revert values not satisfying the x-drop criteria
        H = _mm256_blendv_epi8(H_orig, H, xdrop_cmp);
        update_score = _mm256_blendv_epi8(gap_orig, update_score, xdrop_cmp);

        // check which updates should be stores
        __m256i both_cmp = _mm256_cmpgt_epi32(H, H_orig);
        __m256i gap_cmp = _mm256_cmpgt_epi32(update_score, gap_orig);

        // update scores in cells of DP table

        // update gap count
        _mm256_maskstore_epi32(&update_gap_count[i], gap_cmp, _mm256_blendv_epi8(
            ones,
            incoming_count,
            _mm256_cmpeq_epi32(update_score, update_score_extend)
        ));

        // update scores
        _mm256_maskstore_epi32(&update_scores[i], both_cmp, H);
        _mm256_maskstore_epi32(&update_gap_scores[i], gap_cmp, update_score);

        // update traceback operation
        _mm256_maskstore_epi32(&update_ops[i], both_cmp,
            _mm256_blendv_epi8(_mm256_loadu_si256((__m256i*)&profile_ops[i]),
                               insert_p,
                               update_cmp)
        );

        // update traceback node
        const __m128i *cmp = (__m128i*)&both_cmp;
        _mm256_maskstore_epi64(&update_prevs[i], _mm256_cvtepi32_epi64(cmp[0]), prev_node);
        _mm256_maskstore_epi64(&update_prevs[i + 4], _mm256_cvtepi32_epi64(cmp[1]), prev_node);
        cmp = (__m128i*)&gap_cmp;
        _mm256_maskstore_epi64(&update_gap_prevs[i], _mm256_cvtepi32_epi64(cmp[0]), prev_node);
        _mm256_maskstore_epi64(&update_gap_prevs[i + 4], _mm256_cvtepi32_epi64(cmp[1]), prev_node);

        // update mask
        both_cmp = _mm256_or_si256(both_cmp, gap_cmp);
        _mm256_maskstore_epi32(&updated_mask[i], both_cmp, both_cmp);
    }
}

#endif

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t>
inline void compute_updates(const DBGAlignerConfig &config,
                            score_t *update_scores,
                            score_t *update_gap_scores,
                            NodeType *update_prevs,
                            int32_t *update_gap_count,
                            NodeType *update_gap_prevs,
                            Cigar::Operator *update_ops,
                            const NodeType &prev_node,
                            const NodeType &node,
                            const score_t *incoming_scores,
                            const score_t *incoming_gap_scores,
                            const int32_t *incoming_gap_count,
                            const score_t *profile_scores,
                            const Cigar::Operator *profile_ops,
                            AlignedVector<int32_t> &updated_mask,
                            size_t length,
                            score_t xdrop_cutoff) {
    assert(length);

    // handle first element (i.e., no match update possible)
    score_t update_score = std::max(
        incoming_scores[0] + config.gap_opening_penalty,
        incoming_gap_scores[0] + config.gap_extension_penalty
    );

    if (update_score >= xdrop_cutoff) {
        if (update_score > update_gap_scores[0]) {
            update_gap_count[0] = update_score == incoming_scores[0] + config.gap_opening_penalty
                ? 1 : incoming_gap_count[0] + 1;
            update_gap_scores[0] = update_score;
            update_gap_prevs[0] = prev_node;
            updated_mask[0] = 0xFFFFFFFF;
        }

        if (update_score > update_scores[0]) {
            update_scores[0] = update_score;
            update_ops[0] = Cigar::INSERTION;
            update_prevs[0] = prev_node;
            updated_mask[0] = 0xFFFFFFFF;
        }
    }

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
                        update_gap_count,
                        reinterpret_cast<long long int*>(update_gap_prevs),
                        incoming_scores,
                        incoming_gap_scores,
                        reinterpret_cast<const int32_t*>(profile_scores),
                        reinterpret_cast<const int32_t*>(profile_ops),
                        incoming_gap_count,
                        updated_mask.data(),
                        _mm256_set1_epi32(xdrop_cutoff - 1));
        i = length;
    }

#endif

    for (; i < length; ++i) {
        score_t H_orig = update_scores[i];

        // check match score
        score_t match_score = incoming_scores[i - 1] + score_t(profile_scores[i]);
        score_t H = std::max(H_orig, match_score);

        // check insert score
        score_t update_score_open = incoming_scores[i] + score_t(config.gap_opening_penalty);
        score_t update_score_extend = incoming_gap_scores[i] + score_t(config.gap_extension_penalty);
        score_t update_score = std::max(update_score_open, update_score_extend);

        // pick best one
        bool update_cmp = update_score > H;
        H = std::max(H, update_score);

        // update if scores improve
        if (H >= xdrop_cutoff) {
            if (update_score > update_gap_scores[i]) {
                update_gap_scores[i] = update_score;
                update_gap_count[i] = update_score == update_score_open
                    ? 1 : incoming_gap_count[i] + 1;
                update_gap_prevs[i] = prev_node;
                updated_mask[i] = 0xFFFFFFFF;
            }

            if (H > H_orig) {
                update_scores[i] = H;
                update_ops[i] = update_cmp ? Cigar::INSERTION : profile_ops[i];
                update_prevs[i] = prev_node;
                updated_mask[i] = 0xFFFFFFFF;
            }
        }
    }

    update_del_scores(config,
                      update_scores,
                      update_prevs,
                      update_ops,
                      node,
                      updated_mask.data(),
                      length,
                      xdrop_cutoff);
}


template <typename NodeType>
void DefaultColumnExtender<NodeType>
::operator()(std::function<void(DBGAlignment&&, NodeType)> callback,
             score_t min_path_score) {
    assert(graph_);
    assert(columns_to_update.empty());

    const auto &path = get_seed();

    if (!graph_->outdegree(path.back())) {
        callback(DBGAlignment(), NodeType());
        return;
    }

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
    xdrop_cutoff = score_cutoff - config_.xdrop;

    if (xdrop_cutoff > start_score)
        return;

    begin = 0;
    end = size;

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
        if (c != '$' && (dp_table.size() < max_num_nodes || dp_table.count(next_node))) {
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
        columns_to_update.pop();

        NodeType node = std::get<0>(top);
        score_t best_score_update = std::get<1>(top);
        const auto &cur_col = dp_table.find(node)->second;

        // if this happens, then it means that the column was in the priority
        // queue multiple times, so we don't need to consider it again
        if (best_score_update != cur_col.last_priority_value())
            continue;

        auto out_columns = fork_extension(node, callback, min_path_score);

        assert(std::all_of(out_columns.begin(), out_columns.end(), [&](const auto &pair) {
            return graph_->traverse(node, pair.second) == pair.first;
        }));

        update_columns(node, out_columns, min_path_score);
    }

    logger->trace("Extension completed:\tquery size:\t{}\tseed size:\t{}\texplored nodes:\t{}",
                  query.size(), path_->size(), dp_table.size());

    assert(start_score > config_.min_cell_score);

    // no good path found
    if (start_node == SequenceGraph::npos
            || start_score == get_seed().get_score()
            || score_cutoff > start_score) {
        reset();
        callback(Alignment<NodeType>(), NodeType());
        return;
    }

    // check to make sure that start_node stores the best starting point
    assert(start_score == dp_table.best_score().second);

    if (dp_table.find(start_node)->second.best_op() != Cigar::MATCH)
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
    if (out_columns.empty())
        return;

    // set boundaries for vertical band
    auto *incoming = &dp_table.find(incoming_node).value();

    if (dp_table.size() == 1 && out_columns.size() && out_columns.front().first != incoming_node) {
        update_del_scores(config_,
                          incoming->scores.data(),
                          incoming->prev_nodes.data(),
                          incoming->ops.data(),
                          incoming_node,
                          nullptr,
                          size,
                          xdrop_cutoff);
    }

    std::tie(begin, end) = get_column_boundaries(incoming->scores,
                                                 size,
                                                 incoming->best_pos,
                                                 xdrop_cutoff,
                                                 config_.bandwidth);

    if (begin >= end)
        return;

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
            next_column.gap_count.data() + begin,
            next_column.gap_prev_nodes.data() + begin,
            next_column.ops.data() + begin,
            incoming_node,
            next_node,
            incoming->scores.data() + begin,
            incoming->gap_scores.data() + begin,
            incoming->gap_count.data() + begin,
            profile_score[next_column.last_char].data() + query.size() - size + begin,
            profile_op[next_column.last_char].data() + query.size() - size + begin,
            updated_mask,
            end - begin,
            xdrop_cutoff
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
                xdrop_cutoff = std::max(start_score - config_.xdrop, xdrop_cutoff);
                assert(start_score == dp_table.best_score().second);
                score_cutoff = std::max(start_score, min_path_score);
            }

            check_and_push(ColumnRef(next_node, *best_update, !overlapping_range_));
        }
    }
}


template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
