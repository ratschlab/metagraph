#include "aligner_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include <Eigen/StdVector>

#include "common/algorithms.hpp"
#include "common/bounded_priority_queue.hpp"
#include "common/logger.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"

using mg::common::logger;


template <typename T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;

template <typename NodeType, typename score_t = typename DPTable<NodeType>::score_t>
using ColumnRef = std::pair<NodeType, score_t>;


/*
 * Helpers for DefaultColumnExtender::operator()
 */

// since insertion invalidates references, return a vector of pointers
template <typename NodeType,
          typename Column = typename DPTable<NodeType>::Column,
          typename score_t = typename DPTable<NodeType>::score_t>
inline std::vector<ColumnRef<NodeType>>
get_outgoing_columns(const DeBruijnGraph &graph,
                     DPTable<NodeType> &dp_table,
                     NodeType cur_node,
                     size_t size,
                     size_t best_pos,
                     score_t min_cell_score) {
    std::vector<ColumnRef<NodeType>> out_columns;

    graph.call_outgoing_kmers(
        cur_node,
        [&](auto next_node, char c) {
            auto find = dp_table.find(next_node);
            if (find == dp_table.end()) {
                find = dp_table.emplace(
                    next_node,
                    Column(size,
                           min_cell_score,
                           c,
                           best_pos + 1 != size ? best_pos + 1 : best_pos)
                ).first;
            }

            assert(find != dp_table.end());
            out_columns.emplace_back(next_node, find->second.best_score());
        }
    );

    return out_columns;
}

#ifdef __AVX2__

inline void compute_match_insert_updates_avx2(size_t &i,
                                              size_t length,
                                              int64_t prev_node,
                                              int32_t gap_opening_penalty,
                                              int32_t gap_extension_penalty,
                                              int32_t *&update_scores,
                                              long long int *&update_prevs,
                                              int32_t *&update_ops,
                                              const int32_t *&incoming_scores,
                                              const int32_t *&incoming_ops,
                                              const int8_t *&char_scores,
                                              const int32_t *&match_ops) {
    const __m256i prev_packed = _mm256_set1_epi64x(prev_node);
    const __m256i ins_packed = _mm256_set1_epi32(Cigar::Operator::INSERTION);
    const __m256i gap_open_packed = _mm256_set1_epi32(gap_opening_penalty);
    const __m256i gap_extend_packed = _mm256_set1_epi32(gap_extension_penalty);

    for (; i + 8 <= length; i += 8) {
        __m256i incoming_packed = _mm256_loadu_si256((__m256i*)(incoming_scores - 1));

        // compute match and delete scores
        __m256i match_scores_packed = _mm256_add_epi32(
            expandepi8_epi64(*(uint64_t*)char_scores),
            incoming_packed
        );

        __m256i insert_scores_packed = _mm256_add_epi32(
            rshiftpushback_epi32(incoming_packed, incoming_scores[7]),
            _mm256_blendv_epi8(gap_open_packed,
                               gap_extend_packed,
                               _mm256_cmpeq_epi32(
                                   _mm256_loadu_si256((__m256i*)incoming_ops),
                                   ins_packed
                               ))
        );

        // update scores
        __m256i max_scores = _mm256_max_epi32(match_scores_packed, insert_scores_packed);
        __m256i both_cmp = _mm256_cmpgt_epi32(
            max_scores,
            _mm256_loadu_si256((__m256i*)update_scores)
        );
        _mm256_maskstore_epi32(update_scores, both_cmp, max_scores);

        // update ops
        // pick match operation if match score >= gap score
        // pick delete operation if gap score > match score
        _mm256_maskstore_epi32(
            update_ops,
            both_cmp,
            _mm256_blendv_epi8(_mm256_load_si256((__m256i*)match_ops),
                               ins_packed,
                               _mm256_cmpgt_epi32(insert_scores_packed,
                                                  match_scores_packed))
        );

        // update prev nodes
        // TODO: this can be done with one AVX512 instruction
        __m128i *cmp_array = (__m128i*)&both_cmp;
        _mm256_maskstore_epi64(update_prevs,
                               _mm256_cvtepi32_epi64(cmp_array[0]),
                               prev_packed);
        _mm256_maskstore_epi64(update_prevs + 4,
                               _mm256_cvtepi32_epi64(cmp_array[1]),
                               prev_packed);

        update_scores += 8;
        update_prevs += 8;
        update_ops += 8;
        incoming_scores += 8;
        incoming_ops += 8;
        char_scores += 8;
        match_ops += 8;
    }
}

#endif

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t>
inline void compute_match_insert_updates(const DBGAlignerConfig &config,
                                         score_t *update_scores,
                                         NodeType *update_prevs,
                                         Cigar::Operator *update_ops,
                                         const NodeType &prev_node,
                                         const score_t *incoming_scores,
                                         const Cigar::Operator *incoming_ops,
                                         const AlignedVector<int8_t> &char_scores,
                                         const AlignedVector<Cigar::Operator> &match_ops,
                                         size_t length) {
    static_assert(sizeof(NodeType) == sizeof(long long int));
    static_assert(sizeof(Cigar::Operator) == sizeof(int32_t));
    auto update_prevs_cast = reinterpret_cast<long long int*>(update_prevs);
    auto update_ops_cast = reinterpret_cast<int32_t*>(update_ops);
    auto incoming_ops_cast = reinterpret_cast<const int32_t*>(incoming_ops);
    auto char_scores_data = char_scores.data();
    auto match_ops_cast = reinterpret_cast<const int32_t*>(match_ops.data());

    // handle first element (i.e., no match update possible)
    score_t gap_score = *incoming_scores + (*incoming_ops_cast == Cigar::Operator::INSERTION
        ? config.gap_extension_penalty : config.gap_opening_penalty);

    if (gap_score > *update_scores) {
        *update_ops_cast = Cigar::Operator::INSERTION;
        *update_prevs_cast = prev_node;
        *update_scores = gap_score;
    }

    ++incoming_scores;
    ++incoming_ops_cast;
    ++update_scores;
    ++update_prevs_cast;
    ++update_ops_cast;

    size_t i = 1;

#ifdef __AVX2__
    static_assert(std::is_same<score_t, int32_t>::value);

    // update 8 scores at a time
    compute_match_insert_updates_avx2(
        i,
        length,
        prev_node,
        config.gap_opening_penalty,
        config.gap_extension_penalty,
        update_scores,
        update_prevs_cast,
        update_ops_cast,
        incoming_scores,
        incoming_ops_cast,
        char_scores_data,
        match_ops_cast
    );
#endif

    // handle the residual
    for (; i < length; ++i) {
        if (*(incoming_scores - 1) + *char_scores_data > *update_scores) {
            *update_ops_cast = *match_ops_cast;
            *update_prevs_cast = prev_node;
            *update_scores = *(incoming_scores - 1) + *char_scores_data;
        }

        // TODO: enable check for deletion after insertion?
        gap_score = *incoming_scores + (*incoming_ops_cast == Cigar::Operator::INSERTION
            ? config.gap_extension_penalty : config.gap_opening_penalty);

        if (gap_score > *update_scores) {
            *update_ops_cast = Cigar::Operator::INSERTION;
            *update_prevs_cast = prev_node;
            *update_scores = gap_score;
        }

        ++incoming_scores;
        ++incoming_ops_cast;
        ++char_scores_data;
        ++match_ops_cast;
        ++update_scores;
        ++update_prevs_cast;
        ++update_ops_cast;
    }
}

#ifdef __AVX2__

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t,
          typename Column = typename DPTable<NodeType>::Column>
inline size_t update_column_avx2(bool &updated,
                                 AlignedVector<score_t> &update_scores,
                                 AlignedVector<Cigar::Operator> &update_ops,
                                 AlignedVector<NodeType> &update_prevs,
                                 Column &next_column,
                                 size_t overall_begin) {
    static_assert(sizeof(Cigar::Operator) == sizeof(int32_t));
    static_assert(sizeof(NodeType) == sizeof(long long int));
    static_assert(std::is_same<score_t, int32_t>::value);
    auto next_column_scores = reinterpret_cast<int32_t*>(
        &next_column.scores[overall_begin]
    );

    auto next_column_ops = reinterpret_cast<int32_t*>(
        &next_column.ops[overall_begin]
    );

    auto next_column_prev_nodes = reinterpret_cast<long long int*>(
        &next_column.prev_nodes[overall_begin]
    );

    size_t i = 0;
    for (; i + 8 <= update_scores.size(); i += 8) {
        // load update scores
        __m256i updates = _mm256_load_si256((__m256i*)&update_scores[i]);

        // compare updates to column
        __m256i cmp = _mm256_cmpgt_epi32(
            updates,
            _mm256_loadu_si256((__m256i*)next_column_scores)
        );

        updated |= bool(_mm256_movemask_epi8(cmp));

        // store updates in column
        _mm256_maskstore_epi32(next_column_scores, cmp, updates);
        _mm256_maskstore_epi32(next_column_ops, cmp,
                               _mm256_load_si256((__m256i*)&update_ops[i]));

        // TODO: this can be done with one AVX512 instruction
        __m128i *cmp_array = (__m128i*)&cmp;
        _mm256_maskstore_epi64(next_column_prev_nodes,
                               _mm256_cvtepi32_epi64(cmp_array[0]),
                               _mm256_load_si256((__m256i*)&update_prevs[i]));
        next_column_prev_nodes += 4;

        _mm256_maskstore_epi64(next_column_prev_nodes,
                               _mm256_cvtepi32_epi64(cmp_array[1]),
                               _mm256_load_si256((__m256i*)&update_prevs[i + 4]));
        next_column_prev_nodes += 4;

        next_column_scores += 8;
        next_column_ops += 8;
    }

    return i;
}

#endif


/*
 * DefaultColumnExtender::operator()
 */

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::operator()(const DBGAlignment &path,
             std::string_view query,
             std::function<void(DBGAlignment&&, NodeType)> callback,
             bool orientation,
             score_t min_path_score) {
    assert(graph_);
    size_t start_index = path.get_query_end() - 1 - query.data();
    const score_t *match_score_begin = partial_sums_.data() + start_index;

    // this extender only works if at least one character has been matched
    assert(path.get_query_end() > path.get_query().data());
    assert(query.data() + query.size() >= path.get_query_end());

    const auto *align_start = path.get_query_end();
    size_t size = query.data() + query.size() - align_start + 1;

    assert(config_.match_score(std::string_view(align_start - 1, size))
        == *match_score_begin);
    assert(config_.get_row(query.back())[query.back()] == match_score_begin[size - 1]);

    // stop path early if it can't be better than the min_path_score
    if (path.get_score() + match_score_begin[1] < min_path_score)
        return;

    // if the nodes from this seed were in the previous alignment and had a
    // better score, don't redo the extension
    assert(path.get_query_end() > query.data());
    size_t query_offset = path.get_query_end() - query.data() - 1;
    assert(query_offset + size == query.size());

    if (dp_table.size()
            && query_offset >= dp_table.get_query_offset()
            && std::all_of(path.begin(), path.end(),
                           [&](auto i) { return dp_table.find(i) != dp_table.end(); })) {
        auto find = dp_table.find(path.back());
        if (find != dp_table.end()
                && find->second.scores.at(query_offset - dp_table.get_query_offset())
                       >= path.get_score()) {
            return;
        }
    }

    reset();
    if (!dp_table.add_seed(path.back(),
                           *(align_start - 1),
                           path.get_score(),
                           config_.min_cell_score,
                           size,
                           0,
                           config_.gap_opening_penalty,
                           config_.gap_extension_penalty,
                           query_offset))
        return;

    // for storage of intermediate values
    AlignedVector<int8_t> char_scores;
    AlignedVector<Cigar::Operator> match_ops;

    AlignedVector<score_t> update_scores;
    AlignedVector<Cigar::Operator> update_ops;
    AlignedVector<node_index> update_prevs;

    // keep track of which columns to use next
    BoundedPriorityQueue<ColumnRef<NodeType>,
                         std::vector<ColumnRef<NodeType>>,
                         utils::LessSecond> columns_to_update(
        config_.queue_size
    );

    NodeType start_node = path.back();
    score_t start_score = dp_table.find(start_node)->second.best_score();
    score_t score_cutoff = std::max(start_score, min_path_score);

    columns_to_update.emplace(start_node, start_score);
    while (columns_to_update.size()) {
        const auto next_pair = columns_to_update.pop_top();
        const auto &cur_col = dp_table.find(next_pair.first)->second;

        // if this happens, then it means that the column was in the priority
        // queue multiple times, so we don't need to consider it again
        if (next_pair.second != cur_col.best_score())
            continue;

        // get next columns
        auto out_columns = get_outgoing_columns(*graph_,
                                                dp_table,
                                                next_pair.first,
                                                size,
                                                cur_col.best_pos,
                                                config_.min_cell_score);

        // update columns
        for (auto &[next_node, next_score] : out_columns) {
            auto *iter = &*dp_table.find(next_node);
            auto &next_column = const_cast<typename DPTable<NodeType>::Column&>(
                iter->second
            );

            // the get_outgoing_columns call may have invalidated cur_col, so
            // get a pointer to the column again
            const auto &incoming = dp_table.find(next_pair.first)->second;

            // find incoming nodes to check for alignment extension
            // set boundaries for vertical band
            size_t best_pos = incoming.best_pos;
            size_t begin = best_pos >= config_.bandwidth ? best_pos - config_.bandwidth : 0;
            size_t end = config_.bandwidth <= size - best_pos ? best_pos + config_.bandwidth : size;
            assert(end > begin);

            size_t next_best = next_column.best_pos;
            size_t overall_begin = std::min(
                begin,
                next_best >= config_.bandwidth ? next_best - config_.bandwidth : 0
            );

            size_t overall_end = std::max(
                end,
                config_.bandwidth <= size - next_best ? next_best + config_.bandwidth : size
            );
            assert(overall_begin <= best_pos);
            assert(overall_end > best_pos);


            update_scores.assign(overall_end - overall_begin, config_.min_cell_score);
            update_prevs.assign(overall_end - overall_begin, SequenceGraph::npos);
            update_ops.assign(overall_end - overall_begin, Cigar::Operator::CLIPPED);

            const auto &row = config_.get_row(next_column.last_char);
            const auto &op_row = Cigar::get_op_row(next_column.last_char);
            char_scores.resize(end - begin - 1);
            match_ops.resize(end - begin - 1);
            std::transform(align_start + begin, align_start + end - 1,
                           char_scores.begin(),
                           [&row](char c) { return row[c]; });
            std::transform(align_start + begin, align_start + end - 1,
                           match_ops.begin(),
                           [&op_row](char c) { return op_row[c]; });

            compute_match_insert_updates(config_,
                                         update_scores.data() + (begin - overall_begin),
                                         update_prevs.data() + (begin - overall_begin),
                                         update_ops.data() + (begin - overall_begin),
                                         next_pair.first,
                                         incoming.scores.data() + begin,
                                         incoming.ops.data() + begin,
                                         char_scores,
                                         match_ops,
                                         end - begin);

            // compute deletion scores
            score_t delete_score;
            for (size_t i = 1; i < update_scores.size(); ++i) {
                delete_score = update_scores[i - 1]
                    + (update_ops[i - 1] == Cigar::Operator::DELETION
                        ? config_.gap_extension_penalty
                        : config_.gap_opening_penalty);

                if (delete_score > update_scores[i]) {
                    update_ops[i] = Cigar::Operator::DELETION;
                    update_prevs[i] = next_node;
                    update_scores[i] = delete_score;
                }
            }


            // update column scores
            bool updated = false;
            size_t i = 0;
#ifdef __AVX2__
            i = update_column_avx2(updated,
                                   update_scores,
                                   update_ops,
                                   update_prevs,
                                   next_column,
                                   overall_begin);
#endif
            // handle the residual
            for (; i < update_scores.size(); ++i) {
                if (update_scores[i] > next_column.scores[overall_begin + i]) {
                    next_column.ops[overall_begin + i] = update_ops[i];
                    next_column.prev_nodes[overall_begin + i] = update_prevs[i];
                    next_column.scores[overall_begin + i] = update_scores[i];

                    updated = true;
                }
            }


            // put column back in the priority queue if it's updated
            if (updated) {
                // store max pos
                auto max_pos = std::max_element(
                    next_column.scores.begin() + overall_begin,
                    next_column.scores.begin() + overall_end
                );

                next_column.best_pos = max_pos - next_column.scores.begin();

                if (*max_pos > start_score) {
                    start_node = iter->first;
                    start_score = iter->second.best_score();
                    score_cutoff = std::max(start_score, min_path_score);
                }

                // branch and bound
                // TODO: this cuts off too early (before the scores have converged)
                //       so the code below has to be used to compute correct scores
                if (!std::equal(match_score_begin + overall_begin,
                                match_score_begin + overall_end,
                                next_column.scores.begin() + overall_begin,
                                [&](auto a, auto b) { return a + b < score_cutoff; })
                        && score_cutoff - iter->second.best_score() <= config_.xdrop)
                    columns_to_update.emplace(iter->first, iter->second.best_score());
            }
        }
    }

    assert(start_score > config_.min_cell_score);

    // no good path found
    if (start_node == SequenceGraph::npos || start_score == path.get_score() || score_cutoff > start_score)
        return;

    // check to make sure that start_node stores the best starting point
    assert(start_score
        == std::max_element(dp_table.begin(), dp_table.end(),
                            [](const auto &a, const auto &b) {
                                return a.second < b.second;
                            })->second.best_score());

    if (dp_table.find(start_node)->second.best_op() != Cigar::Operator::MATCH)
        logger->trace("best alignment does not end with a MATCH");

    // get all alignments
    dp_table.extract_alignments(*graph_,
                                config_,
                                std::string_view(align_start, query.data() + query.size() - align_start),
                                callback,
                                orientation,
                                min_path_score,
                                &start_node);
}


template <typename NodeType>
void DefaultColumnExtender<NodeType>
::initialize_query(const std::string_view query) {
    partial_sums_.resize(query.size());
    std::transform(query.begin(), query.end(),
                   partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query) == partial_sums_.front());
    assert(config_.get_row(query.back())[query.back()] == partial_sums_.back());
}


template class DefaultColumnExtender<>;
