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
std::pair<typename DPTable<NodeType>::const_iterator, bool> DefaultColumnExtender<NodeType>
::emplace_node(NodeType node, NodeType, char c, size_t size, size_t best_pos) {
    auto find = dp_table.find(node);
    if (find == dp_table.end()) {
        return dp_table.emplace(node,
                                typename DPTable<NodeType>::Column(
                                    size, config_.min_cell_score, c,
                                    best_pos + 1 != size ? best_pos + 1 : best_pos
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
bool DefaultColumnExtender<NodeType>
::extendable(size_t begin, size_t end, const score_t *scores) const {
    // TODO: this cuts off too early (before the scores have converged)
    //       so path scores have to be recomputed after alignment
    return !std::equal(
        match_score_begin + begin, match_score_begin + end,
        scores + begin,
        [&](auto a, auto b) { return a + b < score_cutoff; }
    );
};

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

inline bool compute_match_insert_updates_avx2(size_t &i,
                                              size_t length,
                                              int64_t prev_node,
                                              int32_t gap_opening_penalty,
                                              int32_t gap_extension_penalty,
                                              int32_t *&update_scores,
                                              long long int *&update_prevs,
                                              int32_t *&update_ops,
                                              const int32_t *&incoming_scores,
                                              const int32_t *&incoming_ops,
                                              const int32_t *&char_scores,
                                              const int32_t *&match_ops) {
    bool updated = false;
    const __m256i prev_packed = _mm256_set1_epi64x(prev_node);
    const __m256i ins_packed = _mm256_set1_epi32(Cigar::Operator::INSERTION);
    const __m256i gap_open_packed = _mm256_set1_epi32(gap_opening_penalty);
    const __m256i gap_extend_packed = _mm256_set1_epi32(gap_extension_penalty);
    const __m256i min_score = _mm256_set1_epi32(std::numeric_limits<int32_t>::min());

    for (; i + 8 <= length; i += 8) {
        __m256i incoming_packed = _mm256_loadu_si256((__m256i*)(incoming_scores - 1));

        // compute match and delete scores
        __m256i match_scores_packed = _mm256_add_epi32(
            _mm256_load_si256((__m256i*)char_scores),
            incoming_packed
        );

        __m256i incoming_ops_packed = _mm256_loadu_si256((__m256i*)incoming_ops);
        __m256i incoming_is_insert = _mm256_cmpeq_epi32(incoming_ops_packed, ins_packed);

        __m256i incoming_shift = rshiftpushback_epi32(incoming_packed, incoming_scores[7]);
        __m256i insert_scores_packed = _mm256_add_epi32(
            incoming_shift,
            _mm256_blendv_epi8(gap_open_packed, gap_extend_packed, incoming_is_insert)
        );

        // disallow insertion after deletion
        static_assert(Cigar::Operator::DELETION > Cigar::Operator::INSERTION);
        insert_scores_packed = _mm256_blendv_epi8(
            insert_scores_packed, min_score,
            _mm256_cmpgt_epi32(incoming_ops_packed, ins_packed) // incoming is delete
        );

        // update scores
        __m256i max_scores = _mm256_max_epi32(match_scores_packed, insert_scores_packed);
        __m256i both_cmp = _mm256_cmpgt_epi32(
            max_scores,
            _mm256_loadu_si256((__m256i*)update_scores)
        );

        updated |= static_cast<bool>(_mm256_movemask_epi8(both_cmp));

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

    return updated;
}

#endif

template <typename NodeType,
          typename score_t = typename Alignment<NodeType>::score_t>
inline bool compute_match_insert_updates(const DBGAlignerConfig &config,
                                         score_t *update_scores,
                                         NodeType *update_prevs,
                                         Cigar::Operator *update_ops,
                                         const NodeType &prev_node,
                                         const score_t *incoming_scores,
                                         const Cigar::Operator *incoming_ops,
                                         const AlignedVector<int32_t> &char_scores,
                                         const AlignedVector<Cigar::Operator> &match_ops,
                                         size_t length) {
    bool updated = false;
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
        updated = true;
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
    updated |= compute_match_insert_updates_avx2(
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
    // TODO: is there a better way to handle this?
    for (; i < length; ++i) {
        if (*(incoming_scores - 1) + *char_scores_data > *update_scores) {
            updated = true;
            *update_ops_cast = *match_ops_cast;
            *update_prevs_cast = prev_node;
            *update_scores = *(incoming_scores - 1) + *char_scores_data;
        }

        // disallow INSERTION after DELETION?
        if (*incoming_ops_cast != Cigar::Operator::DELETION) {
            gap_score = *incoming_scores + (*incoming_ops_cast == Cigar::Operator::INSERTION
                ? config.gap_extension_penalty : config.gap_opening_penalty);

            if (gap_score > *update_scores) {
                updated = true;
                *update_ops_cast = Cigar::Operator::INSERTION;
                *update_prevs_cast = prev_node;
                *update_scores = gap_score;
            }
        }

        ++incoming_scores;
        ++incoming_ops_cast;
        ++char_scores_data;
        ++match_ops_cast;
        ++update_scores;
        ++update_prevs_cast;
        ++update_ops_cast;
    }

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

    columns_to_update.emplace(start_node, start_score);

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
        if (next_pair.second != cur_col.best_score()) {
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
    for (const auto &[next_node, c] : out_columns) {
        const auto &cur_col = dp_table.find(incoming_node)->second;
        auto iter = &*emplace_node(next_node, incoming_node,
                                   c, size, cur_col.best_pos).first;
        auto &next_column = const_cast<typename DPTable<NodeType>::Column&>(
            iter->second
        );

        // the get_outgoing_columns call may have invalidated cur_col, so
        // get a pointer to the column again
        const auto &incoming = dp_table.find(incoming_node)->second;

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


        bool updated = false;

        const auto &row = config_.get_row(next_column.last_char);
        const auto &op_row = Cigar::get_op_row(next_column.last_char);

        // for storage of intermediate values
        AlignedVector<score_t> char_scores(end - begin - 1);
        AlignedVector<Cigar::Operator> match_ops(end - begin - 1);

        // TODO: do this with SIMD gather operations?
        std::transform(align_start + begin, align_start + end - 1,
                       char_scores.begin(),
                       [&row](char c) { return row[c]; });
        std::transform(align_start + begin, align_start + end - 1,
                       match_ops.begin(),
                       [&op_row](char c) { return op_row[c]; });

        updated |= compute_match_insert_updates(
            config_,
            next_column.scores.data() + begin,
            next_column.prev_nodes.data() + begin,
            next_column.ops.data() + begin,
            incoming_node,
            incoming.scores.data() + begin,
            incoming.ops.data() + begin,
            char_scores,
            match_ops,
            end - begin
        );

        // compute deletion scores
        score_t delete_score;
        for (size_t i = overall_begin + 1; i < overall_end; ++i) {
            // disallow DELETION after INSERTION
            if (next_column.ops[i - 1] == Cigar::Operator::INSERTION)
                continue;

            delete_score = next_column.scores[i - 1]
                + (next_column.ops[i - 1] == Cigar::Operator::DELETION
                    ? config_.gap_extension_penalty
                    : config_.gap_opening_penalty);

            if (delete_score > next_column.scores[i]) {
                updated = true;
                next_column.ops[i] = Cigar::Operator::DELETION;
                next_column.prev_nodes[i] = next_node;
                next_column.scores[i] = delete_score;
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
            if (extendable(overall_begin, overall_end, next_column.scores.data()))
                columns_to_update.emplace(iter->first, iter->second.best_score());

            assert(start_score == dp_table.best_score().second);
        }
    }
}


template class DefaultColumnExtender<>;
