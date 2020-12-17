#include "aligner_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "common/logger.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

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

        std::transform(query.begin(), query.end(), profile_score_row.begin(),
                       [&row](char q) { return row[q]; });
        std::transform(query.begin(), query.end(), profile_op_row.begin(),
                       [&op_row](char q) { return op_row[q]; });
    }
}

template <class DPTable, class ColumnIt, class Config>
std::pair<size_t, size_t> get_column_boundaries(DPTable &dp_table,
                                                ColumnIt column_it,
                                                DBGAlignerConfig::score_t xdrop_cutoff,
                                                const Config &config) {
    size_t bandwidth = config.bandwidth;
    auto &column = column_it.value();
    auto &scores = column.scores;
    size_t size = column.size();
    size_t best_pos = column.best_pos;
    assert(best_pos >= column.start_index);
    assert(best_pos - column.start_index < scores.size());

    size_t shift = column.start_index;
    if (xdrop_cutoff > scores[best_pos - shift])
        return std::make_pair(best_pos, best_pos);

    size_t begin = best_pos >= bandwidth ? best_pos - bandwidth : 0;
    size_t end = bandwidth <= size - best_pos ? best_pos + bandwidth : size;

    if (column.min_score_ < xdrop_cutoff) {
        begin = std::max(begin, column.start_index);
        assert(column.start_index + scores.size() >= 8);
        size_t cur_end = std::min(end, column.start_index + scores.size() - 8);
        assert(begin < cur_end);
        if (scores.at(cur_end - 1 - shift) >= xdrop_cutoff) {
            size_t old_size = scores.size();
            dp_table.expand_to_cover(column_it, begin, end);
            assert(shift == column.start_index);
            if (scores.size() > old_size) {
                update_del_scores(config,
                                  scores.data() + begin - shift,
                                  column.prev_nodes.data() + begin - shift,
                                  column.ops.data() + begin - shift,
                                  nullptr,
                                  end - begin,
                                  xdrop_cutoff);
            }
        } else {
            end = cur_end;
        }
    } else {
        size_t old_size = scores.size();
        dp_table.expand_to_cover(column_it, begin, end);
        shift = column.start_index;
        if (scores.size() > old_size) {
            update_del_scores(config,
                              scores.data() + begin - shift,
                              column.prev_nodes.data() + begin - shift,
                              column.ops.data() + begin - shift,
                              nullptr,
                              end - begin,
                              xdrop_cutoff);
        }
    }

    while (begin < end && scores[begin - shift] < xdrop_cutoff) {
        ++begin;
    }

    if (begin == scores.size() + shift)
        return std::make_pair(begin, begin);

    // ensure that the next position is included in the range [begin, end)
    size_t cur_end = best_pos + 2;
    while (cur_end < end && scores[cur_end - shift] >= xdrop_cutoff) {
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
    dp_table.expand_to_cover(column_it, begin, end);

    assert(begin < end);
    assert(begin <= best_pos);
    if (end <= best_pos)
        std::cerr << begin << " " << best_pos << " " << end << std::endl;
    assert(end > best_pos);

    return std::make_pair(begin, end);
}

template <typename NodeType>
std::pair<typename DPTable<NodeType>::iterator, bool> DefaultColumnExtender<NodeType>
::emplace_node(NodeType node, NodeType, char c, size_t size,
               size_t best_pos, size_t last_priority_pos, size_t begin, size_t end) {
    end = std::min(end, size);
    auto find = dp_table.find(node);
    if (find == dp_table.end()) {
        return dp_table.emplace(node,
                                typename DPTable<NodeType>::Column(
                                    size, config_.min_cell_score, c,
                                    best_pos + 1 != size ? best_pos + 1 : best_pos,
                                    last_priority_pos, begin, end
                                ));
    } else {
        dp_table.expand_to_cover(find, begin, end);
        auto [node_begin, node_end] = get_column_boundaries(
            dp_table, find, xdrop_cutoff, config_
        );

        if (node_begin != node_end)
            overlapping_range_ |= (begin < std::min(size, node_end + 9)
                && std::min(size, end + 9)> node_begin);
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
                          column.scores.data() + begin - column.start_index,
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

template <typename score_t>
inline void update_del_scores(const DBGAlignerConfig &config,
                              score_t *update_scores,
                              uint8_t *update_prevs,
                              Cigar::Operator *update_ops,
                              int8_t *updated_mask,
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
                update_prevs[i] = 0xFF;

                if (updated_mask)
                    updated_mask[i] = 0xFF;

                del_score += config.gap_extension_penalty;
                ++i;
            }
            --i;
        }
    }
}

#ifdef __AVX2__

// Drop-in replacement for _mm_loadu_si64
inline __m128i mm_loadu_si64(const void *mem_addr) {
    return _mm_loadl_epi64((const __m128i*)mem_addr);
}

// Drop-in replacement for _mm_storeu_si64
inline void mm_storeu_si64(void *mem_addr, __m128i a) {
    _mm_storel_epi64((__m128i*)mem_addr, a);
}

inline void mm_maskstorel_epi8(int8_t *mem_addr, __m128i mask, __m128i a) {
    __m128i orig = mm_loadu_si64((__m128i*)mem_addr);
    a = _mm_blendv_epi8(orig, a, mask);
    mm_storeu_si64(mem_addr, a);
}

inline void compute_HE_avx2(size_t length,
                            __m128i prev_node,
                            __m256i gap_opening_penalty,
                            __m256i gap_extension_penalty,
                            int32_t *update_scores,
                            int32_t *update_gap_scores,
                            uint8_t *update_prevs,
                            int8_t *update_ops,
                            int32_t *update_gap_count,
                            uint8_t *update_gap_prevs,
                            const int32_t *incoming_scores,
                            const int32_t *incoming_gap_scores,
                            const int8_t *profile_scores,
                            const int8_t *profile_ops,
                            const int32_t *incoming_gap_count,
                            int8_t *updated_mask,
                            __m256i xdrop_cutoff) {
    assert(update_scores != incoming_scores);
    assert(update_gap_scores != incoming_gap_scores);

    __m128i insert_p = _mm_set1_epi8(Cigar::INSERTION);
    for (size_t i = 1; i < length; i += 8) {
        // store score updates
        // load previous values for cells to update
        __m256i H_orig = _mm256_loadu_si256((__m256i*)&update_scores[i]);

        // compute match score
        __m256i incoming_p = _mm256_loadu_si256((__m256i*)&incoming_scores[i - 1]);
        __m256i match_score = _mm256_add_epi32(
            incoming_p, _mm256_cvtepi8_epi32(mm_loadu_si64(&profile_scores[i]))
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
        __m128i update_gap_prev = prev_node;

        // compute updated gap size count
        __m256i is_extend = _mm256_cmpeq_epi32(update_score, update_score_extend);
        __m256i incoming_count = _mm256_add_epi32(
            _mm256_set1_epi32(1),
            _mm256_and_si256(_mm256_loadu_si256((__m256i*)&incoming_gap_count[i]), is_extend)
        );

        __m256i update_gap_scores_orig = _mm256_loadu_si256((__m256i*)&update_gap_scores[i]);
        __m256i update_gap_count_orig = _mm256_loadu_si256((__m256i*)&update_gap_count[i]);
        __m128i update_gap_prevs_orig = mm_loadu_si64(&update_gap_prevs[i]);
        __m256i gap_updated = _mm256_cmpgt_epi32(update_score, update_gap_scores_orig);
        __m128i gap_updated_small = mm256_cvtepi32_epi8(gap_updated);

        update_score = _mm256_blendv_epi8(update_gap_scores_orig, update_score, gap_updated);
        incoming_count = _mm256_blendv_epi8(update_gap_count_orig, incoming_count, gap_updated);
        update_gap_prev = _mm_blendv_epi8(update_gap_prevs_orig, update_gap_prev, gap_updated_small);

        // compute score for cell update. check if inserting a gap improves the update
        __m256i update_cmp = _mm256_cmpgt_epi32(update_score, H);
        H = _mm256_max_epi32(H, update_score);

        // determine which indices satisfy the x-drop criteria
        __m256i xdrop_cmp = _mm256_cmpgt_epi32(H, xdrop_cutoff);

        // revert values not satisfying the x-drop criteria
        H = _mm256_blendv_epi8(H_orig, H, xdrop_cmp);

        __m256i both_cmp = _mm256_cmpgt_epi32(H, H_orig);
        if (!_mm256_movemask_epi8(both_cmp))
            continue;

        __m128i both_cmp_small = mm256_cvtepi32_epi8(both_cmp);
        __m128i update_cmp_small = mm256_cvtepi32_epi8(update_cmp);

        // update scores
        _mm256_maskstore_epi32(&update_scores[i], both_cmp, H);

        _mm256_storeu_si256((__m256i*)&update_gap_scores[i],
                            _mm256_blendv_epi8(update_gap_scores_orig, update_score, both_cmp));
        _mm256_storeu_si256((__m256i*)&update_gap_count[i],
                            _mm256_blendv_epi8(update_gap_count_orig, incoming_count, both_cmp));

        update_gap_prev = _mm_blendv_epi8(update_gap_prevs_orig, update_gap_prev, both_cmp_small);
        mm_storeu_si64(&update_gap_prevs[i], update_gap_prev);

        __m128i update_op = _mm_blendv_epi8(mm_loadu_si64(&profile_ops[i]), insert_p, update_cmp_small);
        mm_maskstorel_epi8(&update_ops[i], both_cmp_small, update_op);

        __m128i updated_mask_orig = mm_loadu_si64(&updated_mask[i]);
        mm_storeu_si64(&updated_mask[i], _mm_or_si128(updated_mask_orig, both_cmp_small));

        __m128i update_prev = _mm_blendv_epi8(prev_node, update_gap_prev, update_cmp_small);
        mm_maskstorel_epi8((int8_t*)&update_prevs[i], both_cmp_small, update_prev);
    }
}

#endif

// direct translation of compute_HE_avx2 to scalar code
inline void compute_HE(size_t length,
                       uint8_t prev_node,
                       int32_t gap_opening_penalty,
                       int32_t gap_extension_penalty,
                       int32_t *update_scores,
                       int32_t *update_gap_scores,
                       uint8_t *update_prevs,
                       Cigar::Operator *update_ops,
                       int32_t *update_gap_count,
                       uint8_t *update_gap_prevs,
                       const int32_t *incoming_scores,
                       const int32_t *incoming_gap_scores,
                       const int8_t *profile_scores,
                       const Cigar::Operator *profile_ops,
                       const int32_t *incoming_gap_count,
                       int8_t *updated_mask,
                       int32_t xdrop_cutoff) {
    // round to cover the same address space as the AVX2 version
    for (size_t j = 1; j < length; j += 8) {
        for (size_t i = j; i < j + 8; ++i) {
            // store score updates
            // load previous values for cells to update
            int32_t H_orig = update_scores[i];

            // compute match score
            int32_t incoming_p = incoming_scores[i - 1];
            int32_t match_score = incoming_p + profile_scores[i];

            // compute score for cell update
            int32_t H = std::max(H_orig, match_score);

            // compute insert score
            int32_t update_score_open = incoming_scores[i] + gap_opening_penalty;
            int32_t update_score_extend = incoming_gap_scores[i] + gap_extension_penalty;
            int32_t update_score = std::max(update_score_open, update_score_extend);
            int8_t update_gap_prev = prev_node;

            // compute updated gap size count
            int32_t incoming_count = 1 + (update_score == update_score_extend
                ? incoming_gap_count[i]
                : 0);

            if (update_score <= update_gap_scores[i]) {
                update_score = update_gap_scores[i];
                incoming_count = update_gap_count[i];
                update_gap_prev = update_gap_prevs[i];
            }

            // compute score for cell update. check if inserting a gap improves the update
            int32_t update_cmp = update_score > H ? 0xFFFFFFFF : 0x0;
            H = std::max(H, update_score);

            // determine which indices satisfy the x-drop criteria
            int32_t xdrop_cmp = H > xdrop_cutoff ? 0xFFFFFFFF : 0x0;

            // revert values not satisfying the x-drop criteria
            if (!xdrop_cmp)
                H = H_orig;

            int32_t both_cmp = H > H_orig ? 0xFFFFFFFF : 0x0;
            if (!both_cmp)
                continue;

            // update scores
            update_scores[i] = H;
            update_gap_scores[i] = update_score;
            update_gap_count[i] = incoming_count;
            update_gap_prevs[i] = update_gap_prev;

            update_ops[i] = update_cmp ? Cigar::INSERTION : profile_ops[i];
            updated_mask[i] = 0xFF;
            update_prevs[i] = update_cmp ? update_gap_prev : prev_node;
        }
    }
}

template <typename NodeType,
          typename Column,
          typename score_t = typename Alignment<NodeType>::score_t>
inline void compute_updates(Column &update_column,
                            size_t begin,
                            const DBGAlignerConfig &config,
                            const NodeType &prev_node,
                            const NodeType &node,
                            const score_t *incoming_scores,
                            const score_t *incoming_gap_scores,
                            const int32_t *incoming_gap_count,
                            const int8_t *profile_scores,
                            const Cigar::Operator *profile_ops,
                            AlignedVector<int8_t> &updated_mask,
                            size_t length,
                            score_t xdrop_cutoff) {
    assert(length);
    size_t shift = update_column.start_index;
    score_t *update_scores = update_column.scores.data() + begin - shift;
    score_t *update_gap_scores = update_column.gap_scores.data() + begin - shift;
    uint8_t *update_prevs = update_column.prev_nodes.data() + begin - shift;
    int32_t *update_gap_count = update_column.gap_count.data() + begin - shift;
    uint8_t *update_gap_prevs = update_column.gap_prev_nodes.data() + begin - shift;
    Cigar::Operator *update_ops = update_column.ops.data() + begin - shift;

    // handle first element (i.e., no match update possible)
    score_t update_score = std::max(
        incoming_scores[0] + config.gap_opening_penalty,
        incoming_gap_scores[0] + config.gap_extension_penalty
    );

    uint8_t prev_node_rank = prev_node != node
        ? update_column.rank_prev_node(prev_node)
        : 0xFF;

    if (update_score >= xdrop_cutoff && update_score > update_scores[0]) {
        update_scores[0] = update_score;
        update_ops[0] = Cigar::INSERTION;
        update_prevs[0] = prev_node_rank;
        updated_mask[0] = 0xFF;

        update_gap_scores[0] = update_score;
        update_gap_prevs[0] = prev_node_rank;
        update_gap_count[0] = update_score == incoming_gap_scores[0] + config.gap_extension_penalty
            ? incoming_gap_count[0] + 1
            : 1;
    }

    // ensure sizes are the same before casting
    static_assert(sizeof(Cigar::Operator) == sizeof(int8_t));

    auto update_block = [&]() {
        compute_HE(length, prev_node_rank,
                   config.gap_opening_penalty, config.gap_extension_penalty,
                   update_scores, update_gap_scores,
                   update_prevs,
                   update_ops, update_gap_count,
                   update_gap_prevs,
                   incoming_scores, incoming_gap_scores,
                   profile_scores, profile_ops,
                   incoming_gap_count, updated_mask.data(), xdrop_cutoff - 1);
    };

#ifdef __AVX2__

    if (prev_node != node) {
        // update 8 scores at a time
        compute_HE_avx2(length,
                        _mm_set1_epi8(prev_node_rank),
                        _mm256_set1_epi32(config.gap_opening_penalty),
                        _mm256_set1_epi32(config.gap_extension_penalty),
                        update_scores, update_gap_scores,
                        update_prevs,
                        reinterpret_cast<int8_t*>(update_ops),
                        update_gap_count,
                        update_gap_prevs,
                        incoming_scores, incoming_gap_scores,
                        profile_scores, reinterpret_cast<const int8_t*>(profile_ops),
                        incoming_gap_count, updated_mask.data(),
                        _mm256_set1_epi32(xdrop_cutoff - 1));
    } else {
        update_block();
    }

#else

    update_block();

#endif

    update_del_scores(config, update_scores, update_prevs, update_ops,
                      updated_mask.data(), length, xdrop_cutoff);
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
        size_t shift = find->second.start_index;
        if (query_clipping - dp_table.get_query_offset() >= shift
                && find->second.scores[query_clipping - dp_table.get_query_offset() - shift]
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

    max_num_nodes = std::numeric_limits<size_t>::max();
    if (config_.max_nodes_per_seq_char < std::numeric_limits<double>::max()) {
        max_num_nodes = std::min(max_num_nodes,
            static_cast<size_t>(std::ceil(config_.max_nodes_per_seq_char
                * static_cast<double>(size)))
        );
    }

    double num_bytes = static_cast<double>(dp_table.num_bytes()) / 1024 / 1024;
    if (num_bytes > config_.max_ram_per_alignment)
        logger->warn("Alignment RAM limit too low. Alignment may be fragmented.");

    extend_main(callback, min_path_score);
}

template <typename NodeType>
std::deque<std::pair<NodeType, char>> DefaultColumnExtender<NodeType>
::fork_extension(NodeType node,
                 std::function<void(DBGAlignment&&, NodeType)>,
                 score_t) {
    overlapping_range_ = false;
    std::deque<std::pair<DeBruijnGraph::node_index, char>> out_columns;

    double num_bytes = static_cast<double>(dp_table.num_bytes()) / 1024 / 1024;
    bool added = false;

    if (dynamic_cast<const DBGSuccinct*>(graph_)) {
        const auto &dbg_succ = *dynamic_cast<const DBGSuccinct*>(graph_);
        const auto &boss = dbg_succ.get_boss();
        graph_->adjacent_outgoing_nodes(node, [&](auto next_node) {
            auto find = dp_table.find(next_node);
            char c = find == dp_table.end()
                ? boss.decode(
                      boss.get_W(dbg_succ.kmer_to_boss_index(next_node)) % boss.alph_size
                  )
                : find->second.last_char;
            if (c != '$' && ((dp_table.size() < max_num_nodes
                                && num_bytes <= config_.max_ram_per_alignment)
                    || find != dp_table.end())) {
                if (next_node == node) {
                    out_columns.emplace_front(next_node, c);
                } else {
                    out_columns.emplace_back(next_node, c);
                }
                added = true;
                num_bytes = static_cast<double>(dp_table.num_bytes()) / 1024 / 1024;
            }
        });
    } else {
        graph_->call_outgoing_kmers(node, [&](auto next_node, char c) {
            if (c != '$' && ((dp_table.size() < max_num_nodes
                                && num_bytes <= config_.max_ram_per_alignment)
                    || dp_table.count(next_node))) {
                if (next_node == node) {
                    out_columns.emplace_front(next_node, c);
                } else {
                    out_columns.emplace_back(next_node, c);
                }
                added = true;
                num_bytes = static_cast<double>(dp_table.num_bytes()) / 1024 / 1024;
            }
        });
    }

    if (added && num_bytes > config_.max_ram_per_alignment)
        logger->warn("Alignment RAM limit too low: {} MB > {} MB. Alignment may be fragmented.",
                     num_bytes, config_.max_ram_per_alignment);

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
    auto incoming_find = dp_table.find(incoming_node);

    if (dp_table.size() == 1 && out_columns.size() && out_columns.front().first != incoming_node) {
        dp_table.expand_to_cover(incoming_find, 0, size);
        update_del_scores(config_,
                          incoming_find.value().scores.data(),
                          incoming_find.value().prev_nodes.data(),
                          incoming_find.value().ops.data(),
                          nullptr,
                          size,
                          xdrop_cutoff);
    }

    std::tie(begin, end) = get_column_boundaries(dp_table, incoming_find,
                                                 xdrop_cutoff, config_);

    if (begin >= end)
        return;

    for (const auto &[next_node, c] : out_columns) {
        auto emplace = emplace_node(next_node, incoming_node, c, size, begin, begin, begin, end);

        // emplace_node may have invalidated incoming, so update the iterator
        if (emplace.second)
            incoming_find = dp_table.find(incoming_node);

        auto &incoming = incoming_find.value();
        auto &next_column = emplace.first.value();

        assert(begin >= next_column.start_index);
        assert(begin >= incoming.start_index);
        assert(next_column.start_index + next_column.scores.size() >= end);
        assert(incoming.start_index + incoming.scores.size() >= end);

        // store the mask indicating which cells were updated
        // this is padded to ensure that the vectorized code doesn't access
        // out of bounds
        AlignedVector<int8_t> updated_mask(end - begin + 9, 0x0);

        assert(next_column.scores.size() == next_column.gap_scores.size());
        assert(incoming.scores.size() == incoming.gap_scores.size());

        size_t shift = incoming.start_index;
        compute_updates(
            next_column,
            begin,
            config_,
            incoming_node,
            next_node,
            incoming.scores.data() + begin - shift,
            incoming.gap_scores.data() + begin - shift,
            incoming.gap_count.data() + begin - shift,
            profile_score[next_column.last_char].data() + query.size() - size + begin,
            profile_op[next_column.last_char].data() + query.size() - size + begin,
            updated_mask,
            end - begin,
            xdrop_cutoff
        );

        shift = next_column.start_index;

        // Find the maximum changed value
        const score_t *best_update = nullptr;
        size_t best_pos = -1;
        for (size_t i = begin, j = 0; i < size && j < updated_mask.size(); ++i, ++j) {
            if (updated_mask[j] && (!best_update || next_column.scores[i - shift] > *best_update)) {
                best_update = &next_column.scores[i - shift];
                best_pos = i;
                assert(i >= shift);
                assert(i < next_column.size());
            }
        }

        // put column back in the priority queue if it's updated
        if (best_update) {
            assert(updated_mask[best_pos - begin]);

            next_column.last_priority_pos = best_pos;

            if (*best_update > next_column.best_score()) {
                next_column.best_pos = next_column.last_priority_pos;

                assert(next_column.best_pos >= begin);
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
