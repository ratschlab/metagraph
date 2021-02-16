#include "aligner_extender_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include <priority_deque.hpp>

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

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

#endif


template <typename NodeType>
DefaultColumnExtender<NodeType>::DefaultColumnExtender(const DeBruijnGraph &graph,
                                                       const DBGAlignerConfig &config,
                                                       std::string_view query)
      : graph_(graph), config_(config), query_(query) {
    assert(config_.check_config_scores());
    partial_sums_.reserve(query_.size() + 1);
    partial_sums_.resize(query_.size(), 0);
    std::transform(query_.begin(), query_.end(),
                   partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query_) == partial_sums_.front());
    assert(config_.get_row(query_.back())[query_.back()] == partial_sums_.back());
    partial_sums_.push_back(0);

    for (char c : graph_.alphabet()) {
        auto &p_score_row = profile_score_.emplace(c, query_.size() + 8).first.value();
        auto &p_op_row = profile_op_.emplace(c, query_.size() + 8).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = Cigar::get_op_row(c);

        std::transform(query_.begin(), query_.end(), p_score_row.begin(),
                       [&row](char q) { return row[q]; });
        std::transform(query_.begin(), query_.end(), p_op_row.begin(),
                       [&op_row](char q) { return op_row[q]; });
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &seed) {
    seed_ = &seed;
    reset();
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::operator()(ExtensionCallback callback,
                                                 score_t min_seed_score) {
    const char *align_start = seed_->get_query().data() + seed_->get_query().size() - 1;
    size_t start = align_start - query_.data();
    size_t size = query_.size() - start + 1;
    assert(start + size == partial_sums_.size());
    match_score_begin_ = partial_sums_.data() + start;

    extend_window_ = { align_start, size - 1 };
    assert(extend_window_[0] == seed_->get_query().back());

    size_t size_per_column = sizeof(size_t) * 2 + sizeof(Column) + size * (
        sizeof(score_t) * 3 + sizeof(AlignNode) * 2 + sizeof(Cigar::Operator) * 3
    );

    constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;

    assert(seed_->get_cigar().back().first == Cigar::MATCH
        || seed_->get_cigar().back().first == Cigar::MISMATCH);

    auto &first_column = table_.emplace(
        graph_.max_index() + 1,
        Column{ { std::make_tuple(
            ScoreVec(1, ninf), ScoreVec(1, ninf), ScoreVec(1, ninf),
            OpVec(1, Cigar::CLIPPED), OpVec(1, Cigar::CLIPPED), OpVec(1, Cigar::CLIPPED),
            AlignNode{}, PrevVec(1, NONE), PrevVec(1, NONE),
            0 /* offset */, 0 /* max_pos */
        )}, false }
    ).first.value().first[0];
    sanitize(first_column);
    auto &[S, E, F, OS, OE, OF, prev_node, PS, PF, offset, max_pos] = first_column;

    size_t num_columns = 1;

    S[0] = seed_->get_score() - profile_score_[seed_->get_sequence().back()][start];

    AlignNode start_node{ graph_.max_index() + 1,
                          seed_->get_sequence()[seed_->get_sequence().size() - 2],
                          0 };

    typedef std::pair<AlignNode, score_t> Ref;
    boost::container::priority_deque<Ref, std::vector<Ref>, utils::LessSecond> best_starts;
    best_starts.emplace(start_node, S[0]);

    std::priority_queue<Ref, std::vector<Ref>, utils::LessSecond> stack;
    stack.emplace(start_node, 0);

    while (stack.size()) {
        AlignNode prev = stack.top().first;
        stack.pop();

        size_t total_size = num_columns * size_per_column
            + table_.size() * sizeof(NodeType);

        if (total_size > config_.max_ram_per_alignment * 1000000)
            break;

        if (num_columns > config_.max_nodes_per_seq_char * extend_window_.size())
            break;

        for (const auto &[next, c] : get_outgoing(prev)) {
            auto &column_pair = table_[next];
            auto &[column, converged] = column_pair;
            if (converged)
                continue;

            auto &column_prev = table_[std::get<0>(prev)].first;

            score_t xdrop_cutoff = (best_starts.size()
                ? best_starts.maximum().second
                : seed_->get_score()) - config_.xdrop;

            // compute bandwidth based on xdrop criterion
            size_t min_i;
            size_t max_i;
            {
                const auto &S_prev = std::get<0>(column_prev[std::get<2>(prev)]);
                size_t offset_prev = std::get<9>(column_prev[std::get<2>(prev)]);
                size_t max_pos_prev = std::get<10>(column_prev[std::get<2>(prev)]);
                assert(max_pos_prev - offset_prev < S_prev.size());
                assert(std::max_element(S_prev.begin(), S_prev.end())
                    == S_prev.begin() + (max_pos_prev - offset_prev));

                if (S_prev[max_pos_prev - offset_prev] < xdrop_cutoff)
                    continue;

                min_i = std::max((size_t)1, max_pos_prev);
                max_i = std::min({ min_i + 2, S_prev.size() + offset_prev + 2, size });
                while (min_i >= std::max((size_t)2, offset_prev)
                        && S_prev[min_i - offset_prev] >= xdrop_cutoff) {
                    --min_i;
                }
                while (max_i - offset_prev < S_prev.size()
                        && S_prev[max_i - offset_prev] >= xdrop_cutoff) {
                    ++max_i;
                }

                if (max_i - offset_prev >= S_prev.size())
                    max_i = size;

                if (min_i >= max_i)
                    continue;
            }

            size_t depth = column.size();
            size_t cur_size = max_i - min_i + 1;

            auto &next_column = column.emplace_back(
                ScoreVec(cur_size, ninf), ScoreVec(cur_size, ninf),
                ScoreVec(cur_size, ninf),
                OpVec(cur_size, Cigar::CLIPPED), OpVec(cur_size, Cigar::CLIPPED),
                OpVec(cur_size, Cigar::CLIPPED),
                prev, PrevVec(cur_size, NONE), PrevVec(cur_size, NONE),
                min_i - 1 /* offset */, 0 /* max_pos */
            );
            sanitize(next_column);
            auto &[S, E, F, OS, OE, OF, prev_node, PS, PF, offset, max_pos] = next_column;
            ++num_columns;

            assert(cur_size + offset <= size);

            AlignNode cur{ next, c, depth };

            auto &[S_prev, E_prev, F_prev, OS_prev, OE_prev, OF_prev,
                   prev_node_prev, PS_prev, PF_prev, offset_prev, max_pos_prev]
                = column_prev[std::get<2>(prev)];
            assert(S_prev.size() + offset_prev <= size);

            // compute column boundaries for updating the match and deletion scores
            // need i + offset - offset_prev - 1 >= 0
            // &&   i + offset - offset_prev - 1 < S_prev.size()
            size_t match_begin = offset_prev + 1 >= offset ? offset_prev + 1 - offset : 0;
            size_t match_end = S_prev.size() + offset_prev + 1 >= offset
                ? std::min(S_prev.size() + offset_prev + 1 - offset, cur_size)
                : 0;
            assert(match_begin + offset);

            // need i + offset - offset_prev >= 0
            // &&   i + offset - offset_prev < S_prev.size()
            size_t del_begin = match_begin ? match_begin - 1 : 0;
            size_t del_end = S_prev.size() + offset_prev >= offset
                ? std::min(S_prev.size() + offset_prev - offset, cur_size)
                : 0;

            assert(del_begin <= match_begin);
            assert(match_begin - del_begin <= 1);
            assert(match_end == std::min(cur_size, del_end + 1));

            // set prev node vector for deletion
            std::fill(PF.begin() + del_begin, PF.begin() + del_end, PREV);

#ifdef __AVX2__
            // compute deletion score vector
            for (size_t i = del_begin; i < del_end; i += 8) {
                __m256i del_open = _mm256_add_epi32(
                    _mm256_loadu_si256((__m256i*)&S_prev[i + offset - offset_prev]),
                    _mm256_set1_epi32(config_.gap_opening_penalty)
                );
                __m256i del_extend = _mm256_add_epi32(
                    _mm256_loadu_si256((__m256i*)&F_prev[i + offset - offset_prev]),
                    _mm256_set1_epi32(config_.gap_extension_penalty)
                );

                _mm256_storeu_si256((__m256i*)&F[i],
                                    _mm256_max_epi32(del_open, del_extend));

                __m128i del_op_v = _mm_blendv_epi8(
                    _mm_set1_epi8(Cigar::MATCH),
                    _mm_set1_epi8(Cigar::DELETION),
                    mm256_cvtepi32_epi8(_mm256_cmpgt_epi32(del_extend, del_open))
                );

                mm_storeu_si64(&OF[i], del_op_v);
            }

            // compute initial match/mismatch scores
            for (size_t i = match_begin; i < match_end; i += 8) {
                __m256i s_prev_v = _mm256_loadu_si256(
                    (__m256i*)&S_prev[i + offset - offset_prev - 1]
                );
                __m256i profile_v = _mm256_cvtepi8_epi32(
                    mm_loadu_si64(&profile_score_[c][start + i + offset - 1])
                );

                _mm256_storeu_si256((__m256i*)&S[i],
                                    _mm256_add_epi32(s_prev_v, profile_v));
            }
#else
            // compute deletion score vector
            for (size_t i = del_begin; i < del_end; ++i) {
                score_t del_open = S_prev[i + offset - offset_prev]
                                        + config_.gap_opening_penalty;
                score_t del_extend = F_prev[i + offset - offset_prev]
                                        + config_.gap_extension_penalty;
                F[i] = std::max(del_open, del_extend);
                OF[i] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;
            }

            // compute initial match/mismatch scores
            for (size_t i = match_begin; i < match_end; ++i) {
                S[i] = S_prev[i + offset - offset_prev - 1]
                        + profile_score_[c][start + i + offset - 1];
            }
#endif

            // compute insert and best scores
            bool updated = false;
            for (size_t i = 0; i < cur_size; ++i) {
                if (i) {
                    score_t ins_open = S[i - 1] + config_.gap_opening_penalty;
                    score_t ins_extend = E[i - 1] + config_.gap_extension_penalty;
                    E[i] = std::max(ins_open, ins_extend);
                    OE[i] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;
                } else {
                    E[i] = ninf;
                    OE[i] = Cigar::CLIPPED;
                }

                S[i] = std::max({ 0, S[i], E[i], F[i] });

                if (S[i] >= xdrop_cutoff) {
                    xdrop_cutoff = std::max(xdrop_cutoff, S[i] - config_.xdrop);
                    updated = true;
                } else {
                    S[i] = ninf;
                }
            }

            // compute traceback vectors
#ifdef __AVX2__
            for (size_t i = 0; i < cur_size; i += 8) {
                __m256i e_v = _mm256_load_si256((__m256i*)&E[i]);
                __m256i f_v = _mm256_load_si256((__m256i*)&F[i]);
                __m256i s_v = _mm256_load_si256((__m256i*)&S[i]);
                __m128i mask = mm256_cvtepi32_epi8(
                    _mm256_cmpgt_epi32(s_v, _mm256_setzero_si256())
                );
                __m128i equal_e = mm256_cvtepi32_epi8(_mm256_cmpeq_epi32(s_v, e_v));
                __m128i equal_f = mm256_cvtepi32_epi8(_mm256_cmpeq_epi32(s_v, f_v));

                __m128i ps_v = _mm_blendv_epi8(_mm_set1_epi8(PREV),
                                               _mm_set1_epi8(CUR),
                                               equal_e);
                mm_maskstorel_epi8((int8_t*)&PS[i], mask, ps_v);

                __m128i os_v = _mm_blendv_epi8(
                    mm_loadu_si64(&profile_op_[c][start + i + offset - 1]),
                    _mm_set1_epi8(Cigar::DELETION),
                    equal_f
                );
                os_v = _mm_blendv_epi8(os_v, _mm_set1_epi8(Cigar::INSERTION), equal_e);
                mm_maskstorel_epi8((int8_t*)&OS[i], mask, os_v);
            }
#else
            for (size_t i = 0; i < cur_size; ++i) {
                if (S[i] > 0) {
                    if (S[i] == E[i]) {
                        PS[i] = CUR;
                        OS[i] = Cigar::INSERTION;
                    } else if (S[i] == F[i]) {
                        PS[i] = PREV;
                        OS[i] = Cigar::DELETION;
                    } else {
                        assert(i + offset >= offset_prev + 1
                            && i + offset - offset_prev - 1 < S_prev.size());
                        PS[i] = PREV;
                        OS[i] = profile_op_[c][start + i + offset - 1];
                    }
                }
            }
#endif

            // extend to the right with insertion scores
            while (offset + S.size() < size && S.back() >= xdrop_cutoff) {
                score_t ins_open = S.back() + config_.gap_opening_penalty;
                score_t ins_extend = E.back() + config_.gap_extension_penalty;
                E.push_back(std::max(ins_open, ins_extend));
                F.push_back(ninf);
                S.push_back(std::max({ 0, E.back(), F.back() }));

                OS.push_back(Cigar::CLIPPED);
                OE.push_back(ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH);
                OF.push_back(Cigar::CLIPPED);

                PS.push_back(NONE);
                PF.push_back(NONE);
                if (S.back() > 0 && S.back() == E.back()) {
                    updated = true;
                    PS.back() = CUR;
                    OS.back() = Cigar::INSERTION;
                }
            }
            sanitize(next_column);

            auto max_it = std::max_element(S.begin(), S.end());
            max_pos = (max_it - S.begin()) + offset;
            assert(max_pos < size);

            converged = !updated || has_converged(column_pair, 1, cur_size);

            if (best_starts.size() < config_.num_alternative_paths) {
                best_starts.emplace(cur, *max_it);
            } else if (*max_it > best_starts.minimum().second) {
                best_starts.update(best_starts.begin(), Ref{ cur, *max_it });
            }

            assert(match_score_begin_[max_pos]
                == config_.match_score(extend_window_.substr(max_pos)));
            score_t score_rest = *max_it + match_score_begin_[max_pos];

            assert(xdrop_cutoff == best_starts.maximum().second - config_.xdrop);

            if (*max_it >= xdrop_cutoff && score_rest >= min_seed_score)
                stack.emplace(Ref{ cur, *max_it });
        };
    }

    while (best_starts.size()) {
        auto [best_node, max_score] = best_starts.maximum();
        best_starts.pop_maximum();

        const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
            = table_[std::get<0>(best_node)].first[std::get<2>(best_node)];

        assert(S[max_pos - offset] == max_score);

        if (max_pos < 2 && std::get<0>(best_node) == seed_->back()
                && !std::get<2>(best_node)) {
            callback(Alignment<NodeType>(), 0);
            return;
        }

        Cigar cigar;
        if (max_pos + 1 < size)
            cigar.append(Cigar::CLIPPED, size - max_pos - 1);

        size_t pos = max_pos;
        std::vector<NodeType> path;
        std::string seq;
        NodeType start_node = 0;
        score_t score = max_score;

        Cigar::Operator last_op = std::get<3>(
            table_[std::get<0>(best_node)].first[std::get<2>(best_node)])[pos - offset];
        assert(last_op == Cigar::MATCH);

        while (true) {
            const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
                = table_[std::get<0>(best_node)].first[std::get<2>(best_node)];

            assert(last_op == Cigar::MATCH || last_op == Cigar::MISMATCH);

            if (pos == 1 && std::get<0>(best_node) == seed_->back()
                    && !std::get<2>(best_node)
                    && OS[pos - offset] == seed_->get_cigar().back().first) {
                assert(std::get<0>(prev) == graph_.max_index() + 1);
                start_node = seed_->back();
                score -= seed_->get_score();
                break;
            } else if (OS[pos - offset] == Cigar::CLIPPED) {
                start_node = DeBruijnGraph::npos;
                break;
            }

            last_op = OS[pos - offset];

            switch (OS[pos - offset]) {
                case Cigar::MATCH:
                case Cigar::MISMATCH: {
                    cigar.append(OS[pos - offset]);
                    path.push_back(std::get<0>(best_node));
                    seq += std::get<1>(best_node);
                    assert((OS[pos - offset] == Cigar::MATCH)
                        == (graph_.get_node_sequence(std::get<0>(best_node)).back()
                            == extend_window_[pos - 1]));
                    switch (PS[pos - offset]) {
                        case NONE: { best_node = {}; } break;
                        case PREV: { best_node = prev; } break;
                        case CUR: {}
                    }
                    --pos;
                } break;
                case Cigar::INSERTION: {
                    assert(PS[pos - offset] == CUR);
                    while (last_op == Cigar::INSERTION) {
                        last_op = OE[pos - offset];
                        assert(last_op == Cigar::MATCH || last_op == Cigar::INSERTION);
                        cigar.append(Cigar::INSERTION);
                        --pos;
                        assert(pos);
                    }
                } break;
                case Cigar::DELETION: {
                    while (last_op == Cigar::DELETION) {
                        const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
                            = table_[std::get<0>(best_node)].first[std::get<2>(best_node)];
                        last_op = OF[pos - offset];
                        assert(last_op == Cigar::MATCH || last_op == Cigar::DELETION);
                        path.push_back(std::get<0>(best_node));
                        seq += std::get<1>(best_node);
                        cigar.append(Cigar::DELETION);
                        switch (PF[pos - offset]) {
                            case NONE: { best_node = {}; } break;
                            case PREV: { best_node = prev; } break;
                            case CUR: {}
                        }
                    }
                } break;
                case Cigar::CLIPPED: { assert(false); }
            }

            assert(pos);
        }

        if (max_score < min_seed_score)
            return;

        if (pos > 1)
            cigar.append(Cigar::CLIPPED, pos - 1);

        std::reverse(cigar.begin(), cigar.end());
        std::reverse(path.begin(), path.end());
        std::reverse(seq.begin(), seq.end());

        Alignment<NodeType> extension({ extend_window_.data() + pos, max_pos - pos },
                                      std::move(path),
                                      std::move(seq),
                                      score,
                                      std::move(cigar),
                                      0,
                                      seed_->get_orientation(),
                                      graph_.get_k() - 1);

        assert(extension.is_valid(graph_, &config_));
        callback(std::move(extension), start_node);
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::call_visited_nodes(const std::function<void(NodeType, size_t, size_t)> &callback) const {
    size_t window_start = extend_window_.data() - query_.data();
    for (const auto &[node, columns] : table_) {
        size_t start = query_.size();
        size_t end = 0;
        for (const auto &column : columns.first) {
            const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos] = column;

            auto it = std::find_if(S.begin(), S.end(), [](score_t s) { return s > 0; });
            size_t start_c = (it - S.begin()) + offset;
            start = std::min(start, start_c ? start_c - 1 : start_c);
            end = std::max(end, max_pos ? max_pos - 1 : max_pos);
        }

        if (start < end)
            callback(node, window_start + start, window_start + end);
    }
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>
::has_converged(const Column &column, size_t begin, size_t end) {
    if (column.second)
        return true;

    if (column.first.size() <= 1)
        return false;

    const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos] = column.first.back();
    const auto &[S_b, E_b, F_b, OS_b, OE_b, OF_b, prev_b, PS_b, PF_b, offset_b, max_pos_b]
        = column.first[column.first.size() - 2];

    return offset == offset_b && max_pos == max_pos_b
        && std::get<0>(prev) == std::get<0>(prev_b) && S.size() == S_b.size()
        && std::equal(S.begin() + begin, S.begin() + end, S_b.begin() + begin)
        && std::equal(E.begin() + begin, E.begin() + end, E_b.begin() + begin)
        && std::equal(F.begin() + begin, F.begin() + end, F_b.begin() + begin)
        && std::equal(OS.begin() + begin, OS.begin() + end, OS_b.begin() + begin)
        && std::equal(OE.begin() + begin, OE.begin() + end, OE_b.begin() + begin)
        && std::equal(OF.begin() + begin, OF.begin() + end, OF_b.begin() + begin)
        && std::equal(PS.begin() + begin, PS.begin() + end, PS_b.begin() + begin)
        && std::equal(PF.begin() + begin, PF.begin() + end, PF_b.begin() + begin);
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::sanitize(Scores &scores) {
    auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos] = scores;

    size_t size = S.size();
    size_t pad_size = ((size + 7) / 8) * 8 + 8;
    size_t size_diff = pad_size - size;

    if (!size_diff)
        return;

    S.reserve(pad_size);
    E.reserve(pad_size);
    F.reserve(pad_size);
    OS.reserve(pad_size);
    OE.reserve(pad_size);
    OF.reserve(pad_size);
    PS.reserve(pad_size);
    PF.reserve(pad_size);

    memset(&S[size], 0, sizeof(typename decltype(S)::value_type) * size_diff);
    memset(&E[size], 0, sizeof(typename decltype(E)::value_type) * size_diff);
    memset(&F[size], 0, sizeof(typename decltype(F)::value_type) * size_diff);
    memset(&OS[size], 0, sizeof(typename decltype(OS)::value_type) * size_diff);
    memset(&OE[size], 0, sizeof(typename decltype(OE)::value_type) * size_diff);
    memset(&OF[size], 0, sizeof(typename decltype(OF)::value_type) * size_diff);
    memset(&PS[size], 0, sizeof(typename decltype(PS)::value_type) * size_diff);
    memset(&PF[size], 0, sizeof(typename decltype(PF)::value_type) * size_diff);
}

template <typename NodeType>
std::vector<std::pair<NodeType, char>> DefaultColumnExtender<NodeType>
::get_outgoing(const AlignNode &node) const {
    std::vector<std::pair<NodeType, char>> outgoing;
    if (std::get<0>(node) == graph_.max_index() + 1) {
        outgoing.emplace_back(seed_->back(), seed_->get_sequence().back());
    } else {
        graph_.call_outgoing_kmers(std::get<0>(node), [&](NodeType next, char c) {
            if (c != boss::BOSS::kSentinel)
                outgoing.emplace_back(next, c);
        });
    }

    return outgoing;
}

template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
