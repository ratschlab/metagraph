#include "aligner_extender_methods.hpp"

#include <tsl/hopscotch_set.h>

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {

typedef DBGAlignerConfig::score_t score_t;
constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;


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
        auto &p_score_row = profile_score_.emplace(c, query_.size() + 9).first.value();
        auto &p_op_row = profile_op_.emplace(c, query_.size() + 9).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = kCharToOp[c];

        // the first cell in a DP table row is one position before the last matched
        // character, so we need to shift the indices of profile_score_ and profile_op_
        std::transform(query_.begin(), query_.end(), p_score_row.begin() + 1,
                       [&row](char q) { return row[q]; });

        std::transform(query_.begin(), query_.end(), p_op_row.begin() + 1,
                       [&op_row](char q) { return op_row[q]; });
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &seed) {
    assert(seed.size());
    assert(seed.get_cigar().size());
    assert(seed.get_cigar().back().first == Cigar::MATCH
        || seed.get_cigar().back().first == Cigar::MISMATCH);

    seed_ = &seed;
    reset();
}

template <typename Node, typename Column>
std::pair<size_t, size_t> get_band(const Node &prev,
                                   const Column &column_prev,
                                   score_t xdrop_cutoff) {
    const auto &S_prev = std::get<0>(column_prev[std::get<2>(prev)]);
    size_t offset_prev = std::get<9>(column_prev[std::get<2>(prev)]);
    size_t max_pos_prev = std::get<10>(column_prev[std::get<2>(prev)]);
    assert(max_pos_prev - offset_prev < S_prev.size());
    assert(std::max_element(S_prev.begin(), S_prev.end())
        == S_prev.begin() + (max_pos_prev - offset_prev));

    size_t start_pos = max_pos_prev - offset_prev;
    if (S_prev[start_pos] < xdrop_cutoff)
        return {};

    auto stop = [cutoff=std::max(xdrop_cutoff, ninf)](score_t s) { return s < cutoff; };
    auto min_rit = std::find_if(std::make_reverse_iterator(S_prev.begin() + start_pos),
                                S_prev.rend(), stop);
    auto max_it = std::find_if(S_prev.begin() + start_pos, S_prev.end(), stop);

    return std::make_pair(S_prev.rend() - min_rit + offset_prev,
                          max_it - S_prev.begin() + offset_prev);
}

template <typename NodeType,
          typename Column,
          typename Scores,
          typename ProfileScore,
          typename ProfileOp>
bool update_column(const DeBruijnGraph &graph_,
                   const DBGAlignerConfig &config_,
                   const Column &column_prev,
                   Scores &next_column,
                   char c,
                   size_t start,
                   size_t size,
                   score_t &xdrop_cutoff,
                   const ProfileScore &profile_score_,
                   const ProfileOp &profile_op_,
                   const Alignment<NodeType> &seed_) {
    typedef DefaultColumnExtender<NodeType> Extender;

    auto &[S, E, F, OS, OE, OF, prev_node, PS, PF, offset, max_pos] = next_column;
    size_t cur_size = S.size();
    assert(cur_size + offset <= size);

    auto &[S_prev, E_prev, F_prev, OS_prev, OE_prev, OF_prev,
           prev_node_prev, PS_prev, PF_prev, offset_prev, max_pos_prev]
        = column_prev[std::get<2>(prev_node)];
    assert(S_prev.size() + offset_prev <= size);

    // compute column boundaries for updating the match and deletion scores
    ssize_t offset_diff = static_cast<ssize_t>(offset_prev) - offset;

    // to define the boundaries for match scores
    // need i + offset - offset_prev - 1 >= 0
    // &&   i + offset - offset_prev - 1 < S_prev.size()
    // so offset_diff + 1 <= i < S_prev.size() + offset_diff + 1
    size_t match_begin = std::max((ssize_t)0, offset_diff + 1);
    size_t match_end = std::max(
        match_begin,
        static_cast<size_t>(std::min(S_prev.size() + offset_diff + 1, cur_size))
    );

    // to define the boundaries for deletion scores
    // need i + offset - offset_prev >= 0
    // &&   i + offset - offset_prev < S_prev.size()
    // so offset_diff <= i < S_prev.size() + offset_diff
    size_t del_begin = std::max((ssize_t)0, offset_diff);
    size_t del_end = std::max(
        del_begin,
        static_cast<size_t>(std::min(S_prev.size() + offset_diff, cur_size))
    );

    assert(del_end <= match_end);
    assert(del_end + 1 >= match_end);


    const score_t *sprev = &S_prev[offset - offset_prev];
    const score_t *fprev = &F_prev[offset - offset_prev];
    const int8_t *profile = &profile_score_.find(c)->second[start + offset];
    const Cigar::Operator *profile_o = &profile_op_.find(c)->second[start + offset];

    std::fill(PS.begin() + match_begin, PS.begin() + match_end, Extender::PREV);
    std::fill(PF.begin() + del_begin, PF.begin() + del_end, Extender::PREV);

    bool updated = false;

    auto update_match = [sprev,profile,profile_o,S=S.data(),OS=OS.data()](ssize_t i) {
        S[i + 1] = *(sprev + i) + profile[i + 1];
        OS[i + 1] = profile_o[i + 1];
    };

    auto update_del = [&config_,sprev,fprev,F=F.data(),OF=OF.data(),
                       S=S.data(),OS=OS.data()](size_t i) {
        score_t del_open = sprev[i] + config_.gap_opening_penalty;
        score_t del_extend = fprev[i] + config_.gap_extension_penalty;
        F[i] = std::max(del_open, del_extend);
        OF[i] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;

        if (F[i] > S[i]) {
            S[i] = F[i];
            OS[i] = Cigar::DELETION;
        }
    };

    if (del_begin < std::min(del_end, match_begin))
        update_del(del_begin);

    if (match_begin < match_end)
        update_match(static_cast<ssize_t>(match_begin) - 1);

#ifndef __AVX2__
    for (size_t i = match_begin; i < del_end; ++i) {
        update_del(i);
        update_match(i);
    }
#else
    static_assert(sizeof(score_t) == sizeof(int32_t));
    for (size_t i = match_begin; i < del_end; i += 8) {
        // vectorized update_match(i)
        __m256i sprev_v = _mm256_loadu_si256((__m256i*)&sprev[i]);
        __m256i profile_v = _mm256_cvtepi8_epi32(mm_loadu_si64(&profile[i + 1]));
        _mm256_storeu_si256((__m256i*)&S[i + 1], _mm256_add_epi32(sprev_v, profile_v));
        *((uint64_t*)&OS[i + 1]) = *((uint64_t*)&profile_o[i + 1]);

        // vectorized update_del(i)
        __m256i gap_open = _mm256_set1_epi32(config_.gap_opening_penalty);
        __m256i del_open = _mm256_add_epi32(sprev_v, gap_open);

        __m256i fprev_v = _mm256_loadu_si256((__m256i*)&fprev[i]);
        __m256i gap_extend = _mm256_set1_epi32(config_.gap_extension_penalty);
        __m256i del_extend = _mm256_add_epi32(fprev_v, gap_extend);

        __m256i del_score = _mm256_max_epi32(del_open, del_extend);
        _mm256_storeu_si256((__m256i*)&F[i], del_score);

        __m128i del_op_v = _mm_blendv_epi8(
            _mm_set1_epi8(Cigar::MATCH),
            _mm_set1_epi8(Cigar::DELETION),
            mm256_cvtepi32_epi8(_mm256_cmpgt_epi32(del_extend, del_open))
        );
        mm_storeu_si64(&OF[i], del_op_v);

        __m256i s_v = _mm256_loadu_si256((__m256i*)&S[i]);
        __m256i score_max = _mm256_max_epi32(s_v, del_score);
        _mm256_storeu_si256((__m256i*)&S[i], score_max);

        __m128i mask = mm256_cvtepi32_epi8(_mm256_cmpgt_epi32(del_score, s_v));
        mm_maskstorel_epi8((int8_t*)&OS[i], mask, _mm_set1_epi8(Cigar::DELETION));
    }
#endif

    auto update_max = [&xdrop_cutoff,&updated,&config_,ninf=ninf,
                       S=S.data(),E=E.data(),F=F.data(),
                       OS=OS.data(),OE=OE.data(),OF=OF.data(),
                       PS=PS.data(),PF=PF.data()](size_t i) {
        if (S[i] < xdrop_cutoff) {
            S[i] = ninf;
            E[i] = ninf;
            F[i] = ninf;
            OS[i] = Cigar::CLIPPED;
            OE[i] = Cigar::CLIPPED;
            OF[i] = Cigar::CLIPPED;
            PS[i] = Extender::NONE;
            PF[i] = Extender::NONE;
        } else if (S[i] <= 0) {
            S[i] = 0;
            OS[i] = Cigar::CLIPPED;
            PS[i] = Extender::NONE;
            updated = true;
        } else {
            xdrop_cutoff = std::max(xdrop_cutoff, S[i] - config_.xdrop);
            updated = true;
        }
    };

    update_max(0);

    size_t i = 1;
    for ( ; i < cur_size; ++i) {
        score_t ins_open = S[i - 1] + config_.gap_opening_penalty;
        score_t ins_extend = E[i - 1] + config_.gap_extension_penalty;
        E[i] = std::max(ins_open, ins_extend);
        OE[i] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;

        if (E[i] > S[i]) {
            S[i] = E[i];
            OS[i] = Cigar::INSERTION;
            PS[i] = Extender::CUR;
        }

        update_max(i);

        if (S[i] == ninf && E[i] == ninf)
            break;
    }

    for ( ; i < cur_size; ++i) {
        update_max(i);
    }

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

        PS.push_back(Extender::NONE);
        PF.push_back(Extender::NONE);
        if (S.back() > 0) {
            assert(S.back() == E.back());
            updated = true;
            PS.back() = Extender::CUR;
            OS.back() = Cigar::INSERTION;
            xdrop_cutoff = std::max(xdrop_cutoff, S.back() - config_.xdrop);
        }
    }

    // make sure that the first operation taken matches the seed
    std::ignore = graph_;
    std::ignore = seed_;
    assert(std::get<0>(prev_node) != graph_.max_index() + 1 || std::get<2>(prev_node)
        || (offset <= 1 && S[1 - offset] == seed_.get_score()
            && OS[1 - offset] == seed_.get_cigar().back().first
            && PS[1 - offset] == Extender::PREV));

    return updated;
}

template <typename NodeType, typename AlignNode, class Table, class StartSet>
void backtrack(const Table &table_,
               const Alignment<NodeType> &seed_,
               const DeBruijnGraph &graph_,
               const DBGAlignerConfig &config_,
               score_t min_path_score,
               AlignNode best_node,
               StartSet &prev_starts,
               size_t size,
               std::string_view extend_window_,
               std::string_view query,
               std::vector<Alignment<NodeType>> &extensions) {
    typedef DefaultColumnExtender<NodeType> Extender;

    Cigar cigar;
    std::vector<NodeType> path;
    std::string seq;
    NodeType start_node = DeBruijnGraph::npos;

    assert(table_.count(std::get<0>(best_node)));
    const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
        = table_.find(std::get<0>(best_node))->second.first.at(std::get<2>(best_node));

    score_t max_score = S[max_pos - offset];
    score_t score = max_score;
    Cigar::Operator last_op = OS[max_pos - offset];
    assert(last_op == Cigar::MATCH);

    if (max_pos + 1 < size)
        cigar.append(Cigar::CLIPPED, size - max_pos - 1);

    size_t pos = max_pos;
    while (true) {
        assert(table_.count(std::get<0>(best_node)));
        const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
            = table_.find(std::get<0>(best_node))->second.first.at(std::get<2>(best_node));

        prev_starts.emplace(best_node);

        assert(last_op == Cigar::MATCH || last_op == Cigar::MISMATCH);
        last_op = OS[pos - offset];

        if (last_op == Cigar::CLIPPED || S[pos - offset] == 0) {
            assert(S[pos - offset] == 0);
            max_score = score;
            break;
        } else if (pos == 1 && last_op != Cigar::DELETION) {
            score -= S[pos - offset];
            if (std::get<0>(prev) != graph_.max_index() + 1
                    || std::get<0>(best_node) != seed_.back()
                    || std::get<2>(best_node)
                    || last_op != seed_.get_cigar().back().first) {
                // last op in the seed was skipped
                // TODO: reconstruct the entire alignment. for now, throw this out
                return;
            } else {
                assert(seed_.get_score() == S[pos - offset]);
                start_node = seed_.back();
            }
            break;
        }

        switch (last_op) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                cigar.append(last_op);
                path.push_back(std::get<0>(best_node));
                seq += std::get<1>(best_node);
                assert((last_op == Cigar::MATCH)
                    == (graph_.get_node_sequence(std::get<0>(best_node)).back()
                        == extend_window_[pos - 1]));
                switch (PS[pos - offset]) {
                    case Extender::PREV: { best_node = prev; } break;
                    case Extender::CUR: {} break;
                    case Extender::NONE: { assert(false); }
                }
                --pos;
                assert(pos);
            } break;
            case Cigar::INSERTION: {
                assert(PS[pos - offset] == Extender::CUR);
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
                    assert(table_.count(std::get<0>(best_node)));
                    const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
                        = table_.find(std::get<0>(best_node))->second.first.at(std::get<2>(best_node));
                    last_op = OF[pos - offset];
                    assert(last_op == Cigar::MATCH || last_op == Cigar::DELETION);
                    path.push_back(std::get<0>(best_node));
                    seq += std::get<1>(best_node);
                    cigar.append(Cigar::DELETION);
                    switch (PF[pos - offset]) {
                        case Extender::PREV: { best_node = prev; } break;
                        case Extender::CUR: {} break;
                        case Extender::NONE: { assert(false); }
                    }
                    prev_starts.emplace(best_node);
                }
            } break;
            case Cigar::CLIPPED: { assert(false); }
        }

        assert(pos);
    }

    if (max_score < min_path_score)
        return;

    if (pos > 1)
        cigar.append(Cigar::CLIPPED, pos - 1);

    std::reverse(cigar.begin(), cigar.end());
    std::reverse(path.begin(), path.end());
    std::reverse(seq.begin(), seq.end());

    Alignment<NodeType> extension(
        { extend_window_.data() + pos, max_pos - pos },
        std::move(path), std::move(seq), score, std::move(cigar),
        0, seed_.get_orientation(), graph_.get_k() - 1
    );

    std::ignore = config_;
    assert(extension.is_valid(graph_, &config_));
    extension.extend_query_end(query.data() + query.size());

    if (start_node) {
        auto next_path = seed_;
        next_path.append(std::move(extension));
        next_path.trim_offset();
        assert(next_path.is_valid(graph_, &config_));

        DEBUG_LOG("Alignment (extended): {}", next_path);
        extensions.emplace_back(std::move(next_path));
    } else {
        extension.extend_query_begin(query.data());
        extension.trim_offset();
        assert(extension.is_valid(graph_, &config_));

        DEBUG_LOG("Alignment (trim seed): {}", extension);
        extensions.emplace_back(std::move(extension));
    }
}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>::get_extensions(score_t min_path_score)
        -> std::vector<DBGAlignment> {
    const char *align_start = seed_->get_query().data() + seed_->get_query().size() - 1;
    size_t start = align_start - query_.data();
    size_t size = query_.size() - start + 1;
    assert(start + size == partial_sums_.size());
    match_score_begin_ = partial_sums_.data() + start;

    extend_window_ = { align_start, size - 1 };
    DEBUG_LOG("Extend query window: {}", extend_window_);
    assert(extend_window_[0] == seed_->get_query().back());

    auto &first_column = table_.emplace(
        graph_.max_index() + 1,
        Column{ { std::make_tuple(
            ScoreVec(1, ninf), ScoreVec(1, ninf), ScoreVec(1, ninf),
            OpVec(1, Cigar::CLIPPED), OpVec(1, Cigar::CLIPPED), OpVec(1, Cigar::CLIPPED),
            AlignNode{}, PrevVec(1, NONE), PrevVec(1, NONE),
            0 /* offset */, 0 /* max_pos */
        ) }, false }
    ).first.value().first[0];
    sanitize(first_column);
    auto &[S, E, F, OS, OE, OF, prev_node, PS, PF, offset, max_pos] = first_column;

    size_t num_columns = 1;
    constexpr size_t column_vector_size = sizeof(std::pair<NodeType, std::pair<Column, bool>>);

    auto get_column_size = [&](const Scores &scores) {
        size_t size = std::get<0>(scores).capacity();
        return sizeof(Scores) + size * (
            sizeof(score_t) * 3 + sizeof(Cigar::Operator) * 3 + sizeof(NodeId) * 2
        );
    };
    size_t total_size = column_vector_size + get_column_size(first_column);

    S[0] = seed_->get_score() - profile_score_[seed_->get_sequence().back()][start + 1];

    AlignNode start_node{ graph_.max_index() + 1,
                          seed_->get_sequence()[seed_->get_sequence().size() - 2],
                          0, 0 };

    typedef std::pair<AlignNode, score_t> Ref;
    Ref best_start{ start_node, S[0] };
    std::vector<Ref> starts;

    std::priority_queue<Ref, std::vector<Ref>, utils::LessSecond> stack;
    stack.emplace(start_node, S[0]);

    while (stack.size()) {
        AlignNode prev = stack.top().first;
        stack.pop();

        if (static_cast<double>(total_size) / 1000000
                > config_.max_ram_per_alignment) {
            DEBUG_LOG("Alignment RAM limit reached, stopping extension");
            break;
        }

        if (static_cast<double>(num_columns) / extend_window_.size()
                > config_.max_nodes_per_seq_char) {
            DEBUG_LOG("Alignment node limit reached, stopping extension");
            break;
        }

        size_t next_distance_from_origin = std::get<3>(prev) + 1;

        for (const auto &[next, c] : get_outgoing(prev)) {
            auto &column_pair = table_[next];
            auto &[column, converged] = column_pair;
            if (converged)
                continue;

            assert(table_.count(std::get<0>(prev)));
            auto &column_prev = table_[std::get<0>(prev)].first;

            score_t xdrop_cutoff = best_start.second - config_.xdrop;

            // compute bandwidth based on xdrop criterion
            auto [min_i, max_i] = get_band(prev, column_prev, xdrop_cutoff);
            if (min_i >= max_i)
                continue;

            max_i = std::min(max_i + 1, size);

            size_t depth = column.size();
            size_t cur_size = max_i - min_i;

            Scores next_column(ScoreVec(cur_size, ninf), ScoreVec(cur_size, ninf),
                               ScoreVec(cur_size, ninf), OpVec(cur_size, Cigar::CLIPPED),
                               OpVec(cur_size, Cigar::CLIPPED),
                               OpVec(cur_size, Cigar::CLIPPED),
                               prev, PrevVec(cur_size, NONE), PrevVec(cur_size, NONE),
                               min_i /* offset */, 0 /* max_pos */);
            sanitize(next_column);

            bool updated = update_column<NodeType>(
                graph_, config_, column_prev, next_column, c, start, size,
                xdrop_cutoff, profile_score_, profile_op_, *seed_
            );
            sanitize(next_column);

            auto &[S, E, F, OS, OE, OF, prev_node, PS, PF, offset, max_pos] = next_column;

            auto max_it = std::max_element(S.begin(), S.end());
            max_pos = (max_it - S.begin()) + offset;
            assert(max_pos < size);

            converged = !updated || has_converged(column_pair, next_column);

            const score_t *match = &match_score_begin_[offset];
            bool extendable = false;
            for (size_t i = 0; i < S.size() && !extendable; ++i) {
                if (S[i] >= 0 && S[i] + match[i] >= min_path_score)
                    extendable = true;
            }

            bool add_to_table = false;
            AlignNode cur{ next, c, depth, next_distance_from_origin };
            if (OS[max_pos - offset] == Cigar::MATCH && *max_it > best_start.second) {
                best_start.first = cur;
                best_start.second = *max_it;
                add_to_table = true;
            }

            assert(xdrop_cutoff == best_start.second - config_.xdrop);

            if (*max_it >= xdrop_cutoff && extendable) {
                stack.emplace(cur, *max_it);
                add_to_table = true;
            }

            if (add_to_table) {
                total_size += get_column_size(next_column) + (!depth * column_vector_size);
                ++num_columns;
                if (OS[max_pos - offset] == Cigar::MATCH)
                    starts.emplace_back(cur, *max_it);

                column.emplace_back(std::move(next_column));
            } else if (!depth) {
                table_.erase(next);
            }
        }
    }

    std::sort(starts.begin(), starts.end(), utils::GreaterSecond());
    assert(starts.empty() || starts[0].second == best_start.second);

    struct AlignNodeHash {
        uint64_t operator()(const AlignNode &x) const {
            uint64_t seed = hasher1(std::get<0>(x));
            return seed ^ (hasher2(std::get<2>(x)) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        }

        std::hash<NodeType> hasher1;
        std::hash<size_t> hasher2;
    };

    tsl::hopscotch_set<AlignNode, AlignNodeHash> prev_starts;

    std::vector<DBGAlignment> extensions;
    for (const auto &[best_node, max_score] : starts) {
        if (prev_starts.count(best_node))
            continue;

        assert(table_.count(std::get<0>(best_node)));
        const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos]
            = table_[std::get<0>(best_node)].first.at(std::get<2>(best_node));

        assert(S[max_pos - offset] == max_score);

        if (max_pos < 2 && std::get<0>(best_node) == seed_->back()
                && !std::get<2>(best_node)) {
            if (seed_->get_score() >= min_path_score) {
                DEBUG_LOG("Alignment (seed): {}", *seed_);
                extensions.emplace_back(*seed_);
                extensions.back().extend_query_end(query_.data() + query_.size());
                extensions.back().trim_offset();
                assert(extensions.back().is_valid(graph_, &config_));
            }
        } else {
            assert(OS[max_pos - offset] == Cigar::MATCH);
            backtrack<NodeType>(table_, *seed_, graph_, config_, min_path_score, best_node,
                                prev_starts, size, extend_window_, query_, extensions);
        }

        assert(extensions.size() < 2
            || extensions.back().get_score() <= extensions[extensions.size() - 2].get_score());

        if (extensions.size() == config_.num_alternative_paths)
            break;
    }

    return extensions;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::call_visited_nodes(const std::function<void(NodeType, size_t, size_t)> &callback) const {
    size_t window_start = extend_window_.data() - query_.data();
    for (const auto &[node, columns] : table_) {
        size_t start = query_.size();
        size_t end = 0;
        size_t start_distance_from_origin = 0;
        for (const auto &column : columns.first) {
            const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos] = column;

            auto it = std::find_if(S.begin(), S.end(), [](score_t s) { return s > 0; });
            auto rit = std::find_if(S.rbegin(), S.rend(), [](score_t s) { return s > 0; });
            size_t start_c = (it - S.begin()) + offset;
            size_t end_c = (S.rend() - rit) + offset;

            size_t prev_distance_from_origin = std::get<3>(prev);

            if (start_c)
                --start_c;

            if (start_c < start) {
                start = start_c;
                start_distance_from_origin = prev_distance_from_origin
                    + (OS[it - S.begin()] != Cigar::INSERTION);
            }

            end = std::max(end, end_c);
        }

        if (start < end) {
            assert(start_distance_from_origin);
            --start_distance_from_origin;

            size_t start_pos = window_start + start + 1 - graph_.get_k();
            if (start_distance_from_origin < seed_->get_offset())
                start_pos += seed_->get_offset() - start_distance_from_origin;

            callback(node, start_pos, window_start + end);
        }
    }
}

template <typename NodeType>
bool DefaultColumnExtender<NodeType>::has_converged(const Column &column,
                                                    const Scores &next) {
    if (column.second)
        return true;

    if (column.first.empty())
        return false;

    const auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos] = next;
    const auto &[S_b, E_b, F_b, OS_b, OE_b, OF_b, prev_b, PS_b, PF_b, offset_b, max_pos_b]
        = column.first.back();

    return offset == offset_b && max_pos == max_pos_b
        && std::get<0>(prev) == std::get<0>(prev_b)
        && S == S_b && E == E_b && F == F_b && OS == OS_b && OE == OE_b && OF == OF_b
        && PS == PS_b && PF == PF_b;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::sanitize(Scores &scores) {
    auto &[S, E, F, OS, OE, OF, prev, PS, PF, offset, max_pos] = scores;

    size_t size = S.size();
    size_t pad_size = ((size + 7) / 8) * 8 + 8;
    size_t size_diff = pad_size - size;

    assert(size_diff);

    S.reserve(pad_size);
    E.reserve(pad_size);
    F.reserve(pad_size);
    OS.reserve(pad_size);
    OE.reserve(pad_size);
    OF.reserve(pad_size);
    PS.reserve(pad_size);
    PF.reserve(pad_size);

    std::fill(&S[size], &S[size] + size_diff, ninf);
    std::fill(&E[size], &E[size] + size_diff, ninf);
    std::fill(&F[size], &F[size] + size_diff, ninf);
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
