#include "aligner_extender_methods.hpp"

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
DefaultColumnExtender<NodeType>::DefaultColumnExtender(const DeBruijnGraph &graph,
                                                       const DBGAlignerConfig &config,
                                                       std::string_view query)
      : graph_(graph), config_(config), query(query) {
    assert(config_.check_config_scores());
    partial_sums_.resize(query.size());
    std::transform(query.begin(), query.end(),
                   partial_sums_.begin(),
                   [&](char c) { return config_.get_row(c)[c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(config_.match_score(query) == partial_sums_.front());
    assert(config_.get_row(query.back())[query.back()] == partial_sums_.back());

    for (char c : graph_.alphabet()) {
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

template <typename NodeType>
void DefaultColumnExtender<NodeType>::initialize(const DBGAlignment &seed) {
    // this extender only works if at least one character has been matched
    assert(seed.get_query().size());
    assert(query.data() <= seed.get_query().data());

    const char *query_end = query.data() + query.size();
    const char *seed_end = seed.get_query().data() + seed.get_query().size();
    assert(query_end >= seed_end);

    // size includes the last character of the seed since the upper-left corner
    // of the score matrix is the seed score
    size = query_end - seed_end + 1;

    // extend_window_ doesn't include this last character
    extend_window_ = { seed_end - 1, size };

    // match_score_begin = partial_sums_.data() + (seed_end - 1 - query.data());
    // assert(config_.get_row(query.back())[query.back()] == match_score_begin[size - 1]);
    // assert(config_.match_score(std::string_view(extend_window_.data() - 1, size))
    //     == *match_score_begin);

    // start_node = seed.back();
    this->path_ = &seed;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::operator()(ExtensionCallback callback,
                                                 score_t min_path_score) {
    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<NodeType> PrevVec;
    typedef AlignedVector<Cigar::Operator> OpVec;
    typedef std::tuple<ScoreVec, ScoreVec, ScoreVec, PrevVec, OpVec> Column;
    tsl::hopscotch_map<NodeType, Column> table;

    constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;

    auto emplace = table.emplace(path_->back(), std::make_tuple(
        ScoreVec(size, 0), ScoreVec(size, ninf), ScoreVec(size, ninf),
        PrevVec(size, DeBruijnGraph::npos), OpVec(size, Cigar::CLIPPED)
    ));

    Cigar::Operator last_op = path_->get_cigar().back().first;

    auto &[S, E, F, P, O] = emplace.first.value();
    S[0] = path_->get_score();

    score_t gap_score = path_->get_score() - config_.get_row(path_->get_query().back())[path_->get_sequence().back()]
                            + 2 * config_.gap_opening_penalty;
    switch (last_op) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            E[0] = gap_score;
            F[0] = gap_score;
        } break;
        case Cigar::INSERTION: {
            E[0] = path_->get_score();
            F[0] = gap_score;
        } break;
        case Cigar::DELETION: {
            E[0] = gap_score;
            F[0] = path_->get_score();
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    for (size_t i = 1; i < size; ++i) {
        E[i] = E[i - 1] + config_.gap_extension_penalty;
        if (E[i] > S[i]) {
            S[i] = E[i];
            P[i] = path_->back();
            O[i] = Cigar::INSERTION;
        }
    }

    if (last_op == Cigar::DELETION) {
        F[1] = path_->get_score() + config_.gap_opening_penalty;
        if (F[1] > S[1]) {
            S[1] = F[1];
            O[1] = Cigar::DELETION;
        }
        for (size_t i = 2; i < size; ++i) {
            F[i] = F[i - 1] + config_.gap_extension_penalty;
            if (F[i] > S[i]) {
                S[i] = F[i];
                O[i] = Cigar::DELETION;
            }
        }
    } else if (last_op != Cigar::INSERTION) {
        for (size_t i = 1; i < size; ++i) {
            F[i] = gap_score + config_.gap_extension_penalty * i;
            if (F[i] > S[i]) {
                S[i] = F[i];
                O[i] = Cigar::DELETION;
            }
        }
    }

    typedef std::pair<score_t, NodeType> ColumnRef;
    std::priority_queue<ColumnRef, std::vector<ColumnRef>, utils::LessFirst> queue;
    queue.emplace(0, path_->back());

    score_t max_score = std::max(min_path_score, path_->get_score());
    size_t max_pos = 0;
    NodeType best_node = 0;

    while (queue.size()) {
        ColumnRef top = queue.top();
        queue.pop();

        NodeType prev = top.second;
        assert(prev);

        bool pushed = false;

        graph_.call_outgoing_kmers(prev, [&](NodeType next, char c) {
            if (c == boss::BOSS::kSentinel)
                return;

            auto find = table.find(next);
            if (find == table.end()) {
                find = table.emplace(next, std::make_tuple(
                    ScoreVec(size, 0), ScoreVec(size, ninf), ScoreVec(size, ninf),
                    PrevVec(size, 0), OpVec(size, Cigar::CLIPPED)
                )).first;
            }

            auto &[S, E, F, P, O] = find.value();
            auto prev_find = table.find(prev);
            assert(prev_find != table.end());
            auto &[S_prev, E_prev, F_prev, P_prev, O_prev] = prev_find.value();

            score_t max_diff = 0;
            size_t max_pos_local = 0;

            for (size_t i = 0; i < F.size(); ++i) {
                score_t new_val = std::max({ F[i],
                                             F_prev[i] + config_.gap_extension_penalty,
                                             S_prev[i] + config_.gap_opening_penalty });
                F[i] = new_val;
            }

            if (F[0] > S[0]) {
                max_diff = std::max(max_diff, F[0] - S[0]);
                S[0] = F[0];
                P[0] = prev;
                O[0] = Cigar::DELETION;
            }

            score_t max_val = S[0];

            for (size_t i = 1; i < S.size(); ++i) {
                score_t ins_val = std::max({ E[i],
                                             E[i - 1] + config_.gap_extension_penalty,
                                             S[i - 1] + config_.gap_opening_penalty });
                E[i] = ins_val;

                score_t match_val = S_prev[i - 1] + config_.get_row(c)[extend_window_[i]];

                score_t best_val = std::max({ S[i], E[i], F[i], match_val });
                if (best_val == S[i]) {
                } else if (best_val == match_val) {
                    P[i] = prev;
                    O[i] = c == extend_window_[i] ? Cigar::MATCH : Cigar::MISMATCH;
                } else if (best_val == ins_val) {
                    P[i] = next;
                    O[i] = Cigar::INSERTION;
                } else {
                    P[i] = prev;
                    O[i] = Cigar::DELETION;
                }

                max_diff = std::max(max_diff, best_val - S[i]);
                S[i] = best_val;

                if (best_val > max_val) {
                    max_val = best_val;
                    max_pos_local = i;
                }
            }

            if (max_val > max_score) {
                max_score = max_val;
                max_pos = max_pos_local;
                best_node = next;
            }

            std::string_view check_window = extend_window_;
            check_window.remove_prefix(max_pos + 1);
            // if (max_diff > 0 && max_val >= max_score - config_.xdrop && S[max_pos] + config_.match_score(check_window) >= max_score) {
            if (max_diff > 0) {
                queue.emplace(max_diff, next);
                pushed = true;
            }
        });

        if (pushed)
            queue.emplace(std::move(top));
    }

    if (!max_pos)
        return;

    assert(std::get<0>(table[best_node])[max_pos] == max_score);
    // std::cout << "max score: " << max_score << "\n";

    // backtrack
    Cigar cigar;
    if (max_pos + 1 < F.size())
        cigar.append(Cigar::CLIPPED, F.size() - max_pos - 1);

    size_t pos = max_pos;
    std::vector<NodeType> path;

#ifndef NDEBUG
    score_t score_track = max_score;
#endif

    while (best_node) {
        auto &[S, E, F, P, O] = table[best_node];

        Cigar::Operator cur_op = O[pos];

#ifndef NDEBUG
        // std::cout << score_track << "\t" << pos << " " << (uint64_t)O[pos] << " " << best_node << " " << P[pos] << " " << S[pos] << "\t" << path.size() << "\t" << cigar.to_string() << "\n";

        if (cigar.size() && cigar.back().first != O[pos] && (cigar.back().first == Cigar::INSERTION || cigar.back().first == Cigar::DELETION)) {
            if (score_track == S[pos]) {
                // indel was replaced

                // score_track -= config_.gap_opening_penalty - config_.gap_extension_penalty;
                // if (cigar.back().first == Cigar::INSERTION) {
                //     if (score_track != E[pos])
                //         std::cout << "fail\t" << score_track << " " << E[pos] << "\n";
                //     assert(score_track == E[pos]);
                // } else {
                //     if (score_track != F[pos])
                //         std::cout << "fail\t" << score_track << " " << F[pos] << "\n";
                //     assert(score_track == F[pos]);
                // }
            } else {
                score_track -= config_.gap_opening_penalty - config_.gap_extension_penalty;
                assert(score_track == S[pos]);
            }

        }
#endif

        switch (cur_op) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                cigar.append(cur_op);
                path.push_back(best_node);
#ifndef NDEBUG
                char last_char = graph_.get_node_sequence(best_node).back();
                assert((last_char == extend_window_[pos]) == (O[pos] == Cigar::MATCH));
                score_track -= config_.get_row(extend_window_[pos])[last_char];
#endif
                best_node = P[pos];
                --pos;
            } break;
            case Cigar::INSERTION: {
                assert(P[pos] == best_node);
                cigar.append(cur_op);
#ifndef NDEBUG
                score_track -= config_.gap_extension_penalty;
#endif
                --pos;
            } break;
            case Cigar::DELETION: {
                path.push_back(best_node);
                cigar.append(O[pos]);
                best_node = P[pos];
#ifndef NDEBUG
                score_track -= config_.gap_extension_penalty;
#endif
            } break;
            case Cigar::CLIPPED: {
                best_node = DeBruijnGraph::npos;
                assert(!pos || S[pos] == 0);
            }
        }
    }

    assert(path.size());

    if (pos)
        cigar.append(Cigar::CLIPPED, pos);

    std::reverse(cigar.begin(), cigar.end());

    std::reverse(path.begin(), path.end());
    NodeType start_node = pos ? path.front() : path_->back();
    std::string seq;
    spell_path(graph_, path, seq, graph_.get_k() - 1);

    Alignment<NodeType> extension({ extend_window_.data() + pos + 1, max_pos - pos },
                                  std::move(path),
                                  std::move(seq),
                                  max_score - (pos ? 0 : path_->get_score()),
                                  std::move(cigar),
                                  0,
                                  path_->get_orientation(),
                                  graph_.get_k() - 1);

    assert(extension.is_valid(graph_, &config_));
    callback(std::move(extension), start_node);
}


template <typename NodeType>
auto DefaultColumnExtender<NodeType>::emplace_node(NodeType, NodeType, char, size_t, size_t, size_t, size_t, size_t)
-> std::pair<typename DPTable<NodeType>::iterator, bool> { return {}; }

template <typename NodeType>
bool DefaultColumnExtender<NodeType>::add_seed(size_t) { return true; }

template <typename NodeType>
void DefaultColumnExtender<NodeType>::check_and_push(ColumnRef&&) {}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::extend_main(ExtensionCallback, score_t) {}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::update_columns(NodeType,
                    const std::deque<std::pair<NodeType, char>> &,
                    score_t) {}

template <typename NodeType>
auto DefaultColumnExtender<NodeType>::fork_extension(NodeType, ExtensionCallback, score_t)
-> std::deque<std::pair<NodeType, char>> { return {}; }


template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
