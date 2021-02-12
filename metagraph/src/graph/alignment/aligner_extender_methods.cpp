#include "aligner_extender_methods.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include <priority_deque.hpp>

#include "common/logger.hpp"
#include "common/hashers/hash.hpp"
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
    path_ = &seed;
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>::operator()(ExtensionCallback callback,
                                                 score_t min_path_score) {
    // std::cout << "START\t" << *path_ << "\n";

    typedef std::pair<NodeType, size_t> AlignNode;
    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<AlignNode> PrevVec;
    typedef AlignedVector<Cigar::Operator> OpVec;
    typedef std::pair<std::vector<std::tuple<ScoreVec, ScoreVec, ScoreVec, PrevVec, OpVec>>, bool> Column;

    tsl::hopscotch_map<NodeType, Column> table;

    const char *align_start = path_->get_query().data() + path_->get_query().size() - 1;
    size_t size = query.data() + query.size() - align_start + 1;
    extend_window_ = { align_start, size - 1 };
    assert(extend_window_[0] == path_->get_query().back());

    size_t size_per_column = sizeof(Column) + size * (
        sizeof(score_t) * 3 + sizeof(AlignNode) + sizeof(Cigar::Operator)
    );

    constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;

    Cigar::Operator last_op = path_->get_cigar().back().first;
    assert(last_op == Cigar::MATCH || last_op == Cigar::MISMATCH);

    auto &[S, E, F, P, O] = table.emplace(
        graph_.max_index() + 1,
        Column{ { std::make_tuple(ScoreVec(size, ninf), ScoreVec(size, ninf),
                                ScoreVec(size, ninf), PrevVec(size, AlignNode{}),
                                OpVec(size, Cigar::CLIPPED)) }, false }
    ).first.value().first[0];

    size_t num_columns = 1;

    S[0] = path_->get_score() - config_.get_row(extend_window_[0])[
        path_->get_sequence().back()
    ];

    AlignNode start_node{ graph_.max_index() + 1, 0 };

    typedef std::pair<AlignNode, score_t> Ref;
    boost::container::priority_deque<Ref, std::vector<Ref>, utils::LessSecond> best_starts;
    best_starts.emplace(start_node, S[0]);

    // std::vector<AlignNode> stack { best_node };
    std::priority_queue<Ref, std::vector<Ref>, utils::LessSecond> stack;
    stack.emplace(start_node, 0);

    // std::cout << "seed\t" << *path_ << "\n";
    while (stack.size()) {
        AlignNode prev = stack.top().first;
        stack.pop();

        size_t total_size = num_columns * size_per_column + table.size() * sizeof(NodeType);
        if (total_size > config_.max_ram_per_alignment * 1000000)
            break;

        if (num_columns > config_.max_nodes_per_seq_char * extend_window_.size())
            break;

        // std::cout << "foo\t" << prev.first << " " << prev.second << "\n";

        std::vector<std::pair<NodeType, char>> outgoing;
        if (prev.first == graph_.max_index() + 1) {
            outgoing.emplace_back(path_->back(), path_->get_sequence().back());
        } else {
            graph_.call_outgoing_kmers(prev.first, [&](NodeType next, char c) {
                if (c != boss::BOSS::kSentinel)
                    outgoing.emplace_back(next, c);
            });
        }

        for (const auto &[next, c] : outgoing) {
            auto &[column, converged] = table[next];
            if (converged)
                continue;

            size_t depth = column.size();
            column.emplace_back(ScoreVec(size, ninf), ScoreVec(size, ninf),
                                ScoreVec(size, ninf), PrevVec(size, AlignNode{}),
                                OpVec(size, Cigar::CLIPPED));
            ++num_columns;

            AlignNode cur{ next, depth };

            auto &[S, E, F, P, O] = column.back();
            auto &[S_prev, E_prev, F_prev, P_prev, O_prev] = table[prev.first].first[prev.second];

            for (size_t i = 1; i < S.size(); ++i) {
                E[i] = std::max(E[i - 1] + config_.gap_extension_penalty,
                                S[i - 1] + config_.gap_opening_penalty);
                F[i] = std::max(F_prev[i] + config_.gap_extension_penalty,
                                S_prev[i] + config_.gap_opening_penalty);
                score_t match_val = S_prev[i - 1] + config_.get_row(c)[extend_window_[i - 1]];

                S[i] = std::max({ E[i], F[i], match_val, 0 });

                if (S[i] > 0) {
                    if (S[i] == E[i]) {
                        P[i] = cur;
                        O[i] = Cigar::INSERTION;
                    } else if (S[i] == F[i]) {
                        P[i] = prev;
                        O[i] = Cigar::DELETION;
                    } else if (S[i] == match_val) {
                        P[i] = prev;
                        O[i] = c == extend_window_[i - 1] ? Cigar::MATCH : Cigar::MISMATCH;
                    }
                }

                assert(prev.first || prev.second || i != 1
                    || (S[1] == path_->get_score() && O[1] == path_->get_cigar().back().first
                        && next == path_->back() && !depth));

                // std::cout << i << " "
                //           << Cigar::opt_to_char(O[1]) << " "
                //           << S[1] << " "
                //           << path_->get_score() << " "
                //           << Cigar::opt_to_char(path_->get_cigar().back().first) << " "
                //           << prev.first << "," << prev.second << " "
                //           << next << "," << depth << " "
                //           << path_->back() << "\n";
                assert(i != 1 || O[1] == Cigar::DELETION || O[1] == Cigar::CLIPPED
                    || (prev.first == graph_.max_index() + 1 && !prev.second && S[1] == path_->get_score() && O[1] == path_->get_cigar().back().first
                        && next == path_->back() && !depth));

                if (best_starts.size() < config_.num_alternative_paths) {
                    best_starts.emplace(cur, S[i]);
                } else if (S[i] > best_starts.minimum().second) {
                    best_starts.update(best_starts.begin(), Ref{ cur, S[i] });
                }
            }

            if (depth) {
                auto &[S_b, E_b, F_b, P_b, O_b] = column[column.size() - 2];
                if (S == S_b && E == E_b && F == F_b && O == O_b) {
                    for (size_t i = 0; i < S.size(); ++i) {
                        if (P[i].first == P_b[i].first && P[i].second == P_b[i].second + 1) {
                            converged = true;
                        } else {
                            converged = false;
                            break;
                        }
                    }
                }
            }

            auto max_it = std::max_element(S.begin(), S.end());

            // std::cout << *max_it << "\t" << (max_it - S.begin()) << "\t"
            //           << (prev.first != graph_.max_index() + 1 ? graph_.get_node_sequence(prev.first) : std::string(graph_.get_k(), '$')) << "\t"
                      // << graph_.get_node_sequence(cur.first) << "\t";
                      // << cur.first << "," << cur.second << "\t"
                      // << prev.first << "," << prev.second << "\t";
            // for (size_t i = 0; i < S.size(); ++i) {
            //     std::cout << std::max(S[i],-9) << "," << std::max(E[i],-9) << "," << std::max(F[i],-9) << "\t";
            // }
            // std::cout << "\n";

            score_t score_rest = *max_it + config_.match_score(
                extend_window_.substr(max_it - S.begin())
            );

            score_t xdrop_cutoff = best_starts.maximum().second - config_.xdrop;

            // std::cout << "bar\t" << *max_it << " " << (max_it - S.begin()) << "\n";
            if (*max_it >= xdrop_cutoff && score_rest >= min_path_score) {
                // std::cout << "bar\t" << cur.first << "," << cur.second << "\n";
                stack.emplace(Ref{ cur, *max_it });
            }
        };
    }

    while (best_starts.size()) {
        auto [best_node, max_score] = best_starts.maximum();
        best_starts.pop_maximum();

        const auto &best_scores = std::get<0>(table[best_node.first].first[best_node.second]);
        size_t max_pos = std::max_element(best_scores.begin(), best_scores.end())
            - best_scores.begin();

        assert(best_scores[max_pos] == max_score);

        if (max_pos < 2 && best_node.first == path_->back() && !best_node.second) {
            callback(Alignment<NodeType>(), 0);
            return;
        }

        Cigar cigar;
        if (max_pos + 1 < F.size())
            cigar.append(Cigar::CLIPPED, F.size() - max_pos - 1);

        size_t pos = max_pos;
        std::vector<NodeType> path;
        NodeType start_node = 0;
        score_t score = max_score;

        while (true) {
            auto &[S, E, F, P, O] = table[best_node.first].first[best_node.second];

            // std::cout << pos << "\t" << Cigar::opt_to_char(O[pos]) << "\t" << best_node.first << "," << best_node.second << "\t";
            // for (size_t i = 0; i < S.size(); ++i) {
            //     std::cout << std::max(S[i],-9) << "," << std::max(E[i],-9) << "," << std::max(F[i],-9) << "\t";
            // }
            // std::cout << "\n";

            // std::cout << "p\t"
            //           << pos << "," << max_pos << "\t"
            //           << best_node.first << "," << best_node.second << "\t"
            //           << P[pos].first << "," << P[pos].second << "\t"
            //           << S[pos] << "," << E[pos] << "," << F[pos] << "\t"
            //           << cigar.to_string() << " "
            //           << Cigar::opt_to_char(O[pos]) << "\t"
            //           << path_->back() << "\n";

            if (pos == 1 && best_node.first == path_->back() && !best_node.second && O[pos] == path_->get_cigar().back().first) {
                assert(P[pos].first == graph_.max_index() + 1);
                assert(!P[pos].second);
                start_node = path_->back();
                score -= path_->get_score();
                // std::cout << "a\t"
                //           << pos << "," << max_pos << "\t"
                //           << best_node.first << "," << best_node.second << "\t"
                //           << P[pos].first << "," << P[pos].second << "\t"
                //           << S[pos] << "\t"
                //           << cigar.to_string() << " "
                //           << Cigar::opt_to_char(O[pos]) << "\t"
                //           << path_->back() << "\n";
                break;
            }

            if (O[pos] == Cigar::CLIPPED || (pos == 1 && O[pos] != Cigar::DELETION)) {
                start_node = DeBruijnGraph::npos;
                // std::cout << "b\t"
                //           << pos << "," << max_pos << "\t"
                //           << best_node.first << "," << best_node.second << "\t"
                //           << P[pos].first << "," << P[pos].second << "\t"
                //           << S[pos] << "\t"
                //           << cigar.to_string() << " "
                //           << Cigar::opt_to_char(O[pos]) << "\t"
                //           << path_->back() << "\n";
                break;
            }

            switch (O[pos]) {
                case Cigar::MATCH:
                case Cigar::MISMATCH: {
                    cigar.append(O[pos]);
                    path.push_back(best_node.first);
                    assert((O[pos] == Cigar::MATCH)
                        == (graph_.get_node_sequence(best_node.first).back() == extend_window_[pos - 1]));
                    best_node = P[pos];
                    --pos;
                } break;
                case Cigar::INSERTION: {
                    assert(P[pos] == best_node);
                    cigar.append(O[pos]);
                    --pos;
                } break;
                case Cigar::DELETION: {
                path.push_back(best_node.first);
                    cigar.append(O[pos]);
                    best_node = P[pos];
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
        std::string seq;
        spell_path(graph_, path, seq, graph_.get_k() - 1);


        Alignment<NodeType> extension({ extend_window_.data() + pos, max_pos - pos },
                                      std::move(path),
                                      std::move(seq),
                                      score,
                                      std::move(cigar),
                                      0,
                                      path_->get_orientation(),
                                      graph_.get_k() - 1);
        // std::cout << *path_ << " " << extension << "\n";

        assert(extension.is_valid(graph_, &config_));
        callback(std::move(extension), start_node);
    }
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
