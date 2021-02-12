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
        auto &profile_score_row = profile_score_.emplace(c, query_.size() + 8).first.value();
        auto &profile_op_row = profile_op_.emplace(c, query_.size() + 8).first.value();

        const auto &row = config_.get_row(c);
        const auto &op_row = Cigar::get_op_row(c);

        std::transform(query_.begin(), query_.end(), profile_score_row.begin(),
                       [&row](char q) { return row[q]; });
        std::transform(query_.begin(), query_.end(), profile_op_row.begin(),
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

    size_t size_per_column = sizeof(Column) + size * (
        sizeof(score_t) * 3 + sizeof(AlignNode) * 2 + sizeof(Cigar::Operator) * 3
    );

    constexpr score_t ninf = std::numeric_limits<score_t>::min() + 100;

    assert(seed_->get_cigar().back().first == Cigar::MATCH
        || seed_->get_cigar().back().first == Cigar::MISMATCH);

    auto &[S, E, F, OS, OE, OF, PS, PF] = table_.emplace(
        graph_.max_index() + 1,
        Column{ { std::make_tuple(
            ScoreVec(size, ninf), ScoreVec(size, ninf), ScoreVec(size, ninf),
            OpVec(size, Cigar::CLIPPED), OpVec(size, Cigar::CLIPPED), OpVec(size, Cigar::CLIPPED),
            PrevVec(size, AlignNode{}), PrevVec(size, AlignNode{})
        )}, false }
    ).first.value().first[0];

    size_t num_columns = 1;

    S[0] = seed_->get_score() - profile_score_[seed_->get_sequence().back()][start];

    AlignNode start_node{ graph_.max_index() + 1, 0 };

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

        std::vector<std::pair<NodeType, char>> outgoing;
        if (prev.first == graph_.max_index() + 1) {
            outgoing.emplace_back(seed_->back(), seed_->get_sequence().back());
        } else {
            graph_.call_outgoing_kmers(prev.first, [&](NodeType next, char c) {
                if (c != boss::BOSS::kSentinel)
                    outgoing.emplace_back(next, c);
            });
        }

        for (const auto &[next, c] : outgoing) {
            auto &[column, converged] = table_[next];
            if (converged)
                continue;

            size_t depth = column.size();
            column.emplace_back(
                ScoreVec(size, ninf), ScoreVec(size, ninf), ScoreVec(size, ninf),
                OpVec(size, Cigar::CLIPPED), OpVec(size, Cigar::CLIPPED), OpVec(size, Cigar::CLIPPED),
                PrevVec(size, AlignNode{}), PrevVec(size, AlignNode{})
            );
            ++num_columns;

            AlignNode cur{ next, depth };

            auto &[S, E, F, OS, OE, OF, PS, PF] = column.back();
            auto &[S_prev, E_prev, F_prev, OS_prev, OE_prev, OF_prev, PS_prev, PF_prev]
                = table_[prev.first].first[prev.second];

            for (size_t i = 1; i < S.size(); ++i) {
                score_t ins_open   = S[i - 1] + config_.gap_opening_penalty;
                score_t ins_extend = E[i - 1] + config_.gap_extension_penalty;
                E[i] = std::max(ins_open, ins_extend);
                OE[i] = ins_open < ins_extend ? Cigar::INSERTION : Cigar::MATCH;

                score_t del_open   = S_prev[i] + config_.gap_opening_penalty;
                score_t del_extend = F_prev[i] + config_.gap_extension_penalty;
                PF[i] = prev;
                F[i] = std::max(del_open, del_extend);
                OF[i] = del_open < del_extend ? Cigar::DELETION : Cigar::MATCH;

                score_t match = S_prev[i - 1] + profile_score_[c][start + i - 1];

                S[i] = std::max({ E[i], F[i], match, 0 });

                if (S[i] > 0) {
                    if (S[i] == E[i]) {
                        PS[i] = cur;
                        OS[i] = Cigar::INSERTION;
                    } else if (S[i] == F[i]) {
                        PS[i] = prev;
                        OS[i] = Cigar::DELETION;
                    } else if (S[i] == match) {
                        PS[i] = prev;
                        OS[i] = profile_op_[c][start + i - 1];
                    }
                }

                // sanity check to ensure that the first operation taken is the
                // last operation from the seed
                assert(i != 1 || OS[1] == Cigar::DELETION || OS[1] == Cigar::CLIPPED
                    || (prev.first == graph_.max_index() + 1
                        && !prev.second && S[1] == seed_->get_score()
                        && OS[1] == seed_->get_cigar().back().first
                        && next == seed_->back() && !depth));
            }

            // TODO
            // if (depth) {
            //     auto &[S_b, E_b, F_b, P_b, O_b] = column[column.size() - 2];
            //     if (S == S_b && E == E_b && F == F_b && O == O_b) {
            //         for (size_t i = 0; i < S.size(); ++i) {
            //             if (P[i].first == P_b[i].first && (P[i].second == P_b[i].second + 1
            //                     || (!P[i].first && !P[i].second && !P_b[i].second))) {
            //                 converged = true;
            //             } else {
            //                 converged = false;
            //                 break;
            //             }
            //         }
            //     }
            // }

            auto max_it = std::max_element(S.begin(), S.end());
            if (best_starts.size() < config_.num_alternative_paths) {
                best_starts.emplace(cur, *max_it);
            } else if (*max_it > best_starts.minimum().second) {
                best_starts.update(best_starts.begin(), Ref{ cur, *max_it });
            }

            assert(match_score_begin_[max_it - S.begin()]
                == config_.match_score(extend_window_.substr(max_it - S.begin())));
            score_t score_rest = *max_it + match_score_begin_[max_it - S.begin()];

            score_t xdrop_cutoff = best_starts.maximum().second - config_.xdrop;

            if (*max_it >= xdrop_cutoff && score_rest >= min_seed_score)
                stack.emplace(Ref{ cur, *max_it });
        };
    }

    while (best_starts.size()) {
        auto [best_node, max_score] = best_starts.maximum();
        best_starts.pop_maximum();

        const auto &best_scores
            = std::get<0>(table_[best_node.first].first[best_node.second]);
        size_t max_pos = std::max_element(best_scores.begin(), best_scores.end())
            - best_scores.begin();

        assert(best_scores[max_pos] == max_score);

        if (max_pos < 2 && best_node.first == seed_->back() && !best_node.second) {
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

        Cigar::Operator last_op
            = std::get<4>(table_[best_node.first].first[best_node.second])[pos];
        assert(last_op == Cigar::MATCH);

        while (true) {
            const auto &[S, E, F, OS, OE, OF, PS, PF]
                = table_[best_node.first].first[best_node.second];

            const auto &O = (last_op == Cigar::MATCH || last_op == Cigar::MISMATCH)
                ? OS
                : (last_op == Cigar::INSERTION ? OE : OF);

            if (last_op == Cigar::MATCH || last_op == Cigar::MISMATCH) {
                if (pos == 1 && best_node.first == seed_->back()
                        && !best_node.second
                        && O[pos] == seed_->get_cigar().back().first) {
                    assert(P[pos].first == graph_.max_index() + 1);
                    assert(!P[pos].second);
                    start_node = seed_->back();
                    score -= seed_->get_score();
                    break;
                } else if (O[pos] == Cigar::CLIPPED) {
                    start_node = DeBruijnGraph::npos;
                    break;
                }
            }

            switch (last_op) {
                case Cigar::MATCH:
                case Cigar::MISMATCH: {
                    last_op = O[pos];
                    switch (O[pos]) {
                        case Cigar::MATCH:
                        case Cigar::MISMATCH: {
                            cigar.append(O[pos]);
                            path.push_back(best_node.first);
                            assert((O[pos] == Cigar::MATCH)
                                == (graph_.get_node_sequence(best_node.first).back()
                                    == extend_window_[pos - 1]));
                            best_node = PS[pos];
                            --pos;
                        } break;
                        case Cigar::INSERTION: { assert(PS[pos] == best_node); } break;
                        case Cigar::DELETION: {} break;
                        case Cigar::CLIPPED: { assert(false); }
                    }
                } break;
                case Cigar::INSERTION: {
                    last_op = O[pos];
                    assert(last_op == Cigar::MATCH || last_op == Cigar::INSERTION);
                    cigar.append(Cigar::INSERTION);
                    --pos;
                } break;
                case Cigar::DELETION: {
                    last_op = O[pos];
                    assert(last_op == Cigar::MATCH || last_op == Cigar::DELETION);
                    path.push_back(best_node.first);
                    cigar.append(Cigar::DELETION);
                    best_node = PF[pos];
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
        std::string seq;
        spell_path(graph_, path, seq, graph_.get_k() - 1);

        Alignment<NodeType> extension({ extend_window_.data() + pos, max_pos - pos },
                                      std::move(path),
                                      std::move(seq),
                                      score,
                                      std::move(cigar),
                                      0,
                                      seed_->get_orientation(),
                                      graph_.get_k() - 1);

        if (!extension.is_valid(graph_, &config_)) {
            std::cerr << "Seed: " << *seed_;
            if (start_node)
                std::cerr << " " << start_node;

            std::cerr << std::endl;
            throw std::runtime_error("Internal error: invalid extension");
        }

        callback(std::move(extension), start_node);
    }
}

template <typename NodeType>
void DefaultColumnExtender<NodeType>
::call_visited_nodes(const std::function<void(NodeType, size_t, size_t)> &callback) const {
    size_t offset = extend_window_.data() - query_.data();
    for (const auto &[node, columns] : table_) {
        size_t start = query_.size();
        size_t end = 0;
        for (const auto &column : columns.first) {
            const auto &[S, E, F, OS, OE, OF, PS, PF] = column;
            assert(S.size() == extend_window_.size() + 1);

            auto it = std::find_if(S.begin(), S.end(), [](score_t s) { return s > 0; });
            start = std::min(start, static_cast<size_t>(it - S.begin()) - 1);

            end = std::max(
                end, static_cast<size_t>(std::max_element(it, S.end()) - S.begin()) - 1
            );
        }
        callback(node, offset + start, offset + end);
    }
}

template class DefaultColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
