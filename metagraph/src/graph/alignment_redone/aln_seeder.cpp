#include "aln_seeder.hpp"

#include <tsl/hopscotch_set.h>

#include "aln_chainer.hpp"
#include "common/logger.hpp"
#include "common/vector_map.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace mtg::graph::align_redone {

std::vector<Anchor> ExactSeeder::get_anchors() const {
    const DeBruijnGraph &graph = query_.get_graph();
    std::vector<Anchor> anchors;

    if (query_.get_query().size() < config_.min_seed_length)
        return anchors;

    if (config_.min_seed_length >= graph.get_k()) {
        for (bool orientation : { false, true }) {
            std::string_view this_query = query_.get_query(orientation);
            size_t begin = 0;
            graph.map_to_nodes_sequentially(this_query,
                [&](Match::node_index node) {
                    if (node != DeBruijnGraph::npos) {
                        anchors.emplace_back(this_query,
                                             begin, begin + graph.get_k(),
                                             orientation,
                                             std::vector<Match::node_index>{ node },
                                             config_);
                    }

                    ++begin;
                }
            );
        }
    } else {
        const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
        if (!dbg_succ)
            return anchors;

        const auto &boss = dbg_succ->get_boss();
        for (bool orientation : { false, true }) {
            std::string_view this_query = query_.get_query(orientation);
            auto encoded = boss.encode(this_query);

            size_t k = graph.get_k();
            for (size_t i = 0; i + k <= this_query.size(); ++i) {
                if (std::find(encoded.begin() + i, encoded.end() + config_.min_seed_length,
                            boss.alph_size) != encoded.end() + config_.min_seed_length) {
                    continue;
                }
                auto [first, last, it] = boss.index_range(
                    encoded.begin() + i,
                    encoded.begin() + i + config_.min_seed_length
                );

                size_t match_size = it - (encoded.begin() + i);
                if (match_size < config_.min_seed_length)
                    continue;

                using edge_index = boss::BOSS::edge_index;
                std::vector<std::tuple<edge_index, edge_index, std::string>> traverse;
                traverse.emplace_back(first, last, "");
                while (traverse.size()) {
                    auto [first, last, suffix] = traverse.back();
                    assert(boss.get_last(last));
                    traverse.pop_back();

                    assert(suffix.size() + match_size < k);
                    if (suffix.size() + match_size == k - 1) {
                        assert(first == last || !boss.get_last(first));
                        assert(boss.succ_last(first) == last);
                        for (edge_index node = first; node <= last; ++node) {
                            if (dbg_succ->in_graph(node)) {
                                char c = boss.decode(boss.get_W(node) % boss.alph_size);
                                if (c != boss::BOSS::kSentinel) {
                                    anchors.emplace_back(this_query,
                                                        i, i + match_size,
                                                        orientation,
                                                        std::vector<Match::node_index>{ node },
                                                        config_,
                                                        suffix + c);
                                }
                            }
                        }

                        continue;
                    }

                    for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
                        boss::BOSS::edge_index next_first = first;
                        boss::BOSS::edge_index next_last = last;
                        if (boss.tighten_range(&next_first, &next_last, s)) {
                            traverse.emplace_back(next_first, next_last, suffix + boss.decode(s));
                        }
                    }
                }
            }
        }
    }

    // if (anchors.size() && (config_.max_seed_length > config_.min_seed_length || config_.min_seed_length > graph.get_k())) {
    //     // common::logger->info("Merging {} anchors", anchors.size());
    //     auto rbegin = anchors.rbegin();
    //     auto rend = anchors.rend();
    //     for (auto it = rbegin; it + 1 != rend; ++it) {
    //         auto &a_j = *it;
    //         auto &a_i = *(it + 1);
    //         if (a_i.empty() || a_j.empty() || a_i.get_orientation() != a_j.get_orientation()
    //                 || a_i.get_label_class() != a_j.get_label_class() || a_i.get_end_trim())
    //             continue;

    //         std::string_view a_i_s = a_i.get_spelling();
    //         std::string_view a_j_s = a_j.get_spelling();
    //         size_t overlap = graph.get_k() - 1;
    //         if (a_i_s.size() >= overlap
    //                 && a_j_s.size() >= overlap
    //                 && a_i.get_clipping() + a_i_s.size() - a_j.get_clipping() >= overlap
    //                 && a_j.get_clipping() + a_j_s.size() - a_i.get_clipping() <= config_.max_seed_length
    //                 && (a_i.get_seed().size() < config_.min_seed_length
    //                         || (graph.has_single_incoming(a_j.get_path()[0])
    //                         && graph.has_single_outgoing(a_i.get_path().back())))
    //                 && std::equal(a_i_s.end() - overlap, a_i_s.end(), a_j_s.begin(), a_j_s.begin() + overlap)) {
    //             // merge them
    //             a_i.append(a_j, config_, &query_.get_graph());
    //             a_j = std::decay_t<decltype(a_j)>();
    //         }
    //     }

    //     auto end = std::remove_if(anchors.begin(), anchors.end(), [&](auto &a) {
    //         return a.empty() || a.get_seed().size() < config_.min_seed_length;
    //     });

    //     anchors.erase(end, anchors.end());
    // }

    return anchors;
}

template <typename T>
class OffsetVector {
  public:
    OffsetVector() : offset_(0) {}

    template <typename... Args>
    OffsetVector(size_t offset, Args&&... args) : offset_(offset), null_(T()), data_(std::forward<Args>(args)...) {}

    size_t size() const { return offset_ + data_.size(); }
    bool empty() const { return offset_ == 0 && data_.empty(); }

    void resize(size_t new_size, const T &val = T()) {
        if (new_size >= offset_) {
            data_.resize(new_size - offset_, val);
            return;
        }

        new_size -= data_.size();
        data_.clear();
        assert(new_size <= offset_);
        offset_ -= new_size;
    }

    size_t offset() const { return offset_; }

    template <typename... Args>
    T& set(size_t i, Args&&... args) {
        get(i) = T(std::forward<Args>(args)...);
        return data_[i - offset_];
    }

    T& get(size_t i) {
        if (i < offset_) {
            data_.insert(data_.begin(), offset_ - i, null_);
            offset_ = i;
        }

        if (i - offset_ >= data_.size())
            data_.resize(i - offset_ + 1);

        return data_[i - offset_];
    }

    T& operator[](size_t i) {
        assert(i >= offset_);
        return data_[i - offset_];
    }

    const T& operator[](size_t i) const {
        return i >= offset_ ? data_[i - offset_] : null_;
    }

    using const_iterator = typename std::vector<T>::const_iterator;
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend() const { return data_.cend(); }

  private:
    size_t offset_;
    const T null_;
    std::vector<T> data_;
};

template <typename Tuple>
using WaveFront = OffsetVector<VectorMap<DeBruijnGraph::node_index, Tuple>>; // S(cost)

template <typename Tuple>
using ScoreTable = std::vector<WaveFront<Tuple>>;

DBGAlignerConfig::score_t cost_to_score(size_t cost,
                                        size_t query_size,
                                        size_t match_size,
                                        DBGAlignerConfig::score_t match_score) {
    DBGAlignerConfig::score_t double_score = match_score * (query_size + match_size) - cost;
    assert(double_score % 2 == 0);
    return double_score / 2;
}

template <typename Tuple>
const Tuple* get_bucket(const ScoreTable<Tuple> &table,
                        size_t cost,
                        size_t query_dist,
                        DeBruijnGraph::node_index node) {
    if (cost >= table.size())
        return nullptr;

    if (query_dist >= table[cost].size() || query_dist < table[cost].offset())
        return nullptr;

    const auto &table_slot = table[cost][query_dist];
    auto it = table_slot.find(node);
    if (it == table_slot.end())
        return nullptr;

    return &it->second;
};

// dist, num_ops, last_node, last_op, last_char_of_this_node, num_matches
using SMap = std::tuple<size_t, size_t, DeBruijnGraph::node_index, Cigar::Operator, char, size_t>;

// dist, num_ops
using EMap = std::tuple<size_t, size_t>;

// dist, num_ops, last_node
using FMap = std::tuple<size_t, size_t, DeBruijnGraph::node_index>;

template <typename StrItr>
void align_impl(const std::function<size_t(DeBruijnGraph::node_index, size_t, size_t)> &has_single_incoming,
                const std::function<DeBruijnGraph::node_index(DeBruijnGraph::node_index, char, size_t, size_t)> &traverse_back,
                const std::function<void(DeBruijnGraph::node_index, const std::function<void(DeBruijnGraph::node_index, char)>&, size_t, size_t)> &call_incoming_kmers,
                const DBGAlignerConfig &config,
                const DeBruijnGraph::node_index start_node,
                size_t start_num_matches,
                const StrItr query_window_begin,
                const StrItr query_window_end,
                size_t max_dist,
                const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t cost)> &callback,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate) {
    using node_index = DeBruijnGraph::node_index;

    if (query_window_begin == query_window_end)
        return;

    size_t query_size = query_window_end - query_window_begin;

    const auto query_window_rbegin = std::make_reverse_iterator(query_window_end);
    const auto query_window_rend = std::make_reverse_iterator(query_window_begin);

    // derived from 2.4.1 in https://doi.org/10.1101/2022.01.12.476087 (v1)
    DBGAlignerConfig::score_t match_score = config.match_score("A");
    DBGAlignerConfig::score_t mismatch_score = config.score_sequences("A", "T");
    assert(config.gap_opening_penalty <= config.gap_extension_penalty);
    assert(config.gap_extension_penalty < 0);
    assert(mismatch_score < match_score);

    ssize_t mismatch_cost = (match_score - mismatch_score) * 2;
    ssize_t gap_ext = match_score - 2 * config.gap_extension_penalty;
    ssize_t gap_opn = match_score - 2 * config.gap_opening_penalty;

    assert(mismatch_cost > 0);
    assert(gap_ext > 0);
    assert(gap_opn >= gap_ext);

    assert(mismatch_cost % 2 == 0);
    assert(gap_ext % 2 == 0);
    assert(gap_opn % 2 == 0);

    common::logger->info("x: {}\to: {}\te: {}", mismatch_cost, gap_opn, gap_ext);

    // S(cost, query_dist, node) = best_dist
    ScoreTable<SMap> S; // best cost
    ScoreTable<EMap> E; // best cost (last operation is insertion)
    ScoreTable<FMap> F; // best cost (last operation is deletion)

    auto set_value = [&](
            auto &table,
            size_t cost,
            size_t query_dist,
            node_index node,
            size_t dist,
            node_index last_node,
            size_t num_ops,
            Cigar::Operator last_op,
            char c,
            size_t num_matches) -> size_t {
        assert(last_op == Cigar::MATCH
            || last_op == Cigar::MISMATCH
            || last_op == Cigar::INSERTION
            || last_op == Cigar::DELETION);
        if (cost >= table.size())
            table.resize(cost + 1);

        if (query_dist >= table[cost].size())
            table[cost].resize(query_dist + 1);

        auto &bucket = table[cost].get(query_dist)[node];
        using table_t = std::decay_t<decltype(bucket)>;
        static_assert(std::is_same_v<table_t, FMap>
                        || std::is_same_v<table_t, EMap>
                        || std::is_same_v<table_t, SMap>);

        bool inserted = dist > std::get<0>(bucket);

        if constexpr(std::is_same_v<table_t, SMap>) {
            inserted |= std::get<3>(bucket) == Cigar::CLIPPED || std::make_pair(dist, num_matches) > std::make_pair(std::get<0>(bucket), std::get<5>(bucket));
        }

        if constexpr(std::is_same_v<table_t, FMap>)
            inserted |= std::get<2>(bucket) == DeBruijnGraph::npos;

        if constexpr(std::is_same_v<table_t, EMap>)
            inserted |= std::get<1>(bucket) == 0;

        if (inserted) {
            std::get<0>(bucket) = dist;
            std::get<1>(bucket) = num_ops;

            if constexpr(std::is_same_v<table_t, FMap>
                            || std::is_same_v<table_t, SMap>) {
                assert(dist >= num_ops);
                std::get<2>(bucket) = last_node;
            }

            if constexpr(std::is_same_v<table_t, SMap>) {
                std::get<3>(bucket) = last_op;
                std::get<4>(bucket) = c;
                std::get<5>(bucket) = num_matches;
            }

            inserted = true;
        }

        if constexpr(std::is_same_v<table_t, SMap>) {
            assert(std::get<3>(bucket) == Cigar::MATCH
                || std::get<3>(bucket) == Cigar::MISMATCH
                || std::get<3>(bucket) == Cigar::INSERTION
                || std::get<3>(bucket) == Cigar::DELETION);
        }

        return inserted;
    };

    size_t num_explored_nodes = 0;
    {
        size_t query_dist = 0;
        DeBruijnGraph::node_index node = start_node;
        SMap data(0, 0, start_node, Cigar::MATCH, *(query_window_rbegin - 1), start_num_matches);
        auto &[best_dist, last_ext, last_node, last_op, c, num_matches] = data;

        for (auto it = query_window_rbegin; it != query_window_rend; ++it) {
            if (terminate_branch(0, data, query_dist, node) || best_dist == max_dist || !has_single_incoming(node, best_dist, query_dist))
                break;

            if (auto prev = traverse_back(node, *it, best_dist, query_dist)) {
                ++num_explored_nodes;
                ++best_dist;
                ++query_dist;
                ++last_ext;
                ++num_matches;
                node = prev;
            } else {
                break;
            }
        }

        set_value(
            S,
            0,
            query_dist,
            node,
            best_dist,
            start_node,
            best_dist,
            Cigar::MATCH,
            best_dist > 0 ? *query_window_rbegin : '\0',
            num_matches
        );
    }

    // size_t max_edit_cost = std::max({ gap_opn, gap_ext, mismatch_cost });

    // for (size_t cost = 0; cost < S.size() + max_edit_cost; ++cost) {
    for (size_t cost = 0; cost < S.size(); ++cost) {
        bool done = false;
        for (size_t query_dist = S[cost].offset(); query_dist < S[cost].size(); ++query_dist) {
            assert(query_dist <= query_size);

            size_t it_dist = 0;
            for (auto it = S[cost][query_dist].begin(); it != S[cost][query_dist].end(); ++it, ++it_dist) {
                node_index node = it->first;
                size_t best_dist = std::get<0>(it->second);
                Cigar::Operator last_op = std::get<3>(it->second);
                size_t num_matches = std::get<5>(it->second);
                assert(last_op != Cigar::CLIPPED);

                std::cerr << "Exp: " << cost << "\t" << node << "\t" << best_dist << "," << query_dist << "\n";

                if (terminate(cost, it->second, query_dist, node)) {
                    done = true;
                    // continue;
                    // common::logger->info("Halted at cost {}", cost);
                    break;
                }

                if (terminate_branch(cost, it->second, query_dist, node))
                    continue;

                // forward creation of insertions
                if (query_dist < query_size) {
                    // extend a previous insertion
                    if (const auto *ins_ext_bucket = get_bucket(E, cost, query_dist, node)) {
                        auto [last_dist, last_num_ops] = *ins_ext_bucket;
                        assert(last_num_ops > 0);
                        assert(query_dist >= last_num_ops);
                        ssize_t next_ext_cost = cost + gap_ext;
                        if (!terminate_branch(next_ext_cost, SMap(last_dist, 0, node, Cigar::INSERTION, '\0', 0), query_dist + 1, node)
                                && set_value(E, next_ext_cost, query_dist + 1, node, last_dist, node, last_num_ops + 1, Cigar::INSERTION, '\0', 0)
                                && set_value(S, next_ext_cost, query_dist + 1, node, last_dist, node, 0, Cigar::INSERTION, '\0', 0)) {
                            it = S[cost][query_dist].begin() + it_dist;
                        }
                    }

                    // open an insertion
                    if (last_op != Cigar::DELETION) {
                        ssize_t next_opn_cost = cost + gap_opn;
                        if (!terminate_branch(next_opn_cost, SMap(best_dist, 0, node, Cigar::INSERTION, '\0', 0), query_dist + 1, node)
                                && set_value(E, next_opn_cost, query_dist + 1, node, best_dist, node, 1, Cigar::INSERTION, '\0', 0)
                                && set_value(S, next_opn_cost, query_dist + 1, node, best_dist, node, 0, Cigar::INSERTION, '\0', 0)) {
                            it = S[cost][query_dist].begin() + it_dist;
                        }
                    }
                }

                if (best_dist < max_dist) {
                    std::vector<std::pair<node_index, char>> prevs;
                    call_incoming_kmers(node, [&](auto prev, char c) {
                        ++num_explored_nodes;
                        if (c != boss::BOSS::kSentinel)
                            prevs.emplace_back(prev, c);
                    }, best_dist, query_dist);

                    if (prevs.empty())
                        continue;

                    // deletion extension
                    if (const auto *del_ext_bucket = get_bucket(F, cost, query_dist, node)) {
                        auto [last_dist, last_num_ops, last_node] = *del_ext_bucket;
                        assert(last_node != DeBruijnGraph::npos);
                        assert(last_num_ops > 0);
                        assert(last_dist >= last_num_ops);
                        if (last_dist < max_dist) {
                            ssize_t next_ext_cost = cost + gap_ext;
                            std::vector<std::pair<DeBruijnGraph::node_index, char>> last_prevs;
                            if (best_dist == last_dist) {
                                last_prevs = prevs;
                            } else {
                                call_incoming_kmers(node, [&](auto prev, char c) {
                                    ++num_explored_nodes;
                                    if (c != boss::BOSS::kSentinel)
                                        last_prevs.emplace_back(prev, c);
                                }, last_dist, query_dist);
                            }

                            for (const auto &[prev, c] : last_prevs) {
                                if (!terminate_branch(next_ext_cost, SMap(last_dist + 1, 0, node, Cigar::DELETION, '\0', 0), query_dist, prev)
                                        && set_value(F, next_ext_cost, query_dist, prev, last_dist + 1, node, last_num_ops + 1, Cigar::DELETION, '\0', 0)
                                        && set_value(S, next_ext_cost, query_dist, prev, last_dist + 1, node, 0, Cigar::DELETION, '\0', 0)) {
                                    it = S[cost][query_dist].begin() + it_dist;
                                }
                            }
                        }
                    }

                    // deletion open
                    if (last_op != Cigar::INSERTION) {
                        ssize_t next_opn_cost = cost + gap_opn;
                        for (const auto &[prev, c] : prevs) {
                            if (!terminate_branch(next_opn_cost, SMap(best_dist + 1, 0, node, Cigar::DELETION, '\0', 0), query_dist, prev)
                                    && set_value(F, next_opn_cost, query_dist, prev, best_dist + 1, node, 1, Cigar::DELETION, '\0', 0)
                                    && set_value(S, next_opn_cost, query_dist, prev, best_dist + 1, node, 0, Cigar::DELETION, '\0', 0)) {
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }
                    }

                    if (query_dist < query_size) {
                        // match
                        auto local_query_window_begin = query_window_begin;
                        auto local_query_window_end = query_window_end;
                        local_query_window_end -= query_dist;

                        auto local_query_window_rbegin = std::make_reverse_iterator(local_query_window_end);
                        auto local_query_window_rend = std::make_reverse_iterator(local_query_window_begin);

                        for (auto [prev, c] : prevs) {
                            SMap data(best_dist + 1, 1, node, Cigar::MATCH, c, num_matches + 1);
                            auto &[cur_best, num_ops, last_node, op, cur_c, cur_num_matches] = data;
                            size_t cur_query_dist = query_dist + 1;
                            size_t cur_cost = cost;
                            if (c != *local_query_window_rbegin) {
                                cur_cost += mismatch_cost;
                                op = Cigar::MISMATCH;
                                cur_num_matches = 0;
                            }

                            for (auto jt = local_query_window_rbegin + 1; jt != local_query_window_rend; ++jt) {
                                if (terminate_branch(cur_cost, data, cur_query_dist, prev) || cur_best == max_dist || !has_single_incoming(prev, cur_best, cur_query_dist))
                                    break;

                                if (auto pprev = traverse_back(prev, *jt, cur_best, cur_query_dist)) {
                                    ++num_explored_nodes;
                                    ++cur_best;
                                    ++cur_query_dist;
                                    ++cur_num_matches;
                                    prev = pprev;
                                } else {
                                    break;
                                }
                            }

                            if (set_value(
                                S,
                                cur_cost,
                                cur_query_dist,
                                prev,
                                cur_best,
                                node,
                                cur_query_dist - query_dist,
                                op,
                                c,
                                cur_num_matches
                            )) {
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }
                    }
                }
            }
            // if (done)
            //     break;
        }

        if (done)
            break;
    }

    // backtrack
    for (size_t start_cost = 0; start_cost < S.size(); ++start_cost) {
        std::vector<std::tuple<size_t, ssize_t, size_t, DeBruijnGraph::node_index>> starts;
        for (size_t start_query_dist = S[start_cost].offset(); start_query_dist < S[start_cost].size(); ++start_query_dist) {
            assert(start_query_dist <= query_size);
            const auto &bucket = S[start_cost][start_query_dist];
            for (const auto &[node, data] : bucket) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                assert(last_op != Cigar::CLIPPED);
                // std::cerr << "BTIter: " << start_cost << "\t" << node << "\t" << dist << "," << start_query_dist << std::endl;
                if (start_backtrack(start_cost, data, start_query_dist, node))
                    starts.emplace_back(dist, static_cast<ssize_t>(dist) - start_query_dist, start_query_dist, node);
            }
        }

        for (auto [dist, diag, query_dist, node] : starts) {
            size_t cost = start_cost;

            const auto *bt = get_bucket(S, cost, query_dist, node);
            assert(bt);

            std::vector<node_index> path;
            Cigar cigar;

            do {
                assert(std::get<0>(*bt) == dist);

                if (std::get<3>(*bt) == Cigar::INSERTION) {
                    assert(std::get<1>(*bt) == 0);
                    assert(cost > 0);
                    const auto *et = get_bucket(E, cost, query_dist, node);
                    assert(et);

                    auto [ins_dist, num_ins] = *et;
                    assert(query_dist >= num_ins);
                    assert(num_ins > 0);
                    assert(ins_dist == dist);

                    size_t prev_cost = cost - gap_opn - (num_ins - 1) * gap_ext;
                    assert(prev_cost < cost);
                    size_t prev_query_dist = query_dist - num_ins;

                    const auto *check_bt = get_bucket(S, prev_cost, prev_query_dist, node);
                    assert(check_bt);
                    assert(std::get<0>(*check_bt) == dist);
                    assert(std::get<3>(*check_bt) != Cigar::CLIPPED);

                    query_dist = prev_query_dist;
                    cost = prev_cost;
                    bt = check_bt;
                    cigar.append(Cigar::INSERTION, num_ins);
                } else if (std::get<3>(*bt) == Cigar::DELETION) {
                    assert(std::get<1>(*bt) == 0);
                    assert(cost > 0);
                    const auto *ft = get_bucket(F, cost, query_dist, node);
                    assert(ft);

                    auto [del_dist, num_del, prev_node] = *ft;
                    assert(prev_node != DeBruijnGraph::npos);
                    assert(del_dist == dist);
                    assert(dist >= num_del);

                    cigar.append(Cigar::DELETION, num_del);

                    size_t prev_cost = cost;
                    size_t prev_dist = dist;
                    while (num_del > 1) {
                        path.emplace_back(node);

                        prev_cost -= gap_ext;
                        --prev_dist;
                        node = prev_node;
                        --num_del;

                        const auto *prev_ft = get_bucket(F, prev_cost, query_dist, node);
                        assert(prev_ft);
                        assert(std::get<0>(*prev_ft) == prev_dist);
                        assert(std::get<1>(*prev_ft) == num_del);
                        ft = prev_ft;
                        prev_node = std::get<2>(*ft);
                    }

                    path.emplace_back(node);

                    prev_cost -= gap_opn;
                    --prev_dist;

                    assert(prev_cost < cost);

                    const auto *check_bt = get_bucket(S, prev_cost, query_dist, prev_node);
                    assert(check_bt);
                    assert(std::get<0>(*check_bt) == prev_dist);
                    assert(std::get<3>(*check_bt) != Cigar::CLIPPED);

                    dist = prev_dist;
                    cost = prev_cost;
                    node = prev_node;
                    bt = check_bt;
                } else if (std::get<3>(*bt) == Cigar::MATCH || std::get<3>(*bt) == Cigar::MISMATCH) {
                    auto [cur_dist, num_match, last_node, last_op, mismatch_char, num_matches] = *bt;

                    if (!num_match) {
                        assert(last_op == Cigar::MATCH);
                        bt = nullptr;
                        break;
                    }

                    assert(query_dist >= num_match);
                    assert(dist >= num_match);
                    assert(cur_dist == dist);

                    node_index traverse_node = last_node;

                    cigar.append(Cigar::MATCH, num_match - 1);
                    cigar.append(last_op);
                    size_t prev_dist = dist - num_match;
                    size_t prev_query = query_dist - num_match;

                    // reconstruct the backwards traversal from last_node to node
                    auto it = query_window_rbegin + prev_query;
                    assert(it != query_window_rend);
                    size_t cur_path_size = path.size();
                    size_t num_traversed = 0;
                    if (std::get<3>(*bt) == Cigar::MISMATCH) {
                        assert(mismatch_char != *it);
                        traverse_node = traverse_back(traverse_node, mismatch_char, prev_dist, prev_query);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        ++it;
                        --num_match;
                        ++num_traversed;
                    } else {
                        assert(mismatch_char == *it);
                    }

                    while (num_match) {
                        assert(it != query_window_rend);
                        traverse_node = traverse_back(traverse_node, *it, prev_dist + num_traversed, prev_query + num_traversed);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        ++it;
                        --num_match;
                        ++num_traversed;
                    }
                    assert(traverse_node == node);

                    std::reverse(path.begin() + cur_path_size, path.end());

                    size_t prev_cost = cost - (last_op == Cigar::MATCH ? 0 : mismatch_cost);
                    assert(prev_cost <= cost);

                    const auto *check_bt = bt;
                    check_bt = nullptr;
                    if (prev_cost || prev_query || prev_dist) {
                        check_bt = get_bucket(S, prev_cost, prev_query, last_node);
                        assert(check_bt);
                        assert(std::get<0>(*check_bt) == prev_dist);
                        assert(std::get<3>(*check_bt) != Cigar::CLIPPED);
                    } else {
                        assert(last_node == start_node);
                    }

                    dist = prev_dist;
                    cost = prev_cost;
                    query_dist = prev_query;
                    node = last_node;

                    if (!check_bt)
                        break;

                    bt = check_bt;
                } else {
                    assert(false && "Different operation found");
                    throw std::runtime_error("Inf. loop");
                }

            } while (bt);

            assert(dist == 0);
            assert(cost == 0);
            assert(query_dist == 0);
            std::cerr << "Ext: " << start_cost << "\t" << cigar.to_string() << "\n";
            callback(std::move(path), std::move(cigar), start_cost);
        }
    }
    common::logger->info("Explored {} nodes", num_explored_nodes);
}

void align_bwd(const DeBruijnGraph &graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate) {
    align_impl(
        [&graph](DeBruijnGraph::node_index node, size_t, size_t) { return graph.has_single_incoming(node); },
        [&graph](DeBruijnGraph::node_index node, char c, size_t, size_t) { return graph.traverse_back(node, c); },
        [&graph](DeBruijnGraph::node_index node, const auto &callback, size_t, size_t) {
            graph.call_incoming_kmers(node, callback);
        },
        config,
        start_node,
        num_matches,
        query_window.begin(), query_window.end(),
        max_dist,
        callback,
        start_backtrack,
        terminate_branch,
        terminate
    );
}

void align_fwd(const DeBruijnGraph &graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate,
               std::string_view suffix = "") {
    align_impl(
        [&graph,&suffix,&start_node](DeBruijnGraph::node_index node, size_t dist, size_t) {
            std::ignore = start_node;
            assert(dist >= suffix.size() || node == start_node);
            return dist < suffix.size() || graph.has_single_outgoing(node);
        },
        [&graph,&suffix,&start_node](DeBruijnGraph::node_index node, char c, size_t dist, size_t) {
            if (dist < suffix.size()) {
                std::ignore = start_node;
                assert(node == start_node);
                return c == suffix[dist] ? node : DeBruijnGraph::npos;
            } else {
                return graph.traverse(node, c);
            }
        },
        [&graph,&suffix,&start_node](DeBruijnGraph::node_index node, const auto &callback, size_t dist, size_t) {
            if (dist < suffix.size()) {
                std::ignore = start_node;
                assert(node == start_node);
                callback(node, suffix[dist]);
            } else {
                graph.call_outgoing_kmers(node, callback);
            }
        },
        config,
        start_node,
        num_matches,
        query_window.rbegin(), query_window.rend(),
        max_dist,
        [&](auto&& path, auto&& cigar, size_t cost) {
            if (path.size() >= suffix.size()) {
                path.resize(path.size() - suffix.size());
                std::reverse(path.begin(), path.end());
                std::reverse(cigar.data().begin(), cigar.data().end());
                callback(std::move(path), std::move(cigar), cost);
            }
        },
        start_backtrack,
        terminate_branch,
        terminate
    );
}

void Extender::extend(const Alignment &aln, const std::function<void(Alignment&&)> &callback, bool no_bwd, bool no_fwd) const {
    std::cerr << "Base " << aln << "\n";
    DBGAlignerConfig::score_t match_score = config_.match_score("A");

    std::vector<Alignment> fwd_exts;
    if (no_fwd || !aln.get_end_clipping()) {
        if (!aln.get_end_trim()) {
            std::cerr << "\tskipfwd\t" << aln << "\n";
            fwd_exts.emplace_back(aln);
        }
    } else {
        std::string_view query_window = aln.get_query();
        query_window.remove_prefix(aln.get_clipping() + aln.get_seed().size());

        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return aln.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.right_end_bonus : 0);
        };

        auto it = std::make_reverse_iterator(aln.get_cigar().data().end());
        auto end = std::make_reverse_iterator(aln.get_cigar().data().begin());
        std::ignore = end;
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;
        size_t max_dist = 0;
        DBGAlignerConfig::score_t best_score = aln.get_score();
        size_t max_query_dist = 0;
        align_fwd(
            query_.get_graph(), config_, aln.get_path().back(), num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& path, auto&& cigar, size_t) {
                auto ext_path = aln.get_path();
                ext_path.insert(ext_path.end(), path.begin(), path.end());

                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                cigar.append(Cigar::CLIPPED, query_window.size() - query_dist);

                Cigar ext_cigar = aln.get_cigar();
                ext_cigar.trim_end_clipping();
                ext_cigar.append(std::move(cigar));

                fwd_exts.emplace_back(
                    query_.get_graph(),
                    aln.get_query(),
                    aln.get_orientation(),
                    std::move(ext_path),
                    config_,
                    std::move(ext_cigar)
                );

                // assert(fwd_exts.back().get_score() == best_score);
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                // size_t num_matches = std::get<5>(data);
                // if (!num_matches || dist < aln.get_end_trim() || (dist == 0 && query_dist == 0))
                // if (!num_matches || (dist == 0 && query_dist == 0))
                if (dist == 0 && query_dist == 0)
                    return false;

                if (query_dist < max_query_dist)
                    return false;

                auto score = get_score(cost, dist, query_dist);
                // std::cerr << "BT?: " << score << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                assert(score <= best_score);
                return score == best_score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate branch
                if (query_dist == query_window.size())
                    return true;

                size_t dist = std::get<0>(data);
                auto score = get_score(cost, dist, query_dist);
                return score <= 0 || config_.xdrop < best_score - score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate
                size_t dist = std::get<0>(data);
                max_dist = std::max(dist, max_dist);
                auto score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    max_query_dist = std::max(max_query_dist, query_dist);
                }
                // if (query_dist == query_window.size())
                //     std::cerr << "Exppp: " << cost << "\n";
                // return query_dist == query_window.size() && dist == max_dist;
                return query_dist == query_window.size();
            },
            aln.get_trim_spelling()
        );

        if (fwd_exts.empty() && !aln.get_end_trim())
            fwd_exts.emplace_back(aln);
    }

    for (auto &fwd_ext : fwd_exts) {
        assert(!fwd_ext.get_end_trim());
        if (no_bwd || !fwd_ext.get_clipping()) {
            std::cerr << "\tskipbwd\t" << fwd_ext << "\n";
            callback(std::move(fwd_ext));
            continue;
        }

        std::string_view query_window = fwd_ext.get_query();
        query_window.remove_suffix(fwd_ext.get_end_clipping() + fwd_ext.get_seed().size());

        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return fwd_ext.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.left_end_bonus : 0);
        };

        DBGAlignerConfig::score_t best_score = fwd_ext.get_score();

        auto it = aln.get_cigar().data().begin();
        auto end = aln.get_cigar().data().end();
        std::ignore = end;
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;

        bool found = false;
        size_t max_dist = 0;
        size_t max_query_dist = 0;
        // std::cerr << "bwd extending\t" << fwd_ext << "\n";
        align_bwd(
            query_.get_graph(), config_, fwd_ext.get_path()[0], num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& path, auto&& cigar, size_t) {
                path.insert(path.end(), fwd_ext.get_path().begin(), fwd_ext.get_path().end());
                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                Cigar ext_cigar(Cigar::CLIPPED, query_window.size() - query_dist);
                ext_cigar.append(std::move(cigar));

                Cigar aln_cigar = fwd_ext.get_cigar();
                aln_cigar.trim_clipping();
                ext_cigar.append(std::move(aln_cigar));

                Alignment next_aln(query_.get_graph(),
                    fwd_ext.get_query(),
                    fwd_ext.get_orientation(),
                    std::move(path),
                    config_,
                    std::move(ext_cigar)
                );

                // assert(next_aln.get_score() == best_score);
                // std::cerr << "found\t" << next_aln << "\n";
                callback(std::move(next_aln));
                found = true;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                // size_t num_matches = std::get<5>(data);
                // if ((dist == 0 && query_dist == 0) || !num_matches)
                //     return false;

                if (dist == 0 && query_dist == 0)
                    return false;

                if (query_dist < max_query_dist)
                    return false;

                auto score = get_score(cost, dist, query_dist);
                assert(score <= best_score);
                return score == best_score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate branch
                if (query_dist == query_window.size())
                    return true;

                size_t dist = std::get<0>(data);
                auto score = get_score(cost, dist, query_dist);
                return score <= 0 || config_.xdrop < best_score - score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate
                size_t dist = std::get<0>(data);
                max_dist = std::max(dist, max_dist);
                auto score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    max_query_dist = std::max(max_query_dist, query_dist);
                }
                // if (query_dist == query_window.size())
                //     std::cerr << "Exppp: " << cost << "\n";
                // return false;
                // return query_dist == query_window.size() && dist == max_dist;
                return query_dist == query_window.size();
            }
        );

        if (!found)
            callback(std::move(fwd_ext));
    }
}

std::vector<Alignment> ExactSeeder::get_inexact_anchors() const {
    std::vector<Anchor> anchors = get_anchors();

    std::vector<Alignment> alignments;

    if (anchors.empty())
        return alignments;

    using AnchorIt = std::vector<Anchor>::iterator;
    std::sort(anchors.begin(), anchors.end(), AnchorLess<Anchor>());

    std::cerr << "Anchors\n";
    for (const auto &a : anchors) {
        std::cerr << "\t" << a << "\t" << a.get_path_spelling() << "\n";
    }

    const auto &graph = query_.get_graph();
    ssize_t k = graph.get_k();

    sdsl::bit_vector selected(anchors.size(), false);

    bool first_chain = true;
    DBGAlignerConfig::score_t first_chain_score = DBGAlignerConfig::ninf;
    size_t num_chains = 0;

    std::string dummy(query_.get_query().size(), '$');

    chain_anchors<AnchorIt>(query_, config_, anchors.begin(), anchors.end(),
        [this,k,&graph,&dummy](
                const Anchor &a_j,
                ssize_t max_dist,
                AnchorIt rbegin,
                AnchorIt rend,
                typename ChainScores<AnchorIt>::iterator rchain_scores,
                const ScoreUpdater<AnchorIt> &update_score) {
            auto rchain_scores_end = rchain_scores + (rend - rbegin);

            auto begin = std::make_reverse_iterator(rend);
            auto end = std::make_reverse_iterator(rbegin);
            auto chain_scores = std::make_reverse_iterator(rchain_scores_end);

            std::string_view query_j = a_j.get_seed();
            // const DBGAlignerConfig::score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));
            const DBGAlignerConfig::score_t &score_j = std::get<0>(*(rchain_scores_end));

            for (auto rit = begin; rit != end; ++rit, ++chain_scores) {
                auto it = rit.base() - 1;
                const Anchor &a_i = *it;
                assert(a_i.get_label_class() == a_j.get_label_class());
                assert(a_i.get_orientation() == a_j.get_orientation());

                std::string_view query_i = a_i.get_seed();

                DBGAlignerConfig::score_t dist = query_j.end() - query_i.end();

                if (dist <= 0 || query_i.begin() >= query_j.begin())
                    continue;

                DBGAlignerConfig::score_t base_score = std::get<0>(*chain_scores);

                auto added_seq_start = std::max(query_j.begin(), query_i.end());
                std::string_view added_seq(added_seq_start, query_j.end() - added_seq_start);

                DBGAlignerConfig::score_t score = base_score + config_.match_score(added_seq);

                // std::cerr << "chh\t" << a_i << " -> " << a_j << "\t" << score << " vs. " << score_j << "\n";
                if (score <= score_j)
                    continue;

                size_t i_seed_size = k - (a_i.get_path_spelling().size() - query_i.size());
                auto last_node_i_query = query_i.end() - i_seed_size;

                if (last_node_i_query >= query_j.begin()) {
                    // there should be overlapping nodes
                    // std::cerr << "olN" << a_i << " -> " << a_j << "\n";
                    size_t j_start = last_node_i_query - query_j.begin();
                    if (a_i.get_path().size() >= j_start
                            && a_j.get_path().size() >= j_start
                            && std::equal(a_j.get_path().begin(), a_j.get_path().begin() + j_start,
                                          a_i.get_path().end() - j_start, a_i.get_path().end())) {
                        // std::cerr << "\t\tworked! " << score << "\n";
                        update_score(score, it, dist);
                    }
                    continue;
                } else if (query_i.end() > query_j.begin()) {
                    // overlapping nucleotides
                    // std::cerr << "oln" << a_i << " -> " << a_j << "\n";
                    const std::string spelling = a_j.get_path_spelling();
                    const char *jt_end = spelling.c_str() + k;
                    const char *jt = jt_end - (query_j.begin() - last_node_i_query);
                    ssize_t num_matches = 0;
                    DeBruijnGraph::node_index last_node = DeBruijnGraph::npos;
                    graph.traverse(a_i.get_path().back(), jt, jt_end, [&](DeBruijnGraph::node_index cur) {
                        ++num_matches;
                        last_node = cur;
                    });
                    // assert(num_matches < jt_end - jt || last_node == a_j.get_path()[0]);
                    // std::cerr << "\t" << std::string_view(jt, jt_end - jt) << "\t" << num_matches << " vs. " << jt_end - jt << "\n";
                    if (num_matches == jt_end - jt && last_node == a_j.get_path()[0]) {
                        // std::cerr << "\t\tworked! " << score << "\n";
                        update_score(score, it, dist);
                    }

                    // std::cerr << "ol" << a_i << " -> " << a_j << "\n";
                    // if (last_node_i_query + 1 >= query_j.begin()) {
                    //     size_t num_chars_in_j = last_node_i_query - query_j.begin() + 1;
                    //     if (num_chars_in_j < a_j.get_path().size()
                    //             && graph.traverse(a_i.get_path().back(), *query_i.end()) == a_j.get_path()[num_chars_in_j]) {
                    //         // traversal worked
                    //         // std::cerr << "\tworked\n";
                    //         update_score(score, it, dist);
                    //     }
                    // // } else {
                    // //     size_t gap = query_j.begin() - (last_node_i_query + 1);
                    // //     // std::cerr << "dl" << a_i << " -> " << a_j << "\t" << gap << "\n";
                    // //     score += config_.gap_opening_penalty + (gap - 1) * config_.gap_extension_penalty;
                    // //     if (score > score_j)
                    // //         update_score(score, it, dist);
                    // }
                    continue;
                // } else if (dist <= k) {
                    // first, check for insertion
                    // if (graph.traverse(a_i.get_path().back(), a_j.get_path_spelling()[k - 1]) == a_j.get_path()[0]) {
                    //     size_t i_seed_size = k - (a_i.get_path_spelling().size() - query_i.size());
                    //     auto last_node_i_query = query_i.end() - i_seed_size;
                    //     size_t ins_length = query_j.begin() - last_node_i_query - 1;
                    //     DBGAlignerConfig::score_t ins_score = score + config_.gap_opening_penalty + (ins_length - 1) * config_.gap_extension_penalty;
                    //     if (ins_score > score_j)
                    //         update_score(ins_score, it, dist);
                    // }
                    // // gap, so we need to count indels
                    // size_t gap = query_j.begin() - query_i.end() + 1;
                    // score += config_.gap_opening_penalty + (gap - 1) * config_.gap_extension_penalty;
                    // if (score > score_j)
                    //     update_score(score, it, dist);
                    // continue;
                } else if (query_j.begin() == query_i.end()) {
                    size_t gap = k;
                    score += config_.gap_opening_penalty + (gap - 1) * config_.gap_extension_penalty;
                    // std::cerr << "gap\t" << a_i << " -> " << a_j << "\t" << gap << "\t" << score << "\n";
                    if (score > score_j) {
                        // std::cerr << "\t\tworked! " << score << "\n";
                        update_score(score, it, dist);
                    }
                    continue;
                }

                // TODO: deal with indels?
                size_t nmismatch = query_j.begin() - query_i.end();
                DBGAlignerConfig::score_t cur_mismatch = config_.score_sequences(std::string_view(&*query_i.end(), nmismatch),
                                                           std::string_view(dummy.c_str(), nmismatch));
                // std::cerr << "uk" << a_i << " -> " << a_j << "\t" << nmismatch << "\t" << cur_mismatch << "\n";
                score += cur_mismatch;
                if (score > score_j) {
                    // std::cerr << "\t\tworked! " << score << "\n";
                    update_score(score, it, dist);
                }
            }
        },
        [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
            assert(chain.size());

            bool ret_val = (!chain[0].first->get_clipping() && !chain[0].first->get_end_clipping())
                            || chain.size() > 1
                            || std::any_of(chain.begin(), chain.end(),
                                           [](const auto &a) { return a.first->get_path().size() > 1; });

            ret_val &= std::all_of(chain.begin(), chain.end(),
                                   [&](const auto &a) { return !selected[a.first - anchors.begin()]; });

            if (ret_val) {
                common::logger->info("Chain\t{}\t{}", score_traceback.back(), ret_val);
                for (const auto &[it, dist] : chain) {
                    std::cerr << "\t" << *it << "\t" << dist << "\n";
                }
            }

            if (!ret_val)
                return false;

            if (first_chain) {
                first_chain_score = score_traceback.back();
                assert(score_traceback.front() <= score_traceback.back());
                first_chain = false;
            }

            ret_val = score_traceback.back() >= first_chain_score * config_.rel_score_cutoff;

            if (ret_val) {
                ++num_chains;
                for (const auto &[it, dist] : chain) {
                    selected[it - anchors.begin()] = true;
                }
            }

            return ret_val;
        },
        [this,match_score=config_.match_score("A")](
                AnchorIt last,
               AnchorIt next,
               Alignment&& aln,
               size_t last_to_next_dist,
               DBGAlignerConfig::score_t score_up_to_now,
               const AlignmentCallback &callback) {
            if (next->get_seed().begin() + 1 == last->get_seed().begin()) {
                // simple extension
                std::cerr << "Simple extension: " << *next << "\t->\t" << aln << std::endl;
                // const auto &graph = query_.get_graph();
                // size_t k = graph.get_k();
                // if (graph.traverse(next->get_path()[0], aln.get_path_spelling()[k-1]) != aln.get_path()[0]) {
                //     std::cerr << "\t" << next->get_path()[0] << " > " << aln.get_path_spelling()[k-1] << "\t" << graph.traverse(next->get_path()[0], aln.get_path_spelling()[k-1]) << " vs. " << aln.get_path()[0] << std::endl;
                //     throw std::runtime_error("FJLSKFDJF");
                // }
                // common::logger->info("Simple extension {} -> {}", *next, aln);
                std::vector<DeBruijnGraph::node_index> path;
                path.resize(aln.get_path().size() + 1);
                path[0] = next->get_path()[0];
                std::copy(aln.get_path().begin(), aln.get_path().end(), path.begin() + 1);

                Cigar cigar(Cigar::CLIPPED, next->get_clipping());

                Cigar aln_cigar = aln.get_cigar();
                aln_cigar.trim_clipping();

                cigar.append(Cigar::MATCH, 1);
                cigar.append(std::move(aln_cigar));

                Alignment next_aln(query_.get_graph(),
                                aln.get_query(),
                                aln.get_orientation(),
                                std::move(path),
                                config_,
                                std::move(cigar),
                                aln.get_end_trim());
                std::cerr << "\tnew aln: " << next_aln << std::endl;
                callback(std::move(next_aln));
                return;
            }

            // common::logger->info("aln extension {} -> {}", *next, aln);
            std::cerr << "Alignment extension: " << *next << "\t->\t" << aln << std::endl;

            std::string_view query_window(
                next->get_seed().data(),
                aln.get_seed().begin() - next->get_seed().begin()
            );
            auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
                return aln.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                        + (query_dist == query_window.size() && !next->get_clipping() ? config_.left_end_bonus : 0);
            };

            bool aln_found = false;

            auto start_backtracking = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                // common::logger->info("Checkbt: n: {}, d: {}, qd: {} / {}, s: {}, nm: {}",
                //                      graph.get_node_sequence(node), dist, query_dist, query_window.size(),
                //                      get_score(cost, query_dist, dist), num_matches);

                if (aln_found || dist == 0 || query_dist == 0 || num_matches < next->get_seed().size())
                    return false;

                return query_dist == query_window.size() && node == next->get_path()[0];
            };

            auto terminate_branch = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                // common::logger->info("Check: n: {}, d: {}, qd: {} / {}, s: {}",
                //                      graph.get_node_sequence(node), dist, query_dist, query_window.size(),
                //                      get_score(cost, query_dist, dist));

                if (dist == 0 && query_dist == 0)
                    return false;

                if (query_dist == query_window.size()) {
                    aln_found |= start_backtracking(cost, data, query_dist, node);
                    return true;
                }

                // if (static_cast<ssize_t>(query_window.size()) - query_dist < static_cast<ssize_t>(next->get_seed().size()) - num_matches)
                if (query_window.size() - query_dist <= next->get_seed().size() && num_matches < next->get_seed().size()
                        && query_window.size() - query_dist < next->get_seed().size() - num_matches) {
                    // std::cerr << "ttt: " << query_window.size() << "," << query_dist << "\t" << next->get_seed().size() << "," << num_matches << "\n";
                    return true;
                }

                DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
                if (score <= 0)
                    return true;

                // if (num_matches < next->get_seed().size())
                //     return false;

                return false;

                // if (score <= 0 || config_.xdrop < local_best_score - score)
                //     return true;



                // // first check if we've hit another anchor
                // auto jt = node_to_anchors.find(node);
                // if (jt == node_to_anchors.end())
                //     return false;

                // for (size_t j : jt->second) {
                //     if (found_next_anchor(j, query_dist, num_matches) && !skipped[j] && is_last_anchor(j))
                //         return true;
                // }

                // return false;
            };

            align_bwd(
                query_.get_graph(),
                config_,
                aln.get_path()[0],
                last->get_seed().size(),
                query_window,
                config_.max_dist_between_seeds,
                [&](auto&& bt_path, auto&& bt_cigar, size_t cost) { // callback
                    if (bt_path[0] != next->get_path()[0])
                        return;

                    // size_t dist = bt_path.size();
                    // size_t query_dist = bt_cigar.get_num_query();
                    // assert(dist == cigar.get_num_ref());
                    // auto added_score = get_score(cost, dist, query_dist) - last->get_score();

                    Cigar cigar(Cigar::CLIPPED, next->get_clipping());
                    cigar.append(std::move(bt_cigar));
                    std::cerr << "\textension: " << *next << "\t" << cigar.to_string() << "\t" << aln << std::endl;

                    Cigar aln_cigar = aln.get_cigar();
                    aln_cigar.trim_clipping();
                    cigar.append(std::move(aln_cigar));
                    bt_path.insert(bt_path.end(), aln.get_path().begin(), aln.get_path().end());

                    Alignment next_aln(query_.get_graph(),
                                    aln.get_query(),
                                    aln.get_orientation(),
                                    std::move(bt_path),
                                    config_,
                                    std::move(cigar),
                                    aln.get_end_trim());
                    std::cerr << "\t\t" << next_aln << std::endl;
                    aln_found = true;
                    callback(std::move(next_aln));
                },
                start_backtracking,
                terminate_branch,
                [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) { // terminate
                    std::ignore = cost;
                    std::ignore = data;
                    std::ignore = query_dist;
                    // const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                    // common::logger->info("Checkt: n: {}, d: {}, qd: {} / {}, s: {}",
                    //                  graph.get_node_sequence(node), dist, query_dist, query_window.size(),
                    //                  get_score(cost, query_dist, dist));
                    return aln_found || (query_dist == query_window.size() && node == next->get_path()[0]);
                    // // terminate
                    // const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                    // DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
                    // best_score = std::max(best_score, score);
                    // local_best_score = std::max(local_best_score, score);

                    // max_dist = std::max(dist, max_dist);
                    // auto jt = all_nodes_to_anchors.find(node);
                    // if (jt == all_nodes_to_anchors.end())
                    //     return false;

                    // bool found = false;
                    // for (const auto &[j, i] : jt->second) {
                    //     auto kt = begin + j;
                    //     if (found_next_anchor(j, query_dist - kt->get_path().size() + 1 + i, num_matches + i)) {
                    //         found = true;
                    //         if (is_last_anchor(j))
                    //             return true;
                    //     }
                    // }

                    // return found && query_dist == query_window.size() && dist == max_dist;
                }
            );
            // callback(std::move(aln));
            // std::ignore = this;
            // // std::cerr << aln << "\t" << *last << "\t" << *next << "\n";
            // std::ignore = aln;
            // std::ignore = last_to_next_dist;
            // std::ignore = score_up_to_now;
        },
        [&alignments](Alignment&& aln) {
            aln.trim_end();
            alignments.emplace_back(std::move(aln));
        }
    );

    size_t num_selected = sdsl::util::cnt_one_bits(selected);
    common::logger->info("Selected {} / {} anchors from {} chains. Produced {} alignments",
                            num_selected,
                            selected.size(),
                            num_chains, alignments.size());

    return alignments;
    // // common::logger->info("Computing alignments between {} anchors", anchors.size());
    // auto connections = get_connections(anchors);
    // // common::logger->info("DONE");
    // std::vector<Alignment> alignments;

    // using AnchorIt = std::vector<Anchor>::iterator;

    // std::sort(anchors.begin(), anchors.end(), AnchorLess<Anchor>());
    // bool first_chain = true;
    // DBGAlignerConfig::score_t first_chain_score = DBGAlignerConfig::ninf;
    // chain_anchors<AnchorIt>(query_, config_, anchors.begin(), anchors.end(),
    //     [this,&connections](const Anchor &a_j,
    //                         ssize_t max_dist,
    //                         AnchorIt begin,
    //                         AnchorIt end,
    //                         typename ChainScores<AnchorIt>::iterator chain_scores,
    //                         const ScoreUpdater<AnchorIt> &update_score) {
    //         // std::cerr << "Anchor: " << a_j << "\n";
    //         auto find_j = connections.find(a_j);
    //         if (find_j == connections.end())
    //             return;

    //         // std::cerr << "End Anchor: " << a_j << "\n";

    //         // const DeBruijnGraph &graph = query_.get_graph();
    //         std::string_view query_j = a_j.get_seed();
    //         // const DBGAlignerConfig::score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));

    //         for (auto it = begin; it != end; ++it) {
    //             const Anchor &a_i = *it;
    //             assert(a_i.get_label_class() == a_j.get_label_class());
    //             assert(a_i.get_orientation() == a_j.get_orientation());

    //             std::string_view query_i = a_i.get_seed();

    //             DBGAlignerConfig::score_t dist = query_j.end() - query_i.end();

    //             if (dist <= 0 || query_i.begin() >= query_j.begin()) {
    //                 ++chain_scores;
    //                 continue;
    //             }

    //             // std::cerr << "checking connect " << a_i << " -> " << a_j << "\n";

    //             auto find_i = find_j->second.find(a_i);
    //             if (find_i == find_j->second.end() || find_i->second.empty()) {
    //                 ++chain_scores;
    //                 continue;
    //             }
    //             // std::cerr << "\tfound!\n";

    //             DBGAlignerConfig::score_t base_score = std::get<0>(*chain_scores);

    //             DBGAlignerConfig::score_t score = DBGAlignerConfig::ninf;
    //             DBGAlignerConfig::score_t coord_dist = dist;
    //             // size_t a_i_chars = graph.get_k() - a_i.get_end_trim();
    //             // size_t a_j_chars = a_j.get_spelling().size();
    //             bool updated = false;
    //             // std::cerr << "Q: " << dist << "\tcd: " << coord_dist << "\n";
    //             for (const auto &[i, added_score, path, cigar] : find_i->second) {
    //                 // std::cerr << "check\t" << a_i << "\t->\t" << cigar.to_string() << "\t->\t" << a_j << "\n";
    //                 // size_t query_dist = cigar.get_num_query();
    //                 // size_t path_dist = path.size();
    //                 // assert(path_dist == cigar.get_num_ref());
    //                 // query_dist is the distance from the front of a_j to the ith node of a_i
    //                 // ssize_t offset = a_j.get_seed().size() + i;
    //                 // offset -= a_i.get_seed().size();
    //                 // DBGAlignerConfig::score_t cur_dist = query_dist + offset;
    //                 // std::cerr << "\t\tp: " << path_dist << "\tq: " << query_dist << "\ti: " << i << "\t-> " << cur_dist;

    //                 std::string_view a_i_front(a_i.get_seed().begin() + i, a_i.get_seed().size() - i);
    //                 DBGAlignerConfig::score_t front_score = config_.match_score(a_i_front);

    //                 if (base_score + added_score - front_score + a_j.get_score() > score) {
    //                     // DBGAlignerConfig::score_t cur_coord_dist = path_dist + offset;
    //                     // std::cerr << "," << cur_coord_dist;

    //                     score = base_score + added_score - front_score + a_j.get_score();
    //                     coord_dist = path.size() + i + a_j.get_seed().size() - a_i.get_seed().size();
    //                     updated = true;
    //                     // std::cerr << "\tworked\t" << coord_dist << "\t" << base_score << "+" << added_score << "-" << front_score << "+" << a_j.get_score() << "=" << score << "\n";
    //                 }
    //                 // std::cerr << "\n";
    //             }

    //             // std::cerr << "connect " << a_i << "\t->\t" << a_j << "\t"
    //             //                         << fmt::format("{}", fmt::join(find_i->second, ",")) << "\t" << score << "\t" << dist << "," << coord_dist << std::endl;

    //             // std::cerr << "checking connect2 " << a_i << " -> " << a_j << "\t" << updated << "\t" << base_score << "->" << score << " vs. " << std::get<0>(*(chain_scores + (end - begin))) << "\t" << dist << "," << coord_dist << "\n";
    //             if (updated)
    //                 update_score(score, it, coord_dist);

    //             ++chain_scores;
    //         }
    //     },
    //     [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
    //         assert(chain.size());

    //         if (first_chain) {
    //             first_chain_score = score_traceback.back();
    //             assert(score_traceback.front() <= score_traceback.back());
    //         }

    //         bool ret_val = (first_chain || score_traceback.back() == first_chain_score || chain.size() > 1 || chain[0].first->get_seed().size() > config_.min_seed_length || (!chain[0].first->get_clipping() && !chain[0].first->get_end_clipping()));
    //         first_chain = false;

    //         // if (ret_val) {
    //             // common::logger->info("Chain\t{}\t{}", score_traceback.back(), ret_val);
    //             // for (const auto &[it, dist] : chain) {
    //             //     std::cerr << "\t" << *it << "\t" << dist;
    //             // }
    //             // std::cerr << std::endl;
    //         // }

    //         return ret_val;
    //     },
    //     [this,&connections](AnchorIt last,
    //                         AnchorIt next,
    //                         Alignment&& aln,
    //                         size_t last_to_next_dist,
    //                         DBGAlignerConfig::score_t score_up_to_now,
    //                         const AlignmentCallback &callback) {
    //         auto find_j = connections.find(*last);
    //         assert(find_j != connections.end());
    //         auto find_i = find_j->second.find(*next);
    //         assert(find_i != find_j->second.end() && find_i->second.size());

    //         Cigar base_cigar(Cigar::CLIPPED, next->get_clipping());

    //         Cigar aln_cigar = aln.get_cigar();
    //         aln_cigar.trim_clipping();

    //         for (const auto &[i, added_score, bt_path, bt_cigar] : find_i->second) {
    //             // std::cerr << "try\t" << *next << "\t->\t" << bt_cigar.to_string() << "\t->\t" << aln << std::endl;
    //             size_t coord_dist = bt_path.size() + i + last->get_seed().size() - next->get_seed().size();
    //             // size_t path_dist = bt_path.size();
    //             // assert(path_dist == bt_cigar.get_num_ref());

    //             // size_t offset = path_dist + last->get_seed().size() - i;
    //             // size_t cur_coord_dist = offset - next->get_seed().size();

    //             if (coord_dist == last_to_next_dist) {
    //                 // std::cerr << "\tconnected\n";
    //                 std::vector<DeBruijnGraph::node_index> path(next->get_path().begin(), next->get_path().begin() + i);
    //                 path.insert(path.end(), bt_path.begin(), bt_path.end());
    //                 path.insert(path.end(), aln.get_path().begin(), aln.get_path().end());

    //                 auto cigar = base_cigar;
    //                 cigar.append(Cigar::MATCH, i);
    //                 cigar.append(Cigar(bt_cigar));
    //                 cigar.append(Cigar(aln_cigar));

    //                 Alignment next_aln(query_.get_graph(),
    //                                 aln.get_query(),
    //                                 aln.get_orientation(),
    //                                 std::move(path),
    //                                 config_,
    //                                 std::move(cigar),
    //                                 aln.get_end_trim());
    //                 // std::cerr << "\tworked!\t" << next_aln << "\n";
    //                 callback(std::move(next_aln));
    //             }
    //         }
    //     },
    //     [&alignments](Alignment&& alignment) {
    //         // std::cerr << "final aln\t" << alignment << "\n";
    //         alignments.emplace_back(std::move(alignment));
    //     }
    // );

    // if (alignments.size()) {
    //     std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
    //         return std::make_pair(a.get_score(), b.get_orientation())
    //              > std::make_pair(b.get_score(), a.get_orientation());
    //     });

    //     alignments.resize(1);
    // }

    // // std::cerr << "done\n";
    // return alignments;
}

// auto ExactSeeder::get_connections(std::vector<Anchor> &anchors) const -> Connections {
//     if (anchors.empty())
//         return {};

//     Connections connections;
//     if (anchors.size() < 2)
//         return connections;

//     std::sort(anchors.begin(), anchors.end(), [&](const auto &a, const auto &b) {
//         return std::make_tuple(b.get_orientation(), b.get_label_class(), a.get_clipping())
//              > std::make_tuple(a.get_orientation(), a.get_label_class(), b.get_clipping());
//     });

//     DBGAlignerConfig::score_t match_score = config_.match_score("A");

//     auto chain = [&](auto begin, auto end) {
//         if (begin == end)
//             return;

//         DBGAlignerConfig::score_t best_score = 0;
//         sdsl::bit_vector skipped(end - begin);
//         skipped[end - begin - 1] = true;
//         size_t smallest_clipping = std::numeric_limits<size_t>::max();
//         tsl::hopscotch_map<DeBruijnGraph::node_index, tsl::hopscotch_set<size_t>> node_to_anchors;
//         tsl::hopscotch_map<DeBruijnGraph::node_index, std::vector<std::pair<size_t, size_t>>> all_nodes_to_anchors;
//         for (auto it = begin; it != end; ++it) {
//             smallest_clipping = std::min(smallest_clipping, it->get_clipping());
//             node_to_anchors[it->get_path().back()].emplace(it - begin);
//             for (size_t i = 0; i < it->get_path().size(); ++i) {
//                 all_nodes_to_anchors[it->get_path()[i]].emplace_back(it - begin, i);
//             }
//         }

//         std::string_view global_query_window = begin->get_query();
//         global_query_window.remove_prefix(smallest_clipping);
//         size_t num_extended = 0;
//         for (auto it = begin; it != end; ++it) {
//             assert(it->get_query().end() == global_query_window.end());
//             if (!it->get_clipping() || skipped[it - begin])
//                 continue;

//             ++num_extended;

//             skipped[it - begin] = true;
//             // std::cerr << "Starting from\t" << *it << "\n";
//             DBGAlignerConfig::score_t local_best_score = it->get_score();
//             size_t next_clipping = std::numeric_limits<size_t>::max();
//             if (it->get_clipping() > (it + 1)->get_clipping())
//                 next_clipping = (it + 1)->get_clipping();

//             DeBruijnGraph::node_index start_node = it->get_path()[0];
//             std::string_view query_window = global_query_window;
//             query_window.remove_suffix(it->get_end_clipping() + it->get_seed().size());
//             size_t max_dist = 0;

//             auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
//                 return it->get_score() + cost_to_score(cost, query_dist, dist, match_score)
//                         + (query_dist == query_window.size() + smallest_clipping ? config_.left_end_bonus : 0);
//             };

//             auto start_backtracking = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
//                 const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
//                 if (dist == 0 || query_dist == 0 || num_matches < config_.min_seed_length)
//                     return false;

//                 // first check if we've hit another anchor
//                 auto jt = all_nodes_to_anchors.find(node);
//                 if (jt == all_nodes_to_anchors.end())
//                     return false;

//                 for (const auto &[j, i] : jt->second) {
//                     auto kt = begin + j;

//                     // std::cerr << "Trying to compute " << *kt << "; " << i << " -> " << *it << "\t" << dist << "," << query_dist << "\n";
//                     assert(kt > it || skipped[j]);

//                     if (num_matches + i < kt->get_seed().size() || kt <= it)
//                         continue;

//                     if (it->get_clipping() != kt->get_clipping() + i + query_dist)
//                         continue;

//                     return true;
//                 }

//                 return false;
//             };

//             auto found_next_anchor = [&](size_t j, size_t query_dist, size_t num_matches) {
//                 if (num_matches < config_.min_seed_length)
//                     return false;

//                 auto kt = begin + j;
//                 assert(kt > it || skipped[j]);

//                 // std::cerr << "Checking " << (kt - begin) << " vs. " << (it - begin) << "\t" << (kt <= it) << "\t" << *kt << " -> " << *it << "\t" << dist << "," << query_dist << "," << num_matches << "\n";
//                 // std::cerr << "\tfoo\t" << it->get_clipping() - query_dist << " vs. " << kt->get_clipping() + kt->get_path().size() - 1 << "\n";

//                 if (num_matches < kt->get_seed().size() - kt->get_path().size() + 1 || kt <= it)
//                     return false;

//                 // std::cerr << "\tbar\t" << (it->get_clipping() - query_dist == kt->get_clipping() + kt->get_path().size() - 1) << "\n";

//                 if (it->get_clipping() != kt->get_clipping() + kt->get_path().size() - 1 + query_dist)
//                     return false;

//                 return true;
//             };

//             auto is_last_anchor = [&](size_t j) {
//                 bool all_prev_skipped = true;
//                 auto kt = begin + j;
//                 for (auto lt = it + 1; lt != kt; ++lt) {
//                     if (lt->get_clipping() > kt->get_clipping() && !skipped[lt - begin]) {
//                         all_prev_skipped = false;
//                         break;
//                     }
//                 }

//                 return all_prev_skipped;
//             };

//             auto terminate_branch = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
//                 const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
//                 if (dist == 0 && query_dist == 0)
//                     return false;

//                 if (query_dist == query_window.size())
//                     return true;

//                 DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
//                 if (score <= 0 || config_.xdrop < local_best_score - score)
//                     return true;

//                 if (num_matches < config_.min_seed_length)
//                     return false;

//                 // first check if we've hit another anchor
//                 auto jt = node_to_anchors.find(node);
//                 if (jt == node_to_anchors.end())
//                     return false;

//                 for (size_t j : jt->second) {
//                     if (found_next_anchor(j, query_dist, num_matches) && !skipped[j] && is_last_anchor(j))
//                         return true;
//                 }

//                 return false;
//             };

//             align_bwd(
//                 query_.get_graph(),
//                 config_,
//                 start_node,
//                 it->get_seed().size(),
//                 query_window,
//                 config_.max_dist_between_seeds,
//                 [&](auto&& path, auto&& cigar, size_t cost) {
//                     size_t dist = path.size();
//                     size_t query_dist = cigar.get_num_query();
//                     assert(dist == cigar.get_num_ref());
//                     DeBruijnGraph::node_index node = path.front();
//                     auto jt = all_nodes_to_anchors.find(node);
//                     assert(jt != all_nodes_to_anchors.end());
//                     bool found = false;
//                     for (const auto &[j, i] : jt->second) {
//                         size_t num_matches = 0;
//                         assert(cigar.size());
//                         if (cigar.data()[0].first == Cigar::MATCH) {
//                             num_matches += cigar.data()[0].second;
//                             if (cigar.size() == 1)
//                                 num_matches += it->get_seed().size();
//                         }

//                         auto kt = begin + j;
//                         if (found_next_anchor(j, query_dist - kt->get_path().size() + 1 + i, num_matches + i)) {
//                             // std::cerr << "foundfo\t" << *kt << "\t->\t" << cigar.to_string() << "\t->\t" << *it << "\n";
//                             connections[*it][*kt].emplace_back(i, get_score(cost, dist, query_dist) - it->get_score(), path, cigar);
//                             found = true;
//                         }
//                     }
//                     assert(found);
//                 },
//                 start_backtracking,
//                 terminate_branch,
//                 [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
//                     // terminate
//                     const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
//                     DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
//                     best_score = std::max(best_score, score);
//                     local_best_score = std::max(local_best_score, score);

//                     max_dist = std::max(dist, max_dist);
//                     auto jt = all_nodes_to_anchors.find(node);
//                     if (jt == all_nodes_to_anchors.end())
//                         return false;

//                     bool found = false;
//                     for (const auto &[j, i] : jt->second) {
//                         auto kt = begin + j;
//                         if (found_next_anchor(j, query_dist - kt->get_path().size() + 1 + i, num_matches + i)) {
//                             found = true;
//                             if (is_last_anchor(j))
//                                 return true;
//                         }
//                     }

//                     return found && query_dist == query_window.size() && dist == max_dist;
//                 }
//             );
//         }
//     };

//     auto begin = anchors.begin();
//     for (auto it = begin; it != anchors.end(); ++it) {
//         if (it != begin && it->get_orientation() != (it - 1)->get_orientation()) {
//             chain(begin, it);
//             begin = it;
//         }
//     }

//     chain(begin, anchors.end());

//     return connections;
// }

} // namespace mtg::graph::align
