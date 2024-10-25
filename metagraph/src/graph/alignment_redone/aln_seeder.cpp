#include "aln_seeder.hpp"

#include <tsl/hopscotch_set.h>

#include "aln_chainer.hpp"
#include "common/logger.hpp"
#include "common/vector_map.hpp"
#include "graph/representation/succinct/boss.hpp"

namespace mtg::graph::align_redone {

std::vector<Alignment> Seeder::get_inexact_anchors() const {
    auto anchors = get_anchors();
    std::vector<Alignment> inexact_anchors;
    inexact_anchors.reserve(anchors.size());
    std::transform(anchors.begin(), anchors.end(), std::back_inserter(inexact_anchors),
                   [&](const Anchor &a) { return Alignment(a); });
    return inexact_anchors;
}

std::vector<Anchor> ExactSeeder::get_anchors() const {
    const DeBruijnGraph &graph = query_.get_graph();
    std::vector<Anchor> anchors;

    if (query_.get_query().size() < config_.min_seed_length)
        return anchors;

    if (config_.min_seed_length == graph.get_k()) {
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

        return anchors;
    }

    return anchors;
}

void call_distances(const DeBruijnGraph &graph,
                    DeBruijnGraph::node_index start,
                    DeBruijnGraph::node_index target,
                    const std::function<void(ssize_t)> &callback,
                    const std::function<bool(ssize_t)> &terminate_branch = [](ssize_t) { return false; },
                    const std::function<bool()> &terminate = []() { return false; }) {
    // TODO: this is a slow DFS-based method, replace with something smarter
    std::vector<std::pair<DeBruijnGraph::node_index, ssize_t>> traversal_stack;
    traversal_stack.emplace_back(start, 0);
    while (traversal_stack.size() && !terminate()) {
        auto [node, coord_dist] = traversal_stack.back();
        traversal_stack.pop_back();

        if (node == target)
            callback(coord_dist);

        if (terminate_branch(coord_dist))
            continue;

        ++coord_dist;
        graph.adjacent_incoming_nodes(node, [&](DeBruijnGraph::node_index next) {
            traversal_stack.emplace_back(next, coord_dist);
        });
    }
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
    return (match_score * (query_size + match_size) - cost) / 2;
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
void align_impl(const std::function<size_t(DeBruijnGraph::node_index)> &has_single_incoming,
                const std::function<DeBruijnGraph::node_index(DeBruijnGraph::node_index, char)> &traverse_back,
                const std::function<void(DeBruijnGraph::node_index, const std::function<void(DeBruijnGraph::node_index, char)>&)> &call_incoming_kmers,
                const DBGAlignerConfig &config,
                const DeBruijnGraph::node_index start_node,
                size_t start_num_matches,
                const StrItr query_window_begin,
                const StrItr query_window_end,
                size_t max_dist,
                const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&)> &callback,
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

        bool inserted = false;

        if constexpr(std::is_same_v<table_t, SMap>)
            inserted = (std::get<3>(bucket) == Cigar::CLIPPED);

        if constexpr(std::is_same_v<table_t, FMap>)
            inserted = (std::get<2>(bucket) == DeBruijnGraph::npos);

        if constexpr(std::is_same_v<table_t, EMap>)
            inserted = (std::get<1>(bucket) == 0);

        if (inserted || dist > std::get<0>(bucket)) {
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
        }

        if constexpr(std::is_same_v<table_t, SMap>) {
            assert(std::get<3>(bucket) == Cigar::MATCH
                || std::get<3>(bucket) == Cigar::MISMATCH
                || std::get<3>(bucket) == Cigar::INSERTION
                || std::get<3>(bucket) == Cigar::DELETION);
        }

        return std::get<0>(bucket);
    };

    auto fill_table = [&]() {
        {
            size_t query_dist = 0;
            DeBruijnGraph::node_index node = start_node;
            SMap data(0, 0, start_node, Cigar::MATCH, *(query_window_rbegin - 1), start_num_matches);
            auto &[best_dist, last_ext, last_node, last_op, c, num_matches] = data;

            for (auto it = query_window_rbegin; it != query_window_rend; ++it) {
                if (terminate_branch(0, data, query_dist, node) || best_dist == max_dist || !has_single_incoming(node))
                    break;

                if (auto prev = traverse_back(node, *it)) {
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

            terminate(0, data, query_dist, node);
        }

        for (size_t cost = 0; cost < std::max({ S.size(), E.size(), F.size() }); ++cost) {
            // switch to matches from long insertions
            if (cost < E.size()) {
                for (size_t query_dist = E[cost].offset(); query_dist < E[cost].size(); ++query_dist) {
                    assert(query_dist <= query_size);
                    if (query_dist == query_size)
                        break;

                    for (const auto &[node, data] : E[cost][query_dist]) {
                        const auto &[best_dist, last_num_ops] = data;

                        if (last_num_ops == 0 || best_dist == max_dist)
                            continue;

                        if (terminate_branch(cost, SMap(best_dist, 0, node, Cigar::INSERTION, '\0', 0), query_dist, node))
                            continue;

                        size_t num_matches = 0;

                        std::vector<std::pair<node_index, char>> prevs;
                        call_incoming_kmers(node, [&](auto prev, char c) {
                            prevs.emplace_back(prev, c);
                        });

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
                                if (terminate_branch(cur_cost, data, cur_query_dist, prev) || cur_best == max_dist || !has_single_incoming(prev))
                                    break;

                                if (auto pprev = traverse_back(prev, *jt)) {
                                    ++cur_best;
                                    ++cur_query_dist;
                                    ++cur_num_matches;
                                    prev = pprev;
                                } else {
                                    break;
                                }
                            }

                            set_value(
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
                            );
                        }
                    }
                }
            }

            // switch to matches from long deletions
            if (cost < F.size()) {
                for (size_t query_dist = F[cost].offset(); query_dist < F[cost].size(); ++query_dist) {
                    assert(query_dist <= query_size);
                    if (query_dist == query_size)
                        break;

                    for (const auto &[node, data] : F[cost][query_dist]) {
                        const auto &[best_dist, last_num_ops, last_node] = data;

                        assert((last_node == DeBruijnGraph::npos) == (last_num_ops == 0));

                        if (last_num_ops == 0 || best_dist == max_dist)
                            continue;

                        if (terminate_branch(cost, SMap(best_dist, 0, last_node, Cigar::DELETION, '\0', 0), query_dist, node))
                            continue;

                        size_t num_matches = 0;

                        std::vector<std::pair<node_index, char>> prevs;
                        call_incoming_kmers(node, [&](auto prev, char c) {
                            prevs.emplace_back(prev, c);
                        });

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
                                if (terminate_branch(cur_cost, data, cur_query_dist, prev) || cur_best == max_dist || !has_single_incoming(prev))
                                    break;

                                if (auto pprev = traverse_back(prev, *jt)) {
                                    ++cur_best;
                                    ++cur_query_dist;
                                    ++cur_num_matches;
                                    prev = pprev;
                                } else {
                                    break;
                                }
                            }

                            set_value(
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
                            );
                        }
                    }
                }
            }

            if (cost >= S.size())
                continue;

            for (size_t query_dist = S[cost].offset(); query_dist < S[cost].size(); ++query_dist) {
                assert(query_dist <= query_size);
                size_t it_dist = 0;
                for (auto it = S[cost][query_dist].begin(); it != S[cost][query_dist].end(); ++it, ++it_dist) {
                    node_index node = it->first;
                    size_t best_dist = std::get<0>(it->second);
                    Cigar::Operator last_op = std::get<3>(it->second);
                    size_t num_matches = std::get<5>(it->second);

                    if (last_op == Cigar::CLIPPED)
                        continue;

                    terminate(cost, it->second, query_dist, node);
                    if (terminate_branch(cost, it->second, query_dist, node))
                        continue;

                    // forward creation of insertions
                    if (query_dist < query_size) {
                        // extend a previous insertion
                        if (const auto *ins_ext_bucket = get_bucket(E, cost, query_dist, node)) {
                            auto [last_dist, last_num_ops] = *ins_ext_bucket;
                            if (last_num_ops > 0) {
                                ssize_t next_ext_cost = cost + gap_ext;
                                size_t stored_dist = set_value(E, next_ext_cost, query_dist + 1, node, last_dist, node, last_num_ops + 1, Cigar::INSERTION, '\0', 0);
                                set_value(S, next_ext_cost, query_dist + 1, node, stored_dist, node, 0, Cigar::INSERTION, '\0', 0);
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }

                        // open an insertion
                        if (last_op != Cigar::DELETION) {
                            ssize_t next_opn_cost = cost + gap_opn;
                            size_t stored_dist = set_value(E, next_opn_cost, query_dist + 1, node, best_dist, node, 1, Cigar::INSERTION, '\0', 0);
                            set_value(S, next_opn_cost, query_dist + 1, node, stored_dist, node, 0, Cigar::INSERTION, '\0', 0);
                            it = S[cost][query_dist].begin() + it_dist;
                        }
                    }

                    if (best_dist < max_dist) {
                        // deletion
                        std::vector<std::pair<node_index, char>> prevs;
                        call_incoming_kmers(node, [&](auto prev, char c) {
                            prevs.emplace_back(prev, c);
                        });
                        if (prevs.size()) {
                            if (const auto *del_ext_bucket = get_bucket(F, cost, query_dist, node)) {
                                // extension
                                auto [last_dist, last_num_ops, last_node] = *del_ext_bucket;
                                if (last_node != DeBruijnGraph::npos) {
                                    assert(last_dist >= last_num_ops);
                                    ssize_t next_ext_cost = cost + gap_ext;
                                    for (const auto &[prev, c] : prevs) {
                                        size_t stored_value = set_value(F, next_ext_cost, query_dist, prev, last_dist + 1, node, last_num_ops + 1, Cigar::DELETION, '\0', 0);
                                        set_value(S, next_ext_cost, query_dist, prev, stored_value, node, 0, Cigar::DELETION, '\0', 0);
                                        it = S[cost][query_dist].begin() + it_dist;
                                    }
                                }
                            }

                            if (last_op != Cigar::INSERTION) {
                                ssize_t next_opn_cost = cost + gap_opn;
                                for (const auto &[prev, c] : prevs) {
                                    size_t stored_value = set_value(F, next_opn_cost, query_dist, prev, best_dist + 1, node, 1, Cigar::DELETION, '\0', 0);
                                    set_value(S, next_opn_cost, query_dist, prev, stored_value, node, 0, Cigar::DELETION, '\0', 0);
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
                                    if (terminate_branch(cur_cost, data, cur_query_dist, prev) || cur_best == max_dist || !has_single_incoming(prev))
                                        break;

                                    if (auto pprev = traverse_back(prev, *jt)) {
                                        ++cur_best;
                                        ++cur_query_dist;
                                        ++cur_num_matches;
                                        prev = pprev;
                                    } else {
                                        break;
                                    }
                                }

                                set_value(
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
                                );
                                it = S[cost][query_dist].begin() + it_dist;

                                // terminate(cur_cost, data, cur_query_dist, prev);
                            }
                        }
                    }

                    if (terminate(cost, it->second, query_dist, node))
                        return;
                }
            }
        }
    };

    fill_table();

    // backtrack
    for (size_t start_cost = 0; start_cost < S.size(); ++start_cost) {
        std::vector<std::tuple<size_t, size_t, size_t, DeBruijnGraph::node_index>> starts;
        for (size_t start_query_dist = S[start_cost].offset(); start_query_dist < S[start_cost].size(); ++start_query_dist) {
            assert(start_query_dist <= query_size);
            for (auto it = S[start_cost][start_query_dist].begin(); it != S[start_cost][start_query_dist].end(); ++it) {
                DeBruijnGraph::node_index node = it->first;
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = it->second;
                if (last_op != Cigar::CLIPPED && start_backtrack(start_cost, it->second, start_query_dist, node))
                    starts.emplace_back(dist, std::abs(static_cast<ssize_t>(dist - start_query_dist)), start_query_dist, node);
            }
        }

        std::sort(starts.begin(), starts.end(), [&](const auto &a, const auto &b) {
            return std::make_pair(std::get<0>(b), std::get<1>(a))
                    < std::make_pair(std::get<0>(a), std::get<1>(b));
        });

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
                    assert(ins_dist == dist);
                    assert(query_dist >= num_ins);
                    assert(num_ins > 0);

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
                    assert(del_dist == dist);
                    assert(dist >= num_del);
                    assert(num_del > 0);
                    assert(prev_node != DeBruijnGraph::npos);

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
                    if (std::get<3>(*bt) == Cigar::MISMATCH) {
                        assert(mismatch_char != *it);
                        traverse_node = traverse_back(traverse_node, mismatch_char);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        ++it;
                        --num_match;
                    } else {
                        assert(mismatch_char == *it);
                    }

                    while (num_match) {
                        assert(it != query_window_rend);
                        traverse_node = traverse_back(traverse_node, *it);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        ++it;
                        --num_match;
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
            callback(std::move(path), std::move(cigar));
        }
    }
}

void align_bwd(const DeBruijnGraph &graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate) {
    align_impl(
        [&graph](DeBruijnGraph::node_index node) { return graph.has_single_incoming(node); },
        [&graph](DeBruijnGraph::node_index node, char c) { return graph.traverse_back(node, c); },
        [&graph](DeBruijnGraph::node_index node, const auto &callback) {
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
               const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate) {
    align_impl(
        [&graph](DeBruijnGraph::node_index node) { return graph.has_single_outgoing(node); },
        [&graph](DeBruijnGraph::node_index node, char c) { return graph.traverse(node, c); },
        [&graph](DeBruijnGraph::node_index node, const auto &callback) {
            graph.call_outgoing_kmers(node, callback);
        },
        config,
        start_node,
        num_matches,
        query_window.rbegin(), query_window.rend(),
        max_dist,
        [&](auto&& path, auto&& cigar) {
            std::reverse(path.begin(), path.end());
            std::reverse(cigar.data().begin(), cigar.data().end());
            callback(std::move(path), std::move(cigar));
        },
        start_backtrack,
        terminate_branch,
        terminate
    );
}

void Extender::extend(const Alignment &aln, const std::function<void(Alignment&&)> &callback) const {
    DBGAlignerConfig::score_t match_score = config_.match_score("A");

    std::vector<Alignment> fwd_exts;
    if (!aln.get_end_clipping()) {
        fwd_exts.emplace_back(aln);
    } else {
        std::string_view query_window = aln.get_query();
        query_window.remove_prefix(aln.get_clipping() + aln.get_seed().size());

        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return aln.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.right_end_bonus : 0);
        };

        auto it = std::make_reverse_iterator(aln.get_cigar().data().end());
        auto end = std::make_reverse_iterator(aln.get_cigar().data().begin());
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;
        DBGAlignerConfig::score_t best_score = aln.get_score();
        align_fwd(
            query_.get_graph(), config_, aln.get_path().back(), num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& path, auto&& cigar) {
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
                size_t num_matches = std::get<5>(data);
                if ((dist == 0 && query_dist == 0) || !num_matches)
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
                auto score = get_score(cost, dist, query_dist);
                best_score = std::max(best_score, score);
                return query_dist == query_window.size();
            }
        );

        if (fwd_exts.empty())
            fwd_exts.emplace_back(aln);
    }

    std::vector<Alignment> alns;
    for (auto &fwd_ext : fwd_exts) {
        alns.emplace_back(fwd_ext);

        if (!fwd_ext.get_clipping())
            continue;

        std::string_view query_window = fwd_ext.get_query();
        query_window.remove_suffix(fwd_ext.get_end_clipping() + fwd_ext.get_seed().size());

        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return fwd_ext.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.left_end_bonus : 0);
        };

        DBGAlignerConfig::score_t best_score = fwd_ext.get_score();

        auto it = aln.get_cigar().data().begin();
        auto end = aln.get_cigar().data().end();
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;

        align_bwd(
            query_.get_graph(), config_, fwd_ext.get_path()[0], num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& path, auto&& cigar) {
                path.insert(path.end(), fwd_ext.get_path().begin(), fwd_ext.get_path().end());
                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                Cigar ext_cigar(Cigar::CLIPPED, query_window.size() - query_dist);
                ext_cigar.append(std::move(cigar));

                Cigar aln_cigar = fwd_ext.get_cigar();
                aln_cigar.trim_clipping();
                ext_cigar.append(std::move(aln_cigar));

                alns.emplace_back(
                    query_.get_graph(),
                    fwd_ext.get_query(),
                    fwd_ext.get_orientation(),
                    std::move(path),
                    config_,
                    std::move(ext_cigar)
                );

                // assert(alns.back().get_score() == best_score);
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                size_t num_matches = std::get<5>(data);
                if ((dist == 0 && query_dist == 0) || !num_matches)
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
                auto score = get_score(cost, dist, query_dist);
                best_score = std::max(best_score, score);
                return query_dist == query_window.size();
            }
        );
    }

    std::sort(alns.begin(), alns.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_score(), a.get_seed().size())
             > std::make_tuple(b.get_score(), b.get_seed().size());
    });

    for (auto &a : alns) {
        callback(std::move(a));
    }
}

std::vector<Alignment> ExactSeeder::get_inexact_anchors() const {
    std::vector<Anchor> anchors = get_anchors();
    std::vector<Alignment> alignments;

    using AnchorIt = std::vector<Anchor>::iterator;

    DBGAlignerConfig::score_t match_score = config_.match_score("A");

    chain_anchors<AnchorIt>(query_, config_, anchors.begin(), anchors.end(),
        [this](const Anchor &a_j,
               ssize_t max_dist,
               AnchorIt begin,
               AnchorIt end,
               typename ChainScores<AnchorIt>::iterator chain_scores,
               const ScoreUpdater<AnchorIt> &update_score) {
            const DeBruijnGraph &graph = query_.get_graph();
            std::string_view query_j = a_j.get_seed();
            const DBGAlignerConfig::score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));

            for (auto it = begin; it != end; ++it) {
                const Anchor &a_i = *it;
                assert(a_i.get_label_class() == a_j.get_label_class());
                assert(a_i.get_orientation() == a_j.get_orientation());

                std::string_view query_i = a_i.get_seed();

                DBGAlignerConfig::score_t dist = query_j.end() - query_i.end();

                if (dist <= 0 || query_i.begin() >= query_j.begin()) {
                    ++chain_scores;
                    return;
                }

                DBGAlignerConfig::score_t base_score = std::get<0>(*chain_scores);

                if (query_i.end() <= query_j.begin()) {
                    base_score += a_j.get_score();
                } else {
                    std::string_view ext(query_i.data() + query_i.size(),
                                         query_j.end() - query_i.end());
                    base_score += config_.match_score(ext);

                    if (!a_j.get_end_clipping())
                        base_score += config_.right_end_bonus;
                }

                if (base_score <= score_j) {
                    ++chain_scores;
                    return;
                }

                DBGAlignerConfig::score_t score = base_score;
                DBGAlignerConfig::score_t coord_dist = dist;
                ssize_t min_diff = std::numeric_limits<ssize_t>::max();
                size_t a_i_chars = graph.get_k() - a_i.get_end_trim();
                size_t a_j_chars = a_j.get_spelling().size();
                call_distances(graph, a_j.get_path()[0], a_i.get_path().back(),
                    [&](ssize_t path_dist) {
                        DBGAlignerConfig::score_t cur_coord_dist
                            = path_dist + a_j_chars - a_i_chars;
                        ssize_t cur_diff = abs(cur_coord_dist - dist);
                        if (cur_diff < min_diff) {
                            min_diff = cur_diff;
                            coord_dist = cur_coord_dist;
                        }
                    },
                    [&](ssize_t cur_dist) {
                        if (cur_dist > dist) {
                            return cur_dist - dist > std::min(min_diff, max_dist);
                        } else {
                            return dist - cur_dist > max_dist;
                        }
                    },
                    [&]() { return min_diff == 0; }
                );

                if (min_diff == std::numeric_limits<ssize_t>::max()) {
                    ++chain_scores;
                    return;
                }

                DBGAlignerConfig::score_t diff = std::abs(coord_dist - dist);
                if (diff != 0) {
                    score += config_.gap_opening_penalty
                            + (diff - 1) * config_.gap_extension_penalty;
                }

                if (query_i.end() < query_j.begin()) {
                    ssize_t mismatched = query_j.begin() - query_i.end();
                    if (dist > coord_dist)
                        mismatched -= dist - coord_dist;

                    if (mismatched > 0) {
                        std::string_view ext(query_i.data() + query_i.size(), mismatched);
                        score += config_.score_sequences(ext, std::string(mismatched, boss::BOSS::kSentinel));
                    }
                }

                update_score(score, it, coord_dist);

                ++chain_scores;
            }
        },
        [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
            return true;
        },
        [this,match_score](AnchorIt last,
                           AnchorIt next,
                           Alignment&& aln,
                           size_t last_to_next_dist,
                           DBGAlignerConfig::score_t score_up_to_now,
                           const AlignmentCallback &callback) {
            DBGAlignerConfig::score_t best_score = std::numeric_limits<DBGAlignerConfig::score_t>::min();
            size_t next_to_last_dist = last_to_next_dist - last->get_spelling().size() + next->get_spelling().size() - next->get_path().size() + 1;
            std::string_view query_window = aln.get_query();
            query_window.remove_suffix(aln.get_end_clipping() + aln.get_seed().size());

            size_t front_trim = next->get_clipping() + next->get_path().size() - 1;
            query_window.remove_prefix(front_trim);

            auto end_branch = [&](size_t, const SMap &data, size_t query_dist, DeBruijnGraph::node_index) {
                return std::get<0>(data) == next_to_last_dist && query_dist == query_window.size();
            };

            auto reached_end = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                return end_branch(cost, data, query_dist, node) && node == next->get_path().back();
            };

            auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
                return aln.get_score() + cost_to_score(cost, query_dist, dist, match_score)
                        + (query_dist == query_window.size() + front_trim ? config_.left_end_bonus : 0);
            };

            auto start_backtrack = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                if (!reached_end(cost, data, query_dist, node))
                    return false;

                size_t dist = std::get<0>(data);
                DBGAlignerConfig::score_t score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    return true;
                } else {
                    return false;
                }
            };

            align_bwd(query_.get_graph(), config_, aln.get_path()[0], last->get_seed().size(), query_window, next_to_last_dist,
                [&](auto&& bt_path, auto&& bt_cigar) {
                    std::vector<DeBruijnGraph::node_index> path(next->get_path().begin(), next->get_path().end() - 1);
                    path.insert(path.end(), bt_path.begin(), bt_path.end());
                    path.insert(path.end(), aln.get_path().begin(), aln.get_path().end());

                    Cigar cigar(Cigar::CLIPPED, next->get_clipping());
                    cigar.append(Cigar::MATCH, next->get_path().size() - 1);
                    cigar.append(std::move(bt_cigar));

                    Cigar aln_cigar = aln.get_cigar();
                    aln_cigar.trim_clipping();
                    cigar.append(std::move(aln_cigar));

                    callback(Alignment(query_.get_graph(),
                                    aln.get_query(),
                                    aln.get_orientation(),
                                    std::move(path),
                                    config_,
                                    std::move(cigar),
                                    aln.get_end_trim()));
                },
                start_backtrack,
                end_branch,
                reached_end);
        },
        [&alignments](Alignment&& alignment) {
            alignments.emplace_back(std::move(alignment));
        }
    );

    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_score(), a.get_seed().size())
             > std::make_tuple(b.get_score(), b.get_seed().size());
    });
    return alignments;
}

std::vector<Alignment> ExactSeeder::get_alignments() const {
    std::vector<Alignment> alns;
    std::vector<Anchor> anchors = get_anchors();
    if (anchors.empty())
        return alns;

    if (anchors.size() == 1) {
        alns.emplace_back(anchors[0]);
        return alns;
    }

    std::sort(anchors.begin(), anchors.end(), [&](const auto &a, const auto &b) {
        return std::make_tuple(b.get_orientation(), a.get_label_class(), a.get_clipping())
             > std::make_tuple(a.get_orientation(), b.get_label_class(), b.get_clipping());
    });

    DBGAlignerConfig::score_t match_score = config_.match_score("A");

    auto chain = [&](auto begin, auto end) {
        if (begin == end)
            return;

        DBGAlignerConfig::score_t best_score = 0;
        sdsl::bit_vector skipped(end - begin);
        skipped[end - begin - 1] = true;
        size_t smallest_clipping = std::numeric_limits<size_t>::max();
        tsl::hopscotch_map<DeBruijnGraph::node_index, tsl::hopscotch_set<size_t>> node_to_anchors;
        for (auto it = begin; it != end; ++it) {
            smallest_clipping = std::min(smallest_clipping, it->get_clipping());
            node_to_anchors[it->get_path().back()].emplace(it - begin);
        }

        std::string_view global_query_window = begin->get_query();
        global_query_window.remove_prefix(smallest_clipping);
        size_t num_extended = 0;
        for (auto it = begin; it != end; ++it) {
            assert(it->get_query().end() == global_query_window.end());
            if (!it->get_clipping() || skipped[it - begin])
                continue;

            ++num_extended;

            skipped[it - begin] = true;
            DBGAlignerConfig::score_t local_best_score = it->get_score();
            size_t next_clipping = std::numeric_limits<size_t>::max();
            if (it->get_clipping() > (it + 1)->get_clipping())
                next_clipping = (it + 1)->get_clipping();

            DeBruijnGraph::node_index start_node = it->get_path()[0];
            std::string_view query_window = global_query_window;
            query_window.remove_suffix(it->get_end_clipping() + it->get_seed().size());

            auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
                return it->get_score() + cost_to_score(cost, query_dist, dist, match_score)
                        + (query_dist == query_window.size() + smallest_clipping ? config_.left_end_bonus : 0);
            };

            auto start_backtracking = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;

                if (query_dist == query_window.size() && num_matches)
                    return true;

                if (dist == 0 || query_dist == 0 || (last_op != Cigar::MATCH && last_op != Cigar::MISMATCH))
                    return false;

                // first check if we've hit another anchor
                auto jt = node_to_anchors.find(node);
                if (jt == node_to_anchors.end())
                    return false;

                bool start_backtrack = false;
                for (size_t j : jt->second) {
                    if (num_matches >= (begin + j)->get_seed().size()) {
                        start_backtrack |= !skipped[j] || j + 1 == skipped.size()
                            || !(begin + j)->get_clipping()
                            || std::none_of(skipped.begin() + j + 1, skipped.end(), [&](bool s) { return s; });
                    }
                }

                if (!start_backtrack) {
                    DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
                    assert(score <= local_best_score);
                    assert(score <= best_score);
                    start_backtrack = (local_best_score == score);
                }

                return start_backtrack;
            };

            auto terminate_branch = [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                if (dist == 0 || query_dist == 0 || (last_op != Cigar::MATCH && last_op != Cigar::MISMATCH))
                    return false;

                // first check if we've hit another anchor
                auto jt = node_to_anchors.find(node);
                if (jt == node_to_anchors.end())
                    return false;

                DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
                if (score <= 0)
                    return true;

                for (size_t j : jt->second) {
                    auto kt = begin + j;
                    assert(kt > it || skipped[j]);

                    if (num_matches < kt->get_seed().size())
                        continue;

                    if (kt > it && !skipped[kt - begin]) {
                        bool all_prev_skipped = true;
                        for (auto lt = it + 1; lt != kt; ++lt) {
                            if (lt->get_clipping() > kt->get_clipping() && !skipped[lt - begin]) {
                                all_prev_skipped = false;
                                break;
                            }
                        }
                        if (all_prev_skipped) {
                            skipped[kt - begin] = true;
                        }
                    }
                }

                return config_.xdrop < best_score - score;
            };

            align_bwd(
                query_.get_graph(),
                config_,
                start_node,
                it->get_seed().size(),
                query_window,
                std::numeric_limits<size_t>::max(),
                [&](auto&& path, auto&& cigar) {
                    size_t num_matches = cigar.get_num_query();
                    Cigar ext_cigar(Cigar::CLIPPED, query_window.size() - num_matches + smallest_clipping);
                    ext_cigar.append(std::move(cigar));
                    ext_cigar.append(Cigar::MATCH, it->get_seed().size());
                    ext_cigar.append(Cigar::CLIPPED, it->get_end_clipping());

                    path.emplace_back(start_node);

                    alns.emplace_back(query_.get_graph(),
                                    it->get_query(),
                                    it->get_orientation(),
                                    std::move(path),
                                    config_,
                                    std::move(ext_cigar));
                },
                start_backtracking,
                terminate_branch,
                [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index) {
                    // terminate
                    const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                    DBGAlignerConfig::score_t score = get_score(cost, query_dist, dist);
                    best_score = std::max(best_score, score);
                    local_best_score = std::max(local_best_score, score);
                    return query_dist == query_window.size();
                }
            );
        }
    };

    auto begin = anchors.begin();
    for (auto it = begin; it != anchors.end(); ++it) {
        if (it != begin && it->get_orientation() != (it - 1)->get_orientation()) {
            chain(begin, it);
            begin = it;
        }
    }

    chain(begin, anchors.end());

    if (alns.empty())
        return alns;

    std::sort(alns.begin(), alns.end(), [&](const auto &a, const auto &b) {
        return a.get_score() > b.get_score();
    });

    return alns;
}

} // namespace mtg::graph::align