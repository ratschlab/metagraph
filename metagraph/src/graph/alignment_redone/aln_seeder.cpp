#include "aln_seeder.hpp"

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

template <typename Map>
using WaveFront = OffsetVector<Map>; // S(cost)

template <typename Map>
using ScoreTable = std::vector<WaveFront<Map>>;

void align(const DeBruijnGraph &graph,
           const DBGAlignerConfig &config,
           const DeBruijnGraph::node_index start_node,
           const Alignment &aln,
           const Anchor &target,
           size_t dist,
           const std::function<void(std::vector<DeBruijnGraph::node_index>&&, Cigar&&)> &callback) {
    assert(dist > 0);
    dist -= target.get_path().size() - 1;

    using node_index = DeBruijnGraph::node_index;

    std::string_view query_window = aln.get_query();
    query_window.remove_suffix(aln.get_end_clipping() + aln.get_seed().size());
    query_window.remove_prefix(target.get_clipping() + target.get_path().size() - 1);
    if (query_window.empty())
        return;

    common::logger->info("Query window: {}", query_window);

    // derived from 2.4.1 in https://doi.org/10.1101/2022.01.12.476087 (v1)
    DBGAlignerConfig::score_t match_score = config.match_score("A");
    DBGAlignerConfig::score_t mismatch_score = config.score_sequences("A", "T");
    assert(config.gap_opening_penalty <= config.gap_extension_penalty);
    assert(mismatch_score < match_score);

    ssize_t mismatch_cost = (match_score - mismatch_score) * 2;
    ssize_t gap_ext = config.gap_extension_penalty * -2 + match_score;
    ssize_t gap_opn = config.gap_opening_penalty * -2 + gap_ext;

    common::logger->info("=: {}, X: {}, O: {}, E: {} -> X: {}, O: {}, E: {}",
                         match_score, mismatch_score, config.gap_opening_penalty, config.gap_extension_penalty,
                         mismatch_cost, gap_opn, gap_ext);

    assert(mismatch_cost > 0);
    assert(gap_ext > 0);
    assert(gap_opn >= gap_ext);

    // S(cost, query_dist, node) = best_dist
    using SMap = VectorMap<node_index, std::tuple<size_t, size_t, node_index, Cigar::Operator, char>>;
    using EMap = VectorMap<node_index, std::tuple<size_t, size_t>>;
    using FMap = VectorMap<node_index, std::tuple<size_t, size_t, node_index>>;

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
            char c = '\0') -> size_t {
        assert(last_op == Cigar::MATCH
            || last_op == Cigar::MISMATCH
            || last_op == Cigar::INSERTION
            || last_op == Cigar::DELETION);
        if (cost >= table.size())
            table.resize(cost + 1);

        auto &table_slot = table[cost].get(query_dist);
        using table_t = std::decay_t<decltype(table_slot)>;

        bool inserted = !table_slot.count(node);
        auto &bucket = table_slot[node];

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

    auto fill_table = [&]() -> size_t {
        size_t query_dist = 0;
        auto node = start_node;
        size_t best_dist = 0;

        for (auto it = query_window.rbegin(); it != query_window.rend(); ++it) {
            if (best_dist == dist || !graph.has_single_incoming(node))
                break;

            if (auto prev = graph.traverse_back(node, *it)) {
                ++best_dist;
                ++query_dist;
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
            best_dist > 0 ? query_window.back() : '\0'
        );

        if (best_dist == dist && query_dist == query_window.size()) {
            common::logger->info("Early stop");
            if (node == target.get_path().back())
                return 0;

            return std::numeric_limits<size_t>::max();
        }

        for (size_t cost = 0; cost < S.size(); ++cost) {
            if (S[cost].empty())
                continue;

            for (size_t query_dist = S[cost].offset(); query_dist < S[cost].size(); ++query_dist) {
                assert(query_dist <= query_window.size());
                if (S[cost][query_dist].empty())
                    continue;

                size_t it_dist = 0;
                for (auto it = S[cost][query_dist].begin(); it != S[cost][query_dist].end(); ++it, ++it_dist) {
                    node_index node = it->first;
                    auto [best_dist, last_ext, last_node, last_op, mismatch_char] = it->second;

                    common::logger->info("Check: c:{}\tn:{} vs. {}\tcd:{} vs. {}\td:{} vs. {}",
                                         cost,
                                         node, target.get_path().back(),
                                         best_dist, dist,
                                         query_dist, query_window.size());
                    if (best_dist == dist && query_dist == query_window.size()) {
                        if (node == target.get_path().back())
                            return cost;

                        common::logger->info("\t\twrong node, skipping");

                        continue;
                    }

                    // forward creation of insertions
                    if (query_dist < query_window.size()) {
                        // extend a previous insertion
                        if (cost < E.size() && query_dist < E[cost].size()) {
                            auto find = E[cost][query_dist].find(node);
                            if (find != E[cost][query_dist].end()) {
                                auto [last_dist, last_num_ops] = find->second;
                                ssize_t next_ext_cost = cost + gap_ext;
                                size_t stored_dist = set_value(E, next_ext_cost, query_dist + 1, node, last_dist, node, last_num_ops + 1, Cigar::INSERTION);
                                set_value(S, next_ext_cost, query_dist + 1, node, stored_dist, node, 0, Cigar::INSERTION);
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }

                        // open an insertion
                        ssize_t next_opn_cost = cost + gap_opn;
                        size_t stored_dist = set_value(E, next_opn_cost, query_dist + 1, node, best_dist, node, 1, Cigar::INSERTION);
                        set_value(S, next_opn_cost, query_dist + 1, node, stored_dist, node, 0, Cigar::INSERTION);
                        it = S[cost][query_dist].begin() + it_dist;
                    }

                    if (best_dist < dist) {
                        // deletion
                        std::vector<std::pair<node_index, char>> prevs;
                        graph.call_incoming_kmers(node, [&](auto prev, char c) {
                            prevs.emplace_back(prev, c);
                        });
                        if (prevs.size()) {
                            if (cost < F.size() && query_dist < F[cost].size()) {
                                // extension
                                auto find = F[cost][query_dist].find(node);
                                if (find != F[cost][query_dist].end()) {
                                    auto [last_dist, last_num_ops, last_node] = find->second;
                                    assert(last_dist >= last_num_ops);
                                    ssize_t next_ext_cost = cost + gap_ext;
                                    for (const auto &[prev, c] : prevs) {
                                        size_t stored_value = set_value(F, next_ext_cost, query_dist, prev, last_dist + 1, node, last_num_ops + 1, Cigar::DELETION);
                                        set_value(S, next_ext_cost, query_dist, prev, stored_value, node, 0, Cigar::DELETION);
                                        it = S[cost][query_dist].begin() + it_dist;
                                    }
                                }
                            }

                            ssize_t next_opn_cost = cost + gap_opn;
                            for (const auto &[prev, c] : prevs) {
                                size_t stored_value = set_value(F, next_opn_cost, query_dist, prev, best_dist + 1, node, 1, Cigar::DELETION);
                                set_value(S, next_opn_cost, query_dist, prev, stored_value, node, 0, Cigar::DELETION);
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }

                        if (query_dist < query_window.size()) {
                            // match
                            std::string_view local_query_window = query_window;
                            local_query_window.remove_suffix(query_dist);

                            for (auto [prev, c] : prevs) {
                                std::cerr << prev;
                                size_t cur_best = best_dist + 1;
                                size_t cur_query_dist = query_dist + 1;
                                size_t cur_cost = cost;
                                Cigar::Operator op = Cigar::MATCH;
                                if (c != local_query_window.back()) {
                                    cur_cost += mismatch_cost;
                                    op = Cigar::MISMATCH;
                                }

                                for (auto jt = local_query_window.rbegin() + 1; jt != local_query_window.rend(); ++jt) {
                                    if (cur_best == dist || !graph.has_single_incoming(prev))
                                        break;

                                    if (auto pprev = graph.traverse_back(prev, *jt)) {
                                        std::cerr << "\t" << pprev;
                                        ++cur_best;
                                        ++cur_query_dist;
                                        prev = pprev;
                                    } else {
                                        break;
                                    }
                                }
                                std::cerr << "\n";

                                common::logger->info("\t\tM: {} -> {}\t{} -> {}\t{}", cost, cur_cost, query_dist, cur_query_dist, Cigar::opt_to_char(op));
                                assert(cur_query_dist > query_dist || cur_cost > cost);
                                set_value(
                                    S,
                                    cur_cost,
                                    cur_query_dist,
                                    prev,
                                    cur_best,
                                    node,
                                    cur_query_dist - query_dist,
                                    op,
                                    c
                                );
                                it = S[cost][query_dist].begin() + it_dist;
                            }
                        }
                    }
                }
            }
        }

        return std::numeric_limits<size_t>::max();
    };

    size_t cost = fill_table();
    if (cost == std::numeric_limits<size_t>::max())
        throw std::runtime_error("No alignment found");

    common::logger->info("FOUND ALIGNMENT AT COST {}", cost);

#ifndef NDEBUG
    for (const auto &wf : S) {
        for (const auto &node_table : wf) {
            for (const auto &[node, data] : node_table) {
                // std::tuple<size_t, size_t, node_index, Cigar::Operator, char>
                const auto &[cur_dist, num_match, last_node, last_op, mismatch_char] = data;
                common::logger->info("S: {}\t{}\t{}\t{}\t{}", cur_dist, num_match, last_node, Cigar::opt_to_char(last_op), mismatch_char);
                assert(last_op == Cigar::MATCH
                    || last_op == Cigar::MISMATCH
                    || last_op == Cigar::INSERTION
                    || last_op == Cigar::DELETION);
            }
        }
    }
#endif

    // backtrack
    common::logger->info("Backtracking");
    node_index node = target.get_path().back();
    size_t query_dist = query_window.size();

    assert(cost < S.size() && query_dist < S[cost].size());
    auto bt = S[cost][query_dist].find(node);
    assert(bt != S[cost][query_dist].end());

    std::vector<node_index> path(target.get_path().begin(), target.get_path().end() - 1);
    Cigar cigar(Cigar::CLIPPED, target.get_clipping());
    cigar.append(Cigar::MATCH, path.size());

    do {
        assert(std::get<0>(bt->second) == dist);
        auto last_bt = bt;

        // check insertion
        if (std::get<3>(bt->second) == Cigar::INSERTION) {
            assert(cost > 0 && cost < E.size() && query_dist < E[cost].size() && query_dist > 0);
            assert(std::get<1>(bt->second) == 0);
            auto et = E[cost][query_dist].find(node);
            assert(et != E[cost][query_dist].end());
            assert(std::get<0>(et->second) == dist);

            auto [ins_dist, num_ins] = et->second;
            assert(query_dist >= num_ins);

            size_t prev_cost = cost - gap_opn - (num_ins - 1) * gap_ext;
            size_t prev_query_dist = query_dist - num_ins;
            assert(prev_query_dist < S[prev_cost].size());

            auto check_bt = S[prev_cost][prev_query_dist].find(node);
            assert(check_bt != S[prev_cost][prev_query_dist].end());
            assert(std::get<0>(check_bt->second) == dist);

            query_dist = prev_query_dist;
            cost = prev_cost;
            bt = check_bt;
            cigar.append(Cigar::INSERTION, num_ins);
        }

        // check deletion
        if (std::get<3>(bt->second) == Cigar::DELETION) {
            assert(cost > 0 && cost < F.size() && query_dist < F[cost].size() && dist > 0);
            assert(std::get<1>(bt->second) == 0);
            auto ft = F[cost][query_dist].find(node);
            assert(ft != F[cost][query_dist].end());

            auto [del_dist, num_del, prev_node] = ft->second;
            common::logger->info("D: {}\t{}\t{}", del_dist, num_del, prev_node);
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

                assert(query_dist < F[prev_cost].size());
                auto prev_ft = F[prev_cost][query_dist].find(node);
                assert(prev_ft != F[prev_cost][query_dist].end());
                assert(std::get<0>(prev_ft->second) == prev_dist);
                assert(std::get<1>(prev_ft->second) == num_del);
                ft = prev_ft;
                prev_node = std::get<2>(ft->second);
            }

            path.emplace_back(node);

            prev_cost -= gap_opn;
            --prev_dist;
            assert(query_dist < S[prev_cost].size());

            auto check_bt = S[prev_cost][query_dist].find(prev_node);
            assert(check_bt != S[prev_cost][query_dist].end());
            assert(std::get<0>(check_bt->second) == prev_dist);

            dist = prev_dist;
            cost = prev_cost;
            node = prev_node;
            bt = check_bt;
        }

        // check match
        if (std::get<3>(bt->second) == Cigar::MATCH || std::get<3>(bt->second) == Cigar::MISMATCH) {
            auto [cur_dist, num_match, last_node, last_op, mismatch_char] = bt->second;
            common::logger->info("M: {}\t{}\t{}\t{}\t{}", cur_dist, num_match, last_node, Cigar::opt_to_char(last_op), mismatch_char);

            if (!num_match) {
                assert(last_op == Cigar::MATCH);
                bt = S[cost][query_dist].end();
                continue;
            }

            assert(query_dist >= num_match);
            assert(dist >= num_match);

            node_index traverse_node = last_node;

            cigar.append(Cigar::MATCH, num_match - 1);
            cigar.append(last_op);
            size_t prev_dist = dist - num_match;
            size_t prev_query = query_dist - num_match;

            // reconstruct the backwards traversal from last_node to node
            auto it = query_window.rbegin() + prev_query;
            assert(it != query_window.rend());
            size_t cur_path_size = path.size();
            if (std::get<3>(bt->second) == Cigar::MISMATCH) {
                assert(mismatch_char != *it);
                traverse_node = graph.traverse_back(traverse_node, mismatch_char);
                assert(traverse_node != DeBruijnGraph::npos);
                std::cerr << "\t" << traverse_node;
                path.emplace_back(traverse_node);
                ++it;
                --num_match;
            } else {
                assert(mismatch_char == *it);
            }

            while (num_match) {
                assert(it != query_window.rend());
                traverse_node = graph.traverse_back(traverse_node, *it);
                std::cerr << "\t" << traverse_node;
                assert(traverse_node != DeBruijnGraph::npos);
                path.emplace_back(traverse_node);
                ++it;
                --num_match;
            }
            std::cerr << std::endl;
            assert(traverse_node == node);

            std::reverse(path.begin() + cur_path_size, path.end());

            size_t prev_cost = cost - (last_op == Cigar::MATCH ? 0 : mismatch_cost);
            assert(prev_query < S[prev_cost].size());

            auto check_bt = S[prev_cost][prev_query].end();
            if (prev_cost || prev_query || prev_dist) {
                check_bt = S[prev_cost][prev_query].find(last_node);
                assert(check_bt != S[prev_cost][prev_query].end());
                assert(std::get<0>(check_bt->second) == prev_dist);
            } else {
                assert(last_node == start_node);
            }

            dist = prev_dist;
            cost = prev_cost;
            query_dist = prev_query;
            node = last_node;
            bt = check_bt;
        }

        if (bt == last_bt)
            throw std::runtime_error("Inf. loop");

    } while (bt != S[cost][query_dist].end());

    assert(dist == 0);
    assert(cost == 0);
    assert(query_dist == 0);
    assert(path.size());
    callback(std::move(path), std::move(cigar));
}

std::vector<Alignment> ExactSeeder::get_inexact_anchors() const {
    std::vector<Anchor> anchors = get_anchors();
    std::vector<Alignment> alignments;

    using AnchorIt = std::vector<Anchor>::iterator;

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

                std::cerr << "Trying to connect "
                          << a_i << " -> " << a_j;

                if (dist <= 0 || query_i.begin() >= query_j.begin()) {
                    std::cerr << "\tbadorder\n";
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
                    std::cerr << "\t" << base_score << " <= " << score_j << "\n";
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
                    std::cerr << "\tnopath\n";
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

                std::cerr << "\td: " << dist << "," << "cd:" << coord_dist
                          << "\tscore: " << std::get<0>(*chain_scores) << "->" << base_score << " -> " << score << " vs. " << score_j << std::endl;
                update_score(score, it, coord_dist);

                ++chain_scores;
            }
        },
        [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
            for (const auto &[it, dist] : chain) {
                std::cerr << "\t" << Alignment(*it) << "->" << dist;
            }
            std::cerr << std::endl;
            return true;
        },
        [this](AnchorIt last,
               AnchorIt next,
               Alignment&& aln,
               size_t last_to_next_dist,
               DBGAlignerConfig::score_t score_up_to_now,
               const AlignmentCallback &callback) {
            size_t next_to_last_dist = last_to_next_dist - last->get_spelling().size() + next->get_spelling().size();
            std::cerr << "Connecting " << Alignment(*next) << " -> " << aln << "\t" << last_to_next_dist << "\n";
            std::cerr << "\ti.e., " << Alignment(*next) << " <- " << aln << "\t" << next_to_last_dist << "\n";
            align(query_.get_graph(), config_, aln.get_path()[0], aln, *next, next_to_last_dist, [&](auto&& path, auto&& cigar) {
                path.insert(path.end(), aln.get_path().begin(), aln.get_path().end());

                Cigar aln_cigar = aln.get_cigar();
                aln_cigar.trim_clipping();
                cigar.append(std::move(aln_cigar));

                common::logger->info("Old path: {}", fmt::join(aln.get_path(), ","));
                common::logger->info("Cigar: {}\tPath: {}", cigar.to_string(), fmt::join(path, ","));

                callback(Alignment(query_.get_graph(),
                                   aln.get_query(),
                                   aln.get_orientation(),
                                   std::move(path),
                                   config_,
                                   std::move(cigar),
                                   aln.get_end_trim()));
            });
        },
        [&alignments](Alignment&& alignment) {
            std::cerr << "Final ALN\t" << alignment << "\n";
            alignments.emplace_back(std::move(alignment));
        }
    );

    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_score(), a.get_seed().size())
             > std::make_tuple(b.get_score(), b.get_seed().size());
    });
    return alignments;
}

} // namespace mtg::graph::align