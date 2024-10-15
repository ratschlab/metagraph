#include "aln_seeder.hpp"

#include "aln_chainer.hpp"
#include "common/logger.hpp"
#include "common/vector_map.hpp"

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
            size_t begin = 0;
            graph.map_to_nodes_sequentially(query_.get_query(orientation),
                [&](Match::node_index node) {
                    if (node != DeBruijnGraph::npos) {
                        anchors.emplace_back(query_.get_query(),
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
            break;

        ++coord_dist;
        graph.adjacent_outgoing_nodes(node, [&](DeBruijnGraph::node_index next) {
            traversal_stack.emplace_back(next, coord_dist);
        });
    }
}

template <typename T>
struct template_parameter;

template <template <typename ...> class C, typename T>
struct template_parameter<C<T>> { using type = T; };

template <typename Map>
using WaveFront = std::vector<Map>; // S(cost)

template <typename Map>
using ScoreTable = std::vector<WaveFront<Map>>;

template <typename T>
using get_table_t = typename template_parameter<std::decay_t<T>>::type;

void global_align(const DeBruijnGraph &graph,
                  const DBGAlignerConfig &config,
                  const Alignment &aln,
                  const Anchor &target,
                  size_t dist,
                  const AlignmentCallback &callback) {
    assert(dist > 0);

    using node_index = DeBruijnGraph::node_index;

    std::string_view query_window = aln.get_query();
    query_window.remove_suffix(aln.get_end_clipping() + aln.get_seed().size());
    query_window.remove_prefix(target.get_clipping());
    if (query_window.empty())
        return;

    // score = 1/2(match_score(query_size + match_spelling_size - cost))
    // => 2*score/match_score - query_size - match_spelling_size = -cost
    // => cost = query_size + match_spelling_size - 2*score/match_score

    DBGAlignerConfig::score_t match_score = config.match_score("A");
    auto score_to_cost = [&,match_score](DBGAlignerConfig::score_t score, ssize_t query=0, ssize_t match=0) -> ssize_t {
        query = query_window.size();
        match = dist;
        return query + match - 2 * (score - aln.get_score()) / match_score;
    };

    auto cost_to_score = [&,match_score](ssize_t cost, ssize_t query=0, ssize_t match=0) -> DBGAlignerConfig::score_t {
        query = query_window.size();
        match = dist;
        return match_score * (query + match - cost) / 2 + aln.get_score();
    };

    size_t mismatch_cost = score_to_cost(-match_score);
    size_t gap_opn = score_to_cost(config.gap_opening_penalty);
    size_t gap_ext = score_to_cost(config.gap_extension_penalty);

    // S(cost, query_dist, node) = best_dist
    // using NodeToBestDists = VectorMap<DeBruijnGraph::node_index, size_t>; // S(cost, query_dist)
    using SMap = VectorMap<node_index, std::tuple<size_t, size_t, node_index, Cigar::Operator, char>>;
    using EMap = VectorMap<node_index, std::tuple<size_t, size_t>>;
    using FMap = VectorMap<node_index, std::tuple<size_t, size_t, node_index>>;

    ScoreTable<SMap> S; // best cost
    ScoreTable<EMap> E; // best cost (last operation is insertion)
    ScoreTable<FMap> F; // best cost (last operation is deletion)

    auto fill_table = [&]() -> size_t {
        size_t query_dist = 0;
        auto node = aln.get_path()[0];
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
        auto &wf = S.emplace_back();
        wf.resize(query_dist + 1);
        wf[query_dist][node] = std::make_tuple(best_dist, best_dist, DeBruijnGraph::npos, Cigar::MATCH, '\0');
        if (best_dist == dist) {
            if (node == target.get_path().back())
                return 0;

            return std::numeric_limits<size_t>::max();
        }

        for (size_t cost = 0; cost < S.size(); ++cost) {
            if (S[cost].empty())
                continue;

            for (size_t query_dist = 0; query_dist < S[cost].size(); ++query_dist) {
                if (S[cost][query_dist].empty())
                    continue;

                size_t it_dist = 0;
                for (auto it = S[cost][query_dist].begin(); it != S[cost][query_dist].end(); ++it, ++it_dist) {
                    assert(query_dist <= query_window.size());

                    auto set_value = [&](auto &table, size_t cost, size_t query_dist, node_index node, size_t dist, node_index last_node, size_t num_ops, Cigar::Operator last_op, char c = '\0') -> size_t {
                        if (cost >= table.size()) {
                            table.resize(cost + 1);
                            if constexpr(std::is_same_v<get_table_t<decltype(table)>, SMap>)
                                it = S[cost][query_dist].begin() + it_dist;
                        }

                        if (query_dist >= table[cost].size()) {
                            table[cost].resize(query_dist + 1);
                            if constexpr(std::is_same_v<get_table_t<decltype(table)>, SMap>)
                                it = S[cost][query_dist].begin() + it_dist;
                        }

                        auto &bucket = table[cost][query_dist][node];
                        if constexpr(std::is_same_v<get_table_t<decltype(table)>, SMap>)
                            it = S[cost][query_dist].begin() + it_dist;

                        if (dist > std::get<0>(bucket)) {
                            std::get<0>(bucket) = dist;
                            std::get<1>(bucket) = num_ops;

                            if constexpr(std::is_same_v<get_table_t<decltype(table)>, FMap>) {
                                std::get<2>(bucket) = last_node;
                            }

                            if constexpr(std::is_same_v<get_table_t<decltype(table)>, SMap>) {
                                std::get<2>(bucket) = last_node;
                                std::get<3>(bucket) = last_op;
                                std::get<4>(bucket) = c;
                            }
                        }

                        return std::get<0>(bucket);
                    };

                    node_index node = it->first;
                    auto [best_dist, last_ext, last_node, last_op, mismatch_char] = it->second;

                    if (best_dist == dist) {
                        if (node == target.get_path().back() && query_dist == query_window.size())
                            return cost;

                        continue;
                    }

                    // DBGAlignerConfig::score_t prev_score = cost_to_score(cost, query_dist, best_dist);

                    // forward creation of insertions
                    if (query_dist + 1 <= query_window.size()) {
                        // extend a previous insertion
                        if (cost < E.size() && E[cost].size() && query_dist < E[cost].size()) {
                            auto find = E[cost][query_dist].find(node);
                            if (find != E[cost][query_dist].end()) {
                                auto [last_dist, last_num_ops] = find->second;
                                ssize_t next_ext_cost = cost + gap_ext;
                                // ssize_t next_ext_cost = score_to_cost(prev_score + config.gap_extension_penalty, query_dist + 1, last_dist);
                                size_t stored_dist = set_value(E, next_ext_cost, query_dist + 1, node, last_dist, node, last_num_ops + 1, Cigar::INSERTION);
                                set_value(S, next_ext_cost, query_dist + 1, node, stored_dist, node, 0, Cigar::INSERTION);
                            }
                        }

                        // open an insertion
                        ssize_t next_opn_cost = cost + gap_opn;
                        // ssize_t next_opn_cost = score_to_cost(prev_score + config.gap_opening_penalty, query_dist + 1, best_dist);
                        size_t stored_dist = set_value(E, next_opn_cost, query_dist + 1, node, best_dist, node, 1, Cigar::INSERTION);
                        set_value(S, next_opn_cost, query_dist + 1, node, stored_dist, node, 0, Cigar::INSERTION);
                    }

                    if (best_dist < dist) {
                        if (query_dist < query_window.size()) {
                            // match
                            std::string_view local_query_window = query_window;
                            local_query_window.remove_suffix(query_dist);

                            // const auto &score_row = config.score_matrix[local_query_window.back()];
                            graph.call_incoming_kmers(node, [&](node_index prev, char c) {
                                size_t cur_best = best_dist + 1;
                                size_t cur_query_dist = query_dist + 1;
                                size_t cur_cost = c == local_query_window.back()
                                    ? cost
                                    : cost + mismatch_cost;
                                    // : score_to_cost(prev_score + score_row[c], cur_query_dist, cur_best);
                                for (auto jt = local_query_window.rbegin() + 1; jt != local_query_window.rend(); ++jt) {
                                    if (cur_best == dist || !graph.has_single_incoming(prev))
                                        break;

                                    if (auto pprev = graph.traverse_back(prev, *jt)) {
                                        ++cur_best;
                                        ++cur_query_dist;
                                        prev = pprev;
                                    } else {
                                        break;
                                    }
                                }

                                assert(cur_query_dist > query_dist || cur_cost > cost);
                                set_value(S, cur_cost, cur_query_dist, prev, cur_best, node, cur_query_dist - query_dist, c == local_query_window.back() ? Cigar::MATCH : Cigar::MISMATCH, c);
                            });

                        } else {
                            // deletion
                            std::vector<node_index> prevs;
                            graph.adjacent_incoming_nodes(node, [&](auto prev) {
                                prevs.emplace_back(prev);
                            });
                            if (prevs.size()) {
                                if (cost < F.size() && F[cost].size() && query_dist < F[cost].size()) {
                                    // extension
                                    auto find = F[cost][query_dist].find(node);
                                    if (find != F[cost][query_dist].end()) {
                                        auto [last_dist, last_num_ops, last_node] = find->second;
                                        ssize_t next_ext_cost = cost + gap_ext;
                                        // ssize_t next_ext_cost = score_to_cost(prev_score + config.gap_extension_penalty, query_dist, last_dist + 1);
                                        for (auto prev : prevs) {
                                            size_t stored_value = set_value(F, next_ext_cost, query_dist, prev, last_dist + 1, node, last_num_ops + 1, Cigar::DELETION);
                                            set_value(S, next_ext_cost, query_dist, prev, stored_value, node, 0, Cigar::DELETION);
                                        }
                                    }
                                }

                                ssize_t next_opn_cost = cost + gap_opn;
                                // ssize_t next_opn_cost = score_to_cost(prev_score + config.gap_opening_penalty, query_dist, best_dist + 1);
                                for (auto prev : prevs) {
                                    size_t stored_value = set_value(F, next_opn_cost, query_dist, prev, best_dist + 1, node, 1, Cigar::DELETION);
                                    set_value(S, next_opn_cost, query_dist, prev, stored_value, node, 0, Cigar::DELETION);
                                }
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

    common::logger->info("FOUND ALIGNMENT AT COST {}\tscore {}", cost, cost_to_score(cost));

    // backtrack
    node_index node = target.get_path().back();
    size_t query_dist = query_window.size();

    assert(cost < S.size() && query_dist < S[cost].size());
    auto bt = S[cost][query_dist].find(node);
    assert(bt != S[cost][query_dist].end());

    std::vector<node_index> path;
    Cigar cigar(Cigar::CLIPPED, target.get_clipping());

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
            assert(num_ins >= query_dist);

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
            assert(std::get<0>(ft->second) == dist);

            auto [del_dist, num_del, prev_node] = ft->second;
            assert(num_del >= dist);

            size_t prev_cost = cost;
            size_t prev_dist = dist;
            while (num_del > 1) {
                prev_cost -= gap_ext;
                --prev_dist;

                path.emplace_back(node);
                cigar.append(Cigar::DELETION);
                node = prev_node;
                --num_del;

                assert(query_dist < F[prev_cost].size());
                auto prev_ft = F[prev_cost][query_dist].find(node);
                assert(ft != F[cost][query_dist].end());
                assert(std::get<0>(ft->second) == prev_dist);
                assert(std::get<1>(ft->second) == num_del);
                ft = prev_ft;
                prev_node = std::get<2>(ft->second);
            }

            prev_cost -= gap_opn;
            cigar.append(Cigar::DELETION);

            --prev_dist;
            assert(query_dist < S[prev_cost].size());

            auto check_bt = S[prev_cost][query_dist].find(prev_node);
            assert(check_bt != S[prev_cost][query_dist].end());
            assert(std::get<0>(check_bt->second) == prev_dist);

            dist = prev_dist;
            cost = prev_cost;
            bt = check_bt;
        }

        // check match
        if (std::get<3>(bt->second) == Cigar::MATCH || std::get<3>(bt->second) == Cigar::MISMATCH) {
            auto [cur_dist, num_match, last_node, last_op, mismatch_char] = bt->second;
            common::logger->info("M: {}\t{}\t{}\t{}\t{}", cur_dist, num_match, last_node, Cigar::opt_to_char(last_op), mismatch_char);
            assert(num_match > 0);
            assert(query_dist >= num_match);
            assert(dist >= num_match);

            node_index traverse_node = last_node;
            if (traverse_node == DeBruijnGraph::npos)
                traverse_node = aln.get_path()[0];

            cigar.append(Cigar::MATCH, num_match - 1);
            cigar.append(last_op);
            size_t prev_dist = dist - num_match;
            size_t prev_query = query_dist - num_match;

            // reconstruct the backwards traversal from last_node to node
            auto it = query_window.rbegin() + prev_query;
            size_t cur_path_size = path.size();
            if (std::get<3>(bt->second) == Cigar::MISMATCH) {
                assert(mismatch_char != '\0');
                assert(it != query_window.rend());
                assert(mismatch_char != *it);
                traverse_node = graph.traverse_back(traverse_node, mismatch_char);
                path.emplace_back(traverse_node);
                ++it;
                --num_match;
            }

            while (num_match) {
                assert(it != query_window.rend());
                traverse_node = graph.traverse_back(traverse_node, *it);
                assert(traverse_node != DeBruijnGraph::npos);
                path.emplace_back(traverse_node);
                ++it;
                --num_match;
            }
            assert(traverse_node == node);

            std::reverse(path.begin() + cur_path_size, path.end());

            size_t prev_cost = cost - (last_op == Cigar::MATCH ? 0 : mismatch_cost);
            assert(prev_query < S[prev_cost].size());

            auto check_bt = S[prev_cost][prev_query].end();
            if (last_node != DeBruijnGraph::npos) {
                check_bt = S[prev_cost][prev_query].find(last_node);
                assert(check_bt != S[prev_cost][prev_query].end());
                assert(std::get<0>(check_bt->second) == prev_dist);
            }

            dist = prev_dist;
            cost = prev_cost;
            query_dist = prev_query;
            bt = check_bt;
        }

        if (bt == last_bt)
            throw std::runtime_error("Inf. loop");

    } while (bt != S[cost][query_dist].end());

    assert(dist == 0);
    assert(cost == 0);
    assert(query_dist == 0);
    assert(path.size());
    path.insert(path.end(), aln.get_path().begin(), aln.get_path().end());

    Cigar aln_cigar = aln.get_cigar();
    aln_cigar.trim_clipping();
    cigar.append(std::move(aln_cigar));

    common::logger->info("Old path: {}", fmt::join(aln.get_path(), ","));
    common::logger->info("Cigar: {}\tPath: {}", cigar.to_string(), fmt::join(path, ","));

    callback(Alignment(graph,
                       aln.get_query(),
                       aln.get_orientation(),
                       std::move(path),
                       config,
                       std::move(cigar),
                       aln.get_end_trim()));
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
                size_t front_dist = a_j.get_clipping() - a_i.get_clipping();
                if (front_dist < a_i.get_path().size()) {
                    if (a_i.get_path()[front_dist] != a_j.get_path()[0]) {
                        ++chain_scores;
                        return;
                    }
                } else if (front_dist < a_i.get_path().size() + graph.get_k() - 1) {
                    std::string path_spelling_i = a_i.get_path_spelling();
                    std::string path_spelling_j = a_j.get_path_spelling();
                    auto [i_mismatch, j_mismatch] = std::mismatch(
                        path_spelling_i.begin() + front_dist, path_spelling_i.end(),
                        path_spelling_j.c_str(), path_spelling_j.c_str() + path_spelling_j.size()
                    );

                    if (i_mismatch != path_spelling_i.end()) {
                        ++chain_scores;
                        return;
                    }

                    if (path_spelling_i.size() - front_dist < graph.get_k() - 1) {
                        // not enough overlap, we need to see if a path exists
                        auto j_check_end = path_spelling_j.c_str() + graph.get_k();
                        assert(j_mismatch < j_check_end);
                        Match::node_index last_node = DeBruijnGraph::npos;
                        graph.traverse(a_i.get_path().back(),
                                       j_mismatch, j_check_end,
                                       [&](Match::node_index node) {
                                           last_node = node;
                                           ++j_mismatch;
                                       });
                        if (j_mismatch != j_check_end || last_node != a_j.get_path()[0]) {
                            ++chain_scores;
                            return;
                        }

                        assert(last_node != DeBruijnGraph::npos);
                    }
                } else {
                    // completely disjoint, we need to explore freely
                    coord_dist = -DBGAlignerConfig::ninf;
                    ssize_t min_diff = std::numeric_limits<ssize_t>::max();
                    call_distances(graph, a_i.get_path().back(), a_j.get_path()[0],
                        [&](ssize_t path_dist) {
                            DBGAlignerConfig::score_t cur_coord_dist
                                = a_i.get_end_trim() + path_dist + a_j.get_path().size() - 1 - a_j.get_end_trim();
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

                    if (coord_dist == -DBGAlignerConfig::ninf) {
                        ++chain_scores;
                        return;
                    }
                }

                DBGAlignerConfig::score_t diff = std::abs(coord_dist - dist);
                score += config_.gap_opening_penalty
                        + (diff - 1) * config_.gap_extension_penalty;

                update_score(score, it, coord_dist);

                ++chain_scores;
            }
        },
        [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
            return true;
        },
        [this](AnchorIt last,
               AnchorIt next,
               Alignment&& aln,
               size_t last_to_next_dist,
               DBGAlignerConfig::score_t score_up_to_now,
               const AlignmentCallback &callback) {
            std::cerr << "Connecting " << Alignment(*next) << " -> " << aln << "\n";
            size_t next_to_last_dist = last_to_next_dist - next->get_seed().size() + last->get_seed().size() - last->get_path().size() + 1;
            global_align(query_.get_graph(), config_, aln, *next, next_to_last_dist, callback);
        },
        [&alignments](Alignment&& alignment) { alignments.emplace_back(std::move(alignment)); }
    );

    return alignments;
}

} // namespace mtg::graph::align