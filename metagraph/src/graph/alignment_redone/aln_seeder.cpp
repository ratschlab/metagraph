#include "aln_seeder.hpp"

#include "aln_chainer.hpp"
#include "common/logger.hpp"
#include "common/vector_map.hpp"

namespace mtg::graph::align {

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

void global_align(const DeBruijnGraph &graph,
                  const DBGAlignerConfig &config,
                  const Alignment &aln,
                  const Anchor &target,
                  size_t dist,
                  const AlignmentCallback &callback) {
    assert(dist > 0);

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
    using NodeToBestDists = VectorMap<DeBruijnGraph::node_index, size_t>; // S(cost, query_dist)
    using WaveFront = std::vector<NodeToBestDists>; // S(cost)
    using ScoreTable = std::vector<WaveFront>;

    ScoreTable S; // best cost
    ScoreTable E; // best cost (last operation is insertion)
    ScoreTable F; // best cost (last operation is deletion)

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
        wf[query_dist][node] = best_dist;
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

                    auto set_value = [&](auto &table, size_t cost, size_t query_dist, auto node, size_t dist) -> size_t {
                        if (cost >= table.size()) {
                            table.resize(cost + 1);
                            if (&table == &S)
                                it = S[cost][query_dist].begin() + it_dist;
                        }

                        if (query_dist >= table[cost].size()) {
                            table[cost].resize(query_dist + 1);
                            if (&table == &S)
                                it = S[cost][query_dist].begin() + it_dist;
                        }

                        auto &bucket = table[cost][query_dist][node];
                        if (&table == &S)
                            it = S[cost][query_dist].begin() + it_dist;

                        bucket = std::max(bucket, dist);

                        return bucket;
                    };

                    DeBruijnGraph::node_index node = it->first;
                    size_t best_dist = it->second;

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
                                size_t last_dist = find->second;
                                ssize_t next_ext_cost = cost + gap_ext;
                                // ssize_t next_ext_cost = score_to_cost(prev_score + config.gap_extension_penalty, query_dist + 1, last_dist);
                                size_t stored_dist = set_value(E, next_ext_cost, query_dist + 1, node, last_dist);
                                set_value(S, next_ext_cost, query_dist + 1, node, stored_dist);
                            }
                        }

                        // open an insertion
                        ssize_t next_opn_cost = cost + gap_opn;
                        // ssize_t next_opn_cost = score_to_cost(prev_score + config.gap_opening_penalty, query_dist + 1, best_dist);
                        size_t stored_dist = set_value(E, next_opn_cost, query_dist + 1, node, best_dist);
                        set_value(S, next_opn_cost, query_dist + 1, node, stored_dist);
                    }

                    if (best_dist < dist) {
                        if (query_dist < query_window.size()) {
                            // match
                            std::string_view local_query_window = query_window;
                            local_query_window.remove_suffix(query_dist);

                            // const auto &score_row = config.score_matrix[local_query_window.back()];
                            graph.call_incoming_kmers(node, [&](DeBruijnGraph::node_index prev, char c) {
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
                                set_value(S, cur_cost, cur_query_dist, prev, cur_best);
                            });

                        } else {
                            // deletion
                            std::vector<DeBruijnGraph::node_index> prevs;
                            graph.adjacent_incoming_nodes(node, [&](auto prev) {
                                prevs.emplace_back(prev);
                            });
                            if (prevs.size()) {
                                if (cost < F.size() && F[cost].size() && query_dist < F[cost].size()) {
                                    // extension
                                    auto find = F[cost][query_dist].find(node);
                                    if (find != F[cost][query_dist].end()) {
                                        size_t last_dist = find->second;
                                        ssize_t next_ext_cost = cost + gap_ext;
                                        // ssize_t next_ext_cost = score_to_cost(prev_score + config.gap_extension_penalty, query_dist, last_dist + 1);
                                        for (auto prev : prevs) {
                                            size_t stored_value = set_value(F, next_ext_cost, query_dist, prev, last_dist + 1);
                                            set_value(S, next_ext_cost, query_dist, prev, stored_value);
                                        }
                                    }
                                }

                                ssize_t next_opn_cost = cost + gap_opn;
                                // ssize_t next_opn_cost = score_to_cost(prev_score + config.gap_opening_penalty, query_dist, best_dist + 1);
                                for (auto prev : prevs) {
                                    size_t stored_value = set_value(F, next_opn_cost, query_dist, prev, best_dist + 1);
                                    set_value(S, next_opn_cost, query_dist, prev, stored_value);
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

    common::logger->trace("FOUND ALIGNMENT AT COST {}", cost);

    // backtrack
    DeBruijnGraph::node_index node = target.get_path().back();
    size_t query_dist = query_window.size();

    assert(cost < S.size() && query_dist < S[cost].size());
    auto bt = S[cost][query_dist].find(node);
    assert(bt != S[cost][query_dist].end());
    assert(bt->second == dist);

//     while (bt != S[cost][query_dist].end()) {
//         assert(bt->second == dist);
//         if (cost < E.size() && query_dist < E[cost].size()) {
//             // check insertion
//             auto et = E[cost][query_dist].find(node);
//             if (et != E[cost][query_dist].end() && et->second == dist) {
//                 // S == E, so we have an insertion. check if it's an extension or an insertion

//                 // backtrack through all
//                 bool found = true;
//                 while (found) {
//                     found = false;
//                     size_t prev_cost = cost - gap_ext;
//                     if (query_dist - 1 < E[prev_cost].size()) {
//                         auto prev_et = E[prev_cost][query_dist - 1].find(node);
//                         if (prev_et != E[prev_cost][query_dist - 1].end() && prev_et->second == dist) {
//                             cost = prev_cost;
//                             --query_dist;
//                             et = prev_et;
//                             found = true;
//                         }
//                     }
//                 }

//                 size_t prev_cost = cost - gap_opn;
//                 if (query_dist - 1 < E[prev_cost].size()) {

//                 } else {

//                 }
//             }
//         }

//         if (cost < F.size() && query_dist < F[cost].size()) {
//             // check deletion
//             auto ft = F[cost][query_dist].find(node);
//             if (ft != F[cost][query_dist].end() && ft->second == dist) {
//                 // S == F, so we may have a deletion
//                 // TODO
//             }
//         }

//         if (query_dist > 0 && dist > 0) {
//             // check match
//             DeBruijnGraph::node_index found_prev = DeBruijnGraph::npos;
//             size_t cur_query_dist = query_dist;
//             for (auto it = S[cost].rbegin(); it != S[cost].rend(); ++it, --cur_query_dist) {
//                 for (const auto &[prev_node, prev_dist] : *it) {
//                     if (query_dist > cur_query_dist && dist > prev_dist && dist - prev_dist == query_dist - cur_query_dist) {
//                         assert(found_prev == DeBruijnGraph::npos);
//                         found_prev = prev_node;
// #ifdef NDEBUG
//                         break;
// #endif
//                     }
//                 }
// #ifdef NDEBUG
//                 if (found_prev != DeBruijnGraph::npos)
//                     break;
// #endif
//             }

//             if (found_prev)
//         }
//     }
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
            assert(a_i.get_label_class() == a_j.get_label_class());
            assert(a_i.get_orientation() == a_j.get_orientation());

            const DeBruijnGraph &graph = query_.get_graph();
            std::string_view query_j = a_j.get_seed();
            const DBGAlignerConfig::score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));

            for (auto it = begin; it != end; ++it) {
                const Anchor &a_i = *it;
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
        [this](AnchorIt next,
               AnchorIt last,
               Alignment&& aln,
               size_t last_to_next_dist,
               DBGAlignerConfig::score_t score_up_to_now,
               const AlignmentCallback &callback) {
            size_t next_to_last_dist = last_to_next_dist - next->get_seed().size() + last->get_seed().size() - last->get_path().size() + 1;
            global_align(query_.get_graph(), config_, aln, *next, next_to_last_dist, callback);
        },
        [&alignments](Alignment&& alignment) { alignments.emplace_back(std::move(alignment)); }
    );

    return alignments;
}

} // namespace mtg::graph::align