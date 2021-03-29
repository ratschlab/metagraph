#include "aligner_labeled.hpp"

#include "graph/annotated_dbg.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/hashers/hash.hpp"
#include "common/vectors/vector_algorithm.hpp"

namespace mtg {
namespace graph {
namespace align {


bool check_targets(const AnnotatedDBG &anno_graph,
                   const Alignment<DeBruijnGraph::node_index> &path) {
    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<uint64_t> ref_targets;
    for (const std::string &label : anno_graph.get_labels(path.get_sequence(), 1.0)) {
        ref_targets.emplace(label_encoder.encode(label));
    }

    for (uint64_t target : path.target_columns) {
        if (!ref_targets.count(target))
            return false;
    }

    return true;
}


void
process_seq_path(const DeBruijnGraph &graph,
                 std::string_view query,
                 const std::vector<DeBruijnGraph::node_index> &query_nodes,
                 const std::function<void(AnnotatedDBG::row_index, size_t)> &callback) {
    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    if (canonical) {
        if (query_nodes.size()) {
            auto first = std::find_if(query_nodes.begin(), query_nodes.end(),
                                      [](auto i) -> bool { return i; });
            if (first == query_nodes.end())
                return;

            size_t start = first - query_nodes.begin();

            if (canonical->get_base_node(*first) == *first) {
                for (size_t i = start; i < query_nodes.size(); ++i) {
                    if (query_nodes[i] != DeBruijnGraph::npos) {
                        callback(AnnotatedDBG::graph_to_anno_index(
                            canonical->get_base_node(query_nodes[i])
                        ), i);
                    }
                }
            } else {
                for (size_t i = query_nodes.size(); i > start; --i) {
                    if (query_nodes[i - 1] != DeBruijnGraph::npos) {
                        callback(AnnotatedDBG::graph_to_anno_index(
                            canonical->get_base_node(query_nodes[i - 1])
                        ), i - 1);
                    }
                }
            }
        }
    } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        for (size_t i = 0; i < query_nodes.size(); ++i) {
            if (query_nodes[i] != DeBruijnGraph::npos)
                callback(AnnotatedDBG::graph_to_anno_index(query_nodes[i]), i);
        }
    } else {
        size_t i = 0;
        if (query.front() == '#') {
            std::string map_query
                = graph.get_node_sequence(query_nodes[0]).substr(0, graph.get_k() - 1);
            map_query += query.substr(graph.get_k() - 1);
            graph.map_to_nodes(map_query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    callback(AnnotatedDBG::graph_to_anno_index(node), i);

                ++i;
            });
        } else {
            graph.map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
                if (node != DeBruijnGraph::npos)
                    callback(AnnotatedDBG::graph_to_anno_index(node), i);

                ++i;
            });
        }
        assert(i == query_nodes.size());
    }
}

std::vector<std::vector<std::pair<uint64_t, size_t>>>
traverse_maximal_by_label(const AnnotatedDBG &anno_graph,
                          const std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> &seq_paths,
                          const Vector<uint64_t> &targets) {
    const DeBruijnGraph &graph = anno_graph.get_graph();
    std::vector<AnnotatedDBG::row_index> rows;

    for (const auto &[seq, nodes] : seq_paths) {
        assert(seq.size() == nodes.size() + graph.get_k() - 1);
        assert(std::find(nodes.begin(), nodes.end(), DeBruijnGraph::npos) == nodes.end());

        assert(std::all_of(nodes.begin(), nodes.end(), [&](auto node) {
            return graph.get_node_sequence(node).back() != boss::BOSS::kSentinel;
        }));

        if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph)) {
            for (size_t i = 0; i < nodes.size(); ++i) {
                rows.push_back(AnnotatedDBG::graph_to_anno_index(canonical->get_base_node(nodes[i])));
            }
        } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
            for (DeBruijnGraph::node_index node : nodes) {
                rows.push_back(AnnotatedDBG::graph_to_anno_index(node));
            }
        } else {
            assert(seq.front() == '#');
            graph.map_to_nodes(graph.get_node_sequence(nodes[0]) + seq.substr(graph.get_k()),
                               [&](DeBruijnGraph::node_index node) {
                assert(node);
                rows.push_back(AnnotatedDBG::graph_to_anno_index(node));
            });
        }
    }

    assert(std::all_of(rows.begin(), rows.end(), [&](auto row) {
        return graph.get_node_sequence(AnnotatedDBG::anno_to_graph_index(row)).back()
            != boss::BOSS::kSentinel;
    }));


    std::vector<std::vector<std::pair<uint64_t, size_t>>> result(seq_paths.size());

    auto masks = anno_graph.get_annotation().get_matrix().has_column(rows, targets);
    for (size_t m = 0; m < masks.size(); ++m) {
        uint64_t target = targets[m];
        size_t l = 0;
        size_t i = 0;
        result[l].emplace_back(target, 1);
        call_ones(masks[m], [&](size_t j) {
            assert(i + seq_paths.at(l).second.size() >= j);
            if (j == i + seq_paths[l].second.size()) {
                i += seq_paths[++l].second.size();
                result[l].emplace_back(target, 1);
            } else {
                result[l].back().second = j - i + 1;
            }
        });
    }

    return result;
}

Vector<uint64_t> LabeledSeedFilter::labels_to_keep(const DBGAlignment &seed) {
    Vector<uint64_t> labels;
    labels.reserve(seed.target_columns.size());

    auto targets = seed.target_columns;

    if (targets.empty())
        targets.push_back(std::numeric_limits<uint64_t>::max());

    for (uint64_t target : targets) {
        std::pair<size_t, size_t> idx_range {
            seed.get_clipping(), seed.get_clipping() + k_ - seed.get_offset()
        };
        size_t found_count = 0;
        for (node_index node : seed) {
            auto emplace = visited_nodes_[target].emplace(node, idx_range);
            auto &range = emplace.first.value();
            if (emplace.second) {
            } else if (range.first > idx_range.first || range.second < idx_range.second) {
                DEBUG_LOG("Node: {}; Prev_range: [{},{})", node, range.first, range.second);
                range.first = std::min(range.first, idx_range.first);
                range.second = std::max(range.second, idx_range.second);
                DEBUG_LOG("Node: {}; cur_range: [{},{})", node, range.first, range.second);
            } else {
                ++found_count;
            }

            if (idx_range.second - idx_range.first == k_)
                ++idx_range.first;

            ++idx_range.second;
        }

        if (found_count != seed.size())
            labels.push_back(target);
    }

    return labels;
}

void LabeledSeedFilter::update_seed_filter(const LabeledNodeRangeGenerator &generator) {
    generator([&](node_index node, uint64_t label, size_t begin, size_t end) {
        auto emplace = visited_nodes_[label].emplace(node, std::make_pair(begin, end));
        auto &range = emplace.first.value();
        if (!emplace.second) {
            range.first = std::min(range.first, begin);
            range.second = std::max(range.second, end);
        }
    });
}

ILabeledDBGAligner::ILabeledDBGAligner(const AnnotatedDBG &anno_graph,
                                       const DBGAlignerConfig &config)
      : anno_graph_(anno_graph),
        graph_(anno_graph_.get_graph()),
        config_(config) {}

auto ILabeledDBGAligner
::map_and_label_query_batch(const QueryGenerator &generate_query) const
        -> std::pair<BatchMapping, BatchLabels> {
    // exact k-mer matchings of each query sequence
    BatchMapping query_nodes;

    // map from Annotation row indices to (i,j), indicating position j in query i
    typedef std::vector<std::pair<size_t, size_t>> RowMapping;
    VectorMap<AnnotatedDBG::row_index, std::pair<RowMapping, RowMapping>> row_to_query_idx;

    // populate maps
    generate_query([&](std::string_view, std::string_view query, bool) {
        size_t i = query_nodes.size();

        // map query sequence to the graph
        query_nodes.emplace_back(map_sequence_to_nodes(graph_, query),
                                 std::vector<node_index>{});

        // update row_to_query_idx
        process_seq_path(graph_, query, query_nodes.back().first,
                         [&](AnnotatedDBG::row_index row, size_t j) {
            row_to_query_idx[row].first.emplace_back(i, j);
        });

        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            query_nodes.back().second = query_nodes.back().first;
            std::string query_rc(query);
            reverse_complement_seq_path(graph_, query_rc, query_nodes.back().second);

            if (graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                process_seq_path(graph_, query_rc, query_nodes.back().second,
                                 [&](AnnotatedDBG::row_index row, size_t j) {
                    row_to_query_idx[row].second.emplace_back(i, j);
                });
            }
        }
    });

    // extract rows from the row index map
    std::vector<AnnotatedDBG::row_index> rows;
    rows.reserve(row_to_query_idx.size());
    for (const auto &[row, mapping] : row_to_query_idx) {
        rows.push_back(row);
    }

    // get annotations for each row
    auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

    tsl::hopscotch_map<Vector<uint64_t>, std::vector<uint64_t>,
                       utils::VectorHash> unique_rows;
    for (size_t i = 0; i < annotation.size(); ++i) {
        unique_rows[annotation[i]].push_back(rows[i]);
    }

    std::vector<VectorMap<Vector<uint64_t>, std::pair<Signature, size_t>,
                          uint64_t, utils::VectorHash>> batch_labels(
        query_nodes.size()
    );
    for (const auto &[targets, rows] : unique_rows) {
        for (auto row : rows) {
            for (const auto &[i, idx] : row_to_query_idx[row].first) {
                if (!batch_labels[i].count(targets)) {
                    batch_labels[i].emplace(targets, std::make_pair(Signature {
                        { query_nodes[i].first.size(), false },
                        config_.forward_and_reverse_complement
                                && graph_.get_mode() != DeBruijnGraph::CANONICAL
                            ? sdsl::bit_vector(query_nodes[i].second.size(), false)
                            : sdsl::bit_vector()
                    }, (size_t)0));
                }

                if (!batch_labels[i][targets].first.first[idx]) {
                    batch_labels[i][targets].first.first[idx] = true;
                    ++batch_labels[i][targets].second;
                }
            }

            if (config_.forward_and_reverse_complement
                    && graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                for (const auto &[i, idx] : row_to_query_idx[row].second) {
                    if (!batch_labels[i].count(targets)) {
                        batch_labels[i].emplace(targets, std::make_pair(Signature {
                            { query_nodes[i].first.size(), false },
                            { query_nodes[i].second.size(), false }
                        }, (size_t)0));
                    }

                    if (!batch_labels[i][targets].first.second[idx]) {
                        batch_labels[i][targets].first.second[idx] = true;
                        ++batch_labels[i][targets].second;
                    }
                }
            }
        }
    }

    if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        for (size_t i = 0; i < batch_labels.size(); ++i) {
            assert(batch_labels[i].size()
                || (std::all_of(query_nodes[i].first.begin(),
                                query_nodes[i].first.end(),
                                [](auto j) { return !j; })
                    && std::all_of(query_nodes[i].second.begin(),
                                   query_nodes[i].second.end(),
                                   [](auto j) { return !j; })));
            for (auto it = batch_labels[i].begin(); it != batch_labels[i].end(); ++it) {
                it.value().first.second = it.value().first.first;
                reverse_bit_vector(it.value().first.second);
                it.value().second *= 2;
            }
        }
    }

    BatchLabels target_columns(query_nodes.size());
    for (size_t i = 0; i < target_columns.size(); ++i) {
        auto batch_labels_vec = const_cast<std::vector<std::pair<Vector<uint64_t>,
                                                                 std::pair<Signature, size_t>>>&&>(
            batch_labels[i].values_container()
        );
        if (batch_labels_vec.size() > config_.num_top_labels) {
            std::sort(batch_labels_vec.begin(), batch_labels_vec.end(),
                      [](const auto &a, const auto &b) {
                return std::make_pair(static_cast<bool>(a.first.size()), a.second.second)
                    < std::make_pair(static_cast<bool>(b.first.size()), b.second.second);
            });
            batch_labels_vec.resize(
                config_.num_top_labels + (batch_labels_vec[0].first.empty())
            );
            if (batch_labels_vec[0].first.empty()) {
                std::rotate(batch_labels_vec.begin(), batch_labels_vec.begin() + 1,
                            batch_labels_vec.end());
            }
        }

        for (auto&& [targets, signature_counts] : batch_labels_vec) {
            target_columns[i].emplace_back(std::move(targets),
                                           std::move(signature_counts.first));
        }
    }

    return { std::move(query_nodes), std::move(target_columns) };
}

template <typename NodeType>
LabeledColumnExtender<NodeType>
::LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                        const DBGAlignerConfig &config,
                        std::string_view query,
                        LabeledSeedFilter &seed_filter,
                        TargetColumnsSet &target_columns,
                        EdgeSetCache &cached_edge_sets)
      : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config, query),
        anno_graph_(anno_graph),
        seed_filter_(std::shared_ptr<LabeledSeedFilter>{}, &seed_filter),
        target_columns_(std::shared_ptr<TargetColumnsSet>{}, &target_columns),
        cached_edge_sets_(std::shared_ptr<EdgeSetCache>{}, &cached_edge_sets) {}

template <typename NodeType>
void LabeledColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    assert(!path.get_offset() || path.target_columns.empty());
    DefaultColumnExtender<NodeType>::initialize(path);
    align_node_to_target_.clear();
    backtrack_start_counter_.clear();

    assert(check_targets(anno_graph_, path));

    auto it = target_columns_->emplace(path.target_columns).first;
    size_t target_column_idx = it - target_columns_->begin();
    assert(target_column_idx < target_columns_->size());
    align_node_to_target_[this->start_node_] = target_column_idx;
}

template <typename NodeType>
auto LabeledColumnExtender<NodeType>
::get_next_align_node(NodeType node, char c, size_t dist_from_origin) const -> AlignNode {
    size_t depth = 0;
    auto find = this->table_.find(node);
    if (find != this->table_.end())
        depth = find->second.first.size();

    return { node, c, depth, dist_from_origin };
}

template <typename NodeType>
auto LabeledColumnExtender<NodeType>::get_outgoing(const AlignNode &node) const -> Edges {
    assert(align_node_to_target_.count(node));

    if (std::get<0>(node) == this->graph_.max_index() + 1)
        return DefaultColumnExtender<NodeType>::get_outgoing(node);

    size_t depth = std::get<3>(node);
    assert(depth);

    uint64_t target_column_idx = align_node_to_target_[node];
    assert(target_column_idx < target_columns_->size());

    if ((target_columns_->begin() + target_column_idx)->empty()) {
        assert(this->seed_->get_offset() >= depth);
        assert(this->seed_->get_offset());
        size_t next_offset = this->seed_->get_offset() - depth;

        assert(this->graph_.get_node_sequence(std::get<0>(node)).substr(next_offset + 1).find(boss::BOSS::kSentinel)
            == std::string::npos);

        Edges edges = DefaultColumnExtender<NodeType>::get_outgoing(node);
        if (next_offset || edges.empty())
            return edges;

        std::string seq(this->graph_.get_k(), '#');
        std::vector<AnnotatedDBG::row_index> rows;
        rows.reserve(edges.size());

        for (const auto &[next_node, c] : edges) {
            seq.back() = c;
            process_seq_path(this->graph_, seq, std::vector<node_index>{ next_node },
                             [&](AnnotatedDBG::row_index row, size_t) {
                rows.emplace_back(row);
            });
        }

        assert(rows.size() == edges.size());

        auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

        Edges out_edges;
        out_edges.reserve(edges.size());
        for (size_t i = 0; i < rows.size(); ++i) {
            DBGAlignment dummy_seed(
                std::string_view(this->seed_->get_query().data(), this->graph_.get_k()),
                std::vector<node_index>{ edges[i].first },
                std::string(seq), 0, Cigar(), this->seed_->get_clipping(),
                this->seed_->get_orientation()
            );
            dummy_seed.target_columns = std::move(annotation[i]);
            annotation[i] = seed_filter_->labels_to_keep(dummy_seed);
            if (annotation[i].empty() || (annotation[i].size() == 1
                    && annotation[i][0] == std::numeric_limits<uint64_t>::max())) {
                DEBUG_LOG("Skipping seed: {}", dummy_seed);
                continue;
            }

            out_edges.emplace_back(std::move(edges[i]));

            AlignNode next = get_next_align_node(edges[i].first, edges[i].second, depth + 1);
            auto it = target_columns_->emplace(annotation[i]).first;
            target_column_idx = it - target_columns_->begin();
            assert(target_column_idx < target_columns_->size());
            assert(!align_node_to_target_.count(next));
            align_node_to_target_[next] = target_column_idx;
        }

        return out_edges;
#ifndef NDEBUG
    } else {
        assert(this->seed_->get_offset() < depth);
        assert(this->graph_.get_node_sequence(std::get<0>(node)).find(boss::BOSS::kSentinel)
            == std::string::npos);
#endif
    }

    if (cached_edge_sets_->count(std::get<0>(node))) {
        const auto &edge_sets = (*cached_edge_sets_)[std::get<0>(node)];
        if (edge_sets.count(target_column_idx)) {
            const LabeledEdges &labeled_edges = edge_sets.find(target_column_idx)->second;
            Edges edges;
            edges.reserve(labeled_edges.size());
            for (auto [node, c, target_column_idx] : labeled_edges) {
                edges.emplace_back(node, c);
                AlignNode next = get_next_align_node(node, c, depth + 1);
                if (!align_node_to_target_.count(next))
                    align_node_to_target_[next] = target_column_idx;
            }
            return edges;
        }
    }

    // prefetch the next unitig
    std::vector<std::pair<std::string, std::vector<NodeType>>> seq_paths;

    const auto &column = this->table_.find(std::get<0>(node))->second;
    auto [min_i, max_i] = this->get_band(node, column, this->xdrop_cutoff_);

    const auto &S_vec = std::get<0>(column.first[std::get<2>(node)]);
    size_t offset = std::get<9>(column.first[std::get<2>(node)]);
    std::vector<score_t> S(S_vec.begin() + min_i - offset, S_vec.begin() + max_i - offset);

    std::string_view q_min = this->query_.substr(this->start_ + min_i);

    std::string seq(this->graph_.get_k(), '#');
    VectorSet<node_index> visited;
    auto &path = const_cast<std::vector<node_index>&>(visited.values_container());
    node_index cur_node = std::get<0>(node);
    visited.emplace(cur_node);
    std::vector<std::pair<node_index, char>> outgoing;

    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_);

    while (cur_node != DeBruijnGraph::npos) {
        outgoing.clear();
        this->graph_.call_outgoing_kmers(cur_node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel)
                outgoing.emplace_back(next, c);
        });

        cur_node = DeBruijnGraph::npos;

        if (outgoing.empty()) {
            if (path.size())
                seq_paths.emplace_back(std::move(seq), std::move(path));
        } else if (outgoing.size() > 1) {
            if (path.size() == 1) {
                seq += '#';
                path.push_back(0);
                for (const auto &[next, c] : outgoing) {
                    if (visited.count(next))
                        continue;

                    seq.back() = c;
                    path.back() = next;
                    seq_paths.emplace_back(seq, path);
                }
            } else {
                seq_paths.emplace_back(std::move(seq), std::move(path));
            }
        } else if (visited.count(outgoing[0].first)) {
            if (path.size() > 1)
                seq_paths.emplace_back(std::move(seq), std::move(path));
        } else if ((path.size() && canonical
                    && (canonical->get_base_node(path.back()) == path.back())
                        != (canonical->get_base_node(outgoing[0].first)
                            == outgoing[0].first))) {
            seq.push_back(outgoing[0].second);
            visited.emplace(outgoing[0].first);
            seq_paths.emplace_back(std::move(seq), std::move(path));
        } else {
            seq.push_back(outgoing[0].second);
            visited.emplace(outgoing[0].first);

            const auto &score_row = this->config_.get_row(seq.back());
            for (size_t i = 0; i < S.size() && i + path.size() >= q_min.size(); ++i) {
                S[i] += score_row[q_min[i + path.size()]];
                if (S[i] >= this->xdrop_cutoff_)
                    cur_node = outgoing[0].first;
            }

            if (cur_node == DeBruijnGraph::npos)
                seq_paths.emplace_back(std::move(seq), std::move(path));
        }
    }

    Edges edges;
    LabeledEdges labeled_edges;
    auto all_extensions = traverse_maximal_by_label(
        anno_graph_, seq_paths, *(target_columns_->begin() + target_column_idx)
    );
    for (size_t j = 0; j < seq_paths.size(); ++j) {
        const auto &[seq, path] = seq_paths[j];
        auto &extensions = all_extensions[j];
        if (path.size() > 1) {
            assert(path[0] == std::get<0>(node));
            std::sort(extensions.begin(), extensions.end(), utils::LessSecond());

            if (extensions.back().second <= 1)
                continue;

            edges.emplace_back(path[1], seq[this->graph_.get_k()]);

            auto it = std::find_if(extensions.begin(), extensions.end(),
                                   [&](const auto &a) { return a.second > 1; });
            Targets targets(extensions.end() - it);
            for (auto jt = it; jt != extensions.end(); ++jt) {
                targets[jt - it] = jt->first;
            }
            std::sort(targets.begin(), targets.end());
            auto jt = target_columns_->emplace(targets).first;
            size_t cur_target_column_idx = jt - target_columns_->begin();
            assert(cur_target_column_idx < target_columns_->size());

            if (it != extensions.begin()) {
                AlignNode next = get_next_align_node(
                    path[1], seq[this->graph_.get_k()], depth + 1);
                assert(!align_node_to_target_.count(next));
                align_node_to_target_[next] = cur_target_column_idx;
            }

            labeled_edges.emplace_back(edges.back().first, edges.back().second, cur_target_column_idx);

            for (size_t i = 2; i < extensions.front().second; ++i) {
                (*cached_edge_sets_)[path[i - 1]][cur_target_column_idx].emplace_back(
                    path[i], seq[i + this->graph_.get_k() - 1], cur_target_column_idx
                );
                if (it->second < i) {
                    auto jt = std::find_if(it + 1, extensions.end(),
                                           [&](const auto &a) { return a.second >= i; });
                    assert(it + targets.size() >= jt);
                    targets.resize(extensions.end() - jt);
                    for (auto kt = jt; kt != extensions.end(); ++kt) {
                        targets[kt - jt] = kt->first;
                    }
                    std::sort(targets.begin(), targets.end());
                    it = jt;
                    auto kt = target_columns_->emplace(targets).first;
                    cur_target_column_idx = kt - target_columns_->begin();
                    assert(cur_target_column_idx < target_columns_->size());
                }

                AlignNode next = get_next_align_node(
                    path[i], seq[i + this->graph_.get_k() - 1], depth + i);
                assert(!align_node_to_target_.count(next));
                align_node_to_target_[next] = cur_target_column_idx;
            }
        }
    }

    auto &edge_sets = (*cached_edge_sets_)[std::get<0>(node)];
    assert(!edge_sets.count(target_column_idx));
    edge_sets[target_column_idx] = labeled_edges;

    return edges;
}

template class LabeledColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
