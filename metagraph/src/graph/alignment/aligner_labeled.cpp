#include "aligner_labeled.hpp"

#include "graph/annotated_dbg.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/hashers/hash.hpp"

namespace mtg {
namespace graph {
namespace align {


inline void reverse_bit_vector(sdsl::bit_vector &v) {
    size_t begin = 0;
    for ( ; begin + begin + 128 <= v.size(); begin += 64) {
        uint64_t a = sdsl::bits::rev(v.get_int(begin));
        uint64_t b = sdsl::bits::rev(v.get_int(v.size() - begin - 64));
        v.set_int(begin, b);
        v.set_int(v.size() - begin - 64, a);
    }

    size_t size = (v.size() % 128) / 2;
    uint64_t a = sdsl::bits::rev(v.get_int(begin, size)) >> (64 - size);
    uint64_t b = sdsl::bits::rev(v.get_int(v.size() - begin - size, size)) >> (64 - size);
    v.set_int(begin, b, size);
    v.set_int(v.size() - begin - size, a, size);
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

std::vector<std::pair<uint64_t, size_t>>
traverse_maximal_by_label(const AnnotatedDBG &anno_graph,
                          const std::string &seq,
                          const std::vector<DeBruijnGraph::node_index> &nodes,
                          const Vector<uint64_t> &targets) {
    const DeBruijnGraph &graph = anno_graph.get_graph();
    assert(seq.size() == nodes.size() + graph.get_k() - 1);
    assert(std::find(nodes.begin(), nodes.end(), DeBruijnGraph::npos) == nodes.end());

    assert(std::all_of(nodes.begin(), nodes.end(), [&](auto node) {
        return graph.get_node_sequence(node).back() != boss::BOSS::kSentinel;
    }));

    std::vector<AnnotatedDBG::row_index> rows(nodes.size());
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph)) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            rows[i] = AnnotatedDBG::graph_to_anno_index(canonical->get_base_node(nodes[i]));
        }
    } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        std::transform(nodes.begin(), nodes.end(), rows.begin(),
                       AnnotatedDBG::graph_to_anno_index);
    } else {
        assert(seq.front() == '#');
        size_t i = 0;
        graph.map_to_nodes(graph.get_node_sequence(nodes[0]) + seq.substr(graph.get_k()),
                           [&](DeBruijnGraph::node_index node) {
            assert(node);
            rows[i++] = AnnotatedDBG::graph_to_anno_index(node);
        });
        assert(i == nodes.size());
    }

    assert(std::all_of(rows.begin(), rows.end(), [&](auto row) {
        return graph.get_node_sequence(AnnotatedDBG::anno_to_graph_index(row)).back()
            != boss::BOSS::kSentinel;
    }));

    auto extensions = anno_graph.get_annotation().get_matrix().extend_maximal(rows, targets);
    std::vector<std::pair<uint64_t, size_t>> result;
    result.reserve(extensions.size());
    for (size_t i = 0; i < extensions.size(); ++i) {
        result.emplace_back(targets[i], extensions[i]);
    }
    return result;
}

Vector<uint64_t> LabeledSeedFilter::labels_to_keep(const DBGAlignment &seed) {
    Vector<uint64_t> labels;
    labels.reserve(seed.target_columns.size());

    std::pair<size_t, size_t> idx_range {
        seed.get_clipping(), seed.get_clipping() + k_ - seed.get_offset()
    };
    for (uint64_t target : seed.target_columns) {
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

    tsl::hopscotch_map<Vector<uint64_t>, std::vector<uint64_t>, utils::VectorHash> unique_rows;
    for (size_t i = 0; i < annotation.size(); ++i) {
        unique_rows[annotation[i]].push_back(rows[i]);
    }

    std::vector<VectorMap<Vector<uint64_t>, Signature,
                          uint64_t, utils::VectorHash>> batch_labels(
        query_nodes.size()
    );
    for (const auto &[targets, rows] : unique_rows) {
        for (auto row : rows) {
            for (const auto &[i, idx] : row_to_query_idx[row].first) {
                if (!batch_labels[i].count(targets)) {
                    batch_labels[i].emplace(targets, Signature {
                        { query_nodes[i].first.size(), false },
                        config_.forward_and_reverse_complement
                                && graph_.get_mode() != DeBruijnGraph::CANONICAL
                            ? sdsl::bit_vector(query_nodes[i].second.size(), false)
                            : sdsl::bit_vector()
                    });
                }

                batch_labels[i][targets].first[idx] = true;
            }

            if (config_.forward_and_reverse_complement
                    && graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                for (const auto &[i, idx] : row_to_query_idx[row].second) {
                    if (!batch_labels[i].count(targets)) {
                        batch_labels[i].emplace(targets, Signature {
                            { query_nodes[i].first.size(), false },
                            { query_nodes[i].second.size(), false }
                        });
                    }

                    batch_labels[i][targets].second[idx] = true;
                }
            }
        }
    }

    if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        for (auto &vmap : batch_labels) {
            for (auto it = vmap.begin(); it != vmap.end(); ++it) {
                it.value().second = it.value().first;
                reverse_bit_vector(it.value().second);
            }
        }
    }

    BatchLabels target_columns(query_nodes.size());
    for (size_t i = 0; i < target_columns.size(); ++i) {
        target_columns[i] = std::move(const_cast<QueryLabels&&>(
            batch_labels[i].values_container()
        ));
    }

    return { std::move(query_nodes), std::move(target_columns) };
}

template <typename NodeType>
LabeledColumnExtender<NodeType>
::LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                        const DBGAlignerConfig &config,
                        std::string_view query)
      : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config, query),
        anno_graph_(anno_graph) {}

template <typename NodeType>
void LabeledColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    assert(!path.get_offset() || path.target_columns.empty());
    DefaultColumnExtender<NodeType>::initialize(path);
    align_node_to_target_.clear();

    auto it = target_columns_.emplace(path.target_columns).first;
    size_t target_column_idx = it - target_columns_.begin();
    assert(target_column_idx < target_columns_.size());
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
    assert(target_column_idx < target_columns_.size());

    if ((target_columns_.begin() + target_column_idx)->empty()) {
        std::cerr << "seed\ta\t" << *this->seed_ << " " << this->seed_->target_columns.size() << " " << depth << "\n";
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

        for (size_t i = 0; i < rows.size(); ++i) {
            AlignNode next = get_next_align_node(edges[i].first, edges[i].second, depth + 1);
            auto it = target_columns_.emplace(annotation[i]).first;
            target_column_idx = it - target_columns_.begin();
            assert(target_column_idx < target_columns_.size());
            assert(!align_node_to_target_.count(next));
            align_node_to_target_[next] = target_column_idx;
        }

        return edges;
#ifndef NDEBUG
    } else {
        std::cerr << "seed\tb\t" << *this->seed_ << " " << this->seed_->target_columns.size() << " " << depth << "\n";
        assert(this->seed_->get_offset() < depth);
        assert(this->graph_.get_node_sequence(std::get<0>(node)).find(boss::BOSS::kSentinel)
            == std::string::npos);
#endif
    }

    if (cached_edge_sets_.count(std::get<0>(node))) {
        const auto &edge_sets = cached_edge_sets_[std::get<0>(node)];
        if (edge_sets.count(target_column_idx))
            return edge_sets.find(target_column_idx)->second;
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
                    if (visited.count(next) || (cached_edge_sets_.count(next)
                            && cached_edge_sets_[next].count(target_column_idx))) {
                        continue;
                    }

                    seq.back() = c;
                    path.back() = next;
                    seq_paths.emplace_back(seq, path);
                }
            } else {
                seq_paths.emplace_back(std::move(seq), std::move(path));
            }
        } else if (visited.count(outgoing[0].first)
                || (cached_edge_sets_.count(outgoing[0].first)
                    && cached_edge_sets_[outgoing[0].first].count(target_column_idx))) {
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
    for (const auto &[seq, path] : seq_paths) {
        if (path.size() > 1) {
            assert(path[0] == std::get<0>(node));
            auto extensions = traverse_maximal_by_label(
                anno_graph_, seq, path,
                *(target_columns_.begin() + target_column_idx)
            );
            std::sort(extensions.begin(), extensions.end(), utils::LessSecond());
            size_t max_ext = extensions.back().second;

            if (max_ext <= 1)
                continue;

            edges.emplace_back(path[1], seq[this->graph_.get_k()]);

            auto it = std::find_if(extensions.begin(), extensions.end(),
                                   [&](const auto &a) { return a.second > 1; });
            Targets targets(extensions.end() - it);
            for (auto jt = it; jt != extensions.end(); ++jt) {
                targets[jt - it] = jt->first;
            }
            Targets sorted_targets = targets;
            std::sort(sorted_targets.begin(), sorted_targets.end());
            size_t cur_target_column_idx;
            if (it != extensions.begin()) {
                auto jt = target_columns_.emplace(sorted_targets).first;
                cur_target_column_idx = jt - target_columns_.begin();
                assert(cur_target_column_idx < target_columns_.size());

                AlignNode next = get_next_align_node(
                    path[1], seq[this->graph_.get_k()], depth + 1);
                assert(!align_node_to_target_.count(next));
                align_node_to_target_[next] = cur_target_column_idx;
            }

            for (size_t i = 2; i < max_ext; ++i) {
                cached_edge_sets_[path[i - 1]][cur_target_column_idx].emplace_back(
                    path[i], seq[i + this->graph_.get_k() - 1]
                );
                if (it->second < i) {
                    auto jt = std::find_if(it + 1, extensions.end(),
                                           [&](const auto &a) { return a.second >= i; });
                    assert(it + targets.size() >= jt);
                    targets.erase(targets.begin(), targets.begin() + (jt - it));
                    sorted_targets = targets;
                    std::sort(sorted_targets.begin(), sorted_targets.end());
                    it = jt;
                    auto kt = target_columns_.emplace(targets).first;
                    cur_target_column_idx = kt - target_columns_.begin();
                    assert(cur_target_column_idx < target_columns_.size());
                }

                AlignNode next = get_next_align_node(
                    path[i], seq[i + this->graph_.get_k() - 1], depth + i);
                assert(!align_node_to_target_.count(next));
                align_node_to_target_[next] = cur_target_column_idx;
            }
        }
    }

    auto &edge_sets = cached_edge_sets_[std::get<0>(node)];
    assert(!edge_sets.count(target_column_idx));
    edge_sets[target_column_idx] = edges;

    return edges;
}

template class LabeledColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
