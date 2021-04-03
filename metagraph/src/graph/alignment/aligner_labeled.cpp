#include "aligner_labeled.hpp"

#include <tsl/hopscotch_set.h>

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
                return a.second.second > b.second.second;
            });
            auto it = batch_labels_vec.begin() + config_.num_top_labels - 1;
            auto end = std::find_if(it + 1, batch_labels_vec.end(),
                                    [&](const auto &a) {
                                        return a.second.second < it->second.second;
                                    });
            batch_labels_vec.erase(end, batch_labels_vec.end());
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
    align_node_to_target_[this->start_node_].first = get_target_id(path.target_columns);
}

template <typename NodeType>
auto LabeledColumnExtender<NodeType>::get_outgoing(const AlignNode &align_node) const -> Edges {
    assert(align_node_to_target_.count(align_node));
    const auto &[node, last_char, count, depth] = align_node;

    Edges base_edges = DefaultColumnExtender<NodeType>::get_outgoing(align_node);

    if (node == this->graph_.max_index() + 1 || base_edges.empty())
        return base_edges;

    assert(depth);

    auto [target_column_idx, lookahead] = align_node_to_target_[align_node];

    if (!target_column_idx) {
        assert(this->seed_->get_offset() >= depth);
        assert(this->seed_->get_offset());
        size_t next_offset = this->seed_->get_offset() - depth;

        assert(this->graph_.get_node_sequence(node).substr(next_offset + 1).find(boss::BOSS::kSentinel)
            == std::string::npos);

        if (next_offset)
            return base_edges;

        std::string seq(this->graph_.get_k(), '#');
        std::vector<AnnotatedDBG::row_index> rows;
        rows.reserve(base_edges.size());

        for (const auto &[next_node, c, next_count, next_depth] : base_edges) {
            seq.back() = c;
            process_seq_path(this->graph_, seq, std::vector<node_index>{ next_node },
                             [&](AnnotatedDBG::row_index row, size_t) {
                rows.emplace_back(row);
            });
        }

        assert(rows.size() == base_edges.size());

        auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

        Edges out_edges;
        out_edges.reserve(base_edges.size());
        for (size_t i = 0; i < rows.size(); ++i) {
            DBGAlignment dummy_seed(
                std::string_view(this->seed_->get_query().data(), this->graph_.get_k()),
                std::vector<node_index>{ std::get<0>(base_edges[i]) },
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

            target_column_idx = get_target_id(annotation[i]);

            update_target_cache(std::get<0>(base_edges[i]), target_column_idx);

            assert(!align_node_to_target_.count(base_edges[i]));
            align_node_to_target_.emplace(base_edges[i],
                                          std::make_pair(target_column_idx, 0));

            out_edges.emplace_back(std::move(base_edges[i]));
        }

        return out_edges;
    }

    assert(this->seed_->get_offset() < depth);
    assert(this->graph_.get_node_sequence(node).find(boss::BOSS::kSentinel)
        == std::string::npos);

    Edges out_edges;

    const Targets &start_targets = get_targets(target_column_idx);

    tsl::hopscotch_set<node_index> visited;
    std::vector<size_t> dist;
    std::vector<AnnotatedDBG::row_index> rows;
    std::vector<std::pair<AlignNode, size_t>> next_nodes;
    std::vector<node_index> nodes;

    for (size_t i = 0; i < base_edges.size(); ++i) {
        const auto &[next_node, c, next_count, next_depth] = base_edges[i];
        visited.emplace(next_node);

        if (next_node == node) {
            assert(!align_node_to_target_.count(base_edges[i]));
            align_node_to_target_.emplace(base_edges[i],
                                          std::make_pair(target_column_idx, lookahead));
            out_edges.emplace_back(base_edges[i]);

        } else {
            if (cached_edge_sets_->count(next_node)) {
                const Targets &next_targets = get_targets((*cached_edge_sets_)[next_node]);
                size_t next_idx = 0;
                if (std::includes(next_targets.begin(), next_targets.end(),
                                  start_targets.begin(), start_targets.end())) {
                    next_idx = target_column_idx;
                } else if (lookahead) {
                    next_idx = get_target_intersection((*cached_edge_sets_)[next_node],
                                                       target_column_idx);
                }

                if (next_idx) {
                    assert(!align_node_to_target_.count(base_edges[i]));
                    align_node_to_target_.emplace(base_edges[i],
                                                  std::make_pair(next_idx, lookahead ? lookahead - 1 : 0));
                    out_edges.emplace_back(base_edges[i]);
                    continue;
                }
            }

            next_nodes.emplace_back(base_edges[i], 0);
        }
    }

    if (next_nodes.size()) {
        for (size_t i = 0; i < next_nodes.size(); ++i) {
            const auto &[next_node, c, next_count, next_depth] = next_nodes[i].first;

            std::string seq(this->graph_.get_k() - 1, '#');
            seq += c;

            std::vector<node_index> path { next_node };
            node_index cur = DeBruijnGraph::npos;

            const auto &column = this->table_.find(node)->second;
            auto [min_i, max_i] = this->get_band(align_node, column, this->xdrop_cutoff_);
            if (this->start_ + min_i + 1 < this->query_.size()) {
                std::string_view q_min = this->query_.substr(this->start_ + min_i + 1);
                size_t offset = std::get<9>(column.first[count]);
                const score_t *S_vec = std::get<0>(column.first[count]).data() + min_i - offset;
                std::vector<score_t> S(S_vec, S_vec + (max_i - min_i));

                const auto &score_row = this->config_.get_row(c);
                for (size_t j = 0; j < S.size() && j + path.size() < q_min.size(); ++j) {
                    S[j] += score_row[q_min[j + path.size()]];
                    if (S[j] >= this->xdrop_cutoff_)
                        cur = next_node;
                }

                while (cur != DeBruijnGraph::npos) {
                    std::vector<std::pair<node_index, char>> outgoing;
                    this->graph_.call_outgoing_kmers(cur, [&](node_index next, char next_c) {
                        outgoing.emplace_back(next, next_c);
                    });
                    cur = DeBruijnGraph::npos;

                    if (outgoing.size() == 1 && visited.emplace(outgoing[0].first).second) {
                        path.push_back(outgoing[0].first);
                        seq += outgoing[0].second;

                        const auto &score_row = this->config_.get_row(seq.back());
                        for (size_t j = 0; j < S.size() && j + path.size() < q_min.size(); ++j) {
                            S[j] += score_row[q_min[j + path.size()]];
                            if (S[j] >= this->xdrop_cutoff_)
                                cur = outgoing[0].first;
                        }
                    }
                }
            }

            next_nodes[i].second = path.size() - 1;

            process_seq_path(this->graph_, seq, path, [&](auto row, size_t d) {
                dist.push_back(d);
                rows.push_back(row);
                nodes.push_back(path[d]);
            });
        }

        assert(static_cast<size_t>(std::count(dist.begin(), dist.end(), 0)) == next_nodes.size());

        auto masks = anno_graph_.get_annotation().get_matrix().has_column(rows, start_targets);

        std::vector<Targets> out_targets(rows.size());
        for (size_t i = 0; i < masks.size(); ++i) {
            call_ones(masks[i], [&](size_t j) { out_targets[j].push_back(start_targets[i]); });
        }

        auto it = next_nodes.begin();
        for (size_t i = 0; i < out_targets.size(); ++i) {
            if (size_t target_idx = get_target_id(out_targets[i])) {
                update_target_cache(nodes[i], target_idx);

                if (!dist[i]) {
                    assert(it != next_nodes.end());
                    assert(std::get<0>(it->first) == nodes[i]);
                    assert(!align_node_to_target_.count(it->first));

                    out_edges.push_back(it->first);
                    align_node_to_target_.emplace(
                        it->first, std::make_pair(target_idx, it->second)
                    );
                }
            }

            if (!dist[i])
                ++it;
        }

        assert(it == next_nodes.end());
    }

    return out_edges;
}

template class LabeledColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
