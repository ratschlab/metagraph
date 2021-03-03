#include "aligner_labeled.hpp"

#include "graph/annotated_dbg.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/vector_map.hpp"
#include "common/utils/template_utils.hpp"

namespace mtg {
namespace graph {
namespace align {

void
process_seq_path(const DeBruijnGraph &graph,
                 std::string_view query,
                 const std::vector<DeBruijnGraph::node_index> &query_nodes,
                 const std::function<void(AnnotatedDBG::row_index, size_t)> &callback) {
    const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    if (canonical) {
        for (size_t i = 0; i < query_nodes.size(); ++i) {
            if (query_nodes[i] != DeBruijnGraph::npos) {
                callback(AnnotatedDBG::graph_to_anno_index(
                    canonical->get_base_node(query_nodes[i])
                ), i);
            }
        }
    } else if (graph.get_mode() != DeBruijnGraph::CANONICAL) {
        for (size_t i = 0; i < query_nodes.size(); ++i) {
            if (query_nodes[i] != DeBruijnGraph::npos)
                callback(AnnotatedDBG::graph_to_anno_index(query_nodes[i]), i);
        }
    } else {
        size_t i = 0;
        graph.map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
            if (node != DeBruijnGraph::npos)
                callback(AnnotatedDBG::graph_to_anno_index(node), i);

            ++i;
        });
        assert(i == query_nodes.size());
    }
}

ILabeledDBGAligner::ILabeledDBGAligner(const AnnotatedDBG &anno_graph,
                                       const DBGAlignerConfig &config,
                                       size_t num_top_labels)
      : anno_graph_(anno_graph),
        graph_(anno_graph_.get_graph()),
        config_(config), num_top_labels_(num_top_labels) {}

auto ILabeledDBGAligner
::map_and_label_query_batch(const QueryGenerator &generate_query) const
        -> std::pair<BatchMapping, BatchLabels> {
    // exact k-mer matchings of each query sequence
    BatchMapping query_nodes;

    // map from Annotation row indices to (i,j), indicating position j in query i
    VectorMap<AnnotatedDBG::row_index,
              std::vector<std::pair<size_t, size_t>>> row_to_query_idx;

    // populate maps
    generate_query([&](std::string_view, std::string_view query, bool) {
        size_t i = query_nodes.size();

        // map query sequence to the graph
        query_nodes.emplace_back(map_sequence_to_nodes(graph_, query));

        // update row_to_query_idx
        process_seq_path(graph_, query, query_nodes.back(),
                         [&](AnnotatedDBG::row_index row, size_t j) {
            row_to_query_idx[row].emplace_back(i, j);
        });
    });

    // extract rows from the row index map
    std::vector<AnnotatedDBG::row_index> rows;
    rows.reserve(row_to_query_idx.size());
    for (const auto &[row, mapping] : row_to_query_idx) {
        rows.push_back(row);
    }

    // get annotations for each row
    auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

    // we now have two maps
    // 1) row indices -> column IDs
    // 2) row indices -> (query i, query_i j) i.e., row -> the k-mer at batch[i][j]
    // and wish to construct the following map
    // query -> ( (col_1, col_1_signature), (col_2, col_2_signature), ... )

    // count labels for each query
    std::vector<VectorMap<uint64_t, uint64_t>> column_counter(query_nodes.size());
    for (size_t i = 0; i < annotation.size(); ++i) {
        for (const auto &[query_id, idx] : row_to_query_idx[rows[i]]) {
            for (uint64_t column : annotation[i]) {
                ++column_counter[query_id][column];
            }
        }
    }

    // compute target columns and initialize signatures for each query
    BatchLabels target_columns(query_nodes.size());
    for (size_t i = 0; i < column_counter.size(); ++i) {
        auto &counter = const_cast<std::vector<std::pair<uint64_t, uint64_t>>&>(
            column_counter[i].values_container()
        );

        // pick the top columns for each query sequence
        size_t num_targets = std::min(counter.size(), num_top_labels_);

        std::sort(counter.begin(), counter.end(), utils::GreaterSecond());
        if (num_targets < counter.size())
            counter.resize(num_targets);

        for (const auto &[column, count] : counter) {
            if (!target_columns[i].count(column)) {
                target_columns[i].emplace(column,
                                          sdsl::bit_vector(query_nodes[i].size(), false));
            }
        }
    }

    // fill signatures for each query
    for (size_t i = 0; i < annotation.size(); ++i) {
        for (const auto &[query_id, idx] : row_to_query_idx[rows[i]]) {
            for (uint64_t column : annotation[i]) {
                if (target_columns[query_id].count(column))
                    target_columns[query_id][column][idx] = true;
            }
        }
    }

    // if any of the queries has no associated labels, mark them as such
    for (size_t i = 0; i < target_columns.size(); ++i) {
        if (target_columns[i].empty()) {
            assert(!query_nodes[i][0]);
            assert(std::equal(query_nodes[i].begin() + 1, query_nodes[i].end(),
                              query_nodes[i].begin()));
            target_columns[i][kNTarget] = sdsl::bit_vector(query_nodes[i].size(), true);
        }
    }

    return std::make_pair(std::move(query_nodes), std::move(target_columns));
}

template <typename NodeType>
LabeledColumnExtender<NodeType>
::LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                        const DBGAlignerConfig &config,
                        std::string_view query)
      : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config, query),
        anno_graph_(anno_graph),
        initial_target_columns_size_(0) {}

template <typename NodeType>
LabeledColumnExtender<NodeType>
::LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                        const DBGAlignerConfig &config,
                        std::string_view query,
                        uint64_t initial_target_column)
      : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config, query),
        anno_graph_(anno_graph),
        target_columns_({ initial_target_column }),
        initial_target_columns_size_(1) {}

template <typename NodeType>
void LabeledColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    target_columns_.resize(initial_target_columns_size_);
    assert(target_columns_.size() <= 1);
    DefaultColumnExtender<NodeType>::initialize(path);
    cached_edge_sets_.clear();
    align_node_to_target_.clear();

    AlignNode start_node{ this->graph_.max_index() + 1, '\0', 0, 0 };

    if (!path.get_offset()) {
        align_node_to_target_[start_node] = 0;
        if (target_columns_.empty())
            target_columns_.push_back(ILabeledDBGAligner::kNTarget);

    } else if (target_columns_.empty()) {
        align_node_to_target_[start_node] = 0;
        target_columns_.push_back(ILabeledDBGAligner::kNTarget);

    } else {
        align_node_to_target_[start_node] = 1;
        target_columns_.push_back(ILabeledDBGAligner::kNTarget);
    }
}

template <typename NodeType>
auto LabeledColumnExtender<NodeType>::get_outgoing(const AlignNode &node) const -> Edges {
    assert(align_node_to_target_.count(node));

    if (std::get<0>(node) == this->graph_.max_index() + 1)
        return DefaultColumnExtender<NodeType>::get_outgoing(node);

    uint64_t target_column_idx = align_node_to_target_[node];
    uint64_t target_column = target_columns_.at(target_column_idx);

    if (target_column == ILabeledDBGAligner::kNTarget) {
        assert(this->seed_->get_offset());
        assert(this->seed_->get_offset() + 1 >= std::get<3>(node));
        size_t next_offset = this->seed_->get_offset() + 1 - std::get<3>(node);

        Edges edges = DefaultColumnExtender<NodeType>::get_outgoing(node);
        if (next_offset)
            return edges;

        std::vector<AnnotatedDBG::row_index> rows;
        rows.reserve(edges.size());

        const CanonicalDBG *canonical = dynamic_cast<const CanonicalDBG*>(&this->graph_);
        for (const auto &[next_node, c] : edges) {
            if (canonical) {
                rows.emplace_back(AnnotatedDBG::graph_to_anno_index(
                    canonical->get_base_node(next_node)
                ));
            } else if (this->graph_.get_mode() != DeBruijnGraph::CANONICAL) {
                rows.emplace_back(AnnotatedDBG::graph_to_anno_index(next_node));
            } else {
                this->graph_.map_to_nodes(this->graph_.get_node_sequence(next_node),
                                          [&](NodeType next_canonical) {
                    rows.emplace_back(AnnotatedDBG::graph_to_anno_index(next_canonical));
                });
            }
        }

        assert(rows.size() == edges.size());

        auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

        Edges out_edges;
        for (size_t i = 0; i < rows.size(); ++i) {
            AlignNode next { edges[i].first, edges[i].second, 0, std::get<3>(node) + 1 };
            auto find = this->table_.find(edges[i].first);
            if (find != this->table_.end())
                std::get<2>(next) = find->second.first.size();

            for (uint64_t target : annotation[i]) {
                target_column_idx = std::find(target_columns_.begin(),
                                              target_columns_.end(),
                                              target) - target_columns_.begin();

                if (target_column_idx == target_columns_.size())
                    target_columns_.emplace_back(target);

                assert(!align_node_to_target_.count(next));
                align_node_to_target_[next] = target_column_idx;

                out_edges.push_back(edges[i]);
                ++std::get<2>(next);
            }
        }

        return out_edges;
    }

    if (cached_edge_sets_.count(std::get<0>(node))) {
        const auto &edge_sets = cached_edge_sets_[std::get<0>(node)];
        if (target_column_idx < edge_sets.size())
            return edge_sets[target_column_idx];
    }

    typedef std::tuple<NodeType /* parent */,
                       NodeType /* child */,
                       char /* edge label */,
                       size_t /* index in query seq */> EdgeDescriptor;
    VectorMap<AnnotatedDBG::row_index, std::vector<EdgeDescriptor>> anno_rows_to_id;

    size_t k = this->graph_.get_k();
    auto push_path = [&](std::string_view seq, const std::vector<NodeType> &path) {
        process_seq_path(this->graph_, seq, path,
                         [&](AnnotatedDBG::row_index row, size_t i) {
            if (seq[i + k - 1] != boss::BOSS::kSentinel) {
                anno_rows_to_id[row].emplace_back(
                    i ? path[i - 1] : std::get<0>(node),
                    path[i], seq[i + k - 1], i
                );
            }
        });
    };

    // prefetch the next unitig
    call_hull_sequences(this->graph_, std::get<0>(node),
        [&](std::string_view seq, const std::vector<NodeType> &path) {
            assert(path.size());
            assert(this->graph_.traverse(std::get<0>(node),
                                         seq[this->graph_.get_k() - 1]) == path[0]);
            push_path(seq, path);
        },
        [&](std::string_view seq, const auto &path, size_t depth, size_t fork_count) {
            bool result = fork_count || depth > this->query_.size()
                || (path.size() && cached_edge_sets_.count(path.back()));

            if (result && depth == 1)
                push_path(seq, path);

            return result;
        }
    );

    std::vector<AnnotatedDBG::row_index> anno_rows;
    anno_rows.reserve(anno_rows_to_id.size());
    for (const auto &pair : anno_rows_to_id) {
        anno_rows.push_back(pair.first);
    }

    auto rows_with_target = anno_graph_.get_annotation().get_matrix().has_column(
        anno_rows, target_column
    );

    anno_rows = std::vector<AnnotatedDBG::row_index>();

    Edges edges;
    for (AnnotatedDBG::row_index row : rows_with_target) {
        for (const auto &[parent_node, child_node, c, i] : anno_rows_to_id[row]) {
            assert(c != boss::BOSS::kSentinel);
            assert(this->graph_.traverse(parent_node, c) == child_node);

            if (!i) {
                assert(parent_node == std::get<0>(node));
                edges.emplace_back(child_node, c);
            } else {
                auto &parent_edge_sets = cached_edge_sets_[parent_node];
                if (target_column_idx >= parent_edge_sets.size())
                    parent_edge_sets.resize(target_column_idx + 1);

                cached_edge_sets_[parent_node][target_column_idx].emplace_back(
                    child_node, c
                );
            }
        }
    }

    auto &edge_sets = cached_edge_sets_[std::get<0>(node)];
    assert(target_column_idx >= edge_sets.size());
    edge_sets.resize(target_column_idx + 1);
    edge_sets[target_column_idx] = edges;
    return edges;
}

template class LabeledColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg