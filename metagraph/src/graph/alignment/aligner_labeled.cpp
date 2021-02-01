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
                callback(
                    AnnotatedDBG::graph_to_anno_index(canonical->get_base_node(query_nodes[i])),
                    i
                );
            }
        }
    } else if (!graph.is_canonical_mode()) {
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
::map_query_batch(const QueryGenerator &generate_query) const
        -> std::pair<BatchMapping, BatchTargets> {
    BatchMapping query_nodes;
    VectorMap<AnnotatedDBG::row_index,
              std::vector<std::pair<size_t, size_t>>> row_to_query_idx;

    size_t num_queries = 0;
    generate_query([&](std::string_view, std::string_view query, bool) {
        query_nodes.emplace_back(map_sequence_to_nodes(graph_, query));
        process_seq_path(graph_, query, query_nodes.back(),
                         [&](AnnotatedDBG::row_index row, size_t i) {
            row_to_query_idx[row].emplace_back(num_queries, i);
        });
        ++num_queries;
    });

    std::vector<AnnotatedDBG::row_index> rows;
    rows.reserve(row_to_query_idx.size());
    for (const auto &[row, mapping] : row_to_query_idx) {
        rows.push_back(row);
    }

    auto annotation = anno_graph_.get_annotation().get_matrix().get_rows(rows);

    // count labels for each query
    std::vector<VectorMap<uint64_t, uint64_t>> column_counter(num_queries);
    for (size_t i = 0; i < annotation.size(); ++i) {
        AnnotatedDBG::row_index row = rows[i];
        for (const auto &[query_id, idx] : row_to_query_idx[row]) {
            for (uint64_t column : annotation[i]) {
                ++column_counter[query_id][column];
            }
        }
    }

    // compute target columns and initialize signatures for each query
    BatchTargets target_columns(num_queries);
    for (size_t j = 0; j < column_counter.size(); ++j) {
        auto &counter = const_cast<std::vector<std::pair<uint64_t, uint64_t>>&>(
            column_counter[j].values_container()
        );

        std::sort(counter.begin(), counter.end(), utils::GreaterSecond());
        size_t num_targets = std::min(counter.size(), num_top_labels_);
        for (size_t i = 0; i < num_targets; ++i) {
            if (!target_columns[j].count(counter[i].first)) {
                target_columns[j].emplace(counter[i].first,
                                          sdsl::bit_vector(query_nodes[j].size(), false));
            }
        }
    }

    // fill signatures for each query
    size_t i = 0;
    for (const auto &[row, mapping] : row_to_query_idx) {
        for (const auto &[query_id, idx] : mapping) {
            for (uint64_t column : annotation[i]) {
                if (target_columns[query_id].count(column))
                    target_columns[query_id][column][idx] = true;
            }
        }

        ++i;
    }

    if (std::all_of(target_columns.begin(), target_columns.end(),
                    [](const auto &a) { return a.empty(); })) {
        for (size_t i = 0; i < target_columns.size(); ++i) {
            target_columns[i][kNTarget] = sdsl::bit_vector(query_nodes[i].size(), true);
        }
    }

    return std::make_pair(std::move(query_nodes), std::move(target_columns));
}


template <typename NodeType>
LabeledColumnExtender<NodeType>
::LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                        const DBGAlignerConfig &config,
                        std::string_view query,
                        uint64_t target_column)
      : DefaultColumnExtender<NodeType>(anno_graph.get_graph(), config, query),
        anno_graph_(anno_graph),
        main_target_column_(target_column) {}


template <typename NodeType>
void LabeledColumnExtender<NodeType>::operator()(ExtensionCallback callback,
                                                 score_t min_path_score) {
    DefaultColumnExtender<NodeType>::operator()([&](DBGAlignment&& extension,
                                                    NodeType start_node) {
        if (seed_extension_.size() && start_node == seed_extension_.back()
                && !extension.get_clipping()) {
            DBGAlignment new_path = seed_extension_;
            new_path.append(std::move(extension));
            callback(std::move(new_path), old_start_);
        } else {
            callback(std::move(extension), start_node);
        }
    }, min_path_score);
}

template <typename NodeType>
void LabeledColumnExtender<NodeType>::initialize(const DBGAlignment &path) {
    alt_seed_ = DBGAlignment();
    seed_extension_ = DBGAlignment();
    old_start_ = NodeType();

    DefaultColumnExtender<NodeType>::initialize(path);

    if (!path.get_offset()) {
        target_column_ = main_target_column_;
        return;
    }

    target_column_ = ILabeledDBGAligner::kNTarget;

    // the path is a suffix match
    assert(path.size() == 1);
    std::string_view subquery(this->query.data(), path.get_clipping() + this->graph_.get_k());

    alt_seed_ = path;
    DefaultColumnExtender<> path_extender(this->graph_, this->config_, subquery);
    path_extender.initialize(alt_seed_);

    bool extended = false;
    path_extender([&](DBGAlignment&& rest, NodeType start_node) {
        if (start_node && !extended && !rest.get_clipping() && !rest.get_end_clipping()) {
            extended = true;
            seed_extension_ = DBGAlignment(rest);
            alt_seed_.append(std::move(rest));
            alt_seed_.trim_offset();
            assert(alt_seed_.is_valid(this->graph_, &this->config_));
            assert(!alt_seed_.get_offset());
        }
    });

    if (extended) {
        auto labels = anno_graph_.get_top_labels(alt_seed_.get_sequence(), 1, 1.0);
        if (labels.size()) {
            old_start_ = path.back();
            DefaultColumnExtender<NodeType>::initialize(alt_seed_);
            target_column_ = anno_graph_.get_annotation().get_label_encoder().encode(
                labels[0].first
            );
        } else {
            alt_seed_ = DBGAlignment();
            seed_extension_ = DBGAlignment();
        }
    }
}

template <typename NodeType>
auto LabeledColumnExtender<NodeType>::fork_extension(NodeType node,
                                                     ExtensionCallback callback,
                                                     score_t min_path_score) -> Edges {
    if (target_column_ == ILabeledDBGAligner::kNTarget)
        return DefaultColumnExtender<NodeType>::fork_extension(node, callback, min_path_score);

    if (cached_edge_sets_.count(node))
        return cached_edge_sets_[node];

    typedef std::tuple<NodeType /* parent */,
                       NodeType /* child */,
                       char /* edge label */,
                       size_t /* index in query seq */> EdgeDescriptor;
    VectorMap<AnnotatedDBG::row_index, std::vector<EdgeDescriptor>> anno_rows_to_id;

    size_t k = this->graph_.get_k();
    auto push_path = [&](std::string_view seq, const std::vector<NodeType> &path) {
        process_seq_path(this->graph_, seq, path, [&](AnnotatedDBG::row_index row, size_t i) {
            if (seq[i + k - 1] != boss::BOSS::kSentinel) {
                anno_rows_to_id[row].emplace_back(
                    i ? path[i - 1] : node, path[i], seq[i + k - 1], i
                );
            }
        });
    };

    // prefetch the next unitig
    call_hull_sequences(this->graph_, node,
        [&](std::string_view seq, const std::vector<NodeType> &path) {
            assert(path.size());
            assert(this->graph_.traverse(node, seq[this->graph_.get_k() - 1]) == path[0]);
            push_path(seq, path);
        },
        [&](std::string_view seq, const auto &path, size_t depth, size_t fork_count) {
            bool result = fork_count || depth > this->query.size()
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
        anno_rows, target_column_
    );

    anno_rows = std::vector<AnnotatedDBG::row_index>();

    Edges edges;
    for (AnnotatedDBG::row_index row : rows_with_target) {
        for (const auto &[parent_node, child_node, c, i] : anno_rows_to_id[row]) {
            assert(c != boss::BOSS::kSentinel);
            assert(this->graph_.traverse(parent_node, c) == child_node);

            if (!this->dp_table.count(child_node) && this->ram_limit_reached())
                continue;

            if (!i) {
                assert(parent_node == node);
                edges.emplace_back(child_node, c);
            } else {
                cached_edge_sets_[parent_node].emplace_back(child_node, c);
            }
        }
    }

    cached_edge_sets_[node] = edges;
    return edges;
}

template class LabeledColumnExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
