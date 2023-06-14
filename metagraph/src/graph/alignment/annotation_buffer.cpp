#include "annotation_buffer.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"

namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef annot::binmat::BinaryMatrix::Row Row;
typedef annot::binmat::BinaryMatrix::Column Column;

// dummy index for an unfetched annotations
static constexpr size_t nannot = std::numeric_limits<size_t>::max();

AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&annotator_.get_matrix())),
        canonical_(dynamic_cast<const CanonicalDBG*>(&graph_)),
        column_sets_({ {} }) {
    if (multi_int_ && graph_.get_mode() != DeBruijnGraph::BASIC) {
        multi_int_ = nullptr;
        logger->warn("Coordinates not supported when aligning to CANONICAL "
                     "or PRIMARY mode graphs");
    }
}

void AnnotationBuffer::fetch_queued_annotations() {
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    std::vector<node_index> queued_nodes;
    std::vector<Row> queued_rows;

    const DeBruijnGraph *base_graph = &graph_;

    if (canonical_)
        base_graph = &canonical_->get_graph();

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    tsl::hopscotch_map<node_index, node_index> dummy_nodes;

    auto get_base_path = [&](const std::vector<node_index> &path) {
        if (base_graph->get_mode() == DeBruijnGraph::CANONICAL) {
            // TODO: avoid this call of spell_path
            std::string query = spell_path(graph_, path);
            return map_to_nodes(*base_graph, query);
        }

        std::vector<node_index> base_path;
        if (canonical_) {
            base_path.reserve(path.size());
            for (node_index node : path) {
                base_path.emplace_back(canonical_->get_base_node(node));
            }

        } else {
            assert(graph_.get_mode() == DeBruijnGraph::BASIC);
            base_path = path;
            if (dynamic_cast<const RCDBG*>(&graph_))
                std::reverse(base_path.begin(), base_path.end());
        }

        return base_path;
    };

    auto queue_node = [&](node_index node, node_index base_node, bool queue_dummy = true) {
        if (base_node == DeBruijnGraph::npos) {
            // this can happen when the base graph is CANONICAL and path[i] is a
            // dummy node
            dummy_nodes.emplace(node, base_node);
            return;
        }

        if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_node))) {
            // skip dummy nodes
            dummy_nodes.emplace(node, base_node);
            return;
        }

        Row row = AnnotatedDBG::graph_to_anno_index(base_node);
        if (canonical_ || graph_.get_mode() == DeBruijnGraph::BASIC) {
            if (node_to_cols_.try_emplace(base_node, nannot).second) {
                queued_rows.push_back(row);
                queued_nodes.push_back(base_node);
            }

            return;
        }

        assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);

        auto find_a = node_to_cols_.find(node);
        auto find_b = node_to_cols_.find(base_node);

        if (find_a == node_to_cols_.end() && find_b == node_to_cols_.end()) {
            node_to_cols_.try_emplace(node, nannot);
            queued_rows.push_back(row);
            queued_nodes.push_back(node);

            if (node != base_node) {
                node_to_cols_.emplace(base_node, nannot);
                queued_rows.push_back(row);
                queued_nodes.push_back(base_node);
            }
        } else if (find_a == node_to_cols_.end() && find_b != node_to_cols_.end()) {
            node_to_cols_.try_emplace(node, find_b->second);
            if (find_b->second == nannot) {
                queued_rows.push_back(row);
                queued_nodes.push_back(node);
            }
        } else if (find_a != node_to_cols_.end() && find_b == node_to_cols_.end()) {
            node_to_cols_.try_emplace(base_node, find_a->second);
        } else {
            size_t label_i = std::min(find_a->second, find_b->second);
            if (label_i != nannot) {
                find_a.value() = label_i;
                find_b.value() = label_i;
            }
        }
    };

    for (const auto &path : queued_paths_) {
        std::vector<node_index> base_path = get_base_path(path);
        assert(base_path.size() == path.size());

        for (size_t i = 0; i < path.size(); ++i) {
            queue_node(path[i], base_path[i]);
        }
    }

    tsl::hopscotch_map<node_index, tsl::hopscotch_map<node_index, std::vector<size_t>>> annotated_to_dummy_nodes;
    for (const auto &[node, base_node] : dummy_nodes) {
        std::vector<std::pair<node_index, std::string>> traversal;
        traversal.emplace_back(node, graph_.get_node_sequence(node));
        assert(traversal.back().second[0] == boss::BOSS::kSentinel);
        while (traversal.size()) {
            auto [cur_node, spelling] = std::move(traversal.back());
            traversal.pop_back();
            if (*(spelling.rbegin() + graph_.get_k() - 1) != boss::BOSS::kSentinel) {
                if (!annotated_to_dummy_nodes[node].count(cur_node)) {
                    annotated_to_dummy_nodes[node][cur_node].emplace_back(spelling.size() - graph_.get_k());
                    node_index base_node = get_base_path({ cur_node })[0];
                    assert(base_node);
                    queue_node(cur_node, base_node);
                }

                continue;
            }

            spelling.push_back(boss::BOSS::kSentinel);
            graph_.call_outgoing_kmers(cur_node, [&](node_index next, char c) {
                auto &[_, next_spelling] = traversal.emplace_back(next, spelling);
                next_spelling.back() = c;
            });
        }
    }

    dummy_nodes.clear();

    queued_paths_.clear();

    if (queued_nodes.empty())
        return;

    auto push_node_labels = [&](node_index node, auto row, auto&& labels, const CoordinateSet &coords = CoordinateSet{}) {
        assert(node_to_cols_.count(node));
        assert(node_to_cols_.count(AnnotatedDBG::anno_to_graph_index()));

        size_t label_i = cache_column_set(std::move(labels));
        node_index base_node = AnnotatedDBG::anno_to_graph_index(row);
        if (graph_.get_mode() == DeBruijnGraph::BASIC) {
            assert(base_node == node);
            node_to_cols_[node] = label_i;
        } else if (canonical_) {
            node_to_cols_[base_node] = label_i;
        } else {
            node_to_cols_[node] = label_i;
            if (base_node != node && node_to_cols_.try_emplace(base_node, label_i).second
                    && has_coordinates()) {
                label_coords_.emplace_back(coords);
            }
        }
    };

    auto node_it = queued_nodes.begin();
    auto row_it = queued_rows.begin();
    if (has_coordinates()) {
        assert(multi_int_);
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(queued_rows)) {
            std::sort(row_tuples.begin(), row_tuples.end(), utils::LessFirst());
            Columns labels;
            CoordinateSet coords;
            labels.reserve(row_tuples.size());
            coords.reserve(row_tuples.size());
            for (auto&& [label, cur_coords] : row_tuples) {
                labels.push_back(label);
                coords.emplace_back(cur_coords.begin(), cur_coords.end());
            }
            assert(node_it != queued_nodes.end());
            auto original_dummy_nodes = annotated_to_dummy_nodes.find(*node_it);
            if (original_dummy_nodes != annotated_to_dummy_nodes.end()) {
                label_coords_.emplace_back(coords);
                push_node_labels(*node_it, *row_it, decltype(labels)(labels), coords);
                for (auto &[dummy_node, dists] : original_dummy_nodes.value()) {
                    assert(dummy_nodes.count(dummy_node));
                    node_index base_dummy_node = dummy_nodes[dummy_node];

                    CoordinateSet cur_coords;
                    for (auto &coord_set : coords) {
                        auto &cur_coord_set = cur_coords.emplace_back();
                        for (ssize_t d : dists) {
                            for (auto coord : coord_set) {
                                cur_coord_set.emplace_back(coord - d);
                            }
                        }
                        std::sort(cur_coord_set.begin(), cur_coord_set.end());
                        cur_coord_set.erase(std::unique(cur_coord_set.begin(),
                                                        cur_coord_set.end()),
                                            cur_coord_set.end());
                    }

                    label_coords_.emplace_back(cur_coords);
                    push_node_labels(dummy_node,
                        AnnotatedDBG::graph_to_anno_index(base_dummy_node ? base_dummy_node : dummy_node),
                        decltype(labels)(labels),
                        cur_coords
                    );
                }
            } else {
                label_coords_.emplace_back(std::move(coords));
                push_node_labels(*node_it, *row_it, std::move(labels), label_coords_.back());
            }
            ++node_it;
            ++row_it;
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
            std::sort(labels.begin(), labels.end());
            assert(node_it != queued_nodes.end());
            auto original_dummy_nodes = annotated_to_dummy_nodes.find(*node_it);
            if (original_dummy_nodes != annotated_to_dummy_nodes.end()) {
                push_node_labels(*node_it, *row_it, decltype(labels)(labels));
                for (const auto &[dummy_node, dists] : original_dummy_nodes->second) {
                    assert(dummy_nodes.count(dummy_node));
                    node_index base_dummy_node = dummy_nodes[dummy_node];
                    push_node_labels(dummy_node,
                        AnnotatedDBG::graph_to_anno_index(base_dummy_node ? base_dummy_node : dummy_node),
                        decltype(labels)(labels)
                    );
                }
            } else {
                push_node_labels(*node_it, *row_it, std::move(labels));
            }
            ++node_it;
            ++row_it;
        }
    }

#ifndef NDEBUG
    for (const auto &[node, val] : node_to_cols_) {
        assert(val != nannot);
    }
#endif
}

auto AnnotationBuffer::get_labels_and_coords(node_index node) const
        -> std::pair<const Columns*, const CoordinateSet*> {
    std::pair<const Columns*, const CoordinateSet*> ret_val { nullptr, nullptr };

    if (canonical_)
        node = canonical_->get_base_node(node);

    auto it = node_to_cols_.find(node);

    // if the node hasn't been seen before, or if its annotations haven't
    // been fetched, return nothing
    if (it == node_to_cols_.end() || it->second == nannot)
        return ret_val;

    ret_val.first = &column_sets_.data()[it->second];

    if (has_coordinates()) {
        assert(static_cast<size_t>(it - node_to_cols_.begin()) < label_coords_.size());
        ret_val.second = &label_coords_[it - node_to_cols_.begin()];
    }

    return ret_val;
}

} // namespace align
} // namespace graph
} // namespace mtg
