#include "annotation_buffer.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/algorithms.hpp"

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

    auto queue_node = [&](node_index node, node_index base_node) {
        if (node_to_cols_.count(node))
            return;

        if (base_node == DeBruijnGraph::npos) {
            // this can happen when the base graph is CANONICAL and path[i] is a
            // dummy node
            dummy_nodes.emplace(node, node);
            return;
        }

        if (boss && (!boss->get_W(dbg_succ->kmer_to_boss_index(base_node))
                    || boss->is_dummy(dbg_succ->kmer_to_boss_index(base_node)))) {
            // skip dummy nodes
            dummy_nodes.emplace(node, base_node);
            return;
        }

        assert(!boss || dbg_succ->get_node_sequence(base_node).find(boss::BOSS::kSentinel) == std::string::npos);

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

    using NodeToDist = tsl::hopscotch_map<node_index, std::vector<size_t>>;
    VectorMap<node_index,
              std::pair<node_index, NodeToDist>> dummy_to_annotated_node;
    for (const auto &[node, base_node] : dummy_nodes) {
        assert(boss);
        assert(base_node);
        assert(!node_to_cols_.count(node));
        assert(!node_to_cols_.count(base_node));

        std::vector<std::pair<node_index, std::string>> traversal;
        traversal.emplace_back(node, graph_.get_node_sequence(node));
        assert(traversal.back().second[0] == boss::BOSS::kSentinel);
        bool discovered = false;
        while (traversal.size()) {
            auto [cur_node, spelling] = std::move(traversal.back());
            traversal.pop_back();

            if (node_to_cols_.count(cur_node)
                    || *(spelling.rbegin() + graph_.get_k() - 1) != boss::BOSS::kSentinel) {
                discovered = true;
                assert(spelling.size() > graph_.get_k());
                auto &mapping = dummy_to_annotated_node.try_emplace(
                    node, std::make_pair(base_node, NodeToDist{})
                ).first.value().second;
                mapping[cur_node].emplace_back(spelling.size() - graph_.get_k());
                node_index cur_base_node = get_base_path({ cur_node })[0];
                assert(cur_base_node);
                assert(!boss->is_dummy(dbg_succ->kmer_to_boss_index(cur_base_node)));
                queue_node(cur_node, cur_base_node);
                continue;
            }

            spelling.push_back(boss::BOSS::kSentinel);
            graph_.call_outgoing_kmers(cur_node, [&,s=std::move(spelling)](node_index next, char c) {
                if (c != boss::BOSS::kSentinel) {
                    auto &[_, next_spelling] = traversal.emplace_back(next, s);
                    next_spelling.back() = c;
                }
            });
        }

        assert(discovered);

        if (base_node != node)
            node_to_cols_.try_emplace(base_node, nannot);

        node_to_cols_.try_emplace(node, nannot);
    }

    dummy_nodes.clear();
    queued_paths_.clear();

    auto push_node_labels = [&](node_index node, auto row, auto&& labels, const CoordinateSet &coords = CoordinateSet{}) {
        node_index base_node = AnnotatedDBG::anno_to_graph_index(row);

        auto node_find = node_to_cols_.find(node);
        auto base_node_find = node_to_cols_.find(base_node);
        assert(node_find != node_to_cols_.end());
        assert(base_node_find != node_to_cols_.end());

        if (has_coordinates()) {
            assert(node_to_cols_.begin() + label_coords_.size() == node_find);
            label_coords_.emplace_back(coords);
        }

        size_t label_i = cache_column_set(std::move(labels));
        if (graph_.get_mode() == DeBruijnGraph::BASIC) {
            assert(base_node == node);
            node_find.value() = label_i;
        } else if (canonical_) {
            node_find.value() = label_i;
            base_node_find.value() = label_i;
        } else {
            node_find.value() = label_i;
            if (base_node != node && base_node_find.value() != label_i) {
                assert(base_node_find->second == nannot);
                base_node_find.value() = label_i;
                if (has_coordinates()) {
                    assert(node_to_cols_.begin() + label_coords_.size() == base_node_find);
                    label_coords_.emplace_back(coords);
                }
            }
        }

        assert(node_find->second != nannot);
        assert(base_node_find->second != nannot);
    };

    if (queued_nodes.size()) {
        auto node_it = queued_nodes.begin();
        auto row_it = queued_rows.begin();
        if (has_coordinates()) {
            assert(multi_int_);
            // extract both labels and coordinates, then store them separately
            for (auto&& row_tuples : multi_int_->get_row_tuples(queued_rows)) {
                assert(row_tuples.size());
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
                push_node_labels(*node_it, *row_it, std::move(labels), coords);
                ++node_it;
                ++row_it;
            }
        } else {
            for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
                assert(labels.size());
                std::sort(labels.begin(), labels.end());
                assert(node_it != queued_nodes.end());
                push_node_labels(*node_it, *row_it, std::move(labels));
                ++node_it;
                ++row_it;
            }
        }
    }

    for (const auto &[dummy_node, mapping_pair] : dummy_to_annotated_node) {
        Columns labels;
        CoordinateSet coords;

        assert(mapping_pair.first != DeBruijnGraph::npos);
        auto row = AnnotatedDBG::graph_to_anno_index(mapping_pair.first);

        const auto &mapping = mapping_pair.second;
        assert(mapping.size());

        for (const auto &[annotated_node, dists] : mapping) {
            auto [cur_labels, cur_coords] = get_labels_and_coords(annotated_node);
            assert(cur_labels);
            assert(!has_coordinates() || cur_coords);
            Columns union_labels;
            if (cur_coords) {
                CoordinateSet union_coords;
                utils::match_indexed_values(labels.begin(), labels.end(), coords.begin(),
                                            cur_labels->begin(), cur_labels->end(), cur_coords->begin(),
                                            [&](const auto label, const auto &c1, const auto &c2) {
                    union_labels.emplace_back(label);
                    auto &merge_coords = union_coords.emplace_back();
                    for (ssize_t d : dists) {
                        utils::set_union(c2.begin(), c2.end(), c1.begin(), c1.end(),
                                         std::back_inserter(merge_coords), -d);
                    }
                });
                std::swap(union_coords, coords);
            } else {
                std::set_union(labels.begin(), labels.end(), cur_labels->begin(), cur_labels->end(),
                               std::back_inserter(union_labels));
            }
            std::swap(union_labels, labels);
        }

        push_node_labels(dummy_node, row, std::move(labels), coords);
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
