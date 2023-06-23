#include "annotation_buffer.hpp"

#include <tsl/hopscotch_set.h>

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector_set.hpp"
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

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(
        canonical_ ? &canonical_->get_graph() : &graph_
    );

    VectorSet<Row> queued_rows;
    std::vector<node_index> queued_nodes;
    tsl::hopscotch_set<node_index> dummy_nodes;

    std::vector<std::pair<node_index, node_index>> add_base_annot;

    std::function<void(node_index, node_index)> queue_node
        = [](node_index, node_index) {};
    if (canonical_) {
        queue_node = [&](node_index node, node_index base_node) {
            assert(base_node);
            auto find_base = node_to_cols_.find(base_node);
            auto row = AnnotatedDBG::graph_to_anno_index(base_node);
            if (find_base != node_to_cols_.end()) {
                assert(find_base->second != nannot || queued_rows.count(row));
                return;
            }

            if (dbg_succ && !dbg_succ->get_mask()
                    && dbg_succ->get_boss().is_dummy(base_node)) {
                dummy_nodes.emplace(node);
                return;
            }

            if (queued_rows.emplace(row).second)
                node_to_cols_.emplace(base_node, nannot);
        };
    } else if (graph_.get_mode() == DeBruijnGraph::BASIC) {
        queue_node = [&](node_index node, node_index = 0) {
            assert(node);
            auto find = node_to_cols_.find(node);
            auto row = AnnotatedDBG::graph_to_anno_index(node);
            if (find != node_to_cols_.end()) {
                assert(find->second != nannot || queued_rows.count(row));
                return;
            }

            if (dbg_succ && !dbg_succ->get_mask()
                    && dbg_succ->get_boss().is_dummy(node)) {
                dummy_nodes.emplace(node);
                return;
            }

            if (queued_rows.emplace(row).second)
                node_to_cols_.emplace(node, nannot);
        };
    } else {
        assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);
        queue_node = [&](node_index node, node_index base_node) {
            assert(node);
            if (base_node) {
                auto find_base = node_to_cols_.find(base_node);
                if (find_base == node_to_cols_.end()) {
                    auto row = AnnotatedDBG::graph_to_anno_index(base_node);
                    if (queued_rows.emplace(row).second) {
                        node_to_cols_.emplace(base_node, nannot);
                        node_to_cols_.emplace(node, nannot);
                        queued_nodes.emplace_back(node);
                    }
                } else if (node != base_node) {
                    node_to_cols_.try_emplace(node, find_base->second);
                    if (find_base->second == nannot)
                        add_base_annot.emplace_back(node, base_node);
                }
            } else {
                assert(dbg_succ);
                assert(!dbg_succ->get_mask());
                assert(dbg_succ->get_boss().is_dummy(node));
                dummy_nodes.emplace(node);
            }
        };
    }

    for (const auto &path : queued_paths_) {
        if (canonical_) {
            for (node_index node : path) {
                queue_node(node, canonical_->get_base_node(node));
            }
        } else if (graph_.get_mode() == DeBruijnGraph::BASIC) {
            for (node_index node : path) {
                queue_node(node, node);
            }
        } else {
            // TODO: avoid this spelling
            std::string spelling = spell_path(graph_, path);
            auto it = path.begin();
            for (node_index base_node : map_to_nodes(graph_, spelling)) {
                assert(it != path.end());
                queue_node(*it, base_node);
                ++it;
            }
            assert(it == path.end());
        }
    }

    tsl::hopscotch_set<node_index> annotated_nodes;
    tsl::hopscotch_map<node_index, std::vector<node_index>> parents;
    for (node_index node : dummy_nodes) {
        assert(dbg_succ);
        assert(!node_to_cols_.count(node));

        std::vector<std::pair<node_index, size_t>> traversal;
        std::string spelling = graph_.get_node_sequence(node);
        traversal.emplace_back(node, spelling.find_last_of(boss::BOSS::kSentinel) + 1);
        assert(traversal.back().second < spelling.size());

        while (traversal.size()) {
            auto [cur_node, num_sentinels_left] = std::move(traversal.back());
            traversal.pop_back();

            node_index cur_base_node = canonical_
                ? canonical_->get_base_node(cur_node)
                : cur_node;
            assert(cur_base_node);

            if (node_to_cols_.count(cur_base_node)) {
                assert(canonical_ || node_to_cols_.count(cur_node));
                annotated_nodes.emplace(cur_node);
                continue;
            }

            if (!num_sentinels_left) {
                assert(!dbg_succ->get_boss().is_dummy(
                    dbg_succ->kmer_to_boss_index(cur_base_node)
                ));
                queue_node(cur_node, cur_base_node);
                assert(node_to_cols_.count(cur_base_node));
                annotated_nodes.emplace(cur_node);
                continue;
            }

            node_to_cols_.try_emplace(cur_base_node, nannot);
            if (!canonical_ && cur_node != cur_base_node) {
                assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);
                node_to_cols_.try_emplace(cur_node, nannot);
            }

            --num_sentinels_left;
            graph_.adjacent_outgoing_nodes(cur_node, [&](node_index next) {
                parents[next].emplace_back(cur_node);
                traversal.emplace_back(next, num_sentinels_left);
            });
        }
    }

    dummy_nodes.clear();

    queued_paths_.clear();

    auto push_node_labels = [&](node_index node,
                                auto row,
                                auto&& labels,
                                const CoordinateSet coords = {}) {
        auto do_push = [&](auto find, size_t labels_i) {
            find.value() = labels_i;
            if (has_coordinates()) {
                assert(coords.size());
                size_t coord_idx = find - node_to_cols_.begin();
                if (coord_idx == label_coords_.size()) {
                    label_coords_.emplace_back(coords);
                } else {
                    label_coords_.resize(std::max(label_coords_.size(), coord_idx + 1));
                    label_coords_[coord_idx] = coords;
                }
            }
        };

        node_index base_node = AnnotatedDBG::anno_to_graph_index(row);
        auto find_base = node_to_cols_.find(base_node);
        assert(find_base != node_to_cols_.end());
        size_t labels_i = cache_column_set(std::move(labels));;
        do_push(find_base, labels_i);

        if (canonical_ || graph_.get_mode() == DeBruijnGraph::BASIC || base_node == node)
            return;

        auto find = node_to_cols_.find(node);
        assert(find != node_to_cols_.end());
        assert(find->second == nannot);
        assert(find_base->second != nannot);
        do_push(find, labels_i);
    };

    auto row_it = queued_rows.begin();
    auto node_it = queued_nodes.begin();
    if (has_coordinates()) {
        assert(multi_int_);
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(queued_rows.values_container())) {
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

            assert(row_it != queued_rows.end());
            if (queued_nodes.size()) {
                assert(node_it != queued_nodes.end());
                push_node_labels(*node_it, *row_it, std::move(labels), coords);
                ++node_it;
            } else {
                push_node_labels(AnnotatedDBG::anno_to_graph_index(*row_it),
                                 *row_it, std::move(labels), coords);
            }
            ++row_it;
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows.values_container())) {
            assert(labels.size());
            std::sort(labels.begin(), labels.end());
            if (queued_nodes.size()) {
                assert(!canonical_ && graph_.get_mode() == DeBruijnGraph::CANONICAL);
                assert(node_it != queued_nodes.end());
                push_node_labels(*node_it, *row_it, std::move(labels));
                ++node_it;
            } else {
                push_node_labels(AnnotatedDBG::anno_to_graph_index(*row_it),
                                 *row_it, std::move(labels));
            }
            ++row_it;
        }
    }

    assert(row_it == queued_rows.end());
    assert(node_it == queued_nodes.end());

    for (const auto &[node, base_node] : add_base_annot) {
        auto find_base = node_to_cols_.find(base_node);
        assert(find_base != node_to_cols_.end());
        assert(find_base->second != nannot);

        auto find = node_to_cols_.find(node);
        assert(find != node_to_cols_.end());
        assert(find->second == nannot || find->second == find_base->second);
        find.value() = find_base->second;
    }

    for (node_index node : annotated_nodes) {
        assert(parents.count(node));
        std::vector<node_index> back_traversal;
        back_traversal.emplace_back(node);
        while (back_traversal.size()) {
            node_index node = back_traversal.back();
            back_traversal.pop_back();
            assert(parents.count(node));

            auto [labels, coords] = get_labels_and_coords(node);
            assert(labels);
            assert(labels->size());

            for (node_index prev : parents[node]) {
                node_index base_node = canonical_ ? canonical_->get_base_node(prev) : prev;
                assert(canonical_ || node_to_cols_.count(prev));
                assert(node_to_cols_.count(base_node));
                auto [prev_labels, prev_coords] = get_labels_and_coords(prev);
                CoordinateSet merged_prev_coords;
                if (!prev_labels) {
                    if (has_coordinates()) {
                        assert(coords);
                        merged_prev_coords.reserve(coords->size());
                        for (auto &tuple : *coords) {
                            auto &prev_tuple = merged_prev_coords.emplace_back();
                            prev_tuple.reserve(tuple.size());
                            for (auto c : tuple) {
                                prev_tuple.emplace_back(c - 1);
                            }
                        }
                    }

                    push_node_labels(prev,
                                     AnnotatedDBG::graph_to_anno_index(base_node),
                                     decltype(*labels)(*labels),
                                     merged_prev_coords);
                } else {
                    Columns merged_columns;
                    if (has_coordinates()) {
                        assert(coords);
                        assert(prev_coords);
                        utils::match_indexed_values(labels->begin(), labels->end(),
                                                    coords->begin(),
                                                    prev_labels->begin(), prev_labels->end(),
                                                    prev_coords->begin(),
                                                    [&](const auto label,
                                                        const auto &c1,
                                                        const auto &c2) {
                            merged_columns.emplace_back(label);
                            auto &merge_coords = merged_prev_coords.emplace_back();
                            utils::set_union(c2.begin(), c2.end(), c1.begin(), c1.end(),
                                             std::back_inserter(merge_coords), -1);
                        });
                    } else {
                        std::set_union(labels->begin(), labels->end(),
                                       prev_labels->begin(), prev_labels->end(),
                                       std::back_inserter(merged_columns));
                    }

                    push_node_labels(prev,
                                     AnnotatedDBG::graph_to_anno_index(base_node),
                                     std::move(merged_columns),
                                     merged_prev_coords);
                }

                if (parents.count(prev))
                    back_traversal.emplace_back(prev);
            }
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
