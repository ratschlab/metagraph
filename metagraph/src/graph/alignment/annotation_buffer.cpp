#include "annotation_buffer.hpp"
#include "graph/annotation_key_mapper.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"

#include <map>

namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef annot::matrix::BinaryMatrix::Row Row;
typedef annot::matrix::BinaryMatrix::Column Column;

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

    const AnnotationKeyMapper key_mapper(graph_);

    for (const auto &path : queued_paths_) {
        const std::string sequence = key_mapper.needs_sequence() ? spell_path(graph_, path) : "";
        auto base_path = key_mapper.keys(path, sequence);
        if (graph_.get_mode() == DeBruijnGraph::BASIC && dynamic_cast<const RCDBG*>(&graph_))
            std::reverse(base_path.begin(), base_path.end());

        assert(base_path.size() == path.size());

        for (size_t i = 0; i < path.size(); ++i) {
            if (base_path[i] == DeBruijnGraph::npos) {
                // this can happen when the base graph is CANONICAL and path[i] is a
                // dummy node
                if (node_to_cols_.try_emplace(path[i], 0).second && has_coordinates())
                    label_coords_.emplace_back();

                continue;
            }

            if (boss && !boss->get_W(base_path[i])) {
                // skip dummy nodes
                if (node_to_cols_.try_emplace(base_path[i], 0).second && has_coordinates())
                    label_coords_.emplace_back();

                if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                        && base_path[i] != path[i]
                        && node_to_cols_.emplace(path[i], 0).second && has_coordinates()) {
                    label_coords_.emplace_back();
                }

                continue;
            }

            Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
            if (canonical_ || graph_.get_mode() == DeBruijnGraph::BASIC) {
                if (node_to_cols_.try_emplace(base_path[i], nannot).second) {
                    queued_rows.push_back(row);
                    queued_nodes.push_back(base_path[i]);
                }

                continue;
            }

            assert(graph_.get_mode() == DeBruijnGraph::CANONICAL);

            auto find_a = node_to_cols_.find(path[i]);
            auto find_b = node_to_cols_.find(base_path[i]);

            if (find_a == node_to_cols_.end() && find_b == node_to_cols_.end()) {
                node_to_cols_.try_emplace(path[i], nannot);
                queued_rows.push_back(row);
                queued_nodes.push_back(path[i]);

                if (path[i] != base_path[i]) {
                    node_to_cols_.emplace(base_path[i], nannot);
                    queued_rows.push_back(row);
                    queued_nodes.push_back(base_path[i]);
                }
            } else if (find_a == node_to_cols_.end() && find_b != node_to_cols_.end()) {
                // copy the value first: inserting may reallocate the map and
                // invalidate find_b
                size_t label_i = find_b->second;
                node_to_cols_.try_emplace(path[i], label_i);
                if (label_i == nannot) {
                    queued_rows.push_back(row);
                    queued_nodes.push_back(path[i]);
                }
            } else if (find_a != node_to_cols_.end() && find_b == node_to_cols_.end()) {
                size_t label_i = find_a->second;
                node_to_cols_.try_emplace(base_path[i], label_i);
            } else {
                size_t label_i = std::min(find_a->second, find_b->second);
                if (label_i != nannot) {
                    find_a.value() = label_i;
                    find_b.value() = label_i;
                }
            }
        }
    }

    queued_paths_.clear();

    if (queued_nodes.empty())
        return;

    auto push_node_labels = [&](auto node_it, auto row_it, auto&& labels) {
        assert(node_it != queued_nodes.end());
        assert(node_to_cols_.count(*node_it));
        assert(node_to_cols_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));

        size_t label_i = cache_column_set(std::move(labels));
        node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);
        if (graph_.get_mode() == DeBruijnGraph::BASIC) {
            assert(base_node == *node_it);
            node_to_cols_[*node_it] = label_i;
        } else if (canonical_) {
            node_to_cols_[base_node] = label_i;
        } else {
            node_to_cols_[*node_it] = label_i;
            if (base_node != *node_it && node_to_cols_.try_emplace(base_node, label_i).second
                    && has_coordinates()) {
                label_coords_.emplace_back(label_coords_.back());
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
            labels.reserve(row_tuples.size());
            label_coords_.emplace_back();
            label_coords_.back().reserve(row_tuples.size());
            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                label_coords_.back().emplace_back(coords.begin(), coords.end());
            }
            push_node_labels(node_it++, row_it++, std::move(labels));
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
            std::sort(labels.begin(), labels.end());
            push_node_labels(node_it++, row_it++, Columns(labels.begin(), labels.end()));
        }
    }

    if (dbg_succ && graph_.get_mode() == DeBruijnGraph::BASIC)
        resolve_source_dummy_labels(*dbg_succ, queued_nodes);

#ifndef NDEBUG
    for (const auto &[node, val] : node_to_cols_) {
        assert(val != nannot);
    }
#endif
}

void AnnotationBuffer::resolve_source_dummy_labels(const DBGSuccinct &dbg_succ,
                                                   const std::vector<node_index> &nodes) {
    // A source dummy k-mer $...$T has no annotation, but a seed shorter than k
    // which occurs only within the first k-1 characters of an input sequence
    // matches only such a k-mer. Give it the labels of the real k-mers that
    // continue it (they all start with T), with their coordinates shifted back
    // to the start of the dummy k-mer.
    for (node_index node : nodes) {
        auto find = node_to_cols_.find(node);
        assert(find != node_to_cols_.end());
        if (find->second)
            continue;  // has labels, so it is a real k-mer

        std::string spelling = dbg_succ.get_node_sequence(node);
        size_t num_dummy = std::find_if(spelling.begin(), spelling.end(),
                                        [](char c) { return c != boss::BOSS::kSentinel; })
                                - spelling.begin();
        if (!num_dummy || num_dummy == spelling.size())
            continue;

        // walk |num_dummy| steps forward to reach the real k-mers
        std::vector<node_index> frontier { node };
        for (size_t step = 0; step < num_dummy && frontier.size(); ++step) {
            std::vector<node_index> next;
            for (node_index n : frontier) {
                dbg_succ.adjacent_outgoing_nodes(n, [&](node_index out) {
                    next.push_back(out);
                });
            }
            frontier = next.size() <= kMaxSourceDummyFanout
                ? std::move(next)
                : std::vector<node_index>{};
        }

        if (frontier.empty())
            continue;

        std::vector<Row> rows;
        rows.reserve(frontier.size());
        for (node_index n : frontier) {
            rows.push_back(AnnotatedDBG::graph_to_anno_index(n));
        }

        // union of the labels (and of their coordinates) of the real k-mers
        std::map<Column, std::vector<int64_t>> label_coords;
        if (has_coordinates()) {
            for (auto&& row_tuples : multi_int_->get_row_tuples(rows)) {
                for (auto&& [label, coords] : row_tuples) {
                    auto &out = label_coords[label];
                    for (int64_t coord : coords) {
                        out.push_back(coord - num_dummy);
                    }
                }
            }
        } else {
            for (auto&& labels : annotator_.get_matrix().get_rows(rows)) {
                for (Column label : labels) {
                    label_coords[label];
                }
            }
        }

        if (label_coords.empty())
            continue;

        Columns labels;
        labels.reserve(label_coords.size());
        for (const auto &[label, coords] : label_coords) {
            labels.push_back(label);
        }
        find.value() = cache_column_set(std::move(labels));

        if (has_coordinates()) {
            auto &coord_set = label_coords_[find - node_to_cols_.begin()];
            coord_set.clear();
            for (auto &[label, coords] : label_coords) {
                std::sort(coords.begin(), coords.end());
                coords.erase(std::unique(coords.begin(), coords.end()), coords.end());
                coord_set.emplace_back(coords.begin(), coords.end());
            }
        }
    }
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

    ret_val.first = &get_cached_column_set(it->second);

    if (has_coordinates()) {
        assert(static_cast<size_t>(it - node_to_cols_.begin()) < label_coords_.size());
        ret_val.second = &label_coords_[it - node_to_cols_.begin()];
    }

    return ret_val;
}

} // namespace align
} // namespace graph
} // namespace mtg
