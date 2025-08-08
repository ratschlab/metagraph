#include "annotation_buffer.hpp"

#include <sstream>

#include "aln_match.hpp"

#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/algorithms.hpp"

namespace mtg {
namespace graph {
namespace align_redone {

using mtg::common::logger;

typedef annot::matrix::BinaryMatrix::Row Row;
typedef annot::matrix::BinaryMatrix::Column Column;

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

    for (const auto &path : queued_paths_) {
        assert(std::find(path.begin(), path.end(), DeBruijnGraph::npos) == path.end());
        assert(path == map_to_nodes_sequentially(graph_, spell_path(graph_, path)));
        std::vector<node_index> base_path;
        if (base_graph->get_mode() == DeBruijnGraph::BASIC) {
            base_path = path;
        } else if (canonical_) {
            base_path.reserve(path.size());
            std::transform(path.begin(), path.end(), std::back_inserter(base_path),
                           [&](node_index node) { return canonical_->get_base_node(node); });
        } else {
            base_path = map_to_nodes(graph_, spell_path(graph_, path));
        }
        assert(base_path.size() == path.size());
        assert(std::find(base_path.begin(), base_path.end(), DeBruijnGraph::npos)
                == base_path.end());

        for (size_t i = 0; i < path.size(); ++i) {
            Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
            auto [it, inserted] = node_to_cols_.try_emplace(base_path[i], nannot);
            if (inserted) {
                assert(path[i] == base_path[i] || !node_to_cols_.count(path[i]));
                if (path[i] != base_path[i])
                    node_to_cols_.emplace(path[i], nannot);

                queued_rows.push_back(row);
                queued_nodes.emplace_back(path[i]);
                if (has_coordinates()) {
                    label_coords_.emplace_back();
                    if (path[i] != base_path[i])
                        label_coords_.emplace_back();
                }
            } else if (base_path[i] != path[i] && it->second != nannot) {
                auto [jt, jnserted] = node_to_cols_.try_emplace(path[i], it->second);
                if (jnserted && has_coordinates())
                    label_coords_.emplace_back(label_coords_[it - node_to_cols_.begin()]);
            }
        }
    }

    queued_paths_.clear();

    if (queued_nodes.empty())
        return;

    auto push_node_labels = [&](auto node_it, auto row_it, auto&& labels, CoordinateSet *coords) {
        assert(node_it != queued_nodes.end());
        assert(row_it != queued_rows.end());
        node_index node = *node_it;
        node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);

        assert(node_to_cols_.count(node));
        assert(node_to_cols_.count(base_node));

        assert(node_to_cols_[node] == nannot);
        assert(node_to_cols_[base_node] == nannot);

        label_class_t label_i = cache_column_set(std::move(labels));
        auto it = node_to_cols_.find(base_node);
        it.value() = label_i;
        if (has_coordinates()) {
            assert(coords);
            label_coords_[it - node_to_cols_.begin()] = *coords;
        }

        if (node != base_node) {
            auto jt = node_to_cols_.find(node);
            jt.value() = label_i;
            if (has_coordinates()) {
                assert(coords);
                label_coords_[jt - node_to_cols_.begin()] = std::move(*coords);
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

            CoordinateSet all_coords;
            all_coords.reserve(row_tuples.size());

            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                all_coords.emplace_back(coords.begin(), coords.end());
            }
            push_node_labels(node_it++, row_it++, std::move(labels), &all_coords);
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
            std::sort(labels.begin(), labels.end());
            push_node_labels(node_it++, row_it++, std::move(labels), nullptr);
        }
    }

#ifndef NDEBUG
    for (const auto &[node, val] : node_to_cols_) {
        assert(val != nannot);
    }
#endif
}

auto AnnotationBuffer::get_labels_id_and_coords(node_index node) const
        -> std::pair<label_class_t, const CoordinateSet*> {
    std::pair<label_class_t, const CoordinateSet*> ret_val { nannot, nullptr };

    if (canonical_)
        node = canonical_->get_base_node(node);

    auto it = node_to_cols_.find(node);

    // if the node hasn't been seen before, or if its annotations haven't
    // been fetched, return nothing
    if (it == node_to_cols_.end() || it->second == nannot)
        return ret_val;

    ret_val.first = it->second;

    if (has_coordinates()) {
        assert(static_cast<size_t>(it - node_to_cols_.begin()) < label_coords_.size());
        ret_val.second = &label_coords_[it - node_to_cols_.begin()];
    }

    return ret_val;
}

auto AnnotationBuffer::get_labels_and_coords(node_index node) const
        -> std::pair<const Columns*, const CoordinateSet*> {
    auto [labels_i, coords] = get_labels_id_and_coords(node);
    return std::make_pair(
        labels_i < column_sets_.size() ? &column_sets_.data()[labels_i] : nullptr,
        coords
    );
}

bool AnnotationBuffer::node_has_labels(node_index node, label_class_t target) const {
    assert(target != nannot);

    auto node_full_columns = get_labels(node);
    assert(node_full_columns);

    const Columns &target_columns = get_cached_column_set(target);
    return utils::count_intersection(node_full_columns->begin(), node_full_columns->end(),
                                     target_columns.begin(), target_columns.end()) == target_columns.size();
}

auto AnnotationBuffer
::get_label_path_coords(const std::vector<node_index> &path,
                        const std::vector<label_class_t> &label_path) const
        -> CoordinateSet {
    assert(path.size() == label_path.size());

    if (!has_coordinates() || path.empty())
        return {};

    label_class_t source = label_path[0];
    assert(source);
    assert(source != nannot);
    assert(std::all_of(label_path.begin() + 1,
                       label_path.end(),
                       [&](label_class_t t) { return t == source; }));
    assert(node_has_labels(path[0], source));

    CoordinateSet ret_val;

    auto [source_full_columns, source_full_coords] = get_labels_and_coords(path[0]);
    assert(source_full_columns);
    assert(source_full_coords);
    assert(source_full_coords->size() == source_full_columns->size());

    auto [target_full_columns, target_full_coords] = get_labels_and_coords(path.back());
    assert(target_full_columns);
    assert(target_full_coords);
    assert(target_full_coords->size() == target_full_columns->size());

    const Columns &source_columns = get_cached_column_set(source);
    assert(source_columns.size());
    auto it_s = source_columns.begin();

    size_t coord_dist = path.size() - 1;
    utils::match_indexed_values(
        source_full_columns->begin(), source_full_columns->end(), source_full_coords->begin(),
        target_full_columns->begin(), target_full_columns->end(), target_full_coords->begin(),
        [&](Column c, const Tuple &source_c, const Tuple &target_c) {
            assert(source_c.size());
            assert(target_c.size());
            if (it_s != source_columns.end() && c == *it_s) {
                Tuple &out = ret_val.emplace_back();
                utils::set_intersection(source_c.begin(), source_c.end(),
                                        target_c.begin(), target_c.end(),
                                        std::back_inserter(out),
                                        coord_dist);
                assert(out.size());
                ++it_s;
            }
        }
    );

    assert(it_s == source_columns.end());

    return ret_val;
}

std::string AnnotationBuffer::generate_column_set_str(label_class_t i, size_t spelling_size) const {
    if (i == nannot)
        return "*";

    std::ostringstream ostd;
    const Columns &columns = get_cached_column_set(i);
    if (columns.empty())
        return "*";

    assert(columns[0] != ncolumn);
    assert(columns[0] < annotator_.get_label_encoder().size());
    ostd << annotator_.get_label_encoder().decode(columns[0]);
    if (label_coords_.size()) {
        assert(label_coords_.size() == node_to_cols_.size());
        for (auto coord : label_coords_[i][0]) {
            ostd << ":" << coord + 1 << "-" << coord + spelling_size;
        }
    }

    for (Column c = 1; c < columns.size(); ++c) {
        assert(columns[c] != ncolumn);
        assert(columns[c] < annotator_.get_label_encoder().size());
        ostd << ";" << annotator_.get_label_encoder().decode(columns[c]);
        if (label_coords_.size()) {
            assert(label_coords_.size() == node_to_cols_.size());
            for (auto coord : label_coords_[i][c]) {
                ostd << ":" << coord + 1 << "-" << coord + spelling_size;
            }
        }
    }

    return ostd.str();
}

} // namespace align
} // namespace graph
} // namespace mtg
