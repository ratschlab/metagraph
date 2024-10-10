#include "annotation_buffer.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/utils/template_utils.hpp"

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
    std::vector<Row> queued_rows;

    const DeBruijnGraph *base_graph = &graph_;

    if (canonical_)
        base_graph = &canonical_->get_graph();

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;
    const DBGSSHash *sshash = dynamic_cast<const DBGSSHash*>(base_graph);

    std::vector<std::tuple<node_index, node_index, int64_t>> to_update;
    for (const auto &path : queued_paths_) {
        std::vector<node_index> base_path;
        std::vector<int64_t> base_path_offsets;
        if (base_graph->get_mode() == DeBruijnGraph::CANONICAL || (sshash && sshash->is_monochromatic())) {
            // TODO: avoid this call of spell_path
            std::string query = spell_path(graph_, path);
            base_path.reserve(path.size());
            call_annotated_nodes_offsets(graph_, query, [&](node_index i, int64_t o) {
                assert(boss || i != DeBruijnGraph::npos);
                base_path.emplace_back(i);
                base_path_offsets.emplace_back(o);
            });
        } else if (canonical_) {
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

        base_path_offsets.resize(base_path.size());

        assert(base_path.size() == path.size());

        for (size_t i = 0; i < path.size(); ++i) {
            if (base_path[i] == DeBruijnGraph::npos
                    || (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i])))) {
                // this can happen when the base graph is CANONICAL and path[i] is a
                // dummy node
                if (node_to_cols_.try_emplace(base_path[i], 0).second && has_coordinates())
                    label_coords_.emplace_back();

                if (node_to_cols_.try_emplace(path[i], 0).second && has_coordinates())
                    label_coords_.emplace_back();

                continue;
            }

            to_update.emplace_back(base_path[i], path[i], base_path_offsets[i]);

            Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
            if (node_to_cols_.try_emplace(base_path[i], nannot).second) {
                queued_rows.push_back(row);
                if (has_coordinates())
                    label_coords_.emplace_back();
            }
        }
    }

    assert(!has_coordinates() || node_to_cols_.size() == label_coords_.size());

    queued_paths_.clear();

    auto row_it = queued_rows.begin();
    if (has_coordinates()) {
        assert(multi_int_);
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(queued_rows)) {
            assert(row_it != queued_rows.end());
            std::sort(row_tuples.begin(), row_tuples.end(), utils::LessFirst());

            node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);
            auto find_base = node_to_cols_.find(base_node);
            assert(find_base != node_to_cols_.end());
            assert(find_base->second == nannot);

            size_t coord_idx = find_base - node_to_cols_.begin();
            assert(coord_idx < label_coords_.size());

            Columns labels;
            labels.reserve(row_tuples.size());
            auto &label_coords = label_coords_[coord_idx];
            label_coords.reserve(row_tuples.size());
            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                label_coords.emplace_back(coords.begin(), coords.end());
            }

            find_base.value() = cache_column_set(std::move(labels));

            ++row_it;
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(queued_rows)) {
            assert(row_it != queued_rows.end());

            std::sort(labels.begin(), labels.end());

            node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);
            auto find_base = node_to_cols_.find(base_node);
            assert(find_base != node_to_cols_.end());
            assert(find_base->second == nannot);

            find_base.value() = cache_column_set(std::move(labels));

            ++row_it;
        }
    }

    for (const auto &[base_node, node, offset] : to_update) {
        auto find_base = node_to_cols_.find(base_node);
        assert(find_base != node_to_cols_.end());
        assert(find_base->second != nannot);

        size_t coord_idx = find_base - node_to_cols_.begin();
        size_t label_i = find_base->second;

        assert(!node_to_cols_.count(node) || node_to_cols_.find(node)->second != nannot);

        if (node_to_cols_.try_emplace(node, label_i).second && has_coordinates()) {
            assert(coord_idx < label_coords_.size());
            label_coords_.emplace_back(label_coords_[coord_idx]);
            for (auto &coords : label_coords_.back()) {
                for (auto &c : coords) {
                    c += offset;
                }
            }
        }
    }

    assert(!has_coordinates() || node_to_cols_.size() == label_coords_.size());

#ifndef NDEBUG
    for (const auto &[node, val] : node_to_cols_) {
        assert(val != nannot);
    }
#endif
}

auto AnnotationBuffer::get_labels_and_coords(node_index node) const
        -> std::pair<const Columns*, const CoordinateSet*> {
    std::pair<const Columns*, const CoordinateSet*> ret_val { nullptr, nullptr };

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
