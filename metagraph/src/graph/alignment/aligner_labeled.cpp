#include "aligner_labeled.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

using MIM = annot::matrix::MultiIntMatrix;


AnnotationBuffer::AnnotationBuffer(const AnnotatedDBG &anno_graph)
      : anno_graph_(anno_graph),
        multi_int_(dynamic_cast<const MIM*>(&anno_graph_.get_annotation().get_matrix())) {
    if (multi_int_ && anno_graph.get_graph().get_mode() == DeBruijnGraph::CANONICAL) {
        multi_int_ = nullptr;
        common::logger->warn("Coordinates not supported when aligning to CANONICAL "
                             "or PRIMARY mode graphs");
    }

    labels_set_.emplace(); // insert empty vector
}

void AnnotationBuffer::flush() {
    auto push_node_labels = [&](auto node_it, auto&& labels) {
        assert(node_it != added_nodes_.end());
        auto label_it = labels_set_.emplace(std::forward<decltype(labels)>(labels)).first;
        assert(labels_[*node_it].first
            == *(added_rows_.begin() + (node_it - added_nodes_.begin())));

        labels_[*node_it].second = label_it - labels_set_.begin();
    };

    auto node_it = added_nodes_.begin();
    if (multi_int_) {
        // extract both labels and coordinates, then store them separately
        for (auto&& row_tuples : multi_int_->get_row_tuples(added_rows_)) {
            Vector<Column> labels;
            labels.reserve(row_tuples.size());
            label_coords_.emplace_back();
            label_coords_.back().reserve(row_tuples.size());
            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                label_coords_.back().emplace_back(std::forward<decltype(coords)>(coords));
            }

            push_node_labels(node_it++, std::move(labels));
        }
    } else {
        for (auto&& labels : anno_graph_.get_annotation().get_matrix().get_rows(added_rows_)) {
            push_node_labels(node_it++, std::forward<decltype(labels)>(labels));
        }
    }

    assert(node_it == added_nodes_.end());

    added_rows_.clear();
    added_nodes_.clear();
}

void LabeledBacktrackingExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char /* last char */)> &callback) {
    std::function<void(node_index, char)> call = callback;

    // TODO: checking for coordinate consistency here can be slow

    if (auto cached_labels = labeled_graph_.get_labels(node)) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        call = [&](node_index next, char c) {
            auto next_labels = labeled_graph_.get_labels(next);
            // If labels at the next node are not cached, always take the edge.
            // In this case, the label consistency will be checked later.
            // If they are cached, the existence of at least one common label is checked.
            if (!next_labels || utils::share_element(cached_labels->get().begin(),
                                                     cached_labels->get().end(),
                                                     next_labels->get().begin(),
                                                     next_labels->get().end())) {
                callback(next, c);
            }
        };
    }

    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance, call);
}

void AnnotationBuffer::add_path(const std::vector<node_index> &path, std::string query) {
    assert(anno_graph_.get_graph().get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    if (path.empty())
        return;

    const DeBruijnGraph &graph = anno_graph_.get_graph();
    const DeBruijnGraph &base_graph = graph.get_base_graph();
    bool base_is_canonical = base_graph.get_mode() == DeBruijnGraph::CANONICAL;

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;

    if (base_is_canonical && query.front() == '#')
        query = graph.get_node_sequence(path[0]) + query.substr(graph.get_k());

    auto call_node = [&](node_index node, node_index base_node) {
        if (base_node != DeBruijnGraph::npos) {
            if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_node)))
                return; // skip dummy nodes

            Row row = AnnotatedDBG::graph_to_anno_index(base_node);
            if (labels_.emplace(node, std::make_pair(row, nannot)).second) {
                added_rows_.push_back(row);
                added_nodes_.push_back(node);
            }
        }
    };
    auto [base_path, reversed] = graph.get_base_path(path, query);
    if (!reversed) {
        for (size_t i = 0; i < base_path.size(); ++i) {
            call_node(path[i], base_path[i]);
        }
    } else {
        auto it = path.rbegin();
        for (size_t i = 0; i < base_path.size(); ++i, ++it) {
            call_node(*it, base_path[i]);
        }
    }
}

void AnnotationBuffer::add_node(node_index node) {
    add_path({ node }, std::string(anno_graph_.get_graph().get_k(), '#'));
}

std::vector<Alignment> LabeledBacktrackingExtender
::extend(score_t min_path_score, bool force_fixed_seed) {
    // the overridden backtrack populates extensions_, so this should return nothing
    DefaultColumnExtender::extend(min_path_score, force_fixed_seed);

    // fetch the alignments from extensions_
    return extensions_.get_alignments();
}

bool LabeledBacktrackingExtender::skip_backtrack_start(size_t i) {
    label_intersection_.clear();
    label_intersection_coords_.clear();

    if (!this->prev_starts.emplace(i).second)
        return false;

    // if backtracking hasn't been started from here yet, get its labels
    node_index node = std::get<3>(this->table[i]);
    auto label_find = diff_label_sets_.find(i);
    if (label_find != diff_label_sets_.end()) {
        // extract a subset of the labels if this node was previously traversed
        label_intersection_ = label_find->second;

        // extract coordinates for these labels
        if (auto fetch = labeled_graph_.get_coordinates(node)) {
            const Vector<Column> &labels = fetch->first.get();
            const Vector<Tuple> &coords = fetch->second.get();
            assert(std::includes(labels.begin(), labels.end(),
                                 label_intersection_.begin(),
                                 label_intersection_.end()));
            auto it = labels.begin();
            for (Column label : label_intersection_) {
                it = std::lower_bound(it, labels.end(), label);
                label_intersection_coords_.push_back(coords[it - labels.begin()]);
            }
        }

    } else if (auto labels = labeled_graph_.get_labels(node)) {
        if (this->seed_->label_columns.size()) {
            // if the seed had labels, intersect with those
            std::set_intersection(labels->get().begin(),
                                  labels->get().end(),
                                  this->seed_->label_columns.begin(),
                                  this->seed_->label_columns.end(),
                                  std::back_inserter(label_intersection_));
            if (auto fetch = labeled_graph_.get_coordinates(node)) {
                const Vector<Column> &coord_labels = fetch->first.get();
                const Vector<Tuple> &coords = fetch->second.get();
                assert(std::includes(coord_labels.begin(), coord_labels.end(),
                                     label_intersection_.begin(),
                                     label_intersection_.end()));
                auto it = coord_labels.begin();
                for (Column label : label_intersection_) {
                    it = std::lower_bound(it, coord_labels.end(), label);
                    label_intersection_coords_.push_back(coords[it - coord_labels.begin()]);
                }
            }
        } else {
            // otherwise take the full label set
            label_intersection_ = *labels;
            if (auto fetch = labeled_graph_.get_coordinates(node)) {
                assert(fetch->first.get() == label_intersection_);
                label_intersection_coords_ = fetch->second.get();
            }
        }
    }

    // we already have the labels for the first node in the path
    last_path_size_ = 1;

    assert(!labeled_graph_.get_coordinate_matrix()
        || label_intersection_.size() == label_intersection_coords_.size());

    // skip backtracking from this node if no labels could be determined for it
    return label_intersection_.empty();
}

struct SetIntersection {
    template <typename... Args>
    void operator()(Args&&... args) const {
        std::set_intersection(std::forward<Args>(args)...);
    }
};

void LabeledBacktrackingExtender
::call_alignments(score_t cur_cell_score,
                  score_t end_score,
                  score_t min_path_score,
                  const std::vector<node_index> &path,
                  const std::vector<size_t> &trace,
                  const Cigar &ops,
                  size_t clipping,
                  size_t offset,
                  std::string_view window,
                  const std::string &match,
                  const std::function<void(Alignment&&)> & /* callback */) {
    assert(path.size());
    assert(ops.size());
    assert(label_intersection_.size());
    assert(trace.size() >= this->graph_->get_k());

    // Normally, we start from the end of the alignment and reconstruct the
    // alignment backwards (shifting coordinates by -1).
    // If we are aligning backwards, then the coordinates need to be shifted by
    // 1 instead.
    ssize_t coord_step = dynamic_cast<const RCDBG*>(graph_) ? 1 : -1;

    size_t label_path_end = trace.size() - this->graph_->get_k() + 1;
    assert(label_path_end <= path.size());
    assert(label_path_end >= last_path_size_);
    constexpr uint64_t ncoord = std::numeric_limits<uint64_t>::max();
    for ( ; last_path_size_ < label_path_end; ++last_path_size_) {
        // update current coordinates
        for (auto &coords : label_intersection_coords_) {
            for (uint64_t &c : coords) {
                if (!c && coord_step == -1) {
                    c = ncoord;
                } else {
                    c += coord_step;
                }
            }

            // if any coordinates went out of bounds, remove them
            coords.erase(std::remove_if(coords.begin(), coords.end(),
                                        [&](const auto &a) { return a == ncoord; }),
                         coords.end());
        }

        const Vector<Column> *label_set = nullptr;

        if (auto label_coords = labeled_graph_.get_coordinates(path[last_path_size_])) {
            // backtrack if coordinates are consistent
            const Vector<Column> &labels = label_coords->first.get();
            const Vector<Tuple> &coords = label_coords->second.get();
            label_set = &labels;
            Vector<Column> label_inter;
            Vector<Tuple> coord_inter;
            utils::indexed_set_op<Tuple, SetIntersection>(
                label_intersection_.begin(),
                label_intersection_.end(),
                label_intersection_coords_.begin(),
                labels.begin(), labels.end(), coords.begin(),
                std::back_inserter(label_inter),
                std::back_inserter(coord_inter)
            );

            std::swap(label_intersection_, label_inter);
            std::swap(label_intersection_coords_, coord_inter);

            if (label_intersection_.empty())
                return;

        } else if (auto labels = labeled_graph_.get_labels(path[last_path_size_])) {
            // backtrack if labels are consistent
            label_set = &labels->get();
            Vector<Column> inter;
            std::set_intersection(label_intersection_.begin(),
                                  label_intersection_.end(),
                                  labels->get().begin(), labels->get().end(),
                                  std::back_inserter(inter));

            std::swap(label_intersection_, inter);

            if (label_intersection_.empty())
                return;
        }

        if (label_set && label_set->size() > label_intersection_.size()
                && this->prev_starts.count(trace[last_path_size_])) {
            // if the new coordinate or label set is a subset, then store
            // the difference to start backtracking from here later
            Vector<Column> diff;
            auto prev_labels = diff_label_sets_.find(trace[last_path_size_]);
            if (prev_labels == diff_label_sets_.end()) {
                std::set_difference(label_set->begin(), label_set->end(),
                                    label_intersection_.begin(),
                                    label_intersection_.end(),
                                    std::back_inserter(diff));
                if (diff.size()) {
                    diff_label_sets_.emplace(trace[last_path_size_], std::move(diff));
                    this->prev_starts.erase(trace[last_path_size_]);
                }
            } else {
                std::set_difference(prev_labels->second.begin(),
                                    prev_labels->second.end(),
                                    label_intersection_.begin(),
                                    label_intersection_.end(),
                                    std::back_inserter(diff));
                std::swap(prev_labels.value(), diff);
                if (diff.size())
                    this->prev_starts.erase(trace[last_path_size_]);
            }
        }
    }

    score_t alignment_score = end_score - cur_cell_score;

    if (ops.data().back().first != Cigar::MATCH
            || window.size() < this->config_.min_seed_length
            || alignment_score < min_path_score)
        return;

    // if the current partial alignment is valid, then add it to the
    // alignment aggregator
    score_t label_score = std::max(
        aggregator_.get_min_path_score(label_intersection_),
        extensions_.get_min_path_score(label_intersection_)
    );

    if (alignment_score < label_score)
        return;

    Alignment alignment = this->construct_alignment(
        ops, clipping, window, path, match, alignment_score, offset
    );

    // alignments with labels should at least cover a full node
    assert(!alignment.get_offset());

    // store the label and coordinate information
    alignment.label_encoder
        = &labeled_graph_.get_anno_graph().get_annotation().get_label_encoder();
    alignment.label_columns = label_intersection_;

    if (labeled_graph_.get_coordinate_matrix()) {
        assert(label_intersection_.size() == label_intersection_coords_.size());
        alignment.label_coordinates = label_intersection_coords_;

        if (coord_step == 1) {
            // if we were aligning backwards, then c represents the
            // end coordinate, so shift it
            for (auto &tuple : alignment.label_coordinates) {
                for (uint64_t &c : tuple) {
                    c -= path.size() - 1;
                }
            }
        }
    }

    extensions_.add_alignment(std::move(alignment));
}

template class ILabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
