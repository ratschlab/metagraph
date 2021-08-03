#include "aligner_labeled.hpp"

#include <unordered_set>

#include <tsl/hopscotch_set.h>

#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/utils/template_utils.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

using MIM = annot::matrix::MultiIntMatrix;


bool check_labels(const DeBruijnGraph &graph,
                  const AnnotatedDBG &anno_graph,
                  const Alignment &path) {
    std::string query = path.get_sequence();
    if (dynamic_cast<const RCDBG*>(&graph))
        ::reverse_complement(query.begin(), query.end());

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<annot::binmat::BinaryMatrix::Column> ref_labels;
    for (const std::string &label : anno_graph.get_labels(query, 1.0)) {
        ref_labels.emplace(label_encoder.encode(label));
    }

    for (annot::binmat::BinaryMatrix::Column label : path.label_columns) {
        if (!ref_labels.count(label))
            return false;
    }

    return true;
}

// Return true if the two sorted ranges share a common element
// adapted from:
// https://www.fluentcpp.com/2020/07/03/how-to-check-if-2-sorted-collections-have-a-common-element/
template <class InputIt1, class InputIt2>
constexpr bool share_element(InputIt1 first1,
                             InputIt1 last1,
                             InputIt2 first2,
                             InputIt2 last2) {
    assert(std::is_sorted(first1, last1));
    assert(std::is_sorted(first2, last2));

    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            first1 = std::lower_bound(first1, last1, *first2);
        } else if (*first2 < *first1) {
            first2 = std::lower_bound(first2, last2, *first1);
        } else {
            return true;
        }
    }

    return false;
}

// Intersect the sorted ranges a1 and b1 with corresponding sorted ranges of
// sorted ranges a2 and b2 (of equal length).
// i.e., For each shared element between a1 and b1, intersect the corresponding
// ranges in a2 and b2.
template <class OutType, class InputIt1, class InputIt2, class InputIt3, class InputIt4,
          class OutputIt1, class OutputIt2>
void indexed_set_intersect(InputIt1 a1_begin,
                           InputIt1 a1_end,
                           InputIt2 a2_begin,
                           InputIt3 b1_begin,
                           InputIt3 b1_end,
                           InputIt4 b2_begin,
                           OutputIt1 out1,
                           OutputIt2 out2) {
    while (a1_begin != a1_end && b1_begin != b1_end) {
        if (*a1_begin < *b1_begin) {
            auto a1_begin_next = std::lower_bound(a1_begin, a1_end, *b1_begin);
            a2_begin += a1_begin_next - a1_begin;
            a1_begin = a1_begin_next;
        } else if (*b1_begin < *a1_begin) {
            auto b1_begin_next = std::lower_bound(b1_begin, b1_end, *a1_begin);
            b2_begin += b1_begin_next - b1_begin;
            b1_begin = b1_begin_next;
        } else {
            OutType merged;
            std::set_intersection(a2_begin->begin(), a2_begin->end(),
                                  b2_begin->begin(), b2_begin->end(),
                                  std::back_inserter(merged));
            if (merged.size()) {
                *out1 = *a1_begin;
                ++out1;
                *out2 = std::move(merged);
                ++out2;
            }
            ++a1_begin;
            ++b1_begin;
            ++a2_begin;
            ++b2_begin;
        }
    }
}

DynamicLabeledGraph::DynamicLabeledGraph(const AnnotatedDBG &anno_graph)
      : anno_graph_(anno_graph),
        multi_int_(dynamic_cast<const MIM*>(&anno_graph_.get_annotation().get_matrix())) {
    if (multi_int_ && anno_graph.get_graph().get_mode() == DeBruijnGraph::CANONICAL) {
        multi_int_ = nullptr;
        common::logger->warn("Coordinates not supported when aligning to CANONICAL "
                             "or PRIMARY mode graphs");
    }

    labels_set_.emplace(); // insert empty vector
}

void DynamicLabeledGraph::flush() {
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
                const std::function<void(node_index, char /* last char */)> &callback,
                size_t table_idx) {
    std::function<void(node_index, char)> call = callback;

    if (auto cached_labels = labeled_graph_.get_labels(node)) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        call = [&](node_index next, char c) {
            auto next_labels = labeled_graph_.get_labels(next);
            // If labels at the next node are not cached, always take the edge.
            // In this case, the label consistency will be checked later.
            // If they are cached, the existence of at least one common label is checked.
            if (!next_labels || share_element(cached_labels->get().begin(),
                                              cached_labels->get().end(),
                                              next_labels->get().begin(),
                                              next_labels->get().end())) {
                callback(next, c);
            }
        };
    }

    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance, call, table_idx);
}

void DynamicLabeledGraph::add_path(const std::vector<node_index> &path,
                                   std::string query) {
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

void DynamicLabeledGraph::add_node(node_index node) {
    add_path({ node }, std::string(anno_graph_.get_graph().get_k(), '#'));
}

bool LabeledBacktrackingExtender::skip_backtrack_start(size_t i) {
    label_intersection_.clear();
    label_intersection_coords_.clear();

    if (this->prev_starts.emplace(i).second) {
        // if backtracking hasn't been started from here yet, get its labels
        node_index node = std::get<3>(this->table[i]);
        auto label_find = diff_label_sets_.find(i);
        if (label_find != diff_label_sets_.end()) {
            // extract a subset of the labels if this node was previously traversed
            label_intersection_ = label_find->second;
            if (auto fetch = labeled_graph_.get_coordinates(node)) {
                const Vector<Column> &labels = fetch->first.get();
                const Vector<Tuple> &coords = fetch->second.get();
                assert(std::includes(labels.begin(), labels.end(),
                                     label_intersection_.begin(), label_intersection_.end()));
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
                                         label_intersection_.begin(), label_intersection_.end()));
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
    }

    assert(!labeled_graph_.get_coordinate_matrix()
        || label_intersection_.size() == label_intersection_coords_.size());

    // skip backtracking from this node if no labels could be determined for it
    return label_intersection_.empty();
}

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

    ssize_t coord_step = dynamic_cast<const RCDBG*>(graph_) ? 1 : -1;

    size_t label_path_end = trace.size() - this->graph_->get_k() + 1;
    if (label_path_end > last_path_size_) {
        assert(label_path_end <= path.size());
        for (size_t i = last_path_size_; i < label_path_end; ++i) {
            for (auto &coords : label_intersection_coords_) {
                for (uint64_t &c : coords) {
                    c += coord_step;
                }
            }

            const Vector<Column> *label_set = nullptr;

            if (auto label_coords = labeled_graph_.get_coordinates(path[i])) {
                const Vector<Column> &labels = label_coords->first.get();
                const Vector<Tuple> &coords = label_coords->second.get();
                label_set = &labels;
                Vector<Column> label_inter;
                Vector<Tuple> coord_inter;
                indexed_set_intersect<Tuple>(
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

            } else if (auto labels = labeled_graph_.get_labels(path[i])) {
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
                    && this->prev_starts.count(trace[i])) {
                Vector<Column> diff;
                auto prev_labels = diff_label_sets_.find(trace[i]);
                if (prev_labels == diff_label_sets_.end()) {
                    std::set_difference(label_set->begin(), label_set->end(),
                                        label_intersection_.begin(),
                                        label_intersection_.end(),
                                        std::back_inserter(diff));
                    if (diff.size()) {
                        diff_label_sets_.emplace(trace[i], std::move(diff));
                        this->prev_starts.erase(trace[i]);
                    }
                } else {
                    std::set_difference(prev_labels->second.begin(),
                                        prev_labels->second.end(),
                                        label_intersection_.begin(),
                                        label_intersection_.end(),
                                        std::back_inserter(diff));
                    std::swap(prev_labels.value(), diff);
                    if (diff.size())
                        this->prev_starts.erase(trace[i]);
                }
            }
        }

        last_path_size_ = label_path_end;
    }

    if (ops.data().back().first == Cigar::MATCH
            && window.size() >= this->config_.min_seed_length
            && end_score - cur_cell_score >= min_path_score) {
        score_t label_score = std::max(
            aggregator_.get_min_path_score(label_intersection_),
            extensions_.get_min_path_score(label_intersection_)
        );

        if (end_score - cur_cell_score >= label_score) {
            Alignment alignment = this->construct_alignment(
                ops, clipping, window, path, match, end_score - cur_cell_score, offset
            );

            assert(!alignment.get_offset());
            alignment.label_columns = label_intersection_;
            assert(check_labels(*this->graph_, labeled_graph_.get_anno_graph(), alignment));

            if (labeled_graph_.get_coordinate_matrix()) {
                assert(label_intersection_.size() == label_intersection_coords_.size());
                alignment.label_coordinates.resize(label_intersection_coords_.size());
                for (size_t i = 0; i < label_intersection_coords_.size(); ++i) {
                    auto &cur_coords = alignment.label_coordinates[i];
                    for (uint64_t c : label_intersection_coords_[i]) {
                        // alignment coordinates are 1-based
                        ++c;
                        if (coord_step == 1)
                            c = c + graph_->get_k() - alignment.get_sequence().size();

                        cur_coords.emplace_back(c, c + alignment.get_sequence().size() - 1);
                    }
                }
            }

            extensions_.add_alignment(std::move(alignment));
        }
    }
}

template class ILabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
