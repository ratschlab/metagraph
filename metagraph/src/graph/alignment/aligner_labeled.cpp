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


bool check_targets(const DeBruijnGraph &graph,
                   const AnnotatedDBG &anno_graph,
                   const Alignment &path) {
    std::string query = path.get_sequence();
    if (dynamic_cast<const RCDBG*>(&graph))
        ::reverse_complement(query.begin(), query.end());

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    tsl::hopscotch_set<annot::binmat::BinaryMatrix::Column> ref_targets;
    for (const std::string &label : anno_graph.get_labels(query, 1.0)) {
        ref_targets.emplace(label_encoder.encode(label));
    }

    for (annot::binmat::BinaryMatrix::Column target : path.target_columns) {
        if (!ref_targets.count(target))
            return false;
    }

    return true;
}

DynamicLabeledGraph::DynamicLabeledGraph(const AnnotatedDBG &anno_graph)
      : anno_graph_(anno_graph), multi_int_(nullptr) {
    if (anno_graph.get_graph().get_mode() != DeBruijnGraph::CANONICAL) {
        multi_int_ = dynamic_cast<const annot::matrix::MultiIntMatrix*>(
            &anno_graph_.get_annotation().get_matrix()
        );
    }

    targets_set_.emplace(); // insert empty vector
}

void DynamicLabeledGraph::flush() {
    auto push_node_labels = [&](auto node_it, auto&& labels) {
        assert(node_it != added_nodes_.end());
        auto label_it = targets_set_.emplace(std::forward<decltype(labels)>(labels)).first;
        assert(targets_[*node_it].first
            == *(added_rows_.begin() + (node_it - added_nodes_.begin())));

        targets_[*node_it].second = label_it - targets_set_.begin();
    };

    auto node_it = added_nodes_.begin();
    if (multi_int_) {
        for (auto&& row_tuples : multi_int_->get_row_tuples(added_rows_)) {
            Vector<Column> labels;
            labels.reserve(row_tuples.size());
            target_coords_.emplace_back();
            target_coords_.back().reserve(row_tuples.size());
            for (auto&& [label, coords] : row_tuples) {
                labels.push_back(label);
                target_coords_.back().emplace_back(std::forward<decltype(coords)>(coords));
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
            if (!next_labels || utils::count_intersection(cached_labels->get().begin(),
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
            if (targets_.emplace(node, std::make_pair(row, nannot)).second) {
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
    target_intersection_.clear();
    target_intersection_coords_.clear();

    if (this->prev_starts.emplace(i).second) {
        // if backtracking hasn't been started from here yet, get its labels
        node_index node = std::get<3>(this->table[i]);
        auto target_find = diff_target_sets_.find(i);
        if (target_find != diff_target_sets_.end()) {
            // extract a subset of the labels if this node was previously traversed
            target_intersection_ = target_find->second;
            if (auto fetch = labeled_graph_.get_coordinates(node)) {
                const Vector<Column> &labels = fetch->first.get();
                const Vector<Tuple> &coords = fetch->second.get();
                assert(std::includes(labels.begin(), labels.end(),
                                     target_intersection_.begin(), target_intersection_.end()));
                auto it = labels.begin();
                for (Column label : target_intersection_) {
                    it = std::lower_bound(it, labels.end(), label);
                    target_intersection_coords_.push_back(coords[it - labels.begin()]);
                }
            }

        } else if (auto labels = labeled_graph_.get_labels(node)) {
            if (this->seed_->target_columns.size()) {
                // if the seed had labels, intersect with those
                std::set_intersection(labels->get().begin(),
                                      labels->get().end(),
                                      this->seed_->target_columns.begin(),
                                      this->seed_->target_columns.end(),
                                      std::back_inserter(target_intersection_));
                if (auto fetch = labeled_graph_.get_coordinates(node)) {
                    const Vector<Column> &coord_labels = fetch->first.get();
                    const Vector<Tuple> &coords = fetch->second.get();
                    assert(std::includes(coord_labels.begin(), coord_labels.end(),
                                         target_intersection_.begin(), target_intersection_.end()));
                    auto it = coord_labels.begin();
                    for (Column label : target_intersection_) {
                        it = std::lower_bound(it, coord_labels.end(), label);
                        target_intersection_coords_.push_back(coords[it - coord_labels.begin()]);
                    }
                }
            } else {
                // otherwise take the full label set
                target_intersection_ = *labels;
                if (auto fetch = labeled_graph_.get_coordinates(node)) {
                    assert(fetch->first.get() == *labels);
                    target_intersection_coords_ = fetch->second.get();
                }
            }
        }

        // we already have the labels for the first node in the path
        last_path_size_ = 1;
    }

    assert(!labeled_graph_.get_coordinate_matrix()
        || target_intersection_.size() == target_intersection_coords_.size());

    // skip backtracking from this node if no labels could be determined for it
    return target_intersection_.empty();
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
    assert(target_intersection_.size());
    assert(trace.size() >= this->graph_->get_k());

    size_t label_path_end = trace.size() - this->graph_->get_k() + 1;
    if (label_path_end > last_path_size_) {
        for (size_t i = last_path_size_; i < label_path_end; ++i) {
            assert(static_cast<size_t>(i) < path.size());

            for (auto &coords : target_intersection_coords_) {
                for (uint64_t &c : coords) {
                    --c;
                }
            }

            const Vector<Column> *label_set = nullptr;

            if (auto label_coords = labeled_graph_.get_coordinates(path[i])) {
                const Vector<Column> &labels = label_coords->first.get();
                const Vector<Tuple> &coords = label_coords->second.get();
                label_set = &labels;
                Vector<Column> label_inter;
                Vector<Tuple> coord_inter;
                auto it = target_intersection_.begin();
                auto jt = labels.begin();
                auto kt = target_intersection_coords_.begin();
                auto lt = coords.begin();
                while (it != target_intersection_.end() && jt != labels.end()) {
                    if (*it < *jt) {
                        auto it_next = std::lower_bound(it, target_intersection_.end(), *jt);
                        kt += it_next - it;
                        it = it_next;
                    } else if (*jt < *it) {
                        auto jt_next = std::lower_bound(jt, labels.end(), *it);
                        lt += jt_next - jt;
                        jt = jt_next;
                    } else {
                        Tuple coord_merge;
                        std::set_intersection(lt->begin(), lt->end(), kt->begin(), kt->end(),
                                              std::back_inserter(coord_merge));
                        if (coord_merge.size()) {
                            label_inter.push_back(*it);
                            coord_inter.push_back(std::move(coord_merge));
                        }
                        ++it;
                        ++jt;
                        ++kt;
                        ++lt;
                    }
                }

                std::swap(target_intersection_, label_inter);
                std::swap(target_intersection_coords_, coord_inter);

                if (target_intersection_.empty())
                    return;

            } else if (auto labels = labeled_graph_.get_labels(path[i])) {
                label_set = &labels->get();
                Vector<Column> inter;
                std::set_intersection(target_intersection_.begin(),
                                      target_intersection_.end(),
                                      labels->get().begin(), labels->get().end(),
                                      std::back_inserter(inter));

                std::swap(target_intersection_, inter);

                if (target_intersection_.empty())
                    return;
            }

            if (label_set && label_set->size() > target_intersection_.size()
                    && this->prev_starts.count(trace[i])) {
                Vector<Column> diff;
                auto prev_labels = diff_target_sets_.find(trace[i]);
                if (prev_labels == diff_target_sets_.end()) {
                    std::set_difference(label_set->begin(), label_set->end(),
                                        target_intersection_.begin(),
                                        target_intersection_.end(),
                                        std::back_inserter(diff));
                    if (diff.size()) {
                        diff_target_sets_.emplace(trace[i], std::move(diff));
                        this->prev_starts.erase(trace[i]);
                    }
                } else {
                    std::set_difference(prev_labels->second.begin(),
                                        prev_labels->second.end(),
                                        target_intersection_.begin(),
                                        target_intersection_.end(),
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
        score_t target_score = std::max(
            aggregator_.get_min_path_score(target_intersection_),
            extensions_.get_min_path_score(target_intersection_)
        );

        if (end_score - cur_cell_score >= target_score) {
            Alignment alignment = this->construct_alignment(
                ops, clipping, window, path, match, end_score - cur_cell_score, offset
            );

            assert(!alignment.get_offset());
            alignment.target_columns = target_intersection_;
            assert(check_targets(*this->graph_, labeled_graph_.get_anno_graph(), alignment));

            if (labeled_graph_.get_coordinate_matrix()) {
                assert(target_intersection_.size() == target_intersection_coords_.size());
                alignment.target_coordinates.resize(target_intersection_coords_.size());
                for (size_t i = 0; i < target_intersection_coords_.size(); ++i) {
                    auto &cur_coords = alignment.target_coordinates[i];
                    for (uint64_t c : target_intersection_coords_[i]) {
                        // alignment coordinates are 1-based
                        cur_coords.emplace_back(
                            c + 1, c + alignment.get_nodes().size() + graph_->get_k() - 1
                        );
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
