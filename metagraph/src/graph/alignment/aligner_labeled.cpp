#include "aligner_labeled.hpp"

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

typedef DeBruijnGraph::node_index node_index;


bool check_targets(const DeBruijnGraph &graph,
                   const AnnotatedDBG &anno_graph,
                   const Alignment<node_index> &path) {
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
      : anno_graph_(anno_graph),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(
            &anno_graph_.get_annotation().get_matrix())
        ) {
    targets_set_.emplace(); // insert empty vector
}

void DynamicLabeledGraph::flush() {
    auto node_it = added_nodes_.begin();
    for (const auto &labels : anno_graph_.get_annotation().get_matrix().get_rows(added_rows_)) {
        assert(node_it != added_nodes_.end());

        auto label_it = targets_set_.emplace(labels).first;
        assert(labels == *label_it);
        assert(targets_[*node_it].first
            == *(added_rows_.begin() + (node_it - added_nodes_.begin())));

        targets_[*node_it].second = label_it - targets_set_.begin();

        ++node_it;
    }

    assert(node_it == added_nodes_.end());

    added_rows_.clear();
    added_nodes_.clear();
}

bool DynamicLabeledGraph::is_coord_consistent(node_index node, node_index next,
                                              std::string sequence) const {
    if (!multi_int_ || !node || node == next)
        return true;

    const DeBruijnGraph &graph = anno_graph_.get_graph();
    bool base_is_canonical = graph.get_base_graph().get_mode() == DeBruijnGraph::CANONICAL;

    if (base_is_canonical && sequence.front() == '#') {
        sequence = graph.get_node_sequence(node) + sequence.back();
    }

    auto [base_path, reversed] = graph.get_base_path({ node, next }, sequence);

    if ((!reversed && !base_path[1]) || (reversed && !base_path[0]))
        return false;

    base_path[0] = AnnotatedDBG::graph_to_anno_index(base_path[0]);
    base_path[1] = AnnotatedDBG::graph_to_anno_index(base_path[1]);

    auto tuples = multi_int_->get_row_tuples(base_path);

    std::vector<uint64_t> coordinates[2];
    for (size_t i = 0; i < 2; ++i) {
        for (const auto &[j, tuple] : tuples[i]) {
            for (uint64_t coord : tuple) {
                // TODO: make sure the offsets are correct (query max_int in multi_int)
                // TODO: if this takes up a significant amount of time, preallocate
                //       the entire vector beforehand
                coordinates[i].push_back(j * 1e15 + coord);
            }
        }
        assert(std::is_sorted(coordinates[i].begin(), coordinates[i].end()));
    }

    // check if coord[0] + 1 == coord[1]
    for (uint64_t &c : coordinates[0]) {
        ++c;
    }

    bool check_intersection =
        utils::count_intersection(coordinates[0].begin(), coordinates[0].end(),
                                  coordinates[1].begin(), coordinates[1].end());

    if (check_intersection || !base_is_canonical)
        return check_intersection;

    // for a base canonical graph, also check
    // coord[0] + 1 == coord[1] + 2 <=> coord[0] - 1 == coord[1]
    for (uint64_t &c : coordinates[1]) {
        c += 2;
    }

    return utils::count_intersection(coordinates[0].begin(), coordinates[0].end(),
                                     coordinates[1].begin(), coordinates[1].end());
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
::call_outgoing(NodeType node,
                size_t max_prefetch_distance,
                const std::function<void(NodeType, char /* last char */)> &callback,
                size_t table_idx) {
    auto cached_labels = labeled_graph_[node];
    if (this->config_.label_every_n && cached_labels) {
        max_prefetch_distance = std::min(max_prefetch_distance, this->config_.label_every_n);
        std::vector<NodeType> nodes { node };
        std::string seq(this->graph_->get_k(), '#');
        for (size_t i = 0; i < max_prefetch_distance; ++i) {
            size_t outdegree = 0;
            this->graph_->call_outgoing_kmers(nodes.back(), [&](NodeType next, char c) {
                if (c != boss::BOSS::kSentinel) {
                    if (!outdegree) {
                        nodes.push_back(next);
                        seq += c;
                    }

                    ++outdegree;
                }
            });

            if (outdegree != 1)
                break;
        }

        labeled_graph_.add_path(nodes, seq);
        labeled_graph_.flush();
        cached_labels = labeled_graph_[node];
    }

    std::function<void(NodeType, char)> call = callback;

    if (labeled_graph_.get_coordinate_matrix()) {
        // coordinate consistency
        call = [&,dummy=std::string(this->graph_->get_k(), '#')](NodeType next, char c) mutable {
            dummy.back() = c;
            if (labeled_graph_.is_coord_consistent(node, next, dummy))
                callback(next, c);
        };

    } else if (cached_labels) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        call = [&](NodeType next, char c) {
            auto next_labels = labeled_graph_[next];
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

    DefaultColumnExtender<NodeType>::call_outgoing(node, max_prefetch_distance, call, table_idx);
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

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::update_seed_filter(node_index node,
                                                               size_t query_start,
                                                               const score_t *s_begin,
                                                               const score_t *s_end) {
    if (SeedFilteringExtender<NodeType>::update_seed_filter(node, query_start, s_begin, s_end)) {
        labeled_graph_.add_node(node);
        return true;
    } else {
        return false;
    }
}

template <typename NodeType>
bool LabeledBacktrackingExtender<NodeType>::skip_backtrack_start(size_t i) {
    target_intersection_.clear();

    if (this->prev_starts.emplace(i).second) {
        // if backtracking hasn't been started from here yet, get its labels

        auto target_find = diff_target_sets_.find(i);
        if (target_find != diff_target_sets_.end()) {
            // extract a subset of the labels if this node was previously traversed
            target_intersection_ = target_find->second;

        } else if (auto labels = labeled_graph_[std::get<3>(this->table[i])]) {
            if (this->seed_->target_columns.size()) {
                // if the seed had labels, intersect with those
                std::set_intersection(labels->get().begin(),
                                      labels->get().end(),
                                      this->seed_->target_columns.begin(),
                                      this->seed_->target_columns.end(),
                                      std::back_inserter(target_intersection_));
            } else {
                // otherwise take the full label set
                target_intersection_ = *labels;
            }
        }

        // we already have the labels for the first node in the path
        last_path_size_ = 1;
    }

    // skip backtracking from this node if no labels could be determined for it
    return target_intersection_.empty();
}

template <class AlignmentCompare>
void ILabeledAligner<AlignmentCompare>
::set_target_coordinates(IDBGAligner::DBGAlignment &alignment) const {
    const auto *multi_int = labeled_graph_.get_coordinate_matrix();

    if (!multi_int)
        return;

    using Column = DynamicLabeledGraph::Column;
    const std::vector<IDBGAligner::node_index> &path = alignment.get_nodes();
    const Vector<Column> &target_columns = alignment.target_columns;
    auto &target_coordinates = alignment.target_coordinates;

    typedef std::vector<std::pair<size_t, std::pair<uint64_t, uint64_t>>> TargetCoords;
    target_coordinates.assign(target_columns.size(), TargetCoords{});

    typedef tsl::hopscotch_map<int64_t, std::vector<size_t>> RelativeCoordsMap;
    tsl::hopscotch_map<Column, RelativeCoordsMap> row_coordinates;
    for (Column target : target_columns) {
        row_coordinates.emplace(target, RelativeCoordsMap{});
    }

    auto anno_rows = labeled_graph_.get_anno_rows(path);
    std::vector<size_t> missing_rows;
    for (size_t i = 0; i < anno_rows.size(); ++i) {
        if (anno_rows[i] == DynamicLabeledGraph::nrow)
            missing_rows.push_back(i);
    }

    anno_rows.erase(std::remove_if(anno_rows.begin(), anno_rows.end(),
                                   [](auto row) { return row == DynamicLabeledGraph::nrow; }),
                    anno_rows.end());

    auto tuples = multi_int->get_row_tuples(anno_rows);
    size_t offset = 0;
    for (size_t i = 0; i < tuples.size(); ++i) {
        while (offset < missing_rows.size() && i + offset == missing_rows[offset]) {
            ++offset;
        }

        for (const auto &[j, coords] : tuples[i]) {
            auto find = row_coordinates.find(j);
            if (find != row_coordinates.end()) {
                for (int64_t coord : coords) {
                    find.value()[coord - i - offset].emplace_back(i + offset);
                }
            }
        }
    }

    for (auto it = row_coordinates.begin(); it != row_coordinates.end(); ++it) {
        TargetCoords &cur_target_coords = target_coordinates[it->first];
        bool full_length_found = false;
        for (auto jt = it.value().begin(); jt != it.value().end(); ++jt) {
            int64_t start = jt->first;
            auto &relative_coords = jt.value();
            if (relative_coords.size()) {
                std::sort(relative_coords.begin(), relative_coords.end());
                relative_coords.erase(std::unique(relative_coords.begin(),
                                                  relative_coords.end()),
                                      relative_coords.end());
                std::pair<size_t, size_t> cur_range(relative_coords[0], relative_coords[0]);
                for (size_t i = 1; i < relative_coords.size(); ++i) {
                    if (relative_coords[i] == cur_range.second + 1) {
                        ++cur_range.second;
                    } else {
                        if (cur_range.second - cur_range.first + 1 == path.size()) {
                            assert(!cur_range.first);
                            if (full_length_found)
                                continue;

                            full_length_found = true;
                        }
                        cur_target_coords.emplace_back(cur_range.first,
                            std::make_pair(cur_range.first + start,
                                           cur_range.second + start)
                        );
                        cur_range.first = relative_coords[i];
                        cur_range.second = relative_coords[i];
                    }
                }

                if (cur_range.second - cur_range.first + 1 == path.size()) {
                    assert(!cur_range.first);
                    if (full_length_found)
                        continue;

                    full_length_found = true;
                }

                cur_target_coords.emplace_back(cur_range.first,
                    std::make_pair(cur_range.first + start,
                                   cur_range.second + start)
                );
            }
        }

        if (cur_target_coords.size()) {
            std::sort(cur_target_coords.begin(), cur_target_coords.end(),
                      [](const auto &a, const auto &b) {
                return std::make_tuple(b.second.second - b.second.first, a.first, a.second.first)
                    < std::make_tuple(a.second.second - a.second.first, b.first, b.second.first);
            });

            if (full_length_found) {
                auto it = std::find_if(cur_target_coords.begin(),
                                       cur_target_coords.end(),
                                       [&path](const auto &a) {
                    return a.second.second - a.second.first + 1 < path.size();
                });
                cur_target_coords.erase(it, cur_target_coords.end());
            }
        }
    }
}

template <typename NodeType>
void LabeledBacktrackingExtender<NodeType>
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
                  const std::function<void(DBGAlignment&&)> & /* callback */) {
    assert(path.size());
    assert(ops.size());
    assert(target_intersection_.size());
    assert(trace.size() >= this->graph_->get_k());

    size_t label_path_end = trace.size() - this->graph_->get_k() + 1;
    if (label_path_end > last_path_size_) {
        for (size_t i = last_path_size_; i < label_path_end; ++i) {
            assert(static_cast<size_t>(i) < path.size());
            if (auto labels = labeled_graph_[path[i]]) {
                Vector<Column> inter;
                std::set_intersection(target_intersection_.begin(),
                                      target_intersection_.end(),
                                      labels->get().begin(), labels->get().end(),
                                      std::back_inserter(inter));

                if (labels->get().size() > inter.size() && this->prev_starts.count(trace[i])) {
                    Vector<Column> diff;
                    auto prev_labels = diff_target_sets_.find(trace[i]);
                    if (prev_labels == diff_target_sets_.end()) {
                        std::set_difference(labels->get().begin(), labels->get().end(),
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

                std::swap(target_intersection_, inter);

                if (target_intersection_.empty())
                    return;
            }
        }

        last_path_size_ = label_path_end;
    }

    if (ops.back().first == Cigar::MATCH
            && window.size() >= this->config_.min_seed_length
            && end_score - cur_cell_score >= min_path_score) {
        score_t target_score = std::max(
            aggregator_.get_min_path_score(target_intersection_),
            extensions_.get_min_path_score(target_intersection_)
        );

        if (end_score - cur_cell_score >= target_score) {
            DBGAlignment alignment = this->construct_alignment(
                ops, clipping, window, path, match, end_score - cur_cell_score, offset
            );

            assert(!alignment.get_offset());
            alignment.target_columns = target_intersection_;
            assert(check_targets(*this->graph_, labeled_graph_.get_anno_graph(), alignment));

            extensions_.add_alignment(std::move(alignment));
        }
    }
}

template <typename NodeType>
auto LabeledBacktrackingExtender<NodeType>
::extend(score_t min_path_score, bool fixed_seed) -> std::vector<DBGAlignment> {
    extensions_.clear();
    DefaultColumnExtender<NodeType>::extend(min_path_score, fixed_seed);

    std::vector<DBGAlignment> extensions;
    extensions_.call_alignments(
        [&](DBGAlignment&& alignment) { extensions.emplace_back(std::move(alignment)); },
        [&]() { return extensions.size() && extensions.back().get_score() < min_path_score; }
    );

    return extensions;
}

template class ILabeledAligner<>;
template class LabeledBacktrackingExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
