#include "aligner_labeled.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

using MIM = annot::matrix::MultiIntMatrix;


AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const MIM*>(&annotator_.get_matrix())),
        labels_set_({ {} }) {
    if (multi_int_ && graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        multi_int_ = nullptr;
        common::logger->warn("Coordinates not supported when aligning to CANONICAL "
                             "or PRIMARY mode graphs");
    }
}

void AnnotationBuffer::flush() {
    if (added_rows_.empty())
        return;

    auto push_node_labels = [&](auto node_it, auto row_it, auto&& labels) {
        assert(node_it != added_nodes_.end());
        assert(row_it != added_rows_.end());
        auto label_it = labels_set_.emplace(std::forward<decltype(labels)>(labels)).first;
        assert(labels_.count(*node_it));
        assert(labels_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));
        assert(labels_[*node_it].first
            == *(added_rows_.begin() + (node_it - added_nodes_.begin())));

        labels_[*node_it].second = label_it - labels_set_.begin();
        labels_[AnnotatedDBG::anno_to_graph_index(*row_it)].second = label_it - labels_set_.begin();
    };

    auto node_it = added_nodes_.begin();
    auto row_it = added_rows_.begin();
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

            push_node_labels(node_it++, row_it++, std::move(labels));
        }
    } else {
        for (auto&& labels : annotator_.get_matrix().get_rows(added_rows_)) {
            push_node_labels(node_it++, row_it++, std::forward<decltype(labels)>(labels));
        }
    }

    assert(node_it == added_nodes_.end());
    assert(row_it == added_rows_.end());

    added_rows_.clear();
    added_nodes_.clear();
}

struct OverlapWithDiff {
    OverlapWithDiff(ssize_t diff) : diff_(diff) {}

    template <class T1, class T2>
    bool operator()(const T1 &tuple1, const T2 &tuple2) const {
        auto a_begin = tuple1.begin();
        const auto a_end = tuple1.end();
        auto b_begin = tuple2.begin();
        const auto b_end = tuple2.end();

        assert(std::is_sorted(a_begin, a_end));
        assert(std::is_sorted(b_begin, b_end));

        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin + diff_ == *b_begin)
                return true;

            if (*a_begin + diff_ < *b_begin) {
                ++a_begin;
            } else {
                ++b_begin;
            }
        }

        return false;
    }

    ssize_t diff_;
};

void LabeledBacktrackingExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char /* last char */, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    bool in_seed = std::get<6>(table[table_i]) + 1 - this->seed_->get_offset()
                        < this->seed_->get_sequence().size();
    size_t next_offset = std::get<6>(table[table_i]) + 1;
    assert(node == std::get<3>(table[table_i]));

    if (in_seed && (next_offset < graph_->get_k() || force_fixed_seed || this->fixed_seed())) {
        DefaultColumnExtender::call_outgoing(node, max_prefetch_distance, callback,
                                             table_i, force_fixed_seed);
        return;
    }

    assert(graph_->get_node_sequence(node).back() == std::get<5>(table[table_i]));
    for ( ; last_buffered_table_i_ < table.size(); ++last_buffered_table_i_) {
        labeled_graph_.add_node(std::get<3>(table[last_buffered_table_i_]));
    }

    std::vector<std::tuple<node_index, char, score_t>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                         [&](node_index next, char c, score_t score) {
        outgoing.emplace_back(next, c, score);
        labeled_graph_.add_node(next);
    }, table_i, force_fixed_seed);

    if (outgoing.size() > 1)
        labeled_graph_.flush();

    const Alignment &seed = *this->seed_;

    auto [base_labels, base_coords] = labeled_graph_.get_labels_and_coordinates(node);

    if (base_coords) {
        // check label and coordinate consistency

        // first, determine a base node from which to compare coordinates
        // by default, node is used (the parent node of next)
        ssize_t dist = 1;
        bool rev_align = dynamic_cast<const RCDBG*>(graph_);
        if (seed.label_coordinates.size()) {
            // if the seed has coordinates, use the seed as the base
            base_labels = std::cref(seed.label_columns);
            base_coords = std::cref(seed.label_coordinates);
            ssize_t offset = std::get<6>(table[table_i]);
            ssize_t dist_from_origin = offset - (seed.get_offset() - 1);
            dist = dist_from_origin - seed.get_offset()
                - seed.get_sequence().size()
                + seed.get_nodes().size()
                - (seed.get_sequence().size() - graph_->get_k()) * rev_align;
        } else {
            // otherwise, check the first node in the traversal
            ssize_t k = graph_->get_k();
            std::tie(base_labels, base_coords)
                = labeled_graph_.get_labels_and_coordinates(std::get<3>(table[0]));
            ssize_t base_offset = std::get<6>(table[0]);
            if (!base_coords || base_coords->get().empty()) {
                // if the first node is not annotated, then check the last fork point
                size_t base_table_i = xdrop_cutoffs_[std::get<9>(table[table_i])].first;
                std::tie(base_labels, base_coords)
                    = labeled_graph_.get_labels_and_coordinates(std::get<3>(table[base_table_i]));
                base_offset = std::get<6>(table[base_table_i]);
            }

            if (base_coords) {
                // determine the distance of each child node to the selected base node
                dist = std::get<6>(table[table_i]) - base_offset + 1;
                if (base_offset < k)
                    dist -= k - seed.get_offset();
            }

            // if the seed has stored labels, but no stored coordinates, check
            // the base node to make sure that it has intersecting labels with the seed
            if (seed.label_columns.size()
                    && !utils::share_element(base_labels->get().begin(),
                                             base_labels->get().end(),
                                             seed.label_columns.begin(),
                                             seed.label_columns.end())) {
                return;
            }
        }

        for (const auto &[next, c, score] : outgoing) {
            if (base_labels->get().empty()) {
                callback(next, c, score);
                continue;
            }

            auto [next_labels, next_coords]
                = labeled_graph_.get_labels_and_coordinates(next);

            // check coordinate consistency later if they are not cached now
            if (!next_coords) {
                callback(next, c, score);
                continue;
            }

            // if we are traversing backwards, then negate the coordinate delta
            ssize_t dist_sign = rev_align ? -1 : 1;

            // check if at least one label has consistent coordinates
            if (!utils::have_consistent(base_labels->get().begin(), base_labels->get().end(),
                                        base_coords->get().begin(),
                                        next_labels->get().begin(), next_labels->get().end(),
                                        next_coords->get().begin(),
                                        OverlapWithDiff(dist * dist_sign)))
                continue;

            callback(next, c, score);
        }
    } else if (base_labels) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        for (const auto &[next, c, score] : outgoing) {
            if (base_labels->get().empty()) {
                callback(next, c, score);
                continue;
            }

            auto next_labels = labeled_graph_.get_labels(next);

            // if not fetched, check later
            if (!next_labels) {
                callback(next, c, score);
                continue;
            }

            if ((seed.label_columns.size()
                    && !utils::share_element(seed.label_columns.begin(),
                                             seed.label_columns.end(),
                                             next_labels->get().begin(),
                                             next_labels->get().end()))
                    || !utils::share_element(base_labels->get().begin(),
                                             base_labels->get().end(),
                                             next_labels->get().begin(),
                                             next_labels->get().end()))
                continue;

            callback(next, c, score);
        }
    } else {
        for (const auto &[next, c, score] : outgoing) {
            callback(next, c, score);
        }
    }
}

auto AnnotationBuffer::add_path(const std::vector<node_index> &path, std::string&& query)
        -> std::pair<std::vector<node_index>, bool> {
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
                && "PRIMARY graphs must be wrapped into CANONICAL");

    if (path.empty())
        return {};

    assert(query.size() >= graph_.get_k());
    assert(path.size() == query.size() - graph_.get_k() + 1);

    // TODO: this cascade of graph unwrapping is ugly, find a cleaner way to do it
    const DeBruijnGraph *base_graph = &graph_;
    bool reverse_path = false;
    if (const auto *rc_dbg = dynamic_cast<const RCDBG*>(base_graph)) {
        base_graph = &rc_dbg->get_graph();
        reverse_path = true;
    }

    const auto *canonical = dynamic_cast<const CanonicalDBG*>(base_graph);
    if (canonical)
        base_graph = &canonical->get_graph();

    bool base_is_canonical = (base_graph->get_mode() == DeBruijnGraph::CANONICAL);

    std::vector<node_index> base_path;
    if (base_is_canonical) {
        if (query.front() == '#')
            query = graph_.get_node_sequence(path[0]) + query.substr(graph_.get_k());

        if (reverse_path)
            reverse_complement(query.begin(), query.end());

        base_path = map_to_nodes(*base_graph, query);
    } else if (canonical) {
        base_path.reserve(path.size());
        for (node_index node : path) {
            base_path.emplace_back(canonical->get_base_node(node));
        }
    } else {
        base_path = path;
    }

    assert(base_path.size() == path.size());

    if (reverse_path)
        std::reverse(base_path.begin(), base_path.end());

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(base_graph);
    const boss::BOSS *boss = dbg_succ ? &dbg_succ->get_boss() : nullptr;
    for (size_t i = 0; i < path.size(); ++i) {
        if (base_path[i] == DeBruijnGraph::npos)
            continue;

        if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i])))
            continue; // skip dummy nodes

        Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
        if (labels_.emplace(path[i], std::make_pair(row, nannot)).second
                && (path[i] == base_path[i]
                    || labels_.emplace(base_path[i], std::make_pair(row, nannot)).second)) {
            added_rows_.push_back(row);
            added_nodes_.push_back(path[i]);
        } else {
            auto find_n = labels_.find(path[i]);
            auto find_b = labels_.find(base_path[i]);
            assert(find_n != labels_.end());
            assert(find_b != labels_.end());
            assert(find_n->second.first == find_b->second.first);
            auto label_i = std::min(find_n->second.second, find_b->second.second);
            find_n.value().second = label_i;
            find_b.value().second = label_i;
        }
    }

    return std::make_pair(std::move(base_path), reverse_path);
}

auto AnnotationBuffer::add_node(node_index node) -> node_index {
    return add_path({ node }, std::string(graph_.get_k(), '#')).first[0];
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
        std::tie(label_intersection_, label_intersection_coords_) = label_find->second;

    } else {
        // otherwise, take the full label set
        auto [fetch_labels, fetch_coords]
            = labeled_graph_.get_labels_and_coordinates(node);
        if (fetch_labels) {
            label_intersection_ = fetch_labels->get();
            if (fetch_coords) {
                label_intersection_coords_ = fetch_coords->get();
                if (seed_label_coordinates_.size()) {
                    const Alignment &seed = *this->seed_;
                    bool rev_align = dynamic_cast<const RCDBG*>(graph_);

                    ssize_t offset = std::get<6>(this->table[i]);
                    ssize_t dist_from_origin = offset - (seed.get_offset() - 1) - 1;
                    ssize_t dist = dist_from_origin - seed.get_offset()
                                    - seed.get_sequence().size()
                                    + seed.get_nodes().size()
                                    - (seed.get_sequence().size() - graph_->get_k()) * rev_align;
                    ssize_t dist_sign = rev_align ? 1 : -1;

                    Vector<Column> labels_inter;
                    Vector<Tuple> coords_inter;
                    utils::indexed_set_op<Tuple>(
                        label_intersection_.begin(), label_intersection_.end(),
                        label_intersection_coords_.begin(),
                        seed_labels_.begin(), seed_labels_.end(),
                        seed_label_coordinates_.begin(),
                        std::back_inserter(labels_inter), std::back_inserter(coords_inter),
                        CoordIntersection(dist * dist_sign)
                    );
                    std::swap(labels_inter, label_intersection_);
                    std::swap(coords_inter, label_intersection_coords_);
                }
            } else if (seed_labels_.size()) {
                Vector<Column> labels_inter;
                std::set_intersection(seed_labels_.begin(),
                                      seed_labels_.end(),
                                      label_intersection_.begin(),
                                      label_intersection_.end(),
                                      std::back_inserter(labels_inter));
                std::swap(labels_inter, label_intersection_);
            }
        }
    }

    // we already have the labels for the first node in the path
    last_path_size_ = 1;

    assert(!labeled_graph_.get_coordinate_matrix()
        || label_intersection_.size() == label_intersection_coords_.size());

    cur_min_path_score_ = label_intersection_.size()
        ? extensions_.get_min_path_score(label_intersection_)
        : std::numeric_limits<score_t>::max();

    // skip backtracking from this node if no labels could be determined for it
    return label_intersection_.empty();
}

struct CoordDiff {
    CoordDiff(size_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end) {
            if (b_begin == b_end || *a_begin + offset_ < *b_begin) {
                *out = *a_begin;
                ++out;
                ++a_begin;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                ++a_begin;
                ++b_begin;
            }
        }
    }

    size_t offset_;
};

void LabeledBacktrackingExtender
::call_alignments(score_t cur_cell_score,
                  score_t end_score,
                  score_t min_path_score,
                  const std::vector<node_index> &path,
                  const std::vector<size_t> &trace,
                  size_t table_i,
                  const Cigar &ops,
                  size_t clipping,
                  size_t offset,
                  std::string_view window,
                  const std::string &match,
                  score_t extra_penalty,
                  const std::function<void(Alignment&&)> & /* callback */) {
    assert(path.size());
    assert(ops.size());
    assert(label_intersection_.size());

    score_t alignment_score = end_score - (clipping ? cur_cell_score : 0);
    if (cur_min_path_score_ > alignment_score)
        return;

    // Normally, we start from the end of the alignment and reconstruct the
    // alignment backwards (shifting coordinates by -1).
    // If we are aligning backwards, then the coordinates need to be shifted by
    // 1 instead.
    ssize_t coord_step = dynamic_cast<const RCDBG*>(graph_) ? 1 : -1;

    size_t label_path_end = trace.size()
        - std::min(this->graph_->get_k(), trace.size()) + 1;
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
                                        [ncoord](const auto &a) { return a == ncoord; }),
                         coords.end());
        }

        size_t old_size = 0;
        size_t cur_size = 0;
        size_t new_size = 0;

        const Vector<Column> *label_set = nullptr;
        const Vector<Tuple> *label_coord = nullptr;

        auto [fetch_labels, fetch_coords]
                = labeled_graph_.get_labels_and_coordinates(path[last_path_size_]);
        if (fetch_coords) {
            // backtrack if coordinates are consistent
            const Vector<Column> &labels = fetch_labels->get();
            const Vector<Tuple> &coords = fetch_coords->get();
            Vector<Column> label_inter;
            Vector<Tuple> coord_inter;
            label_set = &labels;
            label_coord = &coords;
            std::tie(old_size, cur_size, new_size)
                = utils::indexed_set_op<Tuple>(
                    label_intersection_.begin(), label_intersection_.end(),
                    label_intersection_coords_.begin(),
                    labels.begin(), labels.end(), coords.begin(),
                    std::back_inserter(label_inter), std::back_inserter(coord_inter),
                    CoordIntersection()
                );

            std::swap(label_intersection_, label_inter);
            std::swap(label_intersection_coords_, coord_inter);

            if (label_intersection_.empty()) {
                cur_min_path_score_ = std::numeric_limits<score_t>::max();
                return;
            }

        } else if (fetch_labels) {
            // backtrack if labels are consistent
            const Vector<Column> &labels = fetch_labels->get();
            Vector<Column> inter;
            label_set = &labels;
            std::set_intersection(label_intersection_.begin(), label_intersection_.end(),
                                  labels.begin(), labels.end(), std::back_inserter(inter));
            old_size = label_intersection_.size();
            cur_size = labels.size();
            new_size = inter.size();

            std::swap(label_intersection_, inter);

            if (label_intersection_.empty()) {
                cur_min_path_score_ = std::numeric_limits<score_t>::max();
                return;
            }
        }

        if (cur_size > new_size && this->prev_starts.count(trace[last_path_size_])) {
            // if the new coordinate or label set is a subset, then store
            // the difference to start backtracking from here later
            Vector<Column> diff;
            Vector<Tuple> diff_coords;
            auto prev_find = diff_label_sets_.find(trace[last_path_size_]);
            if (prev_find == diff_label_sets_.end()) {
                if (label_intersection_coords_.size()) {
                    utils::indexed_set_op<Tuple>(
                        label_set->begin(), label_set->end(),
                        label_coord->begin(),
                        label_intersection_.begin(), label_intersection_.end(),
                        label_intersection_coords_.begin(),
                        std::back_inserter(diff), std::back_inserter(diff_coords),
                        CoordDiff()
                    );
                } else if (label_intersection_.size()) {
                    std::set_difference(label_set->begin(), label_set->end(),
                                        label_intersection_.begin(),
                                        label_intersection_.end(),
                                        std::back_inserter(diff));
                }
                if (diff.size()) {
                    diff_label_sets_.emplace(trace[last_path_size_],
                                             std::make_pair(std::move(diff),
                                                            std::move(diff_coords)));
                    this->prev_starts.erase(trace[last_path_size_]);
                }
            } else {
                auto &[prev_labels, prev_coords] = prev_find.value();
                if (label_intersection_coords_.size()) {
                    utils::indexed_set_op<Tuple>(
                        prev_labels.begin(), prev_labels.end(),
                        prev_coords.begin(),
                        label_intersection_.begin(), label_intersection_.end(),
                        label_intersection_coords_.begin(),
                        std::back_inserter(diff), std::back_inserter(diff_coords),
                        CoordDiff()
                    );
                } else {
                    std::set_difference(prev_labels.begin(), prev_labels.end(),
                                        label_intersection_.begin(),
                                        label_intersection_.end(),
                                        std::back_inserter(diff));
                }
                std::swap(prev_labels, diff);
                std::swap(prev_coords, diff_coords);
                if (diff.size())
                    this->prev_starts.erase(trace[last_path_size_]);
            }
        }
    }

    cur_min_path_score_ = label_intersection_.size()
        ? extensions_.get_min_path_score(label_intersection_)
        : std::numeric_limits<score_t>::max();

    if ((!config_.allow_left_trim && table_i) || (
            config_.allow_left_trim && ((clipping && ops.data().back().first != Cigar::MATCH)
            || window.size() < this->config_.min_seed_length
            || alignment_score < std::max(min_path_score, cur_min_path_score_))))
        return;

    Alignment alignment = this->construct_alignment(
        ops, clipping, window, path, match, alignment_score, offset, extra_penalty
    );

    // store the label and coordinate information
    alignment.label_encoder = &labeled_graph_.get_annotator().get_label_encoder();

    alignment.label_columns = label_intersection_;

    if (labeled_graph_.get_coordinate_matrix()) {
        assert(label_intersection_.size() == label_intersection_coords_.size());
        alignment.label_coordinates = label_intersection_coords_;

        // if we were aligning backwards, then c represents the
        // end coordinate, so shift it
        if (coord_step == 1) {
            for (auto &tuple : alignment.label_coordinates) {
                for (uint64_t &c : tuple) {
                    c -= path.size() - 1;
                }
            }
        }

        if (alignment.get_sequence().size() < this->graph_->get_k()) {
            size_t residual_offset = this->graph_->get_k() - alignment.get_sequence().size();
            for (auto &tuple : alignment.label_coordinates) {
                for (uint64_t &c : tuple) {
                    c += residual_offset;
                }
            }
        }
    }

    if (!table_i && seed_labels_.size()) {
        Vector<Column> diff;
        if (seed_label_coordinates_.size()) {
            Vector<Tuple> diff_coords;
            utils::indexed_set_op<Tuple>(
                seed_labels_.begin(), seed_labels_.end(),
                seed_label_coordinates_.begin(),
                alignment.label_columns.begin(), alignment.label_columns.end(),
                alignment.label_coordinates.begin(),
                std::back_inserter(diff), std::back_inserter(diff_coords),
                CoordDiff()
            );
            std::swap(diff_coords, seed_label_coordinates_);
        } else {
            std::set_difference(seed_labels_.begin(), seed_labels_.end(),
                                alignment.label_columns.begin(),
                                alignment.label_columns.end(),
                                std::back_inserter(diff));
        }
        std::swap(diff, seed_labels_);
    }

    // instead of calling back the alignment, we add it to the local alignment buffer
    // after extension is done, the best ones are called back
    extensions_.add_alignment(std::move(alignment));
}

} // namespace align
} // namespace graph
} // namespace mtg
