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

template <class T1, class T2>
bool overlap_with_diff(const T1 &tuple1, const T2 &tuple2, size_t diff) {
    auto a_begin = tuple1.begin();
    const auto a_end = tuple1.end();
    auto b_begin = tuple2.begin();
    const auto b_end = tuple2.end();

    assert(std::is_sorted(a_begin, a_end));
    assert(std::is_sorted(b_begin, b_end));

    while (a_begin != a_end && b_begin != b_end) {
        if (*a_begin + diff == *b_begin)
            return true;

        if (*a_begin + diff < *b_begin) {
            ++a_begin;
        } else {
            ++b_begin;
        }
    }

    return false;
}

template <class AlignmentCompare>
void LabeledBacktrackingExtender<AlignmentCompare>
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
        DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                             [&](node_index next, char c, score_t score) {
            node_labels_.emplace_back(node_labels_[table_i]);
            callback(next, c, score);
        }, table_i, force_fixed_seed);
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

    // if (outgoing.size() > 1)
        labeled_graph_.flush();

    const Alignment &seed = *this->seed_;

    auto [base_labels, base_coords] = labeled_graph_.get_labels_and_coordinates(node);
    assert(base_labels);
    assert(base_labels->get().size());
    const Vector<Column> &node_labels = labeled_graph_.get_labels_from_index(node_labels_[table_i]);

    if (!base_coords) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        for (const auto &[next, c, score] : outgoing) {
            auto next_labels = labeled_graph_.get_labels(next);
            assert(next_labels);

            Vector<Column> intersect_labels;
            std::set_intersection(node_labels.begin(),
                                  node_labels.end(),
                                  next_labels->get().begin(),
                                  next_labels->get().end(),
                                  std::back_inserter(intersect_labels));

            if (intersect_labels.size()) {
                node_labels_.emplace_back(labeled_graph_.emplace_label_set(
                    std::move(intersect_labels)
                ));
                callback(next, c, score);
            }
        }

        return;
    }

    // check label and coordinate consistency
    assert(seed.label_coordinates.size());

    // first, determine a base node from which to compare coordinates
    // by default, node is used (the parent node of next)
    bool rev_align = dynamic_cast<const RCDBG*>(graph_);

    // if the seed has coordinates, use the seed as the base
    auto seed_label_it = node_labels.begin();
    auto seed_label_end_it = node_labels.end();
    base_labels = std::cref(seed.label_columns);
    base_coords = std::cref(seed.label_coordinates);
    ssize_t offset = std::get<6>(table[table_i]);
    ssize_t dist_from_origin = offset - (seed.get_offset() - 1);
    ssize_t dist = dist_from_origin - seed.get_offset()
        - seed.get_sequence().size()
        + seed.get_nodes().size()
        - (seed.get_sequence().size() - graph_->get_k()) * rev_align;

    for (const auto &[next, c, score] : outgoing) {
        auto [next_labels, next_coords]
            = labeled_graph_.get_labels_and_coordinates(next);

        assert(next_coords);

        // if we are traversing backwards, then negate the coordinate delta
        ssize_t dist_sign = rev_align ? -1 : 1;

        // check if at least one label has consistent coordinates
        Vector<Column> intersect_labels;
        try {
            utils::match_indexed_values(
                base_labels->get().begin(), base_labels->get().end(),
                base_coords->get().begin(),
                next_labels->get().begin(), next_labels->get().end(),
                next_coords->get().begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    while (seed_label_it != seed_label_end_it && c > *seed_label_it) {
                        ++seed_label_it;
                    }

                    if (seed_label_it == seed_label_end_it)
                        throw std::exception();

                    if (overlap_with_diff(coords, other_coords, dist * dist_sign))
                        intersect_labels.push_back(c);
                }
            );
        } catch (const std::exception&) {}

        if (intersect_labels.size()) {
            // found a consistent pair of coordinates
            node_labels_.emplace_back(labeled_graph_.emplace_label_set(
                std::move(intersect_labels)
            ));
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

template <class AIt, class BIt, class OutIt, class OutIt2>
void set_intersection_difference(AIt a_begin,
                                 AIt a_end,
                                 BIt b_begin,
                                 BIt b_end,
                                 OutIt intersection_out,
                                 OutIt2 diff_out) {
    while (a_begin != a_end) {
        if (b_begin == b_end || *a_begin < *b_begin) {
            *diff_out = *a_begin;
            ++diff_out;
            ++a_begin;
        } else if (*a_begin > *b_begin) {
            ++b_begin;
        } else {
            *intersection_out = *a_begin;
            ++intersection_out;
            ++a_begin;
            ++b_begin;
        }
    }
}

template <class AlignmentCompare>
bool LabeledBacktrackingExtender<AlignmentCompare>::skip_backtrack_start(size_t i) {
    assert(fetched_label_i_ != AnnotationBuffer::nannot);
    assert(node_labels_[i] != AnnotationBuffer::nannot);
    if (!this->prev_starts.emplace(i).second)
        return false;

    const Vector<Column> &end_labels = labeled_graph_.get_labels_from_index(node_labels_[i]);
    const Vector<Column> &left_labels = labeled_graph_.get_labels_from_index(fetched_label_i_);

    label_intersection_ = Vector<Column>{};
    Vector<Column> label_diff;
    set_intersection_difference(end_labels.begin(), end_labels.end(),
                                left_labels.begin(), left_labels.end(),
                                std::back_inserter(label_intersection_),
                                std::back_inserter(label_diff));

    fetched_label_i_ = labeled_graph_.emplace_label_set(std::move(label_diff));
    assert(fetched_label_i_ != AnnotationBuffer::nannot);

    return label_intersection_.empty();
}

template <class AlignmentCompare>
void LabeledBacktrackingExtender<AlignmentCompare>
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
                  const std::function<void(Alignment&&)> &callback) {
    DefaultColumnExtender::call_alignments(cur_cell_score, end_score, min_path_score,
                                           path, trace, table_i, ops, clipping, offset,
                                           window, match, extra_penalty,
                                           [&](Alignment&& alignment) {
        alignment.label_encoder = &labeled_graph_.get_annotator().get_label_encoder();

        auto [base_labels, base_coords] = labeled_graph_.get_labels_and_coordinates(alignment.get_nodes().front());
        assert(base_labels);
        assert(base_labels->get().size());

        if (!base_coords) {
            alignment.label_columns = std::move(label_intersection_);
            callback(std::move(alignment));
            return;
        }

        auto seed_label_it = label_intersection_.begin();
        auto seed_label_end_it = label_intersection_.end();

        if (alignment.get_nodes().size() == 1) {
            auto it = base_labels->get().begin();
            auto end = base_labels->get().end();
            auto c_it = base_coords->get().begin();
            while (seed_label_it != seed_label_end_it && it != end) {
                if (*seed_label_it < *it) {
                    ++seed_label_it;
                } else if (*seed_label_it > *it) {
                    ++it;
                    ++c_it;
                } else {
                    alignment.label_columns.emplace_back(*it);
                    alignment.label_coordinates.emplace_back(*c_it);
                    ++it;
                    ++c_it;
                    ++seed_label_it;
                }
            }
        } else {
            auto [cur_labels, cur_coords] = labeled_graph_.get_labels_and_coordinates(alignment.get_nodes().back());
            assert(cur_labels);
            assert(cur_labels->get().size());
            assert(cur_coords);
            assert(cur_coords->get().size());

            bool rev_align = dynamic_cast<const RCDBG*>(graph_);
            ssize_t offset = std::get<6>(table[table_i]);
            const Alignment &seed = *this->seed_;
            ssize_t dist_from_origin = offset - (seed.get_offset() - 1);
            ssize_t dist = dist_from_origin - seed.get_offset()
                - seed.get_sequence().size()
                + seed.get_nodes().size()
                - (seed.get_sequence().size() - graph_->get_k()) * rev_align;
            ssize_t dist_sign = rev_align ? -1 : 1;
            try {
                CoordIntersection intersect_coords(dist * dist_sign);
                utils::match_indexed_values(
                    base_labels->get().begin(), base_labels->get().end(),
                    base_coords->get().begin(),
                    cur_labels->get().begin(), cur_labels->get().end(),
                    cur_coords->get().begin(),
                    [&](Column c, const auto &coords, const auto &other_coords) {
                        while (seed_label_it != seed_label_end_it && c > *seed_label_it) {
                            ++seed_label_it;
                        }

                        if (seed_label_it == seed_label_end_it)
                            throw std::exception();

                        Tuple overlap;
                        intersect_coords(coords.begin(), coords.end(),
                                         other_coords.begin(), other_coords.end(),
                                         std::back_inserter(overlap));
                        if (overlap.size()) {
                            alignment.label_columns.emplace_back(c);
                            alignment.label_coordinates.emplace_back(std::move(overlap));
                        }
                    }
                );
            } catch (const std::exception&) {}
        }

        if (alignment.label_coordinates.empty())
            return;

        if (alignment.get_sequence().size() < this->graph_->get_k()) {
            size_t residual_offset = this->graph_->get_k() - alignment.get_sequence().size();
            for (auto &tuple : alignment.label_coordinates) {
                for (uint64_t &c : tuple) {
                    c += residual_offset;
                }
            }
        }

        callback(std::move(alignment));
    });
}

template class LabeledBacktrackingExtender<>;

} // namespace align
} // namespace graph
} // namespace mtg
