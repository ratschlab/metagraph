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

AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&annotator_.get_matrix())),
        labels_set_({ {} }) {
    if (multi_int_ && graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        multi_int_ = nullptr;
        common::logger->warn("Coordinates not supported when aligning to CANONICAL "
                             "or PRIMARY mode graphs");
    }
}

auto AnnotationBuffer::add_node(node_index node) -> node_index {
    return add_path({ node }, std::string(graph_.get_k(), '#')).first[0];
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
        auto find_a = labels_.find(path[i]);
        auto find_b = labels_.find(base_path[i]);

        if (base_path[i] == DeBruijnGraph::npos) {
            // this can happen when the base graph is CANONICAL and path[i] is a
            // dummy node
            if (find_a != labels_.end()) {
                find_a.value() = std::make_pair(AnnotatedDBG::graph_to_anno_index(base_path[i]), 0);
            } else {
                labels_[path[i]] = std::make_pair(AnnotatedDBG::graph_to_anno_index(base_path[i]), 0);
                if (multi_int_)
                    label_coords_.emplace_back();
            }
            continue;
        }

        if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i]))) {
            // skip dummy nodes
            if (find_a != labels_.end()) {
                find_a.value() = std::make_pair(AnnotatedDBG::graph_to_anno_index(base_path[i]), 0);
            } else {
                labels_[path[i]] = std::make_pair(AnnotatedDBG::graph_to_anno_index(base_path[i]), 0);
                if (multi_int_)
                    label_coords_.emplace_back();
                find_b = labels_.find(base_path[i]);
            }

            if (find_b != labels_.end()) {
                find_b.value() = std::make_pair(AnnotatedDBG::graph_to_anno_index(base_path[i]), 0);
            } else {
                labels_[base_path[i]] = std::make_pair(AnnotatedDBG::graph_to_anno_index(base_path[i]), 0);
                if (multi_int_)
                    label_coords_.emplace_back();
            }
            continue;
        }

        Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);

        if (find_a == labels_.end() && find_b == labels_.end()) {
            auto val = std::make_pair(row, nannot);
            labels_.emplace(path[i], val);
            added_rows_.push_back(row);
            added_nodes_.push_back(path[i]);

            if (path[i] != base_path[i]) {
                labels_.emplace(base_path[i], val);
                added_rows_.push_back(row);
                added_nodes_.push_back(base_path[i]);
            }
        } else if (find_a == labels_.end() && find_b != labels_.end()) {
            find_a = labels_.emplace(path[i], find_b->second).first;
            if (find_b->second.second == nannot) {
                added_rows_.push_back(row);
                added_nodes_.push_back(path[i]);
            }
        } else if (find_a != labels_.end() && find_b == labels_.end()) {
            find_b = labels_.emplace(base_path[i], find_a->second).first;
        } else {
            size_t label_i = std::min(find_a->second.second, find_b->second.second);
            if (label_i != nannot) {
                find_a.value().second = label_i;
                find_b.value().second = label_i;
            }
        }
    }

    return std::make_pair(std::move(base_path), reverse_path);
}

void AnnotationBuffer::flush() {
    assert(added_rows_.size() == added_nodes_.size());
    if (added_rows_.empty())
        return;

    // TODO: this cascade of unwrapping is ugly, find a better way to do this
    const DeBruijnGraph *base_graph = &graph_;
    if (const auto *rc_dbg = dynamic_cast<const RCDBG*>(base_graph))
        base_graph = &rc_dbg->get_graph();

    const auto *canonical = dynamic_cast<const CanonicalDBG*>(base_graph);

    auto push_node_labels = [&](auto node_it, auto row_it, auto&& labels) {
        assert(node_it != added_nodes_.end());
        assert(row_it != added_rows_.end());
        assert(labels_.count(*node_it));
        assert(labels_.count(AnnotatedDBG::anno_to_graph_index(*row_it)));
        assert(labels_[*node_it].first
            == *(added_rows_.begin() + (node_it - added_nodes_.begin())));

        size_t label_i = emplace_label_set(std::forward<decltype(labels)>(labels));
        node_index base_node = AnnotatedDBG::anno_to_graph_index(*row_it);
        labels_[*node_it].second = label_i;
        if (base_node != *node_it) {
            auto find = labels_.find(base_node);
            if (find == labels_.end()) {
                labels_[base_node] = std::make_pair(*row_it, label_i);
                if (multi_int_)
                    label_coords_.emplace_back(label_coords_.back());
            }
        } else if (canonical) {
            node_index rc_node = canonical->reverse_complement(*node_it);
            auto find = labels_.find(rc_node);
            if (find == labels_.end()) {
                labels_[rc_node] = std::make_pair(*row_it, label_i);
                if (multi_int_)
                    label_coords_.emplace_back(label_coords_.back());
            }
        }
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

#ifndef NDEBUG
    for (const auto &[node, val] : labels_) {
        assert(val.second != nannot);
    }
#endif
}


template <class T1, class T2>
bool overlap_with_diff(const T1 &tuple1, const T2 &tuple2, ssize_t diff) {
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

struct CoordIntersection {
    CoordIntersection(ssize_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin + offset_ < *b_begin) {
                ++a_begin;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                *out = *a_begin;
                ++a_begin;
                ++b_begin;
                ++out;
            }
        }
    }

    ssize_t offset_;
};

void LabeledExtender::flush() {
    labeled_graph_.flush();
    for ( ; last_flushed_table_i_ < table.size(); ++last_flushed_table_i_) {
        size_t parent_i = std::get<4>(table[last_flushed_table_i_]);
        assert(parent_i < last_flushed_table_i_);

        auto clear = [&]() {
            node_labels_[last_flushed_table_i_] = 0;
            auto &S = std::get<0>(table[last_flushed_table_i_]);
            auto &E = std::get<1>(table[last_flushed_table_i_]);
            auto &F = std::get<2>(table[last_flushed_table_i_]);
            std::fill(S.begin(), S.end(), config_.ninf);
            std::fill(E.begin(), E.end(), config_.ninf);
            std::fill(F.begin(), F.end(), config_.ninf);
        };

        if (!node_labels_[parent_i]) {
            clear();
            continue;
        }

        if (node_labels_[parent_i] != node_labels_[last_flushed_table_i_])
            continue;

        node_index node = std::get<3>(table[last_flushed_table_i_]);
        const auto &parent_labels
            = labeled_graph_.get_labels_from_index(node_labels_[parent_i]);

        auto cur_labels = labeled_graph_.get_labels(node);
        assert(cur_labels);
        Vector<Column> intersect_labels;
        std::set_intersection(parent_labels.begin(),
                              parent_labels.end(),
                              cur_labels->get().begin(),
                              cur_labels->get().end(),
                              std::back_inserter(intersect_labels));
        if (intersect_labels.empty()) {
            clear();
            continue;
        }

        node_labels_[last_flushed_table_i_]
            = labeled_graph_.emplace_label_set(std::move(intersect_labels));
    }
}

bool LabeledExtender::set_seed(const Alignment &seed) {
    if (DefaultColumnExtender::set_seed(seed)) {
        assert(labeled_graph_.is_flushed(seed.get_nodes()));
        fetched_label_i_ = labeled_graph_.emplace_label_set(seed.label_columns);
        assert(fetched_label_i_ != AnnotationBuffer::nannot);
        node_labels_.assign(1, fetched_label_i_);
        base_coords_ = seed.label_coordinates;
        if (base_coords_.size()) {
            if (dynamic_cast<const RCDBG*>(graph_)) {
                for (auto &coords : base_coords_) {
                    for (auto &c : coords) {
                        c += seed.get_nodes().size() - 1 - seed.get_offset();
                    }
                }
            } else {
                if (seed.get_offset()) {
                    for (auto &coords : base_coords_) {
                        for (auto &c : coords) {
                            c -= seed.get_offset();
                        }
                    }
                }
            }
        }

        return true;
    }

    return false;
}

void LabeledExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    bool in_seed = std::get<6>(table[table_i]) + 1 - this->seed_->get_offset()
                        < this->seed_->get_sequence().size();
    size_t next_offset = std::get<6>(table[table_i]) + 1;
    assert(node == std::get<3>(table[table_i]));

    if (in_seed && (next_offset < graph_->get_k() || force_fixed_seed || this->fixed_seed())) {
        assert(labeled_graph_.is_flushed(node));
        DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                             [&](node_index next, char c, score_t score) {
            assert(labeled_graph_.is_flushed(next));
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

    if (outgoing.size() == 1) {
        const auto &[next, c, score] = outgoing[0];
        node_labels_.emplace_back(node_labels_[table_i]);
        callback(next, c, score);
        return;
    }

    flush();
    assert(labeled_graph_.is_flushed(node));
#ifndef NDEBUG
    for (const auto &[next, c, score] : outgoing) {
        assert(labeled_graph_.is_flushed(next));
    }
#endif

    const Alignment &seed = *this->seed_;

    auto [base_labels, base_coords] = labeled_graph_.get_labels_and_coordinates(node);
    assert(base_labels);
    assert(base_labels->get().size());
    const auto &node_labels = labeled_graph_.get_labels_from_index(node_labels_[table_i]);

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


    // if the seed has coordinates, use the seed as the base
    base_labels = std::cref(seed.label_columns);
    base_coords = std::cref(base_coords_);
    ssize_t offset = std::get<6>(table[table_i]);
    ssize_t dist_from_origin = offset - (seed.get_offset() - 1);
    ssize_t dist = dist_from_origin - seed.get_sequence().size() + seed.get_nodes().size();

    // if we are traversing backwards, then negate the coordinate delta
    if (dynamic_cast<const RCDBG*>(graph_))
        dist *= -1;

    for (const auto &[next, c, score] : outgoing) {
        auto [next_labels, next_coords]
            = labeled_graph_.get_labels_and_coordinates(next);

        assert(next_coords);

        // check if at least one label has consistent coordinates
        Vector<Column> intersect_labels;
        auto seed_label_it = node_labels.begin();
        auto seed_label_end_it = node_labels.end();
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

                    if (c < *seed_label_it)
                        return;

                    if (overlap_with_diff(coords, other_coords, dist))
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

bool LabeledExtender::skip_backtrack_start(size_t i) {
    assert(fetched_label_i_ != AnnotationBuffer::nannot);
    assert(node_labels_[i] != AnnotationBuffer::nannot);
    if (!this->prev_starts.emplace(i).second)
        return true;

    const auto &end_labels = labeled_graph_.get_labels_from_index(node_labels_[i]);
    const auto &left_labels = labeled_graph_.get_labels_from_index(fetched_label_i_);

    label_intersection_ = Vector<Column>{};
    label_diff_ = Vector<Column>{};
    set_intersection_difference(left_labels.begin(), left_labels.end(),
                                end_labels.begin(), end_labels.end(),
                                std::back_inserter(label_intersection_),
                                std::back_inserter(label_diff_));
    label_diff_.push_back(AnnotationBuffer::nannot);

    return label_intersection_.empty();
}

void LabeledExtender::call_alignments(score_t cur_cell_score,
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

        auto [base_labels, base_coords]
            = labeled_graph_.get_labels_and_coordinates(alignment.get_nodes().front());
        assert(base_labels);
        assert(base_labels->get().size());

        ssize_t dist = alignment.get_nodes().size() - 1;
        if (!clipping) {
            base_labels = std::cref(seed_->label_columns);
            if (base_coords) {
                base_coords = std::cref(seed_->label_coordinates);
                dist -= seed_->get_offset();
                if (dynamic_cast<const RCDBG*>(graph_)) {
                    dist = (alignment.get_sequence().size() - seed_->get_sequence().size()) * -1;
                }
            }
        } else if (dynamic_cast<const RCDBG*>(graph_)) {
            dist *= -1;
        }

        auto update_fetched = [&]() {
            if (label_diff_.size() && label_diff_.back() == AnnotationBuffer::nannot) {
                label_diff_.pop_back();
                fetched_label_i_ = labeled_graph_.emplace_label_set(std::move(label_diff_));
                assert(fetched_label_i_ != AnnotationBuffer::nannot);
                label_diff_ = Vector<Column>{};
            }
        };

        if (!base_coords) {
            alignment.label_columns = std::move(label_intersection_);
            update_fetched();
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
            auto [cur_labels, cur_coords]
                = labeled_graph_.get_labels_and_coordinates(alignment.get_nodes().back());
            assert(cur_labels);
            assert(cur_labels->get().size());
            assert(cur_coords);
            assert(cur_coords->get().size());
            CoordIntersection intersect_coords(dist);
            try {
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

                        if (c < *seed_label_it)
                            return;

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

        update_fetched();

        callback(std::move(alignment));
    });
}

template <class AlignmentCompare>
void ILabeledAligner<AlignmentCompare>::filter_seeds(BatchSeeders &seeders) const {
    common::logger->trace("Filtering seeds by label");
    std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds;
    std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds_rc;

    size_t num_seeds = 0;
    size_t num_seeds_rc = 0;

    for (auto &[seeder, nodes, seeder_rc, nodes_rc] : seeders) {
        counted_seeds.emplace_back(seeder->get_seeds(), seeder->get_num_matches());
        num_seeds += counted_seeds.back().first.size();

        auto add_seeds = [&](const auto &seeds) {
            for (const Alignment &seed : seeds) {
                labeled_graph_.add_path(
                    seed.get_nodes(),
                    std::string(seed.get_offset(), '#') + seed.get_sequence()
                );
            }
        };

        add_seeds(counted_seeds.back().first);

#if ! _PROTEIN_GRAPH
        if (seeder_rc) {
            counted_seeds_rc.emplace_back(seeder_rc->get_seeds(),
                                          seeder_rc->get_num_matches());
            num_seeds_rc += counted_seeds_rc.back().first.size();
            add_seeds(counted_seeds_rc.back().first);
        }
#endif
    }

    labeled_graph_.flush();

    size_t num_seeds_left = 0;
    size_t num_seeds_rc_left = 0;

    for (size_t i = 0; i < counted_seeds.size(); ++i) {
        auto &[seeder, nodes, seeder_rc, nodes_rc] = seeders[i];
        auto &[seeds, num_matching] = counted_seeds[i];
        if (seeds.size()) {
            num_matching = filter_seeds(seeds);
            num_seeds_left += seeds.size();
        }

        seeder = make_shared<ManualSeeder>(std::move(seeds), num_matching);

#if ! _PROTEIN_GRAPH
        if (seeder_rc) {
            assert(seeder_rc);
            auto &[seeds, num_matching] = counted_seeds_rc[i];
            if (seeds.size()) {
                num_matching = filter_seeds(seeds);
                num_seeds_rc_left += seeds.size();
            }

            seeder_rc = make_shared<ManualSeeder>(std::move(seeds), num_matching);
        }
#endif
    }

    common::logger->trace("Kept {}/{} seeds",
                          num_seeds_left + num_seeds_rc_left,
                          num_seeds + num_seeds_rc);

    common::logger->trace("Prefetching labels");
    labeled_graph_.flush();
}

inline size_t get_num_matches(const std::vector<Alignment> &seeds) {
    size_t num_matching = 0;
    size_t last_end = 0;
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (seeds[i].empty())
            continue;

        size_t begin = seeds[i].get_clipping();
        size_t end = begin + seeds[i].get_query().size();
        if (end > last_end) {
            num_matching += end - begin;
            if (begin < last_end)
                num_matching -= last_end - begin;
        }

        if (size_t offset = seeds[i].get_offset()) {
            size_t clipping = seeds[i].get_clipping();
            for (++i; i < seeds.size()
                        && seeds[i].get_offset() == offset
                        && seeds[i].get_clipping() == clipping; ++i) {}
            --i;
        }

        last_end = end;
    }
    return num_matching;
}

template <class AlignmentCompare>
size_t ILabeledAligner<AlignmentCompare>
::filter_seeds(std::vector<Alignment> &seeds) const {
    VectorMap<Column, uint64_t> label_counter;
    for (const Alignment &seed : seeds) {
        assert(labeled_graph_.is_flushed(seed.get_nodes()));
        for (node_index node : seed.get_nodes()) {
            auto labels = labeled_graph_.get_labels(node);
            assert(labels);
            for (uint64_t label : labels->get()) {
                ++label_counter[label];
            }
        }
    }

    if (label_counter.empty())
        return get_num_matches(seeds);

    std::vector<std::pair<Column, uint64_t>> label_counts
        = const_cast<std::vector<std::pair<Column, uint64_t>>&&>(
            label_counter.values_container()
        );

    std::sort(label_counts.begin(), label_counts.end(), utils::GreaterSecond());

    double cutoff = static_cast<double>(label_counts[0].second)
        * this->config_.rel_score_cutoff;
    auto it = std::find_if(label_counts.begin() + 1, label_counts.end(),
                           [cutoff](const auto &a) { return a.second < cutoff; });

    label_counts.erase(it, label_counts.end());

    Vector<Column> labels;
    labels.reserve(label_counts.size());
    for (const auto &[label, count] : label_counts) {
        labels.push_back(label);
    }
    std::sort(labels.begin(), labels.end());

    for (Alignment &seed : seeds) {
        const std::vector<node_index> &nodes = seed.get_nodes();
        auto [fetch_labels, fetch_coords]
            = labeled_graph_.get_labels_and_coordinates(nodes[0]);
        if (!fetch_labels)
            continue;

        if (fetch_coords) {
            auto a_begin = fetch_labels->get().begin();
            auto a_end = fetch_labels->get().end();
            auto a_c_begin = fetch_coords->get().begin();
            auto b_begin = labels.begin();
            auto b_end = labels.end();
            while (a_begin != a_end && b_begin != b_end) {
                if (*a_begin < *b_begin) {
                    ++a_begin;
                    ++a_c_begin;
                } else if (*a_begin > *b_begin) {
                    ++b_begin;
                } else {
                    seed.label_columns.emplace_back(*a_begin);
                    seed.label_coordinates.emplace_back(*a_c_begin);
                    ++a_begin;
                    ++a_c_begin;
                    ++b_begin;
                }
            }
        } else {
            std::set_intersection(fetch_labels->get().begin(),
                                  fetch_labels->get().end(),
                                  labels.begin(), labels.end(),
                                  std::back_inserter(seed.label_columns));
        }

        seed.label_encoder = &labeled_graph_.get_annotator().get_label_encoder();

        for (size_t i = 1; i < nodes.size() && seed.label_columns.size(); ++i) {
            auto [next_fetch_labels, next_fetch_coords]
                = labeled_graph_.get_labels_and_coordinates(nodes[i]);
            if (next_fetch_coords) {
                Vector<Column> label_inter;
                Alignment::CoordinateSet coord_inter;
                CoordIntersection intersect_coords(i);
                utils::match_indexed_values(
                    seed.label_columns.begin(), seed.label_columns.end(),
                    seed.label_coordinates.begin(),
                    next_fetch_labels->get().begin(), next_fetch_labels->get().end(),
                    next_fetch_coords->get().begin(),
                    [&](auto col, const auto &coords, const auto &other_coords) {
                        Alignment::Tuple overlap;
                        intersect_coords(coords.begin(), coords.end(),
                                         other_coords.begin(), other_coords.end(),
                                         std::back_inserter(overlap));
                        if (overlap.size()) {
                            label_inter.push_back(col);
                            coord_inter.push_back(std::move(overlap));
                        }
                    }
                );

                std::swap(seed.label_columns, label_inter);
                std::swap(seed.label_coordinates, coord_inter);
            } else if (next_fetch_labels) {
                Vector<Column> temp;
                std::set_intersection(next_fetch_labels->get().begin(),
                                      next_fetch_labels->get().end(),
                                      seed.label_columns.begin(),
                                      seed.label_columns.end(),
                                      std::back_inserter(temp));
                std::swap(temp, seed.label_columns);
            } else {
                seed.label_columns.clear();
                seed.label_coordinates.clear();
            }
        }

        labeled_graph_.emplace_label_set(seed.label_columns);

        if (seed.get_offset() && seed.label_coordinates.size()) {
            for (auto &tuple : seed.label_coordinates) {
                for (auto &coord : tuple) {
                    coord += seed.get_offset();
                }
            }
        }
    }

    auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [&](const auto &a) {
        assert(labeled_graph_.get_index(a.label_columns) != AnnotationBuffer::nannot);
        return a.label_columns.empty();
    });

    seeds.erase(seed_it, seeds.end());

    return get_num_matches(seeds);
}

template class ILabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg