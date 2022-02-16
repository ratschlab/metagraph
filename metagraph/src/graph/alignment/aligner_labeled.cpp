#include "aligner_labeled.hpp"

#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

AnnotationBuffer::AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator)
      : graph_(graph),
        annotator_(annotator),
        multi_int_(dynamic_cast<const annot::matrix::MultiIntMatrix*>(&annotator_.get_matrix())),
        labels_set_({ {} }) {
    if (multi_int_ && graph_.get_mode() == DeBruijnGraph::CANONICAL) {
        multi_int_ = nullptr;
        logger->warn("Coordinates not supported when aligning to CANONICAL "
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
        Row row = AnnotatedDBG::graph_to_anno_index(base_path[i]);
        std::pair<Row, size_t> val { row, 0 };

        if (base_path[i] == DeBruijnGraph::npos) {
            // this can happen when the base graph is CANONICAL and path[i] is a
            // dummy node
            if (labels_.emplace(path[i], val).second && multi_int_)
                label_coords_.emplace_back();

            continue;
        }

        if (boss && !boss->get_W(dbg_succ->kmer_to_boss_index(base_path[i]))) {
            // skip dummy nodes
            if (labels_.emplace(path[i], val).second && multi_int_)
                label_coords_.emplace_back();

            if (labels_.emplace(base_path[i], val).second && multi_int_)
                label_coords_.emplace_back();

            continue;
        }

        auto find_a = labels_.find(path[i]);
        auto find_b = labels_.find(base_path[i]);

        if (find_a == labels_.end() && find_b == labels_.end()) {
            val.second = nannot;
            labels_.emplace(path[i], val);
            added_rows_.push_back(row);
            added_nodes_.push_back(path[i]);

            if (path[i] != base_path[i]) {
                labels_.emplace(base_path[i], val);
                added_rows_.push_back(row);
                added_nodes_.push_back(base_path[i]);
            }
        } else if (find_a == labels_.end() && find_b != labels_.end()) {
            labels_.emplace(path[i], find_b->second);
            if (find_b->second.second == nannot) {
                added_rows_.push_back(row);
                added_nodes_.push_back(path[i]);
            }
        } else if (find_a != labels_.end() && find_b == labels_.end()) {
            labels_.emplace(base_path[i], find_a->second);
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
            LabelSet labels;
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

LabeledExtender::LabeledExtender(const IDBGAligner &aligner, std::string_view query)
        : LabeledExtender(dynamic_cast<const LabeledAligner<>&>(aligner).get_annotation_buffer(),
                          aligner.get_config(), query) {}

void LabeledExtender::flush() {
    annotation_buffer_.flush();
    for ( ; last_flushed_table_i_ < table.size(); ++last_flushed_table_i_) {
        auto &table_elem = table[last_flushed_table_i_];

        size_t parent_i = table_elem.parent_i;
        assert(parent_i < last_flushed_table_i_);

        auto clear = [&]() {
            node_labels_[last_flushed_table_i_] = 0;
            std::fill(table_elem.S.begin(), table_elem.S.end(), config_.ninf);
            std::fill(table_elem.E.begin(), table_elem.E.end(), config_.ninf);
            std::fill(table_elem.F.begin(), table_elem.F.end(), config_.ninf);
        };

        if (!node_labels_[parent_i]) {
            clear();
            continue;
        }

        if (node_labels_[parent_i] != node_labels_[last_flushed_table_i_])
            continue;

        node_index node = table_elem.node;
        const auto &parent_labels
            = annotation_buffer_.get_labels_from_index(node_labels_[parent_i]);

        auto cur_labels = annotation_buffer_.get_labels(node);
        assert(cur_labels);
        LabelSet intersect_labels;
        std::set_intersection(parent_labels.begin(), parent_labels.end(),
                              cur_labels->begin(), cur_labels->end(),
                              std::back_inserter(intersect_labels));
        if (intersect_labels.empty()) {
            clear();
        } else {
            node_labels_[last_flushed_table_i_]
                = annotation_buffer_.emplace_label_set(std::move(intersect_labels));
        }
    }
}

bool LabeledExtender::set_seed(const Alignment &seed) {
    if (!DefaultColumnExtender::set_seed(seed))
        return false;

    assert(std::all_of(seed.get_nodes().begin(), seed.get_nodes().end(),
                       [&](node_index n) { return annotation_buffer_.get_labels(n); }));

    // the first node of the seed has already been flushed
    last_flushed_table_i_ = 1;

    remaining_labels_i_ = annotation_buffer_.emplace_label_set(seed.label_columns);
    assert(remaining_labels_i_ != AnnotationBuffer::nannot);
    node_labels_.assign(1, remaining_labels_i_);
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

void LabeledExtender
::call_outgoing(node_index node,
                size_t max_prefetch_distance,
                const std::function<void(node_index, char, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    assert(table.size() == node_labels_.size());

    // if we are in the seed and want to force the seed to be fixed, automatically
    // take the next node in the seed
    size_t next_offset = table[table_i].offset + 1;
    bool in_seed = (next_offset - seed_->get_offset() < seed_->get_sequence().size())
        && (next_offset < graph_->get_k() || force_fixed_seed || fixed_seed());

    std::vector<std::tuple<node_index, char, score_t>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
                                         [&](node_index next, char c, score_t score) {
        outgoing.emplace_back(next, c, score);
        if (!in_seed)
            annotation_buffer_.add_node(next);
    }, table_i, force_fixed_seed);

    if (outgoing.empty())
        return;

    if (outgoing.size() == 1) {
        // Assume that annotations are preserved in unitigs. Violations of this
        // assumption are corrected after the next flush
        const auto &[next, c, score] = outgoing[0];
        node_labels_.emplace_back(node_labels_[table_i]);
        callback(next, c, score);
        return;
    }

    // flush the AnnotationBuffer and correct for annotation errors introduced above
    flush();

    if (!node_labels_[table_i])
        return;

    assert(annotation_buffer_.get_labels(node));

    // use the label set of the current node in the alignment tree as the basis
    const auto &node_labels = annotation_buffer_.get_labels_from_index(node_labels_[table_i]);

    // no coordinates are present in the annotation
    if (!annotation_buffer_.get_labels_and_coordinates(node).second) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        for (const auto &[next, c, score] : outgoing) {
            auto next_labels = annotation_buffer_.get_labels(next);
            assert(next_labels);

            LabelSet intersect_labels;
            std::set_intersection(node_labels.begin(), node_labels.end(),
                                  next_labels->begin(), next_labels->end(),
                                  std::back_inserter(intersect_labels));

            if (intersect_labels.size()) {
                node_labels_.emplace_back(annotation_buffer_.emplace_label_set(
                    std::move(intersect_labels)
                ));
                callback(next, c, score);
            }
        }

        return;
    }

    // check label and coordinate consistency
    // use the seed as the basis for labels and coordinates
    assert(seed_->label_coordinates.size());

    // compute the coordinate distance from base_coords
    size_t dist = next_offset - graph_->get_k() + 1;

    for (const auto &[next, c, score] : outgoing) {
        const AnnotationBuffer::LabelSet *base_labels = &seed_->label_columns;
        const AnnotationBuffer::CoordinateSet *base_coords = &base_coords_;
        auto [next_labels, next_coords]
            = annotation_buffer_.get_labels_and_coordinates(next);

        assert(next_coords);

        // if we are traversing backwards, then negate the coordinate delta
        if (dynamic_cast<const RCDBG*>(graph_)) {
            std::swap(base_labels, next_labels);
            std::swap(base_coords, next_coords);
        }

        // check if at least one label has consistent coordinates
        LabelSet intersect_labels;

        auto cur_label_it = node_labels.begin();
        auto cur_label_end_it = node_labels.end();
        try {
            utils::match_indexed_values(
                base_labels->begin(), base_labels->end(), base_coords->begin(),
                next_labels->begin(), next_labels->end(), next_coords->begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    // also, intersect with the label set of node
                    while (cur_label_it != cur_label_end_it && c > *cur_label_it) {
                        ++cur_label_it;
                    }

                    if (cur_label_it == cur_label_end_it)
                        throw std::exception();

                    if (c < *cur_label_it)
                        return;

                    // then check coordinate consistency with the seed
                    if (overlap_with_diff(coords, other_coords, dist))
                        intersect_labels.push_back(c);
                }
            );
        } catch (const std::exception&) {}

        if (intersect_labels.size()) {
            // found a consistent set of coordinates, so assign labels for this next node
            node_labels_.emplace_back(annotation_buffer_.emplace_label_set(
                std::move(intersect_labels)
            ));
            callback(next, c, score);
        }
    }
}

bool LabeledExtender::skip_backtrack_start(size_t i) {
    assert(remaining_labels_i_ != AnnotationBuffer::nannot);
    assert(node_labels_[i] != AnnotationBuffer::nannot);

    // if this alignment tree node has been visited previously, ignore it
    if (!prev_starts.emplace(i).second)
        return true;

    // check if this starting point involves seed labels which have not been considered yet
    const auto &end_labels = annotation_buffer_.get_labels_from_index(node_labels_[i]);
    const auto &left_labels = annotation_buffer_.get_labels_from_index(remaining_labels_i_);

    label_intersection_ = LabelSet{};
    label_diff_ = LabelSet{};
    set_intersection_difference(left_labels.begin(), left_labels.end(),
                                end_labels.begin(), end_labels.end(),
                                std::back_inserter(label_intersection_),
                                std::back_inserter(label_diff_));
    label_diff_.push_back(AnnotationBuffer::nannot);

    return label_intersection_.empty();
}

void LabeledExtender::call_alignments(score_t end_score,
                                      const std::vector<node_index> &path,
                                      const std::vector<size_t> & /* trace */,
                                      const Cigar &ops,
                                      size_t clipping,
                                      size_t offset,
                                      std::string_view window,
                                      const std::string &match,
                                      score_t extra_penalty,
                                      const std::function<void(Alignment&&)> &callback) {
    Alignment alignment = construct_alignment(ops, clipping, window, path, match,
                                              end_score, offset, extra_penalty);
    alignment.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();

    auto [base_labels, base_coords]
        = annotation_buffer_.get_labels_and_coordinates(alignment.get_nodes().front());
    assert(base_labels);
    assert(base_labels->size());

    if (!clipping)
        base_labels = &seed_->label_columns;

    auto call_alignment = [&]() {
        if (label_diff_.size() && label_diff_.back() == AnnotationBuffer::nannot) {
            label_diff_.pop_back();
            remaining_labels_i_ = annotation_buffer_.emplace_label_set(std::move(label_diff_));
            assert(remaining_labels_i_ != AnnotationBuffer::nannot);
            label_diff_ = LabelSet{};
        }

        callback(std::move(alignment));
    };

    if (!base_coords) {
        alignment.label_columns = std::move(label_intersection_);
        call_alignment();
        return;
    }

    ssize_t dist = alignment.get_nodes().size() - 1;
    if (!clipping) {
        base_coords = &seed_->label_coordinates;
        dist -= seed_->get_offset();
        if (dynamic_cast<const RCDBG*>(graph_))
            dist = alignment.get_sequence().size() - seed_->get_sequence().size();
    }

    auto label_it = label_intersection_.begin();
    auto label_end_it = label_intersection_.end();

    if (alignment.get_nodes().size() == 1) {
        auto it = base_labels->begin();
        auto end = base_labels->end();
        auto c_it = base_coords->begin();
        while (label_it != label_end_it && it != end) {
            if (*label_it < *it) {
                ++label_it;
            } else if (*label_it > *it) {
                ++it;
                ++c_it;
            } else {
                alignment.label_columns.emplace_back(*it);
                alignment.label_coordinates.emplace_back(*c_it);
                ++it;
                ++c_it;
                ++label_it;
            }
        }
    } else {
        auto [cur_labels, cur_coords]
            = annotation_buffer_.get_labels_and_coordinates(alignment.get_nodes().back());
        assert(cur_labels);
        assert(cur_labels->size());
        assert(cur_coords);
        assert(cur_coords->size());
        if (dynamic_cast<const RCDBG*>(graph_)) {
            std::swap(cur_labels, base_labels);
            std::swap(cur_coords, base_coords);
        }
        CoordIntersection intersect_coords(dist);
        try {
            utils::match_indexed_values(
                base_labels->begin(), base_labels->end(), base_coords->begin(),
                cur_labels->begin(), cur_labels->end(), cur_coords->begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    while (label_it != label_end_it && c > *label_it) {
                        ++label_it;
                    }

                    if (label_it == label_end_it)
                        throw std::exception();

                    if (c < *label_it)
                        return;

                    Alignment::Tuple overlap;
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

    if (dynamic_cast<const RCDBG*>(graph_) && alignment.get_offset()) {
        for (auto &coords : alignment.label_coordinates) {
            for (auto &c : coords) {
                c += alignment.get_offset();
            }
        }
    }

    call_alignment();
}

template <class Seeder, class Extender, class AlignmentCompare>
LabeledAligner<Seeder, Extender, AlignmentCompare>
::LabeledAligner(const DeBruijnGraph &graph,
                 const DBGAlignerConfig &config,
                 const Annotator &annotator)
      : DBGAligner<Seeder, Extender, AlignmentCompare>(graph, config),
        annotation_buffer_(graph, annotator) {
    if (annotation_buffer_.get_coordinate_matrix()
            && std::is_same_v<Extender, LabeledExtender>) {
        // do not use a global xdrop cutoff since we need separate cutoffs
        // for each label
        this->config_.global_xdrop = false;
    }
}

template <class Seeder, class Extender, class AlignmentCompare>
LabeledAligner<Seeder, Extender, AlignmentCompare>::~LabeledAligner() {
    logger->trace("Buffered {}/{} nodes and {} label combinations",
                  annotation_buffer_.num_nodes_buffered(),
                  this->graph_.num_nodes(),
                  annotation_buffer_.num_label_sets());
}

template <class Seeder, class Extender, class AlignmentCompare>
auto LabeledAligner<Seeder, Extender, AlignmentCompare>
::build_seeders(const std::vector<IDBGAligner::Query> &seq_batch,
                const std::vector<AlignmentResults> &wrapped_seqs) const -> BatchSeeders {
    BatchSeeders seeders
        = DBGAligner<Seeder, Extender, AlignmentCompare>::build_seeders(seq_batch, wrapped_seqs);

    // now we're going to filter the seeds
    logger->trace("Filtering seeds by label");
    std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds;
    std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds_rc;

    size_t num_seeds = 0;
    size_t num_seeds_rc = 0;

    for (auto &[seeder, nodes, seeder_rc, nodes_rc] : seeders) {
        counted_seeds.emplace_back(seeder->get_seeds(), seeder->get_num_matches());
        num_seeds += counted_seeds.back().first.size();

        auto add_seeds = [&](const auto &seeds) {
            for (const Alignment &seed : seeds) {
                annotation_buffer_.add_path(
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

    annotation_buffer_.flush();

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

    logger->trace("Kept {}/{} seeds", num_seeds_left + num_seeds_rc_left,
                                      num_seeds + num_seeds_rc);
    logger->trace("Prefetching labels");
    annotation_buffer_.flush();

    return seeders;
}

inline size_t get_num_matches(const std::vector<Alignment> &seeds) {
    size_t num_matching = 0;
    size_t last_end = 0;
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (seeds[i].empty())
            continue;

        size_t begin = seeds[i].get_clipping();
        size_t end = begin + seeds[i].get_query_view().size();
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

template <class Seeder, class Extender, class AlignmentCompare>
size_t LabeledAligner<Seeder, Extender, AlignmentCompare>
::filter_seeds(std::vector<Alignment> &seeds) const {
    VectorMap<Column, uint64_t> label_counter;
    for (const Alignment &seed : seeds) {
        for (node_index node : seed.get_nodes()) {
            assert(annotation_buffer_.get_labels(node));
            for (uint64_t label : *annotation_buffer_.get_labels(node)) {
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

    LabelSet labels;
    labels.reserve(label_counts.size());
    for (const auto &[label, count] : label_counts) {
        labels.push_back(label);
    }
    std::sort(labels.begin(), labels.end());

    for (Alignment &seed : seeds) {
        const std::vector<node_index> &nodes = seed.get_nodes();
        auto [fetch_labels, fetch_coords]
            = annotation_buffer_.get_labels_and_coordinates(nodes[0]);
        if (!fetch_labels)
            continue;

        if (fetch_coords) {
            auto a_begin = fetch_labels->begin();
            auto a_end = fetch_labels->end();
            auto a_c_begin = fetch_coords->begin();
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
            std::set_intersection(fetch_labels->begin(), fetch_labels->end(),
                                  labels.begin(), labels.end(),
                                  std::back_inserter(seed.label_columns));
        }

        seed.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();

        for (size_t i = 1; i < nodes.size() && seed.label_columns.size(); ++i) {
            auto [next_fetch_labels, next_fetch_coords]
                = annotation_buffer_.get_labels_and_coordinates(nodes[i]);
            if (next_fetch_coords) {
                LabelSet label_inter;
                Alignment::CoordinateSet coord_inter;
                CoordIntersection intersect_coords(i);
                utils::match_indexed_values(
                    seed.label_columns.begin(), seed.label_columns.end(),
                    seed.label_coordinates.begin(),
                    next_fetch_labels->begin(), next_fetch_labels->end(),
                    next_fetch_coords->begin(),
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
                LabelSet temp;
                std::set_intersection(next_fetch_labels->begin(), next_fetch_labels->end(),
                                      seed.label_columns.begin(), seed.label_columns.end(),
                                      std::back_inserter(temp));
                std::swap(temp, seed.label_columns);
            } else {
                seed.label_columns.clear();
                seed.label_coordinates.clear();
            }
        }

        annotation_buffer_.emplace_label_set(seed.label_columns);

        if (seed.get_offset() && seed.label_coordinates.size()) {
            for (auto &tuple : seed.label_coordinates) {
                for (auto &coord : tuple) {
                    coord += seed.get_offset();
                }
            }
        }
    }

    auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [&](const auto &a) {
        assert(annotation_buffer_.get_index(a.label_columns) != AnnotationBuffer::nannot);
        return a.label_columns.empty();
    });

    seeds.erase(seed_it, seeds.end());

    return get_num_matches(seeds);
}

template class LabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
