#include "aligner_labeled.hpp"

#include <progress_bar.hpp>

#include "graph/representation/rc_dbg.hpp"
#include "common/algorithms.hpp"
#include "common/unix_tools.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

typedef annot::binmat::BinaryMatrix::Row Row;
typedef annot::binmat::BinaryMatrix::Column Column;
typedef AnnotationBuffer::Columns Columns;
typedef DeBruijnGraph::node_index node_index;

// dummy index for an unfetched annotations
static constexpr size_t nannot = std::numeric_limits<size_t>::max();

template <class T1, class T2>
bool overlap_with_diff(const T1 &tuple1, const T2 &tuple2, int64_t diff) {
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

template <class AIt, class BIt, class CallbackInter, class CallbackDiff1, class CallbackDiff2>
void call_set_intersection_difference(AIt a_begin, AIt a_end,
                                      BIt b_begin, BIt b_end,
                                      const CallbackInter &callback_inter,
                                      const CallbackDiff1 &callback_diff1,
                                      const CallbackDiff2 &callback_diff2) {
    while (a_begin != a_end || b_begin != b_end) {
        if (b_begin == b_end) {
            callback_diff1(*a_begin);
            ++a_begin;
        } else if (a_begin == a_end || *a_begin > *b_begin) {
            callback_diff2(*b_begin);
            ++b_begin;
        } else {
            if (*a_begin == *b_begin) {
                callback_inter(*a_begin);
                ++b_begin;
            } else {
                callback_diff1(*a_begin);
            }

            ++a_begin;
        }
    }
}

template <class AIt, class BIt, class OutIt, class OutIt2, class OutIt3>
void set_intersection_difference(AIt a_begin, AIt a_end,
                                 BIt b_begin, BIt b_end,
                                 OutIt intersection_out,
                                 OutIt2 diff_out,
                                 OutIt3 diff_out2) {
    call_set_intersection_difference(a_begin, a_end, b_begin, b_end,
        [&](const auto &inter) { *intersection_out = inter; ++intersection_out; },
        [&](const auto &diff1) { *diff_out = diff1; ++diff_out; },
        [&](const auto &diff2) { *diff_out2 = diff2; ++diff_out2; }
    );
}

LabeledExtender::LabeledExtender(const IDBGAligner &aligner, std::string_view query)
    : LabeledExtender(aligner.get_graph(), aligner.get_config(),
                      dynamic_cast<const LabeledAligner<>&>(aligner).annotation_buffer_,
                      query) {}

void LabeledExtender::flush(size_t table_i, const std::vector<node_index> &outnodes) {
    if (!table_i) {
        // flush everything
        for (size_t i = table.size(); i > 1; --i) {
            if (node_labels_[i - 1] == nannot)
                flush(i - 1);
        }

        return;
    }

    assert(node_labels_.size() == table.size());

    // construct path
    std::vector<node_index> path;
    std::vector<size_t> table_path;
    while (node_labels_[table_i] == nannot) {
        assert(table[table_i].node != DeBruijnGraph::npos);
        path.push_back(table[table_i].node);
        table_path.push_back(table_i);
        table_i = table[table_i].parent_i;
    }

    assert(table[table_i].node != DeBruijnGraph::npos);
    path.push_back(table[table_i].node);
    std::reverse(path.begin(), path.end());
    annotation_buffer_.queue_path(std::move(path));
    for (node_index next : outnodes) {
        annotation_buffer_.queue_path(1, next);
    }
    annotation_buffer_.fetch_queued_annotations();

    size_t parent_i = table_i;
    for (auto it = table_path.rbegin(); it != table_path.rend(); ++it) {
        size_t table_i = *it;
        auto &table_elem = table[table_i];

        assert(parent_i < table_i);

        auto clear = [&]() {
            DEBUG_LOG("Removed table element {}", table_i);
            node_labels_[table_i] = 0;
            std::fill(table_elem.S.begin(), table_elem.S.end(), config_.ninf);
            std::fill(table_elem.E.begin(), table_elem.E.end(), config_.ninf);
            std::fill(table_elem.F.begin(), table_elem.F.end(), config_.ninf);
        };

        if (!node_labels_[parent_i]) {
            clear();
            std::swap(table_i, parent_i);
            continue;
        }

        if (table_elem.node == DeBruijnGraph::npos) {
            std::swap(table_i, parent_i);
            continue;
        }

        const auto &parent_labels
            = annotation_buffer_.get_cached_column_set(node_labels_[parent_i]);

        auto cur_labels = annotation_buffer_.get_labels(table_elem.node);
        assert(cur_labels);

#ifndef NDEBUG
        if (table[parent_i].offset >= 0
                && static_cast<size_t>(table[parent_i].offset) >= graph_->get_k() - 1) {
            auto parent_real_labels = annotation_buffer_.get_labels(table[parent_i].node);
            assert(parent_real_labels);
            assert(parent_real_labels->size() >= parent_labels.size());
            Columns diff;
            std::set_difference(parent_labels.begin(), parent_labels.end(),
                                parent_real_labels->begin(), parent_real_labels->end(),
                                std::back_inserter(diff));
            assert(diff.empty());
        }
#endif

        Columns intersect_labels;
        std::set_intersection(parent_labels.begin(), parent_labels.end(),
                              cur_labels->begin(), cur_labels->end(),
                              std::back_inserter(intersect_labels));

        if (intersect_labels.empty()) {
            clear();
        } else {
            node_labels_[table_i]
                = annotation_buffer_.cache_column_set(std::move(intersect_labels));
        }

        std::swap(table_i, parent_i);
    }
}

bool LabeledExtender::set_seed(const Alignment &seed) {
    if (!DefaultColumnExtender::set_seed(seed))
        return false;

    assert(std::all_of(seed.get_nodes().begin(), seed.get_nodes().end(),
                    [&](node_index n) {
                        return n == DeBruijnGraph::npos || annotation_buffer_.get_labels(n);
                    }));

    remaining_labels_i_ = annotation_buffer_.cache_column_set(seed.label_columns);
    assert(remaining_labels_i_);
    assert(remaining_labels_i_ != nannot);
    node_labels_.assign(1, remaining_labels_i_);
    base_coords_ = seed.label_coordinates;
    if (!annotation_buffer_.has_coordinates())
        return true;

    if (!base_coords_.size()) {
        base_coords_ = *annotation_buffer_.get_coords(seed.get_nodes()[0],
                                                      false,
                                                      &seed.label_columns);
        for (auto &tuple : base_coords_) {
            for (auto &c : tuple) {
                c += seed.get_offset();
            }
        }
    }

    if (dynamic_cast<const RCDBG*>(graph_)) {
        for (auto &coords : base_coords_) {
            for (auto &c : coords) {
                c += seed.get_nodes().size() - 1 - seed.get_offset();
            }
        }
    } else if (seed.get_offset()) {
        for (auto &coords : base_coords_) {
            for (auto &c : coords) {
                c -= seed.get_offset();
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

    std::vector<std::tuple<node_index, char, score_t>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
        [&](node_index next, char c, score_t score) {
            outgoing.emplace_back(next, c, score);
        },
        table_i, force_fixed_seed
    );

    if (outgoing.empty())
        return;

    size_t next_offset = table[table_i].offset + 1;
    if (outgoing.size() == 1) {
        // Assume that annotations are preserved in unitigs. Violations of this
        // assumption are corrected after the next flush
        const auto &[next, c, score] = outgoing[0];
        bool in_seed = next_offset - seed_->get_offset() < seed_->get_sequence().size()
                        && (next_offset < graph_->get_k() || force_fixed_seed);
        if (in_seed) {
            node_labels_.emplace_back(node_labels_[table_i]);
        } else {
            auto [column_diff, coord_diff] = annotation_buffer_.get_label_and_coord_diff(node, next);
            if (column_diff.empty()) {
                node_labels_.emplace_back(node_labels_[table_i]);
            } else if (coord_diff) {
                Columns next_columns;
                size_t dist = next_offset - graph_->get_k() + 1;
                if (dynamic_cast<const RCDBG*>(graph_))
                    dist *= -1;

                auto b_it = seed_->label_columns.begin();
                auto b_c_it = base_coords_.begin();
                const auto &columns = annotation_buffer_.get_cached_column_set(node_labels_[table_i]);
                auto c_it = columns.begin();

                auto d_it = column_diff.begin();
                auto d_c_it = coord_diff->begin();

                while (c_it != columns.end()) {
                    assert(b_it != seed_->label_columns.end());
                    while (*b_it < *c_it) {
                        ++b_it;
                        ++b_c_it;
                        assert(b_it != seed_->label_columns.end());
                    }
                    assert(*b_it == *c_it);
                    while (d_it != column_diff.end() && *d_it < *c_it) {
                        ++d_it;
                        ++d_c_it;
                    }

                    if (d_it == column_diff.end() || *d_it > *c_it) {
                        next_columns.push_back(*c_it);

                    } else {
                        assert(*d_it == *c_it);
                        assert(d_c_it->size());
                        Alignment::Tuple next_tuple;
                        utils::set_difference(b_c_it->begin(), b_c_it->end(),
                                              d_c_it->begin(), d_c_it->end(),
                                              std::back_inserter(next_tuple),
                                              dist);
                        if (next_tuple.size())
                            next_columns.push_back(*c_it);

                        ++d_it;
                        ++d_c_it;
                    }
                    ++c_it;
                    ++b_it;
                    ++b_c_it;
                }

                if (next_columns.empty())
                    return;

                size_t next_label_i = annotation_buffer_.cache_column_set(std::move(next_columns));
                node_labels_.emplace_back(next_label_i);

            } else {
                const auto &columns = annotation_buffer_.get_cached_column_set(node_labels_[table_i]);
                Columns next_columns;
                std::set_difference(columns.begin(), columns.end(),
                                    column_diff.begin(), column_diff.end(),
                                    std::back_inserter(next_columns));

                if (next_columns.empty())
                    return;

                size_t next_label_i = annotation_buffer_.cache_column_set(std::move(next_columns));
                node_labels_.emplace_back(next_label_i);
            }
        }
        callback(next, c, score);
        return;
    }

    std::vector<node_index> outnodes;
    outnodes.reserve(outgoing.size());
    for (const auto &[next, c, score] : outgoing) {
        assert(next != DeBruijnGraph::npos);
        outnodes.push_back(next);
    }

    // flush the AnnotationBuffer and correct for annotation errors introduced above
    flush(table_i, outnodes);

    if (!node_labels_[table_i])
        return;

    assert(annotation_buffer_.get_labels(node));

    // use the label set of the current node in the alignment tree as the basis
    const auto &columns = annotation_buffer_.get_cached_column_set(node_labels_[table_i]);

    // no coordinates are present in the annotation
    if (!annotation_buffer_.has_coordinates()) {
        // label consistency (weaker than coordinate consistency):
        // checks if there is at least one label shared between adjacent nodes
        for (const auto &[next, c, score] : outgoing) {
            auto next_labels = annotation_buffer_.get_labels(next);
            assert(next_labels);

            Columns intersect_labels;
            std::set_intersection(columns.begin(), columns.end(),
                                  next_labels->begin(), next_labels->end(),
                                  std::back_inserter(intersect_labels));

            if (intersect_labels.size()) {
                node_labels_.push_back(annotation_buffer_.cache_column_set(
                                                std::move(intersect_labels)));
                callback(next, c, score);
            }
        }

        return;
    }

    annotation_buffer_.prefetch_coords(outnodes);

    // check label and coordinate consistency
    // use the seed as the basis for labels and coordinates
    assert(seed_->label_coordinates.size());

    // compute the coordinate distance from base_coords
    size_t dist = next_offset - graph_->get_k() + 1;

    for (const auto &[next, c, score] : outgoing) {
        const Columns *base_labels = &seed_->label_columns;
        std::shared_ptr<const CoordinateSet> base_coords {
            std::shared_ptr<const CoordinateSet>{}, &base_coords_
        };

        auto [next_labels, next_coords_mut]
            = annotation_buffer_.get_labels_and_coords(next, false);

        assert(next_coords_mut);
        std::shared_ptr<const CoordinateSet> next_coords = next_coords_mut;

        // if we are traversing backwards, then negate the coordinate delta
        if (dynamic_cast<const RCDBG*>(graph_)) {
            std::swap(base_labels, next_labels);
            std::swap(base_coords, next_coords);
        }

        // check if at least one label has consistent coordinates
        Columns intersect_labels;

        try {
            auto col_it = columns.begin();
            auto col_end = columns.end();

            utils::match_indexed_values(
                base_labels->begin(), base_labels->end(), base_coords->begin(),
                next_labels->begin(), next_labels->end(), next_coords->begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    // also, intersect with the label set of node
                    while (col_it != col_end && c > *col_it) {
                        ++col_it;
                    }

                    if (col_it == col_end)
                        throw std::exception();

                    if (c < *col_it)
                        return;

                    // then check coordinate consistency with the seed
                    if (overlap_with_diff(coords, other_coords, dist))
                        intersect_labels.push_back(c);
                }
            );
        } catch (const std::exception&) {}

        if (intersect_labels.size()) {
            // found a consistent set of coordinates, so assign labels for this next node
            node_labels_.push_back(annotation_buffer_.cache_column_set(
                                                std::move(intersect_labels)));
            callback(next, c, score);
        }
    }
}

bool LabeledExtender::skip_backtrack_start(size_t i) {
    assert(remaining_labels_i_ != nannot);
    assert(node_labels_[i] != nannot);

    // if this alignment tree node has been visited previously, ignore it
    assert(remaining_labels_i_);
    if (!prev_starts.emplace(i).second)
        return true;

    // check if this starting point involves seed labels which have not been considered yet
    const auto &end_labels = annotation_buffer_.get_cached_column_set(node_labels_[i]);
    const auto &left_labels = annotation_buffer_.get_cached_column_set(remaining_labels_i_);

    label_intersection_ = Columns{};
    label_diff_ = Columns{};
    utils::set_intersection_difference(left_labels.begin(), left_labels.end(),
                                       end_labels.begin(), end_labels.end(),
                                       std::back_inserter(label_intersection_),
                                       std::back_inserter(label_diff_));
    label_diff_.push_back(nannot);

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
                                      score_t extra_score,
                                      const std::function<void(Alignment&&)> &callback) {
    Alignment alignment = construct_alignment(ops, clipping, window, path, match,
                                              end_score, offset, extra_score);
    alignment.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();

    auto [base_labels, base_coords_mut]
        = annotation_buffer_.get_labels_and_coords(alignment.get_nodes().front(), !clipping);
    assert(base_labels);
    assert(base_labels->size());
    std::shared_ptr<const CoordinateSet> base_coords = base_coords_mut;

    if (!clipping)
        base_labels = &seed_->label_columns;

    auto call_alignment = [&]() {
        assert(alignment.label_columns.size());
        if (label_diff_.size() && label_diff_.back() == nannot) {
            label_diff_.pop_back();
            remaining_labels_i_ = annotation_buffer_.cache_column_set(std::move(label_diff_));
            assert(remaining_labels_i_ != nannot);
            label_diff_ = Columns{};
        }

        callback(std::move(alignment));
    };

    if (!annotation_buffer_.has_coordinates()) {
        alignment.label_columns = std::move(label_intersection_);
        call_alignment();
        return;
    }

    ssize_t dist = alignment.get_nodes().size() - 1;
    if (!clipping) {
        base_coords = std::shared_ptr<const CoordinateSet> {
            std::shared_ptr<const CoordinateSet>{}, &seed_->label_coordinates
        };
        dist -= seed_->get_offset();
        if (dynamic_cast<const RCDBG*>(graph_))
            dist = alignment.get_sequence().size() - seed_->get_sequence().size();
    }
    assert(base_coords);

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
        auto [cur_labels, cur_coords_mut]
            = annotation_buffer_.get_labels_and_coords(alignment.get_nodes().back(), false);
        assert(cur_labels);
        assert(cur_labels->size());
        assert(cur_coords_mut);
        assert(cur_coords_mut->size());
        std::shared_ptr<const CoordinateSet> cur_coords = cur_coords_mut;
        if (dynamic_cast<const RCDBG*>(graph_)) {
            std::swap(cur_labels, base_labels);
            std::swap(cur_coords, base_coords);
        }
        try {
            utils::match_indexed_values(
                base_labels->begin(), base_labels->end(), base_coords->begin(),
                cur_labels->begin(), cur_labels->end(), cur_coords->begin(),
                [&](Column c, const auto &coords, const auto &other_coords) {
                    while (label_it != label_end_it && c > *label_it) {
                        ++label_it;
                    }

                    if (label_it == label_end_it)
                        throw std::bad_function_call();

                    if (c < *label_it)
                        return;

                    Alignment::Tuple overlap;
                    utils::set_intersection(coords.begin(), coords.end(),
                                            other_coords.begin(), other_coords.end(),
                                            std::back_inserter(overlap),
                                            dist);
                    if (overlap.size()) {
                        alignment.label_columns.emplace_back(c);
                        alignment.label_coordinates.emplace_back(std::move(overlap));
                    }
                }
            );
        } catch (const std::bad_function_call&) {}
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
        annotation_buffer_(graph, annotator, config.row_batch_size,
                           config.row_batch_size) {
    // do not use a global xdrop cutoff since we need separate cutoffs for each label
    if (annotation_buffer_.has_coordinates()) {
        logger->trace("Coordinates detected. Enabling seed chaining");
        this->config_.global_xdrop = false;
        this->config_.chain_alignments = true;
    }

    this->config_.min_seed_length = std::min(graph.get_k(), this->config_.min_seed_length);
    this->config_.max_seed_length = std::min(graph.get_k(), this->config_.max_seed_length);
}

template <class Seeder, class Extender, class AlignmentCompare>
LabeledAligner<Seeder, Extender, AlignmentCompare>::~LabeledAligner() {
    logger->trace("Buffered {}/{} nodes and {} label combinations. Mem usage {} MB",
                  annotation_buffer_.num_nodes_buffered(),
                  this->graph_.num_nodes(),
                  annotation_buffer_.num_column_sets(),
                  get_curr_RSS() / 1e6);
}

template <class Seeder, class Extender, class AlignmentCompare>
auto LabeledAligner<Seeder, Extender, AlignmentCompare>
::build_seeders(const std::vector<IDBGAligner::Query> &seq_batch,
                const std::vector<AlignmentResults> &wrapped_seqs) const -> BatchSeeders {
    BatchSeeders seeders
        = DBGAligner<Seeder, Extender, AlignmentCompare>::build_seeders(seq_batch, wrapped_seqs);

    // now we're going to filter the seeds
    logger->trace("Filtering seeds by label. Cur mem usage {} MB", get_curr_RSS() / 1e6);
    std::vector<std::pair<std::vector<Seed>, size_t>> counted_seeds;
    std::vector<std::pair<std::vector<Seed>, size_t>> counted_seeds_rc;

    size_t num_seeds = 0;
    size_t num_seeds_rc = 0;

#if ! _PROTEIN_GRAPH
    std::vector<bool> has_rc;
    has_rc.reserve(seeders.size());
#endif

    for (auto &[seeder, seeder_rc] : seeders) {
        counted_seeds.emplace_back(seeder->get_seeds(), seeder->get_num_matches());
        seeder.reset();
        num_seeds += counted_seeds.back().first.size();

        auto add_seeds = [&](const auto &seeds) {
            for (const Seed &seed : seeds) {
                annotation_buffer_.queue_path(seed.get_nodes());
            }
        };

        add_seeds(counted_seeds.back().first);

#if ! _PROTEIN_GRAPH
        has_rc.emplace_back(seeder_rc);
        if (seeder_rc) {
            counted_seeds_rc.emplace_back(seeder_rc->get_seeds(),
                                          seeder_rc->get_num_matches());
            seeder_rc.reset();
            num_seeds_rc += counted_seeds_rc.back().first.size();
            add_seeds(counted_seeds_rc.back().first);
        }
#endif
    }

    logger->trace("Prefetching labels for {} seeds. Cur mem usage {} MB",
                  num_seeds + num_seeds_rc, get_curr_RSS() / 1e6);
    annotation_buffer_.fetch_queued_annotations();
    logger->trace("Done prefetching. Cur mem usage {} MB", get_curr_RSS() / 1e6);

    size_t num_seeds_left = 0;
    size_t num_seeds_rc_left = 0;

    for (size_t i = 0; i < counted_seeds.size(); ++i) {
        auto &[seeder, seeder_rc] = seeders[i];
        auto &[seeds, num_matching] = counted_seeds[i];
        if (seeds.size()) {
            num_matching = filter_seeds(seeds);
            num_seeds_left += seeds.size();
        }

        seeder = make_shared<ManualMatchingSeeder>(std::move(seeds), num_matching, this->config_);

#if ! _PROTEIN_GRAPH
        if (has_rc[i]) {
            auto &[seeds, num_matching] = counted_seeds_rc[i];
            if (seeds.size()) {
                num_matching = filter_seeds(seeds);
                num_seeds_rc_left += seeds.size();
            }

            seeder_rc = make_shared<ManualMatchingSeeder>(std::move(seeds), num_matching, this->config_);
        }
#endif
    }

    logger->trace("Old seed count: {}\tNew seed count: {}",
                  num_seeds + num_seeds_rc,
                  num_seeds_left + num_seeds_rc_left);

    return seeders;
}

template <class Seeder, class Extender, class AlignmentCompare>
size_t LabeledAligner<Seeder, Extender, AlignmentCompare>
::filter_seeds(std::vector<Seed> &seeds) const {
    if (seeds.empty())
        return 0;

    size_t query_size = seeds[0].get_clipping() + seeds[0].get_end_clipping()
                            + seeds[0].get_query_view().size();

    Columns labels;

    {
        VectorMap<Column, sdsl::bit_vector> label_mapper;
        for (size_t j = 0; j < seeds.size(); ++j) {
            const Seed &seed = seeds[j];
            assert(seed.get_nodes().size() == 1);
            size_t end = seed.get_clipping() + this->graph_.get_k() - seed.get_offset();
            node_index node = seed.get_nodes()[0];

            assert(annotation_buffer_.get_labels(node));

            for (uint64_t label : *annotation_buffer_.get_labels(node)) {
                auto &indicator = label_mapper[label];
                if (indicator.empty())
                    indicator = sdsl::bit_vector(query_size, false);

                for (size_t i = seed.get_clipping(); i < end; ++i) {
                    indicator[i] = true;
                }
            }
        }

        if (label_mapper.empty()) {
            seeds.clear();
            return 0;
        }

        std::vector<std::pair<Column, uint64_t>> label_counts;
        label_counts.reserve(label_mapper.size());
        for (const auto &[c, indicator] : label_mapper) {
            label_counts.emplace_back(c, sdsl::util::cnt_one_bits(indicator));
        }

        std::sort(label_counts.begin(), label_counts.end(), utils::GreaterSecond());

        double cutoff = this->config_.min_exact_match * query_size;
        auto it = std::find_if(label_counts.begin(), label_counts.end(),
                               [cutoff](const auto &a) { return a.second < cutoff; });

        DEBUG_LOG("Keeping {} / {} labels",
                  std::distance(label_counts.begin(), it), label_counts.size());

        label_counts.erase(it, label_counts.end());

        labels.reserve(label_counts.size());
        for (const auto &[label, count] : label_counts) {
            labels.push_back(label);
        }
    }

    if (labels.empty()) {
        seeds.clear();
        return 0;
    }

    std::sort(labels.begin(), labels.end());

    for (size_t j = 0; j < seeds.size(); ++j) {
        Seed &seed = seeds[j];
        const std::vector<node_index> &nodes = seed.get_nodes();
        assert(nodes.size() == 1);
        if (!seed.label_encoder) {
            seed.label_columns.clear();
            auto fetch_labels = annotation_buffer_.get_labels(nodes[0]);
            assert(fetch_labels);
            std::set_intersection(fetch_labels->begin(), fetch_labels->end(),
                                  labels.begin(), labels.end(),
                                  std::back_inserter(seed.label_columns));
            if (seed.label_columns.size())
                seed.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();
        }
    }

    auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [&](const auto &a) {
        return !a.label_encoder || a.label_columns.empty();
    });

    seeds.erase(seed_it, seeds.end());

    assert(std::all_of(seeds.begin(), seeds.end(), [&](const auto &a) {
        return a.get_query_view().size() >= this->config_.min_seed_length;
    }));

    return get_num_char_matches_in_seeds(seeds.begin(), seeds.end());
}

template class LabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
