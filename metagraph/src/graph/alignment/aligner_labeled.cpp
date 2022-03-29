#include "aligner_labeled.hpp"

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

void LabeledExtender::flush() {
    annotation_buffer_.fetch_queued_annotations();
    for ( ; last_flushed_table_i_ < table.size(); ++last_flushed_table_i_) {
        auto &table_elem = table[last_flushed_table_i_];

        size_t parent_i = table_elem.parent_i;
        assert(parent_i < last_flushed_table_i_);

        auto clear = [&]() {
            DEBUG_LOG("Removed table element {}", last_flushed_table_i_);
            node_labels_[last_flushed_table_i_] = 0;
            std::fill(table_elem.S.begin(), table_elem.S.end(), config_.ninf);
            std::fill(table_elem.E.begin(), table_elem.E.end(), config_.ninf);
            std::fill(table_elem.F.begin(), table_elem.F.end(), config_.ninf);
        };

        if (!node_labels_[parent_i]) {
            clear();
            continue;
        }

        if (table_elem.node == DeBruijnGraph::npos)
            continue;

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
            node_labels_[last_flushed_table_i_]
                = annotation_buffer_.cache_column_set(std::move(intersect_labels));
        }
    }
}

bool LabeledExtender::set_seed(const Alignment &seed) {
    if (!DefaultColumnExtender::set_seed(seed))
        return false;

    assert(std::all_of(seed.get_nodes().begin(), seed.get_nodes().end(),
                    [&](node_index n) {
                        return n == DeBruijnGraph::npos || annotation_buffer_.get_labels(n);
                    }));

    // the first node of the seed has already been flushed
    last_flushed_table_i_ = 1;

    remaining_labels_i_ = annotation_buffer_.cache_column_set(seed.label_columns);
    assert(remaining_labels_i_ != nannot);
    node_labels_.assign(1, remaining_labels_i_);
    base_coords_ = seed.label_coordinates;
    if (!base_coords_.size())
        return true;

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

    // if we are in the seed and want to force the seed to be fixed, automatically
    // take the next node in the seed
    size_t next_offset = table[table_i].offset + 1;
    bool in_seed = next_offset - seed_->get_offset() < seed_->get_sequence().size()
                    && (next_offset < graph_->get_k() || force_fixed_seed);

    std::vector<std::tuple<node_index, char, score_t>> outgoing;
    DefaultColumnExtender::call_outgoing(node, max_prefetch_distance,
        [&](node_index next, char c, score_t score) {
            outgoing.emplace_back(next, c, score);
            if (!in_seed)
                annotation_buffer_.queue_path({ next });
        },
        table_i, force_fixed_seed
    );

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
    const auto &columns = annotation_buffer_.get_cached_column_set(node_labels_[table_i]);

    // no coordinates are present in the annotation
    if (!annotation_buffer_.get_labels_and_coords(node).second) {
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
        auto [next_labels, next_coords]
            = annotation_buffer_.get_labels_and_coords(next, false);

        assert(next_coords);

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

    auto [base_labels, base_coords]
        = annotation_buffer_.get_labels_and_coords(alignment.get_nodes().front(), clipping);
    assert(base_labels);
    assert(base_labels->size());

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

    if (!base_coords) {
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
            = annotation_buffer_.get_labels_and_coords(alignment.get_nodes().back(), false);
        assert(cur_labels);
        assert(cur_labels->size());
        assert(cur_coords);
        assert(cur_coords->size());
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
                        throw std::exception();

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
        annotation_buffer_(graph, annotator, config.row_batch_size,
                           config.max_num_seeds_per_locus),
        max_seed_length_(config.max_seed_length) {
    // do not use a global xdrop cutoff since we need separate cutoffs for each label
    if (annotation_buffer_.has_coordinates())
        this->config_.global_xdrop = false;

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
                annotation_buffer_.queue_path(std::vector<node_index>(seed.get_nodes()));
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

inline size_t get_num_matches(const std::vector<Seed> &seeds) {
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

template <class AIt, class BIt, class CIt, class OutIt, class OutIt2>
void matched_intersection(AIt a_begin, AIt a_end, BIt a_c_begin,
                          CIt b_begin, CIt b_end,
                          OutIt out1, OutIt2 out2) {
    while (a_begin != a_end && b_begin != b_end) {
        if (*a_begin < *b_begin) {
            ++a_begin;
            ++a_c_begin;
        } else if (*a_begin > *b_begin) {
            ++b_begin;
        } else {
            *out1 = *a_begin;
            *out2 = *a_c_begin;
            ++out1;
            ++out2;
            ++a_begin;
            ++a_c_begin;
            ++b_begin;
        }
    }
}

template <class Seeder, class Extender, class AlignmentCompare>
size_t LabeledAligner<Seeder, Extender, AlignmentCompare>
::filter_seeds(std::vector<Seed> &seeds) const {
    if (seeds.empty())
        return 0;

    size_t query_size = seeds[0].get_clipping() + seeds[0].get_end_clipping()
                            + seeds[0].get_query_view().size();

    Columns labels;
    tsl::hopscotch_map<node_index, std::vector<size_t>> node_to_seeds;

    {
        VectorMap<Column, sdsl::bit_vector> label_mapper;
        for (size_t j = 0; j < seeds.size(); ++j) {
            const Seed &seed = seeds[j];
            assert(seed.get_nodes().size() == 1);
            size_t end = seed.get_clipping() + this->graph_.get_k() - seed.get_offset();
            node_index node = seed.get_nodes()[0];

            if (max_seed_length_ > this->config_.max_seed_length)
                node_to_seeds[node].emplace_back(j);

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

        if (label_mapper.empty())
            return get_num_matches(seeds);

        std::vector<std::pair<Column, uint64_t>> label_counts;
        label_counts.reserve(label_mapper.size());
        for (const auto &[c, indicator] : label_mapper) {
            label_counts.emplace_back(c, sdsl::util::cnt_one_bits(indicator));
        }

        std::sort(label_counts.begin(), label_counts.end(), utils::GreaterSecond());

        double cutoff = this->config_.min_exact_match * query_size;
        auto it = std::find_if(label_counts.begin(), label_counts.end(),
                               [cutoff](const auto &a) { return a.second < cutoff; });

        logger->trace("Keeping {} / {} labels",
                      std::distance(label_counts.begin(), it), label_counts.size());

        label_counts.erase(it, label_counts.end());

        labels.reserve(label_counts.size());
        for (const auto &[label, count] : label_counts) {
            labels.push_back(label);
        }
    }

    if (labels.empty()) {
        seeds.clear();
        return get_num_matches(seeds);
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

    if (max_seed_length_ > this->config_.max_seed_length) {
        // merge seeds
        for (size_t j = 0; j < seeds.size(); ++j) {
            Seed &seed = seeds[j];
            if (annotation_buffer_.has_coordinates() && !(j % 100)) {
                std::vector<node_index> nodes;
                for (size_t k = j; k < std::min(j + 100, seeds.size()); ++k) {
                    if (seeds[k].size())
                        nodes.emplace_back(seeds[k].get_nodes()[0]);
                }
                annotation_buffer_.prefetch_coords(nodes);
            }

            if (seed.empty() || !seed.get_end_clipping() || !seed.label_encoder)
                continue;

            node_index next_node = this->graph_.traverse(
                seed.get_nodes().back(),
                seed.get_query_view().data()[seed.get_query_view().size()]
            );

            if (!next_node)
                continue;

            auto find = node_to_seeds.find(next_node);
            if (find == node_to_seeds.end())
                continue;

            for (size_t j : find->second) {
                auto &next_seed = seeds[j];
                if (next_seed.empty() || !next_seed.label_encoder
                        || next_seed.get_end_clipping() >= seed.get_end_clipping()
                        || next_seed.get_offset() != seed.get_offset()
                        || next_seed.label_columns != seed.label_columns) {
                    continue;
                }
                assert(next_node == next_seed.get_nodes()[0]);

                if (annotation_buffer_.has_coordinates()) {
                    if (seed.label_coordinates.size() != seed.label_columns.size()) {
                        auto [fetch_labels, fetch_coords]
                            = annotation_buffer_.get_labels_and_coords(seed.get_nodes()[0], false);
                        Columns dummy;
                        matched_intersection(fetch_labels->begin(), fetch_labels->end(),
                                             fetch_coords->begin(),
                                             seed.label_columns.begin(),
                                             seed.label_columns.end(),
                                             std::back_inserter(dummy),
                                             std::back_inserter(seed.label_coordinates));
                        assert(dummy == seed.label_columns);
                    }

                    if (next_seed.label_coordinates.size() != next_seed.label_columns.size()) {
                        auto [fetch_labels, fetch_coords]
                            = annotation_buffer_.get_labels_and_coords(next_seed.get_nodes()[0], false);
                        Columns dummy;
                        matched_intersection(fetch_labels->begin(), fetch_labels->end(),
                                             fetch_coords->begin(),
                                             next_seed.label_columns.begin(),
                                             next_seed.label_columns.end(),
                                             std::back_inserter(dummy),
                                             std::back_inserter(next_seed.label_coordinates));
                        assert(dummy == next_seed.label_columns);
                    }
                }

                size_t k = 0;
                for ( ; k < seed.label_coordinates.size(); ++k) {
                    if (!overlap_with_diff(seed.label_coordinates[k],
                                           next_seed.label_coordinates[k],
                                           seed.get_nodes().size())) {
                        break;
                    }
                }

                if (k != seed.label_coordinates.size())
                    continue;

                // they can be merged
                seed.expand(next_seed.get_nodes());
                next_seed = Seed();
                assert(next_seed.empty());
            }
        }
    }

    for (Seed &seed : seeds) {
        if (seed.label_coordinates.size()) {
            size_t removed_labels = 0;
            for (size_t i = 0; i < seed.label_coordinates.size(); ++i) {
                if (seed.label_coordinates[i].size() > this->config_.max_num_seeds_per_locus) {
                    seed.label_coordinates[i] = Alignment::Tuple();
                    ++removed_labels;
                }
            }

            if (removed_labels == seed.label_coordinates.size())
                seed.label_encoder = nullptr;
        } else if (annotation_buffer_.has_coordinates() && seed.label_columns.size()) {
            seed.label_encoder = nullptr;
        }
    }

    auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [&](const auto &a) {
        return !a.label_encoder;
    });

    seeds.erase(seed_it, seeds.end());

    for (Seed &seed : seeds) {
        if (seed.get_offset() && seed.label_coordinates.size()) {
            for (auto &tuple : seed.label_coordinates) {
                for (auto &coord : tuple) {
                    coord += seed.get_offset();
                }
            }
        }
    }

    assert(std::all_of(seeds.begin(), seeds.end(), [&](const auto &a) {
        return a.get_query_view().size() >= this->config_.min_seed_length;
    }));

    return get_num_matches(seeds);
}

template class LabeledAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
