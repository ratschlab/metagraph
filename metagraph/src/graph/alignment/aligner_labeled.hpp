#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <optional>

#include <tsl/hopscotch_map.h>
#include <tsl/ordered_set.h>

#include "dbg_aligner.hpp"
#include "common/hashers/hash.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename Key, class Hash = std::hash<Key>, typename IndexType = uint64_t,
          class EqualTo = std::equal_to<Key>, class Allocator = std::allocator<Key>,
          class Container = std::vector<Key, Allocator>>
using VectorSet = tsl::ordered_set<Key, Hash, EqualTo, Allocator, Container, IndexType>;


class AnnotationBuffer {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::binmat::BinaryMatrix::Row Row;
    typedef annot::matrix::MultiIntMatrix::Tuple Tuple;

    typedef std::reference_wrapper<const Alignment::LabelSet> LabelSet;
    typedef std::reference_wrapper<const Alignment::CoordinateSet> CoordsSet;

    static constexpr Row nrow = std::numeric_limits<Row>::max();

    AnnotationBuffer(const AnnotatedDBG &anno_graph);

    const AnnotatedDBG& get_anno_graph() const { return anno_graph_; }
    const annot::matrix::MultiIntMatrix* get_coordinate_matrix() const { return multi_int_; }

    // flush the buffer and fetch their annotations from the AnnotatedDBG
    void flush();

    // push (a) node(s) to the buffer
    node_index add_node(node_index node);

    std::pair<std::vector<node_index>, bool>
    add_path(const std::vector<node_index> &path, std::string sequence);

    // get the annotations and coordinates of a node if they have been fetched
    std::pair<std::optional<LabelSet>, std::optional<CoordsSet>>
    get_labels_and_coordinates(node_index node) const {
        std::pair<std::optional<LabelSet>, std::optional<CoordsSet>> ret_val {
            std::nullopt, std::nullopt
        };

        auto it = labels_.find(node);

        // if the node hasn't been seen before, or if its annotations haven't
        // been flushed, return nothing
        if (it == labels_.end() || it->second.second == nannot)
            return ret_val;

        ret_val.first = std::cref(labels_set_.data()[it->second.second]);

        // if no coordinates are present, return just the labels
        if (!multi_int_)
            return ret_val;

        assert(static_cast<size_t>(it - labels_.begin()) < label_coords_.size());
        ret_val.second = std::cref(label_coords_[it - labels_.begin()]);
        return ret_val;
    }

    // get the annotations of a node if they have been fetched
    inline std::optional<LabelSet> get_labels(node_index node) const {
        return get_labels_and_coordinates(node).first;
    }

    size_t num_cached() const { return added_rows_.size(); }

  private:
    const AnnotatedDBG &anno_graph_;
    const annot::matrix::MultiIntMatrix *multi_int_;

    // placeholder index for an unfetched annotation
    static constexpr size_t nannot = std::numeric_limits<size_t>::max();

    // keep a unique set of annotation rows
    VectorSet<Vector<Column>, utils::VectorHash> labels_set_;

    // map nodes to indexes in labels_set_
    VectorMap<node_index, std::pair<Row, size_t>> labels_;

    // map each element in labels_ to a set of coordinates
    std::vector<Vector<Tuple>> label_coords_;

    // buffer of accessed nodes and their corresponding annotation rows
    std::vector<Row> added_rows_;
    std::vector<node_index> added_nodes_;
};

struct CoordIntersection {
    CoordIntersection(size_t offset = 0) : offset_(offset) {}

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

    size_t offset_;
};

class LabeledExtender : public DefaultColumnExtender {
  public:
    typedef AnnotationBuffer::Column Column;
    typedef AnnotationBuffer::Tuple Tuple;
    typedef AlignmentAggregator<LocalAlignmentLess> Aggregator;

    LabeledExtender(AnnotationBuffer &labeled_graph,
                    const DBGAlignerConfig &config,
                    const Aggregator &,
                    std::string_view query)
          : DefaultColumnExtender(labeled_graph.get_anno_graph().get_graph(), config, query),
            labeled_graph_(labeled_graph) {}

    virtual ~LabeledExtender() {}

  protected:
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window) override {
        // extract all labels for explored nodes
        labeled_graph_.flush();

        // run backtracking
        return DefaultColumnExtender::backtrack(min_path_score, window);
    }

    virtual std::vector<Alignment> extend(score_t min_path_score,
                                          bool force_fixed_seed) override {
        last_buffered_table_i_ = 0;
        return DefaultColumnExtender::extend(min_path_score, force_fixed_seed);
    }

    AnnotationBuffer &labeled_graph_;
    size_t last_buffered_table_i_;
};

class LabeledBacktrackingExtender : public LabeledExtender {
  public:
    LabeledBacktrackingExtender(AnnotationBuffer &labeled_graph,
                                const DBGAlignerConfig &config,
                                const Aggregator &aggregator,
                                std::string_view query)
          : LabeledExtender(labeled_graph, config, aggregator, query),
            extensions_(labeled_graph_.get_anno_graph().get_graph(),
                        aggregator.get_query(false),
                        aggregator.get_query(true), this->config_) {}

    virtual ~LabeledBacktrackingExtender() {}

  protected:
    virtual bool set_seed(const Alignment &seed) override {
        if (DefaultColumnExtender::set_seed(seed)) {
            seed_labels_ = seed.label_columns;
            seed_label_coordinates_ = seed.label_coordinates;
            return true;
        }

        return false;
    }

    virtual std::vector<Alignment> extend(score_t min_path_score,
                                          bool force_fixed_seed) override final {
        // the overridden backtrack populates extensions_, so this should return nothing
        LabeledExtender::extend(min_path_score, force_fixed_seed);

        // fetch the alignments from extensions_
        return extensions_.get_alignments();
    }

    // backtrack through the DP table to reconstruct alignments
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window) override final {
        // reset the per-node temporary label storage
        diff_label_sets_.clear();

        // run backtracking
        return LabeledExtender::backtrack(min_path_score, window);
    }

    // overrides for backtracking helpers
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &) const override final {
        return this->seed_->label_columns.size() && seed_labels_.empty();
    }

    virtual bool terminate_backtrack() const override final { return label_intersection_.empty(); }
    virtual bool skip_backtrack_start(size_t i) override final;

    // since multi-node seeds may span across different labels, we no longer
    // want the restriction that the seed must be a prefix of the extended alignment
    virtual bool fixed_seed() const override final { return false; }

    // this override ensures that outgoing nodes are label- and coordinate-consistent
    // (when applicable)
    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char, score_t)> &callback,
                               size_t table_i,
                               bool force_fixed_seed = false) override final;

    // this method calls multiple label-consistent alignments by backtracking
    virtual void call_alignments(score_t cur_cell_score,
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
                                 const std::function<void(Alignment&&)> &callback) override final;

    virtual void pop(size_t i) override { DefaultColumnExtender::pop(i); }

  private:
    // local set of alignments
    Aggregator extensions_;

    // keep track of the label set for the current backtracking
    Vector<Column> label_intersection_;
    score_t cur_min_path_score_;
    size_t last_path_size_;

    Vector<Tuple> label_intersection_coords_;

    Vector<Column> seed_labels_;
    Vector<Tuple> seed_label_coordinates_;

    // After a node has been visited during backtracking, we keep track of which
    // of its labels haven't been considered yet. This way, if backtracking is
    // called from this node, then we can restrict it to these labels.
    tsl::hopscotch_map<size_t, std::pair<Vector<Column>, Vector<Tuple>>> diff_label_sets_;
};

template <class Extender = LabeledBacktrackingExtender,
          class Seeder = UniMEMSeeder,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    typedef Alignment::score_t score_t;
    typedef Alignment::node_index node_index;
    typedef Alignment::Column Column;

    LabeledAligner(const AnnotatedDBG &anno_graph,
                   const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config)
          : ISeedAndExtendAligner<AlignmentCompare>(graph, config),
            labeled_graph_(anno_graph) {
        if (labeled_graph_.get_coordinate_matrix()
                && std::is_same_v<Extender, LabeledBacktrackingExtender>) {
            // do not use a global xdrop cutoff since we need separate cutoffs
            // for each label
            this->config_.global_xdrop = false;
        }
    }

  protected:
    typedef typename ISeedAndExtendAligner<AlignmentCompare>::BatchSeeders BatchSeeders;
    mutable AnnotationBuffer labeled_graph_;

    std::shared_ptr<IExtender>
    build_extender(std::string_view query,
                   const typename Extender::Aggregator &aggregator,
                   const DBGAlignerConfig &config) const override final {
        return std::make_shared<Extender>(labeled_graph_, config, aggregator, query);
    }

    std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<IDBGAligner::node_index> &nodes) const override final {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }

    void filter_seeds(BatchSeeders &seeders) const override final {
        common::logger->trace("Filtering seeds by label");
        std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds;
        std::vector<std::pair<std::vector<Alignment>, size_t>> counted_seeds_rc;

        size_t num_seeds = 0;
        size_t num_seeds_rc = 0;

        for (auto &[seeder, nodes, seeder_rc, nodes_rc] : seeders) {
            counted_seeds.emplace_back(seeder->get_seeds(), seeder->get_num_matches());
            num_seeds += counted_seeds.back().first.size();
            for (const Alignment &seed : counted_seeds.back().first) {
                labeled_graph_.add_path(
                    seed.get_nodes(),
                    std::string(seed.get_nodes().size() + this->graph_.get_k() - 1, '#')
                );
            }

#if ! _PROTEIN_GRAPH
            if (seeder_rc) {
                counted_seeds_rc.emplace_back(seeder_rc->get_seeds(), seeder_rc->get_num_matches());
                num_seeds_rc += counted_seeds_rc.back().first.size();
                for (const Alignment &seed : counted_seeds_rc.back().first) {
                    labeled_graph_.add_path(
                        seed.get_nodes(),
                        std::string(seed.get_nodes().size() + this->graph_.get_k() - 1, '#')
                    );
                }
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

  private:
    static inline size_t get_num_matches(const std::vector<Alignment> &seeds) {
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

    // find the most frequent labels among the seeds and restrict graph traversal
    // to those labeled paths during extension
    size_t filter_seeds(std::vector<Alignment> &seeds) const {
        VectorMap<Column, uint64_t> label_counter;
        for (const Alignment &seed : seeds) {
            for (node_index node : seed.get_nodes()) {
                if (auto labels = labeled_graph_.get_labels(node)) {
                    for (uint64_t label : labels->get()) {
                        ++label_counter[label];
                    }
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

            seed.label_encoder = &labeled_graph_.get_anno_graph().get_annotator().get_label_encoder();

            for (size_t i = 1; i < nodes.size() && seed.label_columns.size(); ++i) {
                auto [next_fetch_labels, next_fetch_coords]
                    = labeled_graph_.get_labels_and_coordinates(nodes[i]);
                if (next_fetch_coords) {
                    Vector<Column> label_inter;
                    Alignment::CoordinateSet coord_inter;
                    utils::indexed_set_op<Alignment::Tuple, CoordIntersection>(
                        seed.label_columns.begin(),
                        seed.label_columns.end(),
                        seed.label_coordinates.begin(),
                        next_fetch_labels->get().begin(),
                        next_fetch_labels->get().end(),
                        next_fetch_coords->get().begin(),
                        std::back_inserter(label_inter),
                        std::back_inserter(coord_inter),
                        i
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

            if (seed.get_offset() && seed.label_coordinates.size()) {
                for (auto &tuple : seed.label_coordinates) {
                    for (auto &coord : tuple) {
                        coord += seed.get_offset();
                    }
                }
            }
        }

        auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [](const auto &a) {
            return a.label_columns.empty();
        });

        seeds.erase(seed_it, seeds.end());

        return get_num_matches(seeds);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
