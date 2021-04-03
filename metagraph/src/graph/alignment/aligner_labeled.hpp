#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_set.hpp"
#include "common/utils/template_utils.hpp"
#include "common/hashers/hash.hpp"
#include "common/vectors/bitmap.hpp"

namespace mtg {
namespace graph {

class AnnotatedDBG;

namespace align {


class ILabeledDBGAligner : public ISeedAndExtendAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef Vector<uint64_t> Targets;

    ILabeledDBGAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config);

    virtual ~ILabeledDBGAligner() {}

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    typedef std::pair<std::vector<node_index> /* forward */,
                      std::vector<node_index> /* reverse complement */ > Mapping;
    typedef std::pair<sdsl::bit_vector /* forward */,
                      sdsl::bit_vector /* reverse complement */ > Signature;
    typedef std::vector<Mapping> BatchMapping;
    typedef std::vector<std::pair<Targets, Signature>> QueryLabels;
    typedef std::vector<QueryLabels> BatchLabels;

    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;

    std::pair<BatchMapping, BatchLabels>
    map_and_label_query_batch(const QueryGenerator &generate_query) const;
};

bool check_targets(const AnnotatedDBG &anno_graph,
                   const Alignment<DeBruijnGraph::node_index> &path);

class LabeledSeedFilter : public ISeedFilter {
  public:
    LabeledSeedFilter(size_t k) : k_(k) {}
    Vector<uint64_t> labels_to_keep(const DBGAlignment &seed);
    void update_seed_filter(const LabeledNodeRangeGenerator &generator);

  private:
    size_t k_;
    tsl::hopscotch_map<node_index, tsl::hopscotch_map<uint64_t, std::pair<size_t, size_t>>> visited_nodes_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class LabeledColumnExtender;

template <class BaseSeeder = ExactSeeder<>,
          class Extender = LabeledColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledDBGAligner : public ILabeledDBGAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef IDBGAligner::AlignmentCallback AlignmentCallback;

    template <typename... Args>
    LabeledDBGAligner(Args&&... args) : ILabeledDBGAligner(std::forward<Args>(args)...) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final;

  protected:
    typedef LabeledSeeder<BaseSeeder> Seeder;

    std::shared_ptr<ISeeder<node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 std::vector<node_index>&& base_nodes,
                 const std::vector<Targets> &targets,
                 const std::vector<std::unique_ptr<bitmap>> &signatures) const;
};

template <typename NodeType>
class LabeledColumnExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef ILabeledDBGAligner::Targets Targets;

    struct TargetColumnsSet {
        typedef std::vector<std::pair<Targets, size_t>> Storage;
        typedef Storage::iterator iterator;
        typedef Storage::const_iterator const_iterator;

        Storage targets_;
        tsl::hopscotch_map<Targets, size_t, utils::VectorHash> map_;
        std::vector<size_t> empty_slots_;

        template <typename... Args>
        std::pair<iterator, bool> emplace(Args&&... args) {
            iterator ret;
            auto [it, inserted] = map_.emplace(
                std::forward<Args>(args)...,
                empty_slots_.size() ? empty_slots_.back() : targets_.size()
            );

            assert(it->second <= targets_.size());
            if (it->second == targets_.size()) {
                assert(inserted);
                targets_.emplace_back(it->first, 0);
                ret = targets_.end() - 1;
            } else {
                assert(it->second < targets_.size());
                ret = targets_.begin() + it->second;

                if (inserted) {
                    assert(empty_slots_.size());
                    assert(it->second == empty_slots_.back());
                    assert(ret->first.empty());

                    ret->first = it->first;
                    empty_slots_.pop_back();
                }
            }

            ++ret->second;

            return std::make_pair(ret, inserted);
        }

        iterator begin() { return targets_.begin(); }
        iterator end() { return targets_.end(); }

        size_t size() const {
            assert(targets_.size() == map_.size() + empty_slots_.size());
            return targets_.size();
        }

        std::pair<iterator, bool> merge(iterator old_it, const Targets &update) {
            assert(old_it != end());

            const Targets &old_targets = old_it->first;
            assert(old_it->second);
            if (std::includes(old_targets.begin(), old_targets.end(),
                              update.begin(), update.end())) {
                return std::make_pair(old_it, false);
            }

            Targets target_union;
            std::set_union(old_targets.begin(), old_targets.end(),
                           update.begin(), update.end(),
                           std::back_inserter(target_union));

            --old_it->second;
            if (!old_it->second) {
                map_.erase(old_it->first);
                auto [it, inserted] = map_.emplace(std::move(target_union),
                                                   old_it - targets_.begin());
                if (inserted) {
                    old_it->first = it->first;
                } else {
                    old_it->first = Targets{};
                    old_it->second = 0;
                    empty_slots_.push_back(old_it - targets_.begin());
                    old_it = targets_.begin() + it->second;
                }
                ++old_it->second;
            } else {
                old_it = emplace(std::move(target_union)).first;
            }

            return std::make_pair(old_it, true);
        }

        void increment(size_t i) {
            assert(i < targets_.size());
            ++targets_[i].second;
        }

        void decrement(size_t i) {
            assert(i < targets_.size());
            --targets_[i].second;
            if (!targets_[i].second) {
                map_.erase(targets_[i].first);
                targets_[i].first = Targets{};
                empty_slots_.push_back(i);
            }
        }
    };

    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::score_t score_t;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef tsl::hopscotch_map<NodeType, size_t> EdgeSetCache;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                          const DBGAlignerConfig &config,
                          std::string_view query,
                          LabeledSeedFilter &seed_filter,
                          TargetColumnsSet &target_columns,
                          EdgeSetCache &cached_edge_sets);

    virtual ~LabeledColumnExtender() {}

    virtual void initialize(const DBGAlignment &path) override;

  protected:
    typedef typename DefaultColumnExtender<NodeType>::Column Column;
    typedef typename DefaultColumnExtender<NodeType>::Scores Scores;
    typedef typename DefaultColumnExtender<NodeType>::AlignNode AlignNode;
    typedef typename DefaultColumnExtender<NodeType>::AlignNodeHash AlignNodeHash;
    typedef std::vector<AlignNode> Edges;

    virtual Edges get_outgoing(const AlignNode &node) const override;

    virtual void add_scores_to_column(Column &column,
                                      Scores&& scores,
                                      const AlignNode &node) override {
        const AlignNode &prev = std::get<6>(scores);
        assert(align_node_to_target_.count(prev));

        if (!align_node_to_target_.count(node)) {
            target_columns_->increment(align_node_to_target_[prev].first);
            align_node_to_target_.emplace(
                node,
                std::make_pair(align_node_to_target_[prev].first,
                               align_node_to_target_[prev].second
                                   ? align_node_to_target_[prev].first - 1
                                   : 0)
            );
        }

        DefaultColumnExtender<NodeType>::add_scores_to_column(
            column, std::move(scores), node
        );
    }

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &,
                                      const AlignNode &node) const override {
        assert(align_node_to_target_.count(node));
        size_t target_columns_idx = align_node_to_target_[node].first;

        auto find = backtrack_start_counter_.find(target_columns_idx);
        assert(find == backtrack_start_counter_.end()
            || find->second <= this->config_.num_alternative_paths);

        return find != backtrack_start_counter_.end()
            && find->second == this->config_.num_alternative_paths;
    }

    virtual void backtrack(score_t min_path_score,
                           AlignNode best_node,
                           tsl::hopscotch_set<AlignNode, AlignNodeHash> &prev_starts,
                           std::vector<DBGAlignment> &extensions) const override {
        size_t old_size = extensions.size();
        DefaultColumnExtender<NodeType>::backtrack(min_path_score, best_node,
                                                   prev_starts, extensions);
        assert(extensions.size() - old_size <= 1);
        if (extensions.size() > old_size) {
            if (extensions.back().get_offset()) {
                extensions.pop_back();
                return;
            }

            size_t target_columns_idx = align_node_to_target_.find(best_node)->second.first;
            assert(target_columns_idx < target_columns_->size());

            ++backtrack_start_counter_[target_columns_idx];
            extensions.back().target_columns = get_targets(target_columns_idx);

            assert(check_targets(anno_graph_, extensions.back()));
        }
    }

    virtual void pop(const AlignNode &node) override {
        if (align_node_to_target_.count(node)) {
            target_columns_->decrement(align_node_to_target_[node].first);
            align_node_to_target_.erase(node);
        }
    }

    const Targets& get_targets(size_t id) const {
        assert(id < target_columns_->size());
        const Targets &targets = (target_columns_->begin() + id)->first;
        assert(std::is_sorted(targets.begin(), targets.end()));
        return targets;
    }

    size_t get_target_id(const Targets &targets) const {
        assert(std::is_sorted(targets.begin(), targets.end()));
        auto find = target_columns_->emplace(targets).first;
        return find - target_columns_->begin();
    }

    size_t get_cached_target_id(NodeType node) const {
        assert(cached_edge_sets_->count(node));
        return (*cached_edge_sets_)[node];
    }

    bool update_target_cache(NodeType node, size_t a) const {
        if (cached_edge_sets_->count(node)) {
            size_t cached_id = get_cached_target_id(node);
            if (cached_id == a)
                return false;

            auto [it, inserted] = target_columns_->merge(
                target_columns_->begin() + cached_id,
                get_targets(a)
            );

            (*cached_edge_sets_)[node] = it - target_columns_->begin();
            return inserted;
        }

        cached_edge_sets_->emplace(node, a);
        return true;
    }

    size_t get_target_union(size_t a, size_t b) const {
        if (a == b)
            return a;

        const Targets &a_t = get_targets(a);
        const Targets &b_t = get_targets(b);

        Targets new_targets;
        std::set_union(a_t.begin(), a_t.end(), b_t.begin(), b_t.end(),
                       std::back_inserter(new_targets));

        return get_target_id(new_targets);
    }

    size_t get_target_intersection(size_t a, size_t b) const {
        if (a == b)
            return a;

        const Targets &a_t = get_targets(a);
        const Targets &b_t = get_targets(b);

        Targets new_targets;
        std::set_intersection(a_t.begin(), a_t.end(), b_t.begin(), b_t.end(),
                              std::back_inserter(new_targets));

        return get_target_id(new_targets);
    }

  private:
    const AnnotatedDBG &anno_graph_;
    std::shared_ptr<LabeledSeedFilter> seed_filter_;
    std::shared_ptr<TargetColumnsSet> target_columns_;
    std::shared_ptr<EdgeSetCache> cached_edge_sets_;
    mutable tsl::hopscotch_map<AlignNode, std::pair<size_t, size_t>, AlignNodeHash> align_node_to_target_;
    mutable tsl::hopscotch_map<size_t, size_t> backtrack_start_counter_;

    AlignNode get_next_align_node(NodeType node, char c, size_t dist_from_origin) const;
};


template <class BaseSeeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    typedef SeedAndExtendAlignerCore<AlignmentCompare> AlignerCore;
    typedef typename Extender::TargetColumnsSet TargetColumnsSet;
    typedef typename Extender::EdgeSetCache EdgeSetCache;
    auto mapped_batch = map_and_label_query_batch(generate_query);

    TargetColumnsSet target_columns_set;
    target_columns_set.emplace(Targets{});

    EdgeSetCache cached_edge_sets;

    size_t i = 0;
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        auto &[query_nodes_pair, target_columns] = mapped_batch;
        assert(config_.num_alternative_paths);

        LabeledSeedFilter seed_filter(this->graph_.get_k());
        AlignerCore aligner_core(graph_, config_, seed_filter, query, is_reverse_complement);
        DBGQueryAlignment &paths = aligner_core.get_paths();

        std::string_view this_query = paths.get_query(is_reverse_complement);
        std::string_view reverse = paths.get_query(true);
        assert(this_query == query);

        Extender extender(anno_graph_, config_, this_query, seed_filter,
                          target_columns_set, cached_edge_sets);
        Extender extender_rc(anno_graph_, config_, reverse, seed_filter,
                             target_columns_set, cached_edge_sets);

        std::vector<Targets> targets;
        std::vector<std::unique_ptr<bitmap>> signatures;
        std::vector<std::unique_ptr<bitmap>> signatures_rc;

        for (auto&& [target_columns, signature_pair] : target_columns[i]) {
            targets.emplace_back(std::move(target_columns));
            signatures.emplace_back(std::make_unique<bitmap_vector>(
                std::move(signature_pair.first)
            ));
            signatures_rc.emplace_back(std::make_unique<bitmap_vector>(
                std::move(signature_pair.second)
            ));
        }

        auto seeder = build_seeder(this_query, // use this_query since paths stores a copy
                                   is_reverse_complement,
                                   std::move(query_nodes_pair[i].first),
                                   targets,
                                   signatures);
        std::shared_ptr<ISeeder<node_index>> seeder_rc;
        if (config_.forward_and_reverse_complement
                || graph_.get_mode() == DeBruijnGraph::CANONICAL) {
            seeder_rc = build_seeder(reverse,
                                     !is_reverse_complement,
                                     std::move(query_nodes_pair[i].second),
                                     targets,
                                     signatures_rc);
        }

        if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
            assert(!is_reverse_complement);

            auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                     const auto &callback) {
                callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
            };

            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            aligner_core.align_both_directions(*seeder, *seeder_rc,
                                               extender, extender_rc,
                                               build_rev_comp_alignment_core);

        } else if (config_.forward_and_reverse_complement) {
            assert(!is_reverse_complement);
            aligner_core.align_best_direction(*seeder, *seeder_rc, extender, extender_rc);

        } else {
            aligner_core.align_one_direction(is_reverse_complement, *seeder, extender);
        }

        auto &aggregator = aligner_core.get_aggregator().data();
        aggregator.erase(std::numeric_limits<uint64_t>::max());
        if (aggregator.size() > config_.num_top_labels) {
            std::vector<std::pair<uint64_t, score_t>> scored_labels;
            scored_labels.reserve(aggregator.size());
            for (const auto &[target, path_queue] : aggregator) {
                scored_labels.emplace_back(target, path_queue.maximum()->get_score());
            }

            std::sort(scored_labels.begin(), scored_labels.end(), utils::GreaterSecond());
            auto start = scored_labels.begin() + config_.num_top_labels - 1;
            auto it = std::find_if(start + 1, scored_labels.end(), [&](const auto &a) {
                return a.second < start->second;
            });
            for ( ; it != scored_labels.end(); ++it) {
                aggregator.erase(it->first);
            }
        }

        aligner_core.flush();
        callback(header, std::move(paths));
        ++i;
    });
}

template <class BaseSeeder, class Extender, class AlignmentCompare>
inline auto LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::build_seeder(std::string_view query,
               bool is_reverse_complement,
               std::vector<node_index>&& base_nodes,
               const std::vector<Targets> &targets,
               const std::vector<std::unique_ptr<bitmap>> &signatures) const
        -> std::shared_ptr<ISeeder<node_index>> {
    if (this->config_.min_seed_length < this->graph_.get_k()
            && SuffixSeeder<Seeder>::get_base_dbg_succ(this->graph_)) {
        return std::make_shared<SuffixSeeder<Seeder>>(
            targets, signatures, graph_, query, is_reverse_complement,
            std::move(base_nodes), config_
        );
    }

    return std::make_shared<Seeder>(targets, signatures, graph_, query,
                                    is_reverse_complement, std::move(base_nodes), config_);
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
