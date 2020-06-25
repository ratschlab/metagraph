#ifndef __MASKED_DBG_ALIGNER_HPP__
#define __MASKED_DBG_ALIGNER_HPP__

#include "dbg_aligner.hpp"

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "common/vectors/vector_algorithm.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/annotated_graph_algorithm.hpp"


template <class BaseSeeder = ExactSeeder<>>
class LabeledSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::Seed Seed;
    typedef typename Seed::node_index node_index;

    LabeledSeeder(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : BaseSeeder(anno_graph.get_graph(), config), anno_graph_(anno_graph) {}

    void call_seeds(std::function<void(Seed&&)> callback) const;

  private:
    const AnnotatedDBG &anno_graph_;
};


template <class BaseExtender = DefaultColumnExtender<>>
class MaskedColumnExtender : public BaseExtender {
  public:
    typedef BaseExtender base_extender;
    typedef typename BaseExtender::DBGAlignment DBGAlignment;
    typedef typename BaseExtender::node_index node_index;
    typedef typename BaseExtender::score_t score_t;

    MaskedColumnExtender(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : BaseExtender(anno_graph.get_graph(), config),
            matrix_(anno_graph.get_annotation().get_matrix()),
            label_encoder_(anno_graph.get_annotation().get_label_encoder()) {}

    bool set_labels(const std::vector<std::string> &labels);
    bool set_labels(std::vector<uint64_t>&& label_codes);

    virtual
    void operator()(std::function<void(DBGAlignment&&, node_index)> callback,
                    score_t min_path_score = std::numeric_limits<score_t>::min()) override;

  protected:
    virtual std::pair<typename DPTable<node_index>::iterator, bool>
    emplace_node(node_index node, node_index incoming_node, char c, size_t size,
                 size_t best_pos = 0, size_t last_priority_pos = 0) override;

    virtual std::deque<std::pair<node_index, char>>
    fork_extension(node_index node,
                   std::function<void(DBGAlignment&&, node_index)> callback,
                   score_t min_path_score) override;

  private:
    const BinaryMatrix &matrix_;
    const annotate::LabelEncoder<std::string> &label_encoder_;
    tsl::hopscotch_set<uint64_t> label_codes_;

    typedef tsl::hopscotch_map<DeBruijnGraph::node_index,
                               DeBruijnGraph::node_index> CanonicalMap;

    CanonicalMap canonical_map_;

    bool graph_is_canonical() const { return this->get_graph().is_canonical_mode(); }

    // return a node's corresponding labeled node in anno_graph
    node_index get_labeled_node(node_index node);
};


template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>,
          template <class BaseSeeder> class LabeledSeeder = LabeledSeeder,
          template <class BaseExtender> class MaskedExtender = MaskedColumnExtender>
class MaskedDBGAligner : public SeedAndExtendAligner<LabeledSeeder<Seeder>,
                                                     MaskedExtender<Extender>> {
  public:
    typedef SeedAndExtendAligner<LabeledSeeder<Seeder>, MaskedExtender<Extender>> Aligner;
    typedef typename Aligner::DBGAlignment DBGAlignment;
    typedef typename Aligner::AlignmentGenerator AlignmentGenerator;

    MaskedDBGAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : anno_graph_(anno_graph), config_(config) {}

    const DeBruijnGraph& get_graph() const { return anno_graph_.get_graph(); }
    const DBGAlignerConfig& get_config() const { return config_; }

  private:
    LabeledSeeder<Seeder> build_seeder() const {
        return LabeledSeeder<Seeder>(anno_graph_, config_);
    }

    MaskedExtender<Extender> build_extender() const {
        return MaskedExtender<Extender>(anno_graph_, config_);
    }

    virtual void
    align_aggregate(const AlignmentGenerator &alignment_generator,
                    const std::function<void(DBGAlignment&&)> &callback) const override;

    const AnnotatedDBG &anno_graph_;
    DBGAlignerConfig config_;
};


template <class BaseSeeder>
void LabeledSeeder<BaseSeeder>::call_seeds(std::function<void(Seed&&)> callback) const {
    const auto &config = this->get_config();
    const auto &graph = this->get_graph();
    size_t k = graph.get_k();

    tsl::hopscotch_map<std::string, std::vector<Seed>> seeds;

    BaseSeeder::call_seeds([&](Seed&& seed) {
        std::string seed_seq = seed.get_offset()
            ? graph.get_node_sequence(seed.front()).substr(0, seed.get_offset())
                + seed.get_sequence()
            : seed.get_sequence();
        assert(seed_seq.size() >= k);
        for (const auto &[label, signature]
                : anno_graph_.get_top_label_signatures(
                      seed_seq,
                      anno_graph_.get_annotation().num_labels())) {
            size_t i = 0;
            while ((i = sdsl::util::next_bit(signature, i)) < signature.size()) {
                size_t next = sdsl::util::next_zero(signature, i);

                const char *begin = seed.get_query().data() + i;
                size_t size = next - i + k - 1 -
                    (seed.get_offset() > i ? seed.get_offset() - i : 0);

                seeds[label].emplace_back(
                    std::string_view(begin, size),
                    std::vector<node_index>(seed.begin() + i,
                                            seed.begin() + next),
                    config.match_score(std::string_view(begin, size)),
                    seed.get_clipping() + i,
                    seed.get_orientation(),
                    i <= seed.get_offset() ? seed.get_offset() - i : 0
                );

                i = next;
            }
        }
    });

    const auto &query = this->get_query();

    for (auto it = seeds.begin(); it != seeds.end(); ++it) {
        size_t kmer_count = std::accumulate(it.value().begin(), it.value().end(),
                                            size_t(0),
                                            [&](size_t sum, const auto &seed) {
                                                return sum + seed.size();
                                            });
        assert(kmer_count <= query.size() - k + 1);
        if (config.exact_kmer_match_fraction * (query.size() - k + 1) > kmer_count)
            continue;

        for (auto&& seed : it.value()) {
            seed.get_labels() = { it->first };
            callback(std::move(seed));
        }
    }
}

template <class BaseExtender>
bool MaskedColumnExtender<BaseExtender>
::set_labels(const std::vector<std::string> &labels) {
    std::vector<uint64_t> label_codes(labels.size());
    std::transform(labels.begin(), labels.end(), label_codes.begin(),
                   [&](const std::string &label) { return label_encoder_.encode(label); });
    return set_labels(std::move(label_codes));
}

template <class BaseExtender>
bool MaskedColumnExtender<BaseExtender>
::set_labels(std::vector<uint64_t>&& label_codes) {
    if (label_codes.size() != label_codes_.size()
            || std::any_of(label_codes.begin(), label_codes.end(), [&](uint64_t code) {
                return !label_codes_.count(code);
            })) {
        label_codes_.clear();
        label_codes_.insert(label_codes.begin(), label_codes.end());
        return true;
    }

    return false;
}

template <class BaseExtender>
void MaskedColumnExtender<BaseExtender>
::operator()(std::function<void(DBGAlignment&&, node_index)> callback,
             score_t min_path_score) {
    const auto &path = this->get_seed();
    assert(path.get_labels().size());

    if (set_labels(path.get_labels()))
        this->reset();

    BaseExtender::operator()([&](DBGAlignment&& alignment, node_index node) {
        auto &labels = alignment.get_labels();
        labels.resize(label_codes_.size());
        std::transform(label_codes_.begin(), label_codes_.end(), labels.begin(),
                       [&](uint64_t code) { return label_encoder_.decode(code); });
        callback(std::move(alignment), node);
    }, min_path_score);
}

template <class BaseExtender>
auto MaskedColumnExtender<BaseExtender>
::emplace_node(node_index node, node_index incoming_node, char c, size_t size,
               size_t best_pos, size_t last_priority_pos)
                -> std::pair<typename DPTable<node_index>::iterator, bool> {
    assert(this->get_graph().traverse(incoming_node, c) == node);

    if (graph_is_canonical()) {
        // make sure that node and incoming_node are in canonical_map_
        auto a_find = canonical_map_.find(incoming_node);
        auto b_find = canonical_map_.find(node);

        if (a_find == canonical_map_.end() || b_find == canonical_map_.end()) {
            std::array<node_index, 2> in_nodes { incoming_node, node};
            auto it = in_nodes.begin();
            this->get_graph().map_to_nodes(
                this->get_graph().get_node_sequence(incoming_node) + c,
                [&](auto i) {
                    assert(i != DeBruijnGraph::npos);
                    assert(it != in_nodes.end());
                    canonical_map_[*it] = i;
                    ++it;
                }
            );
            assert(it == in_nodes.end());
        }
    }

    return BaseExtender::emplace_node(node, incoming_node, c, size,
                                      best_pos, last_priority_pos);
}

template <class BaseExtender>
auto MaskedColumnExtender<BaseExtender>
::fork_extension(node_index node,
                 std::function<void(DBGAlignment&&, node_index)> callback,
                 score_t min_path_score) -> std::deque<std::pair<node_index, char>> {
    assert(node != DeBruijnGraph::npos);
    std::deque<std::pair<node_index, char>> out_columns;
    for (const auto &[next_node, c]
            : BaseExtender::fork_extension(node, callback, min_path_score)) {
        auto labeled_node = get_labeled_node(next_node);
        assert(!graph_is_canonical() || labeled_node != DeBruijnGraph::npos);
        if (labeled_node == DeBruijnGraph::npos)
            continue;

        auto row = matrix_.get_row(AnnotatedDBG::graph_to_anno_index(labeled_node));
        size_t count = std::count_if(row.begin(), row.end(),
                                     [&](uint64_t code) {
                                         return label_codes_.count(code);
                                     });

        if (count == label_codes_.size())
            out_columns.emplace_back(next_node, c);
    }

    return out_columns;
}

template <class BaseExtender>
auto MaskedColumnExtender<BaseExtender>::get_labeled_node(node_index node) -> node_index {
    assert(node != DeBruijnGraph::npos);
    if (!graph_is_canonical())
        return node;

    auto find = canonical_map_.find(node);
    if (find == canonical_map_.end()) {
        this->get_graph().map_to_nodes(
            this->get_graph().get_node_sequence(node),
            [&](auto i) {
                if (i != DeBruijnGraph::npos)
                    find = canonical_map_.emplace(node, i).first;
            }
        );
        if (find == canonical_map_.end())
            return DeBruijnGraph::npos;
    }

#ifndef NDEBUG
    this->get_graph().map_to_nodes(this->get_graph().get_node_sequence(node),
                                   [&](auto i) { assert(find->second == i); });
#endif

    return find->second;
}

template <class Seeder, class Extender, class AlignmentCompare,
          template <class BaseSeeder> class LabeledSeeder,
          template <class BaseExtender> class MaskedExtender>
void MaskedDBGAligner<Seeder, Extender, AlignmentCompare, LabeledSeeder, MaskedExtender>
::align_aggregate(const AlignmentGenerator &alignment_generator,
                  const std::function<void(DBGAlignment&&)> &callback) const {
    typedef boost::container::priority_deque<DBGAlignment,
                                             std::vector<DBGAlignment>,
                                             AlignmentCompare> AlignmentQueue;
    tsl::hopscotch_map<std::string, AlignmentQueue> path_queues;
    size_t num_alternative_paths = config_.num_alternative_paths;

    alignment_generator(
        [&](DBGAlignment&& alignment) {
#ifndef NDEBUG
            // check labels
            if (!alignment.get_offset()) {
                auto seq_labels = anno_graph_.get_labels(alignment.get_sequence(), 1.0);
                const auto &align_labels = alignment.get_labels();
                assert(seq_labels.size() >= align_labels.size());
                assert(std::all_of(align_labels.begin(), align_labels.end(),
                    [&](const auto &label) {
                        return std::find(seq_labels.begin(), seq_labels.end(), label)
                            != seq_labels.end();
                    }
                ));
            }
#endif
            // put this alignment into the path queues corresponding to its labels
            for (const auto &label : alignment.get_labels()) {
                auto single_label = alignment;
                single_label.get_labels() = { label };
                if (path_queues[label].size() < num_alternative_paths) {
                    path_queues[label].emplace(std::move(single_label));
                } else if (single_label.get_score()
                        > path_queues[label].minimum().get_score()) {
                    path_queues[label].update(path_queues[label].begin(),
                                              std::move(single_label));
                }
            }
        },
        [&](const DBGAlignment &alignment) {
            score_t min_path_score = std::numeric_limits<score_t>::max();
            bool found = false;
            for (const auto &label : alignment.get_labels()) {
                if (path_queues.count(label)) {
                    min_path_score = std::min(
                        min_path_score,
                        path_queues[label].minimum().get_score()
                    );
                    found = true;
                }
            }

            if (found)
                return std::max(min_path_score, config_.min_path_score);

            return config_.min_path_score;
        }
    );

    for (auto it = path_queues.begin(); it != path_queues.end(); ++it) {
        while (it->second.size()) {
            callback(DBGAlignment(it.value().maximum()));
            it.value().pop_maximum();
        }
    }
}

#endif // __MASKED_DBG_ALIGNER_HPP__
