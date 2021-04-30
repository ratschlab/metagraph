#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_set.hpp"
#include "common/utils/template_utils.hpp"
#include "common/hashers/hash.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace graph {
namespace align {

template <class AlignmentCompare = LocalAlignmentLess>
class ILabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    ILabeledAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : anno_graph_(anno_graph), graph_(anno_graph_.get_graph()), config_(config) {}

    virtual ~ILabeledAligner() {}

    virtual void align_batch(const IDBGAligner::QueryGenerator &generate_query,
                             const IDBGAligner::AlignmentCallback &callback) const override {
        ISeedAndExtendAligner<AlignmentCompare>::align_batch(generate_query,
            [&](std::string_view header, IDBGAligner::DBGQueryAlignment&& alignments) {
                auto it = std::remove_if(
                    alignments.begin(), alignments.end(),
                    [](const auto &a) { return a.target_columns.empty(); }
                );
                alignments.erase(it, alignments.end());
                callback(header, std::move(alignments));
            }
        );
    }

    virtual const DeBruijnGraph& get_graph() const override final { return graph_; }
    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;
};

template <typename NodeType = DeBruijnGraph::node_index>
class LabeledBacktrackingExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef DefaultColumnExtender<DeBruijnGraph::node_index> BaseExtender;
    typedef typename BaseExtender::score_t score_t;
    typedef typename BaseExtender::node_index node_index;
    typedef typename BaseExtender::AlignNode AlignNode;
    typedef typename BaseExtender::AlignNodeHash AlignNodeHash;
    typedef typename BaseExtender::DBGAlignment DBGAlignment;

    LabeledBacktrackingExtender(const AnnotatedDBG &anno_graph,
                                const DBGAlignerConfig &config,
                                std::string_view query)
          : BaseExtender(anno_graph.get_graph(), config, query), anno_graph_(anno_graph) {
        targets_set_.emplace(Vector<uint64_t>{});
    }

    virtual ~LabeledBacktrackingExtender() {}

  protected:
    virtual std::vector<AlignNode>
    backtrack(score_t min_path_score,
              AlignNode best_node,
              tsl::hopscotch_set<AlignNode, AlignNodeHash> &prev_starts,
              std::vector<DBGAlignment> &extensions) const override;

    virtual void init_backtrack() const override;

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &,
                                      const AlignNode &) const override { return false; }

  private:
    const AnnotatedDBG &anno_graph_;
    mutable VectorSet<Vector<uint64_t>, utils::VectorHash> targets_set_;
    mutable tsl::hopscotch_map<node_index, size_t> targets_;
};

template <class Seeder = ExactSeeder<>,
          class Extender = LabeledBacktrackingExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public ILabeledAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    LabeledAligner(Args&&... args)
          : ILabeledAligner<AlignmentCompare>(std::forward<Args>(args)...) {
        assert(this->config_.num_alternative_paths);
        if (!this->config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

  protected:
    std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query) const override {
        return std::make_shared<Extender>(this->anno_graph_, this->config_, query);
    }

    std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 std::vector<DeBruijnGraph::node_index>&& nodes) const override {
        if (this->config_.min_seed_length < this->graph_.get_k()
                && SuffixSeeder<Seeder>::get_base_dbg_succ(this->graph_)) {
            return std::make_shared<SuffixSeeder<Seeder>>(
                this->graph_, query, is_reverse_complement, std::move(nodes), this->config_
            );
        } else {
            return std::make_shared<Seeder>(
                this->graph_, query, is_reverse_complement, std::move(nodes), this->config_
            );
        }
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
