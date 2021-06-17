#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include "aligner_alignment.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {


class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    typedef std::tuple<std::string /* header */,
                       std::string /* seq */,
                       bool /* orientation of seq */> Query;
    typedef std::function<void(std::string_view /* header */,
                               DBGQueryAlignment&& /* alignments */)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    // Main aligner
    virtual void align_batch(const std::vector<Query> &seq_batch,
                             const AlignmentCallback &callback) const = 0;

    // Convenience method
    DBGQueryAlignment align(std::string_view query, bool is_reverse_complement = false) const;
};

template <class AlignmentCompare = LocalAlignmentLess>
class ISeedAndExtendAligner : public IDBGAligner {
  public:
    ISeedAndExtendAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {
        if (!config_.min_seed_length)
            config_.min_seed_length = graph_.get_k();

        if (!config_.max_seed_length)
            config_.max_seed_length = graph_.get_k();

        assert(config_.max_seed_length >= config_.min_seed_length);
        assert(config_.num_alternative_paths);
        assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
            && "primary graphs must be wrapped into canonical");

        if (!config_.check_config_scores())
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
    }

    virtual ~ISeedAndExtendAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const AlignmentCallback &callback) const override;

    const DBGAlignerConfig& get_config() const { return config_; }

  protected:
    const DeBruijnGraph &graph_;

    virtual std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query,
                   const AlignmentAggregator<IDBGAligner::node_index, AlignmentCompare> &aggregator) const = 0;

    virtual std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<DeBruijnGraph::node_index> &nodes) const = 0;

    template <class Seeder>
    std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder_impl(std::string_view query,
                      bool is_reverse_complement,
                      const std::vector<DeBruijnGraph::node_index> &nodes) const {
        return std::make_shared<SuffixSeeder<Seeder>>(
            graph_, query, is_reverse_complement, nodes, config_
        );
    }

  private:
    DBGAlignerConfig config_;
};

template <class Extender = DefaultColumnExtender<>,
          class Seeder = UniMEMSeeder<>,
          class AlignmentCompare = LocalAlignmentLess>
class DBGAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    DBGAligner(Args&&... args)
          : ISeedAndExtendAligner<AlignmentCompare>(std::forward<Args>(args)...) {}

  private:
    std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query,
                   const AlignmentAggregator<IDBGAligner::node_index, AlignmentCompare>&) const override {
        return std::make_shared<Extender>(this->graph_, this->get_config(), query);
    }

    std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<DeBruijnGraph::node_index> &nodes) const override {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
