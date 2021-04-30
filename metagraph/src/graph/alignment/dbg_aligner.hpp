#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include <tsl/hopscotch_map.h>

#include "aligner_alignment.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "aligner_aggregator.hpp"
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

    typedef std::function<void(std::string_view /* header */,
                               std::string_view /* seq */,
                               bool /* orientation of seq */)> QueryCallback;
    typedef std::function<void(const QueryCallback&)> QueryGenerator;
    typedef std::function<void(std::string_view /* header */,
                               DBGQueryAlignment&& /* alignments */)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    // Main aligner
    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const = 0;

    // Convenience methods
    DBGQueryAlignment align(std::string_view query,
                            bool is_reverse_complement = false) const;
    void align_batch(const std::vector<std::pair<std::string, std::string>> &seq_batch,
                     const AlignmentCallback &callback) const;
};

template <class AlignmentCompare = LocalAlignmentLess>
class ISeedAndExtendAligner : public IDBGAligner {
  public:
    virtual ~ISeedAndExtendAligner() {}
    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override;

  protected:
    virtual std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query) const = 0;

    virtual std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 std::vector<DeBruijnGraph::node_index>&& nodes) const = 0;
};

template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class DBGAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    DBGAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores())
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
    }

    virtual const DeBruijnGraph& get_graph() const override final { return graph_; }
    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;

    std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query) const override {
        return std::make_shared<Extender>(graph_, config_, query);
    }

    std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 std::vector<DeBruijnGraph::node_index>&& nodes) const override {
        if (config_.min_seed_length < graph_.get_k()
                && SuffixSeeder<Seeder>::get_base_dbg_succ(graph_)) {
            return std::make_shared<SuffixSeeder<Seeder>>(
                graph_, query, is_reverse_complement, std::move(nodes), config_
            );
        } else {
            return std::make_shared<Seeder>(
                graph_, query, is_reverse_complement, std::move(nodes), config_
            );
        }
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
