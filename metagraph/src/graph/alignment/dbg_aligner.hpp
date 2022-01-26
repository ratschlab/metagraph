#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include "aligner_alignment.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "aligner_chainer.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {
namespace align {


class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::score_t score_t;

    typedef std::tuple<std::string /* header */,
                       std::string /* seq */,
                       bool /* orientation of seq */> Query;
    typedef std::function<void(std::string_view /* header */,
                               QueryAlignment&& /* alignments */)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    // Main aligner
    virtual void align_batch(const std::vector<Query> &seq_batch,
                             const AlignmentCallback &callback) const = 0;

    // Convenience method
    QueryAlignment align(std::string_view query, bool is_reverse_complement = false) const;
};

template <class AlignmentCompare = LocalAlignmentLess>
class ISeedAndExtendAligner : public IDBGAligner {
  public:
    ISeedAndExtendAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    virtual ~ISeedAndExtendAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const AlignmentCallback &callback) const override;

    const DBGAlignerConfig& get_config() const { return config_; }

  protected:
    typedef AlignmentAggregator<AlignmentCompare> Aggregator;
    typedef std::vector<std::tuple<std::shared_ptr<ISeeder>, std::vector<node_index>,
                                   std::shared_ptr<ISeeder>, std::vector<node_index>>> BatchSeeders;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;

    virtual std::shared_ptr<IExtender>
    build_extender(std::string_view query,
                   const Aggregator &aggregator,
                   const DBGAlignerConfig &config) const = 0;

    virtual void filter_seeds(BatchSeeders &seeders) const = 0;

    virtual std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<node_index> &nodes) const = 0;

    template <class Seeder>
    std::shared_ptr<ISeeder>
    build_seeder_impl(std::string_view query,
                      bool is_reverse_complement,
                      const std::vector<node_index> &nodes) const {
        return std::make_shared<SuffixSeeder<Seeder>>(
            graph_, query, is_reverse_complement, nodes, config_
        );
    }

    // Generates seeds and extends them. If force_fixed_seed is true, then
    // all alignments must have the seed as a prefix. Otherwise, only the first
    // node of the seed is used as an alignment starting node.
    virtual void align_core(const ISeeder &seeder,
                            IExtender &extender,
                            const std::function<void(Alignment&&)> &callback,
                            const std::function<score_t(const Alignment&)> &get_min_path_score,
                            bool force_fixed_seed) const;

  private:
    // Align the forward and reverse complement of the query sequence in both
    // directions and return the overall best alignment. e.g., for the forward query
    // 1. Find all seeds of its reverse complement
    // 2. Given a seed, extend forwards to get alignment A
    // 3. Reverse complement the alignment to get A', treat it like a new seed
    // 4. Extend A' forwards to get the final alignment A''
    std::tuple<size_t, size_t, size_t>
    align_both_directions(std::string_view forward,
                          std::string_view reverse,
                          const ISeeder &forward_seeder,
                          const ISeeder &reverse_seeder,
                          IExtender &forward_extender,
                          IExtender &reverse_extender,
                          const std::function<void(Alignment&&)> &callback,
                          const std::function<score_t(const Alignment&)> &get_min_path_score) const;

    // Construct a full alignment from a chain by aligning the query agaisnt
    // the graph in the regions of the query in between the chain seeds.
    void extend_chain(std::string_view query,
                      std::string_view query_rc,
                      Chain&& chain,
                      score_t score,
                      size_t &num_extensions,
                      size_t &num_explored_nodes,
                      const std::function<void(Alignment&&)> &callback) const;

    BatchSeeders
    build_seeders(const std::vector<Query> &seq_batch,
                  const std::vector<std::pair<QueryAlignment, Aggregator>> &wrapped_seqs) const;
};

template <class Extender = DefaultColumnExtender,
          class Seeder = UniMEMSeeder,
          class AlignmentCompare = LocalAlignmentLess>
class DBGAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    DBGAligner(Args&&... args)
          : ISeedAndExtendAligner<AlignmentCompare>(std::forward<Args>(args)...) {}

  private:
    typedef typename ISeedAndExtendAligner<AlignmentCompare>::Aggregator Aggregator;
    typedef typename ISeedAndExtendAligner<AlignmentCompare>::BatchSeeders BatchSeeders;

    std::shared_ptr<IExtender>
    build_extender(std::string_view query, const Aggregator&, const DBGAlignerConfig &config) const override final {
        return std::make_shared<Extender>(this->graph_, config, query);
    }

    std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<IDBGAligner::node_index> &nodes) const override final {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }

    void filter_seeds(BatchSeeders &) const override final {}
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
