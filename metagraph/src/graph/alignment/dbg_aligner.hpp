#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include "alignment.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "aligner_chainer.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

class GraphTopology;
class GraphSeqAnnotator;

namespace align {

class IDBGAligner {
  public:
    typedef std::pair<std::string /* header */, std::string /* seq */> Query;
    typedef std::function<void(const std::string& /* header */,
                               AlignmentResults&& /* alignments */)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;

    // Main aligner
    virtual void align_batch(const std::vector<Query> &seq_batch,
                             const AlignmentCallback &callback) const = 0;

    // Convenience method
    AlignmentResults align(std::string_view query) const;

    virtual bool has_coordinates() const = 0;
    virtual bool has_local_coordinates() const = 0;

    virtual std::unique_ptr<SeedFilteringExtender> build_extender(std::string_view query) const = 0;

    virtual const GraphSeqAnnotator* get_seq_annotator() const = 0;
};


template <class Seeder = SuffixSeeder<UniMEMSeeder>,
          class Extender = DefaultColumnExtender,
          class AlignmentCompare = LocalAlignmentLess>
class DBGAligner : public IDBGAligner {
  public:
    DBGAligner(const DeBruijnGraph &graph,
               const DBGAlignerConfig &config,
               std::shared_ptr<const GraphSeqAnnotator> seq_annotator = {});

    virtual ~DBGAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const AlignmentCallback &callback) const override;

    const DeBruijnGraph& get_graph() const override { return graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

    virtual bool has_coordinates() const override { return false; }
    virtual bool has_local_coordinates() const override { return false; }

    virtual std::unique_ptr<SeedFilteringExtender> build_extender(std::string_view query) const {
        return std::make_unique<Extender>(*this, query);
    }

    virtual const GraphSeqAnnotator* get_seq_annotator() const override final { return seq_annotator_.get(); }

  protected:
    typedef typename Seeder::node_index node_index;
    typedef Alignment::score_t score_t;

    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;

    std::shared_ptr<const GraphSeqAnnotator> seq_annotator_;

    typedef std::vector<std::pair<std::unique_ptr<ISeeder>, std::unique_ptr<ISeeder>>> BatchSeeders;

    virtual BatchSeeders build_seeders(const std::vector<Query> &seq_batch,
                                       const std::vector<AlignmentResults> &wrapped_seqs,
                                       std::vector<std::pair<std::vector<Seed>, std::vector<Seed>>> &discarded_seeds) const;

  private:
    /**
     * Align the forward and reverse complement of the query sequence in both
     * directions and return the overall best alignment. e.g., for the forward query
     * 1. Find all seeds of its reverse complement
     * 2. Given a seed, extend forwards to get alignment A
     * 3. Reverse complement the alignment to get A', treat it like a new seed
     * 4. Extend A' forwards to get the final alignment A''
     */
    std::tuple<size_t, size_t, size_t>
    align_both_directions(std::string_view forward,
                          std::string_view reverse,
                          const ISeeder &forward_seeder,
                          const std::unique_ptr<ISeeder> &reverse_seeder,
                          Extender &forward_extender,
                          Extender &reverse_extender,
                          const std::function<void(Alignment&&)> &callback,
                          const std::function<void(Alignment&&)> &callback_discarded,
                          const std::function<score_t()> &get_min_path_score) const;

    // Construct a full alignment from a chain by aligning the query agaisnt
    // the graph in the regions of the query in between the chain seeds.
    void extend_chain(std::string_view query,
                      std::string_view query_rc,
                      Extender &forward_extender,
                      Chain&& chain,
                      size_t &num_extensions,
                      size_t &num_explored_nodes,
                      const std::function<void(Alignment&&)> &callback) const;
};

std::pair<Alignment, Alignment> split_seed(const DeBruijnGraph &graph,
                                           const DBGAlignerConfig &config,
                                           const Alignment &alignment);

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
