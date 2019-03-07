#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <vector>
#include <functional>

#include "sequence_graph.hpp"
#include "aligner_helper.hpp"
#include "aligner_methods.hpp"


class DBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef typename DBGAlignment::score_t score_t;

    explicit DBGAligner(const DeBruijnGraph &graph,
                        const DBGAlignerConfig &config = DBGAlignerConfig(),
                        const Seeder &seed = default_seeder,
                        const Extender &extend = default_extender,
                        const PriorityFunction &priority_function = std::less<DBGAlignment>());

    explicit DBGAligner(const DeBruijnGraph &graph,
                        const Config &config,
                        const Seeder &seed = default_seeder,
                        const Extender &extend = default_extender,
                        const PriorityFunction &priority_function = std::less<DBGAlignment>());

    // Align a sequence to the graph
    template <class StringIt>
    void align(StringIt query_begin_it,
               StringIt query_end_it,
               const std::function<void(DBGAlignment&&)> &callback,
               bool orientation = false,
               score_t min_path_score = std::numeric_limits<score_t>::min()) const;

    // A convenience function
    std::vector<DBGAlignment>
    align(const std::string &query,
          bool orientation = false,
          score_t min_path_score = std::numeric_limits<score_t>::min()) const;

    std::vector<DBGAligner::DBGAlignment>
    align_forward_and_reverse_complement(const std::string &query,
                                         const std::string &reverse_complement_query,
                                         score_t min_path_score
                                             = std::numeric_limits<DBGAligner::score_t>::min()) const;

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }
    const Seeder& get_seeder() const { return seed_; }
    const Extender& get_extender() const { return extend_; }
    const PriorityFunction& get_priority_function() const { return priority_function_; }

  private:
    const DeBruijnGraph& graph_;
    DBGAlignerConfig config_;
    const Seeder seed_;
    const Extender extend_;
    const PriorityFunction priority_function_;
};

std::vector<DBGAligner::DBGAlignment>
extend_mapping_forward_and_reverse_complement(const std::string &query,
                                              const std::string &reverse_complement_query,
                                              const DeBruijnGraph &graph,
                                              const DBGAlignerConfig &config = DBGAlignerConfig(),
                                              DBGAligner::score_t min_path_score
                                                  = std::numeric_limits<DBGAligner::score_t>::min(),
                                              const Extender &extend = default_extender,
                                              const PriorityFunction &priority_function
                                                  = std::less<DBGAligner::DBGAlignment>());


#endif // __DBG_ALIGNER_HPP__
