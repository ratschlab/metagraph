#ifndef __GRAPH_CLEANING_HPP__
#define __GRAPH_CLEANING_HPP__

#include "sequence_graph.hpp"
#include "weighted_graph.hpp"


bool is_unreliable_unitig(const std::string &sequence,
                          const DeBruijnGraph &graph,
                          const IWeighted<DeBruijnGraph::node_index> &node_weights,
                          uint64_t min_median_abundance);

uint64_t estimate_min_kmer_abundance(const DeBruijnGraph &graph,
                                     const IWeighted<DeBruijnGraph::node_index> &node_weights,
                                     uint64_t fallback_cutoff);

#endif // __GRAPH_CLEANING_HPP__
