#ifndef __GRAPH_CLEANING_HPP__
#define __GRAPH_CLEANING_HPP__

#include "representation/base/sequence_graph.hpp"
#include "graph_extensions/node_weights.hpp"


bool is_unreliable_unitig(const std::vector<SequenceGraph::node_index> &path,
                          const NodeWeights &node_weights,
                          uint64_t min_median_abundance);

uint64_t estimate_min_kmer_abundance(const DeBruijnGraph &graph,
                                     const NodeWeights &node_weights,
                                     uint64_t num_singleton_kmers = 0);

#endif // __GRAPH_CLEANING_HPP__
