#ifndef __GRAPH_CLEANING_HPP__
#define __GRAPH_CLEANING_HPP__

#include "representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector.hpp"
#include "graph_extensions/node_weights.hpp"


namespace mtg {
namespace graph {

bool is_unreliable_unitig(const std::vector<SequenceGraph::node_index> &path,
                          const NodeWeights &node_weights,
                          uint64_t min_median_abundance);

uint64_t estimate_min_kmer_abundance(const DeBruijnGraph &graph,
                                     const NodeWeights &node_weights,
                                     uint64_t num_singleton_kmers = 0);

uint64_t estimate_min_kmer_abundance_ztp(const DeBruijnGraph &graph,
                                         const NodeWeights &node_weights,
                                         uint64_t num_singleton_kmers = 0);

std::tuple<int64_t,double,double> estimate_min_kmer_abundance(const bit_vector &idx,
                                     const sdsl::int_vector<> &values,
                                     uint64_t num_singleton_kmers = 0,
                                     bool discard_last_count = false);

} // namespace graph
} // namespace mtg

#endif // __GRAPH_CLEANING_HPP__
