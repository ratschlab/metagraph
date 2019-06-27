#include "graph_cleaning.hpp"

#include <cmath>

#include "utils.hpp"


bool is_unreliable_unitig(const std::string &sequence,
                          const DeBruijnGraph &graph,
                          const IWeighted<DeBruijnGraph::node_index> &node_weights,
                          uint64_t min_median_abundance) {
    if (min_median_abundance <= 1)
        return false;

    const uint64_t num_kmers = sequence.size() - graph.get_k() + 1;
    uint64_t num_weak_kmers = 0;
    uint64_t num_reliable_kmers = 0;

    graph.map_to_nodes(sequence,
        [&](auto node) {
            if (node_weights.get_weight(node) < min_median_abundance) {
                num_weak_kmers++;
            } else {
                num_reliable_kmers++;
            }
        },
        [&]() { return num_weak_kmers * 2 > num_kmers
                        || num_reliable_kmers * 2 >= num_kmers; }
    );

    assert(num_weak_kmers + num_reliable_kmers <= num_kmers);

    // check if median is smaller than the threshold
    return num_weak_kmers * 2 > num_kmers;
}


uint64_t estimate_min_kmer_abundance(const bitmap &node_mask,
                                     const IWeighted<DeBruijnGraph::node_index> &node_weights,
                                     uint64_t fallback_cutoff) {
    // TODO: implement
    std::ignore = node_mask;
    std::ignore = node_weights;
    return fallback_cutoff;
}
