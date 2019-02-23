#include "sequence_graph.hpp"

#include <cassert>


DeBruijnGraph::node_index DeBruijnGraph::kmer_to_node(const char *begin) const {
    node_index node = npos;
    map_to_nodes(std::string(begin, get_k()),
        [&node](node_index i) { node = i; });
    return node;
}

DeBruijnGraph::node_index DeBruijnGraph::kmer_to_node(const std::string &kmer) const {
    assert(kmer.size() == get_k());
    return kmer_to_node(kmer.data());
}

// Check whether graph contains fraction of nodes from the sequence
bool DeBruijnGraph::find(const std::string &sequence,
                         double discovery_fraction) const {
    if (sequence.length() < get_k())
        return false;

    const size_t num_kmers = sequence.length() - get_k() + 1;
    const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    map_to_nodes(sequence,
        [&](node_index node) {
            if (node) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }
        },
        [&]() { return num_kmers_missing > max_kmers_missing
                        || num_kmers_discovered >= min_kmers_discovered; }
    );

    return num_kmers_missing <= max_kmers_missing;
}
