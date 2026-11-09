#ifndef __TESTS_GRAPH_TRAVERSAL_TRAVERSAL_TEST_HELPERS_HPP__
#define __TESTS_GRAPH_TRAVERSAL_TRAVERSAL_TEST_HELPERS_HPP__

#include <random>
#include <string>
#include <vector>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace test {

// |n| random nucleotide sequences of length |len|, drawn from one generator seeded with |seed|
inline std::vector<std::string> make_sequences(size_t n, size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::vector<std::string> seqs;
    for (size_t i = 0; i < n; ++i) {
        std::string s(len, 'A');
        for (char &c : s) {
            c = "ACGT"[gen() % 4];
        }
        seqs.push_back(s);
    }
    return seqs;
}

inline std::string random_seq(size_t len, uint32_t seed) {
    return make_sequences(1, len, seed)[0];
}

inline std::string rc(std::string s) {
    ::reverse_complement(s);
    return s;
}

// the graph modes of this build: BASIC only for protein graphs
inline std::vector<graph::DeBruijnGraph::Mode> all_modes() {
    return {
        graph::DeBruijnGraph::BASIC,
#if ! _PROTEIN_GRAPH
        graph::DeBruijnGraph::CANONICAL,
        graph::DeBruijnGraph::PRIMARY,
#endif
    };
}

} // namespace test
} // namespace mtg

#endif // __TESTS_GRAPH_TRAVERSAL_TRAVERSAL_TEST_HELPERS_HPP__
