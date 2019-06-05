#include "test_dbg_helpers.hpp"

#include "gtest/gtest.h"
#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "boss_construct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"


template <class Graph>
std::unique_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences,
            bool canonical) {
    std::unique_ptr<DeBruijnGraph> graph { new Graph(k, canonical) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence);
    }

    return graph;
}

template <>
std::unique_ptr<DeBruijnGraph>
build_graph<DBGHashString>(uint64_t k,
            const std::vector<std::string> &sequences,
            bool) {
    std::unique_ptr<DeBruijnGraph> graph { new DBGHashString(k) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence);
    }

    return graph;
}

template std::unique_ptr<DeBruijnGraph> build_graph<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, bool);

template <>
std::unique_ptr<DeBruijnGraph>
build_graph<DBGBitmap>(uint64_t k,
                   const std::vector<std::string> &sequences,
                   bool canonical) {
    DBGBitmapConstructor constructor(k, canonical);
    for (const auto &sequence : sequences) {
        constructor.add_sequence(std::string(sequence));
    }

    return std::unique_ptr<DeBruijnGraph>(new DBGBitmap(&constructor));
}

template <>
std::unique_ptr<DeBruijnGraph>
build_graph<DBGSuccinct>(uint64_t k,
                         const std::vector<std::string> &sequences,
                         bool canonical) {
    std::unique_ptr<DeBruijnGraph> graph { new DBGSuccinct(k, canonical) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(std::string(sequence));
    }

    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    return graph;
}


template <class Graph>
std::unique_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences,
                  bool canonical) {
    std::unique_ptr<DeBruijnGraph> graph { new Graph(k, canonical) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(std::string(sequence));
    }
    return graph;
}

template <>
std::unique_ptr<DeBruijnGraph>
build_graph_batch<DBGHashString>(uint64_t k,
                  const std::vector<std::string> &sequences,
                  bool) {
    std::unique_ptr<DeBruijnGraph> graph { new DBGHashString(k) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(std::string(sequence));
    }
    return graph;
}

template std::unique_ptr<DeBruijnGraph> build_graph_batch<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, bool);

template <>
std::unique_ptr<DeBruijnGraph>
build_graph_batch<DBGBitmap>(uint64_t k,
                         const std::vector<std::string> &sequences,
                         bool canonical) {
    DBGBitmapConstructor constructor(k, canonical);
    constructor.add_sequences(sequences);
    return std::unique_ptr<DeBruijnGraph>(new DBGBitmap(&constructor));
}

template <>
std::unique_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinct>(uint64_t k,
                               const std::vector<std::string> &sequences,
                               bool canonical) {
    BOSSConstructor constructor(k - 1, canonical);
    constructor.add_sequences(sequences);
    std::unique_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(&constructor)) };
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);
    return graph;
}

template <class Graph>
bool check_graph(const std::string &alphabet, bool canonical) {
    std::vector<std::string> sequences;

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
        }
        sequences.push_back(seq);
    }

    auto graph = build_graph<Graph>(20, sequences, canonical);

    auto it = DeBruijnGraph::npos;
    graph->call_nodes(
        [&](const auto &i) { it = i; },
        [&]() {
            if (it == DeBruijnGraph::npos)
                return false;

            if (it != graph->kmer_to_node(graph->get_node_sequence(it))) {
                it = DeBruijnGraph::npos;
                return true;
            }

            return false;
        }
    );

    return it != DeBruijnGraph::npos;
}

template bool check_graph<DBGSuccinct>(const std::string &, bool);
template bool check_graph<DBGBitmap>(const std::string &, bool);
template bool check_graph<DBGHashOrdered>(const std::string &, bool);
template bool check_graph<DBGHashString>(const std::string &, bool);
