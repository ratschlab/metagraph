#include "test_dbg_helpers.hpp"

#include "gtest/gtest.h"
#include "boss.hpp"
#include "boss_construct.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"


template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences,
            bool canonical) {
    std::shared_ptr<DeBruijnGraph> graph { new Graph(k, canonical) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence);
    }
    return graph;
}

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, bool);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashString>(uint64_t k,
            const std::vector<std::string> &sequences,
            bool) {
    std::shared_ptr<DeBruijnGraph> graph { new DBGHashString(k) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence);
    }
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGBitmap>(uint64_t k,
                       const std::vector<std::string> &sequences,
                       bool canonical) {
    DBGBitmapConstructor constructor(k, canonical);
    for (const auto &sequence : sequences) {
        constructor.add_sequence(std::string(sequence));
    }
    return std::shared_ptr<DeBruijnGraph>(new DBGBitmap(&constructor));
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinct>(uint64_t k,
                         const std::vector<std::string> &sequences,
                         bool canonical) {
    std::shared_ptr<DeBruijnGraph> graph { new DBGSuccinct(k, canonical) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(std::string(sequence));
    }
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);
    return graph;
}


template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences,
                  bool canonical) {
    std::shared_ptr<DeBruijnGraph> graph { new Graph(k, canonical) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(std::string(sequence));
    }
    return graph;
}

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, bool);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashString>(uint64_t k,
                                 const std::vector<std::string> &sequences,
                                 bool) {
    std::shared_ptr<DeBruijnGraph> graph { new DBGHashString(k) };
    for (const auto &sequence : sequences) {
        graph->add_sequence(std::string(sequence));
    }
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGBitmap>(uint64_t k,
                             const std::vector<std::string> &sequences,
                             bool canonical) {
    DBGBitmapConstructor constructor(k, canonical);
    constructor.add_sequences(sequences);
    return std::shared_ptr<DeBruijnGraph>(new DBGBitmap(&constructor));
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinct>(uint64_t k,
                               const std::vector<std::string> &sequences,
                               bool canonical) {
    BOSSConstructor constructor(k - 1, canonical);
    EXPECT_EQ(k - 1, constructor.get_k());
    constructor.add_sequences(sequences);
    std::shared_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(&constructor), canonical) };
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);
    EXPECT_EQ(k, graph->get_k());
    return graph;
}

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_iterative(uint64_t k,
                      std::function<void(std::function<void(const std::string&)>)> generate,
                      bool canonical) {
    std::vector<std::string> sequences;
    generate([&](const auto &sequence) { sequences.push_back(sequence); });
    return build_graph_batch<Graph>(k, sequences, canonical);
}

template
std::shared_ptr<DeBruijnGraph>
build_graph_iterative<DBGHashOrdered>(uint64_t, std::function<void(std::function<void(const std::string&)>)>, bool);

template
std::shared_ptr<DeBruijnGraph>
build_graph_iterative<DBGHashString>(uint64_t, std::function<void(std::function<void(const std::string&)>)>, bool);

template
std::shared_ptr<DeBruijnGraph>
build_graph_iterative<DBGBitmap>(uint64_t, std::function<void(std::function<void(const std::string&)>)>, bool);

template
std::shared_ptr<DeBruijnGraph>
build_graph_iterative<DBGSuccinct>(uint64_t, std::function<void(std::function<void(const std::string&)>)>, bool);



template <class Graph>
bool check_graph(const std::string &alphabet, bool canonical, bool check_sequence) {
    std::vector<std::string> sequences;

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % alphabet.size()];
        }
        sequences.push_back(seq);
    }

    auto graph = build_graph<Graph>(20, sequences, canonical);

    const auto nnodes = graph->num_nodes();
    for (DeBruijnGraph::node_index i = 1; i <= nnodes; ++i) {
        if (graph->kmer_to_node(graph->get_node_sequence(i)) != i)
            return false;
    }

    if (!check_sequence)
        return true;

    for (const auto &seq : sequences) {
        bool stop = false;
        graph->map_to_nodes(
            seq,
            [&](const auto &i) {
                stop = !i || graph->kmer_to_node(graph->get_node_sequence(i)) != i;
            },
            [&]() { return stop; }
        );

        if (stop)
            return false;
    }

    return true;
}

template bool check_graph<DBGSuccinct>(const std::string &, bool, bool);
template bool check_graph<DBGBitmap>(const std::string &, bool, bool);
template bool check_graph<DBGHashOrdered>(const std::string &, bool, bool);
template bool check_graph<DBGHashString>(const std::string &, bool, bool);
