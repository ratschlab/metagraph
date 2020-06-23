#include "test_dbg_helpers.hpp"

#include "gtest/gtest.h"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/bitmap/dbg_bitmap_construct.hpp"


namespace mtg {
namespace test {

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences,
            bool canonical) {
    std::shared_ptr<DeBruijnGraph> graph { new Graph(k, canonical) };

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    return graph;
}

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, bool);

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashFast>(uint64_t, const std::vector<std::string> &, bool);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashString>(uint64_t k,
            const std::vector<std::string> &sequences,
            bool) {
    std::shared_ptr<DeBruijnGraph> graph { new DBGHashString(k) };

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

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

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<1>>(uint64_t k,
                                   const std::vector<std::string> &sequences,
                                   bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).get_boss().index_suffix_ranges(1);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<2>>(uint64_t k,
                                   const std::vector<std::string> &sequences,
                                   bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)2));
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<10>>(uint64_t k,
                                    const std::vector<std::string> &sequences,
                                    bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)10));
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloomFPR<1, 1>>(uint64_t k,
                                       const std::vector<std::string> &sequences,
                                       bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloomFPR<1, 10>>(uint64_t k,
                                        const std::vector<std::string> &sequences,
                                        bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0 / 10);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloom<4, 1>>(uint64_t k,
                                    const std::vector<std::string> &sequences,
                                    bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 1);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloom<4, 50>>(uint64_t k,
                                     const std::vector<std::string> &sequences,
                                     bool canonical) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 50);
    return graph;
}


template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences,
                  bool canonical) {
    return build_graph<Graph>(k, sequences, canonical);
}

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, bool);

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashFast>(uint64_t, const std::vector<std::string> &, bool);

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashString>(uint64_t, const std::vector<std::string> &, bool);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGBitmap>(uint64_t k,
                             const std::vector<std::string> &sequences,
                             bool canonical) {
    DBGBitmapConstructor constructor(k, canonical);
    constructor.add_sequences(std::vector<std::string>(sequences));
    return std::shared_ptr<DeBruijnGraph>(new DBGBitmap(&constructor));
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinct>(uint64_t k,
                               const std::vector<std::string> &sequences,
                               bool canonical) {
    BOSSConstructor constructor(k - 1, canonical);
    EXPECT_EQ(k - 1, constructor.get_k());
    constructor.add_sequences(std::vector<std::string>(sequences));
    std::shared_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(&constructor), canonical) };
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);
    EXPECT_EQ(k, graph->get_k());
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<1>>(uint64_t k,
                                         const std::vector<std::string> &sequences,
                                         bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).get_boss().index_suffix_ranges(1);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<2>>(uint64_t k,
                                         const std::vector<std::string> &sequences,
                                         bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)2));
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<10>>(uint64_t k,
                                          const std::vector<std::string> &sequences,
                                          bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)10));
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloomFPR<1, 1>>(uint64_t k,
                                             const std::vector<std::string> &sequences,
                                             bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloomFPR<1, 10>>(uint64_t k,
                                              const std::vector<std::string> &sequences,
                                              bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0 / 10);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloom<4, 1>>(uint64_t k,
                                          const std::vector<std::string> &sequences,
                                          bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 1);
    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloom<4, 50>>(uint64_t k,
                                           const std::vector<std::string> &sequences,
                                           bool canonical) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, canonical);
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 50);
    return graph;
}


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

    bool node_remap_failed = false;
    graph->call_nodes(
        [&graph, &node_remap_failed](DeBruijnGraph::node_index i) {
            if (graph->kmer_to_node(graph->get_node_sequence(i)) != i)
                node_remap_failed = true;
        },
        [&node_remap_failed]() { return node_remap_failed; }
    );
    if (node_remap_failed)
        return false;

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
template bool check_graph<DBGSuccinctIndexed<1>>(const std::string &, bool, bool);
template bool check_graph<DBGSuccinctIndexed<2>>(const std::string &, bool, bool);
template bool check_graph<DBGSuccinctIndexed<10>>(const std::string &, bool, bool);
template bool check_graph<DBGSuccinctBloomFPR<1, 1>>(const std::string &, bool, bool);
template bool check_graph<DBGSuccinctBloomFPR<1, 10>>(const std::string &, bool, bool);
template bool check_graph<DBGSuccinctBloom<4, 1>>(const std::string &, bool, bool);
template bool check_graph<DBGSuccinctBloom<4, 50>>(const std::string &, bool, bool);
template bool check_graph<DBGBitmap>(const std::string &, bool, bool);
template bool check_graph<DBGHashOrdered>(const std::string &, bool, bool);
template bool check_graph<DBGHashFast>(const std::string &, bool, bool);
template bool check_graph<DBGHashString>(const std::string &, bool, bool);

} // namespace test
} // namespace mtg
