#include "test_dbg_helpers.hpp"

#include "gtest/gtest.h"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/bitmap/dbg_bitmap_construct.hpp"


namespace mtg {
namespace test {

using namespace mtg::graph::boss;


template <class Graph>
std::shared_ptr<DeBruijnGraph> make_graph_primary(std::shared_ptr<DeBruijnGraph> graph) {
    std::vector<std::string> contigs;
    graph->call_sequences([&](const std::string &contig, const auto &) {
        contigs.push_back(contig);
    }, 1, true);

    return std::make_shared<CanonicalDBG>(
        build_graph_batch<Graph>(graph->get_k(), contigs, BuildMode::BASE),
        true,
        2 // make the cache tiny to test out thread safety for all graphs
    );
}

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences,
            BuildMode mode) {
    auto graph = std::make_shared<Graph>(k, mode == BuildMode::CANONICAL);

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    if (mode == 2)
        return make_graph_primary<Graph>(graph);

    return graph;
}

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, BuildMode);

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashFast>(uint64_t, const std::vector<std::string> &, BuildMode);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashString>(uint64_t k,
                           const std::vector<std::string> &sequences,
                           BuildMode mode) {
    auto graph = std::make_shared<DBGHashString>(k);

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    if (mode == 2)
        return make_graph_primary<DBGHashString>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGBitmap>(uint64_t k,
                       const std::vector<std::string> &sequences,
                       BuildMode mode) {
    DBGBitmapConstructor constructor(k, mode == BuildMode::CANONICAL);
    for (const auto &sequence : sequences) {
        constructor.add_sequence(std::string(sequence));
    }

    auto graph = std::make_shared<DBGBitmap>(&constructor);

    if (mode == 2)
        return make_graph_primary<DBGBitmap>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinct>(uint64_t k,
                         const std::vector<std::string> &sequences,
                         BuildMode mode) {
    auto graph = std::make_shared<DBGSuccinct>(k, mode == BuildMode::CANONICAL);

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    graph->mask_dummy_kmers(1, false);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<1>>(uint64_t k,
                                   const std::vector<std::string> &sequences,
                                   BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).get_boss().index_suffix_ranges(1);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<2>>(uint64_t k,
                                   const std::vector<std::string> &sequences,
                                   BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)2));

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<10>>(uint64_t k,
                                    const std::vector<std::string> &sequences,
                                    BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)10));

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloomFPR<1, 1>>(uint64_t k,
                                       const std::vector<std::string> &sequences,
                                       BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloomFPR<1, 10>>(uint64_t k,
                                        const std::vector<std::string> &sequences,
                                        BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0 / 10);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloom<4, 1>>(uint64_t k,
                                    const std::vector<std::string> &sequences,
                                    BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 1);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloom<4, 50>>(uint64_t k,
                                     const std::vector<std::string> &sequences,
                                     BuildMode mode) {
    auto graph = build_graph<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 50);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}


template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences,
                  BuildMode mode) {
    auto graph = build_graph<Graph>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );

    if (mode == 2)
        return make_graph_primary<Graph>(graph);

    return graph;
}

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, BuildMode);

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashFast>(uint64_t, const std::vector<std::string> &, BuildMode);

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashString>(uint64_t, const std::vector<std::string> &, BuildMode);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGBitmap>(uint64_t k,
                             const std::vector<std::string> &sequences,
                             BuildMode mode) {
    DBGBitmapConstructor constructor(k, mode == BuildMode::CANONICAL);
    constructor.add_sequences(std::vector<std::string>(sequences));
    auto graph = std::make_shared<DBGBitmap>(&constructor);

    if (mode == 2)
        return make_graph_primary<DBGBitmap>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinct>(uint64_t k,
                               const std::vector<std::string> &sequences,
                               BuildMode mode) {
    BOSSConstructor constructor(k - 1, mode == BuildMode::CANONICAL);
    EXPECT_EQ(k - 1, constructor.get_k());
    constructor.add_sequences(std::vector<std::string>(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), mode == BuildMode::CANONICAL);
    graph->mask_dummy_kmers(1, false);
    EXPECT_EQ(k, graph->get_k());

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<1>>(uint64_t k,
                                         const std::vector<std::string> &sequences,
                                         BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).get_boss().index_suffix_ranges(1);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<2>>(uint64_t k,
                                         const std::vector<std::string> &sequences,
                                         BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)2));

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<10>>(uint64_t k,
                                          const std::vector<std::string> &sequences,
                                          BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph)
        .get_boss()
        .index_suffix_ranges(std::min(k - 1, (uint64_t)10));

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloomFPR<1, 1>>(uint64_t k,
                                             const std::vector<std::string> &sequences,
                                             BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloomFPR<1, 10>>(uint64_t k,
                                              const std::vector<std::string> &sequences,
                                              BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter_from_fpr(1.0 / 10);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloom<4, 1>>(uint64_t k,
                                          const std::vector<std::string> &sequences,
                                          BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 1);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloom<4, 50>>(uint64_t k,
                                           const std::vector<std::string> &sequences,
                                           BuildMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(
        k, sequences,
        mode == BuildMode::CANONICAL ? BuildMode::CANONICAL : BuildMode::BASE
    );
    dynamic_cast<DBGSuccinct&>(*graph).initialize_bloom_filter(4.0, 50);

    if (mode == 2)
        return make_graph_primary<DBGSuccinct>(graph);

    return graph;
}


template <class Graph>
bool check_graph(const std::string &alphabet, BuildMode mode, bool check_sequence) {
    std::vector<std::string> sequences;

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % alphabet.size()];
        }
        sequences.push_back(seq);
    }

    auto graph = build_graph<Graph>(20, sequences, mode);

    bool node_remap_failed = false;
    graph->call_nodes(
        [&graph, &node_remap_failed](DeBruijnGraph::node_index i) {
            if (graph->kmer_to_node(graph->get_node_sequence(i)) != i) {
                node_remap_failed = true;
                std::cerr << "Node failed\n"
                          << i << " " << graph->get_node_sequence(i) << "\n"
                          << graph->kmer_to_node(graph->get_node_sequence(i)) << " "
                          << graph->get_node_sequence(graph->kmer_to_node(graph->get_node_sequence(i))) << "\n";
            }
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

template bool check_graph<DBGSuccinct>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctIndexed<1>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctIndexed<2>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctIndexed<10>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctBloomFPR<1, 1>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctBloomFPR<1, 10>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctBloom<4, 1>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGSuccinctBloom<4, 50>>(const std::string &, BuildMode, bool);
template bool check_graph<DBGBitmap>(const std::string &, BuildMode, bool);
template bool check_graph<DBGHashOrdered>(const std::string &, BuildMode, bool);
template bool check_graph<DBGHashFast>(const std::string &, BuildMode, bool);
template bool check_graph<DBGHashString>(const std::string &, BuildMode, bool);

} // namespace test
} // namespace mtg
