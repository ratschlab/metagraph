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
size_t max_test_k() {
    return 256 / kmer::KmerExtractorBOSS::bits_per_char;
}
template size_t max_test_k<DBGSuccinct>();
template size_t max_test_k<DBGSuccinctIndexed<1>>();
template size_t max_test_k<DBGSuccinctIndexed<2>>();
template size_t max_test_k<DBGSuccinctIndexed<10>>();
template size_t max_test_k<DBGSuccinctBloomFPR<1, 1>>();
template size_t max_test_k<DBGSuccinctBloomFPR<1, 10>>();
template size_t max_test_k<DBGSuccinctBloom<4, 1>>();
template size_t max_test_k<DBGSuccinctBloom<4, 50>>();

template<> size_t max_test_k<DBGBitmap>() {
    return 63. / kmer::KmerExtractor2Bit().bits_per_char;
}
template<> size_t max_test_k<DBGHashOrdered>() {
    return 256 / kmer::KmerExtractor2Bit::bits_per_char;
}
template<> size_t max_test_k<DBGHashFast>() {
    return 256 / kmer::KmerExtractor2Bit::bits_per_char;
}
template<> size_t max_test_k<DBGHashString>() {
    return 100;
}

template <class Graph>
std::shared_ptr<DeBruijnGraph> make_graph_primary(std::shared_ptr<DeBruijnGraph> graph) {
    std::vector<std::string> contigs;
    graph->call_sequences([&](const std::string &contig, const auto &) {
        contigs.push_back(contig);
    }, 1, true);

    return build_graph_batch<Graph>(graph->get_k(), contigs, DBGMode::NORMAL);
}

template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph(uint64_t k,
            const std::vector<std::string> &sequences,
            DBGMode mode) {
    auto graph = std::make_shared<Graph>(k, mode == DBGMode::CANONICAL
                                            ? DeBruijnGraph::CANONICAL
                                            : DeBruijnGraph::BASIC);

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    if (mode == DBGMode::PRIMARY) {
        return make_graph_primary<Graph>(graph);
    } else if (mode == DBGMode::CANONICAL_WRAPPER) {
        return std::make_shared<CanonicalDBG>(make_graph_primary<Graph>(graph), 2 /* cache size */);
    }

    return graph;
}

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, DBGMode);

template
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashFast>(uint64_t, const std::vector<std::string> &, DBGMode);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGHashString>(uint64_t k,
                           const std::vector<std::string> &sequences,
                           DBGMode mode) {
    auto graph = std::make_shared<DBGHashString>(k);

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    if (mode == DBGMode::PRIMARY) {
        return make_graph_primary<DBGHashString>(graph);
    } else if (mode == DBGMode::CANONICAL_WRAPPER) {
        return std::make_shared<CanonicalDBG>(make_graph_primary<DBGHashString>(graph), 2 /* cache size */);
    }

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGBitmap>(uint64_t k,
                       const std::vector<std::string> &sequences,
                       DBGMode mode) {
    DBGBitmapConstructor constructor(k, mode == DBGMode::CANONICAL
                                        ? DeBruijnGraph::CANONICAL
                                        : DeBruijnGraph::BASIC);
    for (const auto &sequence : sequences) {
        constructor.add_sequence(std::string(sequence));
    }

    auto graph = std::make_shared<DBGBitmap>(&constructor);

    if (mode == DBGMode::PRIMARY) {
        return make_graph_primary<DBGBitmap>(graph);
    } else if (mode == DBGMode::CANONICAL_WRAPPER) {
        return std::make_shared<CanonicalDBG>(make_graph_primary<DBGBitmap>(graph), 2 /* cache size */);
    }

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinct>(uint64_t k,
                         const std::vector<std::string> &sequences,
                         DBGMode mode) {
    auto graph = std::make_shared<DBGSuccinct>(k, mode == DBGMode::CANONICAL
                                                    ? DeBruijnGraph::CANONICAL
                                                    : DeBruijnGraph::BASIC);

    uint64_t max_index = graph->max_index();

    for (const auto &sequence : sequences) {
        graph->add_sequence(sequence, [&](auto i) { ASSERT_TRUE(i <= ++max_index); });
    }

    [&]() { ASSERT_EQ(max_index, graph->max_index()); }();

    graph->mask_dummy_kmers(1, false);

    if (mode == DBGMode::PRIMARY) {
        return make_graph_primary<DBGSuccinct>(graph);
    } else if (mode == DBGMode::CANONICAL_WRAPPER) {
        return std::make_shared<CanonicalDBG>(make_graph_primary<DBGSuccinct>(graph), 2 /* cache size */);
    }

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<1>>(uint64_t k,
                                   const std::vector<std::string> &sequences,
                                   DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    BOSS *boss;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        DeBruijnGraph &mutable_graph = const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        );
        boss = &dynamic_cast<DBGSuccinct&>(mutable_graph).get_boss();
    } else {
        boss = &dynamic_cast<DBGSuccinct&>(*graph).get_boss();
    }
    boss->index_suffix_ranges(1);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<2>>(uint64_t k,
                                   const std::vector<std::string> &sequences,
                                   DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    BOSS *boss;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        DeBruijnGraph &mutable_graph = const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        );
        boss = &dynamic_cast<DBGSuccinct&>(mutable_graph).get_boss();
    } else {
        boss = &dynamic_cast<DBGSuccinct&>(*graph).get_boss();
    }
    boss->index_suffix_ranges(std::min(k - 1, (uint64_t)2));

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctIndexed<10>>(uint64_t k,
                                    const std::vector<std::string> &sequences,
                                    DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    BOSS *boss;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        DeBruijnGraph &mutable_graph = const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        );
        boss = &dynamic_cast<DBGSuccinct&>(mutable_graph).get_boss();
    } else {
        boss = &dynamic_cast<DBGSuccinct&>(*graph).get_boss();
    }
    boss->index_suffix_ranges(std::min(k - 1, (uint64_t)10));

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloomFPR<1, 1>>(uint64_t k,
                                       const std::vector<std::string> &sequences,
                                       DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter_from_fpr(1.0);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloomFPR<1, 10>>(uint64_t k,
                                        const std::vector<std::string> &sequences,
                                        DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter_from_fpr(1.0 / 10);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloom<4, 1>>(uint64_t k,
                                    const std::vector<std::string> &sequences,
                                    DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter(4.0, 1);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph<DBGSuccinctBloom<4, 50>>(uint64_t k,
                                     const std::vector<std::string> &sequences,
                                     DBGMode mode) {
    auto graph = build_graph<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter(4.0, 50);

    return graph;
}


template <class Graph>
std::shared_ptr<DeBruijnGraph>
build_graph_batch(uint64_t k,
                  const std::vector<std::string> &sequences,
                  DBGMode mode) {
    return build_graph<Graph>(k, sequences, mode);
}

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashOrdered>(uint64_t, const std::vector<std::string> &, DBGMode);

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashFast>(uint64_t, const std::vector<std::string> &, DBGMode);

template
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGHashString>(uint64_t, const std::vector<std::string> &, DBGMode);

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGBitmap>(uint64_t k,
                             const std::vector<std::string> &sequences,
                             DBGMode mode) {
    DBGBitmapConstructor constructor(k, mode == DBGMode::CANONICAL
                                        ? DeBruijnGraph::CANONICAL
                                        : DeBruijnGraph::BASIC);
    constructor.add_sequences(std::vector<std::string>(sequences));
    auto graph = std::make_shared<DBGBitmap>(&constructor);

    if (mode == DBGMode::PRIMARY) {
        return make_graph_primary<DBGBitmap>(graph);
    } else if (mode == DBGMode::CANONICAL_WRAPPER) {
        return std::make_shared<CanonicalDBG>(make_graph_primary<DBGBitmap>(graph), 2 /* cache size */);
    }

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinct>(uint64_t k,
                               const std::vector<std::string> &sequences,
                               DBGMode mode) {
    BOSSConstructor constructor(k - 1, mode == DBGMode::CANONICAL
                                        ? DeBruijnGraph::CANONICAL
                                        : DeBruijnGraph::BASIC);
    EXPECT_EQ(k - 1, constructor.get_k());
    constructor.add_sequences(std::vector<std::string>(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), mode == DBGMode::CANONICAL
                                                                        ? DeBruijnGraph::CANONICAL
                                                                        : DeBruijnGraph::BASIC);
    graph->mask_dummy_kmers(1, false);
    EXPECT_EQ(k, graph->get_k());

    if (mode == DBGMode::PRIMARY) {
        return make_graph_primary<DBGSuccinct>(graph);
    } else if (mode == DBGMode::CANONICAL_WRAPPER) {
        return std::make_shared<CanonicalDBG>(make_graph_primary<DBGSuccinct>(graph), 2 /* cache size */);
    }

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<1>>(uint64_t k,
                                         const std::vector<std::string> &sequences,
                                         DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    BOSS *boss;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        DeBruijnGraph &mutable_graph = const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        );
        boss = &dynamic_cast<DBGSuccinct&>(mutable_graph).get_boss();
    } else {
        boss = &dynamic_cast<DBGSuccinct&>(*graph).get_boss();
    }
    boss->index_suffix_ranges(1);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<2>>(uint64_t k,
                                         const std::vector<std::string> &sequences,
                                         DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    BOSS *boss;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        DeBruijnGraph &mutable_graph = const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        );
        boss = &dynamic_cast<DBGSuccinct&>(mutable_graph).get_boss();
    } else {
        boss = &dynamic_cast<DBGSuccinct&>(*graph).get_boss();
    }
    boss->index_suffix_ranges(std::min(k - 1, (uint64_t)2));

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctIndexed<10>>(uint64_t k,
                                          const std::vector<std::string> &sequences,
                                          DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    BOSS *boss;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        DeBruijnGraph &mutable_graph = const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        );
        boss = &dynamic_cast<DBGSuccinct&>(mutable_graph).get_boss();
    } else {
        boss = &dynamic_cast<DBGSuccinct&>(*graph).get_boss();
    }
    boss->index_suffix_ranges(std::min(k - 1, (uint64_t)10));

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloomFPR<1, 1>>(uint64_t k,
                                             const std::vector<std::string> &sequences,
                                             DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter_from_fpr(1.0);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloomFPR<1, 10>>(uint64_t k,
                                              const std::vector<std::string> &sequences,
                                              DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter_from_fpr(1.0 / 10);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloom<4, 1>>(uint64_t k,
                                          const std::vector<std::string> &sequences,
                                          DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter(4.0, 1);

    return graph;
}

template <>
std::shared_ptr<DeBruijnGraph>
build_graph_batch<DBGSuccinctBloom<4, 50>>(uint64_t k,
                                           const std::vector<std::string> &sequences,
                                           DBGMode mode) {
    auto graph = build_graph_batch<DBGSuccinct>(k, sequences, mode);
    DBGSuccinct *dbg_succ;
    if (mode == DBGMode::CANONICAL_WRAPPER) {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(const_cast<DeBruijnGraph&>(
            std::dynamic_pointer_cast<CanonicalDBG>(graph)->get_graph()
        ));
    } else {
        dbg_succ = &dynamic_cast<DBGSuccinct&>(*graph);
    }
    dbg_succ->initialize_bloom_filter(4.0, 50);

    return graph;
}


template <class Graph>
bool check_graph(const std::string &alphabet, DBGMode mode, bool check_sequence) {
    std::vector<std::string> sequences;

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % alphabet.size()];
        }
        sequences.push_back(seq);
    }

#if _PROTEIN_GRAPH
    auto graph = build_graph<Graph>(12, sequences, mode);
#else
    auto graph = build_graph<Graph>(20, sequences, mode);
#endif

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

template bool check_graph<DBGSuccinct>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctIndexed<1>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctIndexed<2>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctIndexed<10>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctBloomFPR<1, 1>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctBloomFPR<1, 10>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctBloom<4, 1>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGSuccinctBloom<4, 50>>(const std::string &, DBGMode, bool);
template bool check_graph<DBGBitmap>(const std::string &, DBGMode, bool);
template bool check_graph<DBGHashOrdered>(const std::string &, DBGMode, bool);
template bool check_graph<DBGHashFast>(const std::string &, DBGMode, bool);
template bool check_graph<DBGHashString>(const std::string &, DBGMode, bool);

} // namespace test
} // namespace mtg
