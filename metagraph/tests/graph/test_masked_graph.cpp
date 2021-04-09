#include "gtest/gtest.h"

#include <set>

#include "../test_helpers.hpp"
#include "all/test_dbg_helpers.hpp"

#include "graph/representation/masked_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

template <typename Graph>
class MaskedDeBruijnGraphTest : public DeBruijnGraphTest<Graph> {};
typedef ::testing::Types<DBGBitmap,
                         DBGHashString,
                         DBGHashOrdered,
                         DBGHashFast,
                         DBGSuccinct> GraphsToMask;
TYPED_TEST_SUITE(MaskedDeBruijnGraphTest, GraphsToMask);

template <typename Graph>
class MaskedStableDeBruijnGraphTest : public DeBruijnGraphTest<Graph> {};
typedef ::testing::Types<DBGBitmap,
                         DBGSuccinct> StableGraphsToMask;
TYPED_TEST_SUITE(MaskedStableDeBruijnGraphTest, StableGraphsToMask);


TYPED_TEST(MaskedStableDeBruijnGraphTest, CallPathsNoMask) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);
                MaskedDeBruijnGraph masked(graph, [](auto) { return true; });
                // EXPECT_TRUE(check_graph_nodes(masked));

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                    masked.call_sequences([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });
                EXPECT_EQ(*graph, *reconstructed)
                    << dynamic_cast<const TypeParam&>(*graph) << std::endl
                    << dynamic_cast<const TypeParam&>(*reconstructed);
            }
        }
    }
}

TYPED_TEST(MaskedStableDeBruijnGraphTest, CallUnitigsNoMask) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);
                MaskedDeBruijnGraph masked(graph, [](auto) { return true; });
                // EXPECT_TRUE(check_graph_nodes(masked));

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                    masked.call_unitigs([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });
                EXPECT_EQ(*graph, *reconstructed)
                    << dynamic_cast<const TypeParam&>(*graph) << std::endl
                    << dynamic_cast<const TypeParam&>(*reconstructed);
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallPathsNoMask) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);
                MaskedDeBruijnGraph masked(graph, [](auto) { return true; });
                // EXPECT_TRUE(check_graph_nodes(masked));

                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences);

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                    masked.call_sequences([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });
                EXPECT_EQ(*stable_graph, *reconstructed)
                    << dynamic_cast<const TypeParam&>(*graph) << std::endl
                    << dynamic_cast<const TypeParam&>(*reconstructed);
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsNoMask) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);
                MaskedDeBruijnGraph masked(graph, [](auto) { return true; });
                // EXPECT_TRUE(check_graph_nodes(masked));

                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences);

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                    masked.call_unitigs([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });
                EXPECT_EQ(*stable_graph, *reconstructed)
                    << dynamic_cast<const TypeParam&>(*graph) << std::endl
                    << dynamic_cast<const TypeParam&>(*reconstructed);
            }
        }
    }
}

TYPED_TEST(MaskedStableDeBruijnGraphTest, CallPathsMaskFirstKmer) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto full_graph = build_graph_batch<TypeParam>(k, sequences);

                auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
                MaskedDeBruijnGraph graph(
                    full_graph,
                    [&](const auto &index) {
                        return index != DeBruijnGraph::npos
                            && full_graph->get_node_sequence(index) != first_kmer;
                    }
                );
                // EXPECT_TRUE(check_graph_nodes(graph));

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                    graph.call_sequences([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });

                auto ref = build_graph_iterative<TypeParam>(
                    k,
                    [&](const auto &callback) {
                        graph.call_nodes([&](const auto &index) {
                            callback(graph.get_node_sequence(index));
                        });
                    }
                );
                EXPECT_EQ(*ref, *reconstructed)
                    << *full_graph << std::endl
                    << first_kmer << std::endl
                    << graph << std::endl
                    << *ref << std::endl
                    << *reconstructed;

                std::set<std::string> ref_nodes;
                full_graph->call_nodes([&](auto i) {
                    ref_nodes.insert(full_graph->get_node_sequence(i));
                });

                std::set<std::string> rec_nodes;
                reconstructed->call_nodes([&](const auto &index) {
                    rec_nodes.insert(reconstructed->get_node_sequence(index));
                });
                if (ref_nodes.size())
                    rec_nodes.insert(first_kmer);
                EXPECT_EQ(ref_nodes, rec_nodes);
            }
        }
    }
}

TYPED_TEST(MaskedStableDeBruijnGraphTest, CallUnitigsMaskFirstKmer) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto full_graph = build_graph_batch<TypeParam>(k, sequences);

                auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
                MaskedDeBruijnGraph graph(
                    full_graph,
                    [&](const auto &index) {
                        return index != DeBruijnGraph::npos
                            && full_graph->get_node_sequence(index) != first_kmer;
                    }
                );
                // EXPECT_TRUE(check_graph_nodes(graph));

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                    graph.call_unitigs([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });

                auto ref = build_graph_iterative<TypeParam>(
                    k,
                    [&](const auto &callback) {
                        graph.call_nodes([&](const auto &index) {
                            callback(graph.get_node_sequence(index));
                        });
                    }
                );
                EXPECT_EQ(*ref, *reconstructed)
                    << *full_graph << std::endl
                    << first_kmer << std::endl
                    << graph << std::endl
                    << *ref << std::endl
                    << *reconstructed;

                std::set<std::string> ref_nodes;
                full_graph->call_nodes([&](auto i) {
                    ref_nodes.insert(full_graph->get_node_sequence(i));
                });

                std::set<std::string> rec_nodes;
                reconstructed->call_nodes([&](const auto &index) {
                    rec_nodes.insert(reconstructed->get_node_sequence(index));
                });
                if (ref_nodes.size())
                    rec_nodes.insert(first_kmer);
                EXPECT_EQ(ref_nodes, rec_nodes);
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallPathsMaskFirstKmer) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto full_graph = build_graph_batch<TypeParam>(k, sequences);

                auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
                MaskedDeBruijnGraph graph(
                    full_graph,
                    [&](const auto &index) {
                        return index != DeBruijnGraph::npos
                            && full_graph->get_node_sequence(index) != first_kmer;
                    }
                );
                // EXPECT_TRUE(check_graph_nodes(graph));

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                    graph.call_sequences([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });

                auto ref = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph.call_nodes([&](const auto &index) {
                            callback(graph.get_node_sequence(index));
                        });
                    }
                );
                EXPECT_EQ(*ref, *reconstructed)
                    << *full_graph << std::endl
                    << first_kmer << std::endl
                    << graph << std::endl
                    << *ref << std::endl
                    << *reconstructed;

                std::set<std::string> ref_nodes;
                full_graph->call_nodes([&](auto i) {
                    ref_nodes.insert(full_graph->get_node_sequence(i));
                });

                std::set<std::string> rec_nodes;
                reconstructed->call_nodes([&](const auto &index) {
                    rec_nodes.insert(reconstructed->get_node_sequence(index));
                });
                if (ref_nodes.size())
                    rec_nodes.insert(first_kmer);
                EXPECT_EQ(ref_nodes, rec_nodes);
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsMaskFirstKmer) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto full_graph = build_graph_batch<TypeParam>(k, sequences);

                auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
                MaskedDeBruijnGraph graph(
                    full_graph,
                    [&](const auto &index) {
                        return index != DeBruijnGraph::npos
                            && full_graph->get_node_sequence(index) != first_kmer;
                    }
                );
                // EXPECT_TRUE(check_graph_nodes(graph));

                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                    graph.call_unitigs([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });

                auto ref = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph.call_nodes([&](const auto &index) {
                            callback(graph.get_node_sequence(index));
                        });
                    }
                );
                EXPECT_EQ(*ref, *reconstructed)
                    << *full_graph << std::endl
                    << first_kmer << std::endl
                    << graph << std::endl
                    << *ref << std::endl
                    << *reconstructed;

                std::set<std::string> ref_nodes;
                full_graph->call_nodes([&](auto i) {
                    ref_nodes.insert(full_graph->get_node_sequence(i));
                });

                std::set<std::string> rec_nodes;
                reconstructed->call_nodes([&](const auto &index) {
                    rec_nodes.insert(reconstructed->get_node_sequence(index));
                });
                if (ref_nodes.size())
                    rec_nodes.insert(first_kmer);
                EXPECT_EQ(ref_nodes, rec_nodes);
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallContigsMaskPath) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 10; ++k) {
            std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                                 "ATGCAGTACTGAG",
                                                 "GGGGGGGGGGGGG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            for (const auto &sequence : sequences) {
                sdsl::bit_vector mask(full_graph->max_index() + 1, false);

                full_graph->call_nodes([&](auto index) { mask[index] = true; });

                full_graph->map_to_nodes(
                    sequence,
                    [&](const auto &index) { mask[index] = false; }
                );

                MaskedDeBruijnGraph graph(full_graph,
                                          std::make_unique<bit_vector_stat>(std::move(mask)), true);

                EXPECT_EQ(*full_graph, graph.get_base_graph());
                EXPECT_TRUE(check_graph_nodes(graph));

                size_t counter = 0;
                graph.map_to_nodes(
                    sequence,
                    [&](const auto &index) {
                        EXPECT_EQ(DeBruijnGraph::npos, index);
                        ++counter;
                    }
                );

                EXPECT_EQ(sequence.size() + 1 - graph.get_k(), counter);

                // check if reconstructed graph matches
                std::mutex seq_mutex;
                auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                    graph.call_sequences([&](const auto &seq, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                        std::lock_guard<std::mutex> lock(seq_mutex);
                        callback(seq);
                    }, num_threads);
                });

                std::multiset<std::string> called_nodes;
                graph.call_nodes([&](const auto &index) {
                    called_nodes.insert(graph.get_node_sequence(index));
                });

                std::multiset<std::string> rec_nodes;
                reconstructed->call_nodes([&](const auto &index) {
                    rec_nodes.insert(reconstructed->get_node_sequence(index));
                });
                EXPECT_EQ(called_nodes, rec_nodes) << num_threads;
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsMaskPath) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 4;

        /*
        GCATGGT\     GTACT
         ACCGGT\    /
                GGTA
          TAGGT/    \
           GGGT/     GTATT (masked out)
        */
        std::vector<std::string> sequences { "GCATGGTACT",
                                             "ACCGGTACT",
                                             "TAGGTACT",
                                             "GGGTACT",
                                             "GGGGGGGGGGGG",
                                             "GTATT" };
        std::vector<std::string> masked_out { "GTATT", "GGGGGGGGGGGG" };
        std::vector<std::multiset<std::string>> unitig_sets {
            { "GGTACT", "GCATGGT", "ACCGGT", "TAGGT", "GGGT" },
            { "GGTACT", "GCATGGT", "ACCGGT", "TAGGT" },
            { "GGTACT", "GCATGGT", "ACCGGT" },
            { "GGTACT", "GCATGGT" },
        };

        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        sdsl::bit_vector mask(full_graph->max_index() + 1, false);

        for (const auto &sequence : sequences) {
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) { mask[index] = true; }
            );
        }

        for (const auto &sequence_out : masked_out) {
            full_graph->map_to_nodes(
                sequence_out,
                [&](const auto &index) { mask[index] = false; }
            );
        }

        MaskedDeBruijnGraph graph(full_graph,
                                  std::make_unique<bit_vector_stat>(std::move(mask)));

        for (const auto &sequence_out : masked_out) {
            graph.map_to_nodes(
                sequence_out,
                [](const auto &index) { ASSERT_EQ(DeBruijnGraph::npos, index); }
            );
        }

        for (size_t min_tip_size = 1; min_tip_size <= 4; ++min_tip_size) {
            std::multiset<std::string> unitigs;
            std::mutex seq_mutex;
            graph.call_unitigs(
                [&](const auto &unitig, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(graph, unitig));
                    graph.map_to_nodes(
                        unitig,
                        [&](const auto &index) {
                            EXPECT_TRUE(graph.in_subgraph(index));
                            EXPECT_NE(DeBruijnGraph::npos, index);
                        }
                    );
                    std::lock_guard<std::mutex> lock(seq_mutex);
                    unitigs.insert(unitig);
                },
                num_threads,
                min_tip_size
            );
            EXPECT_EQ(unitig_sets[min_tip_size - 1], unitigs) << min_tip_size;
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(MaskedStableDeBruijnGraphTest, CallUnitigsSingleKmerFormCanonical) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 4; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {
                std::string rev = sequences[0];
                reverse_complement(rev.begin(), rev.end());

                auto full_graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
                sdsl::bit_vector mask(full_graph->max_index() + 1, true);
                full_graph->map_to_nodes_sequentially(sequences[0], [&](auto i) {
                    mask[i] = false;
                });
                full_graph->map_to_nodes_sequentially(rev, [&](auto i) {
                    mask[i] = false;
                });
                MaskedDeBruijnGraph graph(full_graph,
                                          std::make_unique<bit_vector_stat>(std::move(mask)));

                graph.map_to_nodes(sequences[0], [&](auto node) {
                    EXPECT_EQ(DeBruijnGraph::npos, node);
                });
                graph.map_to_nodes(rev, [&](auto node) {
                    EXPECT_EQ(DeBruijnGraph::npos, node);
                });
                for (size_t i = 1; i < sequences.size(); ++i) {
                    size_t j = 0;
                    graph.map_to_nodes(sequences[i], [&](auto node) {
                        if (node) {
                            EXPECT_TRUE(sequences[0].find(sequences[i].substr(j, k)) == std::string::npos
                                && rev.find(sequences[i].substr(j, k)) == std::string::npos);
                        } else {
                            EXPECT_TRUE(sequences[0].find(sequences[i].substr(j, k)) != std::string::npos
                                || rev.find(sequences[i].substr(j, k)) != std::string::npos);
                        }
                        ++j;
                    });
                }

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto full_stable_graph = build_graph_batch<DBGSuccinct>(k, sequences, DeBruijnGraph::CANONICAL);
                sdsl::bit_vector stable_mask(full_stable_graph->max_index() + 1, true);
                full_stable_graph->map_to_nodes_sequentially(sequences[0], [&](auto i) {
                    stable_mask[i] = false;
                });
                full_stable_graph->map_to_nodes_sequentially(rev, [&](auto i) {
                    stable_mask[i] = false;
                });

                MaskedDeBruijnGraph stable_masked_graph(
                    full_stable_graph,
                    std::make_unique<bit_vector_stat>(std::move(stable_mask))
                );

                std::mutex seq_mutex;
                auto stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        stable_masked_graph.call_unitigs(
                            [&](const auto &sequence, const auto &path) {
                                ASSERT_EQ(path, map_sequence_to_nodes(stable_masked_graph, sequence));
                                std::unique_lock<std::mutex> lock(seq_mutex);
                                callback(sequence);
                            },
                            num_threads,
                            1, // min_tip_size
                            true // kmers_in_single_form
                        );
                    },
                    DeBruijnGraph::CANONICAL
                );
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph.call_unitigs(
                            [&](const auto &sequence, const auto &path) {
                                ASSERT_EQ(path, map_sequence_to_nodes(graph, sequence));
                                std::unique_lock<std::mutex> lock(seq_mutex);
                                callback(sequence);
                            },
                            num_threads,
                            1, // min_tip_size
                            true // kmers_in_single_form
                        );
                    },
                    DeBruijnGraph::CANONICAL
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}
#endif

TEST(MaskedDBGSuccinct, CallUnitigsMaskLastEdges) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 15; k += 2) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<DBGSuccinct>(k, sequences, DeBruijnGraph::BASIC);
                auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);

                dbg_succ->reset_mask();
                const auto &boss = dbg_succ->get_boss();
                sdsl::bit_vector mask(boss.num_edges() + 1, true);
                size_t num_kmers = mask.size() - 1;
                for (size_t i = 1; i < mask.size(); ++i) {
                    if (boss.get_last(i) || boss.get_node_seq(i)[0] == boss.kSentinelCode) {
                        mask[i] = false;
                        ASSERT_LT(0, num_kmers);
                        --num_kmers;
                    }
                }
                MaskedDeBruijnGraph masked_graph(graph, std::make_unique<bit_vector_stat>(std::move(mask)));
                std::atomic<size_t> counted_kmers(0);
                masked_graph.call_unitigs([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(masked_graph, seq));
                    counted_kmers += path.size();
                }, num_threads);
                EXPECT_EQ(num_kmers, counted_kmers);
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CheckNodes) {
    for (size_t k = 3; k <= 10; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        for (const auto &sequence : sequences) {
            sdsl::bit_vector mask(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask[index] = true; });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask[index] = false;
                }
            );

            MaskedDeBruijnGraph graph(full_graph,
                                      std::make_unique<bit_vector_stat>(std::move(mask)));

            EXPECT_TRUE(check_graph_nodes(graph));

            std::multiset<MaskedDeBruijnGraph::node_index> nodes;
            graph.call_nodes([&](const auto &node) { nodes.insert(node); });

            std::multiset<MaskedDeBruijnGraph::node_index> ref_nodes;
            full_graph->call_nodes([&](auto i) {
                if (graph.in_subgraph(i))
                    ref_nodes.insert(i);
            });

            EXPECT_EQ(ref_nodes, nodes);
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CheckNonExistant) {
    for (size_t k = 3; k <= 10; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        std::string nonexistant = "TTTTTTTTTTTT";
        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        for (const auto &sequence : sequences) {
            sdsl::bit_vector mask(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask[index] = true; });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask[index] = false;
                }
            );

            MaskedDeBruijnGraph graph(full_graph,
                                      std::make_unique<bit_vector_stat>(std::move(mask)));

            EXPECT_TRUE(check_graph_nodes(graph));

            graph.map_to_nodes(nonexistant, [&](auto node) {
                EXPECT_EQ(DeBruijnGraph::npos, node);
            });

            graph.map_to_nodes(sequence, [&](auto node) {
                EXPECT_EQ(DeBruijnGraph::npos, node);
            });
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CheckOutgoingNodes) {
    for (size_t k = 3; k <= 10; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        for (const auto &sequence : sequences) {
            sdsl::bit_vector mask(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask[index] = true; });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask[index] = false;
                }
            );

            MaskedDeBruijnGraph graph(full_graph,
                                      std::make_unique<bit_vector_stat>(std::move(mask)));

            EXPECT_TRUE(check_graph_nodes(graph));

            graph.call_nodes(
                [&](const auto &node) {
                    std::vector<MaskedDeBruijnGraph::node_index> outnodes;
                    graph.adjacent_outgoing_nodes(node, [&](auto i) { outnodes.push_back(i); });
                    EXPECT_EQ(outnodes.size(), graph.outdegree(node));
                    std::vector<MaskedDeBruijnGraph::node_index> outnodes_full;
                    full_graph->adjacent_outgoing_nodes(node, [&](auto i) { outnodes_full.push_back(i); });
                    outnodes_full.erase(std::remove_if(outnodes_full.begin(),
                                                       outnodes_full.end(),
                                                       [&](auto i) { return !graph.in_subgraph(i); }),
                                        outnodes_full.end());
                    EXPECT_EQ(convert_to_set(outnodes_full), convert_to_set(outnodes));
                }
            );
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CheckIncomingNodes) {
    for (size_t k = 3; k <= 10; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        for (const auto &sequence : sequences) {
            sdsl::bit_vector mask(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask[index] = true; });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask[index] = false;
                }
            );

            MaskedDeBruijnGraph graph(full_graph,
                                      std::make_unique<bit_vector_stat>(std::move(mask)));

            EXPECT_TRUE(check_graph_nodes(graph));

            graph.call_nodes(
                [&](const auto &node) {
                    std::vector<MaskedDeBruijnGraph::node_index> innodes;
                    graph.adjacent_incoming_nodes(node, [&](auto i) { innodes.push_back(i); });
                    EXPECT_EQ(innodes.size(), graph.indegree(node));
                    std::vector<MaskedDeBruijnGraph::node_index> innodes_full;
                    full_graph->adjacent_incoming_nodes(node, [&](auto i) { innodes_full.push_back(i); });
                    innodes_full.erase(std::remove_if(innodes_full.begin(),
                                                      innodes_full.end(),
                                                      [&](auto i) { return !graph.in_subgraph(i); }),
                                        innodes_full.end());
                    EXPECT_EQ(convert_to_set(innodes_full), convert_to_set(innodes));
                }
            );
        }
    }
}

} // namespace
