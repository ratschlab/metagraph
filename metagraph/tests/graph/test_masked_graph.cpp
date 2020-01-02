#include "gtest/gtest.h"

#include <set>

#include "../test_helpers.hpp"
#include "all/test_dbg_helpers.hpp"

#include "masked_graph.hpp"


template <typename Graph>
class MaskedDeBruijnGraphTest : public DeBruijnGraphTest<Graph> { };
TYPED_TEST_CASE(MaskedDeBruijnGraphTest, GraphTypes);

template <typename Graph>
class MaskedStableDeBruijnGraphTest : public DeBruijnGraphTest<Graph> { };
TYPED_TEST_CASE(MaskedStableDeBruijnGraphTest, StableGraphTypes);


TYPED_TEST(MaskedStableDeBruijnGraphTest, CallPathsNoMask) {
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

            auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                masked.call_sequences([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                    callback(seq);
                });
            });
            EXPECT_EQ(*graph, *reconstructed)
                << dynamic_cast<const TypeParam&>(*graph) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed);
        }
    }
}

TYPED_TEST(MaskedStableDeBruijnGraphTest, CallUnitigsNoMask) {
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

            auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                masked.call_unitigs([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                    callback(seq);
                });
            });
            EXPECT_EQ(*graph, *reconstructed)
                << dynamic_cast<const TypeParam&>(*graph) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed);
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallPathsNoMask) {
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

            auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                masked.call_sequences([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                    callback(seq);
                });
            });
            EXPECT_EQ(*stable_graph, *reconstructed)
                << dynamic_cast<const TypeParam&>(*graph) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed);
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsNoMask) {
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

            auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                masked.call_unitigs([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(masked, seq));
                    callback(seq);
                });
            });
            EXPECT_EQ(*stable_graph, *reconstructed)
                << dynamic_cast<const TypeParam&>(*graph) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed);
        }
    }
}

TYPED_TEST(MaskedStableDeBruijnGraphTest, CallPathsMaskFirstKmer) {
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

            auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                graph.call_sequences([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                    callback(seq);
                });
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

TYPED_TEST(MaskedStableDeBruijnGraphTest, CallUnitigsMaskFirstKmer) {
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

            auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                graph.call_unitigs([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                    callback(seq);
                });
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

TYPED_TEST(MaskedDeBruijnGraphTest, CallPathsMaskFirstKmer) {
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

            auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                graph.call_sequences([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                    callback(seq);
                });
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

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsMaskFirstKmer) {
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

            auto reconstructed = build_graph_iterative<DBGSuccinct>(k, [&](const auto &callback) {
                graph.call_unitigs([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                    callback(seq);
                });
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

TYPED_TEST(MaskedDeBruijnGraphTest, CallContigsMaskPath) {
    for (size_t k = 3; k <= 10; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        for (const auto &sequence : sequences) {
            auto mask = std::make_unique<bit_vector_stat>(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask->set(index, true); });

            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) { mask->set(index, false); }
            );

            MaskedDeBruijnGraph graph(full_graph, std::move(mask), true);
            EXPECT_EQ(*full_graph, graph.get_graph());
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
            auto reconstructed = build_graph_iterative<TypeParam>(k, [&](const auto &callback) {
                graph.call_sequences([&](const auto &seq, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(graph, seq));
                    callback(seq);
                });
            });

            std::multiset<std::string> called_nodes;
            graph.call_nodes([&](const auto &index) {
                called_nodes.insert(graph.get_node_sequence(index));
            });

            std::multiset<std::string> rec_nodes;
            reconstructed->call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed->get_node_sequence(index));
            });
            EXPECT_EQ(called_nodes, rec_nodes);
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsMaskPath) {
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

    auto mask = std::make_unique<bit_vector_stat>(
        full_graph->max_index() + 1, false
    );

    for (const auto &sequence : sequences) {
        full_graph->map_to_nodes(
            sequence,
            [&](const auto &index) { mask->set(index, true); }
        );
    }

    for (const auto &sequence_out : masked_out) {
        full_graph->map_to_nodes(
            sequence_out,
            [&](const auto &index) { mask->set(index, false); }
        );
    }

    MaskedDeBruijnGraph graph(full_graph, std::move(mask));

    for (const auto &sequence_out : masked_out) {
        graph.map_to_nodes(
            sequence_out,
            [](const auto &index) { ASSERT_EQ(DeBruijnGraph::npos, index); }
        );
    }

    for (size_t min_tip_size = 1; min_tip_size <= 4; ++min_tip_size) {
        std::multiset<std::string> unitigs;
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
                unitigs.insert(unitig);
            },
            min_tip_size
        );
        EXPECT_EQ(unitig_sets[min_tip_size - 1], unitigs) << min_tip_size;
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CheckNodes) {
    for (size_t k = 3; k <= 10; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        auto full_graph = build_graph_batch<TypeParam>(k, sequences);

        for (const auto &sequence : sequences) {
            auto mask = std::make_unique<bit_vector_stat>(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask->set(index, true); });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask->set(index, false);
                }
            );

            MaskedDeBruijnGraph graph(full_graph, std::move(mask));
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
            auto mask = std::make_unique<bit_vector_stat>(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask->set(index, true); });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask->set(index, false);
                }
            );

            MaskedDeBruijnGraph graph(full_graph, std::move(mask));
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
            auto mask = std::make_unique<bit_vector_stat>(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask->set(index, true); });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask->set(index, false);
                }
            );

            MaskedDeBruijnGraph graph(full_graph, std::move(mask));
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
            auto mask = std::make_unique<bit_vector_stat>(full_graph->max_index() + 1, false);

            full_graph->call_nodes([&](auto index) { mask->set(index, true); });

            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask->set(index, false);
                }
            );

            MaskedDeBruijnGraph graph(full_graph, std::move(mask));
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
