#include "gtest/gtest.h"

#include "test_dbg_helpers.hpp"

#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "masked_graph.hpp"


template <typename Graph>
class MaskedDeBruijnGraphTest : public DeBruijnGraphTest<Graph> { };
TYPED_TEST_CASE(MaskedDeBruijnGraphTest, GraphTypes);

TYPED_TEST(MaskedDeBruijnGraphTest, CallPathsNoMask) {
    for (size_t k = 3; k <= 10; ++k) {
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACACTAG",
                                                                        "AACGACATG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AGACACTGA"
                                                                        "GACTACGTA"
                                                                        "ACTAACGTA" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AGACACAGT"
                                                                        "GACTTGCAG"
                                                                        "ACTAGTCAG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACTCGTAGC",
                                                                        "AAATGCGTAGC" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACT",
                                                                        "AAATG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsNoMask) {
    for (size_t k = 3; k <= 10; ++k) {
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACACTAG",
                                                                        "AACGACATG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AGACACTGA"
                                                                        "GACTACGTA"
                                                                        "ACTAACGTA" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AGACACAGT"
                                                                        "GACTTGCAG"
                                                                        "ACTAGTCAG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACTCGTAGC",
                                                                        "AAATGCGTAGC" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACT",
                                                                        "AAATG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << *dynamic_cast<TypeParam*>(graph.get_graph_ptr().get()) << std::endl
                << *dynamic_cast<TypeParam*>(reconstructed.get_graph_ptr().get());
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallPathsMaskFirstKmer) {
    for (size_t k = 3; k <= 10; ++k) {
        {
            std::vector<std::string> sequences{ "AAACACTAG",
                                                "AACGACATG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            const std::vector<std::string> sequences{ "AGACACTGA",
                                                      "GACTACGTA",
                                                      "ACTAACGTA" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            std::vector<std::string> sequences{ "AGACACAGT",
                                                "GACTTGCAG",
                                                "ACTAGTCAG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            std::vector<std::string> sequences{ "AAACTCGTAGC",
                                                "AAATGCGTAGC" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            std::vector<std::string> sequences{ "AAACT",
                                                "AAATG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsMaskFirstKmer) {
    for (size_t k = 3; k <= 10; ++k) {
        {
            std::vector<std::string> sequences{ "AAACACTAG",
                                                "AACGACATG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            const std::vector<std::string> sequences{ "AGACACTGA",
                                                      "GACTACGTA",
                                                      "ACTAACGTA" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            std::vector<std::string> sequences{ "AGACACAGT",
                                                "GACTTGCAG",
                                                "ACTAGTCAG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            std::vector<std::string> sequences{ "AAACTCGTAGC",
                                                "AAATGCGTAGC" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            if (ref_nodes.size())
                rec_nodes.insert(first_kmer);
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
        {
            std::vector<std::string> sequences{ "AAACT",
                                                "AAATG" };
            auto full_graph = build_graph_batch<TypeParam>(k, sequences);

            auto first_kmer = sequences[0].substr(0, std::min(k, sequences[0].length()));
            MaskedDeBruijnGraph graph(
                full_graph,
                [&](const auto &index) {
                    return index != DeBruijnGraph::npos
                        && full_graph->get_node_sequence(index) != first_kmer;
                }
            );
            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));

            MaskedDeBruijnGraph ref(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) {
                    graph.call_nodes([&](const auto &index) {
                        callback(graph.get_node_sequence(index));
                    });
                }
            ));
            EXPECT_EQ(ref, reconstructed)
                << *full_graph << std::endl
                << first_kmer << std::endl
                << graph << std::endl
                << ref << std::endl
                << reconstructed;

            std::set<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::set<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
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
            std::unique_ptr<bit_vector> mask(
                new bit_vector_stat(full_graph->num_nodes() + 1, true)
            );
            mask->set(DeBruijnGraph::npos, false);
            std::set<std::string> erased;
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) {
                    erased.insert(full_graph->get_node_sequence(index));
                    mask->set(index, false);
                }
            );

            MaskedDeBruijnGraph graph(full_graph, mask.release());

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            std::multiset<std::string> ref_nodes;
            full_graph->call_nodes([&](const auto &index) {
                ref_nodes.insert(full_graph->get_node_sequence(index));
            });

            std::multiset<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            rec_nodes.insert(erased.begin(), erased.end());
            EXPECT_EQ(ref_nodes, rec_nodes);
        }
    }
}
