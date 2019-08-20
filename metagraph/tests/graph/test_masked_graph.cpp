#include "gtest/gtest.h"

#include "test_dbg_helpers.hpp"

#include <set>

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
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
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
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
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
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACTCGTAGC",
                                                                        "AAATGCGTAGC" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACT",
                                                                        "AAATG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
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
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
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
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
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
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACTCGTAGC",
                                                                        "AAATGCGTAGC" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
        }
        {
            MaskedDeBruijnGraph graph(build_graph_batch<TypeParam>(k, { "AAACT",
                                                                        "AAATG" }));

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_unitigs(callback); }
            ));
            EXPECT_EQ(graph, reconstructed)
                << dynamic_cast<const TypeParam&>(*graph.get_graph_ptr()) << std::endl
                << dynamic_cast<const TypeParam&>(*reconstructed.get_graph_ptr());
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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                ref_nodes.insert(full_graph->get_node_sequence(i));
            }

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
            full_graph->map_to_nodes(
                sequence,
                [&](const auto &index) { mask->set(index, false); }
            );

            MaskedDeBruijnGraph graph(full_graph, mask.release());

            size_t counter = 0;
            graph.map_to_nodes(
                sequence,
                [&](const auto &index) {
                    EXPECT_EQ(DeBruijnGraph::npos, index);
                    ++counter;
                }
            );

            EXPECT_EQ(sequence.size() + 1 - graph.get_k(), counter);

            MaskedDeBruijnGraph reconstructed(build_graph_iterative<TypeParam>(
                k,
                [&](const auto &callback) { graph.call_sequences(callback); }
            ));

            std::multiset<std::string> all_nodes;
            for (DeBruijnGraph::node_index i = 1; i <= full_graph->num_nodes(); ++i) {
                all_nodes.insert(full_graph->get_node_sequence(i));
            }

            std::multiset<std::string> ref_nodes;
            for (DeBruijnGraph::node_index i = 1; i <= graph.num_nodes(); ++i) {
                ref_nodes.insert(graph.get_node_sequence(i));
            }
            EXPECT_EQ(all_nodes, ref_nodes);

            std::multiset<std::string> called_nodes;
            graph.call_nodes([&](const auto &index) {
                called_nodes.insert(graph.get_node_sequence(index));
            });

            std::multiset<std::string> rec_nodes;
            reconstructed.call_nodes([&](const auto &index) {
                rec_nodes.insert(reconstructed.get_node_sequence(index));
            });
            EXPECT_EQ(called_nodes, rec_nodes);
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

            std::multiset<MaskedDeBruijnGraph::node_index> nodes;
            graph.call_nodes([&](const auto &node) { nodes.insert(node); });

            std::multiset<MaskedDeBruijnGraph::node_index> ref_nodes;
            for (size_t i = 1; i <= full_graph->num_nodes(); ++i) {
                if (graph.in_graph(i))
                    ref_nodes.insert(i);
            }

            EXPECT_EQ(ref_nodes, nodes);
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

            graph.call_nodes(
                [&](const auto &node) {
                    std::vector<MaskedDeBruijnGraph::node_index> outnodes;
                    graph.adjacent_outgoing_nodes(node, &outnodes);
                    EXPECT_EQ(outnodes.size(), graph.outdegree(node));
                    std::vector<MaskedDeBruijnGraph::node_index> outnodes_full;
                    full_graph->adjacent_outgoing_nodes(node, &outnodes_full);
                    outnodes_full.erase(std::remove_if(outnodes_full.begin(),
                                                       outnodes_full.end(),
                                                       [&](auto i) {
                                                           return !graph.in_graph(i);
                                                       }),
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

            graph.call_nodes(
                [&](const auto &node) {
                    std::vector<MaskedDeBruijnGraph::node_index> innodes;
                    graph.adjacent_incoming_nodes(node, &innodes);
                    EXPECT_EQ(innodes.size(), graph.indegree(node));
                    std::vector<MaskedDeBruijnGraph::node_index> innodes_full;
                    full_graph->adjacent_incoming_nodes(node, &innodes_full);
                    innodes_full.erase(std::remove_if(innodes_full.begin(),
                                                      innodes_full.end(),
                                                      [&](auto i) {
                                                          return !graph.in_graph(i);
                                                      }),
                                        innodes_full.end());
                    EXPECT_EQ(convert_to_set(innodes_full), convert_to_set(innodes));
                }
            );
        }
    }
}
