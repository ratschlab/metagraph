#include "gtest/gtest.h"

#define private public
#define protected public

#include <set>

#include "../../test_helpers.hpp"
#include "test_dbg_helpers.hpp"

#include "common/threads/threading.hpp"

namespace {

using namespace mtg;
using namespace mtg::test;

TYPED_TEST_SUITE(DeBruijnGraphTest, GraphTypes);
TYPED_TEST_SUITE(StableDeBruijnGraphTest, StableGraphTypes);


TYPED_TEST(DeBruijnGraphTest, CallPathsEmptyGraph) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            auto empty = build_graph<TypeParam>(k);
            std::vector<std::string> sequences;
            std::mutex seq_mutex;
            empty->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*empty, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                sequences.push_back(sequence);
            }, num_threads);
            ASSERT_EQ(0u, sequences.size());

            EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences));
            EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences));
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsEmptyGraph) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            auto empty = build_graph<TypeParam>(k);
            std::vector<std::string> sequences;
            std::mutex seq_mutex;
            empty->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*empty, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                sequences.push_back(sequence);
            }, num_threads);
            ASSERT_EQ(0u, sequences.size());

            EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences));
            EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences));
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsOneSelfLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= max_test_k<TypeParam>(); ++k) {
            std::vector<std::string> sequences { std::string(2 * k, 'A') };
            auto graph = build_graph<TypeParam>(k, sequences);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
            ASSERT_EQ(1u, graph->num_nodes());
            ASSERT_EQ(1u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsOneSelfLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= max_test_k<TypeParam>(); ++k) {
            std::vector<std::string> sequences { std::string(2 * k, 'A') };
            auto graph = build_graph<TypeParam>(k, sequences);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
            ASSERT_EQ(1u, graph->num_nodes());
            ASSERT_EQ(1u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsThreeSelfLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= max_test_k<TypeParam>(); ++k) {
            std::vector<std::string> sequences { std::string(2 * k, 'A'),
                                                 std::string(2 * k, 'G'),
                                                 std::string(2 * k, 'C') };
            auto graph = build_graph<TypeParam>(k, sequences);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
            ASSERT_EQ(3u, graph->num_nodes());
            ASSERT_EQ(3u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsExtractsLongestOneLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 4; k < std::min((size_t)14, max_test_k<TypeParam>()); ++k) {
            std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                                 "GGGGGGGGGGGGG" };
            auto graph = build_graph<TypeParam>(k, sequences);

            std::vector<std::string> contigs;
            std::mutex seq_mutex;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                contigs.push_back(sequence);
            }, num_threads);

            EXPECT_EQ(2u, contigs.size());
            EXPECT_EQ(convert_to_set({ "ATGCAGTACTCAG", std::string(k, 'G') }),
                      convert_to_set(contigs)) << k;
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsExtractsLongestTwoLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 4; k < std::min((size_t)14, max_test_k<TypeParam>()); ++k) {
            std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                                 "ATGCAGTACTGAG",
                                                 "GGGGGGGGGGGGG" };
            auto graph = build_graph<TypeParam>(k, sequences);

            std::vector<std::string> contigs;
            std::mutex seq_mutex;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                contigs.push_back(sequence);
            }, num_threads);

            EXPECT_EQ(3u, contigs.size());
        }

        // TODO: There are no guarantees on extracted sequence length when
        // executing multithreaded
        break;
    }
}

TYPED_TEST(DeBruijnGraphTest, CallContigsUniqueKmers) {
    for (size_t num_threads : { 1, 4 }) {
        std::string sequence = "GCAAATAAC";
        auto graph = build_graph<TypeParam>(3, { sequence });

        std::atomic<size_t> num_kmers = 0;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
            num_kmers += sequence.size() - 2;
        }, num_threads);

        EXPECT_EQ(sequence.size() - 2, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsUniqueKmersCycle) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 4;
        std::string sequence = "AAACCCGGGTTTAA";
        auto graph = build_graph<TypeParam>(k, { sequence });

        std::atomic<size_t> num_unitigs = 0;
        std::atomic<size_t> num_kmers = 0;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
            num_unitigs++;
            num_kmers += sequence.size() - k + 1;
        }, num_threads);

        EXPECT_EQ(1u, num_unitigs);
        EXPECT_EQ(sequence.size() - k + 1, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallContigsUniqueKmersCycle) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 4;
        std::string sequence = "AAACCCGGGTTTAAA";
        auto graph = build_graph<TypeParam>(k, { sequence });

        std::atomic<size_t> num_contigs = 0;
        std::atomic<size_t> num_kmers = 0;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
            num_contigs++;
            num_kmers += sequence.size() - k + 1;
        }, num_threads);

        EXPECT_EQ(1u, num_contigs);
        EXPECT_EQ(sequence.size() - k + 1, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsFourLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= max_test_k<TypeParam>(); ++k) {
            std::vector<std::string> sequences { std::string(2 * k, 'A'),
                                                 std::string(2 * k, 'G'),
                                                 std::string(2 * k, 'C') };
            auto graph = build_graph<TypeParam>(k, sequences);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
            ASSERT_EQ(3u, graph->num_nodes());
            ASSERT_EQ(3u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_to_nodes_sequentially(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallPaths) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);

                std::vector<std::string> reconst;
                std::mutex seq_mutex;
                graph->call_sequences([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconst.push_back(sequence);
                }, num_threads);
                auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst);

                EXPECT_EQ(*graph, *reconstructed_graph);
            }
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallUnitigs) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);

                std::vector<std::string> reconst;
                std::mutex seq_mutex;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconst.push_back(sequence);
                }, num_threads);
                auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst);

                EXPECT_EQ(*graph, *reconstructed_graph);
            }
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(StableDeBruijnGraphTest, CallPathsFromCanonical) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                std::vector<std::string> reconst;
                std::mutex seq_mutex;
                graph->call_sequences([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconst.push_back(sequence);
                }, num_threads);
                auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, DeBruijnGraph::CANONICAL);

                EXPECT_EQ(*graph, *reconstructed_graph);
            }
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallPathsFromCanonicalSingleKmerForm) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                std::vector<std::string> reconst;
                std::mutex seq_mutex;
                graph->call_sequences([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconst.push_back(sequence);
                }, num_threads, true);
                auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, DeBruijnGraph::CANONICAL);

                EXPECT_EQ(*graph, *reconstructed_graph);
            }
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallUnitigsFromCanonical) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                std::vector<std::string> reconst;
                std::mutex seq_mutex;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconst.push_back(sequence);
                }, num_threads);
                auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, DeBruijnGraph::CANONICAL);

                EXPECT_EQ(*graph, *reconstructed_graph);
            }
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallUnitigsFromCanonicalSingleKmerForm) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                std::vector<std::string> reconst;
                std::mutex seq_mutex;
                // TODO: why is min_tip_size == 0?
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconst.push_back(sequence);
                }, num_threads, 0, true);
                auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, DeBruijnGraph::CANONICAL);

                EXPECT_EQ(*graph, *reconstructed_graph);
            }
        }
    }
}
#endif

TYPED_TEST(DeBruijnGraphTest, CallPaths) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences);

                std::mutex seq_mutex;
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph->call_sequences([&](const auto &sequence, const auto &path) {
                            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                            std::unique_lock<std::mutex> lock(seq_mutex);
                            callback(sequence);
                        }, num_threads);
                    }
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigs) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences);

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences);

                std::mutex seq_mutex;
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph->call_unitigs([&](const auto &sequence, const auto &path) {
                            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                            std::unique_lock<std::mutex> lock(seq_mutex);
                            callback(sequence);
                        }, num_threads, 1);
                    }
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsWithoutTips) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 3;
        auto graph = build_graph<TypeParam>(k, { "ACTAAGC",
                                                 "TCTAAGC" });
        ASSERT_EQ(6u, graph->num_nodes());

        std::set<std::string> unitigs;
        std::mutex seq_mutex;
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "ACT", "TCT", "CTAAGC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "ACT", "TCT", "CTAAGC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "CTAAGC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "CTAAGC" }), unitigs);


        graph = build_graph<TypeParam>(k, { "ACTAAGC",
                                            "ACTAAGT" });
        ASSERT_EQ(6u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "ACTAAG", "AGC", "AGT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "ACTAAG", "AGC", "AGT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "ACTAAG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "ACTAAG" }), unitigs);


        graph = build_graph<TypeParam>(k, { "ACTAAGCCC",
                                            "AAAGC",
                                            "TAAGCA" });
        ASSERT_EQ(9u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCCC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCC", "CCC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCCC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCC", "CCC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAGC", "GCC", "CCC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCC", "CCC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 3);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAGC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCC", "CCC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        // EXPECT_EQ(std::set<std::string>({ "AAGC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "AAGC", "AAA", "ACTAA", "GCC", "CCC" }), unitigs);


        graph = build_graph<TypeParam>(k, { "ACGAAGCCT",
                                            "AAGC",
                                            "TAAGCA" });
        ASSERT_EQ(9u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "ACGAA", "TAA", "AAGC", "GCA", "GCCT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "ACGAA", "TAA", "AAGC", "GCA", "GCCT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "ACGAA", "AAGC", "GCCT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 3);
        EXPECT_EQ(std::set<std::string>({ "ACGAA", "AAGC" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "AAGC" }), unitigs);


        graph = build_graph<TypeParam>(k, { "TCTAAGCCG",
                                            "CATAAGCCG",
                                            "CATAACCGA" });
        ASSERT_EQ(12u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "CATA", "TCTA", "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "CATA", "TCTA", "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "CATA", "TCTA", "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 3);
        EXPECT_EQ(std::set<std::string>({ "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsWithoutTips2) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 5;
        auto graph = build_graph<TypeParam>(k, { "ACTATAGCTAGTCTATGCGA",
                                                 "ACTATAGCTAGTCTAA",
                                                 "ACTATAGCTA",
                                                 "ACTATAGCTT",
                                                 "ACTATC", });
        ASSERT_EQ(19u, graph->num_nodes());
        std::set<std::string> unitigs;
        std::mutex seq_mutex;
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTT", "AGCTAGTCTA", "TCTAA", "TCTAT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTT", "AGCTAGTCTA", "TCTAA", "TCTAT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT" }), unitigs);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsCheckDegree) {
    for (size_t num_threads : { 1, 4 }) {
        std::vector<std::string> sequences {
            "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGG",
            "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGC",
            "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATAACTGGAAAA",
            "CCGAACTAGTGAAACTGCAACAGACATACGCTGCTCTGAACTCTAAGGC",
            "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
            "ATCGGAAGAGCACACGTCTGAACTCCAGACACTAAGGCATCTCGTATGC",
            "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAATTGGCTGAAAAGC",
            "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAATCGCTCCACCAA"
        };

        auto graph = build_graph_batch<TypeParam>(9, sequences);

        std::multiset<std::string> unitigs {
            "AGACAAATCGCTCCACCAA",
            "AGACAAATTGGCTGAAAAGC",
            "ATCGGAAGAGCACACGTCTGAACT",
            "CAGACATAACTGGAAAA",
            "CAGACATACGCTGCTCTGAACT",
            "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATA",
            "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAAT",
            "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCG",
            "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
            "CCGAACTAGTGAAACTGCAACAGACATA",
            "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAAT",
            "CTGAACTCCAGACACTAAGGCATCTCGTATGC",
            "CTGAACTCTAAGGC",
            "TCTGAACTC"
        };

        std::multiset<std::string> obs_unitigs;
        std::mutex seq_mutex;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
            std::unique_lock<std::mutex> lock(seq_mutex);
            obs_unitigs.insert(sequence);
        }, num_threads, 2);

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(DeBruijnGraphTest, CallUnitigsIndegreeFirstNodeIsZero) {
    for (size_t num_threads : { 1, 4 }) {
        std::vector<std::string> sequences {
            "AGAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
            "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
            "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAATTAG"
        };

        auto graph = build_graph_batch<TypeParam>(31, sequences);

        std::multiset<std::string> unitigs {
            "GAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
            "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
            "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAAT"
        };

        std::multiset<std::string> obs_unitigs;
        std::mutex seq_mutex;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
            std::unique_lock<std::mutex> lock(seq_mutex);
            obs_unitigs.insert(sequence);
        }, num_threads, 2);

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}
#endif

TYPED_TEST(DeBruijnGraphTest, CallUnitigsCross) {
    for (size_t num_threads : { 1, 4 }) {
        // AATTT - ATTTT           TTTAA - TTAAA
        //               > TTTTA <
        // GGTTT - GTTTT           TTTAG - TTAGG

        // build graph from k-mers added in different order
        for (const auto &sequences : {
            std::vector<std::string>({ "AATTTTAAA",
                                       "GGTTTTAGG", }),
            std::vector<std::string>({ "GGTTTTAGG",
                                       "AATTTTAAA", }),
            std::vector<std::string>({ "TTTTAAA",
                                       "TTTTAGG",
                                       "AATTTTA",
                                       "GGTTTTA", }),
            std::vector<std::string>({ "AATTTTA",
                                       "GGTTTTA",
                                       "TTTTAAA",
                                       "TTTTAGG", }) }) {
            auto graph = build_graph_batch<TypeParam>(5, sequences);

            std::multiset<std::string> unitigs {
                "AATTTT",
                "GGTTTT",
                "TTTTA",
                "TTTAAA",
                "TTTAGG",
            };

            std::mutex seq_mutex;

            for (size_t t = 0; t <= 2; ++t) {
                std::multiset<std::string> obs_unitigs;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    obs_unitigs.insert(sequence);
                }, num_threads, t);
                EXPECT_EQ(unitigs, obs_unitigs) << t;
            }

            std::multiset<std::string> long_unitigs {
                "TTTTA",
            };

            for (size_t t = 3; t <= 10; ++t) {
                std::multiset<std::string> obs_long_unitigs;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_to_nodes_sequentially(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    obs_long_unitigs.insert(sequence);
                }, num_threads, 3);
                EXPECT_EQ(long_unitigs, obs_long_unitigs);
            }
        }
    }
}

} // namespace
