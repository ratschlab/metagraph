#include "gtest/gtest.h"

#include "../graph/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "threading.hpp"
#include "annotated_graph_algorithm.hpp"


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_call_significant_indices(double density_cutoff,
                                   double outlabel_mixture) {
    const std::vector<std::string> ingroup { "B", "C" };
    const std::vector<std::string> outgroup { "A" };

    for (size_t k = 3; k < 15; ++k) {
        const std::vector<std::string> sequences {
            std::string("T") + std::string(k - 1, 'A') + std::string(100, 'T'),
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "G"
        };
        const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        const std::unordered_set<std::string> ref { "B", "C", "D" };

        auto masked_dbg = build_masked_graph(*anno_graph,
                                             ingroup,
                                             outgroup,
                                             density_cutoff,
                                             outlabel_mixture);
        EXPECT_EQ(anno_graph->get_graph().num_nodes(), masked_dbg.num_nodes());

        masked_dbg.call_nodes(
            [&](const auto &index) {
                auto cur_labels = anno_graph->get_labels(index);

                obs_labels.insert(cur_labels.begin(), cur_labels.end());
            }
        );

        EXPECT_EQ(ref, obs_labels) << k << " " << density_cutoff << " " << outlabel_mixture;
    }
}

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_call_significant_indices_lazy(double density_cutoff,
                                        double outlabel_mixture) {
    const std::vector<std::string> ingroup { "B", "C" };
    const std::vector<std::string> outgroup { "A" };

    for (size_t k = 3; k < 15; ++k) {
        const std::vector<std::string> sequences {
            std::string("T") + std::string(k - 1, 'A') + std::string(100, 'T'),
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "G"
        };
        const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        const std::unordered_set<std::string> ref { "B", "C", "D" };

        auto masked_dbg = build_masked_graph_lazy(*anno_graph,
                                                  ingroup,
                                                  outgroup,
                                                  density_cutoff,
                                                  outlabel_mixture);
        EXPECT_EQ(anno_graph->get_graph().num_nodes(), masked_dbg.num_nodes());

        masked_dbg.call_nodes(
            [&](const auto &index) {
                auto cur_labels = anno_graph->get_labels(index);

                obs_labels.insert(cur_labels.begin(), cur_labels.end());
            }
        );

        EXPECT_EQ(ref, obs_labels) << k << " " << density_cutoff << " " << outlabel_mixture;
    }
}



template <typename GraphAnnotationPair>
class MaskedDeBruijnGraphAlgorithm : public ::testing::Test {};

TYPED_TEST_CASE(MaskedDeBruijnGraphAlgorithm, GraphAnnotationPairTypes);

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, CallSignificantIndices) {
    for (double d = 0.0; d <= 1.0; d += 0.05) {
        test_call_significant_indices<typename TypeParam::first_type,
                                      typename TypeParam::second_type>(0.0, d);
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, CallSignificantIndicesLazy) {
    for (double d = 0.0; d <= 1.0; d += 0.05) {
        test_call_significant_indices_lazy<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(0.0, d);
    }
}

bool all_mapped_match_first(const SequenceGraph &graph,
                            const std::string &sequence,
                            const DeBruijnGraph::node_index &index) {
    std::vector<SequenceGraph::node_index> indices;
    graph.map_to_nodes(sequence, [&](const auto &i) { indices.push_back(i); });

    EXPECT_FALSE(indices.empty());
    if (indices.empty()) {
        ADD_FAILURE();
        return false;
    }

    EXPECT_EQ(index, indices.front());
    return std::all_of(indices.begin(), indices.end(),
                       [&](const auto &i) { return i != DeBruijnGraph::npos; });
}


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_find_breakpoints(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };

    for (size_t k = 2; k < 7; ++k) {
        const std::vector<std::string> sequences {
            k == 2 ? "ATGC" : "ATGCAGCTTG",
            k == 2 ? "ATGA" : "ATGCAGCTTA"
        };
        const std::vector<std::string> labels { "A", "B" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        ThreadPool thread_pool(pool_size);

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);

        std::atomic<size_t> counter = 0;
        annotated_graph_algorithm::call_breakpoints(
            masked_dbg,
            *anno_graph,
            [&](auto&& path, const std::string &ref, auto&& vlabels) {
                check_json_dump_load(masked_dbg.get_graph(), path, ref);

                auto index = path.front();
                const auto& var = path.get_sequence();

                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), var, index))
                    << k << " " << ref << " " << var;
                EXPECT_FALSE(all_mapped_match_first(masked_dbg, var, index))
                    << k << " " << ref << " " << var;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k;
                for (const auto &label : vlabels) {
                    EXPECT_EQ(std::string("B"), label) << k;
                }
                ++counter;
            },
            &thread_pool
        );

        EXPECT_NE(0u, counter);
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBreakpoints) {
    test_find_breakpoints<typename TypeParam::first_type,
                          typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBreakpointsParallel) {
    test_find_breakpoints<typename TypeParam::first_type,
                          typename TypeParam::second_type>(3);
}


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_find_bubbles(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };

    for (size_t k = 3; k < 7; ++k) {
        std::unordered_map<char, char> complement { { 'A', 'T' },
                                                    { 'T', 'A' },
                                                    { 'C', 'G' },
                                                    { 'G', 'C' } };
        std::vector<std::string> sequences {
            "ATGCAGTACTCAG",
            "ATGCAGTACTCAG"
        };
        const std::vector<std::string> labels { "A", "B" };
        sequences[1][k] = complement[sequences[1][k]];

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        std::mutex add_mutex;
        ThreadPool thread_pool(pool_size);
        const std::unordered_set<std::string> ref { "B" };

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);

        annotated_graph_algorithm::call_bubbles(
            masked_dbg,
            *anno_graph,
            [&](auto&& path, const std::string &ref, auto&& vlabels) {
                check_json_dump_load(masked_dbg.get_graph(), path, ref);

                auto index = path.front();
                const auto& var = path.get_sequence();

                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), var, index))
                    << k << " " << ref << " " << var;
                EXPECT_FALSE(all_mapped_match_first(masked_dbg, var, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var.substr(0, k), index));
                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k;
                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
            },
            &thread_pool
        );

        EXPECT_EQ(ref, obs_labels) << k;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubbles) {
    test_find_bubbles<typename TypeParam::first_type,
                      typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesParallel) {
    test_find_bubbles<typename TypeParam::first_type,
                      typename TypeParam::second_type>(3);
}


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_find_bubbles_incomplete(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };

    for (size_t k = 3; k < 7; ++k) {
        std::unordered_map<char, char> complement { { 'A', 'T' },
                                                    { 'T', 'A' },
                                                    { 'C', 'G' },
                                                    { 'G', 'C' } };
        std::vector<std::string> sequences {
            "AATTACCGGCAGT",
            "AATTACCGGCAGT"
        };
        const std::vector<std::string> labels { "A", "B" };
        auto it = sequences[1].rbegin() + k - 1;
        *it = complement[*it];

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        std::mutex add_mutex;
        ThreadPool thread_pool(pool_size);
        const std::unordered_set<std::string> ref { "B" };

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);

        size_t counter = 0;
        annotated_graph_algorithm::call_bubbles(
            masked_dbg,
            *anno_graph,
            [&](auto&& path, const std::string &ref, auto&& vlabels) {
                check_json_dump_load(masked_dbg.get_graph(), path, ref);

                auto index = path.front();
                const auto& var = path.get_sequence();

                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), var, index))
                    << k << " " << ref << " " << var;
                EXPECT_FALSE(all_mapped_match_first(masked_dbg, var, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var.substr(0, k), index));
                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k;
                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
                ++counter;
            },
            &thread_pool
        );

        EXPECT_EQ(std::unordered_set<std::string>{}, obs_labels) << k;
        EXPECT_EQ(0u, counter);
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesIncomplete) {
    test_find_bubbles_incomplete<typename TypeParam::first_type,
                                 typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesIncompleteParallel) {
    test_find_bubbles_incomplete<typename TypeParam::first_type,
                                 typename TypeParam::second_type>(3);
}


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_find_bubbles_inner_loop(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A", "B" };
    const std::vector<std::string> outgroup { "D" };

    for (size_t k = 3; k < 4; ++k) {
        std::vector<std::string> sequences {
            "ATGATG",
            "ATGCTATG",
            "ATGTATG",
            "TGT"
        };
        const std::vector<std::string> labels { "A", "B", "C", "D" };
        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        std::mutex add_mutex;
        ThreadPool thread_pool(pool_size);
        const std::unordered_set<std::string> ref { "C" };

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);

        size_t counter = 0;
        annotated_graph_algorithm::call_bubbles(
            masked_dbg,
            *anno_graph,
            [&](auto&& path, const std::string &ref, auto&& vlabels) {
                check_json_dump_load(masked_dbg.get_graph(), path, ref);

                auto index = path.front();
                const auto& var = path.get_sequence();

                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), var, index))
                    << k << " " << ref << " " << var;
                EXPECT_FALSE(all_mapped_match_first(masked_dbg, var, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var.substr(0, k), index));
                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k;
                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
                ++counter;
            },
            &thread_pool
        );

        EXPECT_EQ(std::unordered_set<std::string>{}, obs_labels) << k;
        EXPECT_EQ(0u, counter);
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesInnerLoop) {
    test_find_bubbles_inner_loop<typename TypeParam::first_type,
                                 typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesInnerLoopParallel) {
    test_find_bubbles_inner_loop<typename TypeParam::first_type,
                                 typename TypeParam::second_type>(3);
}

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_align_to_masked_graph(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };

    size_t k = 4;

    // TTGC      GCACGGGTC
    //      TGCA
    // ATGC      GCAGTGGTC
    std::vector<std::string> sequences { "TTGCACGGGTC", "ATGCAGTGGTC" };
    const std::vector<std::string> labels { "A", "B" };

    auto config = DBGAlignerConfig(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    std::vector<std::pair<std::string,
                          std::vector<std::pair<std::string,
                                                score_t>>>> query_best_paths {
        { sequences[0],  { { sequences[0],  22 }, { sequences[0],  22 } } },
        { sequences[1],  { { sequences[1],  22 }, { sequences[1],  22 } } },
        { "TTGCAGTGGTC", { { "TTGCAGTGGTC", 22 }, { sequences[0],  14 } } },
        { "ATGCACGGGTC", { { "ATGCACGGGTC", 22 }, { sequences[1],  14 } } },
        { "TTGCAGTGGT",  { { "TTGCAGTGGT",  20 }, { "TTGCACGGGT",  12 } } },
        { "ATGCACGGGT",  { { "ATGCACGGGT",  20 }, { "ATGCAGTGGT",  12 } } },
        { "TTGCAGTGG",   { { "TTGCAGTGG",   18 }, { "TTGCA",       10 } } },
        { "ATGCACGGG",   { { "ATGCACGGG",   18 }, { "ATGCA",       10 } } },
        { "TTGCAGTG",    { { "TTGCAGTG",    16 }, { "TTGCA",       10 } } },
        { "ATGCACGG",    { { "ATGCACGG",    16 }, { "ATGCA",       10 } } },
        { "TTGCAGT",     { { "TTGCAGT",     14 }, { "TTGCA",       10 } } },
        { "ATGCACG",     { { "ATGCACG",     14 }, { "ATGCA",       10 } } },
        { "TTGCAG",      { { "TTGCAG",      12 }, { "TTGCA",       10 } } },
        { "ATGCAC",      { { "ATGCAC",      12 }, { "ATGCA",       10 } } },
        { "TTGCA",       { { "TTGCA",       10 }, { "TTGCA",       10 } } },
        { "ATGCA",       { { "ATGCA",       10 }, { "ATGCA",       10 } } }
    };

    auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);
    const auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(
        anno_graph->get_graph_ptr()
    );
    ASSERT_TRUE(graph.get());

    std::unordered_set<std::string> obs_labels;
    std::mutex add_mutex;
    ThreadPool thread_pool(pool_size);
    const std::unordered_set<std::string> ref { "B" };

    Cigar::initialize_opt_table(graph->alphabet());

    std::vector<DBGAligner<>> aligners {
        DBGAligner<>(*graph, config),
        DBGAligner<>(*graph,
                     config,
                     suffix_seeder<DeBruijnGraph::node_index>,
                     annotated_graph_algorithm::build_masked_graph_extender(*anno_graph))
    };

    size_t j = 0;
    for (const auto& [query, paths] : query_best_paths) {
        ASSERT_EQ(paths.size(), aligners.size());

        for (size_t i = 0; i < paths.size(); ++i) {
            for (auto&& alignment : aligners[i].align(query)) {
                EXPECT_EQ(paths[i].first, alignment.get_sequence()) << j << " " << i;
                EXPECT_EQ(paths[i].second, alignment.get_score()) << j << " " << i;
            }
        }

        ++j;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, AlignToMaskedGraph) {
    test_align_to_masked_graph<typename TypeParam::first_type,
                               typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, AlignToMaskedGraphParallel) {
    test_align_to_masked_graph<typename TypeParam::first_type,
                               typename TypeParam::second_type>(3);
}


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_find_variants(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };

    for (size_t k = 5; k < 11; ++k) {
        std::unordered_map<char, char> complement { { 'A', 'T' },
                                                    { 'T', 'A' },
                                                    { 'C', 'G' },
                                                    { 'G', 'C' } };
        std::vector<std::string> sequences {
            "ATGCAGTACTCAG",
            "ATGCAGTACTCAG"
        };
        const std::vector<std::string> labels { "A", "B" };
        sequences[1][k] = complement[sequences[1][k]];

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        std::mutex add_mutex;
        ThreadPool thread_pool(pool_size);
        const std::unordered_set<std::string> ref { "B" };

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);

        Cigar::initialize_opt_table(masked_dbg.get_graph().alphabet());

        auto score_matrix = DBGAlignerConfig::unit_scoring_matrix(
            1, masked_dbg.get_graph().alphabet()
        );

        annotated_graph_algorithm::call_variants(
            masked_dbg,
            *anno_graph,
            [&](auto&& path, const std::string &query_sequence, auto&& vlabels) {
                check_json_dump_load(masked_dbg.get_graph(), path, query_sequence);

                auto ref = std::string(path.get_query_begin(), path.get_query_end());
                const auto& var = path.get_sequence();

                ASSERT_EQ(query_sequence.c_str() + path.get_clipping(),
                          path.get_query_begin());

                auto index = path.front();
                auto cigar = path.get_cigar().to_string();

                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), var, index))
                    << k << " " << ref << " " << var;
                EXPECT_FALSE(all_mapped_match_first(masked_dbg, var, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var.substr(0, k), index));

                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k << " " << ref;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k << " " << var;

                ASSERT_FALSE(cigar.empty());

                EXPECT_EQ(std::string(std::to_string(k) + "=1X"),
                          std::string(cigar.begin(), cigar.end() - 2))
                    << k << " " << query_sequence << " " << ref << " " << var << " " << cigar;
                EXPECT_EQ('=', cigar.back())
                    << k << " " << query_sequence << " " << ref << " " << var << " " << cigar;

                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
            },
            DBGAlignerConfig(score_matrix),
            default_extender<DeBruijnGraph::node_index>,
            &thread_pool
        );

        EXPECT_EQ(ref, obs_labels) << k;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindVariants) {
    test_find_variants<typename TypeParam::first_type,
                       typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindVariantsParallel) {
    test_find_variants<typename TypeParam::first_type,
                       typename TypeParam::second_type>(3);
}

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_breakpoint_seeder() {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };
    size_t k = 5;

    {
        std::vector<std::string> sequences {
            "TCATTGCGGACCTTGGCCC",
            "TCATTGCTACCTGGCCC",
            "TCATTGCCACCCTTTGGCCC"
        };
        const std::vector<std::string> labels { "A", "B", "C" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);
        auto& background = dynamic_cast<const DeBruijnGraph&>(anno_graph->get_graph());
        auto foreground = build_masked_graph(*anno_graph, ingroup, outgroup);

        Cigar::initialize_opt_table(foreground.get_graph().alphabet());

        auto config = DBGAlignerConfig(DBGAlignerConfig::unit_scoring_matrix(
            1, foreground.get_graph().alphabet()
        ));

        foreground.call_unitigs([&](auto unitig) {
            EXPECT_EQ(sequences[0], unitig);

            std::vector<DeBruijnGraph::node_index> nodes;
            foreground.map_to_nodes(unitig, [&](auto i) { nodes.push_back(i); });
            auto seeds = annotated_graph_algorithm::build_breakpoint_seeder_builder(background)(
                nodes, foreground
            )(
                foreground,
                config,
                &*unitig.begin(),
                &*unitig.end(),
                0,
                false
            );

            EXPECT_EQ(1u, seeds.size());
            EXPECT_EQ("TCATTGC", seeds.front().get_sequence());
        });
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, BreakpointSeeder) {
    test_breakpoint_seeder<typename TypeParam::first_type,
                           typename TypeParam::second_type>();
}

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_find_variants_indel(size_t pool_size = 0) {
    const std::vector<std::string> ingroup { "A" };
    const std::vector<std::string> outgroup { };
    size_t k = 5;

    {
        std::vector<std::string> sequences {
            "TCATTGC""GGACC"  "TTGGCCC",
            "TCATTGC"  "ACC"  "TTGGCCC", // 7=2I6= (CTTGGCCC masked out)
            "TCATTGC" "CACCCT""TTGGCCC"  // 7=1I1X3=2D6= (GGCCC masked out)
        };
        const std::vector<std::string> labels { "A", "B", "C" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        std::mutex add_mutex;
        ThreadPool thread_pool(pool_size);
        const std::unordered_set<std::string> ref { "B" };

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);

        Cigar::initialize_opt_table(masked_dbg.get_graph().alphabet());

        auto score_matrix = DBGAlignerConfig::unit_scoring_matrix(
            1, masked_dbg.get_graph().alphabet()
        );

        annotated_graph_algorithm::call_variants(
            masked_dbg,
            *anno_graph,
            [&](auto&& path, const std::string &query_sequence, auto&& vlabels) {
                auto ref = std::string(path.get_query_begin(), path.get_query_end());
                const auto& var = path.get_sequence();

                ASSERT_EQ(query_sequence.c_str() + path.get_clipping(),
                          path.get_query_begin());

                auto index = path.front();
                auto cigar = path.get_cigar().to_string();

                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), var, index))
                    << k << " " << ref << " " << var;
                EXPECT_FALSE(all_mapped_match_first(masked_dbg, var, index))
                    << k << " " << ref << " " << var;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var.substr(0, k), index));

                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k << " " << ref;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k << " " << var;

                ASSERT_FALSE(cigar.empty());

                EXPECT_EQ("7=2I6=", cigar)
                    << k << " " << query_sequence << " " << ref << " " << var << " " << cigar;

                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
            },
            DBGAlignerConfig(score_matrix),
            default_extender<DeBruijnGraph::node_index>,
            &thread_pool
        );

        EXPECT_EQ(ref, obs_labels) << k;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindVariantsIndel) {
    test_find_variants_indel<typename TypeParam::first_type,
                             typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindVariantsIndelParallel) {
    test_find_variants_indel<typename TypeParam::first_type,
                             typename TypeParam::second_type>(3);
}

