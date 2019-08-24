#include "gtest/gtest.h"

#include "static_annotators_def.hpp"
#include "threading.hpp"
#include "../graph/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"
#include "annotated_graph_algorithm.hpp"


template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_call_significant_indices() {
    const std::vector<std::string> ingroup { "B", "C" };
    const std::vector<std::string> outgroup { "A" };

    for (size_t k = 3; k < 15; ++k) {
        const std::vector<std::string> sequences {
            std::string("T") + std::string(k - 1, 'A') + "T",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "G"
        };
        const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        const std::unordered_set<std::string> ref { "B", "C", "D" };

        auto masked_dbg = build_masked_graph(*anno_graph, ingroup, outgroup);
        EXPECT_EQ(anno_graph->get_graph().num_nodes(), masked_dbg.num_nodes());

        masked_dbg.call_nodes(
            [&](const auto &index) {
                auto cur_labels = anno_graph->get_labels(index);

                obs_labels.insert(cur_labels.begin(), cur_labels.end());
            }
        );

        EXPECT_EQ(ref, obs_labels) << k;
    }
}

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
void test_call_significant_indices_lazy() {
    const std::vector<std::string> ingroup { "B", "C" };
    const std::vector<std::string> outgroup { "A" };

    for (size_t k = 3; k < 15; ++k) {
        const std::vector<std::string> sequences {
            std::string("T") + std::string(k - 1, 'A') + "T",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "C",
            std::string("T") + std::string(k - 1, 'A') + "G"
        };
        const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

        auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

        std::unordered_set<std::string> obs_labels;
        const std::unordered_set<std::string> ref { "B", "C", "D" };

        auto masked_dbg = build_masked_graph_lazy(*anno_graph, ingroup, outgroup);
        EXPECT_EQ(anno_graph->get_graph().num_nodes(), masked_dbg.num_nodes());

        masked_dbg.call_nodes(
            [&](const auto &index) {
                auto cur_labels = anno_graph->get_labels(index);

                obs_labels.insert(cur_labels.begin(), cur_labels.end());
            }
        );

        EXPECT_EQ(ref, obs_labels) << k;
    }
}



template <typename GraphAnnotationPair>
class MaskedDeBruijnGraphAlgorithm : public ::testing::Test {};

TYPED_TEST_CASE(MaskedDeBruijnGraphAlgorithm, GraphAnnotationPairTypes);

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, CallSignificantIndices) {
    test_call_significant_indices<typename TypeParam::first_type,
                                  typename TypeParam::second_type>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, CallSignificantIndicesLazy) {
    test_call_significant_indices_lazy<typename TypeParam::first_type,
                                       typename TypeParam::second_type>();
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
            [&](auto&& path, std::string&& ref, auto&& vlabels) {
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
            [&](auto&& path, std::string&& ref, auto&& vlabels) {
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
            [&](auto&& path, std::string&& ref, auto&& vlabels) {
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
            [&](auto&& path, std::string&& ref, auto&& vlabels) {
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
