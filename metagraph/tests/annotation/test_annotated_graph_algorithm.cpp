#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "boss_construct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"
#include "annotated_dbg.hpp"
#include "masked_graph.hpp"
#include "annotated_graph_algorithm.hpp"
#include "annotate_column_compressed.hpp"
#include "annotation_converters.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "threading.hpp"
#include "../graph/test_dbg_helpers.hpp"

const double cutoff = 0.0;

template <class Graph, class Annotation = annotate::ColumnCompressed<>>
std::unique_ptr<AnnotatedDBG>
build_anno_graph(uint64_t k,
                 const std::vector<std::string> &sequences,
                 const std::vector<std::string> &labels) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences);

    uint64_t num_nodes = graph->num_nodes();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        std::move(graph),
        std::unique_ptr<AnnotatedDBG::Annotator>(
            new annotate::ColumnCompressed<>(num_nodes)
        )
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(sequences[i], { labels[i] });
    }

    if (!std::is_same<Annotation, annotate::ColumnCompressed<>>::value)
        anno_graph = std::make_unique<AnnotatedDBG>(
            std::move(anno_graph->graph_),
            std::unique_ptr<AnnotatedDBG::Annotator>(
                annotate::convert<Annotation>(
                    std::move(dynamic_cast<annotate::ColumnCompressed<>&>(
                        *anno_graph->annotator_
                    )
                )).release()
            )
        );

    return anno_graph;
}

MaskedDeBruijnGraph build_masked_graph(const AnnotatedDBG &anno_graph,
                                       const std::vector<std::string> &ingroup,
                                       const std::vector<std::string> &outgroup) {
    return MaskedDeBruijnGraph(
        std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()),
        annotated_graph_algorithm::mask_nodes_by_label(
            anno_graph,
            ingroup, outgroup,
            [&](size_t incount, size_t outcount) {
                return incount == ingroup.size()
                    && outcount <= cutoff * (incount + outcount);
            }
        ).release()
    );
}


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


template <typename Graph>
class MaskedDeBruijnGraphAlgorithm : public DeBruijnGraphTest<Graph> {};

TYPED_TEST_CASE(MaskedDeBruijnGraphAlgorithm, GraphTypes);

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, CallSignificantIndices) {
    test_call_significant_indices<TypeParam>();
    test_call_significant_indices<TypeParam, annotate::RowFlatAnnotator>();
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

        annotated_graph_algorithm::call_breakpoints(
            masked_dbg,
            *anno_graph,
            [&](const auto &index, const auto &ref, const auto &var, const auto &vlabels) {
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index)) << k;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, ref + var, index)) << k;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(ref + var, 1.0)) << k;
                for (const auto &label : vlabels) {
                    EXPECT_EQ(std::string("B"), label) << k;
                }
            },
            &thread_pool
        );
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBreakpoints) {
    test_find_breakpoints<TypeParam>();
    test_find_breakpoints<TypeParam, annotate::RowFlatAnnotator>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBreakpointsParallel) {
    test_find_breakpoints<TypeParam>(3);
    test_find_breakpoints<TypeParam, annotate::RowFlatAnnotator>(3);
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
            [&](const auto &index, const auto &ref, const auto &var, const auto &vlabels) {
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index)) << k;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var, index)) << k;
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
    test_find_bubbles<TypeParam>();
    test_find_bubbles<TypeParam, annotate::RowFlatAnnotator>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesParallel) {
    test_find_bubbles<TypeParam>(3);
    test_find_bubbles<TypeParam, annotate::RowFlatAnnotator>(3);
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

        annotated_graph_algorithm::call_bubbles(
            masked_dbg,
            *anno_graph,
            [&](const auto &index, const auto &ref, const auto &var, const auto &vlabels) {
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index)) << k;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var, index)) << k;
                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k;
                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
            },
            &thread_pool
        );

        EXPECT_EQ(std::unordered_set<std::string>{}, obs_labels) << k;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesIncomplete) {
    test_find_bubbles_incomplete<TypeParam>();
    test_find_bubbles_incomplete<TypeParam, annotate::RowFlatAnnotator>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesIncompleteParallel) {
    test_find_bubbles_incomplete<TypeParam>(3);
    test_find_bubbles_incomplete<TypeParam, annotate::RowFlatAnnotator>(3);
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

        annotated_graph_algorithm::call_bubbles(
            masked_dbg,
            *anno_graph,
            [&](const auto &index, const auto &ref, const auto &var, const auto &vlabels) {
                EXPECT_TRUE(all_mapped_match_first(anno_graph->get_graph(), ref, index)) << k;
                EXPECT_TRUE(all_mapped_match_first(masked_dbg, var, index)) << k;
                EXPECT_EQ(std::vector<std::string>{ "A" },
                          anno_graph->get_labels(ref, 1.0)) << k;
                EXPECT_EQ(std::vector<std::string>{ "B" },
                          anno_graph->get_labels(var, 1.0)) << k;
                std::lock_guard<std::mutex> lock(add_mutex);
                obs_labels.insert(vlabels.begin(), vlabels.end());
            },
            &thread_pool
        );

        EXPECT_EQ(std::unordered_set<std::string>{}, obs_labels) << k;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesInnerLoop) {
    test_find_bubbles_inner_loop<TypeParam>();
    test_find_bubbles_inner_loop<TypeParam, annotate::RowFlatAnnotator>();
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, FindBubblesInnerLoopParallel) {
    test_find_bubbles_inner_loop<TypeParam>(3);
    test_find_bubbles_inner_loop<TypeParam, annotate::RowFlatAnnotator>(3);
}
