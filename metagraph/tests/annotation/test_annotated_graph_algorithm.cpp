#include "gtest/gtest.h"

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
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
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
        *anno_graph = AnnotatedDBG(
            anno_graph->get_graph_ptr(),
            std::unique_ptr<AnnotatedDBG::Annotator>(
                annotate::convert<Annotation>(
                    std::move(*dynamic_cast<annotate::ColumnCompressed<>*>(
                        &anno_graph->get_annotation()
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
        std::static_pointer_cast<DeBruijnGraph>(anno_graph.get_graph_ptr()),
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
