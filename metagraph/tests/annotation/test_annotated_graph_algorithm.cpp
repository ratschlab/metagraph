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

const double cutoff = 0.0;

// Build AnnotatedDBG from vector of <label, sequence> pairs
template <class Graph>
std::pair<std::shared_ptr<DeBruijnGraph>, std::unique_ptr<AnnotatedDBG>>
build_anno_graph(
      uint64_t k,
      const std::vector<std::pair<std::string, std::string>> &sequences) {
    std::shared_ptr<DeBruijnGraph> graph { new Graph(k) };
    for (const auto &pair : sequences) {
        graph->add_sequence(std::string(pair.second.begin(), pair.second.end()));
    }
    uint64_t num_nodes = graph->num_nodes();

    std::unique_ptr<AnnotatedDBG> anno_graph {
        new AnnotatedDBG(
            graph,
            std::unique_ptr<AnnotatedDBG::Annotator>(
                new annotate::ColumnCompressed<>(num_nodes)
            )
        )
    };

    for (const auto &pair : sequences) {
        anno_graph->annotate_sequence(pair.second, { pair.first });
    }

    return std::make_pair(std::move(graph), std::move(anno_graph));
}

template <>
std::pair<std::shared_ptr<DeBruijnGraph>, std::unique_ptr<AnnotatedDBG>>
build_anno_graph<DBGBitmap>(
      uint64_t k,
      const std::vector<std::pair<std::string, std::string>> &sequences) {
    DBGBitmapConstructor constructor(k);
    for (const auto &pair : sequences) {
        constructor.add_sequence(std::string(pair.second.begin(), pair.second.end()));
    }
    std::shared_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };
    uint64_t num_nodes = graph->num_nodes();

    std::unique_ptr<AnnotatedDBG> anno_graph {
        new AnnotatedDBG(
            graph,
            std::unique_ptr<AnnotatedDBG::Annotator>(
                new annotate::ColumnCompressed<>(num_nodes)
            )
        )
    };

    for (const auto &pair : sequences) {
        anno_graph->annotate_sequence(pair.second, { pair.first });
    }

    return std::make_pair(std::move(graph), std::move(anno_graph));
}

template <>
std::pair<std::shared_ptr<DeBruijnGraph>, std::unique_ptr<AnnotatedDBG>>
build_anno_graph<DBGSuccinct>(
      uint64_t k,
      const std::vector<std::pair<std::string, std::string>> &sequences) {
    std::unique_ptr<DBGSuccinct> succinct { new DBGSuccinct(k) };
    for (const auto &pair : sequences) {
        succinct->add_sequence(std::string(pair.second.begin(), pair.second.end()));
    }
    succinct->mask_dummy_kmers(1, false);
    std::shared_ptr<DeBruijnGraph> graph { succinct.release() };
    uint64_t num_nodes = graph->num_nodes();

    std::unique_ptr<AnnotatedDBG> anno_graph {
        new AnnotatedDBG(
            graph,
            std::unique_ptr<AnnotatedDBG::Annotator>(
                new annotate::ColumnCompressed<>(num_nodes)
            )
        )
    };

    for (const auto &pair : sequences) {
        anno_graph->annotate_sequence(pair.second, { pair.first });
    }

    return std::make_pair(std::move(graph), std::move(anno_graph));
}


template <class Graph>
void test_call_significant_indices() {
    const std::vector<std::string> ingroup { "B", "C" };
    const std::vector<std::string> outgroup { "A" };

    for (size_t k = 3; k < 15; ++k) {
        const std::vector<std::pair<std::string, std::string>> sequences {
            { "A", std::string("T") + std::string(k - 1, 'A') + "T" },
            { "B", std::string("T") + std::string(k - 1, 'A') + "C" },
            { "C", std::string("T") + std::string(k - 1, 'A') + "C" },
            { "D", std::string("T") + std::string(k - 1, 'A') + "C" },
            { "E", std::string("T") + std::string(k - 1, 'A') + "G" }
        };

        auto graph_pair = build_anno_graph<Graph>(k, sequences);

        std::unordered_set<std::string> labels;
        const std::unordered_set<std::string> ref { "B", "C", "D" };

        auto masked_dbg = annotated_graph_algorithm::mask_insignificant_nodes(
            graph_pair.first,
            *graph_pair.second,
            ingroup, outgroup,
            [&](size_t incount, size_t outcount) {
                return incount == ingroup.size()
                    && outcount <= cutoff * (incount + outcount);
            }
        );
        EXPECT_EQ(graph_pair.second->num_nodes(), masked_dbg.num_nodes());

        masked_dbg.call_nodes(
            [&](const auto &index) {
                auto cur_labels = graph_pair.second->get_labels(index);

                labels.insert(cur_labels.begin(), cur_labels.end());
            }
        );

        EXPECT_EQ(ref, labels) << k;
    }
}

TEST(MaskedDeBruijnGraphAlgorithm, CallSignificantIndices) {
    test_call_significant_indices<DBGSuccinct>();
    test_call_significant_indices<DBGHashString>();
    test_call_significant_indices<DBGHashOrdered>();
    test_call_significant_indices<DBGBitmap>();
}

