#include "test_annotated_dbg_helpers.hpp"

#include "../graph/test_dbg_helpers.hpp"

#include "annotation_converters.hpp"
#include "annotated_graph_algorithm.hpp"


template <class Graph, class Annotation>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences);

    uint64_t num_nodes = graph->num_nodes();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        std::move(graph),
        std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
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
                ))
            )
        );

    return anno_graph;
}

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);


const double cutoff = 0.0;

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
        )
    );
}

MaskedDeBruijnGraph build_masked_graph_lazy(const AnnotatedDBG &anno_graph,
                                            const std::vector<std::string> &ingroup,
                                            const std::vector<std::string> &outgroup) {
    return MaskedDeBruijnGraph(
        std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()),
        annotated_graph_algorithm::mask_nodes_by_label(
            anno_graph,
            ingroup, outgroup,
            [&](UInt64Callback incounter, UInt64Callback outcounter) {
                uint64_t incount = incounter();
                if (incount != ingroup.size())
                    return false;

                uint64_t outcount = outcounter();
                return outcount <= cutoff * (incount + outcount);
            }
        )
    );
}
