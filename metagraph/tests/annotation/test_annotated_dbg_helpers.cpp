#include "test_annotated_dbg_helpers.hpp"

#include "../graph/all/test_dbg_helpers.hpp"

#include "annotation/annotation_converters.hpp"
#include "graph/annotated_graph_algorithm.hpp"


namespace mtg {
namespace test {

using namespace mtg::bitmap_graph;

template <class Graph, class Annotation>
std::unique_ptr<AnnotatedDBG>
build_anno_graph(uint64_t k,
                 const std::vector<std::string> &sequences,
                 const std::vector<std::string> &labels,
                 bool canonical) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences, canonical);

    uint64_t max_index = graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<annotate::ColumnCompressed<>>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    if (!std::is_same<Annotation, annotate::ColumnCompressed<>>::value)
        anno_graph = std::make_unique<AnnotatedDBG>(
            graph,
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

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, annotate::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, annotate::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, bool);


MaskedDeBruijnGraph build_masked_graph(const AnnotatedDBG &anno_graph,
                                       const std::vector<std::string> &ingroup,
                                       const std::vector<std::string> &outgroup,
                                       double mask_in_label_fraction,
                                       double mask_out_label_fraction,
                                       double other_label_fraction,
                                       double lazy_evaluation_density_cutoff) {
    size_t insize = ingroup.size();
    size_t outsize = outgroup.size();
    return MaskedDeBruijnGraph(
        std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()),
        annotated_graph_algorithm::mask_nodes_by_node_label(
            anno_graph,
            ingroup,
            outgroup,
            [=,&anno_graph](auto index, auto get_num_in_labels, auto get_num_out_labels) {
                assert(index != DeBruijnGraph::npos);

                size_t num_in_labels = get_num_in_labels();
                if (num_in_labels < mask_in_label_fraction * insize)
                    return false;

                size_t num_out_labels = get_num_out_labels();
                if (num_out_labels < mask_out_label_fraction * outsize)
                    return false;

                size_t num_total_labels = anno_graph.get_labels(index).size();

                return (num_total_labels - num_in_labels - num_out_labels)
                    <= other_label_fraction * num_total_labels;
            },
            lazy_evaluation_density_cutoff
        )
    );
}

} // namespace test
} // namespace mtg
