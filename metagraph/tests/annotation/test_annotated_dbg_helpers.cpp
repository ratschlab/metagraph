#include "test_annotated_dbg_helpers.hpp"

#include "../graph/all/test_dbg_helpers.hpp"

#include "annotation/annotation_converters.hpp"
#include "graph/annotated_graph_algorithm.hpp"


namespace mtg {
namespace test {

using namespace mtg::graph;

template <class Graph, class Annotation>
std::unique_ptr<AnnotatedDBG>
build_anno_graph(uint64_t k,
                 const std::vector<std::string> &sequences,
                 const std::vector<std::string> &labels) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences);

    uint64_t max_index = graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<annot::ColumnCompressed<>>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    if (!std::is_same<Annotation, annot::ColumnCompressed<>>::value)
        anno_graph = std::make_unique<AnnotatedDBG>(
            graph,
            std::unique_ptr<AnnotatedDBG::Annotator>(
                annot::convert<Annotation>(
                    std::move(dynamic_cast<annot::ColumnCompressed<>&>(
                        *anno_graph->annotator_
                    )
                ))
            )
        );

    return anno_graph;
}

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, annot::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, annot::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, annot::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, annot::ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, annot::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, annot::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, annot::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, annot::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, annot::RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&);


} // namespace test
} // namespace mtg
