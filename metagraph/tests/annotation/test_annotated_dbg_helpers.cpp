#include "test_annotated_dbg_helpers.hpp"

#include "../graph/all/test_dbg_helpers.hpp"

#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/canonical_dbg.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"

// this next #include includes AnnotatedDBG. we need access to its protected
// members to modify the underlying annotator
#define protected public
#include "annotation/annotation_converters.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

namespace mtg {
namespace test {

using namespace mtg::graph;
using namespace mtg::annot;

template <class Graph, class Annotation>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels,
                                               DeBruijnGraph::Mode mode,
                                               bool coordinates) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences, mode);

    // TODO: what if CanonicalDBG is not the highest level? find a better way to do this
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(graph.get());
    uint64_t max_index = canonical
        ? canonical->get_graph().max_index()
        : graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<ColumnCompressed<>>(max_index)
    );

    if (coordinates) {
        if constexpr(std::is_same_v<Annotation, ColumnCompressed<>>) {
            std::vector<std::tuple<std::string, std::vector<std::string>, uint64_t>> data;
            for (size_t i = 0; i < sequences.size(); ++i) {
                data.emplace_back(sequences[i],
                                  std::vector<std::string>{ labels[i] },
                                  0);
            }
            anno_graph->annotate_kmer_coords(data);
            std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col_coord");
            anno_graph->get_annotator().serialize(tmp_dir/"test_col_coord");
            anno_graph = std::make_unique<AnnotatedDBG>(
                graph,
                std::unique_ptr<AnnotatedDBG::Annotator>(new ColumnCoordAnnotator(
                    load_coords(std::move(dynamic_cast<ColumnCompressed<>&>(*anno_graph->annotator_)),
                                std::vector<std::string>{
                                    std::string(tmp_dir/"test_col_coord")
                                        + ColumnCompressed<>::kExtension
                                })
                ))
            );
        } else {
            throw std::runtime_error("Coordinate unit tests only implemented for ColumnCompressed");
        }
    } else {
        for (size_t i = 0; i < sequences.size(); ++i) {
            anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
        }

        if (!std::is_same_v<Annotation, ColumnCompressed<>>) {
            anno_graph = std::make_unique<AnnotatedDBG>(
                graph,
                std::unique_ptr<AnnotatedDBG::Annotator>(
                    convert<Annotation>(
                        std::move(dynamic_cast<ColumnCompressed<>&>(*anno_graph->annotator_))
                    )
                )
            );
        }
    }

    return anno_graph;
}

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, ColumnCompressed<>>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGBitmap, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashOrdered, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashFast, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGHashString, RowFlatAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);


} // namespace test
} // namespace mtg
