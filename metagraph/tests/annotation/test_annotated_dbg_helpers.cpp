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
using namespace mtg::annot::matrix;
using common::logger;

typedef typename RowDiffCoordAnnotator::binary_matrix_type CoordRowDiff;

template <class Graph, class Annotation>
std::unique_ptr<AnnotatedDBG> build_anno_graph(uint64_t k,
                                               const std::vector<std::string> &sequences,
                                               const std::vector<std::string> &labels,
                                               DeBruijnGraph::Mode mode,
                                               bool coordinates) {
    assert(sequences.size() == labels.size());
    auto graph = build_graph_batch<Graph>(k, sequences, mode);

    // TODO: what if CanonicalDBG is not the highest level? find a better way to do this
    auto canonical = dynamic_pointer_cast<const CanonicalDBG>(graph);
    const DeBruijnGraph *base_graph = canonical ? &canonical->get_graph() : graph.get();
    uint64_t max_index = base_graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<ColumnCompressed<>>(max_index)
    );

    if (coordinates) {
        std::vector<std::tuple<std::string, std::vector<std::string>, uint64_t>> data;
        for (size_t i = 0; i < sequences.size(); ++i) {
            data.emplace_back(sequences[i], std::vector<std::string>{ labels[i] }, 0);
        }
        anno_graph->annotate_kmer_coords(data);
        if constexpr(std::is_same_v<Annotation, ColumnCompressed<>>) {
            std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
            std::string out_path = tmp_dir/"test_col";
            anno_graph->get_annotator().serialize(out_path);
            return std::make_unique<AnnotatedDBG>(
                graph,
                std::unique_ptr<AnnotatedDBG::Annotator>(new ColumnCoordAnnotator(
                    load_coords(std::move(dynamic_cast<ColumnCompressed<>&>(*anno_graph->annotator_)),
                                std::vector<std::string>{out_path + ColumnCompressed<>::kExtension})
                ))
            );
        } else if constexpr(!std::is_same_v<Annotation, RowDiffColumnAnnotator>) {
            static_assert("Coordinates only supported for ColumnCompressed and RowDiff");
        }
    } else {
        for (size_t i = 0; i < sequences.size(); ++i) {
            anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
        }

        if constexpr(std::is_same_v<Annotation, ColumnCompressed<>>)
            return anno_graph;

        if constexpr(!std::is_same_v<Annotation, RowDiffColumnAnnotator>
                && !std::is_same_v<Annotation, RowDiffCoordAnnotator>) {
            return std::make_unique<AnnotatedDBG>(
                graph,
                std::unique_ptr<AnnotatedDBG::Annotator>(
                    convert<Annotation>(
                        std::move(dynamic_cast<ColumnCompressed<>&>(*anno_graph->annotator_))
                    )
                )
            );
        }
    }

    if constexpr(std::is_same_v<Annotation, RowDiffColumnAnnotator>) {
        static_assert(std::is_same_v<Graph, DBGSuccinct>);
        assert(dynamic_cast<const DBGSuccinct*>(base_graph));

        size_t num_columns = anno_graph->get_annotator().get_label_encoder().size();
        std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
        std::string out_path = tmp_dir/"test_col";
        anno_graph->get_annotator().serialize(out_path);
        std::vector<std::string> files { out_path + ColumnCompressed<>::kExtension };
        base_graph->serialize(tmp_dir/"test_col");
        std::string graph_fname = std::string(tmp_dir/"test_col") + base_graph->file_extension();
        auto annotator = std::make_unique<annot::ColumnCompressed<>>(0);
        convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, tmp_dir,
                            static_cast<annot::RowDiffStage>(0), "", false, coordinates);
        convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, tmp_dir,
                            static_cast<annot::RowDiffStage>(1), "", false, coordinates);
        convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, tmp_dir,
                            static_cast<annot::RowDiffStage>(2), "", false, coordinates);

        if (!annotator->merge_load(files)) {
            logger->error("Cannot load annotations");
            exit(1);
        }

        const std::string anchors_file = graph_fname + annot::binmat::kRowDiffAnchorExt;
        const std::string fork_succ_file = graph_fname + annot::binmat::kRowDiffForkSuccExt;
        if (!std::filesystem::exists(anchors_file)) {
            logger->error("Anchor bitmap {} does not exist.", anchors_file);
            std::exit(1);
        }
        if (!std::filesystem::exists(fork_succ_file)) {
            logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
            std::exit(1);
        }

        assert(annotator->num_labels() == num_columns);

        if (coordinates) {
            std::vector<bit_vector_smart> delimiters;
            std::vector<sdsl::int_vector<>> column_values;

            typedef ColumnCoordAnnotator::binary_matrix_type CoordDiff;
            auto coords_fname = utils::remove_suffix(files[0], annot::ColumnCompressed<>::kExtension)
                                                            + annot::ColumnCompressed<>::kCoordExtension;
            std::ifstream in(coords_fname, std::ios::binary);
            try {
                CoordDiff::load_tuples(in, num_columns, [&](auto&& delims, auto&& values) {
                    delimiters.emplace_back(std::move(delims));
                    column_values.emplace_back(std::move(values));
                });
            } catch (const std::exception &e) {
                logger->error("Couldn't load coordinates from {}\nException: {}", coords_fname, e.what());
                exit(1);
            } catch (...) {
                logger->error("Couldn't load coordinates from {}", coords_fname);
                exit(1);
            }

            auto row_diff = std::make_unique<CoordRowDiff>(
                static_cast<const DBGSuccinct*>(base_graph),
                CoordDiff(std::move(*annotator->release_matrix()),
                          std::move(delimiters), std::move(column_values))
            );
            row_diff->load_anchor(anchors_file);
            row_diff->load_fork_succ(fork_succ_file);
            auto annotator = std::make_unique<RowDiffCoordAnnotator>(
                std::move(row_diff),
                anno_graph->get_annotator().get_label_encoder()
            );
            return std::make_unique<AnnotatedDBG>(graph, std::move(annotator));
        } else {
            auto row_diff = std::make_unique<RowDiffColumnAnnotator::binary_matrix_type>(
                static_cast<const DBGSuccinct*>(base_graph), std::move(*annotator->release_matrix())
            );
            row_diff->load_anchor(anchors_file);
            row_diff->load_fork_succ(fork_succ_file);
            auto annotator = std::make_unique<RowDiffColumnAnnotator>(
                std::move(row_diff),
                anno_graph->get_annotator().get_label_encoder()
            );
            return std::make_unique<AnnotatedDBG>(graph, std::move(annotator));
        }
    } else {
        throw std::runtime_error("Unsupported combination");
    }
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

template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, RowDiffColumnAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);


} // namespace test
} // namespace mtg
