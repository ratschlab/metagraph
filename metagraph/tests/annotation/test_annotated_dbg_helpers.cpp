#include "test_annotated_dbg_helpers.hpp"

#include "test_matrix_helpers.hpp"

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
using namespace mtg::annot::binmat;
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
    auto base_graph = canonical ? canonical->get_graph_ptr() : graph;

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
                && !std::is_same_v<Annotation, RowDiffDiskAnnotator>
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
        assert(dynamic_cast<const DBGSuccinct*>(base_graph.get()));

        std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
        auto out_fs_path = tmp_dir/"test_col";
        std::string out_path = out_fs_path;
        anno_graph->get_annotator().serialize(out_path);
        std::vector<std::string> files { out_path + ColumnCompressed<>::kExtension };
        base_graph->serialize(out_path);
        std::string graph_fname = out_path + base_graph->file_extension();
        if (!std::filesystem::exists(files[0])) {
            logger->error("Failed to serialize annotation to {}.", files[0]);
            std::exit(1);
        }
        if (!std::filesystem::exists(graph_fname)) {
            logger->error("Failed to serialize graph to {}.", graph_fname);
            std::exit(1);
        }

        {
            std::filesystem::path swap_dir = utils::create_temp_dir("", "swap_col");
            convert_to_row_diff(files, graph_fname, 1e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(0), out_path + ".row_count", false, coordinates);
            convert_to_row_diff(files, graph_fname, 1e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(1), out_path + ".row_reduction", false, coordinates);
            convert_to_row_diff(files, graph_fname, 1e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(2), out_path, false, coordinates);
        }

        const std::string anchors_file = graph_fname + kRowDiffAnchorExt;
        const std::string fork_succ_file = graph_fname + kRowDiffForkSuccExt;
        if (!std::filesystem::exists(anchors_file)) {
            logger->error("Anchor bitmap {} does not exist.", anchors_file);
            std::exit(1);
        }
        if (!std::filesystem::exists(fork_succ_file)) {
            logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
            std::exit(1);
        }

        std::unique_ptr<AnnotatedDBG::Annotator> annotator;
        if (coordinates) {
            size_t num_columns = anno_graph->get_annotator().get_label_encoder().size();
            auto diff_annotator = std::make_unique<ColumnCompressed<>>(0);
            if (!diff_annotator->merge_load(files)) {
                logger->error("Cannot load annotations from {}", files[0]);
                exit(1);
            }
            assert(diff_annotator->num_labels() == num_columns);
            std::vector<bit_vector_smart> delimiters;
            std::vector<sdsl::int_vector<>> column_values;

            typedef ColumnCoordAnnotator::binary_matrix_type CoordDiff;
            auto coords_fname = utils::remove_suffix(files[0], ColumnCompressed<>::kExtension)
                                                            + ColumnCompressed<>::kCoordExtension;
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

            annotator.reset(new RowDiffCoordAnnotator(
                anno_graph->get_annotator().get_label_encoder(),
                static_cast<const DBGSuccinct*>(base_graph.get()),
                std::move(*diff_annotator->release_matrix()),
                std::move(delimiters), std::move(column_values)
            ));

        } else {
            auto rd_path = out_fs_path.replace_extension(RowDiffColumnAnnotator::kExtension);
            annotator.reset(new RowDiffColumnAnnotator({}, static_cast<const DBGSuccinct*>(base_graph.get())));
            if (!annotator->load(rd_path)) {
                logger->error("Cannot load annotations from {}", rd_path);
                exit(1);
            }
        }

        IRowDiff &row_diff = const_cast<IRowDiff&>(dynamic_cast<const IRowDiff&>(
            annotator->get_matrix()
        ));
        row_diff.load_anchor(anchors_file);
        row_diff.load_fork_succ(fork_succ_file);
        return std::make_unique<AnnotatedDBG>(graph, std::move(annotator));

    } else if constexpr(std::is_same_v<Annotation, RowDiffDiskAnnotator>) {
        static_assert(std::is_same_v<Graph, DBGSuccinct>);
        assert(dynamic_cast<const DBGSuccinct*>(base_graph.get()));

        std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
        auto out_fs_path = tmp_dir/"test_col";
        std::string out_path = out_fs_path;
        anno_graph->get_annotator().serialize(out_path);
        std::vector<std::string> files { out_path + ColumnCompressed<>::kExtension };
        base_graph->serialize(out_path);
        std::string graph_fname = out_path + base_graph->file_extension();
        if (!std::filesystem::exists(files[0])) {
            logger->error("Failed to serialize annotation to {}.", files[0]);
            std::exit(1);
        }
        if (!std::filesystem::exists(graph_fname)) {
            logger->error("Failed to serialize graph to {}.", graph_fname);
            std::exit(1);
        }

        {
            std::filesystem::path swap_dir = utils::create_temp_dir("", "swap_col");
            convert_to_row_diff(files, graph_fname, 1e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(0), out_path + ".row_count", false, coordinates);
            convert_to_row_diff(files, graph_fname, 1e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(1), out_path + ".row_reduction", false, coordinates);
            convert_to_row_diff(files, graph_fname, 1e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(2), out_path, false, coordinates);
        }

        convert_to_row_diff<RowDiffDiskAnnotator>({ out_path + RowDiffColumnAnnotator::kExtension },
                                                  graph_fname, out_fs_path,
                                                  get_num_threads(), 1e9, tmp_dir);

        auto rd_path = out_path + RowDiffDiskAnnotator::kExtension;
        auto annotator = std::make_unique<RowDiffDiskAnnotator>(
                annot::LabelEncoder(), static_cast<const DBGSuccinct*>(base_graph.get()));
        if (!annotator->load(rd_path)) {
            logger->error("Cannot load annotations from {}", rd_path);
            exit(1);
        }

        return std::make_unique<AnnotatedDBG>(graph, std::move(annotator));

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
template std::unique_ptr<AnnotatedDBG> build_anno_graph<DBGSuccinct, RowDiffDiskAnnotator>(uint64_t, const std::vector<std::string> &, const std::vector<std::string>&, DeBruijnGraph::Mode, bool);


} // namespace test
} // namespace mtg
