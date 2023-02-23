#include "path_index.hpp"

#include <mutex>

#include "graph/annotated_dbg.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "annotation/annotation_converters.hpp"
#include "graph/alignment/alignment.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg::graph {
using namespace annot;
using namespace annot::matrix;
using namespace annot::binmat;

using common::logger;

using Label = AnnotatedDBG::Label;
using Row = MultiIntMatrix::Row;
using Tuple = MultiIntMatrix::Tuple;

static const std::vector<Label> DUMMY { Label(1, 1) };
static const size_t MAX_SIZE = 1000;

template <class PathStorage, class PathBoundaries>
bool PathIndex<PathStorage, PathBoundaries>
::load(const std::string &filename_base) {
    auto in = utils::open_ifstream(filename_base + kPathIndexExtension);
    if (!in->good())
        return false;

    if (!paths_indices_.load(*in))
        return false;

    return path_boundaries_.load(*in);
}

template <class PathStorage, class PathBoundaries>
void PathIndex<PathStorage, PathBoundaries>
::serialize(const std::string &filename_base) const {
    std::ofstream fout(filename_base + kPathIndexExtension);
    paths_indices_.serialize(fout);
    path_boundaries_.serialize(fout);
}

template <class PathStorage, class PathBoundaries>
void PathIndex<PathStorage, PathBoundaries>
::set_graph(std::shared_ptr<const DBGSuccinct> graph) {
    dbg_succ_ = graph;
    if constexpr(std::is_base_of_v<IRowDiff, PathStorage>) {
        static_cast<IRowDiff&>(paths_indices_).set_graph(dbg_succ_.get());
    }
}

template <class PathStorage, class PathBoundaries>
PathIndex<PathStorage, PathBoundaries>
::PathIndex(std::shared_ptr<const DBGSuccinct> graph,
            const std::string &graph_name,
            const std::function<void(const std::function<void(std::string_view)>)> &generate_sequences) {
    if (graph->num_nodes() <= 1)
        return;

    const DBGSuccinct &dbg_succ = *graph;

    LabelEncoder<Label> label_encoder;
    label_encoder.insert_and_encode(DUMMY[0]);

    AnnotatedDBG anno_graph(std::const_pointer_cast<DBGSuccinct>(graph),
                            std::make_unique<ColumnCompressed<>>(dbg_succ.max_index()));

    auto &annotator = const_cast<ColumnCompressed<>&>(
        static_cast<const ColumnCompressed<>&>(anno_graph.get_annotator())
    );

    std::vector<uint64_t> boundaries { 0 };

    std::mutex mu;
    dbg_succ.call_sequences([&](const auto &seq, const auto &path) {
        auto rows = path;
        std::transform(rows.begin(), rows.end(), rows.begin(), AnnotatedDBG::graph_to_anno_index);

        std::lock_guard<std::mutex> lock(mu);
        uint64_t coord = boundaries.back() + dbg_succ.get_k() - 1;
        annotator.add_labels(rows, DUMMY);
        for (auto row : rows) {
            annotator.add_label_coord(row, DUMMY, coord++);
        }

        boundaries.emplace_back(coord);

        if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
            uint64_t coord = boundaries.back() + dbg_succ.get_k() - 1;
            for (auto it = rows.rbegin(); it != rows.rend(); ++it) {
                annotator.add_label_coord(*it, DUMMY, coord++);
            }
        } else if (dbg_succ.get_mode() == DeBruijnGraph::CANONICAL) {
            auto rc_rows = path;
            std::string rc_seq = seq;
            reverse_complement_seq_path(dbg_succ, rc_seq, rc_rows);
            std::transform(rc_rows.begin(), rc_rows.end(), rc_rows.begin(), AnnotatedDBG::graph_to_anno_index);
            uint64_t coord = boundaries.back() + dbg_succ.get_k() - 1;
            annotator.add_labels(rc_rows, DUMMY);
            for (auto row : rc_rows) {
                annotator.add_label_coord(row, DUMMY, coord++);
            }
            boundaries.emplace_back(coord);
        }
    }, get_num_threads());

    size_t seq_count = 0;
    size_t total_seq_count = 0;
    std::shared_ptr<const DeBruijnGraph> check_graph = graph;
    std::shared_ptr<const CanonicalDBG> canonical;
    if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
        canonical = std::make_shared<CanonicalDBG>(graph);
        check_graph = canonical;
    }

    generate_sequences([&](std::string_view seq) {
        total_seq_count += 1 + (dbg_succ.get_mode() != DeBruijnGraph::BASIC);
        auto nodes = map_to_nodes_sequentially(*check_graph, seq);

        if (std::any_of(nodes.begin(), nodes.end(), [&](node_index node) {
            return !node || check_graph->has_multiple_outgoing(node) || check_graph->indegree(node) > 1;
        })) {
            ++seq_count;
            uint64_t coord = boundaries.back() + dbg_succ.get_k() - 1;
            for (node_index &node : nodes) {
                if (node) {
                    if (canonical)
                        node = canonical->get_base_node(node);

                    annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(node), DUMMY, coord);
                }
                ++coord;
            }
            boundaries.emplace_back(boundaries.back() + seq.size());

            if (canonical) {
                ++seq_count;
                uint64_t coord = boundaries.back() + dbg_succ.get_k() - 1;
                boundaries.emplace_back(boundaries.back() + seq.size());
                for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
                    if (*it)
                        annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(*it), DUMMY, coord);

                    ++coord;
                }
            } else if (dbg_succ.get_mode() == DeBruijnGraph::CANONICAL) {
                ++seq_count;
                uint64_t coord = boundaries.back() + dbg_succ.get_k() - 1;
                boundaries.emplace_back(boundaries.back() + seq.size());
                std::string seq_rc(seq);
                reverse_complement_seq_path(dbg_succ, seq_rc, nodes);
                for (node_index node : nodes) {
                    if (node)
                        annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(node), DUMMY, coord);

                    ++coord;
                }
            }
        }
    });

    if (total_seq_count)
        logger->info("Indexed {} / {} sequences", seq_count, total_seq_count);

    assert(annotator.num_labels() == 1);
    assert(std::adjacent_find(boundaries.begin(), boundaries.end()) == boundaries.end());

    path_boundaries_ = bit_vector_smart([&](const auto &callback) {
        std::for_each(boundaries.begin(), boundaries.end() - 1, callback);
    }, boundaries.back(), boundaries.size() - 1);
    boundaries = decltype(boundaries)();

    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
    std::string out_path = tmp_dir/"test_col";
    annotator.serialize(out_path);

    std::vector<std::string> files { out_path + ColumnCompressed<>::kExtension };
    if (!std::filesystem::exists(files[0])) {
        logger->error("Failed to serialize annotation to {}.", files[0]);
        std::exit(1);
    }

    if constexpr(std::is_same_v<PathStorage, ColumnCoordAnnotator::binary_matrix_type>) {
        paths_indices_ = const_cast<PathStorage&&>(load_coords(
            std::move(annotator),
            files
        ).get_matrix());

    } else if constexpr(std::is_same_v<PathStorage, RowDiffCoordAnnotator::binary_matrix_type>) {
        std::string graph_fname = graph_name;
        if (graph_fname.empty()) {
            graph->serialize(out_path);
            graph_fname = out_path + graph->file_extension();
        }

        if (!std::filesystem::exists(graph_fname)) {
            logger->error("Graph path incorrect: {}.", graph_fname);
            std::exit(1);
        }

        {
            std::filesystem::path swap_dir = utils::create_temp_dir("", "swap_col");
            convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(0), out_path + ".row_count", false, true);
            convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(1), out_path + ".row_reduction", false, true);
            convert_to_row_diff(files, graph_fname, 100e9, 100, tmp_dir, swap_dir,
                                static_cast<RowDiffStage>(2), out_path, false, true);
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
        size_t num_columns = anno_graph.get_annotator().get_label_encoder().size();
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
            label_encoder,
            graph.get(),
            std::move(*diff_annotator->release_matrix()),
            std::move(delimiters), std::move(column_values)
        ));

        auto &row_diff = const_cast<PathStorage&>(dynamic_cast<const PathStorage&>(
            annotator->get_matrix()
        ));
        row_diff.load_anchor(anchors_file);
        row_diff.load_fork_succ(fork_succ_file);
        paths_indices_ = const_cast<PathStorage&&>(row_diff);

    } else {
        throw std::runtime_error("Only ColumnCoord and RowDiffCoord annotators supported");
    }

    set_graph(graph);
}

auto IPathIndex
::get_coords(const std::vector<node_index> &nodes) const -> std::vector<RowTuples> {
    std::vector<Row> rows;
    rows.reserve(nodes.size());
    std::transform(nodes.begin(), nodes.end(), std::back_inserter(rows),
                   AnnotatedDBG::graph_to_anno_index);

    auto row_tuples = get_row_tuples(rows);
    std::vector<RowTuples> ret_val;
    ret_val.reserve(row_tuples.size());
    for (auto &tuples : row_tuples) {
        VectorMap<size_t, Tuple> out_tuples;
        assert(tuples.size() <= 1);
        for (const auto &[c, tuple] : tuples) {
            assert(!c);
            assert(std::adjacent_find(tuple.begin(), tuple.end()) == tuple.end());
            for (auto coord : tuple) {
                out_tuples[coord_to_path_id(coord)].emplace_back(coord);
            }
        }
        ret_val.emplace_back(out_tuples.values_container().begin(),
                             out_tuples.values_container().end());
    }

    return ret_val;
}

// template <>
// auto PathIndex<>
// ::get_coords(const std::vector<node_index> &nodes) const -> std::vector<RowTuples> {
//     assert(dbg_succ_);
//     assert(std::all_of(nodes.begin(), nodes.end(),
//         [this](const auto &a) {
//             return dbg_succ_->get_node_sequence(a).find(boss::BOSS::kSentinel)
//                 == std::string::npos;
//         }
//     ));

//     return IPathIndex::get_coords(nodes);

//     // size_t max_path_id = path_boundaries_.num_set_bits() + 1;
//     // tsl::hopscotch_map<Row, std::vector<size_t>> row_inv_map;
//     // for (size_t i = 0; i < rows.size(); ++i) {
//     //     row_inv_map[rows[i]].emplace_back(i);
//     // }

//     // auto [rd_ids, rd_paths_trunc] = paths_indices_.get_rd_ids(rows);
//     // for (size_t i = 0; i < rd_paths_trunc.size(); ++i) {
//     //     if (rd_paths_trunc[i].size() > 1) {
//     //         for (size_t j = 0; j < rd_paths_trunc[i].size(); ++j) {
//     //             for (size_t idx : row_inv_map[rd_ids[rd_paths_trunc[i][j]]]) {
//     //                 ret_val[idx].emplace_back(max_path_id + i, Tuple{ j });
//     //             }
//     //         }
//     //     }
//     // }

//     // return ret_val;
// }

template class PathIndex<>;

} // namespace mtg::graph
