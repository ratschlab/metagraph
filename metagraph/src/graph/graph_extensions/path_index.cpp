#include "path_index.hpp"

#include <mutex>

#include <tsl/hopscotch_set.h>
#include <progress_bar.hpp>

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

constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;

static const std::vector<Label> DUMMY { Label(1, 1) };
static const size_t MAX_SIZE = 1000;

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::load(const std::string &filename_base) {
    auto in = utils::open_ifstream(filename_base + kPathIndexExtension);
    if (!in->good())
        return false;

    if (!paths_indices_.load(*in))
        return false;

    if (!path_boundaries_.load(*in))
        return false;

    logger->trace("Loaded {} paths", path_boundaries_.num_set_bits());

    if (!is_superbubble_start_.load(*in))
        return false;

    logger->trace("Loaded {} superbubbles", is_superbubble_start_.num_set_bits());

    try {
        superbubble_termini_.load(*in);
    } catch (...) {
        return false;
    }

    return true;
}

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::serialize(const std::string &filename_base) const {
    std::ofstream fout(filename_base + kPathIndexExtension);
    paths_indices_.serialize(fout);
    path_boundaries_.serialize(fout);
    is_superbubble_start_.serialize(fout);
    superbubble_termini_.serialize(fout);
}

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::set_graph(std::shared_ptr<const DBGSuccinct> graph) {
    dbg_succ_ = graph;
    if constexpr(std::is_base_of_v<IRowDiff, PathStorage>) {
        static_cast<IRowDiff&>(paths_indices_).set_graph(dbg_succ_.get());
    }
}

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::PathIndex(std::shared_ptr<const DBGSuccinct> graph,
            const std::string &graph_name,
            const std::function<void(const std::function<void(std::string_view)>)> &generate_sequences) {
    const DBGSuccinct &dbg_succ = *graph;

    LabelEncoder<Label> label_encoder;
    label_encoder.insert_and_encode(DUMMY[0]);

    AnnotatedDBG anno_graph(std::const_pointer_cast<DBGSuccinct>(graph),
                            std::make_unique<ColumnCompressed<>>(dbg_succ.max_index()));

    auto &annotator = const_cast<ColumnCompressed<>&>(
        static_cast<const ColumnCompressed<>&>(anno_graph.get_annotator())
    );

    std::shared_ptr<const DeBruijnGraph> check_graph = graph;
    std::shared_ptr<const CanonicalDBG> canonical;
    if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
        canonical = std::make_shared<CanonicalDBG>(graph);
        check_graph = canonical;
    }

    std::vector<uint64_t> boundaries { 0 };
    std::vector<node_index> unitig_fronts;
    std::vector<node_index> unitig_backs;
    tsl::hopscotch_map<node_index, size_t> front_to_unitig_id;
    tsl::hopscotch_map<node_index, size_t> back_to_unitig_id;

    std::mutex mu;
    dbg_succ.call_unitigs([&](const auto &seq, const auto &path) {
        std::ignore = seq;
        auto rows = path;
        std::transform(rows.begin(), rows.end(), rows.begin(), AnnotatedDBG::graph_to_anno_index);

        std::lock_guard<std::mutex> lock(mu);
        front_to_unitig_id[path.front()] = unitig_fronts.size();
        back_to_unitig_id[path.back()] = unitig_fronts.size();
        unitig_fronts.emplace_back(path.front());
        unitig_backs.emplace_back(path.back());
        uint64_t coord = boundaries.back();
        annotator.add_labels(rows, DUMMY);
        for (auto row : rows) {
            annotator.add_label_coord(row, DUMMY, coord++);
        }

        boundaries.emplace_back(coord);
    }, get_num_threads());

    size_t num_unitigs = boundaries.size() - 1;

    size_t seq_count = 0;
    size_t total_seq_count = 0;

    generate_sequences([&](std::string_view seq) {
        total_seq_count += 1 + (dbg_succ.get_mode() != DeBruijnGraph::BASIC);
        auto nodes = map_to_nodes_sequentially(*check_graph, seq);

        if (std::any_of(nodes.begin(), nodes.end(), [&](node_index node) {
            return !node || check_graph->has_multiple_outgoing(node) || check_graph->indegree(node) > 1;
        })) {
            ++seq_count;
            uint64_t coord = boundaries.back();
            for (node_index &node : nodes) {
                if (node) {
                    if (canonical)
                        node = canonical->get_base_node(node);

                    annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(node), DUMMY, coord);
                }
                ++coord;
            }
            boundaries.emplace_back(coord);

            if (canonical) {
                ++seq_count;
                uint64_t coord = boundaries.back();
                for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
                    if (*it)
                        annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(*it), DUMMY, coord);

                    ++coord;
                }
                boundaries.emplace_back(coord);
            } else if (dbg_succ.get_mode() == DeBruijnGraph::CANONICAL) {
                ++seq_count;
                uint64_t coord = boundaries.back();
                std::string seq_rc(seq);
                reverse_complement_seq_path(dbg_succ, seq_rc, nodes);
                for (node_index node : nodes) {
                    if (node)
                        annotator.add_label_coord(AnnotatedDBG::graph_to_anno_index(node), DUMMY, coord);

                    ++coord;
                }
                boundaries.emplace_back(coord);
            }
        }
    });

    if (total_seq_count)
        logger->info("Indexed {} / {} sequences", seq_count, total_seq_count);

    assert(annotator.num_labels() <= 1);
    assert(std::adjacent_find(boundaries.begin(), boundaries.end()) == boundaries.end());

    path_boundaries_ = bit_vector_smart([&](const auto &callback) {
        std::for_each(boundaries.begin(), boundaries.end() - 1, callback);
    }, boundaries.back(), boundaries.size() - 1);

    logger->info("Indexed a total of {} paths", path_boundaries_.num_set_bits());

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

    // enumerate superbubbles
    sdsl::bit_vector is_superbubble_start(num_unitigs, false);

    // TODO: why does this get a bus error?
    auto superbubble_ends = aligned_int_vector(
        // num_unitigs * 2, 0, std::min(uint32_t(64), sdsl::bits::hi(check_graph->max_index()) + 1), 16
        num_unitigs * 2, 0, 64, 16
    );


    std::atomic_thread_fence(std::memory_order_release);

    ProgressBar progress_bar(num_unitigs, "Indexing superbubbles",
                             std::cerr, !common::get_verbose());
    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < num_unitigs; ++i) {
        ++progress_bar;
        tsl::hopscotch_set<size_t> visited;
        tsl::hopscotch_map<size_t, tsl::hopscotch_set<size_t>> seen;
        std::vector<std::pair<size_t, size_t>> traversal_stack;
        traversal_stack.emplace_back(i, 0);
        seen[i].emplace(0);
        while (traversal_stack.size()) {
            auto [unitig_id, dist] = traversal_stack.back();
            traversal_stack.pop_back();
            assert(!visited.count(unitig_id));

            visited.insert(unitig_id);
            bool has_children = false;
            bool has_cycle = false;
            size_t length = boundaries[unitig_id + 1] - boundaries[unitig_id];
            dbg_succ.call_outgoing_kmers(unitig_backs[unitig_id], [&](node_index next, char c) {
                if (c == boss::BOSS::kSentinel)
                    return;

                has_children = true;
                if (has_cycle)
                    return;

                if (next == unitig_fronts[i]) {
                    has_cycle = true;
                    return;
                }

                assert(front_to_unitig_id.count(next));
                size_t next_id = front_to_unitig_id[next];

                seen[next_id].emplace(dist + length);
                bool all_visited = true;
                dbg_succ.call_incoming_kmers(next, [&](node_index sibling, char c) {
                    if (c != boss::BOSS::kSentinel) {
                        assert(back_to_unitig_id.count(sibling));
                        size_t sibling_id = back_to_unitig_id[sibling];
                        if (all_visited && !visited.count(sibling_id))
                            all_visited = false;
                    }
                });

                if (all_visited)
                    traversal_stack.emplace_back(next_id, dist + length);
            });

            if (!has_children || has_cycle)
                break;

            if (traversal_stack.size() == 1 && visited.size() + 1 == seen.size()) {
                auto [unitig_id, dist] = traversal_stack.back();
                traversal_stack.pop_back();

                bool is_cycle = false;
                dbg_succ.adjacent_outgoing_nodes(unitig_backs[unitig_id], [&](node_index next) {
                    if (next == unitig_fronts[i])
                        is_cycle = true;
                });

                if (is_cycle)
                    continue;

                if (std::any_of(seen.begin(), seen.end(),
                                [&](const auto &a) { return a.second.size() != 1; })) {
                    // superbubble with paths of different lengths (i.e., complex)
                    // TODO: handle later
                    continue;
                }

                set_bit(is_superbubble_start.data(), i, true, MO_RELAXED);
                for (const auto &[u_id, d] : seen) {
                    atomic_fetch_and_add(superbubble_ends, u_id * 2, i + 1, mu, MO_RELAXED);
                    atomic_fetch_and_add(superbubble_ends, u_id * 2 + 1, *d.begin(), mu, MO_RELAXED);
                }

                atomic_exchange(superbubble_ends, i * 2, unitig_id + 1, mu, MO_RELAXED);
                atomic_exchange(superbubble_ends, i * 2 + 1, dist, mu, MO_RELAXED);
            }
        }
    }

    std::atomic_thread_fence(std::memory_order_acquire);

    is_superbubble_start_ = SuperbubbleIndicator(std::move(is_superbubble_start));

    logger->info("Indexed {} simple superbubbles", is_superbubble_start_.num_set_bits());

    superbubble_termini_ = SuperbubbleStorage(std::move(superbubble_ends));
}

auto IPathIndex
::get_coords(const std::vector<node_index> &nodes) const -> std::vector<RowTuples> {
    sdsl::bit_vector picked(nodes.size(), true);

    std::vector<Row> rows;
    rows.reserve(nodes.size());
    tsl::hopscotch_map<size_t, std::vector<size_t>> path_id_to_nodes;

    for (size_t i = 0; i < nodes.size(); ++i) {
        if (!has_coord(nodes[i])) {
            picked[i] = false;
            continue;
        }

        rows.emplace_back(AnnotatedDBG::graph_to_anno_index(nodes[i]));
    }

    auto it = picked.begin();
    auto row_tuples = get_row_tuples(rows);
    std::vector<RowTuples> ret_val;
    ret_val.reserve(nodes.size());
    while (it != picked.end() && !*it) {
        ret_val.emplace_back();
        ++it;
    }

    for (auto &tuples : row_tuples) {
        VectorMap<size_t, Tuple> out_tuples;
        assert(tuples.size() <= 1);
        for (const auto &[c, tuple] : tuples) {
            assert(!c);
            assert(std::adjacent_find(tuple.begin(), tuple.end()) == tuple.end());
            for (auto coord : tuple) {
                size_t path_id = coord_to_path_id(coord);
                out_tuples[path_id].emplace_back(coord);
                path_id_to_nodes[path_id].emplace_back(it - picked.begin());
            }
        }
        ret_val.emplace_back(out_tuples.values_container().begin(),
                             out_tuples.values_container().end());
        ++it;
        while (it != picked.end() && !*it) {
            ret_val.emplace_back();
            ++it;
        }
    }

    return ret_val;
}

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
std::pair<size_t, size_t> PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_terminus(size_t path_id) const {
    if (--path_id < is_superbubble_start_.size() && is_superbubble_start_[path_id]) {
        return std::make_pair(superbubble_termini_[path_id * 2],
                              superbubble_termini_[path_id * 2 + 1]);
    }

    return {};
}

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
std::pair<size_t, size_t> PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_and_dist(size_t path_id) const {
    if (--path_id < is_superbubble_start_.size()) {
        if (is_superbubble_start_[path_id])
            return std::make_pair(path_id + 1, 0);

        return std::make_pair(superbubble_termini_[path_id * 2],
                              superbubble_termini_[path_id * 2 + 1]);
    }

    return {};
}

template <class PathStorage, class PathBoundaries, class SuperbubbleIndicator, class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::has_coord(node_index node) const {
    assert(dbg_succ_);
    return node != DeBruijnGraph::npos
        && dbg_succ_->get_node_sequence(node).find(boss::BOSS::kSentinel) == std::string::npos;
}

template class PathIndex<>;

} // namespace mtg::graph
