#include "unitigs.hpp"

#include <mutex>

#include <progress_bar.hpp>

#include "cli/config/config.hpp"
#include "cli/load/load_graph.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "common/unix_tools.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "annotation/annotation_converters.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/threads/threading.hpp"


namespace mtg {
namespace graph {
namespace align {

using common::logger;

Unitigs::Unitigs(const cli::Config &config) : graph_(load_graph_impl(config.fnames[0])) {
    std::filesystem::path tmp_dir = utils::create_temp_dir(config.tmp_dir, "unitigs");
    std::string out_path = tmp_dir/"unitigs";
    std::vector<std::string> files;
    size_t width = sdsl::bits::hi(graph_->num_nodes()) + 1;
    std::vector<size_t> unitigs;
    std::vector<std::unique_ptr<bit_vector>> cols;

    if (!graph_->get_mask()) {
        logger->trace("Marking dummy k-mers");
        DBGSuccinct &ncgraph = const_cast<DBGSuccinct&>(*graph_);
        ncgraph.mask_dummy_kmers(get_num_threads(), false);
        valid_nodes_.reset(ncgraph.release_mask());
        ProgressBar progress_bar(valid_nodes_->num_set_bits(), "Transforming mask",
                                 std::cerr, !common::get_verbose());
        cols.emplace_back(std::make_unique<bit_vector_smart>([&](const auto &callback) {
            valid_nodes_->call_ones([&](node_index i) {
                assert(i);
                callback(AnnotatedDBG::graph_to_anno_index(i));
                ++progress_bar;
            });
        }, ncgraph.num_nodes(), valid_nodes_->num_set_bits()));
    } else {
        logger->error("Dummy mask already present in the graph. "
                      "Support for this not implemented yet.");
        exit(1);
    }

    annot::LabelEncoder<> encoder;
    encoder.insert_and_encode("");
    auto colcomp = std::make_unique<annot::ColumnCompressed<>>(
        std::move(cols), std::move(encoder), 1, tmp_dir, 1e9, width
    );

    std::vector<std::string> labels;
    labels.push_back("");

    logger->trace("Annotating unitigs");
    std::mutex mu;
    size_t counter = 0;
    size_t max_unitig = 0;

    ThreadPool pool(1);
    graph_->call_unitigs([&](const std::string&, const auto &path) {
        if (path.size() == 1)
            return;

#ifndef NDEBUG
        if (graph_->has_single_outgoing(path.back())
                && graph_->has_single_incoming(path.front())) {
            graph_->adjacent_outgoing_nodes(path.back(), [&](node_index next) {
                assert(next == path.front() || !(*valid_nodes_)[next]
                        || graph_->indegree(next) > 1);
            });
            graph_->adjacent_incoming_nodes(path.front(), [&](node_index prev) {
                assert(prev == path.back() || !(*valid_nodes_)[prev]
                        || graph_->outdegree(prev) > 1);
            });
        }
#endif

        {
            std::lock_guard<std::mutex> lock(mu);
            unitigs.emplace_back(path.front());
            unitigs.emplace_back(path.back());
            unitigs.emplace_back(counter);
            max_unitig = std::max({ path.front(), path.back(), max_unitig });
            counter += path.size();
        }

        std::vector<std::pair<annot::ColumnCompressed<>::Index, uint64_t>> coords;
        for (size_t i = 0; i < path.size(); ++i) {
            coords.emplace_back(AnnotatedDBG::graph_to_anno_index(path[i]), i + unitigs.back());
        }

        pool.enqueue([&colcomp,&labels,c=std::move(coords)]() {
            colcomp->add_label_coords(c, labels);
        });
    }, get_num_threads());

    pool.join();

    logger->trace("Initializing unitig vector");
    size_t num_unitigs = unitigs.size() / 3;
    // TODO: replace with int_vector_buffer
    sdsl::int_vector<> boundaries(num_unitigs * 2, 0, sdsl::bits::hi(max_unitig + 1) + 1);
    indicator_ = Indicator([&](const auto &callback) {
        ProgressBar progress_bar(num_unitigs, "Packing unitigs",
                                 std::cerr, !common::get_verbose());
        for (size_t i = 0, j = 0; i < unitigs.size(); i += 3, j += 2) {
            boundaries[j] = unitigs[i];
            boundaries[j + 1] = unitigs[i + 1];
            callback(unitigs[i + 2]);
            ++progress_bar;
        }
    }, counter, num_unitigs);
    boundaries_ = IDVector(std::move(boundaries));
    unitigs = std::vector<size_t>();

    logger->trace("Serializing initial annotation");
    files.push_back(out_path + colcomp->file_extension());
    colcomp->serialize(files[0]);
    colcomp.reset();
    graph_.reset();
    logger->trace("Compressing unitig index");
    logger->trace("Step 0");
    convert_to_row_diff(files,
                        config.fnames[0],
                        config.memory_available * 1e9,
                        config.max_path_length,
                        tmp_dir,
                        tmp_dir,
                        static_cast<annot::RowDiffStage>(0),
                        out_path + ".row_count", false, true);
    logger->trace("Step 1");
    convert_to_row_diff(files,
                        config.fnames[0],
                        config.memory_available * 1e9,
                        config.max_path_length,
                        tmp_dir,
                        tmp_dir,
                        static_cast<annot::RowDiffStage>(1),
                        out_path + ".row_reduction", false, true);
    logger->trace("Step 2");
    convert_to_row_diff(files,
                        config.fnames[0],
                        config.memory_available * 1e9,
                        config.max_path_length,
                        tmp_dir,
                        tmp_dir,
                        static_cast<annot::RowDiffStage>(2),
                        out_path + ".row_reduction", false, true);
    logger->trace("done");
    const std::string anchors_file = config.fnames[0] + annot::binmat::kRowDiffAnchorExt;
    if (!std::filesystem::exists(anchors_file)) {
        logger->error("Anchor bitmap {} does not exist.", anchors_file);
        std::exit(1);
    }
    const std::string fork_succ_file = config.fnames[0] + annot::binmat::kRowDiffForkSuccExt;
    if (!std::filesystem::exists(fork_succ_file)) {
        logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
        std::exit(1);
    }

    common::logger->trace("Loading column");
    auto annotator = std::make_unique<annot::ColumnCompressed<>>(0);
    if (!annotator->merge_load(files)) {
        logger->error("Cannot load annotations");
        exit(1);
    }

    logger->trace("Wrapping as TupleRowDiff");

    std::vector<bit_vector_smart> delimiters(1);
    std::vector<sdsl::int_vector<>> column_values(1);

    auto coords_fname = utils::remove_suffix(files[0], annot::ColumnCompressed<>::kExtension)
                                                    + annot::ColumnCompressed<>::kCoordExtension;
    typedef annot::ColumnCoordAnnotator::binary_matrix_type UnitigCoordDiff;
    std::ifstream in(coords_fname, std::ios::binary);
    try {
        UnitigCoordDiff::load_tuples(in, 1, [&](auto&& delims, auto&& values) {
            delimiters[0] = std::move(delims);
            column_values[0] = std::move(values);
        });
    } catch (const std::exception &e) {
        logger->error("Couldn't load coordinates from {}\nException: {}", coords_fname, e.what());
        exit(1);
    } catch (...) {
        logger->error("Couldn't load coordinates from {}", coords_fname);
        exit(1);
    }


    unitigs_ = UnitigCoords(nullptr,
                            UnitigCoordDiff(std::move(*annotator->release_matrix()),
                                            std::move(delimiters),
                                            std::move(column_values)));

    unitigs_.load_anchor(anchors_file);
    unitigs_.load_fork_succ(fork_succ_file);
    logger->trace("RowDiff support bitmaps loaded");
    load_graph(config.fnames[0]);
    unitigs_.set_graph(graph_.get());
}

bool Unitigs::load(const std::string &filename_base) {
    std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
    std::ifstream fin(fname, std::ios::binary);
    if (!fin.good())
        return false;

    if (!unitigs_.load(fin))
        return false;

    unitigs_.set_graph(graph_.get());

    switch (graph_->get_state()) {
        case boss::BOSS::State::STAT: {
            valid_nodes_.reset(new bit_vector_small());
            break;
        }
        case boss::BOSS::State::FAST: {
            valid_nodes_.reset(new bit_vector_stat());
            break;
        }
        case boss::BOSS::State::DYN: {
            valid_nodes_.reset(new bit_vector_dyn());
            break;
        }
        case boss::BOSS::State::SMALL: {
            valid_nodes_.reset(new bit_vector_small());
            break;
        }
    }

    if (!valid_nodes_->load(fin))
        return false;

    try {
        boundaries_.load(fin);
        return true;
    } catch (...) {
        return false;
    }

    return indicator_.load(fin);
}

void Unitigs::serialize(const std::string &filename_base) const {
    std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
    std::ofstream fout(fname, std::ios::binary);
    unitigs_.serialize(fout);
    valid_nodes_->serialize(fout);
    boundaries_.serialize(fout);
    indicator_.serialize(fout);
}

auto Unitigs::get_unitig(size_t unitig_id) const -> std::pair<node_index, node_index> {
    size_t unitig_id_offset = get_unitig_id_offset();
    if (unitig_id <= unitig_id_offset)
        return std::make_pair(unitig_id, unitig_id);

    unitig_id -= unitig_id_offset;
    return std::make_pair(boundaries_[(unitig_id - 1) * 2],
                          boundaries_[(unitig_id - 1) * 2 + 1]);
}

auto Unitigs::get_unitig_bounds(size_t unitig_id) const
        -> std::pair<std::pair<node_index, node_index>, std::pair<Coord, Coord>> {
    size_t unitig_id_offset = get_unitig_id_offset();
    if (unitig_id <= unitig_id_offset) {
        return std::make_pair(std::make_pair(unitig_id, unitig_id),
                              std::make_pair(0, 1));
    }

    unitig_id -= unitig_id_offset;
    return std::make_pair(std::make_pair(boundaries_[(unitig_id - 1) * 2],
                                         boundaries_[(unitig_id - 1) * 2 + 1]),
                          std::make_pair(indicator_.select1(unitig_id),
                                         indicator_.select1(unitig_id + 1))
    );
}

std::vector<size_t> Unitigs::get_unitig_ids(const std::vector<node_index> &nodes) const {
    auto [indicator, rows] = nodes_to_rows(nodes);
    logger->trace("Fetching unitig IDs");
    auto seed_tuples = unitigs_.get_row_tuples(rows);
    std::vector<size_t> results;
    results.reserve(nodes.size());

    size_t j = 0;
    size_t unitig_id_offset = get_unitig_id_offset();
    for (size_t i = 0; i < nodes.size(); ++i) {
        size_t unitig_id = std::numeric_limits<size_t>::max();
        if (indicator[i]) {
            const auto &coordinates = seed_tuples[j++];
            if (coordinates.size()) {
                assert(coordinates.size() == 1);
                assert(!coordinates[0].first);
                if (coordinates[0].second.size()) {
                    assert(coordinates[0].second.size() == 1);
                    unitig_id = indicator_.rank1(coordinates[0].second[0]);
                    results.emplace_back(unitig_id + unitig_id_offset);
                }
            }
        }

        if (unitig_id == std::numeric_limits<size_t>::max())
            results.emplace_back(get_base_node(nodes[i]));
    }

    return results;
}

auto Unitigs::get_unitig_ids_and_coordinates(const std::vector<node_index> &nodes) const
        -> std::vector<std::pair<size_t, Coord>> {
    auto [indicator, rows] = nodes_to_rows(nodes);
    logger->trace("Fetching unitig IDs");
    auto seed_tuples = unitigs_.get_row_tuples(rows);
    std::vector<std::pair<size_t, Coord>> results;
    results.reserve(nodes.size());

    size_t j = 0;
    size_t unitig_id_offset = graph_->max_index() + 1;
    for (size_t i = 0; i < nodes.size(); ++i) {
        size_t unitig_id = std::numeric_limits<size_t>::max();
        if (indicator[i]) {
            const auto &coordinates = seed_tuples[j++];
            if (coordinates.size()) {
                assert(coordinates.size() == 1);
                assert(!coordinates[0].first);
                if (coordinates[0].second.size()) {
                    assert(coordinates[0].second.size() == 1);
                    unitig_id = indicator_.rank1(coordinates[0].second[0]);
                    results.emplace_back(unitig_id + unitig_id_offset,
                                         coordinates[0].second[0]);
                }
            }
        }

        if (unitig_id == std::numeric_limits<size_t>::max())
            results.emplace_back(get_base_node(nodes[i]), 0);
    }

    return results;
}

std::shared_ptr<const DBGSuccinct> Unitigs::load_graph_impl(const std::string &fname) {
    common::logger->trace("Graph loading...");
    Timer timer;
    return std::dynamic_pointer_cast<const DBGSuccinct>(cli::load_critical_dbg(fname));
    common::logger->trace("Graph loaded in {} sec", timer.elapsed());
}

std::pair<sdsl::bit_vector, std::vector<uint64_t>>
Unitigs::nodes_to_rows(const std::vector<node_index> &nodes) const {
    sdsl::bit_vector indicator(nodes.size(), false);
    std::vector<uint64_t> rows;
    for (size_t i = 0; i < nodes.size(); ++i) {
        node_index node = get_base_node(nodes[i]);

        if ((*valid_nodes_)[graph_->kmer_to_boss_index(node)]) {
            indicator[i] = true;
            rows.emplace_back(AnnotatedDBG::graph_to_anno_index(node));
        }
    }

    return std::make_pair(std::move(indicator), std::move(rows));
}

auto Unitigs::get_base_node(node_index node) const -> node_index {
    if (node > graph_->max_index())
        return node - graph_->max_index();

    return node;
}

} // namespace align
} // namespace graph
} // namespace mtg
