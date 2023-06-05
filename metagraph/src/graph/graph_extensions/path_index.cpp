#include "path_index.hpp"

#include <mutex>

#include <progress_bar.hpp>

#include "graph/annotated_dbg.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "annotation/annotation_converters.hpp"
#include "graph/alignment/alignment.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "common/threads/threading.hpp"
#include "common/vector_set.hpp"
#include "common/batch_accumulator.hpp"


namespace mtg::graph {
using namespace annot;
using namespace annot::matrix;
using namespace annot::binmat;

using common::logger;

using Label = AnnotatedDBG::Label;
using Row = MultiIntMatrix::Row;
using Tuple = MultiIntMatrix::Tuple;
using Column = BinaryMatrix::Column;

const std::string ColumnPathIndex::UNITIG_FRONT_TAG(1, '\1');
const std::string ColumnPathIndex::UNITIG_BACK_TAG(1, '\4');
const std::string ColumnPathIndex::SUPERBUBBLE_TAG(1, '\2');
const std::string ColumnPathIndex::CHAIN_TAG(1, '\3');

const IntMatrix& ColumnPathIndex::get_topo_matrix() const {
    assert(dynamic_cast<const IntMatrix*>(&topo_annotator_->get_matrix()));
    return static_cast<const IntMatrix&>(topo_annotator_->get_matrix());
}

const MultiIntMatrix& ColumnPathIndex::get_coord_matrix() const {
    assert(dynamic_cast<const MultiIntMatrix*>(&anno_graph_->get_annotator().get_matrix()));
    return static_cast<const MultiIntMatrix&>(anno_graph_->get_annotator().get_matrix());
}

std::pair<std::shared_ptr<const AnnotatedDBG>,
          std::shared_ptr<const AnnotatedDBG::Annotator>>
merge_load(std::shared_ptr<DeBruijnGraph> graph, const std::vector<std::string> &prefixes);

ColumnPathIndex::ColumnPathIndex(std::shared_ptr<DeBruijnGraph> graph,
                                 const std::vector<std::string> &prefixes)
      : ColumnPathIndex(merge_load(graph, prefixes)) {}

ColumnPathIndex::ColumnPathIndex(std::shared_ptr<const AnnotatedDBG> anno_graph,
                                 std::shared_ptr<const AnnotatedDBG::Annotator> topo_annotator)
      : anno_graph_(anno_graph), topo_annotator_(topo_annotator) {
    if (!dynamic_cast<const MultiIntMatrix*>(&anno_graph_->get_annotator().get_matrix())) {
        throw std::runtime_error("Annotator must contain coordinates");
    }

    const auto *int_mat = dynamic_cast<const IntMatrix*>(&topo_annotator_->get_matrix());
    if (!int_mat) {
        throw std::runtime_error("Topology annotator must contain counts");
    }

    if (!dynamic_cast<const ColumnMajor*>(&int_mat->get_binary_matrix())) {
        throw std::runtime_error("Topology annotator must be column major");
    }
}

void ColumnPathIndex::annotate_columns(std::shared_ptr<DeBruijnGraph> graph,
                                       const std::string &out_prefix,
                                       const LabeledSeqGenerator &generator,
                                       size_t num_columns_cached,
                                       const std::string &tmp_dir_prefix,
                                       double memory_available_gb,
                                       size_t max_chunks_open) {
    if (graph->get_mode() != DeBruijnGraph::BASIC) {
        throw std::runtime_error("Only implemented for BASIC graphs");
    }
    const size_t k = graph->get_k();
    std::filesystem::path tmp_dir = utils::create_temp_dir(tmp_dir_prefix, "test_col");

    auto anno_graph = std::make_shared<AnnotatedDBG>(graph,
        std::make_unique<ColumnCompressed<>>(graph->max_index(), num_columns_cached,
                                             tmp_dir, memory_available_gb * 1000000000,
                                             8, max_chunks_open)
    );

    ThreadPool thread_pool(get_num_threads() > 1 ? get_num_threads() : 0);
    auto topo_annotator_ptr = std::make_shared<ColumnCompressed<>>(
        graph->max_index() + 1, num_columns_cached, tmp_dir, memory_available_gb * 1000000000,
        64, max_chunks_open
    );

    auto &annotator = const_cast<ColumnCompressed<>&>(static_cast<const ColumnCompressed<>&>(anno_graph->get_annotator()));
    auto &topo_annotator = *topo_annotator_ptr;

    const size_t batch_size = 1'000;
    const size_t batch_length = 100'000;

    BatchAccumulator<std::tuple<std::string, size_t, size_t, size_t, std::string, std::vector<std::string>, uint64_t>> batcher(
        [&](auto&& all_data) {
            std::vector<std::tuple<std::string, std::vector<std::string>, uint64_t>> data;
            data.reserve(all_data.size());

            std::vector<std::pair<bool, bool>> unitig_data;
            unitig_data.reserve(all_data.size());

            std::vector<std::tuple<std::string, std::vector<std::string>, uint64_t>> seq_data;
            seq_data.reserve(all_data.size());

            for (auto &[in_seq, unitig_id, sb_id, chain_id, local_coords, in_labels, coord] : all_data) {
                if (!sb_id)
                    sb_id = unitig_id;

                if (!chain_id)
                    chain_id = sb_id;

                // logger->info("Checking\t{}\t{},{},{}\t{}\t{}", in_seq, unitig_id, sb_id, chain_id, coord, local_coords);
                assert(sb_id >= unitig_id);
                assert(chain_id >= unitig_id);

                assert(coord < graph->max_index());
                const auto &[seq, labels, _] = data.emplace_back(
                    std::move(in_seq), std::move(in_labels), coord
                );
                seq_data.emplace_back(seq, labels, coord);
                unitig_data.emplace_back(sb_id == unitig_id, chain_id == unitig_id);

                if (local_coords.size()) {
                    Labels sec_labels;
                    sec_labels.reserve(labels.size());
                    for (const auto &label : labels) {
                        sec_labels.emplace_back(UNITIG_FRONT_TAG + label);
                    }
                    auto coord_split = utils::split_string(local_coords, ";");
                    assert(coord_split.size() == 2);

                    for (const auto &c : utils::split_string(coord_split[0], ",")) {
                        data.emplace_back(seq, sec_labels, atoll(c.c_str()));
                    }
                    for (const auto &c : utils::split_string(coord_split[1], ",")) {
                        data.emplace_back(seq, sec_labels, static_cast<uint64_t>(atoll(c.c_str())));
                    }
                }
            }

            thread_pool.enqueue([&anno_graph,d=std::move(data)]() {
                anno_graph->annotate_kmer_coords(std::move(d));
            });

            thread_pool.enqueue([&](auto &unitig_data, auto &seq_data) {
                assert(seq_data.size() == unitig_data.size());

                std::vector<std::vector<std::string>> u_labels_v;
                u_labels_v.reserve(unitig_data.size());
                for (size_t i = 0; i < unitig_data.size(); ++i) {
                    auto &u_labels = u_labels_v.emplace_back();
                    const auto &[seq, labels, coord] = seq_data[i];
                    u_labels.reserve(labels.size() * (2 + unitig_data[i].first + unitig_data[i].second));
                    for (const auto &label : labels) {
                        u_labels.emplace_back(graph::ColumnPathIndex::UNITIG_FRONT_TAG + label);
                        u_labels.emplace_back(graph::ColumnPathIndex::UNITIG_BACK_TAG + label);

                        if (unitig_data[i].first)
                            u_labels.emplace_back(graph::ColumnPathIndex::SUPERBUBBLE_TAG + label);

                        if (unitig_data[i].second)
                            u_labels.emplace_back(graph::ColumnPathIndex::CHAIN_TAG + label);
                    }

                    topo_annotator.add_labels({ coord }, u_labels);
                }
            }, std::move(unitig_data), std::move(seq_data));
        },
        batch_size, batch_length, batch_size
    );

    uint64_t coord = 0;
    VectorSet<std::string> label_set;

    generator([&](std::string&& sequence, size_t unitig_id, size_t superbubble_id, size_t chain_id, std::string&& local_coords, auto&& labels) {
        if (sequence.size() >= k) {
            uint64_t num_kmers = sequence.size() - k + 1;
            assert(coord + num_kmers <= graph->max_index());
            for (const auto &label : labels) {
                label_set.insert(graph::ColumnPathIndex::UNITIG_FRONT_TAG + label);
            }
            batcher.push_and_pay(sequence.size(),
                                 std::move(sequence),
                                 unitig_id,
                                 superbubble_id,
                                 chain_id,
                                 std::move(local_coords),
                                 std::move(labels),
                                 coord);
            coord += num_kmers;
        }
    });

    topo_annotator.add_labels({ coord }, label_set.values_container());

    batcher.process_all_buffered();
    thread_pool.join();

    BatchAccumulator<std::tuple<std::string, std::vector<std::string>, uint64_t>> unitig_batcher(
        [&](auto&& data) {
            for (auto &[seq, labels, coord] : data) {
                thread_pool.enqueue([&](std::string seq, auto labels, uint64_t coord) {
                    std::vector<std::string> u_labels_1;
                    std::vector<std::string> u_labels_2;
                    u_labels_1.reserve(labels.size());
                    u_labels_2.reserve(labels.size());
                    for (const auto &label : labels) {
                        u_labels_1.emplace_back(graph::ColumnPathIndex::UNITIG_FRONT_TAG + label);
                        u_labels_2.emplace_back(graph::ColumnPathIndex::UNITIG_BACK_TAG + label);
                    }
                    const auto &graph = anno_graph->get_graph();
                    std::string_view first_kmer(seq.c_str(), k);
                    std::string_view last_kmer(seq.c_str() + seq.size() - k, k);
                    auto first_node = map_to_nodes_sequentially(graph, first_kmer);
                    auto last_node = map_to_nodes_sequentially(graph, last_kmer);
                    assert(first_node[0]);
                    assert(last_node[0]);

                    topo_annotator.add_label_counts({ coord }, u_labels_1, first_node);
                    topo_annotator.add_label_counts({ coord }, u_labels_2, last_node);
                }, std::move(seq), std::move(labels), coord);
            }
        },
        batch_size, batch_length, batch_size
    );

    coord = 0;
    generator([&](std::string&& sequence, size_t, size_t, size_t, std::string&&, auto&& labels) {
        if (sequence.size() >= k) {
            uint64_t num_kmers = sequence.size() - k + 1;
            assert(coord + num_kmers <= graph->max_index());
            unitig_batcher.push_and_pay(sequence.size(), std::move(sequence), std::move(labels), coord);
            coord += num_kmers;
        }
    });

    unitig_batcher.process_all_buffered();
    thread_pool.join();

    logger->trace("Serializing annotations");
    annotator.serialize(out_prefix + ".global");
    topo_annotator.serialize(out_prefix + ".topo");
}

std::pair<std::shared_ptr<const AnnotatedDBG>,
          std::shared_ptr<const AnnotatedDBG::Annotator>>
merge_load(std::shared_ptr<DeBruijnGraph> graph, const std::vector<std::string> &prefixes) {
    std::vector<std::string> global_files;
    std::vector<std::string> topo_files;

    global_files.reserve(prefixes.size());
    topo_files.reserve(prefixes.size());

    for (const auto &prefix : prefixes) {
        global_files.emplace_back(prefix + ".global" + ColumnCompressed<>::kExtension);
        topo_files.emplace_back(prefix + ".topo" + ColumnCompressed<>::kExtension);
    }

    logger->trace("Reloading coordinates");
    ColumnCompressed<> global_annotator(0);
    if (!global_annotator.merge_load(global_files)) {
        logger->error("Failed to load coordinate annotator");
        exit(1);
    }

    auto merged_global_annotator = std::make_unique<ColumnCoordAnnotator>(load_coords(
        std::move(global_annotator), global_files
    ));
    logger->trace("Done");

    logger->trace("Reloading topology");
    ColumnCompressed<> topo_annotator(0);
    if (!topo_annotator.merge_load(topo_files)) {
        logger->error("Failed to load topology annotator");
        exit(1);
    }
    auto merged_topo_annotator = std::make_shared<IntColumnAnnotator>(
        load_counts(std::move(topo_annotator), topo_files)
    );
    logger->trace("Done");

    return std::make_pair(
        std::make_shared<AnnotatedDBG>(graph, std::move(merged_global_annotator)),
        merged_topo_annotator
    );
}

auto ColumnPathIndex::get_chain_info(const Label &check_label, node_index node) const -> InfoPair {
    // TODO: make this more efficient
    for (auto&& [label, ret_val] : get_chain_info({ node })) {
        if (label == check_label)
            return ret_val[0];
    }

    return {};
}

auto ColumnPathIndex::get_chain_info(const std::vector<node_index> &nodes) const
        -> std::vector<LabeledNodesInfo> {
    auto kmer_coords = anno_graph_->get_kmer_coordinates(
        nodes, anno_graph_->get_annotator().get_label_encoder().size(), 0.0, 0.0
    );

    VectorMap<Label, NodesInfo> ret_val;

    const auto &label_encoder = topo_annotator_->get_label_encoder();
    const auto &col_int_mat = get_topo_matrix();
    assert(dynamic_cast<const ColumnMajor*>(&col_int_mat.get_binary_matrix()));

    const auto &col_mat = static_cast<const ColumnMajor&>(col_int_mat.get_binary_matrix());

    for (auto &[stored_label, num_kmer_matches, node_coords] : kmer_coords) {
        assert(node_coords.size() == nodes.size());

        bool is_local_coord = (stored_label.size() && stored_label[0] == UNITIG_FRONT_TAG[0]);
        std::string label = is_local_coord ? stored_label.substr(1) : stored_label;

        auto &nodes_info = ret_val.try_emplace(label, NodesInfo{}).first.value();
        nodes_info.resize(node_coords.size());
        for (size_t i = 0; i < nodes_info.size(); ++i) {
            auto &coords = node_coords[i];
            auto &[chain_info, coord_info] = nodes_info[i];

            if (coords.empty())
                continue;

            if (!is_local_coord) {
                assert(coords.size() == 1);

                Column unitig_col = label_encoder.encode(UNITIG_FRONT_TAG + label);
                Column sb_col = label_encoder.encode(SUPERBUBBLE_TAG + label);
                Column chain_col = label_encoder.encode(CHAIN_TAG + label);

                uint64_t global_coord = coords[0];
                std::get<0>(coord_info) = global_coord;
                uint64_t start_coord = col_mat.data()[unitig_col]->prev1(global_coord);
                assert(col_mat.data()[unitig_col]->rank1(global_coord)
                    == col_mat.data()[unitig_col]->rank1(start_coord));
                chain_info = ChainInfo(
                    col_mat.data()[unitig_col]->rank1(global_coord),
                    col_mat.data()[unitig_col]->rank1(col_mat.data()[sb_col]->next1(start_coord)),
                    col_mat.data()[unitig_col]->rank1(col_mat.data()[chain_col]->next1(start_coord))
                );
                assert(std::get<0>(chain_info));
                assert(std::get<1>(chain_info));
                assert(std::get<2>(chain_info));
                assert(std::get<1>(chain_info) >= std::get<0>(chain_info));
                assert(std::get<2>(chain_info) >= std::get<0>(chain_info));
            } else {
                for (size_t i = 0; i < coords.size(); ++i) {
                    int64_t coord = static_cast<int64_t>(coords[i]);
                    if (coord >= 0) {
                        std::get<1>(coord_info).emplace_back(coord);
                    } else {
                        std::get<2>(coord_info).emplace_back(-(coord + 1));
                    }
                }
            }
        }
    }

    for (const auto &[label, coords] : ret_val) {
        for (const auto &[chain_info, coord_info] : coords) {
            if (!std::get<0>(chain_info))
                continue;

            assert((std::get<0>(chain_info) == std::get<1>(chain_info)
                    && std::get<0>(chain_info) == std::get<2>(chain_info))
                    || (std::get<1>(coord_info).size() && std::get<2>(coord_info).size()));
        }
    }

    return const_cast<std::vector<LabeledNodesInfo>&&>(ret_val.values_container());
}

void ColumnPathIndex::call_distances(const Label &label,
                                     const InfoPair &info_a,
                                     const InfoPair &info_b,
                                     const std::function<void(size_t)> &callback,
                                     int64_t max_distance,
                                     size_t max_steps) const {
    const auto &[unitig_id_a, sb_id_a, chain_id_a] = info_a.first;
    const auto &[unitig_id_b, sb_id_b, chain_id_b] = info_b.first;

    const auto &[global_coord_a, ds_from_start_a, ds_to_end_a] = info_a.second;
    const auto &[global_coord_b, ds_from_start_b, ds_to_end_b] = info_b.second;

    // logger->info("\t\tA: {},{},{},{}", unitig_id_a, sb_id_a, chain_id_a, global_coord_a);
    // logger->info("\t\tB: {},{},{},{}", unitig_id_b, sb_id_b, chain_id_b, global_coord_b);

    std::vector<int64_t> dists;

    if (unitig_id_a == unitig_id_b) {
        // same unitig
        dists.emplace_back(global_coord_b - global_coord_a);
    } else if (unitig_id_b == chain_id_a) {
        // target is at the end of the superbubble or chain

        if (ds_to_end_a.size()) {
            int64_t dist_in_unitig_a = global_coord_a - get_global_coord(label, unitig_id_a);
            int64_t dist_in_unitig_b = global_coord_b - get_global_coord(label, unitig_id_b);
            int64_t dist_offset = dist_in_unitig_b - dist_in_unitig_a;
            for (int64_t dend : ds_to_end_a) {
                dists.emplace_back(dend + dist_offset);
            }
        }
    } else if (chain_id_a == chain_id_b && sb_id_a != sb_id_b) {
        // same chain, different superbubble

        // if the source reaches the end, but the target doesn't, then the source can't
        // reach the target
        if (ds_to_end_a.size() && ds_to_end_b.empty())
            return;

        if (ds_from_start_a.size() > 1 && ds_from_start_b.size() > 1) {
            assert(false);
            throw std::runtime_error("This case in traversal not implemented");
        }

        assert(ds_from_start_a.size());
        assert(ds_from_start_b.size());

        int64_t dist_in_unitig_a = global_coord_a - get_global_coord(label, unitig_id_a);
        int64_t dist_in_unitig_b = global_coord_b - get_global_coord(label, unitig_id_b);
        int64_t dist_offset = dist_in_unitig_b - dist_in_unitig_a;
        for (int64_t d_start_a : ds_from_start_a) {
            for (int64_t d_start_b : ds_from_start_b) {
                dists.emplace_back(d_start_b - d_start_a + dist_offset);
            }
        }
    } else if (sb_id_a == sb_id_b) {
        // same superbubble

        // this was covered in the first case
        assert(chain_id_a || sb_id_a != unitig_id_b);

        // if the source reaches the end, but the target doesn't, then the source can't
        // reach the target
        if (ds_to_end_a.size() && ds_to_end_b.empty())
            return;

        assert(ds_from_start_a.size());

        // traverse in the same superbubble to get the distance
        std::vector<std::tuple<size_t, int64_t, Vector<int64_t>::const_iterator>> traversal_stack;
        int64_t coord_offset = get_global_coord(label, unitig_id_b + 1) - global_coord_b;
        int64_t dist_offset = ds_from_start_a.front() - get_global_coord(label, unitig_id_a) + global_coord_a + coord_offset;
        traversal_stack.emplace_back(unitig_id_a,
                                     get_global_coord(label, unitig_id_a + 1) - global_coord_a - coord_offset,
                                     ds_from_start_b.begin());
        while (traversal_stack.size()) {
            auto [u_id, dist, it] = traversal_stack.back();
            traversal_stack.pop_back();

            if (u_id == sb_id_a)
                continue;

            it = std::find_if(it, ds_from_start_b.end(),
                              [&](int64_t d) { return dist_offset + dist <= d; });

            if (it == ds_from_start_b.end() || dist > max_distance)
                continue;

            adjacent_outgoing_unitigs(label, u_id, [&](size_t next_u_id, int64_t len) {
                if (next_u_id == unitig_id_b) {
                    dists.emplace_back(dist + len);
                } else if (dist + len <= max_distance) {
                    traversal_stack.emplace_back(next_u_id, dist + len, it);
                }
            });
        }
    } else {
        // traverse normally
        std::vector<std::tuple<size_t, int64_t, size_t>> traversal_stack;
        if (chain_id_a == unitig_id_a) {
            int64_t coord_offset_b = get_global_coord(label, unitig_id_b + 1) - global_coord_b;
            traversal_stack.emplace_back(
                unitig_id_a,
                get_global_coord(label, unitig_id_a + 1) - global_coord_a - coord_offset_b,
                0
            );
        } else if (ds_to_end_a.size()) {
            int64_t coord_offset_b = get_global_coord(label, unitig_id_b + 1) - global_coord_b;
            traversal_stack.emplace_back(
                chain_id_a,
                get_global_coord(label, unitig_id_a) - global_coord_a - coord_offset_b,
                0
            );
        }

        while (traversal_stack.size()) {
            auto [u_id, dist, n_steps] = traversal_stack.back();
            traversal_stack.pop_back();

            if (n_steps >= max_steps || dist > max_distance)
                continue;

            adjacent_outgoing_unitigs(label, u_id, [&](size_t next_u_id, int64_t len) {
                if (next_u_id == unitig_id_b) {
                    if (ds_to_end_a.size()) {
                        for (int64_t dend : ds_to_end_a) {
                            dists.emplace_back(dend + dist + len);
                        }
                    } else {
                        dists.emplace_back(dist + len);
                    }
                } else if (dist + len <= max_distance) {
                    traversal_stack.emplace_back(next_u_id, dist + len, n_steps + 1);
                }
            });
        }
    }

    for (int64_t d : dists) {
        if (d > 0)
            callback(d);
    }

    if (dists.empty() || !max_steps)
        return;

    int64_t min_dist = *std::min_element(dists.begin(), dists.end());

    std::vector<int64_t> cycle_lengths;

    std::vector<std::tuple<size_t, size_t, std::vector<int64_t>>> traversal_stack;
    auto &init_dists = std::get<2>(traversal_stack.emplace_back(chain_id_a, 0, std::vector<int64_t>{}));
    if (ds_to_end_a.size()) {
        for (int64_t dstart : ds_from_start_a) {
            for (int64_t dend : ds_to_end_a) {
                if (dstart + dend + min_dist <= max_distance)
                    init_dists.emplace_back(dstart + dend);
            }
        }
    } else if (ds_from_start_a.empty()) {
        int64_t init_dist = get_global_coord(label, unitig_id_a + 1) - get_global_coord(label, unitig_id_a);
        if (init_dist + min_dist <= max_distance)
            init_dists.emplace_back(init_dist);
    } else {
        traversal_stack.clear();
    }

    while (traversal_stack.size()) {
        auto [u_id, n_steps, dists] = std::move(traversal_stack.back());
        traversal_stack.pop_back();

        adjacent_outgoing_unitigs(label, u_id, [&](size_t next_u_id, int64_t len) {
            if (next_u_id == unitig_id_a) {
                for (int64_t d : dists) {
                    cycle_lengths.emplace_back(d);
                }
                return;
            }

            auto next_info = get_chain_info(label, get_unitig_back(label, next_u_id));
            const auto &[next_unitig_id, next_sb_id, next_chain_id] = next_info.first;
            const auto &[next_global_coord, next_ds_from_start, next_ds_to_end] = next_info.second;

            if (next_chain_id == chain_id_a) {
                for (int64_t d : dists) {
                    cycle_lengths.emplace_back(d);
                }
                return;
            }

            if (n_steps + 1 >= max_steps)
                return;

            auto &next_dists = std::get<2>(traversal_stack.emplace_back(next_chain_id, n_steps + 1, std::vector<int64_t>{}));
            if (next_ds_to_end.size()) {
                for (int64_t dend : next_ds_to_end) {
                    for (int64_t d : dists) {
                        if (dend + d + min_dist <= max_distance)
                            next_dists.emplace_back(dend + d);
                    }
                }
            } else if (next_ds_from_start.empty()) {
                for (int64_t d : dists) {
                    if (len + d + min_dist <= max_distance)
                        next_dists.emplace_back(len + d);
                }
            } else {
                traversal_stack.pop_back();
            }
        });
    }

    if (cycle_lengths.empty())
        return;

    tsl::hopscotch_set<int64_t> seen_lengths;
    seen_lengths.emplace(0);
    bool inserted = true;
    while (inserted) {
        inserted = false;
        for (int64_t cycle : cycle_lengths) {
            for (int64_t prev_cycle : seen_lengths) {
                const auto &[it, try_inserted] = seen_lengths.emplace(prev_cycle + cycle);
                if (try_inserted) {
                    for (int64_t d : dists) {
                        if (d + *it > 0 && d + *it <= max_distance) {
                            callback(d + *it);
                            inserted = true;
                        }
                    }
                }
            }
        }
    }
}

int64_t ColumnPathIndex::get_global_coord(const Label &label, size_t unitig_id) const {
    assert(unitig_id);
    const auto &label_encoder = topo_annotator_->get_label_encoder();
    const auto &col_int_mat = get_topo_matrix();
    assert(dynamic_cast<const ColumnMajor*>(&col_int_mat.get_binary_matrix()));

    const auto &col_mat = static_cast<const ColumnMajor&>(col_int_mat.get_binary_matrix());
    Column unitig_col = label_encoder.encode(UNITIG_FRONT_TAG + label);
    return col_mat.data()[unitig_col]->select1(unitig_id);
}

ColumnPathIndex::node_index ColumnPathIndex
::get_unitig_back(const Label &label, size_t unitig_id) const {
    assert(unitig_id);
    const auto &label_encoder = topo_annotator_->get_label_encoder();
    Column unitig_back_col = label_encoder.encode(UNITIG_BACK_TAG + label);
    const auto &col_int_mat = get_topo_matrix();
    assert(dynamic_cast<const ColumnMajor*>(&col_int_mat.get_binary_matrix()));

    const auto &col_mat = static_cast<const ColumnMajor&>(col_int_mat.get_binary_matrix());
    int64_t global_coord = col_mat.data()[unitig_back_col]->select1(unitig_id);
    node_index ret_val = col_int_mat.get_value(global_coord, unitig_back_col);
    assert(ret_val);

    return ret_val;
}

void ColumnPathIndex
::adjacent_outgoing_unitigs(const Label &label,
                            size_t unitig_id,
                            const std::function<void(size_t, int64_t)> &callback) const {
    assert(unitig_id);
    const auto &label_encoder = anno_graph_->get_annotator().get_label_encoder();
    const auto &topo_label_encoder = topo_annotator_->get_label_encoder();
    Column unitig_col = label_encoder.encode(label);
    Column unitig_front_col = topo_label_encoder.encode(UNITIG_FRONT_TAG + label);
    node_index unitig_back = get_unitig_back(label, unitig_id);

    const auto &col_int_mat = get_topo_matrix();
    assert(dynamic_cast<const ColumnMajor*>(&col_int_mat.get_binary_matrix()));

    const auto &col_mat = static_cast<const ColumnMajor&>(col_int_mat.get_binary_matrix());
    const auto &front_col = *col_mat.data()[unitig_front_col];

    const auto &global_coord_mat = get_coord_matrix();
    anno_graph_->get_graph().adjacent_outgoing_nodes(unitig_back, [&](node_index next) {
        auto next_global_coords = global_coord_mat.get_tuple(
            AnnotatedDBG::graph_to_anno_index(next), unitig_col
        );
        if (next_global_coords.size()) {
            size_t next_unitig_id = front_col.rank1(next_global_coords[0]);
            assert(next_unitig_id);

            int64_t length = front_col.next1(next_global_coords[0] + 1) - next_global_coords[0];
            callback(next_unitig_id, length);
        }
    });
}

} // namespace mtg::graph
