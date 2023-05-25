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

            // std::cerr << "foo\t" << std::get<0>(chain_info) << ","
            //                      << std::get<1>(chain_info) << ","
            //                      << std::get<2>(chain_info) << "\t"
            //                      << std::get<1>(coord_info).size() << ","
            //                      << std::get<2>(coord_info).size() << std::endl;

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

    std::vector<int64_t> cycle_lengths;
    adjacent_outgoing_unitigs(label, chain_id_a, [&](size_t next_u_id, int64_t len) {
        if (cycle_lengths.size())
            return;
        auto next_info = get_chain_info(label, get_unitig_back(label, next_u_id));
        if (chain_id_a == std::get<0>(next_info.first)) {
            if (ds_from_start_a.size() && ds_to_end_a.size()) {
                for (int64_t dstart : ds_from_start_a) {
                    for (int64_t dend : ds_to_end_a) {
                        cycle_lengths.emplace_back(dstart + dend + len);
                    }
                }
            } else if (next_u_id == unitig_id_a) {
                cycle_lengths.emplace_back(len);
            }
        }
    });

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
        if ((chain_id_a || sb_id_a) && ds_to_end_a.empty())
            return;

        size_t cycle_target = chain_id_a ? chain_id_a : (sb_id_a ? sb_id_a : unitig_id_a);
        int64_t coord_offset_a = get_global_coord(label, cycle_target + 1)
            - (chain_id_a || sb_id_a ? get_global_coord(label, cycle_target) : global_coord_a);
        int64_t coord_offset_b = get_global_coord(label, unitig_id_b + 1) - global_coord_b;
        std::vector<std::tuple<size_t, int64_t, size_t>> traversal_stack;
        traversal_stack.emplace_back(cycle_target,
                                     coord_offset_a - coord_offset_b,
                                     0);
        while (traversal_stack.size()) {
            auto [u_id, dist, n_steps] = traversal_stack.back();
            traversal_stack.pop_back();

            if (u_id == sb_id_b || n_steps >= max_steps || dist > max_distance)
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

    if (cycle_lengths.empty())
        return;

    tsl::hopscotch_set<int64_t> seen_lengths;
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



constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;

static const std::vector<Label> DUMMY { Label(1, 1) };

static const IPathIndex::SuperbubbleStorage dummy { sdsl::int_vector(1, 0) };

void traverse(const IPathIndex &index,
              size_t start_id,
              size_t target,
              const std::function<void(size_t)> &callback,
              size_t max_dist,
              size_t max_search_depth,
              size_t stop_at = 0) {
    if (start_id == target) {
        callback(0);

        if (!max_dist)
            return;

        bool has_self_loop = false;
        index.adjacent_outgoing_unitigs(start_id, [&](size_t next) {
            if (next == start_id)
                has_self_loop = true;
        });

        if (!has_self_loop)
            return;

        size_t length_1 = index.path_length(start_id);
        for (size_t d = length_1; d <= max_dist; d += length_1) {
            callback(d);
        }

        return;
    }

    std::vector<std::tuple<size_t, size_t, size_t>> search;
    search.emplace_back(start_id, 0, max_search_depth);
    while (search.size()) {
        auto [cur_id, dist, depth] = search.back();
        search.pop_back();

        if (cur_id == target) {
            callback(dist);
            continue;
        }

        if (dist >= max_dist)
            continue;

        if (std::get<0>(index.get_superbubble_terminus(cur_id))) {
            index.call_dists(cur_id, target, [&](size_t d) { callback(dist + d); },
                             max_dist - dist, depth);
            continue;
        }

        if (depth) {
            dist += index.path_length(cur_id);
            index.adjacent_outgoing_unitigs(cur_id, [&](size_t next) {
                if (next != stop_at)
                    search.emplace_back(next, dist, depth - 1);
            });
        }
    }
};

void IPathIndex::call_dists(size_t path_id_1,
                            size_t path_id_2,
                            const std::function<void(size_t)> &callback,
                            size_t max_dist,
                            size_t max_search_depth) const {
    if (path_id_1 == path_id_2) {
        if (!is_unitig(path_id_1)) {
            callback(0);
        } else {
            traverse(*this, path_id_1, path_id_2, callback, max_dist, max_search_depth);
        }

        return;
    }

    if (!max_dist || !is_unitig(path_id_1) || !is_unitig(path_id_2))
        return;

    size_t length_1 = path_length(path_id_1);
    if (length_1 > max_dist)
        return;

    assert(dummy.size() == 1);
    size_t sb1 = path_id_1;
    size_t sb2 = path_id_2;
    auto d1_begin = dummy.begin();
    auto d1_end = dummy.end();
    auto d2_begin = dummy.begin();
    auto d2_end = dummy.end();

    bool is_source1 = std::get<0>(get_superbubble_terminus(path_id_1));
    bool is_source2 = std::get<0>(get_superbubble_terminus(path_id_2));

    if (!is_source1)
        std::tie(sb1, d1_begin, d1_end) = get_superbubble_and_dist(path_id_1);

    if (!is_source2)
        std::tie(sb2, d2_begin, d2_end) = get_superbubble_and_dist(path_id_2);

    assert(std::is_sorted(d1_begin, d1_end));
    assert(std::is_sorted(d2_begin, d2_end));

    if (!sb1) {
        traverse(*this, path_id_1, path_id_2, callback, max_dist, max_search_depth);
        return;
    }

    auto find_in_superbubble = [&](size_t path_id_1, bool is_source1, auto d1_begin, auto d1_end, const auto &callback) {
        if (is_source1) {
            std::for_each(d2_begin, d2_end, callback);
            return;
        }

        auto [t, d_begin, d_end] = get_superbubble_terminus(sb2);
        if (!can_reach_superbubble_terminus(path_id_1)) {
            if (path_id_2 == t || can_reach_superbubble_terminus(path_id_2))
                return;

            traverse(*this, path_id_1, path_id_2, callback, max_dist, max_search_depth, t);
            return;
        }

        if (path_id_2 == t && d_end - d_begin == 1) {
            size_t sb2_to_path_id_2 = *d_begin;
            std::for_each(d1_begin, d1_end, [&](size_t sb1_to_path_id_1) {
                callback(sb2_to_path_id_2 - sb1_to_path_id_1);
            });
            return;
        }

        traverse(*this, path_id_1, path_id_2, callback, max_dist, max_search_depth, t);
    };

    if (sb1 == sb2) {
        find_in_superbubble(path_id_1, is_source1, d1_begin, d1_end, callback);
        return;
    }

    if (!can_reach_superbubble_terminus(path_id_1))
        return;

    auto [t, d_begin, d_end] = get_superbubble_terminus(sb1);
    assert(t);
    assert(d_begin != d_end);
    assert(std::is_sorted(d_begin, d_end));
    assert(*(d_end - 1) >= *(d1_end - 1));

    std::vector<size_t> path_id_1_to_ts;
    if (is_source1) {
        assert(path_id_1 == sb1);
        auto it = *d_begin ? d_begin : d_begin + 1;
        path_id_1_to_ts.reserve(d_end - it);
        std::copy(it, d_end, std::back_inserter(path_id_1_to_ts));
    } else if (d_end == d_begin + 1) {
        size_t sb1_to_t = *d_begin;
        assert(sb1_to_t >= *(d1_end - 1));
        path_id_1_to_ts.resize(d1_end - d1_begin);
        std::transform(d1_begin, d1_end, path_id_1_to_ts.rbegin(),
                       [&](size_t sb1_to_path_id_1) { return sb1_to_t - sb1_to_path_id_1; });
    } else {
        traverse(*this, path_id_1, t, [&](size_t path_id_1_to_t) {
            if (path_id_1_to_t)
                path_id_1_to_ts.emplace_back(path_id_1_to_t);
        }, *(d_end - 1), max_search_depth);
        std::sort(path_id_1_to_ts.begin(), path_id_1_to_ts.end());
        path_id_1_to_ts.erase(std::unique(path_id_1_to_ts.begin(), path_id_1_to_ts.end()),
                              path_id_1_to_ts.end());
    }

    if (path_id_1_to_ts.empty())
        return;

    assert(std::is_sorted(path_id_1_to_ts.begin(), path_id_1_to_ts.end()));
    size_t path_id_1_to_t_min = path_id_1_to_ts.front();

    auto call_path_id_1_to_path_id_2 = [&](size_t t_to_path_id_2) {
        for (size_t path_id_1_to_t : path_id_1_to_ts) {
            callback(path_id_1_to_t + t_to_path_id_2);
        }
    };

    if (t == path_id_2) {
        call_path_id_1_to_path_id_2(0);
        return;
    }

    if (path_id_1_to_t_min > max_dist)
        return;

    if (t == sb2) {
        find_in_superbubble(t, true, dummy.begin(), dummy.end(),
                            call_path_id_1_to_path_id_2);

        return;
    }

    // now we need the distance from t to sb2
    auto [ct1, w1_begin, w1_end] = get_superbubble_chain(t);
    auto [ct2, w2_begin, w2_end] = get_superbubble_chain(sb2);
    assert(!ct1 || !std::get<0>(get_superbubble_terminus(ct1)));
    assert(!ct2 || !std::get<0>(get_superbubble_terminus(ct2)));

    auto call_to_path_id_2 = [&](size_t t_to_sb2) {
        std::for_each(d2_begin, d2_end, [&](size_t sb2_to_path_id_2) {
            size_t t_to_path_id_2 = t_to_sb2 + sb2_to_path_id_2;
            call_path_id_1_to_path_id_2(t_to_path_id_2);
        });
    };

    if (ct1) {
        size_t path_id_1_to_ct1_min = path_id_1_to_t_min + *w1_begin;

        if (ct1 == ct2 && w2_end == w2_begin + 1) {
            std::for_each(w1_begin, w1_end, [&](size_t t_to_ct1) {
                std::for_each(w2_begin, w2_end, [&](size_t sb2_to_ct1) {
                    if (t_to_ct1 >= sb2_to_ct1) {
                        size_t t_to_sb2 = t_to_ct1 - sb2_to_ct1;
                        call_to_path_id_2(t_to_sb2);
                    }
                });
            });
        }

        if (max_dist < path_id_1_to_ct1_min)
            return;

        if (!sb2 && !ct2) {
            traverse(*this, ct1, path_id_2, [&](size_t ct1_to_path_id_2) {
                std::for_each(w1_begin, w1_end, [&](size_t t_to_ct1) {
                    size_t t_to_path_id_2 = t_to_ct1 + ct1_to_path_id_2;
                    call_path_id_1_to_path_id_2(t_to_path_id_2);
                });
            }, max_dist - path_id_1_to_ct1_min, max_search_depth);

            return;
        }

        traverse(*this, ct1, sb2, [&](size_t ct1_to_sb2) {
            std::for_each(w1_begin, w1_end, [&](size_t t_to_ct1) {
                size_t t_to_sb2 = t_to_ct1 + ct1_to_sb2;
                call_to_path_id_2(t_to_sb2);
            });
        }, max_dist - path_id_1_to_ct1_min, max_search_depth);

        return;
    }

    if (sb2) {
        traverse(*this, t, sb2, call_to_path_id_2,
                 max_dist - path_id_1_to_t_min, max_search_depth);
    } else {
        traverse(*this, t, path_id_2, call_path_id_1_to_path_id_2,
                 max_dist - path_id_1_to_t_min, max_search_depth);
    }
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::adjacent_outgoing_unitigs(size_t path_id,
                            const std::function<void(size_t)> &callback) const {
    if (--path_id >= num_unitigs_)
        return;

    dbg_succ_->adjacent_outgoing_nodes(unitig_backs_[path_id], [&](node_index next) {
        if (size_t next_id = unitig_fronts_[next])
            callback(next_id);
    });
}

template <class Indicator, class Storage>
IPathIndex::SuperbubbleInfo get_range(const Indicator &indicator,
                                      const Storage &storage,
                                      size_t i) {
    assert(i);
    size_t begin = indicator.select1(i);
    size_t end = indicator.select1(i + 1);
    assert(begin < storage.size());
    assert(end <= storage.size());
    return IPathIndex::SuperbubbleInfo(storage[begin],
                                       storage.begin() + begin + 1,
                                       storage.begin() + end);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
IPathIndex::SuperbubbleInfo
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_terminus(size_t path_id) const {
    if (path_id > num_unitigs_ || !path_id)
        return SuperbubbleInfo(0, superbubble_termini_.end(), superbubble_termini_.end());

    return get_range(superbubble_termini_b_, superbubble_termini_, path_id);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
IPathIndex::SuperbubbleInfo
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_and_dist(size_t path_id) const {
    if (path_id > num_unitigs_ || !path_id)
        return SuperbubbleInfo(0, superbubble_sources_.end(), superbubble_sources_.end());

    return get_range(superbubble_sources_b_, superbubble_sources_, path_id);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::can_reach_superbubble_terminus(size_t path_id) const {
    return path_id && --path_id < num_unitigs_ && can_reach_terminus_[path_id];
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
IPathIndex::SuperbubbleInfo
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::get_superbubble_chain(size_t path_id) const {
    if (!path_id || --path_id >= num_unitigs_)
        return SuperbubbleInfo(0, unitig_chain_sizes_.end(), unitig_chain_sizes_.end());

    assert(path_id * 2 < unitig_chain_.size());
    if (size_t chain_term_id = unitig_chain_[path_id * 2]) {
        size_t dists_id = unitig_chain_[path_id * 2 + 1];

        size_t begin = unitig_chain_b_.select1(dists_id);
        size_t end = unitig_chain_b_.select1(dists_id + 1);
        assert(begin < unitig_chain_sizes_.size());
        assert(end <= unitig_chain_sizes_.size());

        return SuperbubbleInfo(chain_term_id,
            unitig_chain_sizes_.begin() + begin,
            unitig_chain_sizes_.begin() + end
        );
    }

    return SuperbubbleInfo(0, unitig_chain_sizes_.end(), unitig_chain_sizes_.end());
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
bool PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::load(const std::string &filename_base) {
    auto in = utils::open_ifstream(filename_base + kPathIndexExtension);
    if (!in->good())
        return false;

    num_unitigs_ = load_number(*in);
    num_superbubbles_ = load_number(*in);

    if (!dummy_indicator_.load(*in))
        return false;

    if (!paths_indices_.load(*in))
        return false;

    if (!path_boundaries_.load(*in))
        return false;

    if (!read_indices_.load(*in))
        return false;

    if (!read_boundaries_.load(*in))
        return false;

    logger->trace("Loaded {} paths", path_boundaries_.num_set_bits());

    try {
        unitig_backs_.load(*in);
        unitig_fronts_.load(*in);
        superbubble_sources_.load(*in);
        superbubble_termini_.load(*in);
    } catch (...) {
        return false;
    }

    if (!superbubble_sources_b_.load(*in))
        return false;

    if (!superbubble_termini_b_.load(*in))
        return false;

    if (!can_reach_terminus_.load(*in))
        return false;

    try {
        unitig_chain_.load(*in);
        unitig_chain_sizes_.load(*in);
    } catch (...) {
        return false;
    }

    if (!unitig_chain_b_.load(*in))
        return false;

    logger->trace("Loaded {} superbubbles", num_superbubbles_);

    return true;
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::serialize(const std::string &filename_base) const {
    std::ofstream fout(filename_base + kPathIndexExtension);
    serialize_number(fout, num_unitigs_);
    serialize_number(fout, num_superbubbles_);
    dummy_indicator_.serialize(fout);
    paths_indices_.serialize(fout);
    path_boundaries_.serialize(fout);
    read_indices_.serialize(fout);
    read_boundaries_.serialize(fout);
    unitig_backs_.serialize(fout);
    unitig_fronts_.serialize(fout);
    superbubble_sources_.serialize(fout);
    superbubble_termini_.serialize(fout);
    superbubble_sources_b_.serialize(fout);
    superbubble_termini_b_.serialize(fout);
    can_reach_terminus_.serialize(fout);
    unitig_chain_.serialize(fout);
    unitig_chain_sizes_.serialize(fout);
    unitig_chain_b_.serialize(fout);
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
void PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::set_graph(std::shared_ptr<const DBGSuccinct> graph) {
    dbg_succ_ = graph;
    if constexpr(std::is_base_of_v<IRowDiff, PathStorage>) {
        static_cast<IRowDiff&>(paths_indices_).set_graph(dbg_succ_.get());
        static_cast<IRowDiff&>(read_indices_).set_graph(dbg_succ_.get());
    }
}

template <class PathStorage>
PathStorage compress_path_storage(const DBGSuccinct *graph,
                                  ColumnCompressed<>&& annotator,
                                  const LabelEncoder<Label> &label_encoder,
                                  size_t num_columns,
                                  const std::string &graph_fname,
                                  const std::filesystem::path &tmp_dir,
                                  const std::string &out_path) {
    annotator.serialize(out_path);

    std::vector<std::string> files { out_path + ColumnCompressed<>::kExtension };
    if (!std::filesystem::exists(files[0])) {
        logger->error("Failed to serialize annotation to {}.", files[0]);
        std::exit(1);
    }

    if constexpr(std::is_same_v<PathStorage, ColumnCoordAnnotator::binary_matrix_type>) {
        return const_cast<PathStorage&&>(load_coords(
            std::move(annotator),
            files
        ).get_matrix());
    }

    if constexpr(std::is_same_v<PathStorage, RowDiffCoordAnnotator::binary_matrix_type>) {
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
            graph,
            std::move(*diff_annotator->release_matrix()),
            std::move(delimiters), std::move(column_values)
        ));

        auto &row_diff = const_cast<PathStorage&>(dynamic_cast<const PathStorage&>(
            annotator->get_matrix()
        ));
        row_diff.load_anchor(anchors_file);
        row_diff.load_fork_succ(fork_succ_file);
        return const_cast<PathStorage&&>(row_diff);
    }

    throw std::runtime_error("Only ColumnCoord and RowDiffCoord annotators supported");
}

template <class PathStorage,
          class PathBoundaries,
          class SuperbubbleIndicator,
          class SuperbubbleStorage>
PathIndex<PathStorage, PathBoundaries, SuperbubbleIndicator, SuperbubbleStorage>
::PathIndex(std::shared_ptr<const DBGSuccinct> graph,
            const std::string &graph_name,
            const std::function<void(const std::function<void(std::string_view)>)> &generate_sequences) {
    if (graph->get_mode() != DeBruijnGraph::BASIC) {
        throw std::runtime_error("Only BASIC graphs supported");
    }

    const DBGSuccinct &dbg_succ = *graph;

    LabelEncoder<Label> label_encoder;
    label_encoder.insert_and_encode(DUMMY[0]);

    AnnotatedDBG anno_graph(std::const_pointer_cast<DBGSuccinct>(graph),
                            std::make_unique<ColumnCompressed<>>(dbg_succ.max_index()));

    auto &annotator = const_cast<ColumnCompressed<>&>(
        static_cast<const ColumnCompressed<>&>(anno_graph.get_annotator())
    );

    std::shared_ptr<const DeBruijnGraph> check_graph = graph;

#if 0
    std::shared_ptr<const CanonicalDBG> canonical;
    if (dbg_succ.get_mode() == DeBruijnGraph::PRIMARY) {
        canonical = std::make_shared<CanonicalDBG>(graph);
        check_graph = canonical;
    }
#endif

    std::vector<uint64_t> boundaries { 0 };
    std::vector<node_index> unitig_fronts;
    std::vector<node_index> unitig_backs;
    tsl::hopscotch_map<node_index, size_t> front_to_unitig_id;
    tsl::hopscotch_map<node_index, size_t> back_to_unitig_id;

    if (!dbg_succ.get_mask()) {
        dummy_indicator_ = PathBoundaries(
            dbg_succ.get_boss().mark_all_dummy_edges(get_num_threads())
        );
    } else {
        dummy_indicator_ = PathBoundaries(dbg_succ.max_index() + 1, false);
    }

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

    num_unitigs_ = boundaries.size() - 1;
    logger->info("Indexed {} unitigs", num_unitigs_);

    // enumerate superbubbles
    {
        sdsl::bit_vector can_reach_terminus(num_unitigs_, false);
        std::vector<std::pair<uint64_t, std::vector<size_t>>> superbubble_starts(num_unitigs_);
        std::atomic<size_t> superbubble_start_size { num_unitigs_ };
        std::vector<std::pair<uint64_t, std::vector<size_t>>> superbubble_termini(num_unitigs_);
        std::atomic<size_t> superbubble_termini_size { num_unitigs_ };

        std::atomic<size_t> num_terminal_superbubbles { 0 };
        sdsl::bit_vector not_in_superbubble(num_unitigs_, false);

        std::atomic_thread_fence(std::memory_order_release);

        ProgressBar progress_bar(num_unitigs_, "Indexing superbubbles",
                                 std::cerr, !common::get_verbose());
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < num_unitigs_; ++i) {
            ++progress_bar;
            if (fetch_bit(not_in_superbubble.data(), i, MO_RELAXED))
                continue;

            tsl::hopscotch_set<size_t> visited;
            VectorMap<size_t, tsl::hopscotch_set<size_t>> seen;
            tsl::hopscotch_map<size_t, std::vector<size_t>> parents;
            std::vector<std::pair<size_t, size_t>> traversal_stack;
            traversal_stack.emplace_back(i, 0);
            seen[i].emplace(0);
            bool is_terminal_superbubble = false;
            size_t terminus = 0;
            size_t term_dist = 0;
            while (traversal_stack.size()) {
                auto [unitig_id, dist] = traversal_stack.back();
                traversal_stack.pop_back();
                assert(!visited.count(unitig_id));

                if (fetch_bit(not_in_superbubble.data(), unitig_id, MO_RELAXED)) {
                    is_terminal_superbubble = false;
                    break;
                }

                bool has_cycle = false;
                visited.insert(unitig_id);
                size_t num_children = 0;
                size_t length = boundaries[unitig_id + 1] - boundaries[unitig_id];
                dbg_succ.call_outgoing_kmers(unitig_backs[unitig_id], [&](node_index next, char c) {
                    if (c == boss::BOSS::kSentinel || next == unitig_fronts[unitig_id])
                        return;

                    ++num_children;

                    if (has_cycle)
                        return;

                    if (next == unitig_fronts[i]) {
                        has_cycle = true;
                        return;
                    }

                    assert(front_to_unitig_id.count(next));
                    size_t next_id = front_to_unitig_id[next];

                    bool add_parents = !seen.count(next_id);

                    seen[next_id].emplace(dist + length);
                    bool all_visited = true;
                    dbg_succ.call_incoming_kmers(next, [&](node_index sibling, char c) {
                        if (c != boss::BOSS::kSentinel) {
                            assert(back_to_unitig_id.count(sibling));
                            size_t sibling_id = back_to_unitig_id[sibling];
                            if (sibling_id == next_id)
                                return;

                            if (add_parents)
                                parents[next_id].emplace_back(sibling_id);

                            if (all_visited && !visited.count(sibling_id))
                                all_visited = false;
                        }
                    });

                    if (all_visited)
                        traversal_stack.emplace_back(next_id, dist + length);
                });

                bool reached_end = (traversal_stack.size() == 1 && visited.size() + 1 == seen.size());
                if (has_cycle) {
                    for (const auto &[u_id, d] : seen) {
                        set_bit(not_in_superbubble.data(), u_id, MO_RELAXED);
                    }
                    is_terminal_superbubble = false;
                    break;
                }

                if (!num_children) {
                    is_terminal_superbubble = true;
                }

                if (reached_end) {
                    auto [unitig_id, dist] = traversal_stack.back();
                    traversal_stack.pop_back();

                    terminus = unitig_id;
                    term_dist = dist;
                    for (const auto &[u_id, d] : seen) {
                        for (size_t dist : d) {
                            if (dist > term_dist) {
                                terminus = u_id;
                                term_dist = dist;
                            }
                        }
                    }

                    const auto &d = seen[terminus];

                    superbubble_termini_size.fetch_add(
                        static_cast<ssize_t>(d.size()) - superbubble_termini[i].second.size(),
                        MO_RELAXED
                    );
                    superbubble_termini[i] = std::make_pair(terminus + 1,
                        std::vector<size_t>(d.begin(), d.end()));
                    std::sort(superbubble_termini[i].second.begin(),
                              superbubble_termini[i].second.end());

                    for (const auto &[u_id, d] : seen) {
                        if (u_id == i)
                            continue;

                        std::vector<size_t> cur_d(d.begin(), d.end());
                        std::sort(cur_d.begin(), cur_d.end());

                        #pragma omp critical
                        {
                        if (!superbubble_starts[u_id].first
                                || cur_d.back() < superbubble_starts[u_id].second.back()) {
                            superbubble_start_size -= superbubble_starts[u_id].second.size();
                            superbubble_start_size += cur_d.size();
                            superbubble_starts[u_id] = std::make_pair(i + 1, std::move(cur_d));
                            can_reach_terminus[u_id] = !is_terminal_superbubble;
                            // set_bit(can_reach_terminus.data(), u_id, !is_terminal_superbubble, MO_RELAXED);
                        }
                        }
                    }

                    #pragma omp critical
                    can_reach_terminus[i] = true;
                    // set_bit(can_reach_terminus.data(), i, true, MO_RELAXED);
                }
            }

            if (is_terminal_superbubble && seen.size() > 1 && terminus) {
                sdsl::bit_vector found_map(seen.size(), false);
                std::vector<std::pair<size_t, std::vector<size_t>>> back_traversal_stack;
                back_traversal_stack.reserve(seen.size());
                back_traversal_stack.emplace_back(terminus,
                    std::vector<size_t>(seen[terminus].begin(), seen[terminus].end()));
                while (back_traversal_stack.size()) {
                    auto [cur_id, d] = back_traversal_stack.back();
                    back_traversal_stack.pop_back();

                    found_map[seen.find(cur_id) - seen.begin()] = true;
                    if (parents[cur_id].empty())
                        continue;

                    for (auto &dd : d) {
                        --dd;
                    }

                    for (size_t parent : parents[cur_id]) {
                        auto &[p, p_d] = back_traversal_stack.emplace_back(parent, std::vector<size_t>{});
                        for (auto dd : d) {
                            if (seen[parent].count(dd))
                                p_d.emplace_back(dd);
                        }

                        if (p_d.empty())
                            back_traversal_stack.pop_back();
                    }
                }

                if (!superbubble_termini[i].first) {
                    auto it = found_map.begin();
                    for (const auto &[cur_id, stuff] : seen) {
                        #pragma omp critical
                        {
                        if (superbubble_starts[cur_id].first == i + 1)
                            can_reach_terminus[cur_id] = *it;
                        }
                        ++it;
                    }
                }

                num_terminal_superbubbles.fetch_add(1, MO_RELAXED);
            }
        }

        std::atomic_thread_fence(std::memory_order_acquire);

        {
            sdsl::int_vector<> sources(superbubble_start_size);
            sdsl::bit_vector source_indicator(superbubble_start_size + 1);
            auto it = sources.begin();
            auto jt = source_indicator.begin();
            size_t j = 0;
            for (const auto &[start, d] : superbubble_starts) {
                *jt = true;
                *it = start;
                ++it;
                ++jt;
                std::copy(d.begin(), d.end(), it);
                it += d.size();
                jt += d.size();
                ++j;
            }

            assert(jt + 1 == source_indicator.end());
            *jt = true;
            superbubble_sources_ = SuperbubbleStorage(std::move(sources));
            superbubble_sources_b_ = SuperbubbleIndicator(std::move(source_indicator));
        }

        num_superbubbles_ = 0;
        size_t num_multiple_sizes = 0;
        size_t num_in_superbubble = num_unitigs_;
        {
            sdsl::int_vector<> termini(superbubble_termini_size);
            sdsl::bit_vector termini_indicator(superbubble_termini_size + 1);
            auto it = termini.begin();
            auto jt = termini_indicator.begin();
            size_t i = 1;
            for (const auto &[start, d] : superbubble_termini) {
                *jt = true;
                *it = start;
                ++it;
                ++jt;
                num_superbubbles_ += !d.empty();
                num_multiple_sizes += (d.size() > 1);
                if (!start && !std::get<0>(get_superbubble_and_dist(i)))
                    --num_in_superbubble;

                std::copy(d.begin(), d.end(), it);
                it += d.size();
                jt += d.size();
                ++i;
            }

            assert(jt + 1 == termini_indicator.end());
            *jt = true;
            superbubble_termini_ = SuperbubbleStorage(std::move(termini));
            superbubble_termini_b_ = SuperbubbleIndicator(std::move(termini_indicator));
        }

        can_reach_terminus_ = SuperbubbleIndicator(std::move(can_reach_terminus));

        logger->info("{} / {} unitigs are in a superbubble. Indexed {} superbubbles, of which {} have dead ends and {} have multiple paths to the terminus.",
                     num_in_superbubble, num_unitigs_,
                     num_superbubbles_,
                     num_terminal_superbubbles,
                     num_multiple_sizes);

        sdsl::int_vector<> unitig_backs_sdsl(num_unitigs_);
        std::copy(unitig_backs.begin(), unitig_backs.end(), unitig_backs_sdsl.begin());
        unitig_backs_ = SuperbubbleStorage(std::move(unitig_backs_sdsl));

        sdsl::int_vector<> unitig_fronts_sdsl(dbg_succ.max_index() + 1);
        for (size_t i = 0; i < unitig_fronts.size(); ++i) {
            unitig_fronts_sdsl[unitig_fronts[i]] = i + 1;
        }
        unitig_fronts_ = SuperbubbleStorage(std::move(unitig_fronts_sdsl));
    }

    // chain superbubbles
    {
        sdsl::int_vector<> chain_parent(num_unitigs_ + 1);
        for (size_t i = 1; i <= num_unitigs_; ++i) {
            auto [t, d_begin, d_end] = get_superbubble_terminus(i);
            if (!t || t == i)
                continue;

            // we are at a superbubble start
            auto [t2, d2_begin, d2_end] = get_superbubble_terminus(t);
            if (!t2)
                continue;

            // the superbubble starting at i is chained to the one starting at t
            assert(!chain_parent[t2]);
            chain_parent[t2] = t;
            chain_parent[t] = i;
        }

        sdsl::int_vector<> superbubble_chain(num_unitigs_ * 2);
        std::vector<bool> chain_bounds;
        std::vector<uint64_t> chain_widths;
        size_t chain_i = 0;
        size_t chain_widths_id = 0;
        for (size_t i = 0; i < num_unitigs_; ++i) {
            if (chain_parent[i + 1] && !std::get<0>(get_superbubble_terminus(i + 1))) {
                assert(!superbubble_chain[i * 2]);

                // end of a superbubble chain
                size_t j = i + 1;
                size_t chain_term = i + 1;
                superbubble_chain[(j - 1) * 2] = chain_term;
                superbubble_chain[(j - 1) * 2 + 1] = ++chain_widths_id;
                chain_widths.emplace_back(0);
                chain_bounds.emplace_back(true);

                tsl::hopscotch_set<size_t> widths;
                widths.emplace(0);
                while (chain_parent[j]) {
                    size_t old_j = j;
                    std::ignore = old_j;
                    j = chain_parent[j];
                    size_t widths_start = chain_widths.size();
                    superbubble_chain[(j - 1) * 2] = chain_term;
                    superbubble_chain[(j - 1) * 2 + 1] = ++chain_widths_id;
                    auto [t, d_begin, d_end] = get_superbubble_terminus(j);
                    assert(t == old_j);
                    tsl::hopscotch_set<size_t> next_widths;
                    for (size_t width : widths) {
                        std::for_each(d_begin, d_end, [&](size_t dd) {
                            next_widths.emplace(dd + width);
                        });
                    }
                    std::swap(next_widths, widths);
                    assert(widths.size());
                    for (size_t width : widths) {
                        chain_widths.emplace_back(width);
                        chain_bounds.emplace_back(false);
                    }
                    std::sort(chain_widths.end() - widths.size(), chain_widths.end());
                    chain_bounds[widths_start] = true;
                }

                ++chain_i;
            }
        }

        chain_bounds.emplace_back(true);
        assert(chain_bounds.size() == chain_widths.size() + 1);

        logger->info("Indexed {} superbubble chains containing {} superbubbles",
                     chain_i, chain_widths_id);
        unitig_chain_ = SuperbubbleStorage(std::move(superbubble_chain));

        sdsl::int_vector<> chain_widths_sd(chain_widths.size());
        std::copy(chain_widths.begin(), chain_widths.end(), chain_widths_sd.begin());
        unitig_chain_sizes_ = SuperbubbleStorage(std::move(chain_widths_sd));

        unitig_chain_b_ = SuperbubbleIndicator(to_sdsl(std::move(chain_bounds)));
        assert(unitig_chain_b_.num_set_bits() == chain_widths_id + 1);
    }

    assert(annotator.num_labels() <= 1);
    assert(std::adjacent_find(boundaries.begin(), boundaries.end()) == boundaries.end());

    path_boundaries_ = bit_vector_smart([&](const auto &callback) {
        std::for_each(boundaries.begin(), boundaries.end(), callback);
    }, boundaries.back() + 1, boundaries.size());

    std::string graph_fname = graph_name;
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_col");
    std::string out_path = tmp_dir/"test_col";

    if constexpr(std::is_same_v<PathStorage, RowDiffCoordAnnotator::binary_matrix_type>) {
        if (graph_fname.empty()) {
            graph->serialize(out_path);
            graph_fname = out_path + graph->file_extension();
        }
    }

    size_t num_columns = anno_graph.get_annotator().get_label_encoder().size();
    paths_indices_ = compress_path_storage<PathStorage>(
        graph.get(),
        std::move(annotator),
        label_encoder,
        num_columns,
        graph_fname,
        tmp_dir,
        out_path
    );

    set_graph(graph);

    if (num_unitigs_ > 1) {
        AnnotatedDBG anno_graph(std::const_pointer_cast<DBGSuccinct>(graph),
                                std::make_unique<ColumnCompressed<>>(dbg_succ.max_index()));

        auto &annotator = const_cast<ColumnCompressed<>&>(
            static_cast<const ColumnCompressed<>&>(anno_graph.get_annotator())
        );
        std::atomic<size_t> seq_count { 0 };
        size_t total_seq_count = 0;
        std::vector<uint64_t> boundaries { 0 };

        logger->trace("Indexing extra sequences");
        ThreadPool thread_pool(get_num_threads());
        std::mutex mu;

        generate_sequences([&](std::string_view seq) {
            ++total_seq_count;
            if (total_seq_count && !(total_seq_count % 1000)) {
                logger->trace("Processed {} sequences, indexed {} of them",
                              total_seq_count, seq_count);
            }

            thread_pool.enqueue([&,s=std::string(seq)]() {
                auto nodes = map_to_nodes(*check_graph, s);
                std::vector<Row> rows;
                sdsl::bit_vector picked(nodes.size(), false);
                bool any_picked = false;
                for (size_t i = 0; i < nodes.size(); ++i) {
                    node_index node = nodes[i];
                    if (has_coord(node)) {
                        picked[i] = true;
                        any_picked = true;
                        rows.emplace_back(AnnotatedDBG::graph_to_anno_index(node));
                    }
                }

                if (any_picked) {
                    bool index_read = false;

                    tsl::hopscotch_set<size_t> unitig_ids;
                    for (const auto &tuples : get_path_row_tuples(rows)) {
                        for (const auto &[dummy, tuple] : tuples) {
                            assert(dummy == 0);
                            assert(tuple.size());
                            size_t unitig_id = coord_to_path_id(tuple[0]);
                            if (!std::get<0>(get_superbubble_terminus(unitig_id))
                                    && !std::get<0>(get_superbubble_and_dist(unitig_id))) {
                                index_read = true;
                            }

                            if (unitig_id)
                                unitig_ids.insert(unitig_id);
                        }
                    }

                    if (index_read && unitig_ids.size() > 1) {
                        ++seq_count;
                        std::lock_guard<std::mutex> lock(mu);
                        uint64_t coord = boundaries.back();
                        annotator.add_labels(rows, DUMMY);
                        auto it = rows.begin();
                        for (size_t i = 0; i < nodes.size(); ++i) {
                            if (picked[i])
                                annotator.add_label_coord(*it, DUMMY, coord);

                            ++coord;
                            ++it;
                        }
                        boundaries.emplace_back(coord);
                    }
                }
            });
        });

        thread_pool.join();

        if (!total_seq_count)
            return;

        logger->info("Indexed {} / {} sequences", seq_count, total_seq_count);

        if (!seq_count)
            return;

        read_boundaries_ = bit_vector_smart([&](const auto &callback) {
            std::for_each(boundaries.begin(), boundaries.end(), callback);
        }, boundaries.back() + 1, boundaries.size());

        size_t num_columns = anno_graph.get_annotator().get_label_encoder().size();
        read_indices_ = compress_path_storage<PathStorage>(
            graph.get(),
            std::move(annotator),
            label_encoder,
            num_columns,
            graph_fname,
            tmp_dir,
            out_path
        );
        set_graph(graph);
    }
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
    auto row_tuples = get_path_row_tuples(rows);
    auto read_tuples = get_read_row_tuples(rows);
    std::vector<RowTuples> ret_val;
    ret_val.reserve(nodes.size());
    while (it != picked.end() && !*it) {
        ret_val.emplace_back();
        ++it;
    }

    auto jt = read_tuples.begin();
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
        for (const auto &[c, tuple] : *jt) {
            assert(!c);
            assert(std::adjacent_find(tuple.begin(), tuple.end()) == tuple.end());
            for (auto coord : tuple) {
                size_t read_id = coord_to_read_id(coord);
                out_tuples[read_id].emplace_back(coord);
                path_id_to_nodes[read_id].emplace_back(it - picked.begin());
            }
        }
        ret_val.emplace_back(out_tuples.values_container().begin(),
                             out_tuples.values_container().end());
        ++it;
        while (it != picked.end() && !*it) {
            ret_val.emplace_back();
            ++it;
        }
        ++jt;
    }

    return ret_val;
}

template class PathIndex<>;

} // namespace mtg::graph
