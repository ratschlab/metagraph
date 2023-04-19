#ifndef __TUPLE_ROW_DIFF_HPP__
#define __TUPLE_ROW_DIFF_HPP__

#include <algorithm>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/vector.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

template <class BaseMatrix>
class TupleRowDiff : public binmat::IRowDiff, public MultiIntMatrix {
  public:
    static_assert(std::is_convertible<BaseMatrix*, MultiIntMatrix*>::value);
    static const int SHIFT = 1; // coordinates increase by 1 at each edge

    TupleRowDiff() {}

    TupleRowDiff(const graph::DBGSuccinct *graph, BaseMatrix&& diff)
        : diffs_(std::move(diff)) { graph_ = graph; }

    bool get(Row i, Column j) const override;
    std::vector<Row> get_column(Column j) const override;
    RowTuples get_row_tuples(Row i) const override;
    std::vector<RowTuples> get_row_tuples(const std::vector<Row> &rows) const override;

    /** Returns all labeled traces that pass through a given row.
     * 
     * @param i Index of the row.
     * 
     * @return Vector of pairs (path, column), where path is
     * a vector of Row indices and column is a corresponding Label index.
     */
    std::vector<std::tuple<std::vector<Row>, Column, uint64_t>> get_traces_with_row(Row i) const;

    std::vector<std::tuple<std::vector<Row>, Column, uint64_t>> get_traces_with_row(Row i, std::unordered_map<Row, RowTuples> &rows_annotations, 
    std::vector<Row> &rd_ids, VectorMap<Row, size_t> &node_to_rd, std::vector<RowTuples> &rd_rows, std::unordered_map<Column, std::vector<Row>> &reconstructed_reads) const;
    std::vector<std::vector<std::tuple<std::vector<Row>, Column, uint64_t>>> get_traces_with_row(std::vector<Row> i) const;

    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_attributes() const override { return diffs_.num_attributes(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

  private:
    static void decode_diffs(RowTuples *diffs);
    static void add_diff(const RowTuples &diff, RowTuples *row);

    BaseMatrix diffs_;
};


template <class BaseMatrix>
bool TupleRowDiff<BaseMatrix>::get(Row i, Column j) const {
    SetBitPositions set_bits = get_row(i);
    auto v = std::lower_bound(set_bits.begin(), set_bits.end(), j);
    return v != set_bits.end() && *v == j;
}

template <class BaseMatrix>
std::vector<MultiIntMatrix::Row> TupleRowDiff<BaseMatrix>::get_column(Column j) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    const graph::boss::BOSS &boss = graph_->get_boss();
    assert(!fork_succ_.size() || fork_succ_.size() == boss.get_last().size());

    // TODO: implement a more efficient algorithm
    std::vector<Row> result;
    for (Row i = 0; i < num_rows(); ++i) {
        auto edge = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(i)
        );

        if (boss.get_W(edge) && get(i, j))
            result.push_back(i);
    }
    return result;
}

template <class BaseMatrix>
MultiIntMatrix::RowTuples TupleRowDiff<BaseMatrix>::get_row_tuples(Row row) const {
    return get_row_tuples(std::vector<Row>{ row })[0];
}

template <class BaseMatrix>
std::vector<MultiIntMatrix::RowTuples>
TupleRowDiff<BaseMatrix>::get_row_tuples(const std::vector<Row> &row_ids) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->get_boss().get_last().size());

    // get row-diff paths
    auto [rd_ids, rd_paths_trunc] = get_rd_ids(row_ids);

    std::vector<RowTuples> rd_rows = diffs_.get_row_tuples(rd_ids);
    for (auto &row : rd_rows) {
        decode_diffs(&row);
    }

    rd_ids = std::vector<Row>();

    // reconstruct annotation rows from row-diff
    std::vector<RowTuples> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        RowTuples &result = rows[i];

        auto it = rd_paths_trunc[i].rbegin();
        std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
        result = rd_rows[*it];
        // propagate back and reconstruct full annotations for predecessors
        for (++it ; it != rd_paths_trunc[i].rend(); ++it) {
            std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
            add_diff(rd_rows[*it], &result);
            // replace diff row with full reconstructed annotation
            rd_rows[*it] = result;
        }
        assert(std::all_of(result.begin(), result.end(),
                           [](auto &p) { return p.second.size(); }));
    }

    return rows;
}

template <class BaseMatrix>
bool TupleRowDiff<BaseMatrix>::load(std::istream &in) {
    std::string version(4, '\0');
    in.read(version.data(), 4);
    return anchor_.load(in) && fork_succ_.load(in) && diffs_.load(in);
}

template <class BaseMatrix>
void TupleRowDiff<BaseMatrix>::serialize(std::ostream &out) const {
    out.write("v2.0", 4);
    anchor_.serialize(out);
    fork_succ_.serialize(out);
    diffs_.serialize(out);
}

template <class BaseMatrix>
void TupleRowDiff<BaseMatrix>::decode_diffs(RowTuples *diffs) {
    std::ignore = diffs;
    // no encoding
}

template <class BaseMatrix>
void TupleRowDiff<BaseMatrix>::add_diff(const RowTuples &diff, RowTuples *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));

    if (diff.size()) {
        RowTuples result;
        result.reserve(row->size() + diff.size());

        auto it = row->begin();
        auto it2 = diff.begin();
        while (it != row->end() && it2 != diff.end()) {
            if (it->first < it2->first) {
                result.push_back(*it);
                ++it;
            } else if (it->first > it2->first) {
                result.push_back(*it2);
                ++it2;
            } else {
                if (it2->second.size()) {
                    result.emplace_back(it->first, Tuple{});
                    std::set_symmetric_difference(it->second.begin(), it->second.end(),
                                                  it2->second.begin(), it2->second.end(),
                                                  std::back_inserter(result.back().second));
                }
                ++it;
                ++it2;
            }
        }
        std::copy(it, row->end(), std::back_inserter(result));
        std::copy(it2, diff.end(), std::back_inserter(result));

        row->swap(result);
    }

    assert(std::is_sorted(row->begin(), row->end()));
    for (auto &[j, tuple] : *row) {
        assert(std::is_sorted(tuple.begin(), tuple.end()));
        for (uint64_t &c : tuple) {
            c -= SHIFT;
        }
    }
}

template <class BaseMatrix>
std::vector<std::vector<std::tuple<std::vector<MultiIntMatrix::Row>, MultiIntMatrix::Column, uint64_t>>> TupleRowDiff<BaseMatrix>
::get_traces_with_row(std::vector<Row> i) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->get_boss().get_last().size());

    // a map that stores decompressed annotations for the rows
    std::unordered_map<Row, RowTuples> rows_annotations;

    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;

    std::vector<RowTuples> rd_rows;
    std::unordered_map<Column, std::vector<Row>> reconstructed_reads;

    std::vector<std::vector<std::tuple<std::vector<Row>, Column, uint64_t>>> result;

    for (Row &row_i : i) {
        auto curr_res = get_traces_with_row(row_i, rows_annotations, rd_ids, node_to_rd, rd_rows, reconstructed_reads);
        result.push_back(curr_res);
    }

    return result;
}


template <class BaseMatrix>
std::vector<std::tuple<std::vector<MultiIntMatrix::Row>, MultiIntMatrix::Column, uint64_t>> TupleRowDiff<BaseMatrix>
::get_traces_with_row(Row i) const {

    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->get_boss().get_last().size());

    // a map that stores decompressed annotations for the rows
    std::unordered_map<Row, RowTuples> rows_annotations;

    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;

    std::vector<RowTuples> rd_rows;
    std::unordered_map<Column, std::vector<Row>> reconstructed_reads;

    return get_traces_with_row(i, rows_annotations, rd_ids, node_to_rd, rd_rows, reconstructed_reads);
}

template <class BaseMatrix>
std::vector<std::tuple<std::vector<MultiIntMatrix::Row>, MultiIntMatrix::Column, uint64_t>> TupleRowDiff<BaseMatrix>
::get_traces_with_row(Row i, std::unordered_map<Row, RowTuples> &rows_annotations, std::vector<Row> &rd_ids, VectorMap<Row, size_t> &node_to_rd, 
std::vector<RowTuples> &rd_rows, std::unordered_map<Column, std::vector<Row>> &reconstructed_reads) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->get_boss().get_last().size());    

    std::ignore = reconstructed_reads;

    const graph::boss::BOSS &boss = graph_->get_boss();
    const bit_vector &rd_succ = fork_succ_.size() ? fork_succ_ : boss.get_last();

    std::unordered_map<Column, std::set<uint64_t>> discovered_coordinates;

    std::unordered_map<Column, std::map<uint64_t, Row>> test_new_path_reconstruction;

    // part 1. Graph traversal

    // step 1. Build row-diff path for the input Row
    RowTuples input_row_row_tuples;

    if (!rows_annotations.count(i)) {
        std::vector<size_t> rd_path_trunc; 
        std::vector<Row> rd_ids_start;

        graph::boss::BOSS::edge_index boss_edge = graph_->kmer_to_boss_index(
                    graph::AnnotatedSequenceGraph::anno_to_graph_index(i));

        while (true) {
            Row row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                    graph_->boss_to_kmer_index(boss_edge));

            auto [it, is_new] = node_to_rd.try_emplace(row, rd_ids.size());
            rd_path_trunc.push_back(it.value());

            if (!is_new)
                break;

            rd_ids.push_back(row);
            rd_ids_start.push_back(row); 

            if (anchor_[row])
                break;

            boss_edge = boss.row_diff_successor(boss_edge, rd_succ);
        }

        // step 2.1 Decompress annotations for the input Row 
        std::vector<RowTuples> rd_rows_start = diffs_.get_row_tuples(rd_ids_start);
        for (auto &row : rd_rows_start) {
            // no encoding
            decode_diffs(&row);
        }

        rd_rows.insert(rd_rows.end(), rd_rows_start.begin(), rd_rows_start.end());

        RowTuples result_of_decompression;     
        auto it = rd_path_trunc.rbegin();
        std::sort(rd_rows[*it].begin(), rd_rows[*it].end());

        result_of_decompression = rd_rows[*it];
        rows_annotations[rd_ids[*it]] = result_of_decompression;

        // propagate back and reconstruct full annotations for predecessors
        for (++it ; it != rd_path_trunc.rend(); ++it) {
            std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
            add_diff(rd_rows[*it], &result_of_decompression);
            // replace diff row with full reconstructed annotation
            rd_rows[*it] = result_of_decompression;

            // keep the decompressed annotations for the Rows
            // along that row-diff path
            rows_annotations[rd_ids[*it]] = result_of_decompression;
        }
        assert(std::all_of(result_of_decompression.begin(), result_of_decompression.end(),
                            [](auto &p) { return p.second.size(); }));
        
        
        
        input_row_row_tuples = result_of_decompression;
    } else {
        // if we already decompressed the annotations for this node
        // then simply retrieve them from the map
        input_row_row_tuples = rows_annotations[i];
    }

    

    // step 2.2 Save labels of the input Row in a set
    std::unordered_set<Column> labels_of_the_start_node;

    // step 2.3 Set up adjacency lists for each label of the input Row

    // the adjacency lists for the subgraph traversed forwards
    // std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_parents_map;
    // std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_children_map;

    // the adjacency lists for the subgraph traversed backwards
    // std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_parents_map_backward;
    // std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_children_map_backward;

    std::vector<std::tuple<std::vector<Row>, Column, uint64_t>> result;

    std::unordered_map<Column, uint64_t> input_row_position_in_ref_seq;

    for (auto &[j, tuple] : input_row_row_tuples) {

        // if (!reconstructed_reads.count(j)) {
        //     labels_of_the_start_node.insert(j);

        //     // initialize adjacency lists for each label
        //     labeled_parents_map[j];
        //     labeled_children_map[j];

        //     labeled_parents_map_backward[j];
        //     labeled_children_map_backward[j];

        //     input_row_position_in_ref_seq[j] = tuple.front();
        // } else {
        //     // return this as well
        //     result.push_back(std::make_tuple(reconstructed_reads[j], j, tuple.front()));
        // }
        

        labels_of_the_start_node.insert(j);

        // initialize adjacency lists for each label
        // labeled_parents_map[j];
        // labeled_children_map[j];

        // labeled_parents_map_backward[j];
        // labeled_children_map_backward[j];

        input_row_position_in_ref_seq[j] = tuple.front();
    }

    // std::cout << "labels of the start node:";
    // for (auto & lab_st : labels_of_the_start_node) {
    //     std::cout << lab_st << ", ";
    // }

    // std::cout << '\n';
    // if (labels_of_the_start_node.empty()) {
    //     return result;
    // }


    // step 3. Traverse the graph using DFS 
    // and fill adjacency lists along the way
    std::vector<Row> to_visit;
    to_visit.push_back(i);

    std::unordered_set<Row> visited_nodes_forward;
    while (!to_visit.empty()) {
        // pick a node from the stack
        Row current_parent = to_visit.back();
        to_visit.pop_back();
    
        // if it was already visited, continue
        if (visited_nodes_forward.count(current_parent))
            continue;

        // collect the current node's labels into a set
        std::unordered_set<Column> current_parent_labels;
        for (auto &rowt_parent : rows_annotations[current_parent]) {
            current_parent_labels.insert(rowt_parent.first);
        }
        
        visited_nodes_forward.insert(current_parent);

        graph::AnnotatedSequenceGraph::node_index current_parent_to_graph = graph::AnnotatedSequenceGraph::anno_to_graph_index(current_parent);
        
        // for each adjacent outgoing k-mer 
        graph_->call_outgoing_kmers(current_parent_to_graph, [&](auto next, char c) {
            // std::ignore = c;
            if (c == graph::boss::BOSS::kSentinel)
                return;
            
            Row next_to_anno = graph::AnnotatedSequenceGraph::graph_to_anno_index(next);
            RowTuples next_annotations;

            // if the annotations for the child node are not decompressed,
            // build rd-path and decompress the annotations
            // (repeat steps 1-2.1)
            if (!rows_annotations.count(next_to_anno)) {
                
                // commented last
                std::vector<size_t> rd_path_next; 
                std::vector<Row> rd_ids_next;

                graph::boss::BOSS::edge_index boss_edge = graph_->kmer_to_boss_index(next);

                while (true) {
                    Row row_w = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                            graph_->boss_to_kmer_index(boss_edge));

                    auto [it, is_new] = node_to_rd.try_emplace(row_w, rd_ids.size());
                    rd_path_next.push_back(it.value());

                    if (!is_new)
                        break;                        

                    rd_ids.push_back(row_w);
                    rd_ids_next.push_back(row_w);  

                    if (anchor_[row_w])
                        break;

                    boss_edge = boss.row_diff_successor(boss_edge, rd_succ);
                }

                // get the RowTuples for the new rows along the rd-path
                std::vector<RowTuples> rd_rows_next = diffs_.get_row_tuples(rd_ids_next);
              
                for (auto &row : rd_rows_next) {
                    decode_diffs(&row);
                }

                // add the new rows to the vector of all rows
                rd_rows.insert(rd_rows.end(), rd_rows_next.begin(), rd_rows_next.end());
                
                RowTuples result_of_decompression_next; 

                auto it = rd_path_next.rbegin();
                std::sort(rd_rows[*it].begin(), rd_rows[*it].end());

                result_of_decompression_next = rd_rows[*it];
                rows_annotations[rd_ids[*it]] = result_of_decompression_next;

                // propagate back and reconstruct full annotations for predecessors
                for (++it ; it != rd_path_next.rend(); ++it) {
                    std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
                    add_diff(rd_rows[*it], &result_of_decompression_next);
                    // replace diff row with full reconstructed annotation
                    rd_rows[*it] = result_of_decompression_next;
                    
                    rows_annotations[rd_ids[*it]] = result_of_decompression_next;
                }

                next_annotations = result_of_decompression_next;
            } else {
                // if we already decompressed the annotations for this node
                // then simply retrieve them from the map
                next_annotations = rows_annotations[next_to_anno];
            }

            // acceptable node is a node whose labels match
            // the labels of the start node
            bool acceptable_node = false;
            
            for (auto &[j_next, tuple_next] : next_annotations) {
                // check if the label matches the start row labels
                // and add the child-parent relation to the corresponding map
                if (labels_of_the_start_node.count(j_next) && current_parent_labels.count(j_next)) {
                    acceptable_node = true;
                    // labeled_parents_map[j_next][next_to_anno].push_back(current_parent);
                    // labeled_children_map[j_next][current_parent].push_back(next_to_anno);
                }
            }

            // if the row has at least one label that matches
            // the start row's labels then add it to the stack
            if (acceptable_node) {
                to_visit.push_back(next_to_anno);
            }

        } );
    }

    // DFS backwards is similar to the DFS forwards
    to_visit.clear();
    to_visit.push_back(i);


    std::unordered_set<Row> visited_nodes_backward;
    while (!to_visit.empty()) {
        Row current_child = to_visit.back();
        to_visit.pop_back();

        if (visited_nodes_backward.count(current_child))
            continue;

        std::unordered_set<Column> current_child_labels;
        for (auto &rowt_parent : rows_annotations[current_child]) {
            current_child_labels.insert(rowt_parent.first);
        }
        
        visited_nodes_backward.insert(current_child);

        graph::AnnotatedSequenceGraph::node_index current_child_to_graph = graph::AnnotatedSequenceGraph::anno_to_graph_index(current_child);

        graph_->call_incoming_kmers(current_child_to_graph, [&](auto previous, char c) {

            if (c == graph::boss::BOSS::kSentinel)
                return;
            
            Row previous_to_anno = graph::AnnotatedSequenceGraph::graph_to_anno_index(previous);
            RowTuples previous_annotations;

            if (!rows_annotations.count(previous_to_anno)) {

                // commented last
                std::vector<size_t> rd_path_previous; 
                std::vector<Row> rd_ids_previous;

                graph::boss::BOSS::edge_index boss_edge = graph_->kmer_to_boss_index(previous);

                while (true) {
                    Row row_w = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                            graph_->boss_to_kmer_index(boss_edge));

                    auto [it, is_new] = node_to_rd.try_emplace(row_w, rd_ids.size());
                    rd_path_previous.push_back(it.value());

                    if (!is_new)
                        break;

                    rd_ids.push_back(row_w);
                    rd_ids_previous.push_back(row_w);  

                    if (anchor_[row_w])
                        break;

                    boss_edge = boss.row_diff_successor(boss_edge, rd_succ);
                }

                std::vector<RowTuples> rd_rows_previous = diffs_.get_row_tuples(rd_ids_previous);

                for (auto &row : rd_rows_previous) {
                    decode_diffs(&row);
                }
                
                rd_rows.insert(rd_rows.end(), rd_rows_previous.begin(), rd_rows_previous.end());
                
                RowTuples result_of_decompression_previous;
                auto it = rd_path_previous.rbegin();

                std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
                result_of_decompression_previous = rd_rows[*it];
                
                rows_annotations[rd_ids[*it]] = result_of_decompression_previous;

                for (++it ; it != rd_path_previous.rend(); ++it) {
                    std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
                    add_diff(rd_rows[*it], &result_of_decompression_previous);
                    rd_rows[*it] = result_of_decompression_previous;
                    
                    rows_annotations[rd_ids[*it]] = result_of_decompression_previous;
                }

                previous_annotations = result_of_decompression_previous;
            } else {
                previous_annotations = rows_annotations[previous_to_anno];
            }

            // aceptable node is a node whose labels match
            // the labels of the start node
            bool acceptable_node = false;

            for (auto &[j_prev, tuple_prev] : previous_annotations) {
                if (labels_of_the_start_node.count(j_prev) && current_child_labels.count(j_prev)) {
                    acceptable_node = true;
                    // labeled_parents_map_backward[j_prev][previous_to_anno].push_back(current_child);
                    // labeled_children_map_backward[j_prev][current_child].push_back(previous_to_anno);
                }
            }
            if (acceptable_node) {
                to_visit.push_back(previous_to_anno);
            }
                        
        } );

    }

    // separate rows' annotations by labels
    // std::unordered_map<Column, std::unordered_map<Row, std::vector<uint64_t>>> annotations_map_sep_by_labels;
    for (auto & [row, row_tuples] : rows_annotations) {
        for (auto & [j, tuple] : row_tuples) {
            if (!labels_of_the_start_node.count(j)) {
                continue;
            }
            for (uint64_t &c : tuple) {
                discovered_coordinates[j].insert(c);
                test_new_path_reconstruction[j][c] = row;

                // annotations_map_sep_by_labels[j][row].push_back(c);
            }
        }
    }
    
    // part 2. Paths reconstruction

    // // each labeled sequence can have only one target (that is the last kmer in the sequence)
    // std::unordered_map<Column, Row> labeled_targets;

    // // find the ends of the sequences
    // for (auto & [j, children_map] : labeled_children_map) {
    //     if (!children_map.size()) {
    //         continue;
    //     }        

    //     // we can pick any key node from children_map to follow the edges and find the target
    //     uint64_t final_coord = annotations_map_sep_by_labels[j][children_map.begin()->first].back();
    //     Row final_target = children_map.begin()->first;
    //     Row current_row = final_target;

    //     // while current_row has children
    //     while (children_map.count(current_row)) {
    //         std::vector<Row> current_children = children_map[current_row];

    //         // the next coordinate must increase by 1
    //         uint64_t final_coord_candidate = final_coord + SHIFT;
    //         bool coordinates_can_increase_more = false;

    //         // find an adjacent outgoing node with the next coordinate
    //         for (Row & current_child_candidate : current_children) {
    //             std::vector<uint64_t> curr_cand_coordinates = annotations_map_sep_by_labels[j][current_child_candidate];

    //             // convert the vector of coordinates into set to check if it contains final_coord_candindate
    //             std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(),
    //                                                                    curr_cand_coordinates.end());

    //             // check if set of coordinates for the current child node contains increased coordinate
    //             if (curr_cand_coordinates_set.count(final_coord_candidate)) {
    //                 final_coord = final_coord_candidate;
    //                 current_row = current_child_candidate;
    //                 coordinates_can_increase_more = true;
    //                 final_target = current_row;
    //                 break;
    //             }
    //         }

    //         // if haven't found the child node with the increased coordinate, then current node is the target
    //         if (!coordinates_can_increase_more) {
    //             final_target = current_row;
    //             break;
    //         }
    //     }
    //     labeled_targets[j] = final_target;
    // }

    // // for each label trace the forward path 
    // std::unordered_map<Column, std::vector<std::pair<Row, uint64_t>>> traces_forward_sep_by_labels;
    // for (auto & [target_label, target_row] : labeled_targets) {
    //     Row current_row = target_row;
    //     std::vector<std::pair<Row, uint64_t>> current_trace;
        
    //     // get the very last coordinate of the target node (that is the biggest one)
    //     // question: is the biggest coordinate in annotations always the last one?
    //     uint64_t current_coord = annotations_map_sep_by_labels[target_label][current_row].back(); 
    //     current_trace.push_back(std::make_pair(current_row, current_coord));

    //     // start traversing the incoming edges and trace the path by checking the coordinates decrease
    //     while (labeled_parents_map[target_label].count(current_row)) {                
    //         std::vector<Row> current_parents = labeled_parents_map[target_label][current_row];

    //         uint64_t current_coord_candidate = current_coord - SHIFT;
    //         bool coordinates_can_decrease_more = false;

    //         // find the previous row in path with the decreased coordinate
    //         for (Row & current_row_candidate : current_parents) {
    //             std::vector<uint64_t> curr_cand_coordinates = annotations_map_sep_by_labels[target_label][current_row_candidate];
    //             std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(),
    //                                                                    curr_cand_coordinates.end());

    //             if (curr_cand_coordinates_set.count(current_coord_candidate)) {
    //                 current_coord = current_coord_candidate;
    //                 current_row = current_row_candidate;
    //                 coordinates_can_decrease_more = true;
    //                 break;
    //             }
    //         }

    //         // if haven't found the parent node with the decreased coordinate, then current node is the start of the path
    //         if (!coordinates_can_decrease_more)   
    //             break;                
            
    //         current_trace.push_back(std::make_pair(current_row, current_coord));                     
    //     }

    //     std::reverse(current_trace.begin(), current_trace.end());
    //     traces_forward_sep_by_labels[target_label] = current_trace;
    // }

    // // similarly, the backward traversal targets are found
    // std::unordered_map<Column, Row> labeled_targets_backward;
    // for (auto & [j, children_map] : labeled_children_map_backward) {
    //     if (!children_map.size()) {
    //         continue;
    //     }

    //     uint64_t final_coord = annotations_map_sep_by_labels[j][children_map.begin()->first].front();
    //     Row final_target = children_map.begin()->first;
    //     Row current_row = final_target;

    //     while (children_map.count(current_row)) {
    //         std::vector<Row> current_children = children_map[current_row];

    //         uint64_t final_coord_candidate = final_coord - SHIFT;            
    //         bool coordinates_can_decrease_more = false;

    //         for (Row & current_child_candidate : current_children) {
    //             std::vector<uint64_t> curr_cand_coordinates = annotations_map_sep_by_labels[j][current_child_candidate];
    //             std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(),
    //                                                                    curr_cand_coordinates.end());

    //             if (curr_cand_coordinates_set.count(final_coord_candidate)) {
    //                 final_coord = final_coord_candidate;
    //                 current_row = current_child_candidate;
    //                 coordinates_can_decrease_more = true;
    //                 final_target = current_row;
    //                 break;
    //             }
    //         }

    //         if (!coordinates_can_decrease_more) {
    //             final_target = current_row;
    //             break;
    //         }
    //     }
    //     labeled_targets_backward[j] = final_target;
    // }   

    // // similarly, for each label the backward path is traced
    // std::unordered_map<Column, std::vector<std::pair<Row, uint64_t>>> traces_backward_sep_by_labels;
    // for (auto & [target_label, target_row] : labeled_targets_backward) {
    //     Row current_row = target_row;
    //     std::vector<std::pair<Row, uint64_t>> current_trace;

    //     uint64_t current_coord = annotations_map_sep_by_labels[target_label][current_row].front(); 
    //     current_trace.push_back(std::make_pair(current_row, current_coord));
        
    //     while (labeled_parents_map_backward[target_label].count(current_row)) { 
    //         std::vector<Row> current_parents = labeled_parents_map_backward[target_label][current_row];

    //         uint64_t curr_coord_candidate = current_coord + SHIFT;
    //         bool coordinates_can_increase_more = false;

    //         for (Row & current_row_candidate : current_parents) {
    //             std::vector<uint64_t> curr_cand_coordinates = annotations_map_sep_by_labels[target_label][current_row_candidate];
    //             std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(),
    //                                                                    curr_cand_coordinates.end());

    //             if (curr_cand_coordinates_set.count(curr_coord_candidate)) {
    //                 current_coord = curr_coord_candidate;
    //                 current_row = current_row_candidate;
    //                 coordinates_can_increase_more = true;
    //                 break;
    //             }
    //         }

    //         if (!coordinates_can_increase_more)
    //             break;            
            
    //         current_trace.push_back(std::make_pair(current_row, current_coord));          
    //     }

    //     traces_backward_sep_by_labels[target_label] = current_trace;

    // }

    // new path reconstruction for graphs with unlabeled reads

    std::unordered_map<Column, std::vector<std::pair<Row, uint64_t>>> traces_unlabeled;

    for (auto & [j, coords_map] : test_new_path_reconstruction) {
        for (auto & [ccoord, rrow] : coords_map) {
            traces_unlabeled[j].push_back(std::make_pair(rrow, ccoord));
        }
    }

    // extact reads from walk that has jumps in coordinates

    std::vector<std::tuple<std::vector<std::pair<Row, uint64_t>>, Column, uint64_t>> result_coords;

    for (auto & [j, trace_with_jumps] : traces_unlabeled) {

        uint64_t curr_read_start_coord = trace_with_jumps.front().second;
        uint64_t curr_read_curr_coord = curr_read_start_coord;

        std::vector<Row> curr_read_trace;
        curr_read_trace.push_back(trace_with_jumps.front().first);

        auto curr_edge_node = graph::AnnotatedSequenceGraph::anno_to_graph_index(trace_with_jumps.front().first);


        std::vector<std::pair<Row, uint64_t>> curr_read_trace_coords;
        curr_read_trace_coords.push_back(trace_with_jumps.front());

        for (size_t cur_pos = 1; cur_pos < trace_with_jumps.size(); ++cur_pos) {

            auto curr_next_edge_node = graph::AnnotatedSequenceGraph::anno_to_graph_index(trace_with_jumps[cur_pos].first);
            
            bool edge_exists = false;

            graph_->adjacent_outgoing_nodes(curr_edge_node, [&](auto adj_outg_node) {
                if (adj_outg_node == curr_next_edge_node) {
                    edge_exists = true;
                    return;
                }
            });


            if ((trace_with_jumps[cur_pos].second == (curr_read_curr_coord + SHIFT)) && edge_exists) {
                curr_read_trace.push_back(trace_with_jumps[cur_pos].first);
                curr_read_trace_coords.push_back(trace_with_jumps[cur_pos]);

                curr_read_curr_coord += SHIFT;                

            } else {
                result.push_back(std::make_tuple(curr_read_trace, j, curr_read_start_coord));

                result_coords.push_back(std::make_tuple(curr_read_trace_coords, j, curr_read_start_coord));

                curr_read_start_coord = trace_with_jumps[cur_pos].second;
                curr_read_curr_coord = curr_read_start_coord;

                curr_read_trace.clear();
                curr_read_trace.push_back(trace_with_jumps[cur_pos].first);

                curr_read_trace_coords.clear();
                curr_read_trace_coords.push_back(trace_with_jumps[cur_pos]);
            }

            curr_edge_node = curr_next_edge_node;
        }
    }


    std::vector<std::tuple<std::vector<std::pair<Row, uint64_t>>, Column, uint64_t>> final_result;
    
    for (auto & [read_trace, read_label, read_start_coord] : result_coords) { 
        for (size_t cur_pos = 0; cur_pos < read_trace.size(); ++cur_pos) {
            if (read_trace[cur_pos].first == i) {
                final_result.push_back(std::make_tuple(read_trace, read_label, cur_pos));
                break;
            }
        }
    }

    std::vector<std::tuple<std::vector<Row>, Column, uint64_t>> final_final_result;


    for (auto & [read_trace, read_label, cur_pos] : final_result) {
        std::vector<Row> cur_trace_rows;
        for (auto & el : read_trace) {
            cur_trace_rows.push_back(el.first);
        }

        final_final_result.push_back(std::make_tuple(cur_trace_rows, read_label, cur_pos));
    }

    return final_final_result;

    // // concatenate all labeled paths
    // for (Column const &j : labels_of_the_start_node) {

    //     std::vector<std::pair<Row, uint64_t>> path_to_start = traces_backward_sep_by_labels[j];
    //     std::vector<std::pair<Row, uint64_t>> path_to_end = traces_forward_sep_by_labels[j];

    //     if (path_to_start.empty()) {
            
    //         std::vector<Row> path_to_add_to_result;
    //         path_to_add_to_result.reserve(path_to_end.size());
    //         for (auto & [row, c] : path_to_end) { 
    //             path_to_add_to_result.push_back(row);
    //         }

    //         reconstructed_reads[j] = path_to_add_to_result;
    //         result.push_back(std::make_tuple(path_to_add_to_result, j, input_row_position_in_ref_seq[j]));

    //     } else if (path_to_end.empty()) {

    //         std::vector<Row> path_to_add_to_result;
    //         path_to_add_to_result.reserve(path_to_start.size());
    //         for (auto & [row, c] : path_to_start) { 
    //             path_to_add_to_result.push_back(row);
    //         }

    //         reconstructed_reads[j] = path_to_add_to_result;
    //         result.push_back(std::make_tuple(path_to_add_to_result, j, input_row_position_in_ref_seq[j]));

    //     } else {

    //         // if graph contains cycles then forward and backward paths may overlap
    //         std::vector<std::pair<Row, uint64_t>> starting_path;
    //         std::vector<std::pair<Row, uint64_t>> ending_path;

    //         // find which path has a starting node (i.e. first node with lower coordinate)
    //         if (path_to_start.front().second < path_to_end.front().second)
    //             starting_path = path_to_start;
    //         else
    //             starting_path = path_to_end;

    //         // find which path has an starting node (i.e. last node with bigger coordinate)
    //         if (path_to_start.back().second > path_to_end.back().second)
    //             ending_path = path_to_start;
    //         else
    //             ending_path = path_to_end;

    //         uint64_t last_coord_start_path = starting_path.back().second;
    //         std::vector<Row> path_to_add_to_result;

    //         for (auto it = ending_path.rbegin(); it < ending_path.rend(); ++it) {
    //             if (it->second > last_coord_start_path)
    //                 path_to_add_to_result.insert(path_to_add_to_result.begin(), it->first);
    //             else
    //                 break;
    //         }

    //         std::vector<Row> rows_start_path;
    //         rows_start_path.reserve(starting_path.size());
    //         for (auto & [row, c] : starting_path) {
    //             rows_start_path.push_back(row);
    //         }            

    //         path_to_add_to_result.insert(path_to_add_to_result.begin(),
    //                                      rows_start_path.begin(),
    //                                      rows_start_path.end());

    //         reconstructed_reads[j] = path_to_add_to_result;
    //         result.push_back(std::make_tuple(path_to_add_to_result, j, input_row_position_in_ref_seq[j]));
    //     }
    // }    

    // return result;
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __TUPLE_ROW_DIFF_HPP__
