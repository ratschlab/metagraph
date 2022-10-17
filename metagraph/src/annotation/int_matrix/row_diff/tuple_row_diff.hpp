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
     * @param i Index of the row.
     * @return Vector of pairs (path, column), where path is
     * a vector of Row indices and column is a corresponding Label index.
     */
    std::vector<std::pair<std::vector<Row>, Column>> get_traces_with_row(Row i) const;

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
std::vector<std::pair<std::vector<MultiIntMatrix::Row>, MultiIntMatrix::Column>> TupleRowDiff<BaseMatrix>
::get_traces_with_row(Row i) const {
   
    std::vector<Row> to_visit;
    to_visit.push_back(i);

    // set of rows to fetch annotations for
    std::unordered_set<Row> visited_nodes_to_bottom;

    // during DFS to the bottom of the graph two maps are being constucted: parents_map and children_map
    // children_map is used later to find target nodes (ends of sequences) which are not sink vertices
    // parents_map is used later to trace the paths from start node (Row i) to the found targets
    std::unordered_map<Row, std::vector<Row>> parents_map_to_bottom;
    std::unordered_map<Row, std::vector<Row>> children_map_to_bottom;

    // sink vertices (i.e. vertices with outdegree zero) are the targets which can be found without annotations
    std::vector<Row> sink_targets;
    
    while (!to_visit.empty()) {
        Row current_parent = to_visit.back();
        to_visit.pop_back();
 
        if (visited_nodes_to_bottom.count(current_parent))
            continue;
        
        visited_nodes_to_bottom.insert(current_parent);

        bool node_has_outgoing_kmers = false;

        graph::AnnotatedSequenceGraph::node_index current_parent_to_graph = graph::AnnotatedSequenceGraph::anno_to_graph_index(current_parent);
        
        graph_->adjacent_outgoing_nodes(current_parent_to_graph, [&](auto next) {
            
            Row next_to_anno = graph::AnnotatedSequenceGraph::graph_to_anno_index(next);
            to_visit.push_back(next_to_anno);
            
            if (parents_map_to_bottom.count(next_to_anno)) {
                parents_map_to_bottom[next_to_anno].push_back(current_parent);
            } else {
                parents_map_to_bottom[next_to_anno] = { current_parent };
            }

            if (children_map_to_bottom.count(current_parent)) {
                children_map_to_bottom[current_parent].push_back(next_to_anno);
            } else {
                children_map_to_bottom[current_parent] = { next_to_anno };
            }

            node_has_outgoing_kmers = true;
        } );

        // check if node is a sink vertex (target)
        if (!node_has_outgoing_kmers)
            sink_targets.push_back(current_parent);
    }


    // DFS to the top is similar to dfs to the bottom
    to_visit.clear();
    to_visit.push_back(i);

    // set of nodes to fetch annotations for 
    std::unordered_set<Row> visited_nodes_to_top;

    std::unordered_map<Row, std::vector<Row>> parents_map_to_top;
    std::unordered_map<Row, std::vector<Row>> children_map_to_top;
    
    // source vertices (i.e. vertices with indegree zero) are the targets which can be found without annotations
    std::vector<Row> source_targets;

    while (!to_visit.empty()) {
        Row current_child = to_visit.back();
        to_visit.pop_back();

        if (visited_nodes_to_top.count(current_child))
            continue;
        
        visited_nodes_to_top.insert(current_child);

        bool node_has_incoming_kmers = false;

        graph::AnnotatedSequenceGraph::node_index current_child_to_graph = graph::AnnotatedSequenceGraph::anno_to_graph_index(current_child);

        graph_->adjacent_incoming_nodes(current_child_to_graph, [&](auto previous) {
            
            Row previous_to_anno = graph::AnnotatedSequenceGraph::graph_to_anno_index(previous);
            to_visit.push_back(previous_to_anno);
            
            if (parents_map_to_top.count(previous_to_anno)) {
                parents_map_to_top[previous_to_anno].push_back(current_child);
            } else {
                parents_map_to_top[previous_to_anno] = { current_child };
            }

            if (children_map_to_top.count(current_child)) {
                children_map_to_top[current_child].push_back(previous_to_anno);
            } else {
                children_map_to_top[current_child] = { previous_to_anno };
            }

            node_has_incoming_kmers = true;
        } );

        // check if node is a source vertex (target)
        if (!node_has_incoming_kmers)
            source_targets.push_back(current_child);
    }

    // get set of all nodes for which annotation must be retrieved
    std::unordered_set<Row> all_rows_to_get_annotations(visited_nodes_to_bottom);
    all_rows_to_get_annotations.insert(visited_nodes_to_top.begin(),
                                       visited_nodes_to_top.end());

    // convert this set to vector in order to call method get_traces_with_row()
    std::vector<Row> rows_to_get_annotations_vector;
    rows_to_get_annotations_vector.insert(rows_to_get_annotations_vector.end(),
                                          all_rows_to_get_annotations.begin(),
                                          all_rows_to_get_annotations.end());


    // get annotations for visited nodes
    std::vector<RowTuples> fetched_annotations = get_row_tuples(rows_to_get_annotations_vector);


    // consctruct two map from retrieved annotations
    // fetched_annotations_map is used to construct labeled_parents maps, labeled_children maps, and labeled_targets map
    // fetched_annotations_map_sep_by_labels is used to find target nodes (ends of sequences) which are not sink or source nodes in graph
    std::unordered_map<Row, RowTuples> fetched_annotations_map;
    std::unordered_map<Column, std::unordered_map<Row, std::vector<uint64_t>>> fetched_annotations_map_sep_by_labels;

    // set of labels of the start node is used to trace the paths only with labels that match with the start node (Row i) labels
    std::unordered_set<Column> labels_of_the_start_node;

    for (size_t j = 0; j < rows_to_get_annotations_vector.size(); ++j) {
        fetched_annotations_map[rows_to_get_annotations_vector[j]] = fetched_annotations[j];
        
        // save labels of the start node
        if (rows_to_get_annotations_vector[j] == i) {

            for (auto & rowtuple : fetched_annotations[j]) {
                labels_of_the_start_node.insert(rowtuple.first);
                for (auto & coordinate : rowtuple.second) {
                    fetched_annotations_map_sep_by_labels[rowtuple.first][rows_to_get_annotations_vector[j]].push_back(coordinate);
                }
            }
        } else {
            for (auto & rowtuple : fetched_annotations[j]) {
                for (auto & coordinate : rowtuple.second) {
                    fetched_annotations_map_sep_by_labels[rowtuple.first][rows_to_get_annotations_vector[j]].push_back(coordinate);
                }
            }
        }        
    }

    // create separate parents maps for each label (to trace the paths separately for each label later)
    std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_parents_map;

    for (auto & [child, parents_vector] : parents_map_to_bottom) {

        RowTuples child_annotations = fetched_annotations_map[child];
        
        std::unordered_set<Column> child_annotations_set;

        for (auto & [j, tuple] : child_annotations) {
            child_annotations_set.insert(j);
        }

        for (Row & parent_node : parents_vector) {
            RowTuples parent_annotations = fetched_annotations_map[parent_node];

            for (auto & [parent_j, parent_tuple] : parent_annotations) {

                if (child_annotations_set.count(parent_j) && labels_of_the_start_node.count(parent_j))
                    labeled_parents_map[parent_j][child].push_back(parent_node);
            }
        }
    }


    // same for children map
    std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_children_map;

    for (auto & [parent, children_vector] : children_map_to_bottom) {
        RowTuples parent_annotations = fetched_annotations_map[parent];
        
        std::unordered_set<Column> parent_labels;

        for (auto & [j, tuple] : parent_annotations) {
            parent_labels.insert(j);
        }

        for (Row & child_node : children_vector) {
            RowTuples child_annotations = fetched_annotations_map[child_node];

            for (auto & [child_j, child_tuple] : child_annotations) {
                
                if (parent_labels.count(child_j) && labels_of_the_start_node.count(child_j))
                    labeled_children_map[child_j][parent].push_back(child_node);
            }
        }
    }


    // each labeled sequence can have only one target (that is the last kmer in the sequence)
    std::unordered_map<Column, Row> labeled_targets;

    for (auto & target : sink_targets) {
        RowTuples target_annotations = fetched_annotations_map[target];

        for (auto & [j, tuple] : target_annotations) {
            if (labels_of_the_start_node.count(j))
                labeled_targets[j] = target;
        }
            
    }


    // find addional targets for each label. Those are non-sink nodes that are the ends of the sequences
    std::unordered_map<Column, Row> additional_targets_labeled;

    for (auto & [j, children_map] : labeled_children_map) {

        // no need to find targets for labels that are already found
        if (labeled_targets.count(j))
            continue;

        // we can pick any key node from children_map just to descend to the very bottom and find the target
        uint64_t final_coord = fetched_annotations_map_sep_by_labels[j][children_map.begin()->first].back();
        Row final_target = children_map.begin()->first;

        Row current_row = final_target;

        while (children_map.count(current_row)) {
            std::vector<Row> current_children = children_map[current_row];

            // we are descending in the graph so coordinates must increase
            uint64_t final_coord_candidate = final_coord + SHIFT;
            bool coordinates_can_increase_more = false;

            // find a child node with increased coordinate by one
            for (Row & current_child_candidate : current_children) {
                std::vector<uint64_t> curr_cand_coordinates = fetched_annotations_map_sep_by_labels[j][current_child_candidate];

                // convert to set to check if contains final_coord_candindate
                std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(), curr_cand_coordinates.end());

                if (curr_cand_coordinates_set.count(final_coord_candidate)) {
                    final_coord = final_coord_candidate;
                    current_row = current_child_candidate;
                    coordinates_can_increase_more = true;
                    final_target = current_row;
                    break;
                }
            }

            // if haven't found the children node with next coordinate, then current node is the target
            if (!coordinates_can_increase_more) {
                final_target = current_row;
                break;
            }
        }

        additional_targets_labeled[j] = final_target;  

    }


    // add additional targets to final labeled targets map
    for (auto & [target_label, target_row] : additional_targets_labeled) {
        labeled_targets[target_label] = target_row;
    }


    // for each label trace the path
    std::unordered_map<Column, std::vector<std::pair<Row, uint64_t>>> traces_to_bottom_sep_by_labels;

    for (auto & [target_label, target_row] : labeled_targets) {

        Row current_row = target_row;
        std::vector<std::pair<Row, uint64_t>> current_trace;
        
        // get the very last coordinate of the target node
        uint64_t current_coord = fetched_annotations_map_sep_by_labels[target_label][current_row].back(); 

        current_trace.push_back(std::make_pair(current_row, current_coord));

        // start ascending to the top and trace the path by checking the coordinates decrease
        while (labeled_parents_map[target_label].count(current_row)) {                
            
            std::vector<Row> current_parents = labeled_parents_map[target_label][current_row];

            uint64_t current_coord_candidate = current_coord - SHIFT;
            bool coordinates_can_decrease_more = false;

            // find the previous node in path by checking coordinates decrease
            for (Row & current_row_candidate : current_parents) {
                std::vector<uint64_t> curr_cand_coordinates = fetched_annotations_map_sep_by_labels[target_label][current_row_candidate];
                std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(), curr_cand_coordinates.end());

                if (curr_cand_coordinates_set.count(current_coord_candidate)) {
                    current_coord = current_coord_candidate;
                    current_row = current_row_candidate;
                    coordinates_can_decrease_more = true;
                    break;
                }
            }

            // if haven't found then we reached the top of the path
            if (!coordinates_can_decrease_more)   
                break;                
            
            current_trace.push_back(std::make_pair(current_row, current_coord));
                            
        }

        std::reverse(current_trace.begin(), current_trace.end());
        traces_to_bottom_sep_by_labels[target_label] = current_trace;
    }



    // same operations are being done for parents_map and children_map from the upward traversal of the graph
    std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_parents_map_to_top;

    for (auto & [child, parents_vector] : parents_map_to_top) {
        RowTuples child_annotations = fetched_annotations_map[child];
        
        std::unordered_set<Column> child_annotations_set;

        for (auto & child_anot : child_annotations) {
            child_annotations_set.insert(child_anot.first);
        }

        for (Row & parent_node : parents_vector) {
            RowTuples parent_annotations = fetched_annotations_map[parent_node];

            for (auto & parent_anot : parent_annotations) {
                if (child_annotations_set.count(parent_anot.first) && labels_of_the_start_node.count(parent_anot.first)) 
                    labeled_parents_map_to_top[parent_anot.first][child].push_back(parent_node);
            }
            
        }
    }


    std::unordered_map<Column, std::unordered_map<Row, std::vector<Row>>> labeled_children_map_to_top;

    for (auto & [parent, children_vector] : children_map_to_top) {
        RowTuples parent_annotations = fetched_annotations_map[parent];
        
        std::unordered_set<Column> parent_annotations_set;

        for (auto & parent_anot : parent_annotations) {
            parent_annotations_set.insert(parent_anot.first);
        }

        for (Row & child_node : children_vector) {
            RowTuples child_annotations = fetched_annotations_map[child_node];

            for (auto & child_anot : child_annotations) {
                if (parent_annotations_set.count(child_anot.first) && labels_of_the_start_node.count(child_anot.first)) 
                    labeled_children_map_to_top[child_anot.first][parent].push_back(child_node);
            }
            
        }
    }


    // each Label can have only one target
    std::unordered_map<Column, Row> labeled_targets_to_top;

    for (auto & target : source_targets) {
        RowTuples target_annotations = fetched_annotations_map[target];

        for (auto & [j, tuple] : target_annotations) {
            if (labels_of_the_start_node.count(j))
                labeled_targets_to_top[j] = target;
        }
    }




    std::unordered_map<Column, Row> additional_targets_labeled_to_top;

    for (auto & labd_child_map : labeled_children_map_to_top) {
        Column curr_Label = labd_child_map.first;

        if (labeled_targets_to_top.count(curr_Label)) 
            continue;

        auto curr_child_map = labd_child_map.second;

        uint64_t final_coord = fetched_annotations_map_sep_by_labels[curr_Label][curr_child_map.begin()->first].back();

        Row target_cand = curr_child_map.begin()->first;

        Row curr = target_cand;

        while (curr_child_map.count(curr)) {
            std::vector<Row> currs = curr_child_map[curr];

            uint64_t final_coord_candidate = final_coord - SHIFT;            
            bool coordinates_can_decrease_more = false;

            for (auto & curr_cand : currs) {
                std::vector<uint64_t> curr_cand_coordinates = fetched_annotations_map_sep_by_labels[curr_Label][curr_cand];
                std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(),
                                                                       curr_cand_coordinates.end());

                if (curr_cand_coordinates_set.count(final_coord_candidate)) {
                    final_coord = final_coord_candidate;
                    curr = curr_cand;
                    coordinates_can_decrease_more = true;
                    target_cand = curr;
                    break;
                }
            }

            if (!coordinates_can_decrease_more) {
                target_cand = curr;
                break;
            }
        }

        additional_targets_labeled_to_top[curr_Label] = target_cand;  

    }

    // get final targets
    for (auto & addit_target : additional_targets_labeled_to_top) {
        labeled_targets_to_top[addit_target.first] = addit_target.second;
    }


    // for each Label there is only one trace
    std::unordered_map<Column, std::vector<std::pair<Row, uint64_t>>> traces_to_top_sep_by_labels;

    for (auto & target : labeled_targets_to_top) {
        
        Row target_row = target.second;
        Column curr_Label = target.first;

        Row curr = target_row;
        std::vector<std::pair<Row, uint64_t>> curr_trace;

        uint64_t curr_coord = fetched_annotations_map_sep_by_labels[curr_Label][curr].back(); 

        curr_trace.push_back(std::make_pair(curr, curr_coord));
        
        while (labeled_parents_map_to_top[curr_Label].count(curr)) { 
            
            std::vector<Row> currs = labeled_parents_map_to_top[curr_Label][curr];

            uint64_t curr_coord_candidate = curr_coord + SHIFT;
            bool coordinates_can_increase_more = false;

            for (auto & curr_cand : currs) {
                std::vector<uint64_t> curr_cand_coordinates = fetched_annotations_map_sep_by_labels[curr_Label][curr_cand];
                std::unordered_set<uint64_t> curr_cand_coordinates_set(curr_cand_coordinates.begin(),
                                                                       curr_cand_coordinates.end());

                if (curr_cand_coordinates_set.count(curr_coord_candidate)) {
                    curr_coord = curr_coord_candidate;
                    curr = curr_cand;
                    coordinates_can_increase_more = true;
                    break;
                }
            }

            if (!coordinates_can_increase_more)
                break;            
            
            curr_trace.push_back(std::make_pair(curr, curr_coord));          
        }

        traces_to_top_sep_by_labels[curr_Label] = curr_trace;

    }


    // concatenate all labeled paths from the top to the bottom of the graph
    std::vector<std::pair<std::vector<Row>, Column>> result;

    for (auto & label_of_the_st_node : labels_of_the_start_node) {

        std::vector<std::pair<Row, uint64_t>> path_to_top = traces_to_top_sep_by_labels[label_of_the_st_node];
        std::vector<std::pair<Row, uint64_t>> path_to_bottom = traces_to_bottom_sep_by_labels[label_of_the_st_node];

        if (path_to_top.empty()) {
            
            std::vector<Row> path_to_add_to_result;
            path_to_add_to_result.reserve(path_to_bottom.size());
            for (auto & [row, coordinate] : path_to_bottom) { 
                path_to_add_to_result.push_back(row);
            }

            result.push_back(std::make_pair(path_to_add_to_result, label_of_the_st_node));


        } else if (path_to_bottom.empty()) {

            std::vector<Row> path_to_add_to_result;
            path_to_add_to_result.reserve(path_to_top.size());
            for (auto & [row, coordinate] : path_to_top) { 
                path_to_add_to_result.push_back(row);
            }

            result.push_back(std::make_pair(path_to_add_to_result, label_of_the_st_node));

        } else {

            // if graph contains cycles then paths to the top and to the bottom can overlap
            // so here they are combined
            std::vector<std::pair<Row, uint64_t>> starting_path;
            std::vector<std::pair<Row, uint64_t>> ending_path;

            // find which path has a starting node (i.e. first node with lower coordinate)
            if (path_to_top.front().second < path_to_bottom.front().second)
                starting_path = path_to_top;
            else
                starting_path = path_to_bottom;

            // find which path has an starting node (i.e. last node with bigger coordinate)
            if (path_to_top.back().second > path_to_bottom.back().second)
                ending_path = path_to_top;
            else
                ending_path = path_to_bottom;

            uint64_t last_coord_start_path = starting_path.back().second;
            std::vector<Row> path_to_add_to_result;

            for (auto it = ending_path.rbegin(); it < ending_path.rend(); ++it) {
                if (it->second > last_coord_start_path)
                    path_to_add_to_result.insert(path_to_add_to_result.begin(), it->first);
                else
                    break;
            }

            std::vector<Row> rows_start_path;
            rows_start_path.reserve(starting_path.size());
            for (auto & [row, coordinate] : starting_path) {
                rows_start_path.push_back(row);
            }            

            path_to_add_to_result.insert(path_to_add_to_result.begin(),
                                         rows_start_path.begin(),
                                         rows_start_path.end());

            result.push_back(std::make_pair(path_to_add_to_result, label_of_the_st_node));
        }
    }    

    return result;
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __TUPLE_ROW_DIFF_HPP__
