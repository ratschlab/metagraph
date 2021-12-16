#ifndef __TUPLE_ROW_DIFF_HPP__
#define __TUPLE_ROW_DIFF_HPP__

#include <algorithm>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/vector.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/annotated_dbg.hpp"
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
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    const graph::boss::BOSS &boss = graph_->get_boss();

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
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    const size_t RD_PATH_RESERVE_SIZE = 2;

    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    rd_ids.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;
    node_to_rd.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // Truncated row-diff paths, indexes to |rd_rows|.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<std::pair<size_t, size_t>>> rd_paths_trunc(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        std::vector<std::pair<size_t, size_t>> &rd_path = rd_paths_trunc[i];

        std::vector<size_t> path;
        Vector<std::pair<size_t, Row>> queue;
        queue.emplace_back(0, row_ids[i]);

        while (queue.size()) {
            size_t depth = queue.back().first;
            Row row = queue.back().second;
            queue.pop_back();
            while (depth < path.size()) {
                assert(path.size() > 1);
                rd_path.emplace_back(*(path.rbegin() + 1), *path.rbegin());
                path.pop_back();
            }
            auto [it, is_new] = node_to_rd.try_emplace(row, rd_ids.size());
            path.push_back(it.value());
            // If a node had been reached before, we interrupt the diff path.
            // The annotation for that node will have been reconstructed earlier
            // than for other nodes in this path as well. Thus, we will start
            // reconstruction from that node and don't need its successors.
            if (!is_new)
                continue;

            rd_ids.push_back(row);

            if (anchor_[row])
                continue;

            auto node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row);
            graph_->call_row_diff_successors(node, fork_succ_, [&](auto succ) {
                queue.emplace_back(depth + 1, graph::AnnotatedSequenceGraph::graph_to_anno_index(succ));
            });
        }

        while (path.size() > 1) {
            rd_path.emplace_back(*(path.rbegin() + 1), *path.rbegin());
            path.pop_back();
        }
        assert(path.size());
        rd_path.emplace_back(-1, path[0]);
    }

    node_to_rd = VectorMap<Row, size_t>();

    std::vector<RowTuples> rd_rows = diffs_.get_row_tuples(rd_ids);
    for (auto &row : rd_rows) {
        decode_diffs(&row);
        std::sort(row.begin(), row.end());
    }

    rd_ids = std::vector<Row>();

    // reconstruct annotation rows from row-diff
    std::vector<RowTuples> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        const auto &rd_path = rd_paths_trunc[i];
        // propagate back and reconstruct full annotations for predecessors
        for (size_t j = 0; j + 1 < rd_path.size(); ++j) {
            auto [node, succ] = rd_path[j];
            // reconstruct annotation by adding the diff (full succ + diff)
            add_diff(rd_rows[succ], &rd_rows[node]);
        }
        rows[i] = rd_rows[rd_path.back().second];
        assert(std::all_of(rows[i].begin(), rows[i].end(),
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
                    if (result.back().second.empty())
                        result.pop_back();
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

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __TUPLE_ROW_DIFF_HPP__
